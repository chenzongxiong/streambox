/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2002-2011 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */


#include <iostream>
#include <memory>
#include <string>
#include <vector>

// available under unix only
//#include <sys/utsname.h>

#include <boost/program_options.hpp>
namespace po = boost::program_options;
#include <boost/property_tree/json_parser.hpp>

#include "utilities/kaskopt.hh"

namespace Kaskade
{
  template <class T>
  Options::Options(char const* option, T& parameter, T const& defaultValue, char const* comment)
  : desc(new po::options_description("Allowed options"))
  {
    desc->add_options()("help", "produce help message");
    (*this)(option,parameter,defaultValue,comment);
  }
  
  Options::Options(char const* option, std::string& parameter, char const* defaultValue, char const* comment)
  : Options(option,parameter,std::string(defaultValue),comment)
  {}
  
  template Options::Options(char const* option, bool  & parameter, bool   const& defaultValue, char const* comment);
  template Options::Options(char const* option, int   & parameter, int    const& defaultValue, char const* comment);
  template Options::Options(char const* option, double& parameter, double const& defaultValue, char const* comment);
  template Options::Options(char const* option, std::string& parameter, std::string const& defaultValue, char const* comment);

  template <class T>
  Options& Options::operator()(char const* option, T& parameter, T const& defaultValue, char const* comment)
  {
    desc->add_options()(option,po::value(&parameter)->default_value(defaultValue),comment);
    return *this;
  }
  
  Options& Options::operator()(char const* option, std::string& parameter, char const* defaultValue, char const* comment)
  {
    desc->add_options()(option,po::value(&parameter)->default_value(defaultValue),comment);
    return *this;
  }
  
  template Options& Options::operator()(char const* option, bool  & parameter, bool   const& defaultValue, char const* comment);
  template Options& Options::operator()(char const* option, int   & parameter, int    const& defaultValue, char const* comment);
  template Options& Options::operator()(char const* option, double& parameter, double const& defaultValue, char const* comment);
  template Options& Options::operator()(char const* option, std::string& parameter, std::string const& defaultValue, char const* comment);
  
  boost::program_options::options_description const& Options::descriptions() const
  {
    return *desc;
  }
  
  Options::~Options()
  {}
  
  bool getKaskadeOptions(int argc, char** argv, Options const& options)
  {
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, options.descriptions()), vm);
    po::notify(vm);
    
    if (vm.count("help")) {
      std::cout << options.descriptions() << "\n";
      return true;
    }
    
    return false;
  }
}

// ------------------------------------------------------------------


using namespace std;
using namespace boost::property_tree;

static void copyOver(const ptree *pt, ptree *base, string path) {
  // Value or object or array
  if (!path.empty() && pt->empty())
  {
    // overwrite or add data
    string data = pt->get_value<string>();
    base->put(path, data);
  }
  else if (!path.empty() && pt->count(string()) == pt->size())
  {	
    // overwrite or add array
    ptree::const_iterator it = pt->begin();
    int cnt =  0;
    for (; it != pt->end(); ++it)
    {
      string name;
      name = path + "[" + "cnt" + "]";
      copyOver(&(it->second), base, name);
      cnt++;
    }
  }
  else
  {
    // overwrite or add subtree
    ptree::const_iterator it = pt->begin();
    for (; it != pt->end(); ++it)
    {
      string name;
      if (path.empty())
        name = it->first;
      else
        name  = path + "." + it->first;
      copyOver(&(it->second), base, name);
    }
  }
}

static void sanityCheck(std::unique_ptr<ptree>& pt, int verbosity) {
  string empty;
  string solverType = "names.type." + getParameter(pt, "solver.type", empty), nameType;
  int typeNo = getParameter(pt, solverType, -1), nameNo;
  if (typeNo==-1)
  {
    cout << "Solver type " << getParameter(pt, "solver.type", empty) << " not known" << endl;
  }

  nameType = "names.direct." + getParameter(pt, "solver.direct", empty);
  nameNo = getParameter(pt, nameType, -1);
  if ((nameNo<0)&&(verbosity>0))
  {
    cout << "Direct Solver " << getParameter(pt, "solver.direct", empty) << " not known" << endl;
  }

  nameType = "names.iterate." + getParameter(pt, "solver.iterate", empty);
  nameNo = getParameter(pt, nameType, -1);
  if ((nameNo<0)&&(verbosity>0))
  {
    cout << "Iterate Solver " << getParameter(pt, "solver.iterate", empty) << " not known" << endl;
  }

  string propName = "names.property." + getParameter(pt, "solver.property", empty);
  int propNo = getParameter(pt, propName, -1);
  if ((propNo<0)&&(verbosity>0))
  {
    cout << "Matrix property " << getParameter(pt, "solver.property", empty) << " not known" << endl;
  }

  string preconditioner = getParameter(pt, "solver.preconditioner", empty);
  if (!preconditioner.empty())
  {
    string precondName = "names.preconditioner." + preconditioner;
    int precondNo = getParameter(pt, precondName, -1);
    if ((precondNo<0)&&(verbosity>0))
    {
      cout << "Preconditioner " << getParameter(pt, "solver.preconditioner", empty) << " not known" << endl;
    }
  }
}

static void setRunParameters(ptree& pt) {
  pt.put("run.version", VERSION);
  pt.put("run.kaskade7", KASKADE);

  time_t rawtime;
  struct tm * timeinfo;
  char timeBuffer[256], pidBuffer[256];

  time ( &rawtime );
  timeinfo = localtime ( &rawtime );

  strftime (timeBuffer,256,"%Y-%m-%d-%H-%M",timeinfo);
  string s1(timeBuffer);
  pt.put("run.starttime", s1);
  sprintf(pidBuffer,"%d", getpid());
  string s2(pidBuffer);
  pt.put("run.pid", s2);

  // the following lines only work under unix systems !!!!!!!!!!!!
  /*    struct utsname name;
   * int rc = uname(&name);
   * if (rc==0)
   *   {
   *     string sysname(name.sysname);
   *     pt.put("run.sysname", sysname);
   *     string nodename(name.nodename);
   *     pt.put("run.nodename", nodename);
   *     string release(name.release);
   *     pt.put("run.release", release);
   *     string machine(name.machine);
   *     pt.put("run.machine", machine);
   }*/
}

static bool checkFile(const string &s) {
  FILE *f = fopen(s.c_str(), "r");
  if (f) fclose(f);
  return f;
}

std::unique_ptr<ptree> getKaskadeOptions(int argc, char *argv[], int verbosity, bool dump)
{
  using boost::program_options::value;

  std::unique_ptr<ptree> pt(new ptree);
  string optionFile(KASKADE), empty;
  if (checkFile("default.json"))
    optionFile = "default.json";
  else
    optionFile += "/default.json";
  try
  {
    boost::program_options::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("default", value<string>(), "select default options")
        ("addons", value<string>(), "addional options")
        ("verbose", value<int>(), "verbosity level")
        ;

    boost::program_options::variables_map vm;
    boost::program_options::parsed_options parsed = boost::program_options::command_line_parser(argc, argv).options(desc).allow_unregistered().run();
    store(parsed, vm);
    notify(vm);

    if (vm.count("help"))
    {
      cout << desc << "\n";
      exit(0);
    }
    if (vm.count("default"))
    {
      string personalOptionFile(vm["default"].as<string>());
      read_json(personalOptionFile, *pt);
      if (verbosity>0)
        cout << "Default parameters read from " << personalOptionFile << endl;
    }
    else
    {
      read_json(optionFile, *pt);
      if (verbosity>0)
        cout << "Default parameters read from " << optionFile << endl;
    }
    if (vm.count("addons"))
    {
      string addionalOptionFile(vm["addons"].as<string>());
      ptree new_pt;
      if (!checkFile(addionalOptionFile))
      {
        string dummy(KASKADE);
        addionalOptionFile = dummy + "/" + addionalOptionFile;
        if (verbosity>0)
          cout << "Addional parameters read from " << addionalOptionFile << endl;
      }
      read_json(addionalOptionFile, new_pt);
      if (verbosity>0)
        cout << "Addional parameters read from " << addionalOptionFile << endl;
      string path;
      copyOver(&new_pt, pt.get(), path);
    }
    if (vm.count("verbose"))
    {
      verbosity = vm["verbose"].as<int>();
    }
    pt->put("verbose", verbosity);
    vector<string> addional_parameters = collect_unrecognized(parsed.options, boost::program_options::include_positional);
    int k, lng = addional_parameters.size();
    string val;
    for (k=0; k<lng; k++)
    {
      if ((addional_parameters[k][0]=='-')&&(addional_parameters[k][1]=='-'))
      {
        addional_parameters[k].erase(0,2);
        val = getParameter(pt, addional_parameters[k], empty);
        if (k==lng-1)
        {
          if (verbosity>0)
            cout << "Named parameter("  << k+1 << ") --" << addional_parameters[k] <<
            " needs value" << endl;
          break;
        }
        if (val.empty())
        {
          if (verbosity>0)
            cout << "New parameter " << addional_parameters[k] << " added " <<
            addional_parameters[k+1] << endl;
          pt->put(addional_parameters[k], addional_parameters[k+1]);
          k++;
        }
        else
        {
          if (verbosity>0)
            cout << "Parameter " << addional_parameters[k] << "=" << val <<
            " changed to " << addional_parameters[k+1] << endl;
          pt->put(addional_parameters[k], addional_parameters[k+1]);
          k++;
        }
      }
      else
      {
        if (verbosity>0)
          cout << "Unknown posional parameter("  << k << ") " << addional_parameters[k] << endl;
      }
    }
  }
  catch(exception& e)
  {
    cerr << "error: " << e.what() << "\n";
    exit(0);
  }
  catch(...)
  {
    cerr << "Exception of unknown type!\n";
    exit(0);
  }

  setRunParameters(*pt);
  sanityCheck(pt, verbosity);

  if (dump)  {
    string outFile("run-");
    outFile += getParameter(pt, "run.starttime", empty)+"-"+ getParameter(pt, "run.pid", empty) + ".json";
    write_json(outFile, *pt);
  }
  return pt;
}
