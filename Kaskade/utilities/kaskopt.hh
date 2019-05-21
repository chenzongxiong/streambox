/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2010-2012 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef KASKOPT_HH
#define KASKOPT_HH

#include <memory>
#include <string>

#include <boost/property_tree/ptree.hpp>

/// \internal
namespace boost { namespace program_options { class options_description; } }
/// \endinternal

namespace Kaskade
{
  class Options
  {
  public:
    template <class T>
    Options(char const* option, T& parameter, T const& defaultValue, char const* comment);
    Options(char const* option, std::string& parameter, char const* defaultValue, char const* comment);
    
    template <class T>
    Options& operator()(char const* option, T& parameter, T const& defaultValue, char const* comment);
    Options& operator()(char const* option, std::string& parameter, char const* defaultValue, char const* comment);
    
    boost::program_options::options_description const& descriptions() const;
    
    ~Options();
    
  private:
    std::unique_ptr<boost::program_options::options_description> desc;
  };
  
  /**
   * \ingroup utilities
   * \brief Supports the extraction of program option values from command line arguments.
   * 
   * Based on boost::program_options, command line arguments (or their default values) can be stored in 
   * variables as follows:
   * \code
   * bool a;
   * int b;
   * getKaskadeOptions(argc,argv,Options
   * ("a",a,false,"defines boolean parameter a")
   * ("b",b,42,   "defines integer parameter b"));
   * \endcode
   * 
   * If an argument "--help" is given, the list of allowed options is printed.
   * 
   * \return true if an argument "--help" has been provided, otherwise false.
   */
  bool getKaskadeOptions(int argc, char** argv, Options const& options);
}




/**
 * \brief 
 * 
 * The parameters can be written to a log file named "run-<starttime>-<pid>.json" in JSON format.
 * This is done if \arg dumpParameter is true. starttime and pid have to be provided by the
 * run.starttime and run.pid parameters.
 * 
 * \param argc first argument of main to be passed on
 * \param argv second argument of main to be passed on
 * \param verbosity
 * \param dumpParameter if true, the parameters are written to a log file
 */
std::unique_ptr<boost::property_tree::ptree> getKaskadeOptions(int argc, char *argv[], int verbosity, bool dumpParameter);


/**
 * \brief Extracts parameters from the parameter set.
 * \param pt the parameter set as obtained from getKaskadeOptions.
 * \param s the name of the parameter
 * \param defaultValue is returned if the parameter is not found in the parameter set
 */
template <class TYPE>
TYPE getParameter(std::unique_ptr<boost::property_tree::ptree> const& pt, const std::string s, TYPE defaultValue) 
{
    return (*pt).get(s,defaultValue);
}

#endif
