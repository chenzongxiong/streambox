#include <time.h>
#include <boost/lexical_cast.hpp>

#include "date_and_time.hh"

namespace Kaskade
{
//  std::vector<const char[4]> months = {"Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"};
//  //std::vector<const char[3]> numericalMonths = { "01", "02" };
//
//
//  std::string numericalMonth(std::string const& month)
//  {
//    std::string result;
//    for(int i=0; i<months.size(); ++i) if(month == months[i]) result = boost::lexical_cast<std::string>(i);
//    if( result.size()==1 ) result = std::string("0").append(result);
//    return result;
//  }
//
//  std::string replaceMonth(std::string const& str)
//  {
//    std::string tmp;
//    for(int i=0; i<month.size(); ++i)
//    {
//      tmp = months[i]; tmp.append("_");
//      if(str.find(tmp) != std::string::npos) str.replace(str.find(tmp), 3, numericalMonth(months[i]))
//    }
//  }

  std::string getDateAndTime()
  {
    time_t rawtime;   time(&rawtime);
    std::string str(ctime(&rawtime));
    str.erase(str.rfind('\n'),1);
    return str;
  }

  std::string getDate()
  {
    time_t rawtime;     time(&rawtime);
    std::string str(ctime(&rawtime));
    str.replace(str.find(" "), 1, "_");
    if(str.find("  ") != std::string::npos) str.replace(str.find("  "), 2, "_");
    else str.replace(str.find(" "), 1, "_");
    str.erase(str.find(" "));
    return str;
  }

//  std::string getNumDate()
//  {
//    time_t rawtime;     time(&rawtime);
//    std::string str(ctime(&rawtime));
//    str.replace(str.find(" "), 4, "");
//    if(str.find("  ") != std::string::npos) str.replace(str.find("  "), 2, "_");
//    else str.replace(str.find(" "), 1, "_");
//    str.erase(str.find(" "));
//    return str;
//  }

  /// append current time and date to string str
  std::string& appendDateAndTime(std::string& str)
  {
    str.append( getDateAndTime() );
    return str;
  }

  /// append current date to string str
  std::string& appendDate(std::string& str)
  {
    str.append( getDate() );
    return str;
  }
}
