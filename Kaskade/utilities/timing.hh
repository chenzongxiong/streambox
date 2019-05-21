/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2016-2016 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef TIMING_HH
#define TIMING_HH

#include <memory>
#include <ostream>
#include <stack>
#include <vector>


namespace Kaskade
{
  /**
   * \ingroup utilities
   * \brief Supports gathering and reporting execution times information for nested program parts.
   */
  class Timings
  {
  public:
    /**
     * \brief Returns a reference to a singleton instance.
     */
    static Timings& instance();
    
    /**
     * \brief Constructor
     */
    Timings();
    
    /**
     * \brief Prints a timing report to the given stream.
     */
    std::ostream& report(std::ostream& out) const;
    
    /**
     * \brief Starts or continues the timing of given section.
     * 
     * If the section has not been started before, it is created.
     */
    void start(std::string const& name);
    
    /**
     * \brief Stops the timing of given section.
     */
    void stop(std::string const& name);
    
    /**
     * \brief Resets the timer to an empty state.
     */
    void clear();
    
  private:
    struct Times;
    
    std::unique_ptr<Times> all;
    std::stack<Times*> stack;
  };
  
  // Common interface for stream output
  std::ostream& operator<<(std::ostream& out, Timings const& timings);
}

#endif
