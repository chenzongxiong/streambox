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

#include <cassert>
#include <iomanip>

#include <boost/timer/timer.hpp>

#include "utilities/timing.hh"

namespace Kaskade
{

  struct Timings::Times
  {
    std::string             name;
    size_t                  calls;
    std::vector<Times>      children;
    boost::timer::cpu_timer timer;

    int maxNameLength(int indent) const
    {
      int len = name.size();
      for (auto const& c: children)
        len = std::max(len,indent+c.maxNameLength(2));
      return len;
    }

    void report(std::ostream& out, int indent, int pos, boost::timer::nanosecond_type full) const
    {
      out << std::setw(indent) << ""
          << std::setw(3) << std::right << 100*timer.elapsed().wall/full << "% "
          << std::setw(pos-indent) << std::left << name
          << std::setw(6) << std::right << calls << "  ";
      out << timer.format();

      full = timer.elapsed().wall;
      for (auto const& c: children)
        c.report(out,indent+2,pos,full);
    }
  };

  Timings& Timings::instance()
  {
    static Timings globalTimings;
    return globalTimings;
  }

  Timings::Timings()
  : all(new Times{"all",0,{},{}})
  {
    stack.push(all.get());
    all->timer.start();
  }

  std::ostream& Timings::report(std::ostream& out) const
  {
    int pos = all->maxNameLength(2)+4;
    out << std::setw(pos+4) << " " << "  calls   " << "\n";
    //~ out << std::setw(pos+30) << " " << "\n";
    all->report(out,0,pos,all->timer.elapsed().wall);
    return out;
  }

  void Timings::start(std::string const& name)
  {
    Times* section = nullptr;
    for (auto& c: stack.top()->children)
      if (name == c.name)
      {
        section = &c;
        section->timer.resume();
        break;
      }

    if (! section)
    {
      // section not found - create a new one. Note that this can invalidate pointers into
      // the children container. Fortunately, there are no pointers/iterators around: The
      // stack contains only pointers into higher level containers.
      stack.top()->children.push_back(Times{name,0,{},{}});
      section = &(stack.top()->children.back());
      section->timer.start();
    }

    ++section->calls;
    stack.push(section);
  }

  void Timings::stop(std::string const& name)
  {
    assert(stack.top()->name == name);
    stack.top()->timer.stop();
    stack.pop();
  }

  void Timings::clear()
  {
    while (stack.size()>1)  // remove everything (except the root)
      stack.pop();
    all->timer.start();      // reset timer to zero
    all->calls = 0;
    all->children.clear();
  }


  std::ostream& operator<<(std::ostream& out, Timings const& timings)
  {
    return timings.report(out);
  }

}
