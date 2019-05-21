/*
 * ParDo.cpp
 *
 *  Created on: Jun 30, 2016
 *      Author: xzl
 *
 *  see the documents/comments in beam's ParDo.java
 */
#include <string>
#include <vector>
#include "boost/date_time/posix_time/posix_time.hpp"
#include "Values.h"
#include "core/ExecutionContext.h"
#include "DoFn.h"
#include "WindowFn.h"

using namespace std;


/////////////////////////////////////////////////////////////////////

class DoFnRunners {
	// xzl: unclear whether we need this ...

	template <class T>
	class OutputManager {
	public:
		// Outputs a single element to the receiver (??)
		virtual void output(WindowedValue<T>* output) = 0;
		virtual ~OutputManager() { }
	};
}

/////////////////////////////////////////////////////////////////////

// used in evaluating ParDo
template <class InputT, class OutputT>
class DoFnRunnerBase {
private:

  //////////////////////////////////////////

  template <class T>
  class ListOutputManager : public DoFnRunners::OutputManager {
  private:
    vector<WindowedValue<T>> _outputs; // saved outputs
  public:
    // store the output value in an internal list
    void output(WindowedValue<T> o) {
      _outputs.push_back(o);
    }
    vector<WindowedValue<T>>* getOutput() { return &_outputs; }
  };

  //////////////////////////////////////////

	// default implementation of dofn::Context
	template <class InputT, class OutputT>
	class DoFnContext : public Context {
	private:
		PipelineOptions* _opt;
		DoFn<InputT, OutputT>* _fn;
		WindowFn * _windowfn;
		ListOutputManager * _outputManager;
		StepContext * _stepContext;

	public:
		DoFnContext(PipelineOptions *opt, DoFn<InputT, OutputT>* fn,
		    ListOutputManager *outputmanager, StepContext *stepcontext,
				WindowFn *windowfn) :
					_opt(opt), _fn(fn), _outputManager(outputmanager),
					_stepContext(stepcontext), _windowfn(windowfn) { }

		WindowedValue<InputT>* makeWindowedValue(OutputT output, ptime ts,
			vector<W *>* windows, PaneInfo* pane) {
			if (!windows) {
				assert(0); // unimplemented
			}
			return WindowedValue<InputT>::of(output, ts, windows, *pane);
		}

		void outputWindowedValue(WindowedValue<OutputT> elem) {
		  _outputManager->output(elem);
		  if (_stepContext)
		    _stepContext->noteOutput(elem);
		}
	};

  //////////////////////////////////////////

	// for a single value only?
	template <class InputT, class OutputT>
	class DoFnProcessContext
	    : public DoFn<InputT, OutputT>::ProcessContext {

	  DoFn<InputT, OutputT> _fn;
	  DoFnContext<InputT, OutputT> _context;
	  WindowedValue<InputT> _windowedValue;

	};

  DoFn<InputT, OutputT>* _fn;
};

/////////////////////////////////////////////////////////////////////
