/*
 * PipelineRunner.h
 *
 *  Created on: Jun 21, 2016
 *      Author: xzl
 */

#ifndef PIPELINERUNNER_H_
#define PIPELINERUNNER_H_

#include "Transforms.h"
#include "TransformTreeNode.h"

// an interface
class PipelineRunner {
public:
	static PipelineRunner* fromOptions(PipelineOptions* opt) {
		// avoid InstanceBuilder complication...
		// by design, this should dispatch to concrete pipeline runners
		assert(0);
		return NULL;
	}

	virtual PipelineResult* run (Pipeline* p) = 0;

	PValue* apply(PTransform* p, PValue* input) {
		return p->apply(input);
	}

	virtual ~PipelineRunner() { }
};

////////////////////////////////////////////////////////////////////////


// a minimalist runner. why do we need a separate "evaluator"??
class MinimalPipelineRunner : public PipelineRunner {

public:
	static MinimalPipelineRunner *fromOptions(PipelineOptions *opt) {
		// TODO validate options
		return new MinimalPipelineRunner();
	}

	PipelineResult* run(Pipeline *p) { return NULL; }

private:
	MinimalPipelineRunner() {
		// initialize IO, random seed, etc...
//	  printf("Initializing a %s...\n", __func__);
	}
};


#endif /* PIPELINERUNNER_H_ */
