/*
 * Pipeline.h
 *
 *  Created on: Jun 19, 2016
 *      Author: xzl
 *
 *  NOTE: to avoid circular dep, this header should contain no implementation.
 */

#ifndef PIPELINE_H
#define PIPELINE_H

#include <set>
#include <map>
#include <list>
#include <string>
#include <malloc.h>
#include <assert.h>

class PValue;
class PBegin;

//template <class Derived>
class PTransform;

class PipelineRunner;

class TransformHierarchy;
class TransformTreeNode;
class AppliedPTransform;

class PipelineVisitor;

using namespace std;

enum {
	UNKNOWN,
	STOPPED,
	RUNNING,
	DONE,
	FAILED,
	CANCELLED,
	UPDATED
};

typedef struct state_t {
	bool terminal, hasReplacement;
} State;

static const State states[] = {
	[UNKNOWN] 	= {false, false},
	[STOPPED] 	= {false, false},
	[RUNNING] 	= {false, false},
	[DONE] 		= {true, false},
	[FAILED] 	= {true, false},
	[CANCELLED] = {true, false},
	[UPDATED] 	= {true, true}
};

class PipelineResult {
public:
	State *getState() {
		State* s = (State *)malloc(sizeof(State));
		assert(s);
		s->terminal = _terminal;
		s->hasReplacement = _hasReplacement;
		return s;
	}

private:
	bool _terminal, _hasReplacement;
};

////////////////////////////////////////////////////////////////////////

class  PipelineOptions {
	// TODO
};

////////////////////////////////////////////////////////////////////////
#include <vector>

class PCollection;

class Pipeline {
public:
	static PValue* applyTransform(PValue* input, PTransform* transform);

	static PValue* applyTransform(string name, PValue* input, PTransform* transform);

	static Pipeline* create(PipelineOptions* opt);

	PBegin* begin();

	PValue *apply(PTransform* root);

	PValue *apply1(PTransform* root);

	std::vector<PCollection*> apply2(PTransform* root);

	PipelineRunner *getRunner();

	PipelineOptions *getOptions();

	 /**
	   * Adds the given {@link PValue} to this {@link Pipeline}.
	   *
	   * <p>For internal use only.
	   */
	void addValueInternal(PValue *value);

	PipelineResult* run();

	void traverseTopologically(PipelineVisitor* visitor);

	string toString();

	PBegin*  saved_begin = nullptr;


private:
	PipelineRunner* 			_runner;
	PipelineOptions* 			_options;
	TransformHierarchy*			_transforms;
  set<PValue *>         _values;
	set<string>					_usedFullnames;

	void verifyOutputState(PValue *output,
			TransformTreeNode *node);

	PValue* applyInternal(string name, PValue *input,
			PTransform* transform);

	list<AppliedPTransform* >
		getProducingTransforms(PValue* output);

	string uniquifyInternal(string namePrefix, string origName);

	string buildName(string namePrefix, string name);

protected:
	// ctor
	Pipeline(PipelineRunner* runner, PipelineOptions* options);

};

////////////////////////////////////////////////////////////////////////


#endif // PIPELINE_H
