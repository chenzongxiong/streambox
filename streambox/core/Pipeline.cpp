/*
 * Pipeline.cpp
 *
 *  Created on: Jun 21, 2016
 *      Author: xzl
 *
 *  separate it from the header because otherwise we will have circular dependency with Values.h
 */

#define K2_DEBUG_INFO

#include <assert.h>
#include <list>
#include <algorithm>

#include "Pipeline.h"
#include "PipelineRunner.h"
#include "TransformTreeNode.h"
#include "Values.h"
#include "Transforms.h"
#include "log.h"

using namespace std;

PValue* Pipeline::applyTransform(PValue* input, PTransform* transform) {
	return input->getPipeline()->applyInternal(transform->getName(), input, transform);
}

PValue* Pipeline::applyTransform(string name, PValue* input, PTransform* transform) {
	return input->getPipeline()->applyInternal(name, input, transform);
}

Pipeline* Pipeline::create(PipelineOptions* opt) {
	Pipeline *p;
	if (!opt) // xzl: dirty hacking: shortcut to our minimal runner
		p = new Pipeline(MinimalPipelineRunner::fromOptions(opt), opt);
	else
		p = new Pipeline(PipelineRunner::fromOptions(opt), opt);
	return p;
}

PBegin* Pipeline::begin() {
	return PBegin::in(this);
}

PValue* Pipeline::apply(PTransform* root) {
	return begin()->apply(root);
}

PValue* Pipeline::apply1(PTransform* root) {
  assert(saved_begin == nullptr);
  saved_begin = begin();
  return saved_begin->apply1(root);
}

vector<PCollection*> Pipeline::apply2(PTransform* root) {
  assert(saved_begin == nullptr);
  saved_begin = begin();
  return saved_begin->apply2(root);
}

PipelineRunner* Pipeline::getRunner() { return _runner; }

PipelineOptions* Pipeline::getOptions() { return _options; }


 /**
   * Adds the given {@link PValue} to this {@link Pipeline}.
   *
   * <p>For internal use only.
   */
void Pipeline::addValueInternal(PValue *value) {
	this->_values.insert(value);
}

string Pipeline::uniquifyInternal(string namePrefix, string origName) {
	string name = origName;
	int suffixNum = 2;
	while (true) {
		string candidate = buildName(namePrefix, name);

		VV("   candidate is %s\n", candidate.c_str());

		auto ret = _usedFullnames.insert(candidate);
		if (ret.second)
			return candidate;
		else if (suffixNum >= 100) {
			printf("warning: bug? duplicate type name.\n");
			assert(0);
		}
		name = origName + to_string(suffixNum++);
	}

	return name;
}

string Pipeline::buildName(string namePrefix, string name) {
//	cout << "name is " << name << endl;
	return namePrefix.empty() ? name : namePrefix + "/" + name;
}

/**
 * Applies a {@link PTransform} to the given {@link PInput}.
 *
 * @see Pipeline#apply
 */
PValue* Pipeline::applyInternal(string name, PValue *input, PTransform* transform) {

	VV("%s: PTransform 0x%lx actual type %s \n",
	    __func__, (unsigned long)(transform), typeid(*transform).name());

	input->finishSpecifying();

	TransformTreeNode *parent = _transforms->getCurrent();
	string namePrefix = parent->getFullName();
	string fullName = uniquifyInternal(namePrefix, name);

	// TODO -- check if name is unique...

	TransformTreeNode* child = new TransformTreeNode(parent, transform, fullName, input);
	parent->addComposite(child); // xzl: does this mean that we at least have one node (even for primitive?)

	_transforms->addInput(child, input);

	// TODO -- trace

	_transforms->pushNode(child);
	transform->validate(input);
	// xzl: the following can go recursive
	PValue *output = _runner->apply(transform, input);
	_transforms->setOutput(child, output);

	AppliedPTransform* applied = AppliedPTransform::of(child->getFullName(), input, output, transform);
	// recordAsOutput is a NOOP if already called;
	output->recordAsOutput(applied);
	verifyOutputState(output, child);

	_transforms->popNode();

	return output;
}

/**
* Returns all producing transforms for the {@link PValue PValues} contained
* in {@code output}.
*
* xzl: return a copy of the list... okay?
*/
list<AppliedPTransform* > Pipeline::getProducingTransforms(PValue* output) {
	list<AppliedPTransform* > producingTransforms;
	for (auto const & value : *output->expand()) {
		AppliedPTransform *transform =
				value->getProducingTransformInternal();
		if (transform)
			producingTransforms.push_back(transform);
	}
	return producingTransforms; // return a copy
}

  /**
   * Verifies that the output of a {@link PTransform} is correctly configured in its
   * {@link TransformTreeNode} in the {@link Pipeline} graph.
   *
   * <p>A non-composite {@link PTransform} must have all
   * of its outputs registered as produced by that {@link PTransform}.
   *
   * <p>A composite {@link PTransform} must have all of its outputs
   * registered as produced by the contained primitive {@link PTransform PTransforms}.
   * They have each had the above check performed already, when
   * they were applied, so the only possible failure state is
   * that the composite {@link PTransform} has returned a primitive output.
   */
void Pipeline::verifyOutputState(PValue *output, TransformTreeNode *node) {
	if (!node->isCompositeNode()) {
		PTransform* thisTransform = node->getTransform();
		list<AppliedPTransform* > producingTransforms = getProducingTransforms(output);
		for (auto const& producingTransform : producingTransforms) {
			//assert (thisTransform == producingTransform->getTransform());
			if(thisTransform != producingTransform->getTransform()){
				std::cout << "ERROR: thisTransform != producingTransform->getTransform()" << std::endl;
				abort();
			}
		}
	} else { // composite node
		PTransform* thisTransform = node->getTransform();
		list<AppliedPTransform* > producingTransforms = getProducingTransforms(output);
		for (auto const& producingTransform : producingTransforms) {
			//assert (thisTransform != producingTransform->getTransform());
			if(thisTransform == producingTransform->getTransform()){
				std::cout << "thisTransform == producingTransform->getTransform()" << std::endl;
				abort();
			}
		}
	}
}

PipelineResult* Pipeline::run() {
	// XXX tracing
	return _runner->run(this);
}

void Pipeline::traverseTopologically(PipelineVisitor* visitor) {
	set<PValue* > visitedValues; 	// XXX okay on stack?
	// Visit all the transforms, which should implicitly visit all the values.
	// this invokes TransformHierarchy. why??
	_transforms->visit(visitor, &visitedValues);

	if (!includes(visitedValues.begin(), visitedValues.end(),
			_values.begin(), _values.end())) {
		printf("internal error: should have visited all the values"
				"after visiting all the transforms");
		assert(0);
	}
}

string Pipeline::toString() {
	return "Pipeline#";  // TODO
}

// ctor
Pipeline::Pipeline(PipelineRunner* runner, PipelineOptions* options)
		: _runner(runner), _options(options),
		  _transforms(new TransformHierarchy())
		  { }
