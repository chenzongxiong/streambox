/*
 * TransformTreeNode.h
 *
 *  Created on: Jun 21, 2016
 *      Author: xzl
 *
 *  moved out of Transforms.h, which is supposed to be an interface (included in Values.h).
 */

#ifndef TRANSFORMTREENODE_H_
#define TRANSFORMTREENODE_H_

#include <assert.h>
#include <string>
#include <map>

using namespace std;

class TransformTreeNode {
private:
	bool _finishedSpecifying;
	PTransform * _transform;
	string _fullName;
	TransformTreeNode* _enclosingNode;

	// Nodes for sub-transforms of a composite transform.
	set<TransformTreeNode* > _parts;

	// expanded form
	map<PValue*, TransformTreeNode*> _inputs;
	// unexpanded form
	PValue *_input;
	PValue *_output;

public:
	TransformTreeNode(TransformTreeNode* enclosingNode, PTransform * transform,
			string fullName, PValue* input) :
				_finishedSpecifying(false), _transform(transform),
			_fullName(fullName), _enclosingNode(enclosingNode), _input(input), _output(NULL)
	 { assert((!enclosingNode && !transform) || (enclosingNode && transform)); }

	PTransform * getTransform() { return _transform; }

	TransformTreeNode* getEnclosingNode() { return _enclosingNode; }

	// add a composite operation to the transform node
	void addComposite(TransformTreeNode* node) { _parts.insert(node); }

	bool isCompositeNode() {
		return (!_parts.empty() || returnsOthersOutput() || isRootNode());
	}

	bool isRootNode() { return (!_transform); }

	string getFullName() { return _fullName; } // copy

	// add an input to the transform node
	void addInputProducer(PValue* expandedInput, TransformTreeNode* producer) {
		assert(!_finishedSpecifying);
		_inputs.insert(pair<PValue *, TransformTreeNode *>(expandedInput, producer));
	}

	PValue* getInput() { return _input; }

	map<PValue*, TransformTreeNode*>* getInputs()  { return &_inputs; }

	void setOutput (PValue *output) {
		assert(!_finishedSpecifying);
		assert(!(this->_output));
		_output = output;
	}

	// XXX should we do copy return?
	list<PValue *>* getExpandedOutputs() {
		if (_output)
			return (_output->expand());
		else
			return new list<PValue *>(); // nothing to expand
	}

	PValue *getOutput() { return _output; }

  /**
   * Visit the transform node.
   *
   * <p>Provides an ordered visit of the input values, the primitive
   * transform (or child nodes for composite transforms), then the
   * output values.
   */
	void visit(PipelineVisitor* visitor, set<PValue *>* visitedValues) {
		if (!_finishedSpecifying) {
			finishSpecifying();
		}

		// visit inputs
		for (auto const& ent : _inputs) {
			if ((visitedValues->insert(ent.first)).second)  // insertion okay
				visitor->visitValue(ent.first, ent.second);
		}

		// visit transform
		if (isCompositeNode()) {
			PipelineVisitor::CompositeBehavior recurse = visitor->enterCompositeTransform(this);

			if (recurse == PipelineVisitor::CompositeBehavior::ENTER_TRANSFORM) {
				for (auto const& child : _parts) {
					child->visit(visitor, visitedValues);
				}
			}
			visitor->leaveCompositeTransform(this);
		} else
			visitor->visitPrimitiveTransform(this);

		// visit outputs (each of them)
		for (auto const& ent : *getExpandedOutputs()) {
			if ((visitedValues->insert(ent)).second)  // insertion okay
				visitor->visitValue(ent, this);
			// XXX: should we free the result of getExpandedOutputs()?
		}
	}

  /**
   * Finish specifying a transform.
   *
   * <p>All inputs are finished first, then the transform, then
   * all outputs.
   */
	void finishSpecifying() {
		if (_finishedSpecifying)
			return;
		_finishedSpecifying =  true;

		for (auto const& ent1 : _inputs) {
			if (ent1.second)
				ent1.second->finishSpecifying();
		}

		if (_output)
			_output->finishSpecifyingOutput();
	}

private:
	bool returnsOthersOutput() {
		return false; // TODO -- need better unerstanding
	}


};

////////////////////////////////////////////////////////////////////////

/**
 * Captures information about a collection of transformations and their
 * associated {@link PValue}s.
 *
 * xzl: "hierarchy" is the key...
 *
 * this seems for constructing the pipeline.the construction procedure
 * frequently pushes & pops the stack.
 *
 * so when pipeline construction
 * is done, the only node on the stack is the root node.
 *
 */
class TransformHierarchy {
private:
	list<TransformTreeNode *> _transformStack;  //=  new auto();
	map<PValue *, TransformTreeNode *> _producingTransformNode;  //=  new auto();

public:
	TransformHierarchy() {
		_transformStack.push_back(new TransformTreeNode(NULL, NULL, "", NULL));
	}

	TransformTreeNode *getCurrent() {
		return _transformStack.back();
	}

	void pushNode(TransformTreeNode* current) {
		_transformStack.push_back(current);
	}

	void popNode() {
		_transformStack.pop_back();
		assert(!_transformStack.empty());
	}

  /**
   * Adds an input to the given node.
   *
   * <p>This forces the producing node to be finished.
   * xzl: keeps track of all input PValues and their producing TransformTreeNode
   */
	void addInput(TransformTreeNode* node, PValue* input) {
		for (auto const& i : *input->expand()) {
			auto it = _producingTransformNode.find((PValue *)i);
			assert (it != _producingTransformNode.end()); // unknown input
			TransformTreeNode *producer = (*it).second;
			producer->finishSpecifying();
			node->addInputProducer(i, producer);
		}
	}

	void setOutput(TransformTreeNode *producer, PValue *output) {
		producer->setOutput(output);

		for (auto const& o : *output->expand()) {
			// xzl: expand() will return a list of PValues.
			// the typecast may fail if not PValue*
			_producingTransformNode.insert(
					pair<PValue*, TransformTreeNode*>((PValue *)o, producer));
		}
	}

  /**
   * Visits all nodes in the transform hierarchy, in transitive order.
   *
   * xzl: start from the root node, which should be the only node on the
   * stack.
   */
	void visit(PipelineVisitor* visitor, set<PValue *>* visitedNodes) {
		_transformStack.back()->visit(visitor, visitedNodes);
	}
};


#endif /* TRANSFORMTREENODE_H_ */
