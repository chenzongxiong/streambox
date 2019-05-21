/*
 * PipelineVisitor.h
 *
 *  Created on: Jun 21, 2016
 *      Author: xzl
 */

#ifndef PIPELINEVISITOR_H_
#define PIPELINEVISITOR_H_


class TransformTreeNode;
class PValue;

// an interface
class PipelineVisitor {
public:
	// shall the visitor process the contents of a composite transform?
	enum CompositeBehavior {
		ENTER_TRANSFORM,
		DO_NOT_ENTER_TRANSFORM
	};

	// note that we not only visit transforms but also values. see the Beam comments

	virtual CompositeBehavior enterCompositeTransform(TransformTreeNode* node) = 0;

	virtual void leaveCompositeTransform(TransformTreeNode* node) = 0;

	virtual void visitPrimitiveTransform(TransformTreeNode* node) = 0;

	virtual void visitValue(PValue* value, TransformTreeNode* producer) = 0;

	virtual ~PipelineVisitor() { }
};

// nop visitor.
//class Defaults : public PipelineVisitor {
class PipelineVisitorNop : public PipelineVisitor {
	virtual CompositeBehavior enterCompositeTransform(TransformTreeNode* node) { return PipelineVisitor::ENTER_TRANSFORM; }

	virtual void leaveCompositeTransform(TransformTreeNode* node) { }

	virtual void visitPrimitiveTransform(TransformTreeNode* node) { }

	virtual void visitValue(PValue* value, TransformTreeNode* producer) { }

};

#endif /* PIPELINEVISITOR_H_ */
