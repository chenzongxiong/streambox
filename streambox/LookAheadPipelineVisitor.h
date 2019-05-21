#ifndef LOOKAHEAD_PIPELINE_H
#define LOOKAHEAD_PIPELINE_H

#include "core/PipelineVisitor.h"

class EvaluationBundleContext;

/*
 * set up connections among PValue and PTransform
 */
class LookAheadPipelineVisitor : public PipelineVisitorNop {
public:
  // dispatch transforms to their respective evaluators.
  void visitPrimitiveTransform(TransformTreeNode* node) override;
  void visitValue(PValue *value, TransformTreeNode *producer) override { }
  // ctor
  LookAheadPipelineVisitor(EvaluationBundleContext *ctx);
private:
  EvaluationBundleContext* _ctx;
};

#endif // LOOKAHEAD_PIPELINE_H
