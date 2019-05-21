#include "core/Transforms.h"
#include "core/TransformTreeNode.h"
#include "core/EvaluationBundleContext.h"
#include "LookAheadPipelineVisitor.h"

LookAheadPipelineVisitor::LookAheadPipelineVisitor(EvaluationBundleContext *ctx)
	: _ctx (ctx) { }

// build up the pipeline topology.
void LookAheadPipelineVisitor::visitPrimitiveTransform(TransformTreeNode* node) {
	PTransform * t = node->getTransform();
	PValue * in = node->getInput();
	PValue * out = node->getOutput();
	assert(t);

	in->consumer = t;
	out->producer = t;

	t->inputs.push_back(in);
	t->outputs.push_back(out);

#if 1
	if (!_ctx->source) {
		_ctx-> source = t;
	}
#endif

	// XXX we assume only one "begin" node
	if (!_ctx->begin) {
		_ctx->begin = in;
	}

	// XXX assuming the pipeline only has a linear topology
	printf("%s --> %s --> ",
			in->getName().c_str(),
			t->getName().c_str()
	      );
}
