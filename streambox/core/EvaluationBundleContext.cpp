/*
 * EvaluationBundleContext.cpp
 *
 *  Created on: Sep 25, 2016
 *      Author: xzl
 *
 *  Purdue University, 2016
 */


#include "EvaluationBundleContext.h"

#if 0
long EvaluationBundleContext::Executor::nwait_ = 0;
std::mutex EvaluationBundleContext::Executor::cv_mtx_;
std::condition_variable EvaluationBundleContext::Executor::cv_;
bool EvaluationBundleContext::Executor::isStop = false;
#endif

/* create async tasks for evaluating downstream transforms.
 *
 * @nodeid: the NUMA node where the new task will run.
 * -1: the task will choose the node that has most pending bundles
 *
 * */
void EvaluationBundleContext::SpawnConsumer(PValue *v, int nodeid)
{
	assert(false && "obsoleted");

   assert(v);

   PTransform *transform = v->consumer;

    //std::cout << "++ in " << transform->getName() << " SpawnConsumer" << std::endl;

   if (!transform) {
       // could be the last PValue in the pipeline. No downstream transform.
  	 assert(0);
     return;
   }

   assert(nodeid != -1); /* XXX want to be cautious about this as of now */

   /* if node unspecified, spawn consumer on the node with most bundles.
    * caveat: we may spawn too many consumers on one node in a "wave",
    * which may later cause load imbalance. */
   if (nodeid == -1) {
  	 int maxbundlecount = -1;
		 for (int i = 0; i < MAX_NUMA_NODES; i++) {
			 if (maxbundlecount < v->size(i)) {
				 maxbundlecount = v->size(i);
				 nodeid = i;
			 }
		 }
   }

   ASSERT_VALID_NUMA_NODE(nodeid);

   VV("--> Eval: %s on node %d", transform->getName().c_str(), nodeid);

//   spawn_cnt_[nodeid].fetch_add(1);
   spawn_cnt_[nodeid] += 1; /* atomic */

#ifdef USE_NUMA_TP
		this->pool.runOnNode(nodeid, Task([this, transform, nodeid] {
			/* need to catch -- otherwise tasks may silently die? */
			//			transform->ExecEvaluator(nodeid, this);
#if 1
			try {
				transform->ExecEvaluator(nodeid, this);

				/* for testing -- will throw system error */
//				mutex mtx;
//				unique_lock<mutex> lk(mtx);
//				lk.unlock();
//				lk.unlock();
			} catch (const std::exception& e) {
        std::cout << " a standard exception was caught, with message '"
                  << e.what() << "'\n";
        std::cout << "Type:    " << typeid(e).name() << "\n";
        abort();
			}
#endif
		}));
#else
		this->tp.push([this, transform, nodeid](int id) {
			try {
				//std::cout << " ===== call donwstream " << transform->getName() << "'s ExecEvaluator" << std::endl;
				transform->ExecEvaluator(nodeid, this);
			} catch (const std::exception& e) {
        std::cout << " a standard exception was caught, with message '"
                  << e.what() << "'\n";
        abort();
			}
		});
#endif

}

/* executes evaluation synchronously.
 *
 * although executing eval sync, we still pass in the node id so that
 * the eval knows which node it is on.
 * */
void EvaluationBundleContext::RunConsumer(PValue *v, int nodeid = -1)
{
  assert(v);

  PTransform *transform = v->consumer;

  if (!transform) {
      // could be the last PValue in the pipeline. No downstream transform.
      return;
  }

  // Dispatch to evaluators based on the transform types
  // [no idea how this can be done dynamically. CRTP won't work:
  // it will make @PTransform a template, bad. ]
  VV("Eval: %s ", transform->getName().c_str());

  transform->ExecEvaluator(nodeid, this);
}

/* Execute synchronously, we execute from the "current" NUMA node.
 *
 * XXX should dispatch work to NUMA nodes considering locality and load
 * balance?? */
//#include <numa.h>
//void EvaluationBundleContext::OnNewUpstreamWatermark(ptime wm,
//    PTransform *transform)
//{
//	int cpu = sched_getcpu();
//	int nodeid = numa_node_of_cpu(cpu);
//
//	transform->OnNewUpstreamWatermark(nodeid, this, wm);
//}
//
