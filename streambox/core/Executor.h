/* OBSOLETED -- moved to Eval context.
 *
 * Executor.h
 *
 * the customized, thread pool based, execution engine for streaming pipeline.
 *
 * We use pthread in case we may upgrade to NUMA-aware tp later.
 *
 * Part of code comes from ctpl. See below.
 *
 *  Created on: Nov 14, 2016
 *      Author: xzl
 *
 *  Purdue University, 2016
 */

/*********************************************************
*
*  Copyright (C) 2014 by Vitaliy Vitsentiy
*
*  Licensed under the Apache License, Version 2.0 (the "License");
*  you may not use this file except in compliance with the License.
*  You may obtain a copy of the License at
*
*     http://www.apache.org/licenses/LICENSE-2.0
*
*  Unless required by applicable law or agreed to in writing, software
*  distributed under the License is distributed on an "AS IS" BASIS,
*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
*  See the License for the specific language governing permissions and
*  limitations under the License.
*
*********************************************************/

#ifndef EXECUTOR_H_
#define EXECUTOR_H_

#include <assert.h>
#include <pthread.h>

#include <vector>

#include "config.h"
#include "log.h"
#include "EvaluationBundleContext.h"

using namespace std;

class EvaluationBundleContext;

class Executor {

  namespace detail {
      template <typename T>
      class Queue {
      public:
          bool push(T const & value) {
              std::unique_lock<std::mutex> lock(this->mutex);
              this->q.push(value);
              return true;
          }
          // deletes the retrieved element, do not use for non integral types
          bool pop(T & v) {
              std::unique_lock<std::mutex> lock(this->mutex);
              if (this->q.empty())
                  return false;
              v = this->q.front();
              this->q.pop();
              return true;
          }
          bool empty() {
              std::unique_lock<std::mutex> lock(this->mutex);
              return this->q.empty();
          }

          // xzl: debugging
          int size() {
          	std::unique_lock<std::mutex> lock(this->mutex);
          	return this->q.size();
          }

      private:
          std::queue<T> q;
          std::mutex mutex;
      };
  }

public:
	Executor (EvaluationBundleContext* c, int nthreads)
		: c_(c), threads_(nthreads) {
		assert(nthreads > 0 && nthreads < 200);
	}

	/* fire up all threads */
	bool StartThreads() {

		for (unsigned int i = 0; i < threads_.size(); i++) {
			int rc = pthread_create(&threads_[i], NULL, run, c_);
			assert(rc == 0);
		}
		E("%lu threads fired up...", threads_.size());

		return true;
	}

	void StopThreads() {
		assert(false && "not implemented");
	}

	/* after producing work, call this to notify any waiting threads */
	static bool Notify() {
		std::unique_lock<mutex> lck(cv_mtx_);
		/* XXX increase a global work counter? XXX */
		cv_.notify_one();
	}

private:

	static std::condition_variable cv_;
	static std::mutex cv_mtx_; 				/* only for protecting the cv */

	static void * run(void *t) {

		EvaluationBundleContext *c = (EvaluationBundleContext *)t;
		shared_ptr<BundleBase> bundleptr;
		PTransform *trans;

		std::function<void(int id)> * _f;
		bool isPop;

		while (true) {
			/* XXX check a global work counter or task counter? XXX */

			/* first consider the task queue */
      while ((isPop = this->q.pop(_f)) {  // if there is anything in the queue
          std::unique_ptr<std::function<void(int id)>> func(_f); // at return, delete the function even if an exception occurred
          (*_f)(i);
          if (_flag)
              return;  // the thread is wanted to stop, return even if the queue is not empty yet
          else
              isPop = this->q.pop(_f);
      }

			/* retrieve & consume bundle/punc. */
			while ((bundleptr = c->getOneBundle(&trans))) {
				assert(trans);
				try {
					trans->ExecEvaluator(0, /* nodeid, XXX */, c, bundleptr);
				} catch (const std::exception& e) {
	        std::cout << " a standard exception was caught, with message '"
	                  << e.what() << "'\n";
	        abort();
				}
			}

			/* no work to do. wait to be signaled */
			{
				std::unique_lock<mutex> lck(cv_mtx_);
				nwait_ ++;
				cv_.wait(lck, [] {
						/* XXX check a global work counter? do we really have work to do? XXX */
						return true;
					});
				/* woken up and have the lock */
				nwait_ --;
			}
			/* unlocked */
		}
	}

	/* ----------
	 * task queue
	 * ---------- */

public:
  // pops a functional wrapper to the original function
  std::function<void(int)> pop() {
      std::function<void(int id)> * _f = nullptr;
      this->q.pop(_f);
      std::unique_ptr<std::function<void(int id)>> func(_f); // at return, delete the function even if an exception occurred
      std::function<void(int)> f;
      if (_f)
          f = *_f;
      return f;
  }

  template<typename F, typename... Rest>
  auto push(F && f, Rest&&... rest) ->std::future<decltype(f(0, rest...))> {
      auto pck = std::make_shared<std::packaged_task<decltype(f(0, rest...))(int)>>(
          std::bind(std::forward<F>(f), std::placeholders::_1, std::forward<Rest>(rest)...)
          );
      auto _f = new std::function<void(int id)>([pck](int id) {
          (*pck)(id);
      });
      this->q.push(_f);
      std::unique_lock<std::mutex> lock(this->mutex);
      this->cv.notify_one();
      return pck->get_future();
  }

  // run the user's function that excepts argument int - id of the running thread. returned value is templatized
  // operator returns std::future, where the user can get the result and rethrow the catched exceptins
  template<typename F>
  auto push(F && f) ->std::future<decltype(f(0))> {
      auto pck = std::make_shared<std::packaged_task<decltype(f(0))(int)>>(std::forward<F>(f));
      auto _f = new std::function<void(int id)>([pck](int id) {
          (*pck)(id);
      });
      this->q.push(_f);
      std::unique_lock<std::mutex> lock(this->mutex);
      this->cv.notify_one();
      return pck->get_future();
  }

private:
  detail::Queue<std::function<void(int id)> *> q;  /* task queue */
	EvaluationBundleContext *c_;
	vector<pthread_t> threads_;
	static long nwait_;
};


#endif /* EXECUTOR_H_ */
