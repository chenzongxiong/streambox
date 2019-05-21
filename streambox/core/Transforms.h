/*
 * transforms.cpp
 *
 *  Created on: Jun 18, 2016
 *      Author: xzl
 *
 * since all transforms will ultimately include this file, be cautious in
 * including any file here
 */

#ifndef TRANSFORMS_H
#define TRANSFORMS_H

#include <assert.h>
#include <typeinfo>
#include <string>
#include <set>
#include <map>
#include <list>
#include <tuple>
#include <unordered_set>
#include <unordered_map>
#include <utility>

#include "boost/thread/thread.hpp" // uses macro V()
#include "boost/date_time/posix_time/posix_time.hpp"

#include "config.h"
#include "Pipeline.h"
#include "PipelineVisitor.h"
#include "Values.h"

#include "log.h"

#include "tbb/concurrent_unordered_map.h"

using namespace std;

//class PValue;

class EvaluationBundleContext;

// xzl: avoid introducing template args here.
// in cpp, we can't bound the template args. and in Java, that seems to implement polymorphism?
//template <class DerivedTransform>
class PTransform {
public:
	// hym: for dual input transforms, indicate which side the current transform is, L or R?
 	// 0: downstream only has one input container stream, output bundles should be put into container_ directly
	// 1: downstream has 2 input container streams, 1 means output bundles should be put into left_containers_in
	// 2: downstream has 2 input container streams, 2 means output bundles should be put into right_containers_in
	// 3: current transform has 2 input container streams, and 2 output container streams
	//     			 SimpleMapper1: 1
	// UnbundedInMem: 0 ==>                   ==>  Join: 3 ==> Sink: 0
	//     			 SimpleMapper2: 2
	// XXX side_info should be assigned to containers
	int side_info = 0; // default: downstream is a single input transform
	void set_side_info(int i){
		side_info = i;
	}

	int const & get_side_info(){
		return side_info;
	}

	/* hym: for join and its downstream(type 4 trans)
	 * xzl: partial wms for left and right
	 */
	ptime left_wm;
	ptime right_wm;
	// They are moved here from Join.h

	//hym: apply containr to Join transform
	// Join has two input containers streams, and two ouput container streams
	// Upstream will deposit bundles to left_containers_in or right_containers_in
	// Join evaluator will get bundles from those input containers, and
	// deposit bundles to left_containers_out or right_containers_out
	// when Join process a punc to seal a container, it will move some bundles from
	// left_container_out or right_containerP_out to containers_, which is ordered
	// and owned by Join't downstream
	// XXX TODO downstream grab bundles from left_containers_in or right_containers_in
	// directly

	//hym: these lists are for Join transform Only
	list<bundle_container> left_containers_in; //left input stream
	list<bundle_container> left_containers_out; //left outpout stream
	list<bundle_container> right_containers_in; //right input stream
	list<bundle_container> right_containers_out; //right output stream
	// above are only useful to Join

 	///////////////////////////////////
		// hym: for new design
	///////////////////////////////////
	//hym: these lists are for the transform following Join
	// type 4 trans
	list<bundle_container> left_unordered_containers;
	list<bundle_container> right_unordered_containers;
	list<bundle_container> ordered_containers;

	//hym: locks for lists above
	// type 4 trans
	mutex left_mtx_unordered_containers;
	mutex right_mtx_unordered_containers;
	mutex mtx_ordered_containers;
	////////////// end type 4///////////////////////

	//hym: type 5 trans
	//     has 2 types of containers
	// 	-- unordered_containers
	//	-- ordered_containers
	////////////////////////////////////
        list<bundle_container> unordered_containers;
	//list<bundle_container> ordered_containers; //use type 4's
	mutex mtx_unordered_containers;
	//mutex mtx_ordered_containers; //use type 4's
	/////////////end type 5///////////////////////

  // a unified implementation
	virtual PValue* apply(PValue* input) {
    return PCollection::createPrimitiveOutputInternal(input->getPipeline());
	}

	void validate(PValue* input) { }

	std::string getKindString() {
		// TODO: see Beam impl.
		return string(typeid(*this).name());
	}

	std::string getName() {
	    return (!name.empty()) ? name : getKindString();
	}

	// we don't override anything (unlike Java)
	std::string toString() {
		if (name.empty()) {
			return getKindString();
		} else
			return getName() + " [" + getKindString() + "]";
	}

	virtual ~PTransform() { }

	// we want to hide the internal organization of inputs/outputs

	PValue* getFirstInput() {
	  xzl_assert(!inputs.empty());
	  return inputs[0];
	}

	PValue* getSecondInput() {
	    xzl_assert(!inputs.empty());
	    if (inputs.size() < 2)
	      return nullptr;
	    return inputs[1];
	}

	PValue* getFirstOutput() {
	  xzl_assert(!outputs.empty());
	  return outputs[0];
	}

  PValue* getSecondOutput() {
    xzl_assert(!outputs.empty());
    if (outputs.size() < 2)
      return nullptr;
    return outputs[1];
  }

	// visitor will directly push to it..
  vector<PValue *> inputs;
  vector<PValue *> outputs;

  mutex	mtx_output;  /* keep order among output bundles/punc */
  mutex _mutex;  /* internal state */

  // upstream wm that has been observed lastly. need mtx_watermark
  ptime upstream_wm;

  /* tracking inflight bundles and received (but not yet emitted) punc */
  /*  --- protected by mtx_container --- */
  struct inflight_bundle_container {
  	ptime wm = max_date_time; /* wm unassigned */
  	unordered_set<shared_ptr<BundleBase>> bundles;
  };
  list<inflight_bundle_container> inflight_bundle_containers;

  // bundles being worked on by parallel evaluators
  // [ different from the bundles in the input PValue ]
  // XXX move this to SingleInputTransform
  unordered_set<shared_ptr<BundleBase>> inflight_bundles;

  unordered_set<shared_ptr<BundleBase>> dual_inflight_bundles[2];

  boost::shared_mutex mtx_container; /* protecting all containers */

  //hym: for Join
  // for left input and output
  mutex left_mtx_container_in;
  mutex left_mtx_container_out;
  mutex right_mtx_container_in;
  mutex right_mtx_container_out;

  /* upstream will directly deposit bundles */
//  tbb::concurrent_vector<bundle_container> containers_;
//  atomic<long> containers_low_index_;
//  tbb::concurrent_queue<bundle_container> containers_; // can't access front/end
//  vector<bundle_container> containers_;

  /* oldest at back(). so push_front() appends a new container.
   *
   * a bundle may store a pointer to its enclosing container.
   * so list helps here -- the element pointer should always be
   * valid.
   * */
  list<bundle_container> containers_;

  /* a relaxed counter for bundles. for pressure feedback.
   * NB: we count the bundles in container, i.e. not including inflight
   * */
  atomic<unsigned long> bundle_counter_;

  /* helper */
  void dump_containers_legend() {
  	cout << "legend for puncref:";
  }

  /* A container can be empty, e.g. after retrieving the last bundle */
  /* unsafe. need lock */
	void dump_containers(const char * msg = "") {
#if (CONFIG_KAGE_GLOBAL_DEBUG_LEVEL <= 50)  /* VV(...) */

		printf("dump_containers ('%s'): %s total %lu",
				msg, this->name.c_str(),
				this->containers_.size());
		cout << endl;

		for (auto && cont : this->containers_) {
			/* snapshot */
			long refcnt = cont.refcnt.load();
			long punc_refcnt = cont.punc_refcnt.load();
//			auto s = cont.bundles.size();
			auto s = cont.peekBundleCount();

			printf("%lx: ", (unsigned long)(&cont));
			if (cont.downstream) {
				printf("down:%lx ", (unsigned long)(cont.downstream.load()));
			} else
				printf("down:NA ");

//			if (cont.punc)
//				cout << to_simplest_string1(cont.punc->min_ts).str() << "  ";
//			else
//				cout << "X";
//				cout << "(nullptr)";
//			cout << " puncrefcnt: " << punc_refcnt;
			cout << puncref_key(punc_refcnt) << " ";
//			cout << endl;

//			cout << "bundles: total " << refcnt;
//			cout << " unretrieved: " << s;
//			cout << " outstanding: " << (refcnt - s);
////			cout << endl;

			cout << "ref:" << refcnt;
			cout << "(in:" << s << ") ";
			cout << "|  ";
		}
		cout << endl;
#endif
	}

	void dump_containers_left_in(const char * msg = "") {
//#if (CONFIG_KAGE_GLOBAL_DEBUG_LEVEL <= 40)  /* VV(...) */

		printf("dump_containers: %s ('%s') total %lu",
				this->name.c_str(), msg,
				this->left_containers_in.size());
		cout << endl;

		for (auto && cont : this->left_containers_in) {
			/* snapshot */
			long refcnt = cont.refcnt.load();
			long punc_refcnt = cont.punc_refcnt.load();
			auto s = cont.bundles.size();

			printf("%lx: ", (unsigned long)(&cont));
			if (cont.downstream) {
				printf("down:%lx ", (unsigned long)(cont.downstream.load()));
			} else
				printf("down:NA ");

//			if (cont.punc)
//				cout << to_simplest_string1(cont.punc->min_ts).str() << "  ";
//			else
//				cout << "X";
//				cout << "(nullptr)";
//			cout << " puncrefcnt: " << punc_refcnt;
			cout << punc_refcnt << " ";
//			cout << endl;

//			cout << "bundles: total " << refcnt;
//			cout << " unretrieved: " << s;
//			cout << " outstanding: " << (refcnt - s);
////			cout << endl;

			cout << "ref:" << refcnt;
			cout << "(in:" << s << ") ";
			cout << "|  ";
		}
		cout << endl;
//#endif
	}

	void dump_containers_right_in(const char * msg = "") {
//#if (CONFIG_KAGE_GLOBAL_DEBUG_LEVEL <= 40)  /* VV(...) */

		printf("dump_containers: %s ('%s') total %lu",
				this->name.c_str(), msg,
				this->right_containers_in.size());
		cout << endl;

		for (auto && cont : this->right_containers_in) {
			/* snapshot */
			long refcnt = cont.refcnt.load();
			long punc_refcnt = cont.punc_refcnt.load();
			auto s = cont.bundles.size();

			printf("%lx: ", (unsigned long)(&cont));
			if (cont.downstream) {
				printf("down:%lx ", (unsigned long)(cont.downstream.load()));
			} else
				printf("down:NA ");

//			if (cont.punc)
//				cout << to_simplest_string1(cont.punc->min_ts).str() << "  ";
//			else
//				cout << "X";
//				cout << "(nullptr)";
//			cout << " puncrefcnt: " << punc_refcnt;
			cout << punc_refcnt << " ";
//			cout << endl;

//			cout << "bundles: total " << refcnt;
//			cout << " unretrieved: " << s;
//			cout << " outstanding: " << (refcnt - s);
////			cout << endl;

			cout << "ref:" << refcnt;
			cout << "(in:" << s << ") ";
			cout << "|  ";
		}
		cout << endl;
//#endif
	}

//private
public:  /* xzl: XXX workaround, since these are needed by Join.h */
	/*
	 * support rlock --> ulock and ulock --> wlock
	 *
	 * - if called with rlock held, release; acquire ulock. need restart.
	 * - if called with ulock held, upgrade to wlock in place. no restart needed.
	 *
	 * return: true if wlock (exclusive) is acquired. can proceed to modify
	 * data structures
	 * false if rlock is released and ulock is acquired. need to restart.
	 *
	 * since boost::upgrade_to_unique_lock must be init'd once its declared,
	 * we dynamically create it to defer its init. to avoid manual mgmt of the
	 * dyn object, the caller can simply point a unique_ptr to it.
	 *
	 */
	bool upgrade_locks(boost::shared_lock<boost::shared_mutex>* prlock,
			boost::upgrade_lock<boost::shared_mutex> *pulock,
			unique_ptr<boost::upgrade_to_unique_lock<boost::shared_mutex>>* ppwlock)
	{
		/* at least one lock is owned */
		xzl_assert(prlock->owns_lock()
				|| pulock->owns_lock()
				|| (ppwlock && *ppwlock && (*ppwlock)->owns_lock()));

		auto & rlock = *prlock;
		auto & ulock = *pulock;
		auto & pwlock = *ppwlock;

 		if (rlock.owns_lock()) {  /* upgrade rlock -> ulock */
 			xzl_assert(!ulock.owns_lock());
 			rlock.unlock();
 			ulock.lock();
 			return false;
 		}

 		if (ulock.owns_lock()) { /* upgrade ulock -> wlock in place. */
 			xzl_assert(!rlock.owns_lock());
 			xzl_assert(!pwlock); /* must points to nullptr */
// 			(*ppwlock) = std::move(
// 			pwlock = std::move(
//// 					make_unique<boost::upgrade_to_unique_lock<boost::shared_mutex>>(ulock)
// 					unique_ptr<boost::upgrade_to_unique_lock<boost::shared_mutex>>
// 							(new boost::upgrade_to_unique_lock<boost::shared_mutex>(ulock))
// 					);

 			/* also okay, since assign temp object == move */
 			pwlock =
 					unique_ptr<boost::upgrade_to_unique_lock<boost::shared_mutex>>
 							(new boost::upgrade_to_unique_lock<boost::shared_mutex>(ulock));
 			xzl_assert(pwlock);
 			/* no need start over since we atomically upgrade */
 		}

 		xzl_assert(!rlock.owns_lock() && !ulock.owns_lock()
 				&& pwlock->owns_lock());

 		return true;
	}

public:
	/* Maintain the (relaxed) bundle counter.
	 * return: the value before increment
	 * XXX combine the following two? */
	inline long IncBundleCounter(long delta = 1) {
		return bundle_counter_.fetch_add(delta, std::memory_order_relaxed);
	}

	inline long DecBundleCounter(long delta = 1) {
		return bundle_counter_.fetch_sub(delta, std::memory_order_relaxed);
	}

	/* return the # of *unretrieved* bundles (i.e. not including inflight ones)
	 * in this transform.
	 * Does not count punc.
	 * NB: this is frequently invoked by source to measure the pipeline pressure.
	 * So it has to be cheap.
	 * */
	inline long getNumBundles() {
		return bundle_counter_.load(std::memory_order_relaxed);
	}

#if 0 /* too expensive */
	long getNumBundles() {
		/* lock all containers in case some one suddenly go away */
  	unique_lock<mutex> conlock(mtx_container);
  	long ret = 0;
  	for (auto && cont : this->containers_) {
  		ret +=  cont.refcnt.load();
  	}
  	return ret;
	}
#endif

	/* For statistics: get the bundle count for *each* container and return as
	 * a vector. oldest at the back.
	 *
	 * Can be expensive */
	void getNumBundles(vector<long>* bundles) {
	#if 0
		/* lock all containers in case some one suddenly go away */
		unique_lock<mutex> conlock(mtx_container);
		for (auto && cont : this->containers_) {
			bundles->push_back(cont.refcnt.load());
		}

	#endif
		/* lock all containers in case some one suddenly go away */
//		unique_lock<mutex> conlock(mtx_container);
		boost::shared_lock<boost::shared_mutex> rlock(mtx_container);

		for (auto && cont : this->containers_) {
//			bundles->push_back(cont.refcnt.load());  /* in container + inflight */
//			bundles->push_back(cont.bundles.size());   /* in container only */

			/* w/o lock, the following two may be inconsistent */
			auto c = cont.peekBundleCount();
			if (c == 0) { /* encode perfcnt as neg */
				c = cont.punc_refcnt.load() - PUNCREF_MAX;
//				xzl_assert(c < 0);
			}

			bundles->push_back(c);   /* in container only */
		}
	}

	struct cont_info {
		long bundle;
		ptime punc;
		int side_info;
	};

	void getNumBundlesTs(vector<cont_info>* bundles) {
		/*hym:
		 *     Have already gotten lock before call this function
		 *     so should not get lock again here, otherwise deadlock will happen
		 */
		/* lock all containers in case some one suddenly go away */
//		unique_lock<mutex> conlock(mtx_container);
		//boost::shared_lock<boost::shared_mutex> rlock(mtx_container);
		//old trans, SimpleMapper 1, SimpleMapper2
		if(this->get_side_info() == 0 || this->get_side_info() == 1 || this->get_side_info() == 2){
			for (auto && cont : this->containers_) {
				//			bundles->push_back(cont.refcnt.load());  /* in container + inflight */
				//			bundles->push_back(cont.bundles.size());   /* in container only */
				/* w/o lock, the following two may be inconsistent */
				cont_info info;

				info.bundle  = cont.peekBundleCount();
				if (info.bundle == 0) { /* encode perfcnt as neg */
					info.bundle = cont.punc_refcnt.load() - PUNCREF_MAX;
					xzl_assert(info.bundle < 0);
				}

                // zxchen: punc is leaking
				auto p = cont.getPuncSafe();
				if (p) {
					info.punc = p->min_ts;
					info.side_info = p->get_side_info();
				}

				bundles->push_back(info);   /* in container only */
			}
		//join
		}else if(this->get_side_info() == 3){
			//left_containers_in
			for (auto && cont : this->left_containers_in) {
				//			bundles->push_back(cont.refcnt.load());  /* in container + inflight */
				//			bundles->push_back(cont.bundles.size());   /* in container only */
				/* w/o lock, the following two may be inconsistent */
				cont_info info;

				info.bundle  = cont.peekBundleCount();
				if (info.bundle == 0) { /* encode perfcnt as neg */
					info.bundle = cont.punc_refcnt.load() - PUNCREF_MAX;
					xzl_assert(info.bundle < 0);
				}

				auto p = cont.getPuncSafe();
				if (p) {
                    // zxchen: invoke memory error, since `min_ts` is freed.
                    // std::cout << "p->min_ts: " << p->min_ts << std::endl;
					info.punc = p->min_ts;
					info.side_info = p->get_side_info();
				}

				info.side_info = cont.get_side_info();  // debugging
				xzl_assert(info.side_info == 1);

				bundles->push_back(info);   /* in container only */
			}
			//right_containers_in
			for (auto && cont : this->right_containers_in) {
				//			bundles->push_back(cont.refcnt.load());  /* in container + inflight */
				//			bundles->push_back(cont.bundles.size());   /* in container only */
				/* w/o lock, the following two may be inconsistent */
				cont_info info;

				info.bundle  = cont.peekBundleCount();
				if (info.bundle == 0) { /* encode perfcnt as neg */
					info.bundle = cont.punc_refcnt.load() - PUNCREF_MAX;
					xzl_assert(info.bundle < 0);
				}

				auto p = cont.getPuncSafe();
				if (p) {
					info.punc = p->min_ts;
					info.side_info = p->get_side_info();
				}
				info.side_info = cont.get_side_info();  // debugging
				xzl_assert(info.side_info == 2);

				bundles->push_back(info);   /* in container only */
			}
		}else if(this->get_side_info() == 4){
			//left_unordered_containers
			for (auto && cont : this->left_unordered_containers) {
				//			bundles->push_back(cont.refcnt.load());  /* in container + inflight */
				//			bundles->push_back(cont.bundles.size());   /* in container only */
				/* w/o lock, the following two may be inconsistent */
				cont_info info;

				info.bundle  = cont.peekBundleCount();
				if (info.bundle == 0) { /* encode perfcnt as neg */
					info.bundle = cont.punc_refcnt.load() - PUNCREF_MAX;
					xzl_assert(info.bundle < 0);
				}

				auto p = cont.getPuncSafe();
				if (p) {
					info.punc = p->min_ts;
					info.side_info = p->get_side_info();
				}

				bundles->push_back(info);   /* in container only */
			}
			//right_unordered_containers
			for (auto && cont : this->right_unordered_containers) {
				//			bundles->push_back(cont.refcnt.load());  /* in container + inflight */
				//			bundles->push_back(cont.bundles.size());   /* in container only */
				/* w/o lock, the following two may be inconsistent */
				cont_info info;

				info.bundle  = cont.peekBundleCount();
				if (info.bundle == 0) { /* encode perfcnt as neg */
					info.bundle = cont.punc_refcnt.load() - PUNCREF_MAX;
					xzl_assert(info.bundle < 0);
				}

				auto p = cont.getPuncSafe();
				if (p) {
					info.punc = p->min_ts;
					info.side_info = p->get_side_info();
				}

				bundles->push_back(info);   /* in container only */
			}
			//ordered_containers
			for (auto && cont : this->ordered_containers) {
				//			bundles->push_back(cont.refcnt.load());  /* in container + inflight */
				//			bundles->push_back(cont.bundles.size());   /* in container only */
				/* w/o lock, the following two may be inconsistent */
				cont_info info;

				info.bundle  = cont.peekBundleCount();
				if (info.bundle == 0) { /* encode perfcnt as neg */
					info.bundle = cont.punc_refcnt.load() - PUNCREF_MAX;
					xzl_assert(info.bundle < 0);
				}

				auto p = cont.getPuncSafe();
				if (p) {
					info.punc = p->min_ts;
					info.side_info = p->get_side_info();
				}

				bundles->push_back(info);   /* in container only */
			}
		}else if(this->get_side_info() == 5){
			//unordered_containers
			for (auto && cont : this->unordered_containers) {
				//			bundles->push_back(cont.refcnt.load());  /* in container + inflight */
				//			bundles->push_back(cont.bundles.size());   /* in container only */
				/* w/o lock, the following two may be inconsistent */
				cont_info info;

				info.bundle  = cont.peekBundleCount();
				if (info.bundle == 0) { /* encode perfcnt as neg */
					info.bundle = cont.punc_refcnt.load() - PUNCREF_MAX;
					xzl_assert(info.bundle < 0);
				}

				auto p = cont.getPuncSafe();
				if (p) {
					info.punc = p->min_ts;
					info.side_info = p->get_side_info();
				}

				bundles->push_back(info);   /* in container only */
			}
			//ordered_containers
			for (auto && cont : this->ordered_containers) {
				//			bundles->push_back(cont.refcnt.load());  /* in container + inflight */
				//			bundles->push_back(cont.bundles.size());   /* in container only */
				/* w/o lock, the following two may be inconsistent */
				cont_info info;

				info.bundle  = cont.peekBundleCount();
				if (info.bundle == 0) { /* encode perfcnt as neg */
					info.bundle = cont.punc_refcnt.load() - PUNCREF_MAX;
					xzl_assert(info.bundle < 0);
				}

				auto p = cont.getPuncSafe();
				if (p) {
					info.punc = p->min_ts;
					info.side_info = p->get_side_info();
				}

				bundles->push_back(info);   /* in container only */
			}
		}else{
			assert(false && "wrong side info in getNumBundlesTs");
		}
	}

	string getContainerInfo() {
		ostringstream oss;
		vector<cont_info> cnts;
		this->getNumBundlesTs(&cnts);
		//  	oss << "deposit to: " << this->name << "\n (bundles: ";
		oss << this->name << " (bundles: ";
		for (auto & cnt : cnts) {
			if (cnt.bundle >= 0)
				oss << cnt.bundle << " ";
			else {
				auto r = cnt.bundle + PUNCREF_MAX;
				oss << puncref_key(r) << " ";
				if (r >= 0) { /* punc ever assigned. show punc ts */
					oss << to_simple_string(cnt.punc) << " ";
				}
			}
		}
		oss << ")";
		return oss.str();
	}

#if 0 /* obsoleted */
	/* Decide where to send the output bundles: 1) staging in the container that
	 * encloses the input bundle or 2) directly to downstream?
	 *
	 * 1) iff the container is *not* the oldest container in the transform, i.e.
	 * there is a punc ahead of the output bundles to be assigned.
	 *
	 * 2) otherwise
	 *
	 * must be called when @input_bundle is still valid (i.e. not consumed;
	 * refcnt not decremented)
	 *
	 * return: true if staging okay. false if direct output is needed.
	 */
	bool tryStageBundles(shared_ptr<BundleBase> input_bundle,
			vector<shared_ptr<BundleBase>> const & output_bundles)
	{

		xzl_assert(input_bundle->container); /* enclosing container */

		unique_lock<mutex> conlock(mtx_container);

		xzl_assert(!containers_.empty() && "bug? we at least have one bundle here");

		/* since @containers_ is a list, we can safely take addr of its element */

#ifndef NDEBUG /* debug */
		bool flag = false;
		for (auto & c : containers_) {
			if (&c == input_bundle->container) {
				flag = true;
				break;
			}
		}
		xzl_assert(flag && "bug: enclosing container does not exist?");
#endif

		if (&(containers_.back()) == input_bundle->container) {
			/* the enclosing container is already the oldest. no punc in front of
			 * these @output_bundles. should direct deposit them.
			 *
			 * there may be output bundles staged in the container, which will be deposited
			 * to downstream right before we process the punc of this container.
			 */
			return false;
		}

		/* The enclosing container is *not* the oldest. there're punc(s) ahead.
		 * (or, there's a dead container ahead).
		 *
		 * Stage them.
		 */
		input_bundle->container->staged_bundles.insert(
				input_bundle->container->staged_bundles.end(),
				output_bundles.begin(), output_bundles.end());

//		EE("staging bundles okay");

		return true;
	}
#endif

  /* deposit a data bundle. a later punc can only be deposited
   * after this func returns. */
  void depositOneBundle(shared_ptr<BundleBase> bundle, int node = -1) {

  	/* lock all containers until reliably getting the target container.
  	 * (in case concurrent additions of containers) */
//  	unique_lock<mutex> conlock(mtx_container);

  	boost::shared_lock<boost::shared_mutex> rlock(mtx_container);
  	boost::upgrade_lock<boost::shared_mutex> ulock(mtx_container, boost::defer_lock);
  	unique_ptr<boost::upgrade_to_unique_lock<boost::shared_mutex>> pwlock
  							= nullptr; /* will std::move later */

start:
  	if (containers_.empty()) {
   		if (!upgrade_locks(&rlock, &ulock, &pwlock))
   			goto start;
  		containers_.emplace_front(); /* can't copy */
  	}

  	auto current_con = containers_.begin();
	current_con->set_side_info(this->get_side_info());

  	if (current_con->getPuncSafe()) {

  		/* the latest container already has a punc. open a new container */

		if (!upgrade_locks(&rlock, &ulock, &pwlock))
			goto start;

  		xzl_assert(current_con->punc_refcnt != PUNCREF_UNDECIDED
  				&& "invalid punc refcnt");

		/* The latest container is already dead: clean up.
		 * NB: even if we see punc_refcnt == 1 here, it may be dead
		 * (decremented to 0) beneath us any time. It will be cleaned by future
		 * calls.
		 */
		if (current_con->punc_refcnt == PUNCREF_CONSUMED) {
			/* extensive sanity check... */

			//assert(current_con->punc && "dead container but no valid punc?");
	    		xzl_assert(current_con->getPuncSafe() && "dead container but no valid punc?");
	    		xzl_assert(current_con->peekBundleCount() == 0);
	    		xzl_assert(current_con->refcnt == 0);

			current_con = containers_.erase(current_con);

			if (containers_.empty()) {
				/* after cleanup, no container left. start a new one. */
				containers_.emplace_front();
				current_con = containers_.begin();
				goto open_container;
			} else {
				/* current_con now points to the 2nd most recent container.
				 * sanity check it */
				xzl_assert(current_con->punc_refcnt != PUNCREF_CONSUMED
						&& "bug: can't have two dead containers");
			}
		}

		/* careful when reading from an atomic/volatile value ... */
		auto prefcnt = current_con->punc_refcnt.load();
		xzl_assert( prefcnt == PUNCREF_CONSUMED   /* latest container dies just now */
				|| prefcnt == PUNCREF_RETRIEVED		/* outstanding */
				|| prefcnt == PUNCREF_ASSIGNED		/* unretrieved yet */
		      );

		containers_.emplace_front();
		current_con --;
	}

#if 0
	/* now @current_con points to a valid, non-dead container */

	if (current_con->punc) {
		/* latest container already has a punc.
		 * It is sealed. create a new container. */
		xzl_assert(current_con->punc_refcnt != -1 && "must have a valid refcnt");

		/* current_con->punc_refcnt == 0 may happen: it was 1 last time we checked
		 * it; and it wad decremented underneath us.
		 * It seems okay that we leave it here and let future
		 */
		xzl_assert(current_con->punc_refcnt != 0 && "XXX clean up dead container?");

		containers_.emplace_front();
		current_con --;
	}
#endif

	/* Now @current_con points to an open container. Insert the bundle.
	 *
	 * @current_con won't go away since the punc won't be emitted
	 * concurrently.
	 */
open_container:
		/* should have any type of lock */
		xzl_assert(rlock.owns_lock()
						  || ulock.owns_lock() || (pwlock && pwlock->owns_lock()));

	current_con->putBundleSafe(bundle);
	this->IncBundleCounter();

	//#ifdef MEASURE_LATENCY
	//  	bundle->mark("deposited: " + this->name);
	//#endif

  }

//private:
public:  /* xzl: workaround XXX since it is needed by join */
  bundle_container* localBundleToContainer(
  		shared_ptr<BundleBase> const & input_bundle)
  {
  	bundle_container *upcon = input_bundle->container;
//XXX hym
#if 0 //not all contaier lists called_containers_, so here will be a bug if the container list's name is not containers_
/* No need to lock. Just peek. */
#ifdef DEBUG
		bool flag = false;
		for (auto & c : containers_) {
			if (&c == upcon) {
				flag = true;
				break;
			}
		}
		xzl_assert(flag && "bug: enclosing container does not exist?");
#endif
#endif
		return upcon;
  }

public:
  /* deposit a bundle not to *this* transform, but to the given *downstream*
   * trans.
   * Does not lock all containers_.
   * XXX allow specifying @node for each bundle?
   * */
  void  depositBundlesDownstream(PTransform *downt,
  		shared_ptr<BundleBase> input_bundle,
			vector<shared_ptr<BundleBase>> const & output_bundles, int node = -1)
  {
  	bool try_open = false; /* have we tried open downstream containers? */
	xzl_assert(downt);
  	bundle_container *upcon = localBundleToContainer(input_bundle);
  	xzl_assert(upcon); /* enclosing container */

#if 0
	/* reader lock for *this* trans */
	boost::shared_lock<boost::shared_mutex> conlock(mtx_container);
	/* writer lock for downstream trans */
	boost::unique_lock<boost::shared_mutex> down_conlock(downt->mtx_container);
#endif

deposit:
     {
			if (upcon->downstream) {
				/* fast path. NB @downstream is atomic w/ seq consistency so it's
				 * safe. */
				downt->depositBundlesToContainer(upcon->downstream, output_bundles, node);
				downt->IncBundleCounter(output_bundles.size());
				return;
			}

  	xzl_assert(!try_open && "bug: already tried to open");

  	/* it's possible somebody slipped in & opened the downstream container
  	 * already. So ret == 0  from openDownstreamContainers() then.
       NB: openDownstreamContainer must take lock internally.
  	 */
	if(this->get_side_info() == 0 || this->get_side_info() == 3){ // defaut transforms and join
		openDownstreamContainers(upcon, downt);
	}else if(this->get_side_info() == 1){ //Simplemapper 1, left
		//std::cout << __FILE__ << ": " <<  __LINE__ << std::endl;
		openDownstreamContainers_left(upcon, downt); //downt should be join, open a container in join's left_containers_in
		upcon->set_side_info(1);
		upcon->downstream.load()->set_side_info(1);
	}else if(this->get_side_info() == 2){ //Simplemapper 2, right
		//std::cout << __FILE__ << ": " <<  __LINE__ << std::endl;
		openDownstreamContainers_right(upcon, downt); //downt should be join, open a container in join's right_containers_in
		upcon->set_side_info(2);
		upcon->downstream.load()->set_side_info(2);
	}else{
		std::cout << "side info is " << this->get_side_info() << std::endl;
		assert(false && "Bug: Unknown know side info");
	}
     }

	xzl_assert(upcon->downstream);
  	try_open = true;
  	goto deposit;
  }


  /* see commment above */
  void depositOnePuncDownstream(PTransform* downt,
  		shared_ptr<Punc> const & input_punc, shared_ptr<Punc> punc, int node = -1)
  {
  	bool try_open = false;

  	xzl_assert(downt);
  	bundle_container *upcon = localBundleToContainer(input_punc);
  	xzl_assert(upcon); /* enclosing container */

#if 0
	/* reader lock for *this* trans */
	boost::shared_lock<boost::shared_mutex> conlock(mtx_container);
	/* writer lock for downstream trans */
	boost::unique_lock<boost::shared_mutex> down_conlock(downt->mtx_container);
#endif

deposit:
     {
			if (upcon->downstream) { /* fast path. NB @downstream is atomic, implying fence */
				downt->depositOnePuncToContainer(upcon->downstream, punc, node);
				upcon->downstream.load()->has_punc = 1;
				return;
			}

		xzl_assert(!try_open && "bug: already tried to open");

		//openDownstreamContainers(upcon, downt);
		//hym:XXX
		if(this->get_side_info() == 0 || this->get_side_info() == 3){ //defaut transforms and join
			openDownstreamContainers(upcon, downt);
		}else if(this->get_side_info() == 1){ //Simplemapper 1, left
			openDownstreamContainers_left(upcon, downt); //downt should be join, open a container in join's left_containers_in
			std::cout << __FILE__ << ": " <<  __LINE__ << " open container left because of a punc" << std::endl;
		}else if(this->get_side_info() == 2){ //Simplemapper 2, right
			openDownstreamContainers_right(upcon, downt); //downt should be join, open a container in join's right_containers_in
			std::cout << __FILE__ << ": " <<  __LINE__ << " open container right because of a punc" << std::endl;
		}else{
			xzl_bug("Bug: Unknown know side info");
		}
	}
		xzl_assert(upcon->downstream);
		try_open = true;
		goto deposit;
  }


  /* *This* transform consumes a punc (and therefore closes a container) but
   * does not emit flushed bundles or new punc.
   *
   * Since the up container will be gone soon, we need to mark the downstream
   * container, indicating that i) it will never see a punc and ii) its (effective)
   * punc should come from a later downstream container on the same line.
   *
   */
  void cancelPuncDownstream(PTransform* downt, shared_ptr<Punc> input_punc)
  {
  	xzl_assert(downt);
  	bundle_container *upcon = localBundleToContainer(input_punc);
  	xzl_assert(upcon); /* enclosing container */

deposit:
		if (upcon->downstream) { /* fast path */
			long expected = PUNCREF_UNDECIDED;
			if (upcon->downstream.load()->punc_refcnt.compare_exchange_strong(
					expected, PUNCREF_CANCELED)) {
				return;
			} else {
				bug("invalid punc refcnt?");
			}
		}

		/* if down container not existing, still wants to open one to make
		 * code logic clearer
		 */
		//int ret = openDownstreamContainers(upcon, downt);
		int ret;
		if(this->get_side_info() == 0 || this->get_side_info() == 3){ //defaut transforms and join
			ret = openDownstreamContainers(upcon, downt);
		}else if(this->get_side_info() == 1){ //Simplemapper 1, left
			ret = openDownstreamContainers_left(upcon, downt); //downt should be join, open a container in join's left_containers_in
		}else if(this->get_side_info() == 2){ //Simplemapper 2, right
			ret =openDownstreamContainers_right(upcon, downt); //downt should be join, open a container in join's right_containers_in
		}else{
			xzl_bug("Bug: Unknown know side info");
		}
		xzl_assert(ret);
		xzl_assert(upcon->downstream);
		goto deposit;
  }

private:
  void checkContainerLinks(PTransform *downt) {
#ifdef DEBUG
  	/* sanity check: existing containers in up/downstream should be
  	 * properly linked
  	 */
//  	xzl_assert(downt);
//
//		unique_lock<mutex> up_conlock(this->mtx_container);
//		unique_lock<mutex> down_conlock(downt->mtx_container);
//
//		auto it1 = this->containers_.rbegin();
//		auto it2 = downt->containers_.rbegin();
//
//		for (; it1 != this->containers_.rend() && it2 = downt->containers_.rend();
//					it1++, it2++)
//		{
//			if (it1->downstream != )
//		}

#endif
  }

private:
  /*
   * will lock both transform's containers.
   * XXX avoid?
   *
   * return: # of downstream containers opened */
  int openDownstreamContainers(bundle_container * const upcon,
  		PTransform *downt)
  {
  	xzl_assert(downt && upcon);

//  	/* fast path w/o locking */
//  	if (upcon->downstream)
//  		return 0;

  	/* slow path
  	 *
  	 * in *this* trans, walk from the oldest container until @upcon,
  	 * open a downstream container if there isn't one.
  	 *
  	 * (X = no downstream container;  V = has downstream container)
  	 * possible:
  	 *
  	 * X X V V (oldest)
  	 *
  	 * impossible:
  	 * X V X V (oldest)
  	 *
  	 */

		/* reader lock for *this* trans */
		boost::shared_lock<boost::shared_mutex> conlock(mtx_container);
		/* writer lock for downstream trans */
		boost::unique_lock<boost::shared_mutex> down_conlock(downt->mtx_container);

  	if (upcon->downstream) /* check again */
  		return 0;

  	int cnt = 0;
  	bool miss = false;

  	for (auto it = this->containers_.rbegin();
  			it != this->containers_.rend(); it ++) {

  		/* also lock all downstream containers. good idea? */
//  		unique_lock<mutex> down_conlock(downt->mtx_container);

  		if (!it->downstream) {
  			miss = true;
  			/* alloc & link downstream container. NB: it->downstream is an atomic
  			 * pointer, so it implies a fence here. */
  			downt->containers_.emplace_front();
  			it->downstream = &(downt->containers_.front());

			//hym: XXX
			//assign the side_info to this new container
			//it->downstream.load()->side_info = this->get_side_info(); //0->containrs_, 1->left_containers_in, 2->right_containers_in
			it->downstream.load()->set_side_info(this->get_side_info());
			assert(it->downstream.load()->get_side_info() != 1 && it->downstream.load()->get_side_info() != 2); //containers_

  			cnt ++;
  		} else { /* downstream link exists */
  			xzl_assert (!miss && "bug? an older container misses downstream link.");
#ifdef DEBUG
  			bool found = false;
  			for (auto & con : downt->containers_) {
  				if (&con == it->downstream) {
  					found = true;
  					break;
  				}
  			}
  			/* NB: once upstream assigns wm to the downstream (i.e. upstream
  			 * is RETRIEVED or CONSUMED), the downstream may process the punc
  			 * and be cleaned up prior to the upstream container.
  			 * In this case, a dangling downstream pointer is possible.
  			 *
  			 * XXX other similar paths should do same thing. XXX
  			 */
  			auto puncref = it->punc_refcnt.load();
  			if (!found
  					&& puncref != PUNCREF_CONSUMED
  					&& puncref != PUNCREF_RETRIEVED) {
  				EE("bug: %s: downstream (%s) container does not exist",
  						this->name.c_str(), downt->name.c_str());
						cout << " -------------------------------\n";
						dump_containers("bug");
						downt->dump_containers("bug");
						cout << " -------------------------------\n\n";
						abort();
  				xzl_bug("downstream container does not exist");
  			}
#endif
  		}
  	}

#if 0 /* debug */
  	cout << " -------------------------------\n";
  	dump_containers("link containers");
  	downt->dump_containers("link containers");
  	cout << " -------------------------------\n\n";
#endif

  	return cnt;
  }




//private:
public:
  /* deposit a data bundle to a @container of *this* transform.
   * a later punc can only be deposited after this func returns.
   *
   * no need to lock containers_.
   * IncBundleCounter() done by caller.
   */
  void depositBundlesToContainer(bundle_container * const container,
  		vector<shared_ptr<BundleBase>> const & output_bundles, int node = -1)
  {
  	/* sanity check */
  	xzl_assert(container);

//XXX hym
#if 0 //not all container lists are named this->containers_, e.g. right_containers_in. so this check may fail
#ifdef DEBUG
  	/* the given @container must exist in this trans */
  	{
//  		unique_lock<mutex> conlock(mtx_container);
			boost::shared_lock<boost::shared_mutex> rlock(mtx_container);

			for (auto & c : this->containers_) {
				if (&c == container)
					goto okay;
			}
  	}
  	xzl_assert(false && "bug? @container does not exist in this trans");
okay:
#endif
#endif

		xzl_assert(container->punc_refcnt == PUNCREF_UNDECIDED
				&& "punc already assigned?");
		xzl_assert(container->verifyPuncSafe() && "punc already assigned?");

		for (auto & b : output_bundles) {
			xzl_assert(b);
			container->putBundleSafe(b); /* see func comment above */
		}
		VV("%lu bundles to container %lx", output_bundles.size(),
				(unsigned long)container);
  }

  /* See comment above.
   * no need to lock containers_. */
  void depositOnePuncToContainer(bundle_container * const container,
  		shared_ptr<Punc> const & punc, int node = -1)
  {
	/* sanity check */
  	xzl_assert(container);
  	xzl_assert(punc);

#ifdef MEASURE_LATENCY
  	/*
	ostringstream oss;
  	vector<long> cnts;
  	this->getNumBundles(&cnts);
  	oss << "deposit to: " << this->name << " (bundles: ";
  	for (auto & cnt : cnts) {
  		oss << cnt << " ";
  	}
  	oss << ")";
  	punc->mark(oss.str());
	*/
	ostringstream oss;
  	vector<cont_info> cnts;
  	this->getNumBundlesTs(&cnts);
	oss << "deposit to: " << this->name << "\n (bundles: ";
  	for (auto & cnt : cnts) {
  		if (cnt.bundle >= 0)
  			oss << cnt.bundle << " ";
  		else {
  			auto r = cnt.bundle + PUNCREF_MAX;
  			oss << puncref_key(r) << " ";
  			if (r >= 0) { /* punc ever assigned. show punc ts */
  				oss << to_simple_string(cnt.punc) << "\t";
  			}
//  			oss << endl;
  		}
  	}

	oss << ")";
  	punc->mark(oss.str());
#endif
//XXX hym
#if 0 // see comment above
#ifdef DEBUG
  	{
  		//  		unique_lock<mutex> conlock(mtx_container);
			boost::shared_lock<boost::shared_mutex> rlock(mtx_container);

			/* make sure @container exists in this trans */
			for (auto & c : this->containers_) {
				if (&c == container)
					goto okay;
			}
  	}
  	xzl_assert(false && "bug? @container does not exist in this trans");
okay:
#endif
#endif

		xzl_assert(container->punc_refcnt == PUNCREF_UNDECIDED
				&& "punc already assigned?");
		xzl_assert(container->verifyPuncSafe() && "punc already assigned?");
		container->has_punc = 1;
		container->setPuncSafe(punc);
  }

public:
  /* used by source only. other transforms should use the container version.
   * does not deal with cancelled punc.
   *
   * must be called after all prior bundles are deposited.
   * (exception: flushed bundles can increase the refcnt first
   * and emit the actual bundles later)
   *
   */
  void depositOnePunc(shared_ptr<Punc> punc, int node = -1) {

#ifdef MEASURE_LATENCY

/*
  	ostringstream oss;
  	vector<long> cnts;
  	this->getNumBundles(&cnts);
  	oss << "deposit to: " << this->name << " (bundles: ";
  	for (auto & cnt : cnts) {
  		oss << cnt << " ";
  	}
  	oss << ")";
  	punc->mark(oss.str());
*/
  	ostringstream oss;
  	vector<cont_info> cnts;
  	this->getNumBundlesTs(&cnts);
  	oss << "deposit to: " << this->name << "\n (bundles: ";
  	for (auto & cnt : cnts) {
  		if (cnt.bundle >= 0)
  			oss << cnt.bundle << " ";
  		else {
  			auto r = cnt.bundle + PUNCREF_MAX;
  			oss << puncref_key(r) << " ";
  		}
  		//  			if (r >= 0) { /* punc ever assigned. show punc ts */
			if (cnt.punc != max_date_time) {
					oss << to_simple_string(cnt.punc) << "\t";
			}
  	}
  	oss << ")";
  	punc->mark(oss.str());
//  	punc->mark("deposit to: " + this->name);
#endif
    /* fast path: rlock only, check whether containers need maintenance; if not,
      						assign the punc.

      slow path: if maintenance needed, start over by grabbing upgrade_lock.
      						pretend that we still can do fast path (we don't exclude
      						other readers). if we really need do write, atomically
      						upgrade the lock so other readers are excluded.
      come in with rlock. assuming fast path.

      see upgrade_locks()
    */

	boost::shared_lock<boost::shared_mutex> rlock(mtx_container);
	boost::upgrade_lock<boost::shared_mutex> ulock(mtx_container, boost::defer_lock);
	unique_ptr<boost::upgrade_to_unique_lock<boost::shared_mutex>> pwlock
							= nullptr; /* will std::move later */
start:
   	if (containers_.empty()) {
   		/* Whole transform drained. we don't even have pending flushed
   		 * bundles (otherwise there will be containers).
   		 * Create a new container and seal it with @wm.
   		 */


		if (!upgrade_locks(&rlock, &ulock, &pwlock))
			goto start;
   		containers_.emplace_front();
   		containers_.begin()->setPuncSafe(punc);
   		return;
   	}

   	/* assign or update the wm of the most recent container */
   	auto current = containers_.begin();
	if (current->punc_refcnt == PUNCREF_RETRIEVED) {
		/* punc already emitted. start a new container. */

		if (!upgrade_locks(&rlock, &ulock, &pwlock))
			goto start;

		xzl_assert(current->getPuncSafe() && "must have a valid punc ptr");
		xzl_assert(!current->refcnt && "all bundles must be consumed");
		//   		containers_.push_front(bundle_container());
		containers_.emplace_front();
		containers_.begin()->setPuncSafe(punc);
		return;
	}

	/* clean up any dead container. */
	if (current->punc_refcnt == PUNCREF_CONSUMED) {

		if (!upgrade_locks(&rlock, &ulock, &pwlock))
			goto start;

		/* extensive sanity check... */
		xzl_assert(current->getPuncSafe() && "dead container but no valid punc?");
		xzl_assert(current->peekBundleCount() == 0);
		xzl_assert(current->refcnt == 0);

		current = containers_.erase(current);

		if (containers_.empty()) {
			/* after cleanup, no container left. start a new one. */
			containers_.emplace_front();
			containers_.begin()->setPuncSafe(punc);
			return;
		} else {
			/* current now points to the 2nd most recent container. */
			xzl_assert(current->punc_refcnt != 0
					&& "bug: can't have two dead containers");
		}
	}

   	/* should have any type of lock */
   	xzl_assert(rlock.owns_lock()
   						  || ulock.owns_lock() || (pwlock && pwlock->owns_lock()));

   	/* if we have wlock, we don't need to -Safe for individual container.
   	 * but just be conservative...
   	 */
	if (!current->getPuncSafe() || current->punc_refcnt == PUNCREF_ASSIGNED) {
		current->setPuncSafe(punc); /* assign or overwrite punc */
		//current->set_side_info(this->get_side_info());
	} else {
		xzl_bug("bug?");
	}
  }

  /* Any container in this transform needs help?
   *
   * @has_early_punc: OUT. see comment below.
   */
  bool hasWorkOlderThan(ptime const & t, bool* has_early_punc) {
  	/* 1. check tran wm, if more recent (>t), no work. (unemitted punc must>t)
  	 *
  	 * 2. if we've see an unemitted punc in upstream <t, we should do all
  	 * containers in this trans regardless open/close. Since they'll
  	 * block the punc from emitting from this trans.
  	 *
  	 * 3. if no punc <t observed in upstream, we check the oldest
  	 * container in this transform. if open (punc unassigned), no work here.
  	 * if its punc<=t assigned & unretrieved, has work & return punc & stop.
  	 * [ what if punc outstanding? no work? look at next container?]
  	 *
  	 * report any discovered <t punc to caller.
  	 */

  	if (this->GetWatermarkSafe() > t)
  		return false;
  	return true;
  }

  /* Get a bundle/punc from a container whose punc is assigned and <t.
   *
   * Does not consider bundles from open container (should be done
   * by getOneBundle())
   *
   * This also applies to new cascading container design, which may have
   * multiple open containers towards the front of the list, e.g.
   * O-open X-sealed C-(cancelled)
   *
   * possible:
   * O O O X
   * O O C X
   * O X C X
   *
   * impossible:
   * X O O X
   * */

  virtual shared_ptr<BundleBase> getOneBundleOlderThan(ptime const & t,
  		bool* has_pending_punc, int node = -1) {

		*has_pending_punc = false;

  	/* fastest path: if local wm more recent (>t), no work.
  	 * (and any unemitted punc must >t) */
  	if (this->GetWatermarkSafe() > t) {
  	  return nullptr;
  	}

  	/* fast/slow path: there may be some work in this trans. */


  	boost::shared_lock<boost::shared_mutex> rlock(mtx_container);
  	boost::upgrade_lock<boost::shared_mutex> ulock(mtx_container, boost::defer_lock);
  	unique_ptr<boost::upgrade_to_unique_lock<boost::shared_mutex>> pwlock
  							= nullptr; /* will std::move later */

  	shared_ptr<BundleBase> ret;

  	/* if we upgrade rlock->ulock, we start over from here. this is
  	 * because containers_ may have been changed so much during lock upgrade
  	 */
start:

		if (containers_.empty()) {
//			dump_containers("no containers");
			return nullptr;
		}

		/* Start from oldest container (at back()) and look for one available
		 * bundle.
		 * For each container, we first check available data bundles and then
		 * any ready punc. If no luck, move to a newer container.
		 *
		 * When encountering an empty container, check its refcnt.
		 */

		auto it = containers_.rbegin();

		/* remember the 1st container that 1) has cancelled punc and 2) has
		 * some bundles */
		auto itt = containers_.rend();

		while (it != containers_.rend()) {

			if (it->punc_refcnt == PUNCREF_UNDECIDED) {
				/* we don't deal open container (and all later ones) */
				xzl_assert(it->verifyPuncSafe() && "bug?");
				return nullptr;
			}

			/* punc has been canceled, if we see a later container assigned with
			 * a punc, we may return a bundle from *this* container.
			 */
			if (it->punc_refcnt == PUNCREF_CANCELED) {
				if (it->refcnt == 0) { /* dead container, clean it up. */

					if (!upgrade_locks(&rlock, &ulock, &pwlock))
						goto start;

					/* because punc only canceled after bundles are deposited.*/
					xzl_assert(!it->peekBundleCount());
					it ++;
					it = list<bundle_container>::reverse_iterator(containers_.erase(it.base()));
				} else if (it->peekBundleCount() && itt == containers_.rend()) {
					/* bundles available + oldest container w/ punc canceled */
					itt = it;
					it ++;
				} else {
					/* else: could be no bundles in containers and some bundles are
						outstanding... */
					it ++;
				}
				continue;
			}

			/* punc has *ever* been assigned */

			if (it->getPuncSafe()->min_ts > t) {
				/* punc newer than wm. skip all remaining (later) containers */
				return nullptr;
			}

			/* now we see an assigned punc <=wm. */

			*has_pending_punc = true;

  		/* Grab from the oldest container w/ canceled punc. However, when we
  		 * go back, that container may or may not have bundles (other workers
  		 * may have stolen them all). So we check all containers [itt, it) */
			if (itt != containers_.rend()) {
				ret = itt->getBundleSafe();
				xzl_assert(ret);
				if (ret) {
					this->DecBundleCounter();
					return ret;
				}
			}

			/* non-empty container. grab & go */
//			if ((ret = it->getBundleUnsafe())) {
  		if ((ret = it->getBundleSafe())) {  /* see above comment */
				this->DecBundleCounter();
				return ret;
			}

			/* an empty container... */

			/* some data bundles outstanding.
			 * don't wait. move to a newer container. */
			if (it->refcnt != 0) {
#if 0 /* can no longer verify the following due to concurrent readers. unless we lock the container */
//				if (!(it->punc_refcnt != PUNCREF_RETRIEVED
//										&& it->punc_refcnt != PUNCREF_CONSUMED)) {
//					auto refcnt = it->punc_refcnt.load();
//					EE("bug: punc refcnt is %ld", refcnt);
//				}
				xzl_assert(it->punc_refcnt != PUNCREF_RETRIEVED
						&& it->punc_refcnt != PUNCREF_CONSUMED);
#endif
				it ++;
				continue;
			}

			/* an empty container, punc has been assigned (CANCELED case is handled
			 * above). all bundles already consumed.
			 *
			 * we take a peek of punc_refcnt and decide.
			 * */

			auto punc_refcnt = it->punc_refcnt.load();

			/* first, clean any dead container so that it won't block
			 * later punc.
			 */
			if (punc_refcnt == PUNCREF_CONSUMED) { /* bundles+punc consumed */
				/* has to be the oldest container since we can't have more
				 * than one outstanding punc (the second oldest container, whose
				 * punc is assigned, must still holds its punc).
				 */

	   		if (!upgrade_locks(&rlock, &ulock, &pwlock))
	   			goto start;  /* see comment at @start */

				xzl_assert(it == containers_.rbegin());
//				it = containers_.erase(it);   /* erase w/ reverse iterator */
				it ++;
				it = list<bundle_container>::reverse_iterator(containers_.erase(it.base()));
				continue;
			} else if (punc_refcnt == PUNCREF_RETRIEVED) {
				/* punc outstanding. it may become 0 soon but we don't wait. */
				it ++;
				continue;
			} else if (punc_refcnt == PUNCREF_ASSIGNED) { /* punc not emitted */
				if (it == containers_.rbegin()) { /* oldest container. okay to emit */
					long expected = PUNCREF_ASSIGNED;
					if (!it->punc_refcnt.compare_exchange_strong(expected,
							PUNCREF_RETRIEVED)) {
//						bug("bad punc refcnt");
						/* somebody just took & emitted it. move on to next container */
						it ++;
						continue;
					}
//					auto r = it->punc_refcnt.fetch_sub(1);
//					xzl_assert(r == 2); r = r;

					/* XXX: opt: we may check the next newer container and
					 * coalesce punc as needed. benefit unclear though.
					 */
					return it->getPuncSafe();
				} else {
					/* not oldest container. can't emit punc before older ones. */
#if 0 /* not a good idea */
					/* opt: combine unemitted puncs.
					 *
					 * this is dangerous, since punc_refcnt may change under us.
					 * idea: temporarily decrements punc_refcnt to 1 so that
					 * other evaluators believe the punc is outstanding. after
					 * combing, increment to 2.  */
					auto & older = *(std::prev(it));
					if (older->punc_refcnt == 2) {
						older->wm = current->wm;
//						it = containers_.erase(it); /* erase & continue */
						it ++;
						it = list<bundle_container>::reverse_iterator(containers_.erase(it.base()));
					} else
#endif
						it ++;
					continue;
				}
			} else {
				EE("punc refcnt: %ld", punc_refcnt);
				xzl_assert(false && "illegal punc_perfcnt?");
			}

			xzl_assert(false && "bug?");
		} // while

//		dump_containers("no bundles");

		return nullptr;
  }

  /* Retrieve one bundle from @containers_ of this transform.
   * We may get a data bundle or a punc.
   */

  shared_ptr<BundleBase> getOneBundle(int node = -1) {

  	shared_ptr<BundleBase> ret;

  	/* for locking rule, see comment in the above func */

//  	unique_lock<mutex> conlock(mtx_container);
  	boost::shared_lock<boost::shared_mutex> rlock(mtx_container);
  	boost::upgrade_lock<boost::shared_mutex> ulock(mtx_container, boost::defer_lock);
  	unique_ptr<boost::upgrade_to_unique_lock<boost::shared_mutex>> pwlock
  							= nullptr; /* will std::move later */

start:
  	if (containers_.empty()) {
//  		dump_containers("no containers");
  		return nullptr;
  	}

  	/* Start from oldest container (at back()) and look for one available
  	 * bundle.
  	 * For each container, we first check available data bundles and then
  	 * any ready-to-retrieve punc. If no luck, move to a newer container.
  	 *
  	 * when encountering an empty container (no bundles), check its refcnt.
  	 */

  	auto it = containers_.rbegin();

  	while (it != containers_.rend()) {

  		/* non-empty container. grab & go */
  		if ((ret = it->getBundleSafe())) {
  			VV("got one bundle from container %lx", (unsigned long)(&(*it)));
  			this->DecBundleCounter();
  			return ret;
  		}

  		/* Can't get any bundle: an empty container... */

			if (it->punc_refcnt == PUNCREF_UNDECIDED) {
				/* punc to be assigned. since we allow multi "open" containers, we
				 * should continue to the next container... */
				xzl_assert(it->verifyPuncSafe());
				it ++;
				continue;
			}

			/* some data bundles are outstanding.
			 * don't wait. move to a newer container. */
			if (it->refcnt != 0) {
#if 0 /* can no longer verify the following due to concurrent readers. unless we lock the container */
				xzl_assert(it->punc_refcnt != PUNCREF_CONSUMED
							&& it->punc_refcnt != PUNCREF_RETRIEVED);
#endif
				it ++;
				continue;
			}

			/* an empty container, punc has been determined.
			 * all bundles already consumed (refcnt = 0).
			 *
			 * take a peek of punc_refcnt ...
			 * */

			auto punc_refcnt = it->punc_refcnt.load();

			/* first, clean any dead container so that it won't block the next punc.
			 */

			if (punc_refcnt == PUNCREF_CONSUMED) {

				if (!upgrade_locks(&rlock, &ulock, &pwlock))
					goto start;


				/* has to be the oldest container since we can't have more
				 * than one outstanding punc (the second oldest container, whose
				 * punc is assigned, must still holds its punc).
				 */
				xzl_assert(it == containers_.rbegin());
//				it = containers_.erase(it);   /* erase w/ reverse iterator */
				it ++;
				it = list<bundle_container>::reverse_iterator(containers_.erase(it.base()));
				continue;
			} else if (punc_refcnt == PUNCREF_CANCELED) {
				/* erase */

				if (!upgrade_locks(&rlock, &ulock, &pwlock))
					goto start;

				it ++;
				it = list<bundle_container>::reverse_iterator(containers_.erase(it.base()));
			} else if (punc_refcnt == PUNCREF_RETRIEVED) {
				/* punc outstanding. it may become CONSUMED soon but we don't wait. */
				it ++;
				continue;
			} else if (punc_refcnt == PUNCREF_ASSIGNED) { /* punc not emitted */
				if (it == containers_.rbegin()) {
					/* oldest container. (older dead containers should have been
					 * cleaned up now). okay to emit */
//					auto r = it->punc_refcnt.fetch_sub(1);
//					xzl_assert(r == 2);

					long expected = PUNCREF_ASSIGNED;
					if (!it->punc_refcnt.compare_exchange_strong(expected,
							PUNCREF_RETRIEVED)) {
						//bug("bad punc refcnt");
						/* somebody just took & emitted the wm. move to next container */
						it ++;
						continue;
					}
					/* XXX: opt: we may check the next newer container and
					 * coalesce punc as needed. benefit unclear though.
					 */
					return it->getPuncSafe();
				} else {
					/* not oldest container. can't emit punc before older ones doing so. */
#if 0 /* not a good idea */
					/* opt: combine unemitted puncs.
					 *
					 * this is dangerous, since punc_refcnt may change under us.
					 * idea: temporarily decrements punc_refcnt to 1 so that
					 * other evaluators believe the punc is outstanding. after
					 * combing, increment to 2.  */
					auto & older = *(std::prev(it));
					if (older->punc_refcnt == 2) {
						older->wm = current->wm;
//						it = containers_.erase(it); /* erase & continue */
						it ++;
						it = list<bundle_container>::reverse_iterator(containers_.erase(it.base()));
					} else
#endif
						it ++;
					continue;
				}
			} else {
				EE("punc refcnt: %ld", punc_refcnt);
				xzl_assert(false && "illegal punc_perfcnt?");
			}

	  	xzl_bug("bug?");

  	} // while

//  	dump_containers("no bundles");

  	return nullptr;
  }

  /* ----------------------------------------- */

	/*
		recalculate the watermark for a give transform.
		but does not propagate the watermark to downstream.

		default behavior for state-less transform: watermark
		is decided by the upstream transform and the oldest
		pending work (including those in input and being worked on)

		we don't expect this function is called often

		@upstream_wm: the snapshot of the upstream's watermark. If unspecified,
		actively fetch the current wm from the upstream.

		@return: the wm after recalculation.
	*/
  virtual ptime RefreshWatermark(ptime upstream_wm = max_date_time) {
    PValue* v = getFirstInput(); // XXX deal with multiple ins
    xzl_assert(v);
    PTransform *upstream = v->producer;
    string upstream_name;
    xzl_assert(upstream);

    unique_lock<mutex> lock(mtx_watermark);

    // also examine all inflight bundles
    ptime min_ts_flight = max_date_time;
    for (auto && b : inflight_bundles) {
        if (b->min_ts < min_ts_flight)
          min_ts_flight = b->min_ts;
    }

    if (upstream_wm == max_date_time) {
    	upstream_wm = upstream->watermark;  // fetch the current wm
    	upstream_name = upstream->getName() + string("(live)");
    } else {
    	upstream_name = upstream->getName() + string("(snapshot)");
    }

    /* the upstream wm should be earlier than records that are yet to be
     * processed. otherwise, we will observe the wm before processing records
     * that are prior to the wm. */

    if (!((upstream_wm < v->min_ts) && (upstream_wm < min_ts_flight))) {
    	EE("%s: (%lx): "
//				EE("%s: %s (%lx): new watermark: %s \n"
				" \t\t upstream(%s) watermark: %s\n"
				" \t\t input min_ts: %s\n"
				" \t\t min_ts_flight: %s",
//						to_simple_string(boost::posix_time::microsec_clock::local_time()).c_str(), /* current ts */
					name.c_str(), (unsigned long)this,
					upstream_name.c_str(),
					to_simplest_string(upstream_wm).c_str(),
					to_simplest_string(v->min_ts).c_str(),
					to_simplest_string(min_ts_flight).c_str()
			);
    }

    xzl_assert(upstream_wm < v->min_ts);
    xzl_assert(upstream_wm < min_ts_flight);

    ptime pt = min(min_ts_flight, min(v->min_ts, upstream_wm));

    if (watermark > pt) {
        EE("bug: watermark regression. existing %s > new watermark",
            to_simple_string(watermark).c_str());
        EE("%s: new watermark: %s"
           " (input %s upstream %s",
           name.c_str(),
           to_simple_string(pt).c_str(),
           to_simple_string(v->min_ts).c_str(),
           to_simple_string(upstream_wm).c_str());
        xzl_assert(0);
    }

    xzl_assert(watermark <= pt);

    /* it's possible that watermark does not change (=pt) because it is held
     * back by something other than the upstream wm (e.g. input).
     */
    if (watermark < pt) {
//    if (watermark <= pt) {  /* debugging */
      watermark = pt;
#if 1
      /* debugging */
      if (name == string("[wc-mapper]"))
      	EE("%s: (%lx): new watermark: %s \n"
//				EE("%s: %s (%lx): new watermark: %s \n"
					" \t\t upstream(%s) watermark: %s\n"
					" \t\t input min_ts: %s\n"
					" \t\t min_ts_flight: %s",
//						to_simple_string(boost::posix_time::microsec_clock::local_time()).c_str(), /* current ts */
						name.c_str(), (unsigned long)this,
						to_simplest_string(watermark).c_str(),
						upstream_name.c_str(),
						to_simplest_string(upstream_wm).c_str(),
						to_simplest_string(v->min_ts).c_str(),
						to_simplest_string(min_ts_flight).c_str()
				);
#endif
    }
    return pt;
  }

   /* Actions to take when the transform sees a new snapshot of
    * upstream's watermark. Often, the transform reacts by invoking the
    * corresponding evaluator, e.g. for purging internal state.
    *
    * The action does not necessarily depends the local watermark; it
    * should not depends on upstream's current wm, which may have advanced.
    *
    * default: do nothing
    *
    * @up_wm: the upstream watermark
    XXX rename to OnUpstreamWatermarkChange()
    */
//  virtual void OnNewUpstreamWatermark(int nodeid,
//  		EvaluationBundleContext *c, ptime up_wm) { }

  virtual void ExecEvaluator(int nodeid, EvaluationBundleContext *c,
  		shared_ptr<BundleBase> bundle_ptr = nullptr) = 0;

  // xzl: protected ctor -- only subclasses can be instantiated.
  std::string name;

  struct Statstics {
  	const char * name;
  	double lmbps, lmrps; /* in last interval */
  	double mbps, mrps; 	/* total avg */
  };

  /* @return: is stat filled */
  virtual bool ReportStatistics(Statstics* stat) { return false; }  /* to be overridden */
protected:

	PTransform() { }

	PTransform(std::string name) : name(name) {
		bundle_counter_ = 0; /* can't be put in init list? */
	}

//private:
public:
  // watermark. must be updated atomically. can't be protected by
	// container mutex
  ptime watermark = min_date_time;
  mutex mtx_watermark;

public:
  /* return a copy */
  ptime GetWatermark() {
  	unique_lock<mutex> lock(mtx_watermark);
  	return watermark;
  }

  /* advance trans wm atomically.
   * return: false if @new_wm == trans wm. true if update okay.
   * panic if watermark regression.
   */
  bool SetWatermark(ptime const & new_wm) {
//  bool AdvanceWatermarkAtomic(ptime const & new_wm) {
  	unique_lock<mutex> lock(mtx_watermark);
  	//assert(new_wm >= this->watermark);
  	if (new_wm == this->watermark)
  		return false;
	this->watermark = new_wm;
	return true;
  }

  /* return a copy */
  ptime GetWatermarkSafe() {
  	unique_lock<mutex> lock(mtx_watermark);
  	return watermark;
  }

  /* advance trans wm atomically.
   * return: false if @new_wm == trans wm. true if update okay.
   * panic if watermark regression.
   */
  bool SetWatermarkSafe(ptime const & new_wm) {
//  bool AdvanceWatermarkAtomic(ptime const & new_wm) {
  	unique_lock<mutex> lock(mtx_watermark);

  	if (new_wm < this->watermark) {
  		cout << "wm regression: ";
  		cout << "new:" << new_wm;
  		cout << " old:" << this->watermark << endl;
  	}

  	xzl_assert(new_wm >= this->watermark);

  	if (new_wm == this->watermark)
  		return false;

	this->watermark = new_wm;
	return true;
  }

  #include "Transforms-multi.h"
};

////////////////////////////////////////////////////////////////////////

// xzl: do nothing but pass through.
// to construct the PCollection, we must use template. XXX remove template
template<class T>
//class PTransformNop : public PTransform <PTransformNop> {
class PTransformNop : public PTransform  {
public:
	PTransformNop() : PTransform("nop") { }

	virtual ~PTransformNop() { }
};

/////////////////////////////////////////////////////////////////////

#include "DoFn.h"

/////////////////////////////////////////////////////////////////////

/* used by concrete transforms */

#define DECL_EVAL2(T) \
	using TransformT = T<InputT,OutputT>; \
	using EvalT = T##Evaluator<InputT,OutputT>; \
	template<typename I,typename O> class EvalT

class EvaluationBundleContext;

/////////////////////////////////////////////////////////////////////

#endif // TRANSFORMS_H
