/*
 * Transforms-multi.h
 *
 *  Created on: Jan 19, 2017
 *      Author: xzl
 *
 *  hym's implementation; separate from Transforms.h and is included from there
 */

#ifndef TRANSFORMS_MULTI_H_
#define TRANSFORMS_MULTI_H_

  //hym: join deposits bundles to its left_containers_out or right_containers_out
  //hym: here downt is Join itself
  void  depositBundleLR(PTransform *downt,
  		shared_ptr<BundleBase> input_bundle,
			vector<shared_ptr<BundleBase>> const & output_bundles, int node = -1)
  {
  	bool try_open = false; /* have we tried open downstream containers? */
	//std::cout << __FILE__ << ": " <<  __LINE__ << std::endl;
  	//std::cout << "trans is: " << downt->getName() << std::endl;
	// try_open = try_open;
		xzl_assert(downt);
  	bundle_container *upcon = localBundleToContainer(input_bundle);
  	xzl_assert(upcon); /* enclosing container */

	//XXX rw
	//unique_lock<mutex> conlock(mtx_container);
	// write lock
	//boost::unique_lock<boost::shared_mutex> down_conlock(mtx_container);
deposit:
      {
	//unique_lock<mutex> conlock(mtx_container);

		//{
//			unique_lock<mutex> conlock(mtx_container);
			if (upcon->downstream) {
				/* fast path. NB @downstream is atomic w/ seq consistency so it's
				 * safe. */
				downt->depositBundlesToContainer(upcon->downstream, output_bundles, node);
				downt->IncBundleCounter(output_bundles.size());
				//std::cout << __FILE__ << ": " << __LINE__ << "join deposits bundle to its out: " << upcon->get_side_info() << std::endl;
				return;
			}
		//}
#if 0
	if(upcon->get_side_info() == 1){
		unique_lock<mutex> conlock(left_mtx_container_out); //Join itself out containers
		if (upcon->downstream) {
			downt->depositBundlesToContainer(upcon->downstream, output_bundles, node);
			downt->IncBundleCounter(output_bundles.size());
			return;
		}

	}else if(upcon->get_side_info() == 2){
		unique_lock<mutex> conlock(right_mtx_container_out); //Join itself out containers
		if (upcon->downstream) {
			downt->depositBundlesToContainer(upcon->downstream, output_bundles, node);
			downt->IncBundleCounter(output_bundles.size());
			return;
		}
	}else{
		xzl_assert(false && "unknown side info in join");
	}
#endif
	// try_open = try_open;
  	xzl_assert(!try_open && "bug: already tried to open");

  	/* it's possible somebody slipped in & opened the downstream container
  	 * already. So ret == 0 then.
  	 */
/*
//  	int ret =
  			openDownstreamContainers(upcon, downt);
//  	xzl_assert(ret && "no downstream container opened?");
*/

/*
	if(this->get_side_info() == 0 || this->get_side_info() == 3){ //defaut transforms and join
		openDownstreamContainers(upcon, downt);
	}else if(this->get_side_info() == 1){ //Simplemapper 1, left
		openDownstreamContainers_left(upcon, downt); //downt should be join, open a container in join's left_containers_in
	}else if(this->get_side_info() == 2){ //Simplemapper 2, right
		openDownstreamContainers_right(upcon, downt); //downt should be join, open a container in join's right_containers_in
	}else{
		xzl_assert("Bug: Unknown know side info");
	}
*/
	if(upcon->get_side_info() == 1){
		//std::cout << __FILE__ << ": " <<  __LINE__ << std::endl;
		openJoinOutContainer_left(upcon, downt);
		//std::cout << __FILE__ << ": " <<  __LINE__ << std::endl;
	}else if(upcon->get_side_info() == 2){
		//std::cout << __FILE__ << ": " <<  __LINE__ << std::endl;
		openJoinOutContainer_right(upcon, downt);
		//std::cout << __FILE__ << ": " <<  __LINE__ << std::endl;
	}else{
		xzl_assert("Join: wrong side info");
	}
    }//end lock
	xzl_assert(upcon->downstream);
  	try_open = true;
  	goto deposit;
  }

  //hym: Join deposit one punc to its left_containers_out or right_containers_out
  // here downt is join itself
    void depositPuncLR(PTransform* downt,
    		shared_ptr<Punc> const & input_punc, shared_ptr<Punc> punc, int node = -1)
    {
    	bool try_open = false;

    	xzl_assert(downt);
    	bundle_container *upcon = localBundleToContainer(input_punc);
    	xzl_assert(upcon); /* enclosing container */

  	//XXX rw
  	//unique_lock<mutex> conlock(mtx_container);
  	// write lock
  	//boost::unique_lock<boost::shared_mutex> down_conlock(mtx_container);
  deposit:
       {
  //	unique_lock<mutex> conlock(mtx_container);
  	//{
  //			unique_lock<mutex> conlock(mtx_container);
  			if (upcon->downstream) { /* fast path. NB @downstream is atomic, implying fence */
  				downt->depositOnePuncToContainer(upcon->downstream, punc, node);
  				//upcon->has_punc = 1; //the upcon has a container already
  				upcon->downstream.load()->has_punc = 1;
  				//std::cout << __FILE__ << ": " <<  __LINE__ << " " << this->getName() << " deposit a punc to " << downt->getName() << " to side: " << upcon->get_side_info()  << std::endl;
  				return;
  			}
  	//}
  #if 0
  		if(upcon->get_side_info() == 1){
  			unique_lock<mutex> conlock(left_mtx_container_out); //Join itself out containers
  			if (upcon->downstream) {
  				downt->depositOnePuncToContainer(upcon->downstream, punc, node);
  				upcon->downstream.load()->has_punc = 1;
  				return;
  			}

  		}else if(upcon->get_side_info() == 2){
  			unique_lock<mutex> conlock(right_mtx_container_out); //Join itself out containers
  			if (upcon->downstream) {
  				downt->depositOnePuncToContainer(upcon->downstream, punc, node);
  				upcon->downstream.load()->has_punc = 1;
  				return;
  			}
  		}else{
  			xzl_assert(false && "unknown side info in join");
  		}
  #endif
  		// try_open = try_open;
  		xzl_assert(!try_open && "bug: already tried to open");

  		//openDownstreamContainers(upcon, downt);
  		//hym:XXX
  	/*	if(this->get_side_info() == 0 || this->get_side_info() == 3){ //defaut transforms and join
  			openDownstreamContainers(upcon, downt);
  		}else if(this->get_side_info() == 1){ //Simplemapper 1, left
  			openDownstreamContainers_left(upcon, downt); //downt should be join, open a container in join's left_containers_in
  		}else if(this->get_side_info() == 2){ //Simplemapper 2, right
  			openDownstreamContainers_right(upcon, downt); //downt should be join, open a container in join's right_containers_in
  		}else{
  			xzl_assert("Bug: Unknown know side info");
  		}
  	*/

  		if(upcon->get_side_info() == 1){
  			openJoinOutContainer_left(upcon, downt);
  		}else if(upcon->get_side_info() == 2){
  			openJoinOutContainer_right(upcon, downt);
  		}else{
  			xzl_assert("Join: wrong side info");
  		}
  	}//end lock
  		xzl_assert(upcon->downstream);
  		try_open = true;
  		goto deposit;
    }

    //hym: this is used to open a container in Join's left_containers_in
    //hym: So we don't need to modify depositBundlesDownstream() any more
    //hym: simple mapper will call this function
    //hym: downt is Join here
    int left_num = 0;
    int openDownstreamContainers_left(bundle_container * const upcon,
    		PTransform *downt)
    {
    	xzl_assert(downt && upcon);
  	//std::cout << __FILE__ << ": " <<  __LINE__ << std::endl;
  	//std::cout << "this trans is: " << this->getName() << std::endl;
  	//std::cout << "downt is: " << downt->getName() << std::endl;
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
  #if 0
  	unique_lock<mutex> conlock(mtx_container);

  	//here this is SimpleMapper
    	unique_lock<mutex> down_conlock(downt->mtx_container);
    	//unique_lock<mutex> down_conlock(downt->left_mtx_container_in);
  #endif
  	/* reader lock for *this* trans */
	boost::shared_lock<boost::shared_mutex> conlock(mtx_container);
	/* writer lock for downstream trans */
	boost::unique_lock<boost::shared_mutex> down_conlock(downt->mtx_container);

    	if (upcon->downstream) /* check again */
    		return 0;
  	//std::cout << __FILE__ << ": " <<  __LINE__ << std::endl;
    	int cnt = 0;
    	bool miss = false; (void)sizeof(miss);
  	// miss = miss;
  	//std::cout << __FILE__ << ": " <<  __LINE__ << " this trans is: " << this->getName() << " downt is " << downt->getName() <<std::endl;
    	for (auto it = this->containers_.rbegin();
    			it != this->containers_.rend(); it ++) {

  		//std::cout << __FILE__ << ": " <<  __LINE__ << std::endl;
    		/* also lock all downstream containers. good idea? */
  //  		unique_lock<mutex> down_conlock(downt->mtx_container);

    		if (!it->downstream) {
    			miss = true;
  			//std::cout << __FILE__ << ": " <<  __LINE__ << std::endl;
    			/* alloc & link downstream container. NB: it->downstream is an atomic
    			 * pointer, so it implies a fence here. */
    			downt->left_containers_in.emplace_front();
    			it->downstream = &(downt->left_containers_in.front());

  			//hym: XXX
  			//assigne the side_info to this new container
  			//it->downstream->side_info = this->get_side_info(); //0->containrs_, 1->left_containers_in, 2->right_containers_in
    			it->downstream.load()->set_side_info(this->get_side_info());
  			xzl_assert(it->downstream.load()->get_side_info() == 1);	//left

  			cnt ++;
    		} else { /* downstream link exists */
    			assert (!miss && "bug? an older container misses downstream link.");
  #ifdef DEBUG
    			bool found = false;
    			for (auto & con : downt->left_containers_in) {
    				if (&con == it->downstream) {
    					found = true;
    					break;
    				}
    			}
    			/* NB: downstream container, once received & processed the punc
    			 * from the upstream, may be cleaned up prior to the upstream
    			 * container. In this case, a dangling downstream pointer is possible.
    			 */
    			if (!found && it->punc_refcnt != PUNCREF_CONSUMED) {
    				/*
  				E("%s: downstream (%s) container does not exist",
    						this->name.c_str(), downt->name.c_str());
  					cout << " -------------------------------\n";
  					dump_containers("bug");
  					downt->dump_containers("bug");
  					cout << " -------------------------------\n\n";
    				abort();
  				*/
  				xzl_assert(false && "downstream container does not exist");
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


    //hym: this is used to open a container in Join's right_containers_in
    //hym: SimpleMapper will call this function
    //hym: downt is Join here
    int openDownstreamContainers_right(bundle_container * const upcon,
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
  #if 0
  	unique_lock<mutex> conlock(mtx_container);

  	//here this is SimppleMapper
  	// downt is Join
    	unique_lock<mutex> down_conlock(downt->mtx_container);
    	//unique_lock<mutex> down_conlock(downt->right_mtx_container_in);
  #endif
  	/* reader lock for *this* trans */
	boost::shared_lock<boost::shared_mutex> conlock(mtx_container);
	/* writer lock for downstream trans */
	boost::unique_lock<boost::shared_mutex> down_conlock(downt->mtx_container);
    	
	if (upcon->downstream) /* check again */
    		return 0;

    	int cnt = 0;
    	bool miss = false; (void)sizeof(miss);
  	// miss = miss;
    	for (auto it = this->containers_.rbegin();
    			it != this->containers_.rend(); it ++) {

    		/* also lock all downstream containers. good idea? */
  //  		unique_lock<mutex> down_conlock(downt->mtx_container);

    		if (!it->downstream) {
    			miss = true;
    			//std::cout << __FILE__ << ": " << __LINE__ << std::endl;
  			/* alloc & link downstream container. NB: it->downstream is an atomic
    			 * pointer, so it implies a fence here. */
    			downt->right_containers_in.emplace_front();
    			it->downstream = &(downt->right_containers_in.front());

  			//hym: XXX
  			//assigne the side_info to this new container
  			//it->downstream->side_info = this->get_side_info(); //0->containrs_, 1->left_containers_in, 2->right_containers_in
  			it->downstream.load()->set_side_info(this->get_side_info());
  			xzl_assert(it->downstream.load()->get_side_info() == 2); //right

  			cnt ++;
    		} else { /* downstream link exists */
    			assert (!miss && "bug? an older container misses downstream link.");
  #ifdef DEBUG
    			bool found = false;
    			for (auto & con : downt->right_containers_in) {
    				if (&con == it->downstream) {
    					found = true;
    					break;
    				}
    			}
    			/* NB: downstream container, once received & processed the punc
    			 * from the upstream, may be cleaned up prior to the upstream
    			 * container. In this case, a dangling downstream pointer is possible.
    			 */
    			if (!found && it->punc_refcnt != PUNCREF_CONSUMED) {
    				/*
  				E("%s: downstream (%s) container does not exist",
    						this->name.c_str(), downt->name.c_str());
  					cout << " -------------------------------\n";
  					dump_containers("bug");
  					downt->dump_containers("bug");
  					cout << " -------------------------------\n\n";
    				abort();
  				*/

  				xzl_assert(false && "downstream container does not exist");
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


    /*****************************************************************
     * hym:    For type 4 transform: the following transform of Join
     * 			        1
     *        0                SimpleMapper1         3           4
     * UnbundedInMem<long> ==>                 ==> Join ==> type 4 trans ==> type 5 trans
     *             	           SimpleMapper2
     *			        2
     * type 4 trans have:
     *	-- left_unorderded_containers
     *	-- right_unordered_containers
     *	-- ordered_containers
     *	-- left_wm and right_wm
     * Join:
     *	-- deposit bundles/punc to 4's left/right_unordered_containers
     *	-- move containers in left/right_unordered_containers to ordered_containers
     * type 4 trans:
     * 	-- grab bundles/puncs from ordered_containers first
     *	-- else grab bundles from left/right_unordered_containers
     *	-- don't grab punc from left/right_unorderd_containers, should block getting punc..
     *	-- TODO: deposit bundles/puncs to type 5 trans' unordered_list and ordered_list
     ****************************************************************/

      //TODO: get bundles/puncs from orderd_containers
      shared_ptr<BundleBase> getOneBundleOlderThan_ordered_4(ptime const & t,
      	bool* has_pending_punc, int node = -1)
      {
    	*has_pending_punc = false;

    	//std::cout << __FILE__ << ": " <<  __LINE__ << " This trans is: " << this->getName()  << std::endl;

      	/* fast path: if local wm more recent (>t), no work.
      	 * (and any unemitted punc must >t) */
      	if (this->GetWatermark() > t) {
      	  return nullptr;
      	}

      	/* slowpath: there may be some work in this trans. lock containers */

      	shared_ptr<BundleBase> ret;
    		//unique_lock<mutex> conlock(mtx_ordered_containers);

    		//XXX rw
    		//unique_lock<mutex> conlock(mtx_container);
    		boost::shared_lock<boost::shared_mutex> rlock(mtx_container);
    		boost::upgrade_lock<boost::shared_mutex> ulock(mtx_container, boost::defer_lock);
    		unique_ptr<boost::upgrade_to_unique_lock<boost::shared_mutex>> pwlock
    								= nullptr; /* will std::move later */
    start:
    		if (ordered_containers.empty()) {
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

    		auto it = ordered_containers.rbegin();
    		/* to remember the 1st container that 1) has invalidated punc; 2) has
    		 * some bundles
    		 */
    		auto itt = ordered_containers.rend();

    		while (it != ordered_containers.rend()) {

    			if (it->punc_refcnt == PUNCREF_UNDECIDED) { /* we don't deal open container */
    				xzl_assert(it->verifyPuncSafe() && "bug?");
    				return nullptr;
    			}

    			/* punc has been canceled, if we see a later container assigned with
    			 * a punc, we may return a bundle from *this* container.
    			 */
    			if (it->punc_refcnt == PUNCREF_CANCELED) {
    				if (it->refcnt == 0) { /* dead container, clean it up. */
    					/* because punc only canceled after bundles are deposited.*/

    					//XXX rw
    					if (!upgrade_locks(&rlock, &ulock, &pwlock)){
    						goto start;
    					}

    					xzl_assert(it->bundles.empty());
    					it ++;
    					it = list<bundle_container>::reverse_iterator(ordered_containers.erase(it.base()));
    					std::cout << __FILE__ << __LINE__ << "containers_.erase" << std::endl;
    				} else if (!it->bundles.empty() && itt == ordered_containers.rend()) {
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

    			//if (it->punc->min_ts > t) { /* punc newer than wm. */
    			if (it->getPuncSafe()->min_ts > t) { /* punc newer than wm. */
    				return nullptr;
    			}

    			/* now we see a pending punc <=wm. */

    			*has_pending_punc = true;

    			/* grab from the oldest container w/ canceled punc */
    			if (itt != ordered_containers.rend()) {
    				//ret = itt->getBundleUnsafe();
    				ret = itt->getBundleSafe();
    				xzl_assert(ret);
    				this->DecBundleCounter();
    				return ret;
    			}

    			/* non-empty container. grab & go */
    			//if ((ret = it->getBundleUnsafe())) {
    			if ((ret = it->getBundleSafe())) {
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
    //					E("bug: punc refcnt is %ld", refcnt);
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

    				//XXX rw
    				if (!upgrade_locks(&rlock, &ulock, &pwlock)){
    					goto start;
    				}

    				xzl_assert(it == ordered_containers.rbegin());
    //				it = containers_.erase(it);   /* erase w/ reverse iterator */
    				it ++;
    				it = list<bundle_container>::reverse_iterator(ordered_containers.erase(it.base()));
    #ifdef DEBUG
    				std::cout << __FILE__ << __LINE__ << "containers_.erase " << this->getName()  << std::endl;
    #endif
    				continue;
    			} else if (punc_refcnt == PUNCREF_RETRIEVED) {
    				/* punc outstanding. it may become 0 soon but we don't wait. */
    				it ++;
    				continue;
    			} else if (punc_refcnt == PUNCREF_ASSIGNED) { /* punc not emitted */
    				if (it == ordered_containers.rbegin()) { /* oldest container. okay to emit */
    					long expected = PUNCREF_ASSIGNED;
    					if (!it->punc_refcnt.compare_exchange_strong(expected,
    							PUNCREF_RETRIEVED)) {
    						//bug("bad punc refcnt");
    						it ++;
    						continue;
    					}
    //					auto r = it->punc_refcnt.fetch_sub(1);
    //					xzl_assert(r == 2); r = r;

    					/* XXX: opt: we may check the next newer container and
    					 * coalesce punc as needed. benefit unclear though.
    					 */
    					//return it->punc;
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
    				//E("punc refcnt: %ld", punc_refcnt);
    				xzl_assert(false && "illegal punc_perfcnt?");
    			}

    			xzl_assert(false && "bug?");
    		} // while

    //		dump_containers("no bundles");

    		return nullptr;
      }

      shared_ptr<BundleBase> getOneBundleOlderThan_left_unordered_4(ptime const & t,
      	bool* has_pending_punc, int node = -1)
      {
      	//TODO: get bundles from left_unordered_containers
    	//XXX   Should NOT get punc here
    	*has_pending_punc = false;

    	//std::cout << __FILE__ << ": " <<  __LINE__ << " This trans is: " << this->getName()  << std::endl;

      	/* fast path: if local wm more recent (>t), no work.
      	 * (and any unemitted punc must >t) */
      	if (this->GetWatermark() > t) {
      	  return nullptr;
      	}

      	/* slowpath: there may be some work in this trans. lock containers */

      	shared_ptr<BundleBase> ret;
    		//unique_lock<mutex> conlock(left_mtx_unordered_containers);

    		//XXX rw
    		//unique_lock<mutex> conlock(mtx_container);
    		boost::shared_lock<boost::shared_mutex> rlock(mtx_container);
    		boost::upgrade_lock<boost::shared_mutex> ulock(mtx_container, boost::defer_lock);
    		unique_ptr<boost::upgrade_to_unique_lock<boost::shared_mutex>> pwlock
    								= nullptr; /* will std::move later */
    //start:
    		if (left_unordered_containers.empty()) {
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

    		auto it = left_unordered_containers.rbegin();
    		/* to remember the 1st container that 1) has invalidated punc; 2) has
    		 * some bundles
    		 */
    		auto itt = left_unordered_containers.rend();

    		while (it != left_unordered_containers.rend()) {

    			if (it->punc_refcnt == PUNCREF_UNDECIDED) { /* we don't deal open container */
    				xzl_assert(it->verifyPuncSafe() && "bug?");
    				return nullptr;
    			}

    			/* punc has been canceled, if we see a later container assigned with
    			 * a punc, we may return a bundle from *this* container.
    			 */
    			if (it->punc_refcnt == PUNCREF_CANCELED) {
    				if (it->refcnt == 0) { /* dead container, clean it up. */
    					/* because punc only canceled after bundles are deposited.*/
    					xzl_assert(it->bundles.empty());
    					it ++;
    					//XXX should not destroy container here
    					continue;
    					xzl_assert(false && "should not return Punc here");
    					/*
    					it = list<bundle_container>::reverse_iterator(left_unordered_containers.erase(it.base()));
    					std::cout << __FILE__ << __LINE__ << "containers_.erase" << std::endl;
    					*/
    				} else if (!it->bundles.empty() && itt == left_unordered_containers.rend()) {
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

    			//if (it->punc->min_ts > t) { /* punc newer than wm. */
    			if (it->getPuncSafe()->min_ts > t) { /* punc newer than wm. */
    				return nullptr;
    			}

    			/* now we see a pending punc <=wm. */

    			*has_pending_punc = true;

    			/* grab from the oldest container w/ canceled punc */
    			if (itt != left_unordered_containers.rend()) {
    				//ret = itt->getBundleUnsafe();
    				ret = itt->getBundleSafe();
    				xzl_assert(ret);
    				this->DecBundleCounter();
    				return ret;
    			}

    			/* non-empty container. grab & go */
    			//if ((ret = it->getBundleUnsafe())) {
    			if ((ret = it->getBundleSafe())) {
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
    //					E("bug: punc refcnt is %ld", refcnt);
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
    				xzl_assert(it == left_unordered_containers.rbegin());
    //				it = containers_.erase(it);   /* erase w/ reverse iterator */
    				it ++;
    				//XXX Should NOT destroy containers here
    				//xzl_assert(false && "should not return Punc here");
    				/*
    				it = list<bundle_container>::reverse_iterator(left_unordered_containers.erase(it.base()));
    				std::cout << __FILE__ << __LINE__ << "containers_.erase " << this->getName()  << std::endl;
    				*/
    				continue;
    			} else if (punc_refcnt == PUNCREF_RETRIEVED) {
    				/* punc outstanding. it may become 0 soon but we don't wait. */
    				it ++;
    				continue;
    			} else if (punc_refcnt == PUNCREF_ASSIGNED) { /* punc not emitted */
    				if (it == left_unordered_containers.rbegin()) { /* oldest container. okay to emit */
    					//TODO: we should not return punc here
    					it++;
    					continue;
    				#if 0
    					long expected = PUNCREF_ASSIGNED;
    					if (!it->punc_refcnt.compare_exchange_strong(expected,
    							PUNCREF_RETRIEVED)) {
    						bug("bad punc refcnt");
    					}
    //					auto r = it->punc_refcnt.fetch_sub(1);
    //					xzl_assert(r == 2); r = r;

    					/* XXX: opt: we may check the next newer container and
    					 * coalesce punc as needed. benefit unclear though.
    					 */
    					//XXX should not return punc here
    					xzl_assert(false && "should not return Punc here");
    					return it->punc;
    				#endif
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
    				//E("punc refcnt: %ld", punc_refcnt);
    				std::cout << "punc_refcnt is " << punc_refcnt << std::endl;
    				xzl_assert(false && "illegal punc_perfcnt?");
    			}

    			xzl_assert(false && "bug?");
    		} // while

    //		dump_containers("no bundles");

    		return nullptr;
      }

      //hym: Join deposit a bundle to left_containers_out or
       //hym: here downt is Join it self
       int openJoinOutContainer_left(bundle_container * const upcon,
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
     #if 0
     	//unique_lock<mutex> conlock(mtx_container);
       	unique_lock<mutex> conlock(left_mtx_container_in);

     	// here downt is join actually
       	unique_lock<mutex> down_conlock(downt->left_mtx_container_out);
     #endif
       	/* reader lock for *this* trans */
	boost::shared_lock<boost::shared_mutex> conlock(mtx_container);
	/* writer lock for downstream trans */
	boost::unique_lock<boost::shared_mutex> down_conlock(downt->mtx_container);
	
	if (upcon->downstream) /* check again */
       		return 0;

       	int cnt = 0;
       	bool miss = false;
     	// miss = miss;
       	for (auto it = this->left_containers_in.rbegin();
       			it != this->left_containers_in.rend(); it ++) {

       		/* also lock all downstream containers. good idea? */
     //  		unique_lock<mutex> down_conlock(downt->mtx_container);

       		if (!it->downstream) {
       			miss = true;
       			/* alloc & link downstream container. NB: it->downstream is an atomic
       			 * pointer, so it implies a fence here. */
       			downt->left_containers_out.emplace_front();
       			it->downstream = &(downt->left_containers_out.front());

     			//hym: XXX
     			//assigne the side_info to this new container
     			//it->downstream->side_info = this->get_side_info(); //0->containrs_, 1->left_containers_in, 2->right_containers_in
       			it->downstream.load()->set_side_info(this->get_side_info());
     			xzl_assert(it->downstream.load()->get_side_info() == 1);	//left

     			cnt ++;
       		} else { /* downstream link exists */
       			xzl_assert (!miss && "bug? an older container misses downstream link.");
     #ifdef DEBUG
       			bool found = false;
       			for (auto & con : downt->left_containers_out) {
       				if (&con == it->downstream) {
       					found = true;
       					break;
       				}
       			}
       			/* NB: downstream container, once received & processed the punc
       			 * from the upstream, may be cleaned up prior to the upstream
       			 * container. In this case, a dangling downstream pointer is possible.
       			 */
       			if (!found && it->punc_refcnt != PUNCREF_CONSUMED) {
     			        /*
       				E("%s: downstream (%s) container does not exist",
       						this->name.c_str(), downt->name.c_str());
     					cout << " -------------------------------\n";
     					dump_containers("bug");
     					downt->dump_containers("bug");
     					cout << " -------------------------------\n\n";
       				abort();
     				*/
     				xzl_assert(false && "downstream container does not exist");
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

       //hym: Join deposit a bundle to right_containers_out or
       //hym: here downt is Join it self
       int openJoinOutContainer_right(bundle_container * const upcon,
       		PTransform *downt)
       {
     	//std::cout << __FILE__ << ": " <<  __LINE__ << std::endl;
     	xzl_assert(downt && upcon);

     	//std::cout << __FILE__ << ": " <<  __LINE__ << std::endl;
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
     #if 0
     	//unique_lock<mutex> conlock(mtx_container);
       	unique_lock<mutex> conlock(right_mtx_container_in);

     	//std::cout << __FILE__ << ": " <<  __LINE__ << std::endl;

     	//unique_lock<mutex> down_conlock(downt->mtx_container);
     	unique_lock<mutex> down_conlock(downt->right_mtx_container_out);
     #endif
     	//std::cout << __FILE__ << ": " <<  __LINE__ << "this is " << this->getName() << " downt is " << downt->getName() << std::endl;
  	/* reader lock for *this* trans */
	boost::shared_lock<boost::shared_mutex> conlock(mtx_container);
	/* writer lock for downstream trans */
	boost::unique_lock<boost::shared_mutex> down_conlock(downt->mtx_container);
       	
	if (upcon->downstream) /* check again */
       		return 0;

       	int cnt = 0;
       	bool miss = false; (void)sizeof(miss);
     	// miss = miss;
       	for (auto it = this->right_containers_in.rbegin();
       			it != this->right_containers_in.rend(); it ++) {

     		//std::cout << __FILE__ << ": " <<  __LINE__ << std::endl;
       		/* also lock all downstream containers. good idea? */
     //  		unique_lock<mutex> down_conlock(downt->mtx_container);

       		if (!it->downstream) {
       			miss = true;
       			/* alloc & link downstream container. NB: it->downstream is an atomic
       			 * pointer, so it implies a fence here. */
       			downt->right_containers_out.emplace_front();
       			it->downstream = &(downt->right_containers_out.front());

     			//hym: XXX
     			//assigne the side_info to this new container
     			//it->downstream->side_info = this->get_side_info(); //0->containrs_, 1->left_containers_in, 2->right_containers_in
       			it->downstream.load()->set_side_info(this->get_side_info());
     			xzl_assert(it->downstream.load()->get_side_info() == 1);	//left

     			cnt ++;
       		} else { /* downstream link exists */
       			assert (!miss && "bug? an older container misses downstream link.");
     #ifdef DEBUG
       			bool found = false;
       			for (auto & con : downt->right_containers_out) {
       				if (&con == it->downstream) {
       					found = true;
       					break;
       				}
       			}
       			/* NB: downstream container, once received & processed the punc
       			 * from the upstream, may be cleaned up prior to the upstream
       			 * container. In this case, a dangling downstream pointer is possible.
       			 */
       			if (!found && it->punc_refcnt != PUNCREF_CONSUMED) {
       				/*
     				E("%s: downstream (%s) container does not exist",
       						this->name.c_str(), downt->name.c_str());
     					cout << " -------------------------------\n";
     					dump_containers("bug");
     					downt->dump_containers("bug");
     					cout << " -------------------------------\n\n";
       				abort();
     				*/
     				xzl_assert(false && "downstream container does not exist");

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

       shared_ptr<BundleBase> getOneBundleOlderThan_right_unordered_4(ptime const & t,
       	bool* has_pending_punc, int node = -1)
       {
       	//TODO: get bundles from right_unordered_containers
     	//XXX      Should NOT get punc here

     	*has_pending_punc = false;

     	//std::cout << __FILE__ << ": " <<  __LINE__ << " This trans is: " << this->getName()  << std::endl;

       	/* fast path: if local wm more recent (>t), no work.
       	 * (and any unemitted punc must >t) */
       	if (this->GetWatermark() > t) {
       	  return nullptr;
       	}

       	/* slowpath: there may be some work in this trans. lock containers */

       	shared_ptr<BundleBase> ret;
     		//unique_lock<mutex> conlock(right_mtx_unordered_containers);
     		//XXX rw
     		//unique_lock<mutex> conlock(mtx_container);
     		boost::shared_lock<boost::shared_mutex> rlock(mtx_container);
     		boost::upgrade_lock<boost::shared_mutex> ulock(mtx_container, boost::defer_lock);
     		unique_ptr<boost::upgrade_to_unique_lock<boost::shared_mutex>> pwlock
     								= nullptr; /* will std::move later */
     //start:

     		if (right_unordered_containers.empty()) {
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

     		auto it = right_unordered_containers.rbegin();
     		/* to remember the 1st container that 1) has invalidated punc; 2) has
     		 * some bundles
     		 */
     		auto itt = right_unordered_containers.rend();

     		while (it != right_unordered_containers.rend()) {

     			if (it->punc_refcnt == PUNCREF_UNDECIDED) { /* we don't deal open container */
     				xzl_assert(it->verifyPuncSafe() && "bug?");
     				return nullptr;
     			}

     			/* punc has been canceled, if we see a later container assigned with
     			 * a punc, we may return a bundle from *this* container.
     			 */
     			if (it->punc_refcnt == PUNCREF_CANCELED) {
     				if (it->refcnt == 0) { /* dead container, clean it up. */
     					/* because punc only canceled after bundles are deposited.*/
     					xzl_assert(it->bundles.empty());
     					it ++;
     					//XXX should not destroy container here
     					continue;
     					xzl_assert(false && "should not return Punc here");
     					/*
     					it = list<bundle_container>::reverse_iterator(left_unordered_containers.erase(it.base()));
     					std::cout << __FILE__ << __LINE__ << "containers_.erase" << std::endl;
     					*/
     				} else if (!it->bundles.empty() && itt == right_unordered_containers.rend()) {
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

     			//if (it->punc->min_ts > t) { /* punc newer than wm. */
     			if (it->getPuncSafe()->min_ts > t) { /* punc newer than wm. */
     				return nullptr;
     			}

     			/* now we see a pending punc <=wm. */

     			*has_pending_punc = true;

     			/* grab from the oldest container w/ canceled punc */
     			if (itt != right_unordered_containers.rend()) {
     				//ret = itt->getBundleUnsafe();
     				ret = itt->getBundleSafe();
     				xzl_assert(ret);
     				this->DecBundleCounter();
     				return ret;
     			}

     			/* non-empty container. grab & go */
     			//if ((ret = it->getBundleUnsafe())) {
     			if ((ret = it->getBundleSafe())) {
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
     //					E("bug: punc refcnt is %ld", refcnt);
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
     				xzl_assert(it == right_unordered_containers.rbegin());
     //				it = containers_.erase(it);   /* erase w/ reverse iterator */
     				it ++;
     				//XXX Should NOT destroy containers here
     				//xzl_assert(false && "should not return Punc here");
     				/*
     				it = list<bundle_container>::reverse_iterator(left_unordered_containers.erase(it.base()));
     				std::cout << __FILE__ << __LINE__ << "containers_.erase " << this->getName()  << std::endl;
     				*/
     				continue;
     			} else if (punc_refcnt == PUNCREF_RETRIEVED) {
     				/* punc outstanding. it may become 0 soon but we don't wait. */
     				it ++;
     				continue;
     			} else if (punc_refcnt == PUNCREF_ASSIGNED) { /* punc not emitted */
     				if (it == right_unordered_containers.rbegin()) { /* oldest container. okay to emit */
     					//TODO: we should not return punc here
     					it++;
     					continue;
     				#if 0
     					long expected = PUNCREF_ASSIGNED;
     					if (!it->punc_refcnt.compare_exchange_strong(expected,
     							PUNCREF_RETRIEVED)) {
     						bug("bad punc refcnt");
     					}
     //					auto r = it->punc_refcnt.fetch_sub(1);
     //					xzl_assert(r == 2); r = r;

     					/* XXX: opt: we may check the next newer container and
     					 * coalesce punc as needed. benefit unclear though.
     					 */
     					//XXX should not return punc here
     					xzl_assert(false && "should not return Punc here");
     					return it->punc;
     				#endif
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
     				//E("punc refcnt: %ld", punc_refcnt);
     				xzl_assert(false && "illegal punc_perfcnt?");
     			}

     			xzl_assert(false && "bug?");
     		} // while

     //		dump_containers("no bundles");

     		return nullptr;
       }

       shared_ptr<BundleBase> getOneBundleOlderThan_4(ptime const & t,
       	bool* has_pending_punc, int node = -1)
       {
       	/* hym TODO:
     	 * 	-- grab bundles/puncs from orderd_containers first
     	 *	-- if not, than only get bundles from left/right_unordered_contaienrs
     	 *	   according to left_wm and right_wm
     	 */
     	shared_ptr<BundleBase> ret;
     	//1. from orderd_containers
     	ret = getOneBundleOlderThan_ordered_4(t, has_pending_punc, node);
     	if(ret != nullptr){
     		return ret;
     	}else if(this->left_wm < this->right_wm){
     		//get bundls from left
     		ret = getOneBundleOlderThan_left_unordered_4(t, has_pending_punc, node);
     		if(ret != nullptr){
     			return ret;
     		}else{  //get bundle from right
     			ret = getOneBundleOlderThan_right_unordered_4(t, has_pending_punc, node);
     			return ret;
     		}
     	}else{
     		//get bundls from right
     		ret = getOneBundleOlderThan_right_unordered_4(t, has_pending_punc, node);
     		if(ret != nullptr){
     			return ret;
     		}else{  //get bundle from left
     			ret = getOneBundleOlderThan_left_unordered_4(t, has_pending_punc, node);
     			return ret;
     		}
     	}
       }

       /////////////////////////////////////////////////////////////////////////////////////////////
       shared_ptr<BundleBase> getOneBundle_ordered_4(int node = -1) {
       	//TODO
     	//get bundles and puncs from ordered_containers

     	//std::cout << __FILE__ << ": " <<  __LINE__ << " this trans is: " << this->getName() << std::endl;
       	shared_ptr<BundleBase> ret;

       	//unique_lock<mutex> conlock(mtx_ordered_containers);

     	//XXX rw
     	//unique_lock<mutex> conlock(mtx_container);
     	boost::shared_lock<boost::shared_mutex> rlock(mtx_container);
     	boost::upgrade_lock<boost::shared_mutex> ulock(mtx_container, boost::defer_lock);
     	unique_ptr<boost::upgrade_to_unique_lock<boost::shared_mutex>> pwlock
     								= nullptr; /* will std::move later */
     start:

       	if (ordered_containers.empty()) {
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

       	auto it = ordered_containers.rbegin();

       	while (it != ordered_containers.rend()) {

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
     				xzl_assert(it->punc_refcnt != PUNCREF_CONSUMED
     							&& it->punc_refcnt != PUNCREF_RETRIEVED);
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
     				/* has to be the oldest container since we can't have more
     				 * than one outstanding punc (the second oldest container, whose
     				 * punc is assigned, must still holds its punc).
     				 */

     				//XXX rw
     				if (!upgrade_locks(&rlock, &ulock, &pwlock)){
     					goto start;
     				}

     				xzl_assert(it == ordered_containers.rbegin());
     //				it = containers_.erase(it);   /* erase w/ reverse iterator */
     				it ++;
     				it = list<bundle_container>::reverse_iterator(ordered_containers.erase(it.base()));
     #ifdef DEBUG
     				std::cout << __FILE__ << __LINE__ << "containers_.erase " << this->getName() << std::endl;
     #endif
     				continue;
     			} else if (punc_refcnt == PUNCREF_CANCELED) {
     				/* erase */

     				//XXX rw
     				if (!upgrade_locks(&rlock, &ulock, &pwlock)){
     					goto start;
     				}

     				it ++;
     				it = list<bundle_container>::reverse_iterator(ordered_containers.erase(it.base()));
     #ifdef DEBUG
     				std::cout << __FILE__ << __LINE__ << "containers_.erase" << std::endl;
     #endif
      			} else if (punc_refcnt == PUNCREF_RETRIEVED) {
     				/* punc outstanding. it may become CONSUMED soon but we don't wait. */
     				it ++;
     				continue;
     			} else if (punc_refcnt == PUNCREF_ASSIGNED) { /* punc not emitted */
     				if (it == ordered_containers.rbegin()) {
     					/* oldest container. (older dead containers should have been
     					 * cleaned up now). okay to emit */
     //					auto r = it->punc_refcnt.fetch_sub(1);
     //					xzl_assert(r == 2); r = r;

     					long expected = PUNCREF_ASSIGNED;
     					if (!it->punc_refcnt.compare_exchange_strong(expected,
     							PUNCREF_RETRIEVED)) {
     						//bug("bad punc refcnt");
     						it ++;
     						continue;
     					}
     					/* XXX: opt: we may check the next newer container and
     					 * coalesce punc as needed. benefit unclear though.
     					 */
     					//return it->punc;
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
     				//E("punc refcnt: %ld", punc_refcnt);
     				xzl_assert(false && "illegal punc_perfcnt?");
     			}

     	  	xzl_assert(false && "bug?");

       	} // while

     //  	dump_containers("no bundles");

       	return nullptr;
       }

       shared_ptr<BundleBase> getOneBundle_left_unordered_4(int node = -1) {
       	//TODO
       	//only get bundles from left_unordered_containers
     	//XXX DON'T get punc here

     	//std::cout << __FILE__ << ": " <<  __LINE__ << " this trans is: " << this->getName() << std::endl;
       	shared_ptr<BundleBase> ret;

       	//unique_lock<mutex> conlock(left_mtx_unordered_containers);

     	//XXX rw
     	//unique_lock<mutex> conlock(mtx_container);

     	boost::shared_lock<boost::shared_mutex> rlock(mtx_container);
     	boost::upgrade_lock<boost::shared_mutex> ulock(mtx_container, boost::defer_lock);
     	unique_ptr<boost::upgrade_to_unique_lock<boost::shared_mutex>> pwlock
     								= nullptr; /* will std::move later */
     //start:
       	if (left_unordered_containers.empty()) {
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

       	auto it = left_unordered_containers.rbegin();

       	while (it != left_unordered_containers.rend()) {

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
     				xzl_assert(it->punc_refcnt != PUNCREF_CONSUMED
     							&& it->punc_refcnt != PUNCREF_RETRIEVED);
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
     				/* has to be the oldest container since we can't have more
     				 * than one outstanding punc (the second oldest container, whose
     				 * punc is assigned, must still holds its punc).
     				 */

     				//XXX rw
     			/*	if (!upgrade_locks(&rlock, &ulock, &pwlock)){
     					goto start;
     				}
     			*/
     				xzl_assert(it == left_unordered_containers.rbegin());
     //				it = containers_.erase(it);   /* erase w/ reverse iterator */
     				it ++;
     				/* XXX should not erase any containers here
     				it = list<bundle_container>::reverse_iterator(ordered_containers.erase(it.base()));
     				std::cout << __FILE__ << __LINE__ << "containers_.erase " << this->getName() << std::endl;
     				*/
     				continue;
     			} else if (punc_refcnt == PUNCREF_CANCELED) {
     				/* erase */
     				it ++;
     				/* XXX should not erase any containers here
     				it = list<bundle_container>::reverse_iterator(ordered_containers.erase(it.base()));
     				std::cout << __FILE__ << __LINE__ << "containers_.erase" << std::endl;
     				*/
     				continue;
       		} else if (punc_refcnt == PUNCREF_RETRIEVED) {
     				/* punc outstanding. it may become CONSUMED soon but we don't wait. */
     				it ++;
     				continue;
     			} else if (punc_refcnt == PUNCREF_ASSIGNED) { /* punc not emitted */
     				if (it == left_unordered_containers.rbegin()) {
     					/* oldest container. (older dead containers should have been
     					 * cleaned up now). okay to emit */
     //					auto r = it->punc_refcnt.fetch_sub(1);
     //					xzl_assert(r == 2); r = r;
     				#if 0
     					long expected = PUNCREF_ASSIGNED;
     					if (!it->punc_refcnt.compare_exchange_strong(expected,
     							PUNCREF_RETRIEVED)) {
     						bug("bad punc refcnt");
     					}
     					/* XXX: opt: we may check the next newer container and
     					 * coalesce punc as needed. benefit unclear though.
     					 */
     					return it->punc;
     				#endif
     					//XXX should not return punc here
     					it ++;
     					continue;
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
     				//E("punc refcnt: %ld", punc_refcnt);
     				xzl_assert(false && "illegal punc_perfcnt?");
     			}

     	  	xzl_assert(false && "bug?");

       	} // while

     //  	dump_containers("no bundles");

       	return nullptr;

       }

       shared_ptr<BundleBase> getOneBundle_right_unordered_4(int node = -1) {
       	//TODO
     	//XXX DON'T get punc here

     	//std::cout << __FILE__ << ": " <<  __LINE__ << " this trans is: " << this->getName() << std::endl;
       	shared_ptr<BundleBase> ret;

       	//unique_lock<mutex> conlock(right_mtx_unordered_containers);

     	//XXX rw
     	//unique_lock<mutex> conlock(mtx_container);
     	boost::shared_lock<boost::shared_mutex> rlock(mtx_container);
     	boost::upgrade_lock<boost::shared_mutex> ulock(mtx_container, boost::defer_lock);
     	unique_ptr<boost::upgrade_to_unique_lock<boost::shared_mutex>> pwlock
     								= nullptr; /* will std::move later */
     //start:

       	if (right_unordered_containers.empty()) {
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

       	auto it = right_unordered_containers.rbegin();

       	while (it != right_unordered_containers.rend()) {

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
     				xzl_assert(it->punc_refcnt != PUNCREF_CONSUMED
     							&& it->punc_refcnt != PUNCREF_RETRIEVED);
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
     				/* has to be the oldest container since we can't have more
     				 * than one outstanding punc (the second oldest container, whose
     				 * punc is assigned, must still holds its punc).
     				 */
     				xzl_assert(it == right_unordered_containers.rbegin());
     //				it = containers_.erase(it);   /* erase w/ reverse iterator */
     				it ++;
     				/* XXX should not erase any containers here
     				it = list<bundle_container>::reverse_iterator(ordered_containers.erase(it.base()));
     				std::cout << __FILE__ << __LINE__ << "containers_.erase " << this->getName() << std::endl;
     				*/
     				continue;
     			} else if (punc_refcnt == PUNCREF_CANCELED) {
     				/* erase */
     				it ++;
     				/* XXX should not erase any containers here
     				it = list<bundle_container>::reverse_iterator(ordered_containers.erase(it.base()));
     				std::cout << __FILE__ << __LINE__ << "containers_.erase" << std::endl;
     				*/
     				continue;
       		} else if (punc_refcnt == PUNCREF_RETRIEVED) {
     				/* punc outstanding. it may become CONSUMED soon but we don't wait. */
     				it ++;
     				continue;
     			} else if (punc_refcnt == PUNCREF_ASSIGNED) { /* punc not emitted */
     				if (it == right_unordered_containers.rbegin()) {
     					/* oldest container. (older dead containers should have been
     					 * cleaned up now). okay to emit */
     //					auto r = it->punc_refcnt.fetch_sub(1);
     //					xzl_assert(r == 2); r = r;
     				#if 0
     					long expected = PUNCREF_ASSIGNED;
     					if (!it->punc_refcnt.compare_exchange_strong(expected,
     							PUNCREF_RETRIEVED)) {
     						bug("bad punc refcnt");
     					}
     					/* XXX: opt: we may check the next newer container and
     					 * coalesce punc as needed. benefit unclear though.
     					 */
     					return it->punc;
     				#endif
     					//XXX should not return punc here
     					it ++;
     					continue;
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
     				//E("punc refcnt: %ld", punc_refcnt);
     				xzl_assert(false && "illegal punc_perfcnt?");
     			}

     	  	xzl_assert(false && "bug?");

       	} // while

     //  	dump_containers("no bundles");

       	return nullptr;
       }

       shared_ptr<BundleBase> getOneBundle_4(int node = -1) {

        	/* hym TODO:
      	 * 	-- grab bundles/puncs from orderd_containers first
      	 *	-- if not, than only get bundles from left/right_unordered_contaienrs
      	 *	   according to left_wm and right_wm
      	 */
      	shared_ptr<BundleBase> ret;
      	//1. from orderd_containers
      	ret = getOneBundle_ordered_4(node);
      	if(ret != nullptr){
      		return ret;
      	}else if(this->left_wm < this->right_wm){
      		//get bundls from left
      		ret = getOneBundle_left_unordered_4(node);
      		if(ret != nullptr){
      			return ret;
      		}else{  //get bundle from right
      			ret = getOneBundle_right_unordered_4(node);
      			return ret;
      		}
      	}else{
      		//get bundls from right
      		ret = getOneBundle_right_unordered_4(node);
      		if(ret != nullptr){
      			return ret;
      		}else{  //get bundle from left
      			ret = getOneBundle_left_unordered_4(node);
      			return ret;
      		}
      	}
        }

        //hym TODO: deposit bundles/puncs to type 5 ..

       //type4 trans opens a new container in the unordered_containers of type 5
       //here:  downt is type 5 trans
       // upcon is 4's left_unordered_containers, downt is type 5 trans
       int openDownstreamContainer_unordered_4_left(bundle_container * const upcon,
       	PTransform *downt){
       	// TODO
      	// be care for about assigning side_info to new containers
      	// inherent side_info from upcon: should be 1 or 2

        	xzl_assert(downt && upcon);
      	//std::cout << __FILE__ << ": " <<  __LINE__ << std::endl;
      	//std::cout << "this trans is: " << this->getName() << std::endl;
      	//std::cout << "downt is: " << downt->getName() << std::endl;
      //  	/* fast path w/o locking */
        	if (upcon->downstream)
        		return 0;

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
        //HYM	unique_lock<mutex> conlock(left_mtx_unordered_containers);

      	//here this is SimpleMapper
      	// downt is Join
        	//unique_lock<mutex> down_conlock(downt->mtx_container);
        //HYM	unique_lock<mutex> down_conlock(downt->mtx_unordered_containers);

      /*     these locks have been hold in deposit, should not get these locks again
      	unique_lock<mutex> conlock(mtx_container);
        	unique_lock<mutex> down_conlock(downt->mtx_container);
      */
  	
	/* reader lock for *this* trans */
	boost::shared_lock<boost::shared_mutex> conlock(mtx_container);
	/* writer lock for downstream trans */
	boost::unique_lock<boost::shared_mutex> down_conlock(downt->mtx_container);
        
	if (upcon->downstream) /* check again */
        		return 0;
      	//std::cout << __FILE__ << ": " <<  __LINE__ << std::endl;
        	int cnt = 0;
        	bool miss = false; (void)sizeof(miss);
      	// miss = miss;
      	//std::cout << __FILE__ << ": " <<  __LINE__ << " this trans is: " << this->getName() << " downt is " << downt->getName() <<std::endl;
        	for (auto it = this->left_unordered_containers.rbegin();
        			it != this->left_unordered_containers.rend(); it ++) {

      		//std::cout << __FILE__ << ": " <<  __LINE__ << std::endl;
        		/* also lock all downstream containers. good idea? */
      //  		unique_lock<mutex> down_conlock(downt->mtx_container);

        		if (!it->downstream) {
        			miss = true;
      			//std::cout << __FILE__ << ": " <<  __LINE__ << std::endl;
        			/* alloc & link downstream container. NB: it->downstream is an atomic
        			 * pointer, so it implies a fence here. */
        			downt->unordered_containers.emplace_front();
        			it->downstream = &(downt->unordered_containers.front());

      			//hym: XXX
      			//assigne the side_info to this new container
      			//it->downstream->side_info = this->get_side_info(); //0->containrs_, 1->left_containers_in, 2->right_containers_in
        			//it->downstream.load()->set_side_info(this->get_side_info());
        			it->downstream.load()->set_side_info(it->get_side_info());
      			xzl_assert(it->downstream.load()->get_side_info() == 1);	//left to unordere_containers

      			cnt ++;
        		} else { /* downstream link exists */
        			assert (!miss && "bug? an older container misses downstream link.");
      #ifdef DEBUG
        			bool found = false;
        			for (auto & con : downt->unordered_containers) {
        				if (&con == it->downstream) {
        					found = true;
        					break;
        				}
        			}
        			/* NB: downstream container, once received & processed the punc
        			 * from the upstream, may be cleaned up prior to the upstream
        			 * container. In this case, a dangling downstream pointer is possible.
        			 */
        			if (!found && it->punc_refcnt != PUNCREF_CONSUMED) {
        				/*
      				E("%s: downstream (%s) container does not exist",
        						this->name.c_str(), downt->name.c_str());
      					cout << " -------------------------------\n";
      					dump_containers("bug");
      					downt->dump_containers("bug");
      					cout << " -------------------------------\n\n";
        				abort();
      				*/
      				xzl_assert(false && "downstream container does not exist");
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


       //type4 trans opens a new container in the unordered_containers of type 5
       //here:  downt is type 5 trans
       // upcon is 4's right_unordered_containers, downt is type 5 trans
       int openDownstreamContainer_unordered_4_right(bundle_container * const upcon,
       	PTransform *downt){
       	// TODO
      	// be care for about assigning side_info to new containers
      	// inherent side_info from upcon: should be 1 or 2

        	xzl_assert(downt && upcon);
      	//std::cout << __FILE__ << ": " <<  __LINE__ << std::endl;
      	//std::cout << "this trans is: " << this->getName() << std::endl;
      	//std::cout << "downt is: " << downt->getName() << std::endl;
      //  	/* fast path w/o locking */
        	if (upcon->downstream)
        		return 0;

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
      //HYM  	unique_lock<mutex> conlock(right_mtx_unordered_containers);

      	//here this is SimpleMapper
      	// downt is Join
        	//unique_lock<mutex> down_conlock(downt->mtx_container);
        //HYM	unique_lock<mutex> down_conlock(downt->mtx_unordered_containers);

      /*	these locks have been hodl in deposit, should not get them again
      	unique_lock<mutex> conlock(mtx_container);
      	unique_lock<mutex> down_conlock(downt->mtx_container);
      */
  	/* reader lock for *this* trans */
	boost::shared_lock<boost::shared_mutex> conlock(mtx_container);
	/* writer lock for downstream trans */
	boost::unique_lock<boost::shared_mutex> down_conlock(downt->mtx_container);
        
	if (upcon->downstream) /* check again */
        		return 0;
      	//std::cout << __FILE__ << ": " <<  __LINE__ << std::endl;
        	int cnt = 0;
        	bool miss = false; (void)sizeof(miss);
      	// miss = miss;
      	//std::cout << __FILE__ << ": " <<  __LINE__ << " this trans is: " << this->getName() << " downt is " << downt->getName() <<std::endl;
        	for (auto it = this->right_unordered_containers.rbegin();
        			it != this->right_unordered_containers.rend(); it ++) {

      		//std::cout << __FILE__ << ": " <<  __LINE__ << std::endl;
        		/* also lock all downstream containers. good idea? */
      //  		unique_lock<mutex> down_conlock(downt->mtx_container);

        		if (!it->downstream) {
        			miss = true;
      			//std::cout << __FILE__ << ": " <<  __LINE__ << std::endl;
        			/* alloc & link downstream container. NB: it->downstream is an atomic
        			 * pointer, so it implies a fence here. */
        			downt->unordered_containers.emplace_front();
        			it->downstream = &(downt->unordered_containers.front());

      			//hym: XXX
      			//assigne the side_info to this new container
      			//it->downstream->side_info = this->get_side_info(); //0->containrs_, 1->left_containers_in, 2->right_containers_in
        			//it->downstream.load()->set_side_info(this->get_side_info());
        			it->downstream.load()->set_side_info(it->get_side_info());
      			xzl_assert(it->downstream.load()->get_side_info() == 2);	//right to unordere_containers

      			cnt ++;
        		} else { /* downstream link exists */
        			assert (!miss && "bug? an older container misses downstream link.");
      #ifdef DEBUG
        			bool found = false;
        			for (auto & con : downt->unordered_containers) {
        				if (&con == it->downstream) {
        					found = true;
        					break;
        				}
        			}
        			/* NB: downstream container, once received & processed the punc
        			 * from the upstream, may be cleaned up prior to the upstream
        			 * container. In this case, a dangling downstream pointer is possible.
        			 */
        			if (!found && it->punc_refcnt != PUNCREF_CONSUMED) {
        				/*E("%s: downstream (%s) container does not exist",
        						this->name.c_str(), downt->name.c_str());
      					cout << " -------------------------------\n";
      					dump_containers("bug");
      					downt->dump_containers("bug");
      					cout << " -------------------------------\n\n";
        				abort();
      				*/
      				xzl_assert(false && "downstream container does not exist");
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


       //type4 trans opens a new container in the ordered_containers of type 5
       //here:  downt is type 5 trans
       int openDownstreamContainer_ordered_4(bundle_container * const upcon,
       	PTransform *downt){
       	//TODO
      	// be care for about assigning side_info to new containers
      	// set side_info to 0, all containers in oerdered_containers should have side_info 0

        	xzl_assert(downt && upcon);
      	//std::cout << __FILE__ << ": " <<  __LINE__ << std::endl;
      	//std::cout << "this trans is: " << this->getName() << std::endl;
      	//std::cout << "downt is: " << downt->getName() << std::endl;
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
      //HYM  	unique_lock<mutex> conlock(mtx_ordered_containers);

        	//unique_lock<mutex> down_conlock(downt->mtx_container);
       //HYM 	unique_lock<mutex> down_conlock(downt->mtx_ordered_containers);

      /*	these locks have been hold in deposit, should not get them again
      	unique_lock<mutex> conlock(mtx_container);
       	unique_lock<mutex> down_conlock(downt->mtx_container);
      */
  	/* reader lock for *this* trans */
	boost::shared_lock<boost::shared_mutex> conlock(mtx_container);
	/* writer lock for downstream trans */
	boost::unique_lock<boost::shared_mutex> down_conlock(downt->mtx_container);
        	
	if (upcon->downstream) /* check again */
        		return 0;
      	//std::cout << __FILE__ << ": " <<  __LINE__ << std::endl;
        	int cnt = 0;
        	bool miss = false; (void)sizeof(miss);
      	// miss = miss;
      	//std::cout << __FILE__ << ": " <<  __LINE__ << " this trans is: " << this->getName() << " downt is " << downt->getName() <<std::endl;
        	for (auto it = this->ordered_containers.rbegin();
        			it != this->ordered_containers.rend(); it ++) {

      		//std::cout << __FILE__ << ": " <<  __LINE__ << std::endl;
        		/* also lock all downstream containers. good idea? */
      //  		unique_lock<mutex> down_conlock(downt->mtx_container);

        		if (!it->downstream) {
        			miss = true;
      			std::cout << __FILE__ << ": " <<  __LINE__ << std::endl;
        			/* alloc & link downstream container. NB: it->downstream is an atomic
        			 * pointer, so it implies a fence here. */
        			downt->ordered_containers.emplace_front();
        			it->downstream = &(downt->ordered_containers.front());

      			//hym: XXX
      			//assigne the side_info to this new container
      			//it->downstream->side_info = this->get_side_info(); //0->containrs_, 1->left_containers_in, 2->right_containers_in
        			//it->downstream.load()->set_side_info(this->get_side_info());
        			it->downstream.load()->set_side_info(it->get_side_info());
      			xzl_assert(it->downstream.load()->get_side_info() == 0);	//XXX containers in ordered_containers should have side_info 0

      			cnt ++;
        		} else { /* downstream link exists */
      			//XXX comment this??XXX
        //			assert (!miss && "bug? an older container misses downstream link.");
      #ifdef DEBUG
        			bool found = false;
        			for (auto & con : downt->ordered_containers) {
        				if (&con == it->downstream) {
        					found = true;
        					break;
        				}
        			}
        			/* NB: downstream container, once received & processed the punc
        			 * from the upstream, may be cleaned up prior to the upstream
        			 * container. In this case, a dangling downstream pointer is possible.
        			 */
        			if (!found && it->punc_refcnt != PUNCREF_CONSUMED) {
        				/*E("%s: downstream (%s) container does not exist",
        						this->name.c_str(), downt->name.c_str());
      					cout << " -------------------------------\n";
      					dump_containers("bug");
      					downt->dump_containers("bug");
      					cout << " -------------------------------\n\n";
        				abort();
      				*/

      				xzl_assert(false && "downstream container does not exist");
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

       //type 4 trans deposit a  bundle to downstream(type 5 trans' unordered_contaienrs or ordered_containers)
       void  depositBundleDownstream_4(PTransform *downt,
       	shared_ptr<BundleBase> input_bundle,
      	vector<shared_ptr<BundleBase>> const & output_bundles, int node = -1){
      	//TODO

      	xzl_assert(downt->get_side_info() == 5);

        	bool try_open = false; /* have we tried open downstream containers? */
      	//std::cout << __FILE__ << ": " <<  __LINE__ << std::endl;
        	//std::cout << "trans is: " << downt->getName() << std::endl;
      	// try_open = try_open;
      	xzl_assert(downt);
        	bundle_container *upcon = localBundleToContainer(input_bundle);
        	xzl_assert(upcon); /* enclosing container */

      	//XXX rw
      /*
       	unique_lock<mutex> conlock(mtx_container);
      	unique_lock<mutex> down_conlock(downt->mtx_container);
      */
#if 0
      	/* reader lock for *this* trans */
      	boost::shared_lock<boost::shared_mutex> conlock(mtx_container);
      	/* writer lock for downstream trans */
      	boost::unique_lock<boost::shared_mutex> down_conlock(downt->mtx_container);
#endif
      deposit:
         {		//{
      // 	unique_lock<mutex> conlock(mtx_container);
      //	unique_lock<mutex> down_conlock(downt->mtx_container);
      			if (upcon->downstream) {
      				/* fast path. NB @downstream is atomic w/ seq consistency so it's
      				 * safe. */
      				downt->depositBundlesToContainer(upcon->downstream, output_bundles, node);
      				downt->IncBundleCounter(output_bundles.size());
      				//std::cout << __FILE__ << ": " << __LINE__ << "join deposits bundle to its out: " << upcon->get_side_info() << std::endl;
      				return;
      			}
      		//}
      #if 0
      	if(upcon->get_side_info() == 1 || upcon->get_side_info() == 2){
      		unique_lock<mutex> conlock(downt->mtx_unordered_containers);
      		if (upcon->downstream) {
      			downt->depositBundlesToContainer(upcon->downstream, output_bundles, node);
      			downt->IncBundleCounter(output_bundles.size());
      			return;
      		}

      	}else if(upcon->get_side_info() == 0){ //from 4's ordered_containers
      		unique_lock<mutex> conlock(downt->mtx_ordered_containers);
      		if (upcon->downstream) {
      			downt->depositBundlesToContainer(upcon->downstream, output_bundles, node);
      			downt->IncBundleCounter(output_bundles.size());
      			return;
      		}
      	}else{
      		xzl_assert(false && "Unknown side_info");
      	}
      #endif

      	// try_open = try_open;
        	xzl_assert(!try_open && "bug: already tried to open");

        	/* it's possible somebody slipped in & opened the downstream container
        	 * already. So ret == 0 then.
        	 */
      /*
      //  	int ret =
        			openDownstreamContainers(upcon, downt);
      //  	xzl_assert(ret && "no downstream container opened?");
      */

      /*
      	if(this->get_side_info() == 0 || this->get_side_info() == 3){ //defaut transforms and join
      		openDownstreamContainers(upcon, downt);
      	}else if(this->get_side_info() == 1){ //Simplemapper 1, left
      		openDownstreamContainers_left(upcon, downt); //downt should be join, open a container in join's left_containers_in
      	}else if(this->get_side_info() == 2){ //Simplemapper 2, right
      		openDownstreamContainers_right(upcon, downt); //downt should be join, open a container in join's right_containers_in
      	}else{
      		xzl_assert("Bug: Unknown know side info");
      	}
      */
      	int flag =100; (void)sizeof(flag); // flag = flag;
       //{
      /*
      	unique_lock<mutex> conlock1(left_mtx_unordered_containers);
       	unique_lock<mutex> conlock2(right_mtx_unordered_containers);
       	unique_lock<mutex> conlock5(mtx_ordered_containers);
       	unique_lock<mutex> conlock3(downt->mtx_unordered_containers);
       	unique_lock<mutex> conlock4(downt->mtx_ordered_containers);
      */
       	//unique_lock<mutex> conlock(mtx_container);
       	//unique_lock<mutex> conlock(downt->mtx_container);

      	if(upcon->get_side_info() == 1){
      		//std::cout << __FILE__ << ": " <<  __LINE__ << std::endl;
      		//openJoinOutContainer_left(upcon, downt);
      		openDownstreamContainer_unordered_4_left(upcon, downt);
      		//std::cout << __FILE__ << ": " <<  __LINE__ << std::endl;
      		flag = 1;
      	}else if(upcon->get_side_info() == 2){
      		//std::cout << __FILE__ << ": " <<  __LINE__ << std::endl;
      		//openJoinOutContainer_right(upcon, downt);
      		openDownstreamContainer_unordered_4_right(upcon, downt);
      		//std::cout << __FILE__ << ": " <<  __LINE__ << std::endl;
      		flag = 2;
      	}else if(upcon->get_side_info() == 0){ //from 4's ordered_containers
      		//std::cout << "before open===============" << std::endl;
      		openDownstreamContainer_ordered_4(upcon, downt);
      		//std::cout << "after open=================" << std::endl;
      		flag = 0;
      	}else {
      		xzl_assert(false && "Join: wrong side info");
      	}
         }//end lock
         /*	if(!upcon->downstream){
      		std::cout << "side info is " << upcon->get_side_info() << " flag is " << flag <<  std::endl;
      	}
         */
      	xzl_assert(upcon->downstream);
       //}
        	try_open = true;
        	goto deposit;
       }

       //type 4 trans deposit a punc to downstream(type 5's trans ordered_containers)
       //becase type 4 cannot get punc from left/right_unordered_containers, it can
       //only get punc from its ordered_containers, so the punc should be deposited to
       // type 5 trans' ordered_containers
       void depositPuncDownstream_4(PTransform* downt,
       	shared_ptr<Punc> const & input_punc, shared_ptr<Punc> punc, int node = -1){
       	//TODO

        	bool try_open = false;

        	xzl_assert(downt);
        	bundle_container *upcon = localBundleToContainer(input_punc);
        	xzl_assert(upcon); /* enclosing container */

      	//XXX rw
      /*	unique_lock<mutex> conlock(mtx_container);
      	unique_lock<mutex> down_conlock(downt->mtx_container);
      */
#if 0
      	/* reader lock for *this* trans */
      	boost::shared_lock<boost::shared_mutex> conlock(mtx_container);
      	/* writer lock for downstream trans */
      	boost::unique_lock<boost::shared_mutex> down_conlock(downt->mtx_container);
#endif
      deposit:
           {		//{
      //	unique_lock<mutex> conlock(mtx_container);
      //	unique_lock<mutex> down_conlock(downt->mtx_container);
      			if (upcon->downstream) { /* fast path. NB @downstream is atomic, implying fence */
      				downt->depositOnePuncToContainer(upcon->downstream, punc, node);
      				//upcon->has_punc = 1; //the upcon has a container already
      				upcon->downstream.load()->has_punc = 1;
      				//std::cout << __FILE__ << ": " <<  __LINE__ << " " << this->getName() << " deposit a punc to " << downt->getName() << " to side: " << upcon->get_side_info()  << std::endl;
      				return;
      			}
      		//}
      #if 0
      		if(upcon->get_side_info() == 1 || upcon->get_side_info() == 2){
      			unique_lock<mutex> conlock(downt->mtx_unordered_containers);
      			if (upcon->downstream) {
      				downt->depositOnePuncToContainer(upcon->downstream, punc, node);
      				upcon->downstream.load()->has_punc = 1;
      				return;
      			}

      		}else if(upcon->get_side_info() == 0){ //from 4's ordered_containers
      			unique_lock<mutex> conlock(downt->mtx_ordered_containers);
      			if (upcon->downstream) {
      				downt->depositOnePuncToContainer(upcon->downstream, punc, node);
      				upcon->downstream.load()->has_punc = 1;
      				return;
      			}
      		}else{
      			xzl_assert(false && "Unknown side_info");
      		}
      #endif

      		// try_open = try_open;
      		xzl_assert(!try_open && "bug: already tried to open");

      		//openDownstreamContainers(upcon, downt);
      		//hym:XXX
      	/*	if(this->get_side_info() == 0 || this->get_side_info() == 3){ //defaut transforms and join
      			openDownstreamContainers(upcon, downt);
      		}else if(this->get_side_info() == 1){ //Simplemapper 1, left
      			openDownstreamContainers_left(upcon, downt); //downt should be join, open a container in join's left_containers_in
      		}else if(this->get_side_info() == 2){ //Simplemapper 2, right
      			openDownstreamContainers_right(upcon, downt); //downt should be join, open a container in join's right_containers_in
      		}else{
      			xzl_assert("Bug: Unknown know side info");
      		}
      	*/

      		if(upcon->get_side_info() == 1){
      			//openJoinOutContainer_left(upcon, downt);
      			//openJoinDownstreamContainer_left(upcon, downt);
      			openDownstreamContainer_unordered_4_left(upcon, downt);
      		}else if(upcon->get_side_info() == 2){
      			//openJoinOutContainer_right(upcon, downt);
      			//openJoinDownstreamContainer_right(upcon, downt);
      			openDownstreamContainer_unordered_4_right(upcon, downt);
      		}else if(upcon->get_side_info() == 0){
      			openDownstreamContainer_ordered_4(upcon, downt);
      		}else{
      			xzl_assert(false && "Join: wrong side info");
      		}
      	}//end lock
      		xzl_assert(upcon->downstream);
      		try_open = true;
      		goto deposit;

       }

       /***************************************************************
                   hym: end of type 4 trans
       ****************************************************************/

       /***************************************************************
                   hym: type 5 trans
       ****************************************************************/
       //TODO
       shared_ptr<BundleBase> getOneBundleOlderThan_ordered_5(ptime const & t,
       	bool* has_pending_punc, int node = -1){
       	//TODO

      	*has_pending_punc = false;

      	//std::cout << __FILE__ << ": " <<  __LINE__ << " This trans is: " << this->getName()  << std::endl;

        	/* fast path: if local wm more recent (>t), no work.
        	 * (and any unemitted punc must >t) */
        	if (this->GetWatermark() > t) {
        	  return nullptr;
        	}

        	/* slowpath: there may be some work in this trans. lock containers */

        	shared_ptr<BundleBase> ret;
      		//unique_lock<mutex> conlock(mtx_ordered_containers);
      		//XXX rw
      		//unique_lock<mutex> conlock(mtx_container);

      		boost::shared_lock<boost::shared_mutex> rlock(mtx_container);
      		boost::upgrade_lock<boost::shared_mutex> ulock(mtx_container, boost::defer_lock);
      		unique_ptr<boost::upgrade_to_unique_lock<boost::shared_mutex>> pwlock
      								= nullptr; /* will std::move later */
      start:
      		if (ordered_containers.empty()) {
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

      		auto it = ordered_containers.rbegin();
      		/* to remember the 1st container that 1) has invalidated punc; 2) has
      		 * some bundles
      		 */
      		auto itt = ordered_containers.rend();

      		while (it != ordered_containers.rend()) {

      			if (it->punc_refcnt == PUNCREF_UNDECIDED) { /* we don't deal open container */
      				xzl_assert(it->verifyPuncSafe() && "bug?");
      				return nullptr;
      			}

      			/* punc has been canceled, if we see a later container assigned with
      			 * a punc, we may return a bundle from *this* container.
      			 */
      			if (it->punc_refcnt == PUNCREF_CANCELED) {
      				if (it->refcnt == 0) { /* dead container, clean it up. */
      					/* because punc only canceled after bundles are deposited.*/

      					//XXX rw
      					if (!upgrade_locks(&rlock, &ulock, &pwlock)){
      						goto start;
      					}

      					xzl_assert(it->bundles.empty());
      					it ++;
      					it = list<bundle_container>::reverse_iterator(ordered_containers.erase(it.base()));
      					std::cout << __FILE__ << __LINE__ << "containers_.erase" << std::endl;
      				} else if (!it->bundles.empty() && itt == ordered_containers.rend()) {
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

      			//if (it->punc->min_ts > t) { /* punc newer than wm. */
      			if (it->getPuncSafe()->min_ts > t) { /* punc newer than wm. */
      				return nullptr;
      			}

      			/* now we see a pending punc <=wm. */

      			*has_pending_punc = true;

      			/* grab from the oldest container w/ canceled punc */
      			if (itt != ordered_containers.rend()) {
      				//ret = itt->getBundleUnsafe();
      				ret = itt->getBundleSafe();
      				xzl_assert(ret);
      				this->DecBundleCounter();
      				return ret;
      			}

      			/* non-empty container. grab & go */
      			//if ((ret = it->getBundleUnsafe())) {
      			if ((ret = it->getBundleSafe())) {
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
      //					E("bug: punc refcnt is %ld", refcnt);
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

      				//XXX rw
      				if (!upgrade_locks(&rlock, &ulock, &pwlock)){
      					goto start;
      				}

      				xzl_assert(it == ordered_containers.rbegin());
      //				it = containers_.erase(it);   /* erase w/ reverse iterator */
      				it ++;
      				it = list<bundle_container>::reverse_iterator(ordered_containers.erase(it.base()));
      #ifdef DEBUG
      				std::cout << __FILE__ << __LINE__ << "containers_.erase " << this->getName()  << std::endl;
      #endif
      				continue;
      			} else if (punc_refcnt == PUNCREF_RETRIEVED) {
      				/* punc outstanding. it may become 0 soon but we don't wait. */
      				it ++;
      				continue;
      			} else if (punc_refcnt == PUNCREF_ASSIGNED) { /* punc not emitted */
      				if (it == ordered_containers.rbegin()) { /* oldest container. okay to emit */
      					long expected = PUNCREF_ASSIGNED;
      					if (!it->punc_refcnt.compare_exchange_strong(expected,
      							PUNCREF_RETRIEVED)) {
      						//bug("bad punc refcnt");
      						it ++;
      						continue;
      					}
      //					auto r = it->punc_refcnt.fetch_sub(1);
      //					xzl_assert(r == 2); r = r;

      					/* XXX: opt: we may check the next newer container and
      					 * coalesce punc as needed. benefit unclear though.
      					 */
      					//return it->punc;
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
      				//E("punc refcnt: %ld", punc_refcnt);
      				xzl_assert(false && "illegal punc_perfcnt?");
      			}

      			xzl_assert(false && "bug?");
      		} // while

      //		dump_containers("no bundles");

      		return nullptr;
       }

       shared_ptr<BundleBase> getOneBundleOlderThan_unordered_5(ptime const & t,
       	bool* has_pending_punc, int node = -1){
       	//TODO
      	//should not get punc here
      	// ref getOneBundleOlderThan_left_unordered_4

      	*has_pending_punc = false;

      	//std::cout << __FILE__ << ": " <<  __LINE__ << " This trans is: " << this->getName()  << std::endl;

        	/* fast path: if local wm more recent (>t), no work.
        	 * (and any unemitted punc must >t) */
        	if (this->GetWatermark() > t) {
        	  return nullptr;
        	}

        	/* slowpath: there may be some work in this trans. lock containers */

        	shared_ptr<BundleBase> ret;
      		//unique_lock<mutex> conlock(mtx_unordered_containers);

      		//unique_lock<mutex> conlock(mtx_container);

      		boost::shared_lock<boost::shared_mutex> rlock(mtx_container);
      		boost::upgrade_lock<boost::shared_mutex> ulock(mtx_container, boost::defer_lock);
      		unique_ptr<boost::upgrade_to_unique_lock<boost::shared_mutex>> pwlock
      								= nullptr; /* will std::move later */
      //start:
      		if (unordered_containers.empty()) {
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

      		auto it = unordered_containers.rbegin();
      		/* to remember the 1st container that 1) has invalidated punc; 2) has
      		 * some bundles
      		 */
      		auto itt = unordered_containers.rend();

      		while (it != unordered_containers.rend()) {

      			if (it->punc_refcnt == PUNCREF_UNDECIDED) { /* we don't deal open container */
      				xzl_assert(it->verifyPuncSafe() && "bug?");
      				return nullptr;
      			}

      			/* punc has been canceled, if we see a later container assigned with
      			 * a punc, we may return a bundle from *this* container.
      			 */
      			if (it->punc_refcnt == PUNCREF_CANCELED) {
      				if (it->refcnt == 0) { /* dead container, clean it up. */
      					/* because punc only canceled after bundles are deposited.*/
      					xzl_assert(it->bundles.empty());
      					it ++;
      					//XXX should not destroy container here
      					continue;
      					xzl_assert(false && "should not return Punc here");
      					/*
      					it = list<bundle_container>::reverse_iterator(left_unordered_containers.erase(it.base()));
      					std::cout << __FILE__ << __LINE__ << "containers_.erase" << std::endl;
      					*/
      				} else if (!it->bundles.empty() && itt == unordered_containers.rend()) {
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

      			//if (it->punc->min_ts > t) { /* punc newer than wm. */
      			if (it->getPuncSafe()->min_ts > t) { /* punc newer than wm. */
      				return nullptr;
      			}

      			/* now we see a pending punc <=wm. */

      			*has_pending_punc = true;

      			/* grab from the oldest container w/ canceled punc */
      			if (itt != unordered_containers.rend()) {
      				//ret = itt->getBundleUnsafe();
      				ret = itt->getBundleSafe();
      				xzl_assert(ret);
      				this->DecBundleCounter();
      				return ret;
      			}

      			/* non-empty container. grab & go */
      			//if ((ret = it->getBundleUnsafe())) {
      			if ((ret = it->getBundleSafe())) {
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
      //					E("bug: punc refcnt is %ld", refcnt);
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
      				xzl_assert(it == unordered_containers.rbegin());
      //				it = containers_.erase(it);   /* erase w/ reverse iterator */
      				it ++;
      				//XXX Should NOT destroy containers here
      				//xzl_assert(false && "should not return Punc here");
      				/*
      				it = list<bundle_container>::reverse_iterator(left_unordered_containers.erase(it.base()));
      				std::cout << __FILE__ << __LINE__ << "containers_.erase " << this->getName()  << std::endl;
      				*/
      				continue;
      			} else if (punc_refcnt == PUNCREF_RETRIEVED) {
      				/* punc outstanding. it may become 0 soon but we don't wait. */
      				it ++;
      				continue;
      			} else if (punc_refcnt == PUNCREF_ASSIGNED) { /* punc not emitted */
      				if (it == unordered_containers.rbegin()) { /* oldest container. okay to emit */
      					//TODO: we should not return punc here
      					it++;
      					continue;
      				#if 0
      					long expected = PUNCREF_ASSIGNED;
      					if (!it->punc_refcnt.compare_exchange_strong(expected,
      							PUNCREF_RETRIEVED)) {
      						bug("bad punc refcnt");
      					}
      //					auto r = it->punc_refcnt.fetch_sub(1);
      //					xzl_assert(r == 2); r = r;

      					/* XXX: opt: we may check the next newer container and
      					 * coalesce punc as needed. benefit unclear though.
      					 */
      					//XXX should not return punc here
      					xzl_assert(false && "should not return Punc here");
      					return it->punc;
      				#endif
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
      				//E("punc refcnt: %ld", punc_refcnt);
      				xzl_assert(false && "illegal punc_perfcnt?");
      			}

      			xzl_assert(false && "bug?");
      		} // while

      //		dump_containers("no bundles");

      		return nullptr;
       }

       shared_ptr<BundleBase> getOneBundleOlderThan_5(ptime const & t,
       	bool* has_pending_punc, int node = -1){
      	//TODO
      	shared_ptr<BundleBase> ret;
      	ret = getOneBundleOlderThan_ordered_5(t, has_pending_punc, node);
      	if(ret != nullptr){
      		return ret;
      	}else{
      		ret = getOneBundleOlderThan_unordered_5(t, has_pending_punc, node);
      		return ret;
      	}

       }
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
       shared_ptr<BundleBase> getOneBundle_ordered_5(int node = -1) {
      	//TODO

      	//std::cout << __FILE__ << ": " <<  __LINE__ << " this trans is: " << this->getName() << std::endl;
        	shared_ptr<BundleBase> ret;

        	//unique_lock<mutex> conlock(mtx_ordered_containers);
        	//XXX rw
      	//unique_lock<mutex> conlock(mtx_container);

      	boost::shared_lock<boost::shared_mutex> rlock(mtx_container);
      	boost::upgrade_lock<boost::shared_mutex> ulock(mtx_container, boost::defer_lock);
      	unique_ptr<boost::upgrade_to_unique_lock<boost::shared_mutex>> pwlock
      								= nullptr; /* will std::move later */
      start:
        	if (ordered_containers.empty()) {
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

        	auto it = ordered_containers.rbegin();

        	while (it != ordered_containers.rend()) {

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
      				xzl_assert(it->punc_refcnt != PUNCREF_CONSUMED
      							&& it->punc_refcnt != PUNCREF_RETRIEVED);
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
      				/* has to be the oldest container since we can't have more
      				 * than one outstanding punc (the second oldest container, whose
      				 * punc is assigned, must still holds its punc).
      				 */

      				//XXX rw
      				if (!upgrade_locks(&rlock, &ulock, &pwlock)){
      					goto start;
      				}

      				xzl_assert(it == ordered_containers.rbegin());
      //				it = containers_.erase(it);   /* erase w/ reverse iterator */
      				it ++;
      				it = list<bundle_container>::reverse_iterator(ordered_containers.erase(it.base()));
      #ifdef DEBUG
      				std::cout << __FILE__ << __LINE__ << "containers_.erase " << this->getName() << std::endl;
      #endif
      				continue;
      			} else if (punc_refcnt == PUNCREF_CANCELED) {
      				/* erase */

      				//XXX rw
      				if (!upgrade_locks(&rlock, &ulock, &pwlock)){
      					goto start;
      				}

      				it ++;
      				it = list<bundle_container>::reverse_iterator(ordered_containers.erase(it.base()));
      #ifdef DEBUG
      				std::cout << __FILE__ << __LINE__ << "containers_.erase" << std::endl;
      #endif
       			} else if (punc_refcnt == PUNCREF_RETRIEVED) {
      				/* punc outstanding. it may become CONSUMED soon but we don't wait. */
      				it ++;
      				continue;
      			} else if (punc_refcnt == PUNCREF_ASSIGNED) { /* punc not emitted */
      				if (it == ordered_containers.rbegin()) {
      					/* oldest container. (older dead containers should have been
      					 * cleaned up now). okay to emit */
      //					auto r = it->punc_refcnt.fetch_sub(1);
      //					xzl_assert(r == 2); r = r;

      					long expected = PUNCREF_ASSIGNED;
      					if (!it->punc_refcnt.compare_exchange_strong(expected,
      							PUNCREF_RETRIEVED)) {
      						//bug("bad punc refcnt");
      						it ++;
      						continue;
      					}
      					/* XXX: opt: we may check the next newer container and
      					 * coalesce punc as needed. benefit unclear though.
      					 */
      					//return it->punc;
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
      				//E("punc refcnt: %ld", punc_refcnt);
      				xzl_assert(false && "illegal punc_perfcnt?");
      			}

      	  	xzl_assert(false && "bug?");

        	} // while

      //  	dump_containers("no bundles");

        	return nullptr;
       }

       shared_ptr<BundleBase> getOneBundle_unordered_5(int node = -1) {
       	//TODO
      	//should not get punc here
      	//ref getOneBundle_left_unordered_4

      	//std::cout << __FILE__ << ": " <<  __LINE__ << " this trans is: " << this->getName() << std::endl;
        	shared_ptr<BundleBase> ret;

        	//unique_lock<mutex> conlock(mtx_unordered_containers);
        	//XXX rw
      	//unique_lock<mutex> conlock(mtx_container);
      	boost::shared_lock<boost::shared_mutex> rlock(mtx_container);
      	boost::upgrade_lock<boost::shared_mutex> ulock(mtx_container, boost::defer_lock);
      	unique_ptr<boost::upgrade_to_unique_lock<boost::shared_mutex>> pwlock
      								= nullptr; /* will std::move later */
      //start:

        	if (unordered_containers.empty()) {
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

        	auto it = unordered_containers.rbegin();

        	while (it != unordered_containers.rend()) {

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
      				/* has to be the oldest container since we can't have more
      				 * than one outstanding punc (the second oldest container, whose
      				 * punc is assigned, must still holds its punc).
      				 */
      				xzl_assert(it == unordered_containers.rbegin());
      //				it = containers_.erase(it);   /* erase w/ reverse iterator */
      				it ++;
      				/* XXX should not erase any containers here
      				it = list<bundle_container>::reverse_iterator(ordered_containers.erase(it.base()));
      				std::cout << __FILE__ << __LINE__ << "containers_.erase " << this->getName() << std::endl;
      				*/
      				continue;
      			} else if (punc_refcnt == PUNCREF_CANCELED) {
      				/* erase */
      				it ++;
      				/* XXX should not erase any containers here
      				it = list<bundle_container>::reverse_iterator(ordered_containers.erase(it.base()));
      				std::cout << __FILE__ << __LINE__ << "containers_.erase" << std::endl;
      				*/
      				continue;
        		} else if (punc_refcnt == PUNCREF_RETRIEVED) {
      				/* punc outstanding. it may become CONSUMED soon but we don't wait. */
      				it ++;
      				continue;
      			} else if (punc_refcnt == PUNCREF_ASSIGNED) { /* punc not emitted */
      				if (it == unordered_containers.rbegin()) {
      					/* oldest container. (older dead containers should have been
      					 * cleaned up now). okay to emit */
      //					auto r = it->punc_refcnt.fetch_sub(1);
      //					xzl_assert(r == 2); r = r;
      				#if 0
      					long expected = PUNCREF_ASSIGNED;
      					if (!it->punc_refcnt.compare_exchange_strong(expected,
      							PUNCREF_RETRIEVED)) {
      						bug("bad punc refcnt");
      					}
      					/* XXX: opt: we may check the next newer container and
      					 * coalesce punc as needed. benefit unclear though.
      					 */
      					return it->punc;
      				#endif
      					//XXX should not return punc here
      					it ++;
      					continue;
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
      				//E("punc refcnt: %ld", punc_refcnt);
      				xzl_assert(false && "illegal punc_perfcnt?");
      			}

      	  	xzl_assert(false && "bug?");

        	} // while

      //  	dump_containers("no bundles");

        	return nullptr;
       }

       shared_ptr<BundleBase> getOneBundle_5(int node = -1) {
       	//TODO

      	shared_ptr<BundleBase> ret;
      	ret = getOneBundle_ordered_5(node);
      	if(ret != nullptr){
      		return ret;
      	}else{
      		ret = getOneBundle_unordered_5(node);
      		return ret;
      	}

       }

       /***************************************************************
                   hym: end of type 5 trans
       ****************************************************************/

////////////////////////////////////////////////////
// hym: Following funcs are for Join's new source UnboundedInMemEvaluator_Join
///////////////////////////////////////////////////

//hym: depositOneBundleToJoin_L and depositOneBundleToJoin_R will be called by Source
//hym: Source deposit one bundle to Join's left_containers_in
  void depositOneBundleToJoin_L(shared_ptr<BundleBase> bundle, int node = -1) {

  	/* lock all containers until reliably getting the target container.
  	 * (in case concurrent additions of containers) */
  	//XXX rw
	//unique_lock<mutex> conlock(mtx_container);
	boost::shared_lock<boost::shared_mutex> rlock(mtx_container);
	boost::upgrade_lock<boost::shared_mutex> ulock(mtx_container, boost::defer_lock);
	unique_ptr<boost::upgrade_to_unique_lock<boost::shared_mutex>> pwlock
							= nullptr; /* will std::move later */
start:
  	if (left_containers_in.empty()) {
		//XXX rw
		if (!upgrade_locks(&rlock, &ulock, &pwlock)){
			goto start;
		}
  		left_containers_in.emplace_front( 1 /* side_info */ ); /* can't copy */
  	}

  	auto current_con = left_containers_in.begin();
	//current_con->set_side_info(this->get_side_info());
//	current_con->set_side_info(1); //left
  	
	//if (current_con->punc) { /* the latest container already has a punc */
  	if (current_con->getPuncSafe()) { /* the latest container already has a punc */
		//XXX rw
		if (!upgrade_locks(&rlock, &ulock, &pwlock)){
			goto start;
		}
  		
		xzl_assert(current_con->punc_refcnt != PUNCREF_UNDECIDED && "invalid punc refcnt");

		/* The latest container is already dead: clean up.
		 * NB: even if we see punc_refcnt == 1 here, it may be dead
		 * (decremented to 0) beneath us any time. It will be cleaned by future
		 * calls.
		 */
		if (current_con->punc_refcnt == PUNCREF_CONSUMED) {
			/* extensive sanity check... */

			//xzl_assert(current_con->punc && "dead container but no valid punc?");
			xzl_assert(current_con->getPuncSafe() && "dead container but no valid punc?");
			xzl_assert(current_con->bundles.size() == 0);
			xzl_assert(current_con->refcnt == 0);

			current_con = left_containers_in.erase(current_con);
			//std::cout << __FILE__ << __LINE__ << "containers_.erase" << std::endl;
			if (left_containers_in.empty()) {
				/* after cleanup, no container left. start a new one. */
				left_containers_in.emplace_front(1 /* side_info */ );
				current_con = left_containers_in.begin();
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
//		prefcnt = prefcnt;
		xzl_assert( prefcnt == PUNCREF_CONSUMED   /* latest container dies just now */
				|| prefcnt == PUNCREF_RETRIEVED		/* outstanding */
				|| prefcnt == PUNCREF_ASSIGNED		/* unretrieved yet */
		      );

		left_containers_in.emplace_front(1 /* side_info */ );
		current_con --;
	}

	/* Now @current_con points to an open container. Insert the bundle.
	 *
	 * @current_con won't go away since the punc won't be emitted
	 * concurrently.
	 */
open_container:
	/* should have any type of lock */
	//XXX rw
	xzl_assert(rlock.owns_lock() || ulock.owns_lock() || (pwlock && pwlock->owns_lock()));

	//current_con->putBundleUnsafe(bundle);
	current_con->putBundleSafe(bundle);
	this->IncBundleCounter();

  } //end depositOneBundleToJoin_L


//hym: Source deposit one bundle to Join's right_containers_in
  void depositOneBundleToJoin_R(shared_ptr<BundleBase> bundle, int node = -1) {

  	/* lock all containers until reliably getting the target container.
  	 * (in case concurrent additions of containers) */
  	//XXX rw
//  EE(" XXXXXXXXXX here ");
	//unique_lock<mutex> conlock(mtx_container);
	boost::shared_lock<boost::shared_mutex> rlock(mtx_container);
	boost::upgrade_lock<boost::shared_mutex> ulock(mtx_container, boost::defer_lock);
	unique_ptr<boost::upgrade_to_unique_lock<boost::shared_mutex>> pwlock
							= nullptr; /* will std::move later */
start:
  	if (right_containers_in.empty()) {
		//XXX rw
		if (!upgrade_locks(&rlock, &ulock, &pwlock)){
			goto start;
		}
  		right_containers_in.emplace_front(2 /* side_info */); /* can't copy */
  	}

  	auto current_con = right_containers_in.begin();
	//current_con->set_side_info(this->get_side_info());
//	current_con->set_side_info(2); //right
  	//if (current_con->punc) { /* the latest container already has a punc */
  	if (current_con->getPuncSafe()) { /* the latest container already has a punc */
		//XXX rw
		if (!upgrade_locks(&rlock, &ulock, &pwlock)){
			goto start;
		}

		xzl_assert(current_con->punc_refcnt != PUNCREF_UNDECIDED && "invalid punc refcnt");

		/* The latest container is already dead: clean up.
		 * NB: even if we see punc_refcnt == 1 here, it may be dead
		 * (decremented to 0) beneath us any time. It will be cleaned by future
		 * calls.
		 */
		if (current_con->punc_refcnt == PUNCREF_CONSUMED) {
			/* extensive sanity check... */

			//xzl_assert(current_con->punc && "dead container but no valid punc?");
			xzl_assert(current_con->getPuncSafe() && "dead container but no valid punc?");
			xzl_assert(current_con->bundles.size() == 0);
			xzl_assert(current_con->refcnt == 0);

			current_con = right_containers_in.erase(current_con);
			//std::cout << __FILE__ << __LINE__ << "containers_.erase" << std::endl;
			if (right_containers_in.empty()) {
				/* after cleanup, no container left. start a new one. */
				right_containers_in.emplace_front(2 /* side_info */);
				current_con = right_containers_in.begin();
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
//		prefcnt = prefcnt;
		xzl_assert( prefcnt == PUNCREF_CONSUMED   /* latest container dies just now */
				|| prefcnt == PUNCREF_RETRIEVED		/* outstanding */
				|| prefcnt == PUNCREF_ASSIGNED		/* unretrieved yet */
		      );

		right_containers_in.emplace_front(2 /* side_info */);
		current_con --;
	}

	/* Now @current_con points to an open container. Insert the bundle.
	 *
	 * @current_con won't go away since the punc won't be emitted
	 * concurrently.
	 */
open_container:
	/* should have any type of lock */
	//XXX rw
	xzl_assert(rlock.owns_lock() || ulock.owns_lock() || (pwlock && pwlock->owns_lock()));

	//current_con->putBundleUnsafe(bundle);
	current_con->putBundleSafe(bundle);
	this->IncBundleCounter();

  } //end depositOneBundleToJoin_R

// hym: this func is only called by source
//      deposit a punc to Join's left_containers_in
  void depositOnePuncToJoin_L(shared_ptr<Punc> punc, int node = -1) {

	//std::cout << "depositOnePunc() in trans: " << this->getName() << std::endl;
#ifdef MEASURE_LATENCY

//  	char buf[50];
//  	long cnt = this->getNumBundles();
//  	snprintf(buf, 50, "deposit to: %s (%ld bundles) ",
//  			this->name.c_str(), cnt);
//  	punc->mark(buf);
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
  	oss << "deposit to: " << this->name <<  " " << boost::posix_time::microsec_clock::local_time() <<" " << "\n (bundles: ";
  	for (auto & cnt : cnts) {
  		if (cnt.bundle >= 0)
  			oss << cnt.bundle << " ";
  		else {
  			auto r = cnt.bundle + PUNCREF_MAX;
  			oss << puncref_key(r) << " ";
  		}
		oss << cnt.side_info << " ";
  		//  			if (r >= 0) { /* punc ever assigned. show punc ts */
			if (cnt.punc != max_date_time) {
					oss << to_simple_string(cnt.punc) << "\t";
			}
  	}
  	oss << ")";
  	punc->mark(oss.str());

//  	punc->mark("deposit to: " + this->name);
#endif
	//XXX rw
   	//unique_lock<mutex> conlock(mtx_container);
	boost::shared_lock<boost::shared_mutex> rlock(mtx_container);
	boost::upgrade_lock<boost::shared_mutex> ulock(mtx_container, boost::defer_lock);
	unique_ptr<boost::upgrade_to_unique_lock<boost::shared_mutex>> pwlock
							= nullptr; /* will std::move later */
start:
   	if (left_containers_in.empty()) {
		//std::cout << __FILE__ << ":" << __LINE__ << " containers_.empty() wrong???" << std::endl;
   		/* Whole transform drained. we don't even have pending flushed
   		 * bundles (otherwise there will be containers).
   		 * Create a new container and seal it with @wm.
   		 */
//   		containers_.push_front(bundle_container());
		
		//XXX rw
		if (!upgrade_locks(&rlock, &ulock, &pwlock)){
			goto start;
		}
   		left_containers_in.emplace_front(1 /* side_info */ );
   		//containers_.begin()->setPunc(punc);
   		left_containers_in.begin()->setPuncSafe(punc);
   		return;
   	}
	//std::cout << __FILE__ << ":" << __LINE__ << std::endl;
   	/* assign or update the wm of the most recent container */
   	auto current = left_containers_in.begin();
	if (current->punc_refcnt == PUNCREF_RETRIEVED) {
		/* punc already emitted. start a new container. */

		//XXX rw
		if (!upgrade_locks(&rlock, &ulock, &pwlock)){
			goto start;
		}

		xzl_assert(current->punc && "must have a valid punc ptr");
		xzl_assert(!current->refcnt && "all bundles must be consumed");
		//   		containers_.push_front(bundle_container());
		xzl_bug("shouldn't happen");
		left_containers_in.emplace_front(1 /* side_info */);
		//containers_.begin()->setPunc(punc);
		left_containers_in.begin()->setPuncSafe(punc);
		return;
	}
	//std::cout << __FILE__ << ":" << __LINE__ << std::endl;
	/* clean up any dead container. */
	if (current->punc_refcnt == PUNCREF_CONSUMED) {
		/* extensive sanity check... */

		//XXX rw
		if (!upgrade_locks(&rlock, &ulock, &pwlock)){
			goto start;
		}

		//xzl_assert(current->punc && "dead container but no valid punc?");
		xzl_assert(current->getPuncSafe() && "dead container but no valid punc?");
		xzl_assert(current->bundles.size() == 0);
		xzl_assert(current->refcnt == 0);

		current = left_containers_in.erase(current);
		//std::cout << __FILE__ << __LINE__ << "containers_.erase" << std::endl;
		if (left_containers_in.empty()) {
			/* after cleanup, no container left. start a new one. */
			left_containers_in.emplace_front();
			xzl_bug("shouldn't happen");
			//containers_.begin()->setPunc(punc);
			left_containers_in.begin()->setPuncSafe(punc);
			return;
		} else {
			/* current now points to the 2nd most recent container. */
			xzl_assert(current->punc_refcnt != 0
					&& "bug: can't have two dead containers");
		}
	}
	//std::cout << __FILE__ << ":" << __LINE__ << std::endl;
	//if (!current->punc || current->punc_refcnt == PUNCREF_ASSIGNED) {
	if (!current->getPuncSafe() || current->punc_refcnt == PUNCREF_ASSIGNED) {
		//std::cout << __FILE__ << ":" << __LINE__ << std::endl;
		//current->setPunc(punc); /* assign or overwrite punc */
		current->setPuncSafe(punc); /* assign or overwrite punc */
		//current->set_side_info(this->get_side_info());
	} else {
		xzl_assert(false && "bug?");
	}
  }


// hym: this func is only called by source
//      deposit a punc to Join's left_containers_in
  void depositOnePuncToJoin_R(shared_ptr<Punc> punc, int node = -1) {

	//std::cout << "depositOnePunc() in trans: " << this->getName() << std::endl;
#ifdef MEASURE_LATENCY

//  	char buf[50];
//  	long cnt = this->getNumBundles();
//  	snprintf(buf, 50, "deposit to: %s (%ld bundles) ",
//  			this->name.c_str(), cnt);
//  	punc->mark(buf);
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
  	oss << "deposit to: " << this->name <<  " " << boost::posix_time::microsec_clock::local_time() <<" " << "\n (bundles: ";
  	for (auto & cnt : cnts) {
  		if (cnt.bundle >= 0)
  			oss << cnt.bundle << " ";
  		else {
  			auto r = cnt.bundle + PUNCREF_MAX;
  			oss << puncref_key(r) << " ";
  		}
		oss << cnt.side_info << " ";
  		//  			if (r >= 0) { /* punc ever assigned. show punc ts */
			if (cnt.punc != max_date_time) {
					oss << to_simple_string(cnt.punc) << "\t";
			}
  	}
  	oss << ")";
  	punc->mark(oss.str());

//  	punc->mark("deposit to: " + this->name);
#endif
	//XXX rw
   	//unique_lock<mutex> conlock(mtx_container);
	boost::shared_lock<boost::shared_mutex> rlock(mtx_container);
	boost::upgrade_lock<boost::shared_mutex> ulock(mtx_container, boost::defer_lock);
	unique_ptr<boost::upgrade_to_unique_lock<boost::shared_mutex>> pwlock
							= nullptr; /* will std::move later */
start:
   	if (right_containers_in.empty()) {
		//std::cout << __FILE__ << ":" << __LINE__ << " containers_.empty() wrong???" << std::endl;
   		/* Whole transform drained. we don't even have pending flushed
   		 * bundles (otherwise there will be containers).
   		 * Create a new container and seal it with @wm.
   		 */
//   		containers_.push_front(bundle_container());
		
		//XXX rw
		if (!upgrade_locks(&rlock, &ulock, &pwlock)){
			goto start;
		}
		xzl_bug("shouldn't happen");
   		right_containers_in.emplace_front();
   		//containers_.begin()->setPunc(punc);
   		right_containers_in.begin()->setPuncSafe(punc);
   		return;
   	}
	//std::cout << __FILE__ << ":" << __LINE__ << std::endl;
   	/* assign or update the wm of the most recent container */
   	auto current = right_containers_in.begin();
	if (current->punc_refcnt == PUNCREF_RETRIEVED) {
		/* punc already emitted. start a new container. */

		//XXX rw
		if (!upgrade_locks(&rlock, &ulock, &pwlock)){
			goto start;
		}

		xzl_assert(current->punc && "must have a valid punc ptr");
		xzl_assert(!current->refcnt && "all bundles must be consumed");
		//   		containers_.push_front(bundle_container());
		right_containers_in.emplace_front();
		xzl_bug("shouldn't happen");
		//containers_.begin()->setPunc(punc);
		right_containers_in.begin()->setPuncSafe(punc);
		return;
	}
	//std::cout << __FILE__ << ":" << __LINE__ << std::endl;
	/* clean up any dead container. */
	if (current->punc_refcnt == PUNCREF_CONSUMED) {
		/* extensive sanity check... */

		//XXX rw
		if (!upgrade_locks(&rlock, &ulock, &pwlock)){
			goto start;
		}

		//xzl_assert(current->punc && "dead container but no valid punc?");
		xzl_assert(current->getPuncSafe() && "dead container but no valid punc?");
		xzl_assert(current->bundles.size() == 0);
		xzl_assert(current->refcnt == 0);

		current = right_containers_in.erase(current);
		//std::cout << __FILE__ << __LINE__ << "containers_.erase" << std::endl;
		if (right_containers_in.empty()) {
			xzl_bug("shouldn't happen");
			/* after cleanup, no container left. start a new one. */
			right_containers_in.emplace_front();
			//containers_.begin()->setPunc(punc);
			right_containers_in.begin()->setPuncSafe(punc);
			return;
		} else {
			/* current now points to the 2nd most recent container. */
			xzl_assert(current->punc_refcnt != 0
					&& "bug: can't have two dead containers");
		}
	}
	//std::cout << __FILE__ << ":" << __LINE__ << std::endl;
	//if (!current->punc || current->punc_refcnt == PUNCREF_ASSIGNED) {
	if (!current->getPuncSafe() || current->punc_refcnt == PUNCREF_ASSIGNED) {
		//std::cout << __FILE__ << ":" << __LINE__ << std::endl;
		//current->setPunc(punc); /* assign or overwrite punc */
		current->setPuncSafe(punc); /* assign or overwrite punc */
		//current->set_side_info(this->get_side_info());
	} else {
		xzl_assert(false && "bug?");
	}
  }//end depositOnePuncToJoin_R










////////////////////END For new Join////////////////////////////////////////////
#endif /* TRANSFORMS_MULTI_H_ */
