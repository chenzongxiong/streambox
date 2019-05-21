#ifndef VALUECONTAINER_H_
#define VALUECONTAINER_H_

#include "tbb/concurrent_vector.h"


/*
 * the container for values that may be discontig in memory.
 *
 * this does not contain the actual values.
 *
 * this d/s does not require ts to be saved alongside with each v.
 * If want to save timestamp for each value, V can be Record<>.
 * @V: value element type.
 *
 * NB:
 * by default, the container is MT safe; it can be optionally unsafe.
 * basic container is not MT safe. caller ensures MT safety unless for -safe methods.
 */

/*
 * the container for values in contig memory.
 *
 * xzl: vector<> has not virtual dtor. So this risks resource leakage
 * if the object is deleted through the ptr to vector<>. But just
 * don't do that.
 *
 * If reserving memory, it appears all allocators (tbb, jemalloc, std...) see heavy contention.
 * Is this because reserving a specific memory size stresses a particular bucket and therefore
 * leads to contention?
*/

#ifdef USE_FOLLY_VECTOR /* using folly */
#include "folly/FBVector.h"
template<typename V>
class BasicValueContainer : public folly::fbvector<V> {

  /* we have tons of memory and hate growing */
  const unsigned long init_capacity = (16 * 1024) / sizeof(V);

public:
  ptime min_ts = max_date_time; // user has to maintain this.

    ptime max_ts = min_date_time; // measure ysb latency

  BasicValueContainer() {
  	// detrimental to performance... why
//  	folly::fbvector<V>::reserve(init_capacity);
  }
};
#else /* fallback to std */
template<typename V>
class BasicValueContainer : public std::vector<V> {

  /* we have tons of memory and hate growing */
  const unsigned long init_capacity =
  		(16 * 1024) / sizeof(V);

public:
  ptime min_ts = max_date_time; // user has to maintain this.
  ptime max_ts = min_date_time; // user has to maintain this.
  BasicValueContainer() {
  	// detrimental to performance... cause a lot of wait in new(). this even happens to jemalloc
//  	std::vector<V>::reserve(init_capacity);
  }
};
#endif

#if 0  /* reservation (allocation) leads to lots of cpu spin time, which also kills performance */
#include "tbb/scalable_allocator.h"
template<typename V>
class BasicValueContainer : public std::vector<V, tbb::scalable_allocator<V>> {

  /* we have tons of memory and hate growing */
  const unsigned long init_capacity =
  		(16 * 1024) / sizeof(V);

public:
  ptime min_ts = max_date_time; // user has to maintain this.

//  BasicValueContainer(unsigned long capa = 0) {
//  	if (capa)
//  		std::vector<V,tbb::scalable_allocator<V>>::reserve(capa);  // destroy performance. why? -- stateful transform readers all wait for long
//  }

  BasicValueContainer() {
  	// destroy performance. why? -- stateful transform readers all wait for long
//		std::vector<V,tbb::scalable_allocator<V>>::reserve(init_capacity);
	}
};
#endif


#if 0
template<typename V>
class BasicValueContainer {
  std::vector<V> vec;
public:
  ptime min_ts = max_date_time; // user has to maintain this.
  void push_back(V const & v) {
    vec.push_back(v);
  }
};
#endif

/* good for holding many values */
template<typename V,
				class Container=tbb::concurrent_vector<shared_ptr<BasicValueContainer<V>>>>
struct ValueContainer
{
  using BasicValueContainerT = BasicValueContainer<V>;

  // the basic container that is memory contig, as produced by one
  // thread.

//  vector<shared_ptr<BasicValueContainerT>> vals;
//  tbb::concurrent_vector<shared_ptr<BasicValueContainerT>> vals;
  Container vals;

  // there's a min_ts for the entire container (multiple basic containers).
  // the reason is that we expect that the container (per window per key)
  // will not be further splitted among workers.
  // (although multiple containers can merge).
  // therefore, we don't need to keep finer-grained ts for each basic container.
  //
  // We keep min_ts here, instead of the KVS, so that operation on
  // the value container can update teh min_ts.
  // Besides, key does not have ts.cp
  ptime min_ts = max_date_time;

    // ptime max_ts_tuple_in_win = min_date_time;
    ptime max_ts = min_date_time;

//  mutex mtx_; /* protect the entire d/s */

  void merge_safe(ValueContainer const & other) {
//  	unique_lock<mutex> lock(mtx_);
//  	vals.insert(vals.end(), other.vals.begin(), other.vals.end());
//       zxchen
//       std::cout << "merge safe is called" << std::endl;
  	vals.grow_by(other.vals.begin(), other.vals.end());
  	min_ts = std::min(min_ts, other.min_ts); /* this is unsafe...XXX */

    max_ts = std::max(max_ts, other.max_ts);
  }

  // add a value to the container, which by default, save to the
  // first basic container.
  void add_v(V const & v, ptime const & ts) {
      // std::cout << "add_v is called: " << std::endl;
    if (vals.size() == 0) {
        assert(!_begin); // cached "begin" iterator has to be null
        vals.push_back(make_shared<BasicValueContainerT>());
        /* BasicValueContainerT itself should reserve enough capacity */
    }
    assert(vals[0]);
    vals[0]->push_back(v);   /* hot */

    // update v basic container's min_ts
    if (ts < vals[0]->min_ts)
      vals[0]->min_ts = ts;

    if (ts > vals[0]->max_ts) {
        vals[0]->max_ts = ts;
    }
    // update v container's min_ts
    if (ts < min_ts)
      min_ts = ts;

    if (ts > max_ts) {
        max_ts = ts;
    }
    // destroy & invalidate the "end" iterator but not the "begin" one.
    // since if the container was empty before add_v(), the "begin"
    // iterator is already invalid.
    if (_end)
//      _end.reset(nullptr);
//      _end.reset(); // will become nullptr
      _end = nullptr;
  }

  void add_basic_container(shared_ptr<BasicValueContainerT> const & container) {
      cout << "add basic container: " << endl;

    vals.push_back(container);
    if (min_ts > container->min_ts)
      min_ts = container->min_ts;

    // if (min_ts > container->min_ts)
    //   min_ts = container->min_ts;

    if (_end)
//      _end.reset(nullptr);
//      _end.reset();
      _end = nullptr;
  }

  uint64_t size() {
    uint64_t ret = 0;
    for(auto & basic_container_ptr : vals) {
        ret += basic_container_ptr->size();
    }
    return ret;
  }

  //////////////////////////////////////////////////////////////

  class iterator {
  private:
    ValueContainer const * vcontainer;
    // index of the current basic container. we don't save b container ptr
    uint64_t basic_container_index;
    // iterator within the current basic container
    typename BasicValueContainerT::const_iterator basic_container_iterator;

  public:
    iterator(ValueContainer const * vcontainer, int64_t index)
      : vcontainer(vcontainer) {
      assert(vcontainer);
      assert(vcontainer->vals.size()); // must have at least one basic container

      auto & bcontainers = vcontainer->vals;

      if (index == 0) {
          basic_container_index = 0;
          basic_container_iterator = bcontainers[basic_container_index]->begin();
      } else if (index == -1) {
          basic_container_index = bcontainers.size()-1;
          basic_container_iterator = bcontainers[basic_container_index]->end();
      }else {
          // force iterator can only be constructed as begin() and end()
          EE("not implemented");
          abort();
      }
    }

    bool operator !=(iterator const& other) const {
      return (!((basic_container_index == other.basic_container_index)
              && (basic_container_iterator == other.basic_container_iterator)));
    }

    iterator& operator++() {
      auto & bcontainers = vcontainer->vals;
      auto & bcontainer = bcontainers[basic_container_index];

      assert(basic_container_index < vcontainer->vals.size());

      // 1. if we point to the end() of the last b container, we shouldn't ++;
      // 2. we shouldn't point to the end() of any other b containers
      //    (when that happens, we should automatically move to the begin()
      //     of the next b container)
      assert(basic_container_iterator != bcontainer->end());

      basic_container_iterator ++;
      // hit the end of the current b container
      if (basic_container_iterator == bcontainer->end()) {
          if (basic_container_index == vcontainer->vals.size() - 1) {
              // we are at the end of the last b container. do nothing.
              ;
          } else {
              // move to the beginning of the next b container.
              basic_container_index ++;
              basic_container_iterator = \
                  bcontainers[basic_container_index]->begin();
              // what if the new b container is empty? have to keep moving
              // to the next b container. to be implemented...
              assert(bcontainers[basic_container_index]->size());
          }
      }
      return *this;
    }

    V const & operator* () const {
      return *basic_container_iterator;
    }

  };

  //////////////////////////////////////////////////////////////

#if 1
  /* Shared ptr version. This caches begin() and end() iterators
   * inside the v container. It, however, breaks the const-ness of
   * iterator.
   */
  // will nullptr by default.
  shared_ptr<iterator> _begin;
  shared_ptr<iterator> _end;

  // the following two are *not* const as they will modify the
  // cached begin() and end()
  iterator const & begin () {
    if (!_begin)
      _begin = make_shared<iterator>(this, 0);
    return *_begin;
  }

  iterator const &  end() {
    if (!_end)
        _end = make_shared<iterator>(this, -1);
    return *_end;
  }

  /* naive: construct iterators on the fly. don't cache them inside.
   * this keeps the constness of the iterators.
   * the caller should be cautious not to call cend() in each iteration.
   */
  iterator const cbegin () const {
    return iterator(this, 0);
  }

  iterator const cend() const {
    return iterator(this, -1);
  }

#endif

  using UnsafeValueContainerT = ValueContainer <V,
  		std::vector<shared_ptr<BasicValueContainer<V>>>>;
};

/* unsafe variant */
template<typename V>
using ValueContainerUnsafe = ValueContainer <V,
		std::vector<shared_ptr<BasicValueContainer<V>>>>;

/*
 * same interface as ValueContainer, but only wraps a single value to avoid
 * all the overhead
 */
template<typename V>
struct SingleValueContainer
{
	V v;
	ptime ts;

	/* add a value, which overwrites existing value */
	void add_v(V const & v, ptime const & ts) {
		this->v = v;
		this->ts = ts;
	}

};

/* container for a small number of values. must have same interface as ValueContainer */
template<typename V>
struct SimpleValueContainerUnsafe
{
	vector<V> vals;
  ptime min_ts = max_date_time;  /* we may get ride of this */
    ptime max_ts = min_date_time;  /* we may get ride of this */

  using const_iterator = typename vector<V>::const_iterator;
  using iterator = typename vector<V>::iterator;

  void merge_unsafe(SimpleValueContainerUnsafe const & other) {
  	vals.insert(vals.end(), other.vals.begin(), other.vals.end());
  	min_ts = std::min(min_ts, other.min_ts); /* this is unsafe...XXX */
  	max_ts = std::max(max_ts, other.max_ts); /* this is unsafe...XXX */
  }

	/* add a value, which overwrites existing value */
	void add_v(V const & v, ptime const & ts) {
        cout << "simple value container unsafe add_v is called" << endl;
		vals.push_back(v);
    if (ts < min_ts)
      min_ts = ts;
    if (ts > max_ts)
        max_ts = ts;
	}

  uint64_t size() {
  	return vals.size();
  }

  // SimpleValueContainerUnsafe(const SimpleValueContainerUnsafe & other)
  // 	: vals(other.vals), min_ts(other.min_ts) {
  // }
  SimpleValueContainerUnsafe(const SimpleValueContainerUnsafe & other)
      : vals(other.vals), min_ts(other.min_ts), max_ts(other.max_ts) {
  }

  SimpleValueContainerUnsafe() { }

  const_iterator cbegin() const {
    return vals.cbegin();
  }

  const_iterator cend() const {
    return vals.cend();
  }

  iterator begin() {
  	return vals.begin();
  }

  iterator end() {
  	return vals.end();
  }

//  using UnsafeValueContainerT = SimpleValueContainerUnsafe<V>;
};

/* same as above, but w lock */
template<typename V>
struct SimpleValueContainer
{
	vector<V> vals;
  ptime min_ts = max_date_time;  /* we may get ride of this */
    ptime max_ts = min_date_time;  /* we may get ride of this */
  mutex mtx_; /* protect the entire d/s */

  using const_iterator = typename vector<V>::const_iterator;
  using iterator = typename vector<V>::iterator;
  using UnsafeValueContainerT = SimpleValueContainerUnsafe<V>;

  /* @others is safe */
  void merge_safe(SimpleValueContainer const & other) {
  	unique_lock<mutex> lock(mtx_);
  	vals.insert(vals.end(), other.vals.begin(), other.vals.end());
  	min_ts = std::min(min_ts, other.min_ts); /* this is unsafe...XXX */

  	max_ts = std::max(max_ts, other.max_ts); /* this is unsafe...XXX */
  }

  /* @others is unsafe */
  void merge_safe(UnsafeValueContainerT const & other) {
  	unique_lock<mutex> lock(mtx_);
  	vals.insert(vals.end(), other.vals.begin(), other.vals.end());
  	min_ts = std::min(min_ts, other.min_ts); /* this is unsafe...XXX */

  	max_ts = std::max(max_ts, other.max_ts); /* this is unsafe...XXX */

  }

	/* add a value, which overwrites existing value */
	void add_v(V const & v, ptime const & ts) {
        cout << "simple value container add_v is called" << endl;
		vals.push_back(v);
    if (ts < min_ts)
      min_ts = ts;
    if (ts > max_ts)
        max_ts = ts;
	}


  uint64_t size() {
  	return vals.size();
  }

  // SimpleValueContainer(const SimpleValueContainer & other)
  //     : vals(other.vals), min_ts(other.min_ts) {
  // }

  SimpleValueContainer(const SimpleValueContainer & other)
      : vals(other.vals), min_ts(other.min_ts), max_ts(other.max_ts){
  }

  SimpleValueContainer() { }

  const_iterator cbegin() const {
    return vals.cbegin();
  }

  const_iterator cend() const {
    return vals.cend();
  }

  iterator begin() {
  	return vals.begin();
  }

  iterator end() {
  	return vals.end();
  }

};

#endif /* VALUECONTAINER_H_ */
