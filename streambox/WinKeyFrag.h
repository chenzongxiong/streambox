#ifndef WINDOWKEYEDFRAGMENT_H_
#define WINDOWKEYEDFRAGMENT_H_

#include "tbb/concurrent_unordered_map.h"

/*
 * A window and a vector of Records that are associated with the window
 * (note that these records are contig in memory)
 *
 * "fragment" means that the records here are not necessarily the window's all
 * records.
 *
 * @T: element type, e.g. long
*/
template <typename T>
struct WindowFragment
{
	using RecordT = Record<T>;

	Window w;

	// the ts of the oldest Record in the output. used to compute watermark
	ptime min_ts = max_date_time;
    // the ts of the latest Record in the output, use to compute latency of yahoo streaming benchmark
    ptime max_ts = min_date_time;

	vector<RecordT> vals;

#ifndef NDEBUG
	int id; // for debugging/tracking
#endif

	// add a record (assuming that it falls into the window)
	void add(const RecordT& val) {
		if (val.ts < min_ts)
			min_ts = val.ts;
        if (val.ts > max_ts) {
            max_ts = val.ts;
        }

		vals.push_back(val);
	}

	WindowFragment(const Window& w) : w(w) {
#ifndef NDEBUG
		id = atomic_fetch_add(&cnt, 1);
#endif
	}

	void print(ostream &os) const {
		os << "WindowFragment min_ts: " << min_ts << " data: {";
		for (auto & rec : vals) {
			os << rec.data << ",";
		}
		os << "} end" << endl;
	}

	WindowFragment & merge(const WindowFragment & other) {
		assert(this->w.start == other.w.start
				&& this->w.duration == other.w.duration);
		/* merge two vectors */
		this->vals.insert(this->vals.end(), other.vals.begin(), other.vals.end());
		this->min_ts = (this->min_ts < other.min_ts) ? this->min_ts : other.min_ts;


		return *this;
	}
};

/* One window, multiple keys, and the value containers associated with
 * the keys.
 * Note this does not save KV pairs but K,ValueContainer pairs.
 */

#if 0 /* the default one (good for wordcount) */
template <class KVPair,
	 class UnorderedMap =
	 tbb::concurrent_unordered_map<decltype(KVPair::first),
	 ValueContainer<decltype(KVPair::second)>>
	 /* we defined hasher for fbstring */
	 >
#endif
#if 0
	 /* map is creek's partitioned HT */
	 template <class KVPair,
	 class UnorderedMap =
	 creek::concurrent_unordered_map<decltype(KVPair::first),
	 ValueContainer<decltype(KVPair::second)>>
	 >
#endif
#if 0
	 /* creek + simple vcontainer */
	 template <class KVPair,
	 class UnorderedMap =
	 creek::concurrent_unordered_map<decltype(KVPair::first),
	 SimpleValueContainer<decltype(KVPair::second)>>
	 >
#endif
	 /* must specify both args. */
	 template <class KVPair, class UnorderedMap>
	 struct WinKeyFrag
{
	using K = decltype(KVPair::first);
	using V = decltype(KVPair::second);
	using RecordKV = Record<KVPair>;
	//  using ValueContainerT = ValueContainer<V>;
	using ValueContainerT = typename UnorderedMap::mapped_type;
	using UnsafeValueContainerT = typename ValueContainerT::UnsafeValueContainerT;

	/* keys and their value containers (for a given window).
	 * Can be unordered_map but need extra fine-grained lock. */
	//  using KeyMap = tbb::concurrent_unordered_map<K, ValueContainerT>;
	using KeyMap = UnorderedMap;

	Window w;

	/* each ValueCotainer (contains no actual values; should be small) is
	 * instantiated in place
	 *
	 * tbb::concurrent_unordered_map does not help since ValueContainer also
	 * needs concurrency protection...
	 *
	 */
	//  tbb::concurrent_unordered_map<K, ValueContainerT> vals;
	KeyMap vals;

	// // the ts of the oldest Record in the output. used to compute watermark
	ptime min_ts = max_date_time;
    // the ts of the latest Record in the output, use to compute latency of yahoo streaming benchmark
    ptime max_ts = min_date_time;


	boost::shared_mutex mtx_; /* coarse grained mtx to protect the entire d/s */

	/* explict ctor needed to avoid copying mtx */
	WinKeyFrag(WinKeyFrag const & other)
		: w(other.w), vals(other.vals) { }

	WinKeyFrag() { }

	WinKeyFrag & operator= (WinKeyFrag const & other) {
		w = other.w;
		vals = other.vals;
		/* do nothing with mtx */
		return *this;
	}

#if 0
	void add_basic_container_safe(K const & key,
			shared_ptr<BasicValueContainerT> const & bc_ptr) {
		unique_lock<mutex> lck(mtx_);
		vals[key].add_basic_container(bc_ptr);
	}
#endif

	/* merge by copy. take rlock if no window addition needed; otherwise take
	 * wlock
	 */
#if 0
	void add_vcontainer_safe(K const & key,
			ValueContainerT const & others) {

		boost::shared_lock<boost::shared_mutex> rlock(mtx_);

		if (vals.count(key) == 0) { /* needs to insert new key */
			rlock.unlock();

			boost::unique_lock<boost::shared_mutex> wlock(mtx_);
			if (vals.count(key) == 0) { /* with wlock, check again. */
				vals.emplace(key, ValueContainerT{} /* ctor and move */);
			}
			wlock.unlock();
			/* once set, the key won't be gone. so feel safe to read. */
			rlock.lock();
		}

		vals[key].merge_safe(others);
	}
#endif


	/* this is hot. is it because the copy of vcontainers are constructed? */
	void add_vcontainer_safe(K const & key,
			ValueContainerT const & others) {
		/* @vals is MT safe. the following inserts @key if it is not present */
		/* ver 1 */
		//  	vals[key].merge_safe(others);

		/* ver 2 */
        // zxchen
        // std::cout << "add vcontainer safe is called" << std::endl;
		auto it = vals.find(key);  /* try find -- avoid ctor vcontainer */
		if (it != vals.end()) {
			it->second.merge_safe(others);
		} else {
			/* no found, insert. it may fail, but we have valid iter anyway */
			auto ret = vals.insert(make_pair(key, ValueContainerT{}));
			//  		auto ret = vals.emplace(key, ValueContainerT{});
			/* we always have a valid iterator regardless of insertion fails/succeeds */
			ret.first->second.merge_safe(others); /* ret.first: iterator. ret.first->second: vcontainer */
		}
	}

	void add_vcontainer_safe(K const & key,
			UnsafeValueContainerT const & others) {
		/* @vals is MT safe. the following inserts @key if it is not present */
		/* ver 1 */
		//  	vals[key].merge_safe(others);

		/* ver 2 */
		//  	auto it = vals.find(key);  /* try find -- avoid ctor vcontainer */
		//  	if (it != vals.end()) {
		//  		it->second.merge_safe(others);
		//  	} else {
		//			/* no found, insert. it may fail, but we have valid iter anyway */
		//			auto ret = vals.insert(make_pair(key, ValueContainerT{}));
		////  		auto ret = vals.emplace(key, ValueContainerT{});
		//			/* we always have a valid iterator regardless of insertion fails/succeeds */
		//			ret.first->second.merge_safe(others); /* ret.first: iterator. ret.first->second: vcontainer */
		//  	}

		/*v3 no []. direct insert (does not help much)*/
		auto ret = vals.insert(make_pair(key, ValueContainerT{}));
		/* we always have a valid iterator regardless of insertion fails/succeeds */
		ret.first->second.merge_safe(others); /* ret.first: iterator. ret.first->second: vcontainer */
	}

	// should we keep the Record's ts?
	void add_record_safe(const RecordKV& rec) {
		unique_lock<mutex> lck(mtx_);

		/* v1 */
		// only store value (but not key) in the value container
		//    vals[rec.data.first].add_v(rec.data.second, rec.ts);

		/* v2 */
		auto ret = vals.insert(make_pair(rec.data.first, ValueContainerT{}));
		ret.first->add_v(rec.data.second, rec.ts);
	}

	/* hot -- 1. write to hashtable? 2. instantiate k/v? */
	void add_record_unsafe(const RecordKV& rec) {
		// only store value (but not key) in the value container
		// if @vals is MT safe, this is in fact safe (no lock needed)

		/* ver 1 */
		//		vals[rec.data.first].add_v(rec.data.second, rec.ts);

		/* ver 2 */
		//  	auto & kv = vals[rec.data.first];
		//  	kv.add_v(rec.data.second, rec.ts);

		/* ver 3: avoid using [] */
		auto it = vals.find(rec.data.first);  /* try find -- avoid ctor vcontainer */
		if (it != vals.end()) {
			it->second.add_v(rec.data.second, rec.ts);
		} else {
			/* no found, insert. it may fail, but we have valid iter anyway */
			auto ret = vals.insert(make_pair(rec.data.first, ValueContainerT{}));
			//  		auto ret = vals.emplace(rec.data.first, ValueContainerT{});
			/* ret.first: iterator. ret.first->second: vcontainer */
			ret.first->second.add_v(rec.data.second, rec.ts);
		}

		/* v4: no []. direct insert */
		/* no found, insert. it may fail, but we have valid iter anyway */
		//		auto ret = vals.insert(make_pair(rec.data.first, ValueContainerT{}));
		//		/* ret.first: iterator. ret.first->second: vcontainer */
		//		ret.first->second.add_v(rec.data.second, rec.ts);
	}

	void add_kv_unsafe(KVPair const & kv, ptime const & ts) {
		/* ver1 */
		//  	vals[kv.first].add_v(kv.second, ts);

		/* ver2: avoid using [] */
        // cout << "add_kv_unsafe: " << endl;
		auto ret = vals.insert(make_pair(kv.first, ValueContainerT{}));
		ret.first->second.add_v(kv.second, ts);
	}
    void add_kv_unsafe(KVPair const &kv, ptime const &ts1, ptime const &ts2) {
		auto ret = vals.insert(make_pair(kv.first, ValueContainerT{}));
		ret.first->second.add_v(kv.second, ts1);
        // ret.first->second.add_v(kv.second, ts2);
        ret.first->second.add_v(0, ts2);

        if (min_ts > ts1) {
            min_ts = ts1;
        }
        if (max_ts < ts2) {
            max_ts = ts2;
        }
    }

	// go through each key (the min_ts is saved with the corresponding v
	// container)
	// should we save min_ts in this object?
	// ptime min_ts() {
	// 	ptime min_ts = max_date_time;
	// 	for (auto & kv : vals) {
	// 		if (kv.second.min_ts < min_ts)
	// 			min_ts = kv.second.min_ts;
	// 	}
	// 	return min_ts;
	// }
};

/////////////////////////////////////////////////////////////////////////////////////

/* for eval local use. does not have some merge methods. but can be combined into the shared version */
template <class KVPair, class UnorderedMap>
struct WinKeyFragLocal
{
	using K = decltype(KVPair::first);
	using V = decltype(KVPair::second);
	using RecordKV = Record<KVPair>;
	using ValueContainerT = typename UnorderedMap::mapped_type;

	/* keys and their value containers (for a given window).
	 * Can be unordered_map but need extra fine-grained lock. */
	//  using KeyMap = tbb::concurrent_unordered_map<K, ValueContainerT>;
	using KeyMap = UnorderedMap;

	Window w;

	/* each ValueCotainer (contains no actual values; should be small) is
	 * instantiated in place
	 *
	 * tbb::concurrent_unordered_map does not help since ValueContainer also
	 * needs concurrency protection...
	 *
	 */
	//  tbb::concurrent_unordered_map<K, ValueContainerT> vals;
	KeyMap vals;

	//  boost::shared_mutex mtx_; /* coarse grained mtx to protect the entire d/s */
	//
	//  /* explict ctor needed to avoid copying mtx */
	//  WindowKeyedFragment(WindowKeyedFragment const & other)
	//  	: w(other.w), vals(other.vals) { }
	//
	//  WindowKeyedFragment() { }

	WinKeyFragLocal & operator= (WinKeyFragLocal const & other) {
		w = other.w;
		vals = other.vals;
		/* do nothing with mtx */
		return *this;
	}

	/* this is hot. is it because the copy of vcontainers are constructed? */
	void add_vcontainer_safe(K const & key,
			ValueContainerT const & others) {
		/* @vals is MT safe. the following inserts @key if it is not present */
		/* ver 1 */
		//  	vals[key].merge_safe(others);

		/* ver 2 */
        // std::cout << "add_vcontainer_safe is called" << std::endl;
		auto it = vals.find(key);  /* try find -- avoid ctor vcontainer */
		if (it != vals.end()) {
			it->second.merge_safe(others);
		} else {
			/* no found, insert. it may fail, but we have valid iter anyway */
			auto ret = vals.insert(make_pair(key, ValueContainerT{}));
			//  		auto ret = vals.emplace(key, ValueContainerT{});
			/* we always have a valid iterator regardless of insertion fails/succeeds */
			ret.first->second.merge_safe(others); /* ret.first: iterator. ret.first->second: vcontainer */
		}
	}

	//  // should we keep the Record's ts?
	//  void add_record_safe(const RecordKV& rec) {
	//  	unique_lock<mutex> lck(mtx_);
	//
	//  	/* v1 */
	//  	// only store value (but not key) in the value container
	////    vals[rec.data.first].add_v(rec.data.second, rec.ts);
	//
	//    /* v2 */
	//  	auto ret = vals.insert(make_pair(rec.data.first, ValueContainerT{}));
	//  	ret.first->add_v(rec.data.second, rec.ts);
	//  }

	/* hot -- 1. write to hashtable? 2. instantiate k/v? */
	void add_record_unsafe(const RecordKV& rec) {
		// only store value (but not key) in the value container
		// if @vals is MT safe, this is in fact safe (no lock needed)

		/* ver 1 */
		//		vals[rec.data.first].add_v(rec.data.second, rec.ts);

		/* ver 2 */
		//  	auto & kv = vals[rec.data.first];
		//  	kv.add_v(rec.data.second, rec.ts);

		/* ver 3: avoid using [] */
        // zxchen
        // cout << "add record unsafe is called" << endl;
		auto it = vals.find(rec.data.first);  /* try find -- avoid ctor vcontainer */
		if (it != vals.end()) {
			it->second.add_v(rec.data.second, rec.ts);
		} else {
			/* no found, insert. it may fail, but we have valid iter anyway */
			auto ret = vals.insert(make_pair(rec.data.first, ValueContainerT{}));
			//  		auto ret = vals.emplace(rec.data.first, ValueContainerT{});
			/* ret.first: iterator. ret.first->second: vcontainer */
			ret.first->second.add_v(rec.data.second, rec.ts);
		}
	}

	void add_kv_unsafe(KVPair const & kv, ptime const & ts) {
		/* ver1 */
		//  	vals[kv.first].add_v(kv.second, ts);

		/* ver2: avoid using [] */
		auto ret = vals.insert(make_pair(kv.first, ValueContainerT{}));
		ret.first->second.add_v(kv.second, ts);
	}

	// go through each key (the min_ts is saved with the corresponding v
	// container)
	// should we save min_ts in this object?
	ptime min_ts() {
		ptime min_ts = max_date_time;
		for (auto & kv : vals) {
			if (kv.second.min_ts < min_ts)
				min_ts = kv.second.min_ts;
		}
		return min_ts;
	}
};

/////////////////////////////////////////////////////////////////////////////////////

/* All variants of WindowKeyedFragment
 * By default, it's safe, using tbb's map, and value container
 * It can be:
 * 	- Unsafe (using unsafe d/s)
 * 	- Simple (using simple vcontainer)
 * 	- MT (using internal internally partitioned HT)
 */

/* standard; most conservative */
template<class KVPair>
using WinKeyFrag_Std = WinKeyFrag<KVPair,
      tbb::concurrent_unordered_map<decltype(KVPair::first),
      ValueContainer<decltype(KVPair::second)>>>;
/* we defined hasher for fbstring */

template<class KVPair>
using WinKeyFragLocal_Std = WinKeyFragLocal<KVPair,
      std::unordered_map<decltype(KVPair::first),
      ValueContainerUnsafe<decltype(KVPair::second)>>>;

/* -- for netmon: many keys, few values -- */

template<class KVPair>
using WinKeyFragLocal_Simple = WinKeyFragLocal<KVPair,
      std::unordered_map<decltype(KVPair::first),
      SimpleValueContainerUnsafe<decltype(KVPair::second)>>>;

/* the HT internally partitioned */
template<class KVPair>
using WinKeyFrag_SimpleMT = WinKeyFrag<KVPair,
      creek::concurrent_unordered_map<CONFIG_NETMON_HT_PARTITIONS,
      decltype(KVPair::first),
      SimpleValueContainer<decltype(KVPair::second)>>>;


#endif /* WINDOWKEYEDFRAGMENT_H_ */
