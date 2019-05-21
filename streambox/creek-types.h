#ifndef CREEK_TYPES_H_
#define CREEK_TYPES_H_
#ifdef USE_FOLLY_STRING
#include "folly/FBString.h"
namespace creek {
	using string = folly::fbstring;
}

/* hasher for fbstring. inject into the tbb namespace, must be done before including concurrent ds. otherwise
 * they will fall back to default tbb hasher. */
namespace tbb {
	size_t tbb_hasher(const folly::fbstring & key);
}

/* for dbg -- ideally should provide tbb_hasher() */
//	struct FBStringHash {
//		size_t operator()(const folly::fbstring & key)  const {
//			return tbb::tbb_hasher(key);
//		}
//	};
#else
namespace creek {
	using string = std::string;
}
#endif

namespace creek {
#ifdef USE_NUMA_ALLOC
#include "utilities/threading.hh"
	template<typename T>
		using allocator = Kaskade::NumaAllocator<T>;
#else
	template<typename T>
		using allocator = std::allocator<T>;
#endif
}

namespace creek {
	//	using ippair = uint64_t;
	using ippair = long;
	using tvpair = std::pair<long, long>;
}

#include "creek-map.h"

/* -- concurrent vector -- */
#include "tbb/concurrent_unordered_map.h"
#include "tbb/concurrent_vector.h"
#include "tbb/concurrent_unordered_set.h"

namespace creek {
	template <class T>
		using concurrent_vector = tbb::concurrent_vector<T>;

	template <class T>
		using concurrent_vector_ptr = shared_ptr<creek::concurrent_vector<T>>;

	template <class T>
		using concurrent_unordered_set = tbb::concurrent_vector<T>;
	//    using concurrent_vector_ptr = shared_ptr<creek::concurrent_vector>;
}

#endif /* CREEK_TYPES_H_ */
