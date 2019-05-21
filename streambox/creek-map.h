#ifndef CREEK_MAP_H_
#define CREEK_MAP_H_

#include <atomic>
#include <memory>

#include "tbb/concurrent_unordered_map.h"
//#include "folly/AtomicHashMap.h"
#include "log.h"

using namespace std;

/* a wrapper over multi concurrent maps (to over scalability bottleneck?) */
namespace creek {

//#define MAX_MAPS 1024

//#ifdef USE_TBB_HASHMAP
#if 1
/* wraps around tbb */
template <int MAX_MAPS, typename Key, typename T, typename Hasher = tbb::tbb_hash<Key>, typename Key_equality = std::equal_to<Key>,
         typename Allocator = tbb::tbb_allocator<std::pair<const Key, T> > >
	struct concurrent_unordered_map {

		using Map = tbb::concurrent_unordered_map<Key,T,Hasher,Key_equality,Allocator>;
//		using Map = folly::AtomicHashMap<Key,T>;

    // Type definitions
    typedef Key key_type;
//    typedef typename base_type::value_type value_type;
    typedef std::pair<const Key, T> value_type;
    typedef T mapped_type;
    typedef Hasher hasher;
    typedef Key_equality key_equal;
//    typedef hash_compare key_compare;

//    typedef typename base_type::allocator_type allocator_type;
//    typedef typename base_type::pointer pointer;
//    typedef typename base_type::const_pointer const_pointer;
//    typedef typename base_type::reference reference;
//    typedef typename base_type::const_reference const_reference;
//
//    typedef typename base_type::size_type size_type;
//    typedef typename base_type::difference_type difference_type;
//
//    typedef typename base_type::iterator iterator;
//    typedef typename base_type::const_iterator const_iterator;
//    typedef typename base_type::iterator local_iterator;
//    typedef typename base_type::const_iterator const_local_iterator;

    using iterator = typename Map::iterator;
    using const_iterator = typename Map::const_iterator;
    using size_type = typename Map::size_type;

		std::array<Map, MAX_MAPS> maps;
		std::atomic<int> num_maps;

		concurrent_unordered_map(const int size = MAX_MAPS) : num_maps(size) {
			xzl_bug_on(size > MAX_MAPS || size <= 0);
		}

		/* key to index in the map array */
		inline int k2i(const key_type& k) {
			return (tbb::tbb_hasher(k) % num_maps);
		}

		/* many HT does not implement this */
    mapped_type& operator[](const key_type& k) {
    	return maps[k2i(k)][k];
    }

    mapped_type& at( const key_type& k ) {
    	return maps[k2i(k)].at(k);
    }

    /* this is unsafe */
    size_type size() const {
    	size_type s = 0;
    	for (auto & m : maps)
    		s += m.size();
    }

    /* the caller may compare returned iterator with end(), so
     * we need to provide one.
     * Here we just use the last map's end() as end().
     */
    iterator begin() {
        return maps[0].begin();
    }
    iterator end() {
    	return maps[num_maps - 1].end();
    }

    /* if no found, return the global end() */
    iterator find(const key_type& k) {
    	auto && map = maps[k2i(k)];
			auto it = map.find(k);
			if (it == map.end())
				return this->end(); /* the global end */
			else
				return it;
    }

		Map & get_map(int i) {
			xzl_bug_on(i >= num_maps);
			return maps[i];
		}

    /* not very useful, since other maps don't have this interface.
     * also return whether find succeeds, since caller may
     * not see individual map
     */
    std::pair<iterator, bool> myfind(const key_type& k) {
    	auto && map = maps[k2i(k)];
    	auto it = map.find(k);
    	return make_pair(it, it != map.end());
    }

    // modifiers
    std::pair<iterator, bool> insert(const value_type& x) {
    	auto && map = maps[k2i(x.first)];
    	return map.insert(x);
    }

    template<typename... Args> std::pair<iterator, bool> emplace(Args&&... args) {

    }

    // modifiers
		iterator insert(const_iterator hint, const value_type& x);
		template<class InputIterator>
		void insert(InputIterator first, InputIterator last);

	};
#endif

//#elif defined(USE_FOLLY_HASHMAP)

#if 0
template <typename Key, typename T>
	struct concurrent_unordered_map {

		using Map = folly::AtomicHashMap<Key,T>;

    // Type definitions
    typedef Key key_type;
    typedef std::pair<const Key, T> value_type;
    typedef T mapped_type;
//    typedef Hasher hasher;
//    typedef Key_equality key_equal;

    using iterator = typename Map::iterator;
    using const_iterator = typename Map::const_iterator;
    using size_type = typename Map::size_type;

//		std::vector<Map> maps;  // does not work since stl requires copy ctor
    std::vector<shared_ptr<Map>> pmaps;
		std::atomic<int> num_maps;

		concurrent_unordered_map(const int size = MAX_MAPS)
			: num_maps(size) {
			xzl_bug_on(size > MAX_MAPS || size <= 0);
			for (int i = 0; i < size; i++)
				pmaps.push_back(make_shared<Map>(1000 * 1000));
		}

		/* key to index in the map array */
		inline int k2i(const key_type& k) {
			return (tbb::tbb_hasher(k) % num_maps);
		}

		Map & get_map(int i) {
			xzl_bug_on(i >= num_maps);
			return *pmaps[i];
		}

    /* this is unsafe */
    size_type size() const {
    	size_type s = 0;
    	for (auto & m : pmaps)
    		s += m->size();
    }

    /* the caller may compare returned iterator with end(), so
     * we need to provide one.
     * Here we just use the last map's end() as end().
     */
    iterator end() {
    	return pmaps[num_maps - 1]->end();
    }

    /* if no found, return the global end() */
    iterator find(const key_type& k) {
    	auto && map = pmaps[k2i(k)];
			auto it = map->find(k);
			if (it == map->end())
				return this->end(); /* the global end */
			else
				return it;
    }

    // modifiers
    std::pair<iterator, bool> insert(const value_type& x) {
    	auto && map = pmaps[k2i(x.first)];
    	return map->insert(x);
    }
	};
#endif // if 1

//#endif
};  // namespace

#endif /* CREEK_MAP_H_ */
