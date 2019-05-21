
/* xzl. 
 * g++-5 -std=c++11 alloc.cpp
 * g++-4.8 will fail to create allocator. why?
 */

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <map>
#include <set>
#include <vector>
#include <string>
#include <unordered_map>

using namespace std;

/* 
 * http://www.drdobbs.com/the-standard-librarian-what-are-allocato/184403759?pgno=2 
 */

template <class T> class malloc_allocator
{
public:
  typedef T                 value_type;
  typedef value_type*       pointer;
  typedef const value_type* const_pointer;
  typedef value_type&       reference;
  typedef const value_type& const_reference;
  typedef std::size_t       size_type;
  typedef std::ptrdiff_t    difference_type;
  
  int node;  // xzl
  
  template <class U> 
  struct rebind { typedef malloc_allocator<U> other; };

  malloc_allocator(int node) : node (node) {}
  //malloc_allocator() {}
  malloc_allocator(const malloc_allocator&) {}
  template <class U> 
  malloc_allocator(const malloc_allocator<U>&) {}
  ~malloc_allocator() {}

  pointer address(reference x) const { return &x; }
  const_pointer address(const_reference x) const { 
    return x;
  }

  pointer allocate(size_type n, const_pointer = 0) {
    void* p = std::malloc(n * sizeof(T));
    if (!p)
      throw std::bad_alloc();
    return static_cast<pointer>(p);
  }

  void deallocate(pointer p, size_type) { std::free(p); }

  size_type max_size() const { 
    return static_cast<size_type>(-1) / sizeof(T);
  }

  void construct(pointer p, const value_type& x) { 
    new(p) value_type(x); 
  }
  void destroy(pointer p) { p->~value_type(); }

private:
  void operator=(const malloc_allocator&);
};

template<> class malloc_allocator<void>
{
  typedef void        value_type;
  typedef void*       pointer;
  typedef const void* const_pointer;

  template <class U> 
  struct rebind { typedef malloc_allocator<U> other; };
};


template <class T>
inline bool operator==(const malloc_allocator<T>&, 
                       const malloc_allocator<T>&) {
  return true;
}

template <class T>
inline bool operator!=(const malloc_allocator<T>&, 
                       const malloc_allocator<T>&) {
  return false;
}

int main(int argc, char *argv[]) {
  //malloc_allocator<int> alloc;
  
  malloc_allocator<int> alloc(0);
  
  vector<bool, malloc_allocator<bool>> a1(malloc_allocator<bool>());
  
  //set<int, std::less<int>, malloc_allocator<int>> a(0, 100, alloc);
  set<int, std::less<int>, malloc_allocator<int>> a(alloc);
}