#ifndef TYPE_2_TYPE_HH
#define TYPE_2_TYPE_HH

/**
 * Type2Type template for the efficient emulation of template member function
 * specialication using overloading (see Alexandrescu: Modern C++ Design)
 * Maybe its worth thinking about incorporating Alexandrescu's Loki
 * in order to use more tools like this
 */
template <typename T>
struct Type2Type
{
  typedef T Result;
};

#endif
