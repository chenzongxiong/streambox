#ifndef SELECT_HH
#define SELECT_HH

/// Select a type on compile time.
/**
 * if
 * decide=true: Select::type=T1
 * decide=false: Select::type=T2
 */

template <bool decide, class T1, class T2>
struct Select{
  typedef T1 type;
};

template <class T1, class T2>
struct Select<false,T1,T2>{
  typedef T2 type;
};

#endif
