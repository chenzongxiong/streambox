#ifndef PARDOCOL2_H
#define PARDOCOL2_H

// if we produce any new column, all the old columns are dropped
template <typename T1, typename T2,     // input type of each column
  typename O1, typename O2,   // output type of each column
//  class DoFn1 = DoFn<T1, O1>,
//  class DoFn2 = DoFn<T2, O2>,
//  class InputElement = Element2<T1, T2>,
//  class OutputElement = Element2<O1, O2>,
  int N = 2>
class ParDoCol2 : public PTransform {

  using DoFn1 = DoFn<T1, O1>;
  using DoFn2 = DoFn<T2, O2>;
  using InputElement = Element2<T1, T2>;
  using OutputElement = Element2<O1, O2>;

public:
  string _name;
  DoFn1 * fn1;
  DoFn2 * fn2;

  ParDoCol2(string name, DoFn1* fn1, DoFn2* fn2) :
    _name(name), fn1(fn1), fn2(fn2) { }

};

#endif /* PARDOCOL2_H */
