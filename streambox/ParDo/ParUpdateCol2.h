#ifndef PARUPDATECOL2_H
#define PARUPDATECOL2_H

// apply a set of DoFns to each column.
//  -- (optional) modify individual columns (type won't change)
//  -- (optional) drop individual rows
// T_i: types of column i
template <class T1, class T2,     // input type of each column
//  class DoFn1 = DoFnSimple<T1>,
//  class DoFn2 = DoFnSimple<T2>,
//  class Tuple = tuple<T1, T2>,
  int N = 2>
class ParUpdateCol2 : public PTransform {

  using DoFn1 = DoFnSimple<T1>;
  using DoFn2 = DoFnSimple<T2>;
  using Tuple = tuple<T1, T2>;

private:

public:

  DoFn1 * fn1;
  DoFn2 * fn2;

  ParUpdateCol2(string name, DoFn1* fn1, DoFn2* fn2) :
    PTransform(name), fn1(fn1), fn2(fn2) { }

};

#endif /* PARUPDATECOL2_H */
