/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2002-2011 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef EXTRAPOLATION_HH
#define EXTRAPOLATION_HH

#include <vector>

namespace Kaskade
{
  /**
   * \ingroup timestepping
   * \brief Polynomial extrapolation following Aitken-Neville.
   *
   * An Aitken-Neville extrapolation tableau for polynomial
   * interpolation/extrapolation.
   */
  template <class T>
  class ExtrapolationTableau
  {
  public:
    typedef T value_type;

    explicit ExtrapolationTableau(double hTarget_): hTarget(hTarget_) {}


    /**
     * Removes all previously inserted polynomial values.
     */
    void clear() { data.clear(); }

    /**
     * Returns the current extrapolation value, which takes all inserted
     * values into account. This is supposed to be the most accurate
     * one.
     */
    T const& back() const { return data.back(); }

    /**
     * Returns the order i extrapolation value based on the i+1 most
     * recently pushed values.
     */
    T const& operator[](int i) const {
      assert(0<=i && i<data.size());
      return data[i];
    }

    /**
     * Returns the number of available extrapolation values in the
     * current line of the tableau. This is one more than the maximal
     * order.
     */
    int size() const { return data.size(); }


    /**
     * Inserts polynomial value t(hnew) into the tableau.
     */
    void push_back(T t, double hnew)
    {
      T tmp(t);

      for (int i=0; i<data.size(); ++i) {
        // at loop start: t = t(n+1,i), data[i] = t(n,i)

        // compute coefficients, see Deuflhard/Hohmann 7.1.2
        int n = data.size()-1;
        double a = - (hTarget-hnew)/(hnew-h[n-i]);
        double b = (hTarget-h[n-i])/(hnew-h[n-i]);
        // compute tmp = t(n+1,i+1) = a*t(n,i)+b*t(n+1,i)
        tmp = t; tmp *= b;
        data[i] *= a; tmp += data[i];
        // store data[i] = t(n+1,i)
        data[i] = t;
        // maintain loop invariant: t = t(n+1,i+1)
        t = tmp;
      }
      // store new t(n+1,n+1)
      data.push_back(t);
      h.push_back(hnew);
    }



  private:
    // The array data always contains a complete row of the extrapolation tableau.
    std::vector<T>      data;
    std::vector<double> h;
    double              hTarget;
  };
} // namespace Kaskade

#endif
