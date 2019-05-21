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

#ifndef NORMS_HH
#define NORMS_HH

/**
 * @file
 * @brief  Some norms for FunctionSpaceElement s or FunctionViews
 * @author Anton Schiela
 */

#include <boost/fusion/include/at_c.hpp>
#include <fem/integration.hh>
#include <fem/functionviews.hh>

namespace Kaskade
{
/// L_2-norms
struct L2Norm
{
/// Evaluation of square norm
  template<typename Function> 
  double square(Function const& f) const
  {
    Integral<typename Function::Space> integral;
    FunctionViews::AbsSquare<Function> asf(f);
    return integral(asf);
  }

/// Evaluation of norm
  template<typename Function> 
  double operator()(Function const& f) const
  {
    return std::sqrt(square(f));
  }
};

/// L_2-norms
// struct L2ScalarProduct
// {
// /// Evaluation of square norm
//   template<typename Function> 
//   double operator()(Function f,Function g) 
//   {
//     Integral<typename Function::Space> integral;
//     FunctionViews::Dot<Function> dotfg(f,g);
//     return integral(dotfg);
//   }
// 
// /// Evaluation of norm
//   template<typename Function> 
//   double operator()(Function f) 
//   {
//     return std::sqrt(this->operator()(f,f));
//   }
// };

/// H1-semi-norms
struct H1SemiNorm
{
/// Evaluation of square norm
  template<typename Function> 
  double square(Function const& f) const
  {
    Integral<typename Function::Space> integral;
    FunctionViews::GradientAbsSquare<Function> asf(f);
    return integral(asf);
//    FunctionViews::Gradient<Function> gf(f);
//    L2Norm l2Norm;
//    return l2Norm.square(gf);
  }  

/// Evaluation of norm
  template<typename Function> 
  double operator()(Function const& f) const
  {
    return std::sqrt(square(f));
  }
};


// THE FOLLOWING CODE IS UNDOCUMENTED, APPEARS TO BE RATHER PROJECT-SPECIFIC, AND IS NOWHERE USED IN
// THE SVN REPOSITORY. UNLESS SOME UNEXPECTED PROBLEMS ARISE, THIS CODE WILL BE DELETED AS OF 
// 2016-03-31.
// 
// template <int val>
// struct Decrement
// {
//   static constexpr int value = val-1;
// };
// 
// 
// 
// template <int variableId, int numberOfVariables>
// struct MyNormEvaluator
// {
//   template <class DomainElement, class Norm>
//   static double square(DomainElement const& x, Norm const& norm)
//   {
//     return norm.square( boost::fusion::at_c<variableId>(x.data) ) + MyNormEvaluator<variableId-1,numberOfVariables>::square(x,norm);
//   }
// };
// 
// template <int numberOfVariables>
// struct MyNormEvaluator<0,numberOfVariables>
// {
//   template <class DomainElement, class Norm>
//   static double square(DomainElement const& x, Norm const& norm)
//   {
//     return norm.square( boost::fusion::at_c<0>(x.data) );
//   }
// };
// 
// struct MyH1SemiNorm
// {
//   template <typename DomainElement>
//   static double square(DomainElement const& x)
//   {
//     return MyNormEvaluator<DomainElement::Descriptions::noOfVariables-1, DomainElement::Descriptions::noOfVariables>::square(x,H1SemiNorm());
//   }
// 
//   template <typename DomainElement>
//   double operator()(DomainElement const& x) const
//   {
//     return sqrt(square(x));
//   }
// };

// struct MyL2Norm
// {
//   template <typename DomainElement>
//   static double square(DomainElement const& x)
//   {
//     return MyNormEvaluator<DomainElement::Descriptions::noOfVariables-1, DomainElement::Descriptions::noOfVariables>::square(x,L2Norm());
//   }
// 
//   template <typename DomainElement>
//   double operator()(DomainElement const& x) const
//   {
//     return sqrt(square(x));
//   }
// };


/// H1-norms
struct H1Norm
{
/// Evaluation of square norm
  template<typename Function> 
  double square(Function f) 
  {
    H1SemiNorm h1Norm;
    L2Norm l2Norm;
    return h1Norm.square(f)+l2Norm.square(f);
  }  

/// Evaluation of norm
  template<typename Function> 
  double operator()(Function f) 
  {
    return std::sqrt(square(f));
  }
};

template<class Space> class LocalIntegral;
template<class Grid, class T> class CellData;

/// local (cellwise) H1-semi-norms
template<class Function>
typename CellData<typename Function::Space::Grid, typename Function::ValueType>::CellDataVector localH1SemiNorm(Function const& f)
{
  typedef typename Function::Space Space;
  typedef typename Space::Grid Grid;
  LocalIntegral<Space> localIntegral;
  typename CellData<Grid,typename Function::ValueType>::CellDataVector
    errorIndicator(localIntegral(
                     makeView<FunctionViews::AbsSquare>(makeView<FunctionViews::Gradient>(f))));
  return errorIndicator;
}

/// local (cellwise) L2-norms
template<class Function>
typename CellData<typename Function::Space::Grid, typename Function::ValueType>::CellDataVector localL2Norm(Function const& f)
{
  typedef typename Function::Space Space;
  typedef typename Space::Grid Grid;
  LocalIntegral<Space> localIntegral;
  typename CellData<Grid,typename Function::ValueType>::CellDataVector
    errorIndicator(localIntegral(
                     makeView<FunctionViews::AbsSquare>(f)));
  return errorIndicator;
}


/**
 * @brief boundaryL2Norm computes the L2-norm of an FE function on the whole boundary of the underlying grid.
 *
 * @param function is the FE function to be integrated.
 * @tparam FEFunction is the type of the integrand.
 * @return value of integral
 *
 */
template <class FEFunction>
typename FEFunction::Space::Scalar boundaryL2Norm(FEFunction const& function) {
  FunctionViews::AbsSquare<FEFunction> functionSquared(function);
  return std::sqrt(integrateOverBoundary(functionSquared));
}
} // end of namespace Kaskade
#endif
