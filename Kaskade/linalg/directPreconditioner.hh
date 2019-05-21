/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2002-2016 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef DIRECT_PRECONDITIONER_HH
#define DIRECT_PRECONDITIONER_HH

#include <memory>

#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvercategory.hh>

#include "linalg/crsutil.hh"
#include "linalg/matrixTraits.hh"
#include "linalg/symmetricOperators.hh"
#include "utilities/duneInterface.hh"

namespace Kaskade
{
  /// \internal
  // forward declarations
  template <class,class> class Factorization;
  template <class,class> class NumaBCRSMatrix;
  /// \endinternal

  /**
   * \ingroup direct
   * \brief A direct preconditioner applicable to assembled
   * operators.
   *
   * \tparam Op the assembled operator. Note that Op::matrix_type must be MatrixAsTriplet<> and
   *         Op::Domain must be LinearProductSpace<>. This is satisfied for assembled Galerkin
   *         operators coming from a VariationalFunctionalAssembler.
   */
  template <class Op>
  class DirectPreconditioner: public Dune::Preconditioner<typename Op::Domain, typename Op::Range>
  {
    typedef typename Op::Domain Domain;
    typedef typename Op::Domain Range;
    typedef typename GetScalar<Domain>::type Scalar;

  public:
    static int const category = Dune::SolverCategory::sequential;

    /**
     * \arg op the assembled operator
     * \arg first,last a half-open index range denoting the range where the direct solver will be applied.
     * \arg directType the type of direct solver to be used
     * \arg property the properties of the matrix (MatrixProperties::SYMMETRIC,MatrixProperties::GENERAL)
     */
    DirectPreconditioner(Op const& op, DirectType directType=DirectType::UMFPACK, MatrixProperties property=MatrixProperties::SYMMETRIC)
    : factorization(getFactorization(directType,property,op.template get<MatrixAsTriplet<Scalar> >()))
    {}

    void pre(Domain&, Range&) {}
    void post (Domain&) {}

    void apply (Domain& x, Range const& y)
    {
      std::vector<Scalar> b(y.dim());
      vectorToSequence(y,begin(b));
      factorization->solve(b);
      vectorFromSequence(x,begin(b));
    }

  private:
    std::unique_ptr<Factorization<Scalar>> factorization;
  };
  


} // namespace Kaskade
#endif
