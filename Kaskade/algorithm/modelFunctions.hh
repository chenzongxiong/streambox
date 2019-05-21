#ifndef MODELFUNCTIONS_HH
#define MODELFUNCTIONS_HH

#include <cmath>
#include <memory>

#include "dune/common/fvector.hh"
#include "dune/istl/bvector.hh"
#include "dune/istl/matrix.hh"

#include "algorithm/abstract_interface.hh"

namespace Kaskade
{
  class QuadraticFunction : public DifferentiableScalarFunction
  {
    typedef Dune::FieldVector<double,1> FieldType;
  public:
    QuadraticFunction(double constantPart_, Dune::BlockVector<FieldType> linearPart_, Dune::Matrix<FieldType> quadraticPart_);

    QuadraticFunction(QuadraticFunction const&) = default;

    double d0(std::vector<double> const& argument) const;

    std::vector<double> d1(std::vector<double> const& argument) const;

  private:
    FieldType constantPart;
    Dune::BlockVector<FieldType> linearPart;
  public:
    Dune::Matrix<FieldType> quadraticPart;
  };

  class IsotropicCubicRegularization : public DifferentiableScalarFunction
  {
  public:
    IsotropicCubicRegularization(QuadraticFunction const& normBlf_, double scale_ = 1);
    virtual ~IsotropicCubicRegularization();

    double d0(std::vector<double> const& iterate) const;
    std::vector<double> d1(std::vector<double> const& iterate) const;

  private:
    double scale;
    QuadraticFunction normBlf;
  };


  class CSCubicRegularization : public DifferentiableScalarFunction
  {
  public:
    CSCubicRegularization(double normN_, QuadraticFunction const& norm, double scale_ = 1);
    virtual ~CSCubicRegularization();

    double d0(std::vector<double> const& iterate) const;
    std::vector<double> d1(std::vector<double> const& iterate) const;

  private:
    double scale, normN3;
    QuadraticFunction normT;
  };


  class RegularizedQuadraticFunction : public DifferentiableScalarFunction
  {
  public:
    RegularizedQuadraticFunction(QuadraticFunction quadraticModel_, DifferentiableScalarFunction& regularization_);

    virtual ~RegularizedQuadraticFunction();

    double evalQuadraticModel(std::vector<double>const & iterate) const;

    double d0(std::vector<double> const& iterate) const;
    std::vector<double> d1(std::vector<double> const& iterate) const;

  private:
    QuadraticFunction quadraticModel;
    DifferentiableScalarFunction const& regularization;
  };


  class CubicModel1dForFmin
  {
  public:
    CubicModel1dForFmin(ContinuousScalarFunction const& cb_);

    double operator()(double tau) const;

  private:
    ContinuousScalarFunction const& cb;
  };

  class ContractionModelFunction : public DifferentiableScalarFunction
  {
  public:
    ContractionModelFunction(QuadraticFunction normBlf_, double omegaC_);

    virtual ~ContractionModelFunction();

    double d0(std::vector<double>const & iterate) const;
    std::vector<double> d1(std::vector<double>const & iterate) const;

  private:
    double omegaC;
    QuadraticFunction normBlf;
  };

}  // namespace Kaskade

#endif // MODELFUNCTIONS_HH
