#include "modelFunctions.hh"

namespace Kaskade
{
  namespace // anonymous namespace not visible outside translation unit
  {
    typedef Dune::FieldVector<double,1> FieldType;

    Dune::BlockVector<FieldType> operator*(Dune::Matrix<FieldType> const& A, Dune::BlockVector<FieldType> const& x)
    {
      Dune::BlockVector<FieldType> y(x.size());
      y *= 0;
      for(typename Dune::Matrix<FieldType>::size_type j=0; j<A.M(); ++j)
        for(typename Dune::Matrix<FieldType>::size_type i=0; i<A.N(); ++i)
          y[i] += A[i][j]*x[j];
      return y;
    }

    Dune::BlockVector<FieldType> operator*(double const& alpha, Dune::BlockVector<FieldType> const& x)
    {
      Dune::BlockVector<FieldType> y(x);
      y *= alpha;
      return y;
    }

    Dune::BlockVector<FieldType> operator+(Dune::BlockVector<FieldType> const& x, Dune::BlockVector<FieldType> const& y)
    {
      Dune::BlockVector<FieldType> z(x);
      z += y;
      return z;
    }

    std::vector<double> operator*(double a, std::vector<double> const& v)
    {
      std::vector<double> result(v);
      for(double& d : result) d *= a;
      return result;
    }

    std::vector<double> operator+(std::vector<double> const& v, std::vector<double> const& w)
    {
      assert(v.size() == w.size());
      std::vector<double> result(v);
      for(size_t i=0; i<v.size(); ++i) result[i] += w[i];
      return result;
    }
  }


  QuadraticFunction::QuadraticFunction(double constantPart_, Dune::BlockVector<FieldType> linearPart_, Dune::Matrix<FieldType> quadraticPart_)
    : constantPart(constantPart_), linearPart(linearPart_), quadraticPart(quadraticPart_)
  {}

  double QuadraticFunction::d0(std::vector<double> const& argument) const
  {
    Dune::BlockVector<FieldType> arg(argument.size());
    for(int i=0; i<argument.size(); ++i) arg[i]=argument[i];
    Dune::BlockVector<FieldType> tmp = quadraticPart*arg;
    return constantPart+linearPart*arg+arg*tmp;
  }

  std::vector<double> QuadraticFunction::d1(std::vector<double> const& argument) const
  {
    Dune::BlockVector<FieldType> arg(argument.size());
    for(int i=0; i<argument.size(); ++i) arg[i]=argument[i];
    auto tmp = linearPart+2.0*(quadraticPart*arg);
    std::vector<double> result(tmp.size(),0);
    for(size_t i=0; i<tmp.size(); ++i) result[i] = tmp[i];
    return result;
  }


  IsotropicCubicRegularization::IsotropicCubicRegularization(QuadraticFunction const& normBlf_, double scale_)
    : scale(scale_), normBlf(normBlf_)
  {}

  IsotropicCubicRegularization::~IsotropicCubicRegularization(){}

  double IsotropicCubicRegularization::d0(std::vector<double> const& iterate) const { return scale*pow(normBlf.d0(iterate),1.5); }

  std::vector<double> IsotropicCubicRegularization::d1(std::vector<double> const& iterate) const
  {
    return scale*1.5*pow(normBlf.d0(iterate),0.5)*normBlf.d1(iterate);
  }


  CSCubicRegularization::CSCubicRegularization(double normN_, QuadraticFunction const& norm, double scale_)
    : scale(scale_), normN3(normN_*normN_*normN_), normT(norm)
  {}

  CSCubicRegularization::~CSCubicRegularization(){}

  double CSCubicRegularization::d0(std::vector<double> const& iterate) const { return scale * ( normN3 + normT.d0(iterate) ); }

  std::vector<double> CSCubicRegularization::d1(std::vector<double> const& iterate) const { return scale * normT.d1(iterate); }


  RegularizedQuadraticFunction::RegularizedQuadraticFunction(QuadraticFunction quadraticModel_, DifferentiableScalarFunction& regularization_)
    : quadraticModel(quadraticModel_), regularization(regularization_)
  {}

  RegularizedQuadraticFunction::~RegularizedQuadraticFunction(){}

  double RegularizedQuadraticFunction::evalQuadraticModel(std::vector<double>const & iterate) const { return quadraticModel.d0(iterate); }

  double RegularizedQuadraticFunction::d0(std::vector<double> const & iterate) const { return quadraticModel.d0(iterate)+regularization.d0(iterate); }

  std::vector<double> RegularizedQuadraticFunction::d1(std::vector<double>const & iterate) const { return quadraticModel.d1(iterate)+regularization.d1(iterate); }


  CubicModel1dForFmin::CubicModel1dForFmin(ContinuousScalarFunction const& cb_) : cb(cb_) {}

  double CubicModel1dForFmin::operator()(double tau) const
  {
    std::vector<double> iterate(1);
    iterate[0]=tau;
    return cb.d0(iterate);
  }


  ContractionModelFunction::ContractionModelFunction(QuadraticFunction normBlf_, double omegaC_)
    : omegaC(omegaC_), normBlf(normBlf_)
  {}

  ContractionModelFunction::~ContractionModelFunction(){}

  double ContractionModelFunction::d0(std::vector<double>const & iterate) const { return omegaC/2*pow(normBlf.d0(iterate),0.5); }

  std::vector<double> ContractionModelFunction::d1(std::vector<double>const & iterate) const
  {
    return omegaC/4*pow(normBlf.d0(iterate),-0.5)*normBlf.d1(iterate);
  }
}
