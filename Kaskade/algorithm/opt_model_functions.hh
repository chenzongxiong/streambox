#ifndef OPT_MODEL_FUNCTIONS_HH
#define OPT_MODEL_FUNCTIONS_HH

#include <cmath>
#include <memory>

#include "dune/common/fvector.hh"
#include "dune/istl/bvector.hh"
#include "dune/istl/matrix.hh"

namespace Kaskade
{
  namespace OPT_MODEL_FUNCTIONS_Detail
  {
    template <class FieldType>
    inline Dune::BlockVector<FieldType> operator*(Dune::Matrix<FieldType> const& A, Dune::BlockVector<FieldType> const& x)
    {
      Dune::BlockVector<FieldType> y(x.size());
      y *= 0;
      for(typename Dune::Matrix<FieldType>::size_type j=0; j<A.M(); ++j)
        for(typename Dune::Matrix<FieldType>::size_type i=0; i<A.N(); ++i)
          y[i] += A[i][j]*x[j];
      return y;
    }

    template <class FieldType>
    inline Dune::BlockVector<FieldType> operator*(typename FieldType::value_type const& alpha, Dune::BlockVector<FieldType> const& x)
    {
      Dune::BlockVector<FieldType> y(x);
      y *= alpha;
      return y;
    }

    template <class FieldType>
    inline Dune::BlockVector<FieldType> operator+(Dune::BlockVector<FieldType> const& x, Dune::BlockVector<FieldType> const& y)
    {
      Dune::BlockVector<FieldType> z(x);
      z += y;
      return z;
    }
  }

  /// TODO !!!!! Entscheidung über Vektortyp !!!!!!
  class QuadraticFunction
  {
  public:
    typedef Dune::FieldVector<double,1> FieldType;

    QuadraticFunction(double constantPart_, Dune::BlockVector<FieldType> linearPart_, Dune::Matrix<FieldType> quadraticPart_)
    : constantPart(constantPart_), linearPart(linearPart_), quadraticPart(quadraticPart_)
    {}

    QuadraticFunction(QuadraticFunction const&) = default;

    double d0(std::vector<double> const& argument) const
    {
      using namespace OPT_MODEL_FUNCTIONS_Detail;
      Dune::BlockVector<FieldType> arg(argument.size());
      for(int i=0; i<argument.size(); ++i) arg[i]=argument[i];
      Dune::BlockVector<FieldType> tmp = quadraticPart*arg;
      return constantPart+linearPart*arg+arg*tmp;
    }

    Dune::BlockVector<FieldType> d1(std::vector<double> const& argument) const
    {
      using namespace OPT_MODEL_FUNCTIONS_Detail;
      Dune::BlockVector<FieldType> arg(argument.size());
      for(int i=0; i<argument.size(); ++i) arg[i]=argument[i];
      return linearPart+2.0*(quadraticPart*arg);
    }

  private:
    FieldType constantPart;
    Dune::BlockVector<FieldType> linearPart;
  public:
    Dune::Matrix<FieldType> quadraticPart;
  };

  class CubicModelFunction
  {
  public:
    typedef Dune::FieldVector<double,1> FieldType;

    CubicModelFunction(QuadraticFunction quadraticModel_, QuadraticFunction normBlf_, double omegaL_):
      omegaL(omegaL_), quadraticModel(quadraticModel_), normBlf(normBlf_), regWeight(0)
    {}

    CubicModelFunction(QuadraticFunction quadraticModel_, QuadraticFunction normBlf_,  double omegaL_, double regWeight_)
    : omegaL(omegaL_), quadraticModel(quadraticModel_), normBlf(normBlf_), regWeight(regWeight_)
    {}

    double evalQuadraticModel(std::vector<double>const & iterate) const
    {
      return quadraticModel.d0(iterate);
    }

    double d0(std::vector<double> const & iterate) const
    {
      return quadraticModel.d0(iterate)+omegaL/6*pow(normBlf.d0(iterate),1.5);
    }

    std::vector<double> d1(std::vector<double>const & iterate) const
    {
      using namespace OPT_MODEL_FUNCTIONS_Detail;
      Dune::BlockVector<FieldType> v = quadraticModel.d1(iterate)+omegaL/4*pow(normBlf.d0(iterate),0.5)*normBlf.d1(iterate);
      std::vector<double> w(v.size());
      for(int i=0; i<w.size(); ++i) w[i]=v[i][0];
      return w;
    }

  private:
    double omegaL;
    QuadraticFunction quadraticModel, normBlf;
    double regWeight;
  };

  class CubicModel1dForFmin
  {
  public:
    CubicModel1dForFmin(CubicModelFunction const& cb_) : cb(cb_) {}

    double operator()(double tau)
    {
      std::vector<double> iterate(1);
      iterate[0]=tau;
      return cb.d0(iterate);
    }

  private:
    CubicModelFunction const& cb;
  };

  class ContractionModelFunction
  {
  public:
    typedef Dune::FieldVector<double,1> FieldType;

    ContractionModelFunction(QuadraticFunction normBlf_, double omegaC_):
      omegaC(omegaC_), normBlf(normBlf_)
    {}

    double d0(std::vector<double>const & iterate) const
    {

      return omegaC/2*pow(normBlf.d0(iterate),0.5);
    }

    std::vector<double> d1(std::vector<double>const & iterate) const
    {
      using namespace OPT_MODEL_FUNCTIONS_Detail;
      Dune::BlockVector<FieldType> v =  omegaC/4*pow(normBlf.d0(iterate),-0.5)*normBlf.d1(iterate);
      std::vector<double> w(v.size());
      for(int i=0; i<w.size(); ++i) w[i]=v[i][0];
      return w;
    }

  private:
    double omegaC;
    QuadraticFunction normBlf;
  };

}  // namespace Kaskade
#endif


