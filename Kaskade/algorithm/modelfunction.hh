#ifndef MODELFUNCTION_HH
#define MODELFUNCTION_HH

#include "abstract_interface.hh"
#include "algorithm_base.hh"
#include "istl/matrix.hh"
#include "linalg/simpleLAPmatrix.hh"
#include <vector>
#include "dune/common/fmatrix.hh"

/**
 * @file 
 * @brief Some model functions for various purposes
 * @author Anton Schiela
 *
 */

namespace Kaskade
{

/// Interface for a scalar model function on a scalar domain that may fit data
class ScalarModelFunction
{
public:
  ScalarModelFunction() : data("SamplingData"), modified(false), datastart(0) {}

/// Add new value to the sample points, will be overwritten by the next update, until fixUpdate() is called
virtual void update(double abscissa, double value) { data = std::make_pair(abscissa,value); modified=true; }
  
  /// Compute value of the model function at abscissa
  virtual double getValue(double abscissa) = 0;

  /// Compute integral of the model function from low to high
  virtual double getIntegral(double high, double low) { assert(!"Computation of integrals not implemented"); return 0.0; }

  /// Compute abscissa corresponding to a value. This is of course only possible for monotone functions
  virtual double getAbscissa(double value) { assert(!"Computation of abscissa not implemented"); return 0.0; }
  
  /// fix the last updated value and thus increase the number of sampling points
  virtual void fixUpdate() { if(data.isValid()) data.logValue(); }
  
  virtual ~ScalarModelFunction() {};

  void setDataStart(int ds) {datastart=ds;};
protected:  
  virtual double getDataVal(int i) { return data[i].second;}
  virtual double getDataAbsc(int i) { return data[i].first;}
  virtual double getScaling(int i) { return 1.0;}
  

  LoggedQuantity<std::pair<double, double> > data;
  bool modified;
  int datastart;
};

class PolynomialModel : public ScalarModelFunction
{
public:
  PolynomialModel(std::vector<double> defaults)
    : ScalarModelFunction(), degree(defaults.size()-1), coefficients(defaults), maxels(1000000000)
      {}

  PolynomialModel(int degree_, int maxels_=10000000)
    : ScalarModelFunction(), degree(degree_), coefficients(degree+1,0.0), maxels(maxels_)
  {}

  virtual double getValue(double abscissa)
  {
    computeCoefficients();
    double value=0.0;
    double absc=1.0;
    for(int i=0; i<=degree; ++i)
    {
      value +=absc*coefficients[i];
      absc*=abscissa;
    } 
    return value;
  }

  double getAbscissa(double value)
  {
    assert(degree == 1);
    return (value-coefficients[0])/coefficients[1];
  }

  std::vector<double> const& getCoefficients() { return coefficients;}

protected:

  virtual void computeCoefficients()
  {
    if(modified)
    {
      int enddata=data.size();
      int rows=std::min(enddata-datastart,maxels);
      int startdata=enddata-rows;
      int cols = std::min(degree+1,rows);
      if(rows==0 || cols==0) return; 
      double absc;
      Dune::Matrix<Dune::FieldMatrix<double,1,1> > A(rows,cols);
      std::vector<double> b(rows), x(cols);
      for(int i=0; i<rows; ++i)
      {
        absc=1.0;
        b[i]=getDataVal(startdata+i)*getScaling(startdata+i);
        for(int j=0; j<cols; ++j)
        {
          A[i][j]=absc*getScaling(startdata+i);
          absc*=getDataAbsc(startdata+i);
        }
      }
      LeastSquares(SLAPMatrix<double>(A), b, x);
      for(int i=0; i<x.size(); ++i) coefficients[i]=x[i];
      modified=false;
    }
  }

  int degree;
  std::vector<double> coefficients;
  int maxels;
};

class TransformedPolynomialModel : public PolynomialModel
{
public:
  TransformedPolynomialModel(std::vector<double> defaults)
   : PolynomialModel(defaults)
      {}

  TransformedPolynomialModel(int degree_)
    : PolynomialModel(degree_)
  {}


  virtual double getValue(double abscissa) { return valBwd(PolynomialModel::getValue(abscFwd(abscissa))); }
  virtual double getAbscissa(double value) { return abscBwd(PolynomialModel::getAbscissa(valFwd(value))); }
protected:
  virtual double abscFwd(double a) {return a;}
  virtual double abscBwd(double a) {return a;}
  virtual double valFwd(double v)  {return v;}
  virtual double valBwd(double v)  {return v;}
  virtual double getDataVal(int i) { return valFwd(data[i].second);}
  virtual double getDataAbsc(int i) { return abscFwd(data[i].first);}

};

/// Polynomial model in a loglog scale
class LogLogModel : public TransformedPolynomialModel
{ 
public:
  LogLogModel(std::vector<double> defaults)
    : TransformedPolynomialModel(defaults)
  {}
  LogLogModel(int degree_)
    : TransformedPolynomialModel(degree_)
  {}

  virtual double getIntegral(double high, double low)
  {
    computeCoefficients();
    assert(degree == 1);
    double a=coefficients[1];
    double b=std::exp(coefficients[0]);
    assert(a>-1.0 || high*low > 0);
    return b/(1+a)*(std::pow(high,1+a)-std::pow(low,1+a));
  }

protected:
  virtual double abscFwd(double a) { assert(a>0); return std::log(a); }
  virtual double abscBwd(double a)  { return std::exp(a); }
  virtual double valFwd(double v)  { assert(v>0); return std::log(v); }
  virtual double valBwd(double v)  { return std::exp(v); }
};

class OmegaModel : public LogLogModel
{
public:
  OmegaModel() : LogLogModel(std::vector<double>(2,-2.0)) {};
protected:
  virtual void computeCoefficients()
  {
    if(modified)
    {
      LogLogModel::computeCoefficients();
      assert(coefficients[0] > 0);
      coefficients[1]=std::min(coefficients[1],0.0);
    }
  }
};

class EtaModel : public LogLogModel
{
public:
  EtaModel() : LogLogModel(std::vector<double>(2,-1.0)) {};
protected:
  virtual void computeCoefficients()
  {
    if(modified)
    {
      LogLogModel::computeCoefficients();
      assert(coefficients[0] > 0);
      coefficients[1]=std::min(coefficients[1],0.0);
    }
  }
};

class JModel : public LogLogModel
{
public:
  JModel() : LogLogModel(std::vector<double>(3,1.0)) { setDataStart(1); };
  virtual void update(double abscissa, double value) { LogLogModel::update(abscissa,value); }

protected:
  virtual double getDataVal(int i) { return valFwd(fabs(data[i-1].second-data.value().second));}
  virtual double getDataAbsc(int i) { return abscFwd(data[i-1].first-data.value().first);}
  virtual double getScaling(int i) 
  {
    if(data[i-1].second-data.value().second > 0 && i < data.size()-2)
      return abscBwd(getDataAbsc(data.size()-1))/abscBwd(getDataAbsc(i));
    else return 0.0;
  }


  virtual void computeCoefficients()
  {
    if(modified)
    {
      LogLogModel::computeCoefficients();
    }
  }
};

 class JModelLin : public PolynomialModel
 {
 public:
   JModelLin() : PolynomialModel(1,2) {};
 };

} // namespace Kaskade
#endif
