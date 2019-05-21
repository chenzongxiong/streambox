#ifndef DISCRETE_BRIDGE_HH
#define DISCRETE_BRIDGE_HH

#include <memory>
#include "newton_bridge.hh"
#include "dune/common/fmatrix.hh"
#include "dune/istl/bvector.hh"


namespace Bridge
{

//--------------------------------------------------------------------------------
/// Linearization Implementation for a fixed system, i.e. an inherently finite dimensional system, may be a template parameter Functional
  /** Functional has to implement:
   * typename RT : real type
   * RT d0()
   * d1(std::vector<RT>* rhs)
   * d2(Dune::Matrix<Dune::FieldMatrix<RT, 1, 1> > m, Dune::BlockVector<Dune::FieldVector<RT,1> >* x)
   * int size()
   */
template<class Functional>
class FixedSystemLinearization
{
public:
  typedef typename Functional::RT RT;
  typedef typename Functional::Vector DomainElement;
  typedef typename Functional::Vector ImageElement;
  typedef Functional Implementation;
  class Empty{ public: typedef ImageElement range_type; typedef DomainElement domain_type; }; typedef Empty OperatorType;

//  typedef Functional Linearization;

  FixedSystemLinearization(Functional& fu_, AbstractFunctionSpaceElement const& x_) 
    : x(x_.clone()), fu(fu_), mat(rows(0,nColBlocks()),cols(0,nColBlocks())), mflushed(true)
  {
    fu.setOrigin(getImpl<DomainElement>(*x));
    rhs.resize(0);
  }

  virtual ~FixedSystemLinearization() {};

  int cols(int cbegin, int cend) const { if(cbegin<cend) return std::min(cend,nColBlocks())-std::max(cbegin,0); return 0;}
  int rows(int rbegin, int rend) const { if(rbegin<rend) return std::min(rend,nRowBlocks())-std::max(rbegin,0); return 0;}

  void precompute() 
  {
    createLinearizedSystem();
  }

  virtual void getMatrixBlocks(MatrixAsTriplet<RT>& mat_, int rbegin, int rend, int cbegin, int cend) const 
  { 
    createLinearizedSystem();
    for(int i=rbegin; i<rend; i++)
      for(int j=cbegin; j<cend; j++)
        mat_.addEntry(i-rbegin,j-cbegin,mat[i][j]);
  }
  virtual void getRHSBlocks(std::vector<RT>& rhs_, int rbegin, int rend) const 
  { 
    createRHS();
    rhs_.resize(rend-rbegin);
    for(int i=0; i<rend-rbegin;++i)
      rhs_[i] = rhs[i+rbegin];
  }

  virtual int nColBlocks() const { return fu.size();}

  virtual int nRowBlocks() const { return fu.size();}

  double getValue() const { return fu.d0();}

  AbstractFunctionSpaceElement const& getOrigin() const {return *x;}

  void flush() { mflushed=true; rhs.resize(0); mat*=0.0;}

/// return the implementation
Implementation const& getLinImpl() const {return fu; }

private:
  void createLinearizedSystem() const
  {
    createRHS();
    if(!mflushed)
    {
      return;
    }
    fu.d2(mat);
    mflushed=false;
  } 

  void createRHS() const
  {
    rhs.resize(cols(0,nColBlocks()));
    fu.d1(rhs);
  } 

  mutable std::vector<RT> rhs;
  std::std::unique_ptr<AbstractFunctionSpaceElement > x;
  Functional& fu;
  mutable typename Functional::Matrix mat;
  mutable bool mflushed;
};

  template<class T, class Functional>
  class LinearizationTraits<Dune::BlockVector<T>, Functional>
  {
  public:
    typedef FixedSystemLinearization<Functional> Linearization;
  };

}


#endif
