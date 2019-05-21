#include "abstract_interface.hh"
#include "lagrangeLinearization.hh"

namespace Kaskade
{
  LagrangeLinearization::LagrangeLinearization(AbstractLinearization* Nlin_, AbstractLinearization* Tlin_, AbstractFunctionSpaceElement const& adjointCorrection, int stateId_, int controlId_, bool hasTlin_uu_)
    : Nlin(Nlin_), Tlin(Tlin_), origin(Tlin->getOrigin().clone()), stateId(stateId_), controlId(controlId_), hasTlin_uu(hasTlin_uu_)
  {
    lagrangeCorrection = origin->clone();
    *lagrangeCorrection *= 0.0;
//      existsLagrangian=false;
//      for(int i=0; i<origin.nComponents();++i)
//        if(origin.getRole(i) == "dual") existsLagrangian = true;

    for(int i=0; i<origin->nComponents(); ++i)
      for(int k=0; k<origin->nComponents(); ++k)
        if(origin->getRole(i)!="dual" && origin->getRole(k) == "dual")
          Nlin->ddxpy(*lagrangeCorrection, *origin/*adjointCorrection*/,i,i+1,k,k+1);
  }

 LagrangeLinearization::~LagrangeLinearization() {}

  /// Evaluate f(origin)
 double LagrangeLinearization::eval() const { return Nlin->eval(); }

 double LagrangeLinearization::evalL1norm() const { return Nlin->eval(); }

  /// Evaluate f'(origin)(.) + <Lag. Multiplier, c'(origin)(.)  > this is dual element
 void LagrangeLinearization::evald(AbstractFunctionSpaceElement &g, int rbegin, int rend) const
  {
    // Right hand side is kept in Nlin, but possibly not in Tlin
    if(rend==-1) rend=g.nComponents();
    Nlin->evald(g,rbegin,rend);
  }

  void LagrangeLinearization::getMatrixBlocks(MatrixAsTriplet<double> &mat, int rbegin, int rend, int cbegin, int cend) const
  {
    if(rend==-1) rend = getOrigin().nComponents();
    if(cend==-1) cend = getOrigin().nComponents();
    MatrixAsTriplet<double> tmp;
    size_t rowOffset = 0, colOffset = 0;
    for(int i=rbegin; i<rend; ++i)
    {
      if(i==rbegin) rowOffset = 0;

      for(int j=cbegin; j<cend; ++j)
      {
        if(j==cbegin) colOffset = 0;

        if(origin->getRole(i) != "dual" && origin->getRole(j) != "dual" && !(hasTlin_uu && i==controlId && j==controlId) )
        {
          if( i >= j) Tlin->getMatrixBlocks(tmp,i,i+1,j,j+1);
          else
          {
            Tlin->getMatrixBlocks(tmp,j,j+1,i,i+1);
            tmp.transpose();
          }

          tmp.shiftIndices(rowOffset,colOffset);
          mat += tmp;
        }
        else
        {
          if( i >= j ) Nlin->getMatrixBlocks(tmp,i,i+1,j,j+1);
          else
          {
            Nlin->getMatrixBlocks(tmp,j,j+1,i,i+1);
            tmp.transpose();
          }
          tmp.shiftIndices(rowOffset,colOffset);
          mat += tmp;
        }

        colOffset += tmp.ncols();
      }

      rowOffset += tmp.nrows();
    }
  }

   /// Evaluate hessian (of Lagrangian) times second argument
  void LagrangeLinearization::d2axpy(double a, AbstractFunctionSpaceElement& y, AbstractFunctionSpaceElement const& x, int rbegin, int rend, int cbegin, int cend) const
   {
     if(rend==-1) rend = y.nComponents();
     if(cend==-1) cend = x.nComponents();

     for(int i=rbegin; i<rend; ++i)
       for(int j=cbegin; j<cend; ++j)
       {
         if(origin->getRole(i) != "dual" && origin->getRole(j) != "dual" && !(hasTlin_uu && i==controlId && j==controlId) )
         {
                 if(i >= j) Tlin->d2axpy(a,y,x,i,i+1,j,j+1);
                 else Tlin->d2taxpy(a,y,x,j,j+1,i,i+1);
         }
         else
         {
           if(i >= j) Nlin->d2axpy(a,y,x,i,i+1,j,j+1);
           else Nlin->d2taxpy(a,y,x,j,j+1,i,i+1);
         }
       }
   }

  void LagrangeLinearization::d2taxpy(double a, AbstractFunctionSpaceElement& y, AbstractFunctionSpaceElement const& x, int rbegin, int rend, int cbegin, int cend) const
  {
    assert("not implemented");
  }

  /// Get point of linearization
 AbstractFunctionSpaceElement const& LagrangeLinearization::getOrigin() const { return *origin; }

  /// Precompute data
  // be careful, are Lagrange-Multipliers (needed for Tlin) up to date?
 void LagrangeLinearization::precompute()
  {
    Nlin->precompute();
    if(Tlin!=nullptr) Tlin->precompute();
  }

 void LagrangeLinearization::flush()
  {
    Nlin->flush();
    if(Tlin!=nullptr) Tlin->flush();
  }

//   void touch()
//    {
//      Nlin->touch();
//      Tlin->touch();
//    }
//
//   void connectToSignalForFlush(boost::signals2::signal0<void>& sig)
//    {
//      Nlin->connectToSignalForFlush(sig);
//      Tlin->connectToSignalForFlush(sig);
//    }

  AbstractLinearization const& LagrangeLinearization::getTangentialLinearization() const
  {
    //assert(Tlin);
    if(Tlin!=nullptr) return *Tlin;
    return *Nlin;
  }

  AbstractLinearization& LagrangeLinearization::getTangentialLinearization()
  {
    //assert(Tlin);
    if(Tlin!=nullptr) return *Tlin;
    return *Nlin;
  }

  AbstractLinearization& LagrangeLinearization::getNormalLinearization() { return *Nlin; }
  AbstractLinearization const& LagrangeLinearization::getNormalLinearization() const { return *Nlin; }
}
