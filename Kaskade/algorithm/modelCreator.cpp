#include "modelCreator.hh"

#include "abstract_interface.hh"
#include "lagrangeLinearization.hh"
#include "opt_model_functions.hh"

namespace Kaskade
{
  /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
  /* * * * * * * * * * * QuadraticModelCreator * * * * * * * * * */
  /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
  QuadraticModelCreator::QuadraticModelCreator(AbstractFunctionSpaceElement const& normalStep, std::vector<std::shared_ptr<AbstractFunctionSpaceElement> > const& tangentialBasis, LagrangeLinearization const& lagrangeLinearization, double residualCorrection)
  : p0c0(0), Lxxdndn(0), Lxdt(tangentialBasis.size()),Lxxdndt(tangentialBasis.size()), Lxxdtdt(tangentialBasis.size(),tangentialBasis.size())
  {
    AbstractFunctionSpaceElement const& origin=lagrangeLinearization.getOrigin();
    std::unique_ptr<AbstractFunctionSpaceElement> tmp(origin.clone());
    lagrangeLinearization.evald(*tmp);

    // p0 c'(x_0)
    p0c0 = origin.applyAsDualTo_role(*tmp,"dual");
    p0c0 -= residualCorrection;

    // L_x(x_0,p_0) dt
    for(int k=0; k<tangentialBasis.size();++k)
      Lxdt[k] = tmp->applyAsDualTo_role(*tangentialBasis[k],"primal");

    // tmp = L_xx(x_0,p_0)dn
    *tmp *= 0.0;

    for(int i=0; i<normalStep.nComponents(); ++i)
      for(int k=0; k<normalStep.nComponents(); ++k)
      {
        if(normalStep.getRole(i)=="primal" && tmp->getRole(k)=="primal")
          lagrangeLinearization.ddxpy(*tmp, normalStep, i, i+1, k, k+1);
      }

    // L_xx(x_0,p_0)(dn,dn)
    Lxxdndn = tmp->applyAsDualTo_role(normalStep,"primal");

    // L_xx(x_0,p_0)(dn,dt)
    for(int k=0; k<tangentialBasis.size();++k)
      Lxxdndt[k] = tmp->applyAsDualTo_role(*tangentialBasis[k],"primal");

    // L_xx(x_0,p_0)(dt,dt)
    for(int i=0; i<tangentialBasis.size();++i)
    {
      *tmp *= 0.0;
      // tmp = L_xx(x_0,p_0)dt[i]
      for(int i0=0; i0<normalStep.nComponents(); ++i0)
        for(int k0=0; k0<normalStep.nComponents(); ++k0)
        {
          if(tangentialBasis[i]->getRole(i0)=="primal" && tmp->getRole(k0)=="primal")
            lagrangeLinearization.ddxpy(*tmp, *tangentialBasis[i], i0, i0+1, k0, k0+1);
        }

      // Lxxdtdt[i][k]=L_xx(x_0,p_0)(dt[i],dt[k])
      for(int k=0; k<tangentialBasis.size();++k)
      {
        Lxxdtdt[i][k]=tmp->applyAsDualTo_role(*tangentialBasis[k],"primal");
      }
    }
  }

  QuadraticFunction QuadraticModelCreator::create(double nu) const
  {
    Dune::BlockVector<Dune::FieldVector<double,1>> dndt(Lxxdndt);   dndt *= nu;
    Dune::Matrix<Dune::FieldVector<double,1>> dtdt(Lxxdtdt);        dtdt *= 0.5;
    dndt+=Lxdt;
    return QuadraticFunction(-p0c0*nu+Lxxdndn*0.5*nu*nu,dndt,dtdt);
  }

  /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
  /* * * * * * * * * * * * NormModelCreator  * * * * * * * * * * */
  /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
  NormModelCreator::NormModelCreator(AbstractFunctionSpaceElement const& normalStep, std::vector<std::shared_ptr<AbstractFunctionSpaceElement> > const& tangentialBasis, AbstractScalarProduct const& sp)
  : dndn(0), dndt(tangentialBasis.size()), dtdt(tangentialBasis.size(),tangentialBasis.size())
  {
    // (dn,dn)
    dndn =sp(normalStep,normalStep);

    // (dn,dt)
    for(int k=0; k<tangentialBasis.size();++k)
      dndt[k] =sp(normalStep,*tangentialBasis[k]);

    // (dt,dt)
    for(int i=0; i<tangentialBasis.size();++i)
    {
      for(int k=0; k<tangentialBasis.size();++k)
      {
        dtdt[i][k]=sp(*tangentialBasis[i],*tangentialBasis[k]);
      }
    }
  }


  QuadraticFunction NormModelCreator::create(double nu) const
  {
    Dune::BlockVector<Dune::FieldVector<double,1>> dndt_(dndt);
    dndt_*=2*nu;
    return QuadraticFunction(dndn*nu*nu,dndt_,dtdt);
  }

}
