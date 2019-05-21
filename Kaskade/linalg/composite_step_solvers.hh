#ifndef COMPOSITE_STEP_SOLVERS
#define COMPOSITE_STEP_SOLVERS

#include <memory> // std::unique_ptr
#include <vector>

#include <boost/timer/timer.hpp>

#include "algorithm/opt_interface.hh"
#include "algorithm/newton_bridge.hh"
#include "linalg/triplet.hh"
#include "linalg/tcg.hh"

namespace Kaskade
{
  namespace Bridge
  {
    template<class DirectSolver>
    class DirectInnerSolver
    {
    public:
      DirectInnerSolver(int numberOfBlocks_,bool stabilization_=true) : numberOfBlocks(numberOfBlocks_), stabilization(stabilization_) {}

      // on exit: sol=(dx, p), where dx is normal step, p is least squares update of Lagrange multiplier
      void solveAdjAndNormal(std::vector<double>& sol, SparseLinearSystem const& lin)
      {
        MatrixAsTriplet<double> mat;
        lin.getMatrix(mat);
        //      mat.print();
        if(stabilization)
        {
          // stabilization: cf. Conn/Gould/Toint pg. 110, (5.4.7)
          AT.flush();
          lin.getMatrixBlocks(AT,0,numberOfBlocks,numberOfBlocks,lin.nColBlocks());
          // shift c'(x)^* to block matrix (0 c'(x)^*)
          AT.shiftIndices(0,lin.cols(0,numberOfBlocks));
          // AT is (0 c'(x)^*)
          ATx.resize(0);
        }
        solver.reset(new DirectSolver(mat.nrows(),
            2,
            mat.ridx,
            mat.cidx,
            mat.data,
            MatrixProperties::GENERAL));
        std::vector<double> r,s, soladj;
        sz=lin.size();
        lin.getRHSBlocks(r,0,numberOfBlocks);
        lin.getRHSBlocks(s,numberOfBlocks,lin.nRowBlocks());

        std::vector<double> rg(sz,0.0),rc(sz,0.0);
        for(int i=0;i<r.size();++i) rg[i]=r[i];
        for(int i=0;i<s.size();++i) rc[i+r.size()]=s[i];

        solver->solve(rg,soladj);
        solver->solve(rc,sol);
        for(int i=r.size(); i<sol.size();++i) sol[i]=soladj[i];

      }

      // Application of "preconditioner"
      // stabilization: cf. Conn/Gould/Toint pg. 110, (5.4.7)
      void ax(std::vector<double>& sol, std::vector<double>const &r) const
      {
        std::vector<double> r1(sz,0.0),x;
        if(stabilization && ATx.size()!=0)
          for(int i=0;i<r.size();++i)
            r1[i]=r[i]-ATx[i];
        else
          for(int i=0;i<r.size();++i)
            r1[i]=r[i];
        solver->solve(r1,x);
        if(stabilization) AT.ax(ATx,x);
        for(int i=0; i<sol.size(); ++i) sol[i]=x[i];
      }

      void resolveNormal(std::vector<double>& sol, SparseLinearSystem const& lin)
      {
        std::vector<double> s,r2(sz,0.0);
        lin.getRHSBlocks(s,numberOfBlocks,lin.nRowBlocks());
        for(int i=0;i<s.size();++i) r2[i+lin.rows(0,numberOfBlocks)]=s[i];
        solver->solve(r2,sol);
        for(int i=lin.rows(0,numberOfBlocks); i<sol.size();++i) sol[i]=0.0;

      }
    private:
      MatrixAsTriplet<double> AT;
      mutable std::vector<double> ATx;
      std::unique_ptr<DirectSolver> solver;
      int numberOfBlocks, sz;
      bool stabilization;
    };

    template <class VectorImpl,class InnerSolver>
    class ProjTCGSolver : public AbstractTangentialSpace
    {
    public:
      virtual void setRelativeAccuracy(double accuracy) { acc=accuracy; };
      virtual   double getRelativeAccuracy() { return acc; }
      virtual   double getAbsoluteAccuracy() { return acc; }
      virtual   bool improvementPossible() { return true; }
      virtual   int nSolutionVectors() const { return 2; }
      ProjTCGSolver(InnerSolver& solver_, int numberOfBlocks_) : solver(solver_), numberOfBlocks(numberOfBlocks_) {}
      virtual ~ProjTCGSolver() {}
    private:

      double acc;

      virtual int doSolve(std::vector<AbstractVector* >& correction,
          AbstractLinearization& linT, AbstractLinearization& linN, int start, double nu0)
      {
        SparseLinearSystem const &lT = dynamic_cast<SparseLinearSystem const &>(linT);
        SparseLinearSystem const &lN = dynamic_cast<SparseLinearSystem const &>(linN);

        boost::timer::cpu_timer timer;

        std::vector<double> r(lN.rows(0,numberOfBlocks),0.0);

        // r=f'(x)
        lN.getRHSBlocks(r,0,numberOfBlocks);

	double rsum(0.0);
	for(int i=0; i<r.size(); ++i)
	  rsum += fabs(r[i]);


	std::cout << "r:" << rsum << " " << r.size() << " " << numberOfBlocks << std::endl;

        MatrixAsTriplet<double> m;
        int dimx=lT.rows(0,numberOfBlocks);
        std::vector<double> x,p,normalstep;
        x.reserve(lT.size());
        p.reserve(lT.size());

        if(start >=1)
          dynamic_cast<Bridge::Vector<VectorImpl>& >(*(correction[start-1])).write(normalstep);

        dynamic_cast<Bridge::Vector<VectorImpl>const & >(linN.getOrigin()).write(x );

        MatrixAsTriplet<double> H;
        // H = c'(x)^*
        lN.getMatrixBlocks(H,0,numberOfBlocks,numberOfBlocks,lN.nColBlocks());
        H.shiftIndices(0,lN.cols(0,numberOfBlocks));
        H.axpy(r,x);
        H.flush();
        // r = f'(x)+c'(x)^*p
        lT.getMatrixBlocks(H,0,numberOfBlocks,0,numberOfBlocks);
        // r = f'(x)+c'(x)^*p+nu_0 L_xx delta n
        H.axpy(r,normalstep,nu0);

        for(int i=0; i< r.size(); ++i)
        {
          r[i] *= -1.0;
          x[i] = 0.0;
        }

        // min <Hx,x>+<r,x>
        int exit=projectedtcg(x,p,H,solver,r,1e-6,100);

        for(int i=0; i<dimx;++i) x[i] *= -1.0;
        for(int i=dimx; i<lT.size();++i) x.push_back(0.0);
        dynamic_cast<Bridge::Vector<VectorImpl>& >(*(correction[start])).read(x);

        std::cout << " ProjTCG:" << (double)(timer.elapsed().user)/1e9 << std::endl;

        if(exit==2)
        {
          for(int i=0; i<dimx;++i) p[i] *= -1.0;
          for(int i=dimx; i<lT.size();++i) p.push_back(0.0);
          dynamic_cast<Bridge::Vector<VectorImpl>& >(*(correction[start+1])).read(p);
          return 2;
        }
        return 1;

      }
      InnerSolver& solver;
      int numberOfBlocks;
    };


    template <class VectorImpl,class InnerSolver>
    class TCGSolver : public AbstractTangentialSpace
    {
    public:
      virtual   void setRelativeAccuracy(double accuracy) { acc=accuracy; };
      virtual   double getRelativeAccuracy() { return acc; }
      virtual   double getAbsoluteAccuracy() { return acc; }
      virtual   bool improvementPossible() { return true; }
      virtual   int nSolutionVectors() const { return 2; }
      TCGSolver(InnerSolver& solver_, int numberOfBlocks_) : solver(solver_), numberOfBlocks(numberOfBlocks_) {}
      virtual ~TCGSolver() {}
    private:
      double acc;

      virtual int doSolve(std::vector<AbstractVector* >& correction,
          AbstractLinearization& linT, AbstractLinearization& linN, int start, double d1, double d2, double d3, double d4, double nu0)
      {
        SparseLinearSystem const &lT = dynamic_cast<SparseLinearSystem const &>(linT);
        SparseLinearSystem const &lN = dynamic_cast<SparseLinearSystem const &>(linN);
        MatrixAsTriplet<double> m;
        int dimx=lT.rows(0,numberOfBlocks);
        std::vector<double> x,p,normalstep;
        if(start >=1)
          dynamic_cast<Bridge::Vector<VectorImpl>& >(*(correction[start-1])).write(normalstep);

        solver.solveTCG(x,p,lT,lN,normalstep,nu0);

        for(int i=0; i<dimx;++i) x[i] *= -1.0;
        for(int i=dimx; i<x.size();++i) x[i] *= 0.0;
        dynamic_cast<Bridge::Vector<VectorImpl>& >(*(correction[start])).read(x);
        if(p.size() != 0)
        {
          for(int i=0; i<dimx;++i) p[i] *= -1.0;
          for(int i=dimx; i<p.size();++i) p[i] *= 0.0;
          dynamic_cast<Bridge::Vector<VectorImpl>& >(*(correction[start+1])).read(p);
          return 2;
        }
        return 1;

      }
      InnerSolver& solver;
      int numberOfBlocks;
    };

    template <class VectorImpl,class InnerSolver>
    class PINVSolver : public AbstractNormalDirection
    {
    public:

      virtual void setRelativeAccuracy(double accuracy) { acc=accuracy; };
      virtual   double getRelativeAccuracy() { return acc; }
      virtual   double getAbsoluteAccuracy() { return acc; }
      virtual   bool improvementPossible() { return true; }
      virtual ~PINVSolver() {}

      PINVSolver(InnerSolver& solver_, int numberOfBlocks_) : solver(solver_), numberOfBlocks(numberOfBlocks_) {};

      void computeCorrectionAndAdjointCorrection(Kaskade::AbstractVector&, Kaskade::AbstractVector&, Kaskade::AbstractLinearization&){assert("not implemented");}
      void computeSimplifiedCorrection(Kaskade::AbstractVector&, const Kaskade::AbstractLinearization&) const {assert("not implemented");}

    private:

      double acc;
      virtual void doSolve(AbstractVector& correction,
          AbstractVector& iterate,
          AbstractLinearization& lin)
      {
        SparseLinearSystem const &l = dynamic_cast<SparseLinearSystem const &>(lin);
        int dimx=l.rows(0,numberOfBlocks);
        if(dimx==l.size()) return;
        std::vector<double> xcor(l.size(),0.0);
        solver.solveAdjAndNormal(xcor,l);
        for(int i=0; i< xcor.size(); ++i) xcor[i]*=-1.0;

        std::vector<double> xiterate;
        dynamic_cast<Bridge::Vector<VectorImpl>& >(iterate).write(xiterate);


        // update dual variables of iterate
        for(int i=dimx; i< xcor.size(); ++i)
        {
          xiterate[i]= xcor[i];
          xcor[i]=0.0;
        }

        // read data into correction and iterate
        dynamic_cast<Bridge::Vector<VectorImpl>& >(correction).read(xcor);
        dynamic_cast<Bridge::Vector<VectorImpl>& >(iterate).read(xiterate);
      }


      virtual void doResolve(AbstractVector& correction,
          AbstractLinearization const& lin) const
      {
        SparseLinearSystem const &l = dynamic_cast<SparseLinearSystem const &>(lin);
        int dimx=l.rows(0,numberOfBlocks);
        if(l.size()==dimx) return;
        std::vector<double> x(l.size(),0.0);
        solver.resolveNormal(x,l);
        for(int i=0; i< dimx; ++i) x[i]*=-1.0;
        for(int i=dimx; i< x.size(); ++i) x[i]=0.0;

        dynamic_cast<Bridge::Vector<VectorImpl>& >(correction).read(x);
      }

      MatrixAsTriplet<double> m;
      mutable MatrixAsTriplet<double> P;
      InnerSolver& solver;
      int numberOfBlocks;
    };

    //   template<class LinImpl, class VectorImpl,class Preconditioner>
    //   class TCGWithPreconditioner : public AbstractTangentialSolver
    //   {
    //   public:
    //     virtual void setRelativeAccuracy(double accuracy) { acc=accuracy; };
    //     virtual   double getRelativeAccuracy() { return acc; }
    //     virtual   double getAbsoluteAccuracy() { return acc; }
    //     virtual   bool improvementPossible() { return true; }
    //     virtual   bool localConvergenceLikely() { return lCl; }

    //     virtual   int nSolutionVectors() const { return 2; }

    //     TCGWithPreconditioner(int primalblocks_) :
    //       Me(3), Ae(4), addreg(0.0), primalblocks(primalblocks_) , tcgfile("tcg.log",std::ios::out)
    //     {}

    //     virtual ~TCGWithPreconditioner() {};

    //     virtual bool getNormInfo(std::vector<double>& M,std::vector<double>& A) const
    //     {
    //       M=Me;
    //       A=Ae;
    //       return true;
    //     }

    //   private:

    //     std::vector<double> Me, Ae;
    //     double acc;
    //     double addreg;
    //     bool lCl;

    //     virtual int doSolve(std::vector<AbstractVector* >& correction,
    //                         AbstractLinearization& linT, AbstractLinearization& linN, int start,
    //                         double ThetaAim, double omegaC, double omegaL, double omegaH, double nu0)
    //     {
    //       SparseLinearSystem &lT = dynamic_cast<SparseLinearSystem &>(linT);
    //       SparseLinearSystem &lN = dynamic_cast<SparseLinearSystem &>(linN);

    //       MatrixAsTriplet<double> NM;

    //       lN.getMatrix(NM);


    //       boost::timer timer;

    //       std::vector<double> r(lN.rows(0,primalblocks),0.0);

    //       // r = f'(x)

    //       lN.getRHSBlocks(r,0,primalblocks);
    //       int dimx=lT.rows(0,primalblocks);

    //       if(primalblocks < lT.nColBlocks())
    //       {
    //         // if equality constraints are present, then r += C'(x)^T p
    //         // this is for numerical stability
    //         std::vector<double> x;
    //         dynamic_cast<Bridge::Vector<VectorImpl>const & >(linN.getOrigin()).write(x);
    //         MatrixAsTriplet<double> CPrimeTransposed;
    //         lT.getMatrixBlocks(CPrimeTransposed,0,primalblocks,primalblocks,lT.nColBlocks());
    //         CPrimeTransposed.shiftIndices(0,lN.cols(0,primalblocks));
    //         CPrimeTransposed.axpy(r,x);
    //       }

    //       MatrixAsTriplet<double> Hessian;
    //       lT.getMatrixBlocks(Hessian,0,primalblocks,0,primalblocks);

    //       if(start >=1)
    //       {
    //         // if normal step was taken, then r += nu_0 H delta n
    //         // this guarantees second order approximation

    //         std::vector<double> normalstep;
    //         dynamic_cast<Bridge::Vector<VectorImpl>& >(*(correction[start-1])).write(normalstep);
    //         Hessian.axpy(r,normalstep,nu0);
    //       }

    //       std::vector<double> x(lT.size(),0.0),p(lT.size(),0.0);

    // // Preconditioner is only created at the beginning of the whole algorithm

    // //      if(!(prec_ptr.get()))
    //         prec_ptr.reset(new Preconditioner(lN));


    //       int maxiter=1000;

    //       int exit(2);

    //       double negm(1e300), posdefalpha(1.0);

    // // for unregularized tcg:

    //       addreg = 0.0;

    //       int ntrial(0);

    //       do
    //       {
    //         PrecWrapperForStdVector<Preconditioner> prec(*prec_ptr);
    //         maxiter=1000;
    //         std::cout << "Regularization " << addreg << std::endl;
    //         exit=tcgRegForCubic(x,p,Hessian,prec,r,acc,omegaL,omegaH, addreg,maxiter, ntrial,Me,Ae,tcgfile);
    //         if(exit==-1)
    //           prec_ptr.reset(new Preconditioner(lN));

    //         std::cout << "Negative Curvature: " << Ae[3] << std::endl;
    //         negm=std::min(negm,Ae[3]);
    //         if(addreg >0 ) posdefalpha=0.0;
    //         ntrial++;
    //         lCl = (addreg == 0.0);
    //       } while(exit !=1);

    //       addreg=std::max(0.0,-negm);

    // //        if(posdefalpha>0.0)
    // //              prec_ptr.reset(new Preconditioner(lT));
    //                       //Preconditioner preconditioner(lN);

    //       std::vector<double> Mx(x.size());

    //       NM.ax(Mx,x);
    //       double sum(0.0);
    //       for(int i=0; i<x.size();++i)
    //         sum += x[i]*Mx[i];

    //       std::cout << "xMx:" << sum << " "  << Me[0] << " ";

    //       sum=0.0;

    //       for(int i=0; i<x.size();++i)
    //         sum += p[i]*Mx[i];

    //       std::cout << "pMx:" << sum << " " << Me[2] << " ";

    //       NM.ax(Mx,p);
    //       sum=0.0;
    //       for(int i=0; i<x.size();++i)
    //         sum += p[i]*Mx[i];

    //       std::cout << "pMp:" << sum << " " << Me[1] << std::endl;


    //       // no update of lagrangian multiplier
    //       for(int i=dimx; i<lT.size();++i) x.push_back(0.0);
    //       dynamic_cast<Bridge::Vector<VectorImpl>& >(*(correction[start])).read(x);

    //       std::cout << " ProjTCG: exit:" << exit << " t:" << timer.elapsed() << " it:" << maxiter << std::endl;

    //       if(exit==2)
    //       {
    //         for(int i=dimx; i<lT.size();++i) p.push_back(0.0);
    //         dynamic_cast<Bridge::Vector<VectorImpl>& >(*(correction[start+1])).read(p);
    //         return 2;
    //       }
    //       return 1;

    //     }
    //     int primalblocks;
    //     std::ofstream tcgfile;
    //     std::unique_ptr<Preconditioner> prec_ptr;
    //   };


  }
} // namespace Kaskade
#endif
