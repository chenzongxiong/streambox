#ifndef CGTERMINATIONCRITERIA_HH
#define CGTERMINATIONCRITERIA_HH

namespace Kaskade
{
  /**
   * \ingroup linalgsolution
   * \brief Interface for IterateType::PCG termination criterion policy classes
   * \tparam R a floating point type for real numbers
   */
  template <class R>
  class PCGTerminationCriterion {
  public:
    /**
     * \brief real field type
     */
    typedef R Real;

    /**
     * \brief re-initializes the termination criterion for a new IterateType::CG run
     */
    virtual void clear() = 0;

    /**
     * \brief set requested tolerance
     *
     * \param tol the requested tolerance (nonnegative)
     */
    virtual void setTolerance(Real tol) = 0;

    /**
     * \brief set requested look-ahead count
     *
     * \param lah the requested look-ahead count (positive integer)
     */
    virtual void setLookAhead(int lah) = 0;
    /**
     * @brief addStepQuantities supplies algorithmic quantities to the termination criterion
     * @param stepLength scaling for the conjugate search direction \f$q\f$
     * @param qAq squared energy norm of the conjugate search direction \f$q\f$
     * @param qPq squared \f$P\f$-norm, i. e. the norm induced by the preconditioner, of the conjugate search direction \f$q\f$
     * @param rPINVr squared \f$P^{-1}\f$-norm of the residual
     */
    virtual void addStepQuantities(Real stepLength, Real qAq, Real qPq, Real rPINVr) = 0;


    /**
     * \brief get the maximum number  of allowed iteration steps
     */
    virtual int getMaxIterationSteps() = 0;
    /**
     * \brief termination decision
     * \return true if the iteration has reached the required accuracy
     */
     
    virtual operator bool() = 0;

    virtual bool minimalDecreaseAchieved() { assert("not implemented"); return false; }
  };


  /**
   * Relative energy error termination criterion according to Strakos, Tichy: Error estimation in preconditioned conjugate gradients.
   * Requires that IterateType::CG starts at \f$ x = 0 \f$. More general starting values might be used, but must be chosen such that
   * the estimate for the energy norm of the solution stays positive (see the above mentioned paper for details).
   */
  template <class R>
  class StrakosTichyEnergyErrorTerminationCriterion: public PCGTerminationCriterion<R> {
  public:
    typedef R Real;

    /**
     * \brief constructor
     *
     * The pcg iteration is terminated as soon as either the estimated error satisfies \f$ [\epsilon] \le \mathrm{tol} \f$ or
     * the number of iterations exceeds the limit maxit. Note that the estimate of the relative error requires a look ahead
     * parameter L. Thus, if \f$\mathrm{maxit}\ge\mathrm{L}\f$, then the error is first estimated in the (L+1)-th iteration.
     *
     * \param tol the relative error tolerance for termination
     * \param maxit the maximum number of iterations
     * \param eps maximal attainable accuracy
     */
    StrakosTichyEnergyErrorTerminationCriterion(Real tol_, int maxit_, double eps_ = 1e-12)
      : tol(tol_), maxit(maxit_), eps(eps_)
    {}

    /**
     * \brief constructor
     *
     * The pcg iteration is terminated as soon as either the estimated error satisfies \f$ [\epsilon] \le \mathrm{tol} \f$ or
     * the number of iterations exceeds the limit maxit. Note that the estimate of the relative error requires a look ahead
     * parameter L. Thus, if \f$\mathrm{maxit}\ge\mathrm{L}\f$, then the error is first estimated in the (L+1)-th iteration.
     *
     * \param tol the relative error tolerance for termination
     * \param minTol relative error tolerance to admit truncation in the hybrid cg implementation
     * \param maxit the maximum number of iterations
     * \param eps maximal attainable accuracy
     */
    StrakosTichyEnergyErrorTerminationCriterion(Real tol_, Real minTol_, int maxit_, double eps_ = 1e-12)
      : tol(tol_), maxit(maxit_), eps(eps_), minTol(minTol_)
    {}

    /**
     * \brief re-initializes the termination criterion for a new IterateType::CG run
     */
    virtual void clear()
    {
      scaledGamma2.clear();
      energyNorm = 0;
    }

    /**
     * \brief set requested relative tolerance
     *
     * \param tol the requested tolerance (nonnegative)
     */
    virtual void setTolerance(Real tol_) { tol = tol_; }//minTol = sqrt(tol); }

    /**
     * \brief set requested lookahead value
     *
     * \param lah the requested lookahead (nonnegative)
     *
     * The default value is 50.
     */
    virtual void setLookAhead(int d_) { d = d_; }

    /**
     * @brief addStepQuantities supplies algorithmic quantities to the termination criterion
     * @param stepLength scaling for the conjugate search direction \f$q\f$
     * @param qAq squared energy norm of the conjugate search direction \f$q\f$ (here: unused)
     * @param qPq squared \f$P\f$-norm, i. e. the norm induced by the preconditioner, of the conjugate search direction \f$q\f$ (here: unused)
     * @param rPINVr squared \f$P^{-1}\f$-norm of the residual
     */
    virtual void addStepQuantities(Real stepLength, Real, Real, Real rPINVr)
    {
      scaledGamma2.push_back( stepLength * rPINVr );
      energyNorm += stepLength * rPINVr;
    }

    /**
     * \brief get the maximum number  of allowed iteration steps
     */
    virtual int getMaxIterationSteps()
    {
      return maxit;
    }
    /**
     * \brief termination decision
     * \return true if the iteration has reached the required accuracy
     */
    virtual operator bool()
    {
      return scaledGamma2.size() > d && relativeError() < std::max(eps,tol*tol);
    }

    /**
     * \brief returns the estimated absolute energy error
     */
    Real relativeError()
    {
      if( scaledGamma2.size() < d ) return std::numeric_limits<Real>::max();
      return std::accumulate(scaledGamma2.end() - d, scaledGamma2.end(), 0.) / energyNorm;
    }

    bool minimalDecreaseAchieved() { return relativeError() < minTol; }

  private:
    Real tol;
    int maxit;
    // squared gammas
    std::vector<Real> scaledGamma2;
    Real energyNorm = 0;
    int d = 50;
    double eps = 1e-12, minTol = 1e-2;
  };


  /**
   * Relative error termination criterion based on the norm induced by the preconditioner,
   * according to Strakos, Tichy: Error estimation in preconditioned conjugate gradients.
   * Requires that IterateType::CG starts at \f$ x = 0 \f$. More general starting values might be used, but must be chosen such that
   * the estimate for the energy norm of the solution stays positive (see the above mentioned paper for details).
   */
  template <class R>
  class StrakosTichyPTerminationCriterion: public PCGTerminationCriterion<R> {
  public:
    typedef R Real;

    /**
     * \brief constructor
     *
     * The pcg iteration is terminated as soon as either the estimated error satisfies \f$ [\epsilon] \le \mathrm{tol} \f$ or
     * the number of iterations exceeds the limit maxit. Note that the estimate of the relative error requires a look ahead
     * parameter L. Thus, if \f$\mathrm{maxit}\ge\mathrm{L}\f$, then the error is first estimated in the (L+1)-th iteration.
     *
     * \param tol the relative error tolerance for termination
     * \param maxit the maximum number of iterations
     * \param eps maximal attainable accuracy
     */
    StrakosTichyPTerminationCriterion(Real tol_, int maxit_, double eps_ = 1e-12)
      : tol(tol_), maxit(maxit_), eps(eps_)
    {}

    /**
     * \brief constructor
     *
     * The pcg iteration is terminated as soon as either the estimated error satisfies \f$ [\epsilon] \le \mathrm{tol} \f$ or
     * the number of iterations exceeds the limit maxit. Note that the estimate of the relative error requires a look ahead
     * parameter L. Thus, if \f$\mathrm{maxit}\ge\mathrm{L}\f$, then the error is first estimated in the (L+1)-th iteration.
     *
     * \param tol the relative error tolerance for termination
     * \param minTol relative error tolerance to admit truncation in the hybrid cg implementation
     * \param maxit the maximum number of iterations
     * \param eps maximal attainable accuracy
     */
    StrakosTichyPTerminationCriterion(Real tol_, Real minTol_, int maxit_, double eps_ = 1e-12)
      : tol(tol_), maxit(maxit_), eps(eps_), minTol(minTol_)
    {}

    /**
     * \brief re-initializes the termination criterion for a new IterateType::CG run
     */
    virtual void clear()
    {
      scaledGamma2.clear();
      energyNorm = 0;
    }

    /**
     * \brief set requested relative tolerance
     *
     * \param tol the requested tolerance (nonnegative)
     */
    virtual void setTolerance(Real tol_) { tol = tol_; }

    /**
     * \brief set requested lookahead value
     *
     * \param lah the requested lookahead (nonnegative)
     *
     * The default value is 50.
     */
    virtual void setLookAhead(int d_) { d = d_; }

    /**
     * @brief addStepQuantities supplies algorithmic quantities to the termination criterion
     * @param stepLength scaling for the conjugate search direction \f$q\f$
     * @param qAq squared energy norm of the conjugate search direction \f$q\f$ (here: unused)
     * @param qPq squared \f$P\f$-norm, i. e. the norm induced by the preconditioner, of the conjugate search direction \f$q\f$ (here: unused)
     * @param rPINVr squared \f$P^{-1}\f$-norm of the residual
     */
    virtual void addStepQuantities(Real stepLength, Real qAq, Real qPq, Real rPINVr)
    {
      scaledGamma2.push_back( stepLength * rPINVr );
      energyNorm += stepLength * rPINVr;
      steps2.push_back(qPq/qAq);
    }

    /**
     * \brief get the maximum number  of allowed iteration steps
     */
    virtual int getMaxIterationSteps()
    {
      return maxit;
    }
    /**
     * \brief termination decision
     * \return true if the iteration has reached the required accuracy
     */
    virtual operator bool()
    {
      return scaledGamma2.size() > 2*d && relativeError() < std::max(eps,tol*tol);
    }

    /**
     * \brief returns the estimated absolute energy error
     */
    Real relativeError()
    {
      if( scaledGamma2.size() < 2*d ) return std::numeric_limits<Real>::max();

      size_t j = scaledGamma2.size() - 2*d;
      double tau = 0;
      for(size_t i = j; i < j + d; ++i)
        tau += steps2[i] * ( scaledGamma2[i] + 2 * std::accumulate(scaledGamma2.begin()+i+1, scaledGamma2.end(),0)  );
      return tau;
    }

    bool minimalDecreaseAchieved() { return relativeError() < minTol; }

  private:
    Real tol;
    int maxit;
    // squared gammas
    std::vector<Real> scaledGamma2, steps2;
    Real energyNorm = 0;
    int d = 50;
    double eps = 1e-12, minTol = 1e-2;
  };

}

#endif // CGTERMINATIONCRITERIA_HH
