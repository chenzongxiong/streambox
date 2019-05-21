#ifndef SDC_EULER_HH
#define SDC_EULER_HH

#include <algorithm>
#include <functional>
#include <vector>
#include "timestepping/sdc.hh"
#include "fem/fixdune.hh"

namespace Kaskade
{
  //FIXME:Adapt and move the norms to Kaskade::LinAlg::Norm
  /**
    * \brief A function to compute the partial \f$1-\f$norm or \f$L^{1}-\f$norm of a vector of vectors.
    * If \f$\mathbf{x}\f$ denotes a vector of vectors, where \f$\mathbf{x} = \{ \mathbf{x_1},\ldots , \mathbf{x_n}\}^{T}\f$,
    * then \f$\| \mathbf{x} \|_1 = \sum_{i=1}^N \|\mathbf{x_i}\|$\f. Here the norm \f$\|.\| \f$ is some arbitrary norm and
    * in the current template function it is implemented as a \f$1-\f$norm.
    *
    * \tparam Vector, a vector type, usually Dune::DenseVector
    *
    * \param[in]values a vector of vectors.
    *
    * \return the \f$1-\f$norm of the vector values, usually a double.
    */
  template<class Vector>
  typename Vector::field_type oneNorm(std::vector<Vector> const& values)
  {
    return std::accumulate(values.begin(),values.end(),0.0,[](double init, Vector const& v)
    {
        return init + v.one_norm();
      });
//    auto norm = 0.0;

//    for (auto i = 0u; i < values.size(); ++i)
//      //one_norm: public member function of Dune::DenseVector
//      norm += values[i].one_norm();

//    return norm;
  }

  /**
    *\brief A function to compute the partial \f$\max-\f$norm or \f$L^{\infty}-\f$norm of a vector of vectors.
    * If \f$\mathbf{x}\f$ denotes a vector of vectors, where \f$\mathbf{x} = \{ \mathbf{x_1},\ldots , \mathbf{x_n}\}^{T}\f$,
    * then \f$\| \mathbf{x} \|_{\infty} = \max_{1 \leq i \leq N} \|\mathbf{x_i}\|$\f. Here the norm \f$\|.\| \f$ is some arbitrary norm and
    * in the current template function it is implemented as a \f$\max-\f$norm.
    *
    * \tparam Vector, a vector type, usually Dune::DenseVector
    *
    * \param[in]values a vector of vectors.
    *
    * \return the \f$\infty-\f$norm of the vector values, usually a double.
    */
  template<class Vector>
  typename Vector::field_type maxNorm(std::vector<Vector> const& values)
  {
    typename Vector::field_type max = values[0].infinity_norm();

    for (auto i = 1u; i < values.size(); ++i)
      //infinity_norm: public member function of Dune::DenseVector
      max = std::max(max, values[i].infinity_norm());

    return max;
  }


  //FIXME: make the sdcExEulerIterationStep generic, instead of std::vector<Vector> replace with generic class Vectors.
  /**
    *\brief A single spectral defect correction iteration sweep based on explicit Euler method.
    *
    * This function performs one spectral defect correction (SDC) iteration for a system
    * of ODE's of the form:
    * \f[y'(t) = f(y(t)),\,\, y(0) = y_0 \textrm{ and } t \in [0, T].\f]
    * using the explicit Euler method  on a given time grid. Here SDC is interpreted as fixed point iteration.
    * Given an approximate solution \f$y^{[j]} \in \mathbb{P}_N\f$, the error function
    * \f[\delta^{[j]} = y - y^{[j]}\f]
    * satisfies the defect equation
    * \f[\delta^{[j]}'(t) = y'(t) - y^{[j]}'(t) = f(y(t)) - y^{[j]}'(t).\f]
    * The equivalent Picard equation is
    * \f[\delta^{[j]}(t) = \int_{\tau=0}^t \left( f(y^{[j]}(\tau)+\delta^{[j]}(\tau)) - {y^{[j]}}'(\tau) \right) \,d\tau.\f]
    * Evaluated at the grid nodes \f$ t_i\f$ we obtain
    * \f[
    * \delta^{[j]}_i &= \delta^{[j]}_{i-1} + \int_{\tau=t_{i-1}}^{t_i} \left( f(y^{[j]}(\tau)+\delta^{[j]}(\tau)) - {y^{[j]}}'(\tau) \right) \,d\tau \\
    *                &= \delta^{[j]}_{i-1} + \int_{\tau=t_{i-1}}^{t_i} \left( f(y^{[j]}(\tau)+\delta^{[j]}(\tau)) - f(y^{[j]}(\tau))\right) \,d\tau \\
    *                &\quad + \int_{\tau=t_{i-1}}^{t_i} \hspace{-0.5em}f(y^{[j]}(\tau)) \,d\tau - ( y^{[j]}_i-y^{[j]}_{i-1})
    * \f]
    * starting at \f$\delta_0^{[j]} = 0.\f$ Using left looking rectangular rule for approximating the first integral, and the canonical quadrature by
    * polynomial interpolation on the nodes \f$ t_1,\ldots t_N\f$ for second integral, approximate values \f$ \hat{\delta}_i^{[j]}\f$ for
    * \f$ \delta_i^{[j]}\f$ can be evaluated. With the left looking rectangular rule, we obtain the explicit scheme
    * \f[
    *  \hat\delta^{[j]}_i &= \hat\delta^{[j]}_{i-1}
    * + (t_i-t_{i-1}) \left( f(y^{[j]}_{i-1}+\hat{\delta}^{[j]}_{i-1}) - f(y^{[j]}_{i-1})\right)
    * +  \sum_{k=1}^N S_{ik} f(y^{[j]}_k) - y^{[j]}_i+y^{[j]}_{i-1},
    * \f]
    * where the entries of the spectral quadrature matrix \f$ S \in \mathbb{R}^{N \times N}\f$ are defined in terms of the Lagrange polynomials
    * \f$L_k \in \mathbb{P}_N\f$ satisfying \f$L_k(t_i) = \delta_{ik}\f$ as
    * \f[
    *   S_{ik} = \int_{\tau=t_{i-1}}^{t_i} L_k(\tau) \,d\tau, \quad i,k =1,\dots, N.
    * \f]
    * An improved approximation \f$ y^{[j+1]}\f$ is then obtained by polynomial interpolation of \f$\hat \delta^{[j]}_i \f$,
    * \f[y^{[j+1]} = y^{[j]} + \hat\delta^{[j]}.\f]
    * The norm of the correction \f$\hat\delta^{[j]}\f$ is returned.
    *
    * \tparam Vector, a vector type, usually Dune::DenseVector.
    *
    * \param[in] grid the collocation time grid
    * \param[out] yi stores the current iterate of the solution \f$y\f$.
    * \param[out] dyi is the approximate correction \f$\delta y \f$
    * \param[in] rhsFunc is the function which represents the right hand side \f$f\f$ and is used in th computation of \f$ \delta y\f$.
    * \param[in] normFunc represents the norm that would be used in the computation of norm of \f$ \delta y\f$.
    * \param[in] verbose a boolean used to print out the iterates of \f$y\f$.
    *
    * \return the norm of dyi, the type of the norm is determined by input parameter normFunc.
    */

  //single SDC iteration step, which returns
  template <class Vector>
  typename Vector::field_type sdcExEulerIterationStep(Kaskade::SDCTimeGrid const& grid, std::vector<Vector> & yi, std::vector<Vector> & dyi,
                                                      std::function<Vector(typename Vector::field_type, Vector const&)> rhsFunc,
                                                      std::function<typename Vector::field_type(std::vector<Vector> const&)> normFunc,
                                                      bool verbose = false)
  {
    //extracting time points from the time grid
    auto const& pts = grid.points();
    //extracting values of the integration matrix corresponding to the time grid
    auto const& integ = grid.integrationMatrix();
    //number of subintervals
    int const n = pts.size() - 1;

    //including starting point for yi and dyi
    assert(yi.size() == n+1);
    assert(dyi.size() == n+1);

    //initialize the correction at starting point to zero, where dyi[0] is a vector in general
    //and every component of dyi is set to 0.0 by assignment below.
    dyi[0] = 0.0;

    //declare Vector total to be used inside the for loop
    Vector total(0.0);
    typename Vector::field_type norm = 0.0;

    //The Sdc step for explicit Euler
    //perform n Euler steps
    boost::timer::cpu_timer timer;
    for (int i = 1; i <= n; ++i)
      {
        //Initialize all entries of Vector total to 0.
        for (size_t t = 0; t < total.size(); ++t)
          total[t] = 0.0;
        //Computation of the Lagrange interpolation
        for (int k = 0; k < n; ++k)
          total += integ[i-1][k] * rhsFunc(pts[k], yi[k]);
        //approximate correction computation.
        dyi[i] = dyi[i-1] + (pts[i] - pts[i-1]) *(rhsFunc(pts[i-1],yi[i-1]+dyi[i-1]) - rhsFunc(pts[i-1],yi[i-1])) + total - yi[i] + yi[i-1];

      }

    //compute the next iteration of yi.
    for (int i = 0; i <= n; ++i)
      {
        yi[i] += dyi[i];
        if(verbose) std::cout << "yi[" << i << "] = " << yi[i] << "\n";
      }

    norm = normFunc(dyi);
    return norm;
  }

  /**
   * \brief Base class to perform SDC iterations based on forward Euler method. Total iterations performed depends on the value set for variable maxSDCiterations.
   *
   * \tparam Vector, a vector type, usually Dune::DenseVector.
   * \tparam TimeGrid, a time grid type as implemented in SDCTimeGrid, LobattoTimeGrid is the default choice.
   */


  template <class Vector, class TimeGrid=LobattoTimeGrid>
  class EulerSDC
  {
  public:
    EulerSDC(double t0, double t1, size_t nIntervals, bool verbose_=false)
      : tbegin(t0), tend(t1), verbose(verbose_), grid(nIntervals,tbegin,tend)
    {}

    EulerSDC(double t0, double t1, size_t nIntervals, size_t maxSDCIter, bool verbose_=false)
      : tbegin(t0), tend(t1), maxSDCiterations(maxSDCIter), verbose(verbose_), grid(nIntervals,tbegin,tend)
    {}

    EulerSDC(double t0, double t1, size_t nIntervals, size_t maxSDCIter, double tol, bool verbose_=false)
      : tbegin(t0), tend(t1), maxSDCiterations(maxSDCIter), tolerance(tol), verbose(verbose_), grid(nIntervals, tbegin, tend)
    {}


    void setMaxSDCIterations(size_t maxIter)
    {
      maxSDCiterations = maxIter;
    }

    void setTolerance(double tau)
    {
      tolerance = tau;
    }


    //integrate function to carry out SDC iterations based on explicit Euler method based on the number of iterations(maxSDCiterations) given by user.
    Vector integrate( Vector const& initialValue, std::function<Vector (double,Vector const&)> rightHandSide,
                      std::function<typename Vector::field_type(std::vector<Vector const&>)> normFunc)
    {
      //initialize the vector yi
      std::vector<Vector> yi(grid.points().size(),initialValue), dyi(grid.points().size(),Vector(0.0));

      //print the entries of the vector yi
      if(verbose)
        {
          std::cout << "print entries to the vector yi" << std::endl;
          for (int i = 0; i < yi.size(); ++i)
            std::cout << "[ " << yi[i] << " ]" << std::endl;
        }

      for (int i = 0; i < maxSDCiterations; ++i)
        {
          if(verbose) std::cout << "Iteration : " << i << "\n";
          auto stepNorm = sdcExEulerIterationStep(grid, yi, dyi, rightHandSide, normFunc, verbose);
          if(verbose) std::cout << "norm = " << stepNorm << std::endl;
        }

      return yi.back();
    }

    //integrate function to carry out SDC iterations based on explicit Euler method based on the tolerance given by the user.
    //In case the tolerance is not reached within the default value of maxSDCiterations, the last value of yi at maxSDCiterations is returned.
    Vector integrateTOL( Vector const& initialValue, std::function<Vector(double,Vector const&)> rightHandSide,
                         std::function<typename Vector::field_type(std::vector<Vector> const&)> normFunc)
    {
      //initialize the vector yi with initial value and the correction vector dyi to 0 at all grid points.
      std::vector<Vector> yi(grid.points().size(),initialValue), dyi(grid.points().size(),Vector(0.0));
      //initialize loop variable to keep track of when loop exceeds maxSDCiterations
      int loop = 1;

      //print the entries of the vector yi
      if(verbose)
        {
          std::cout << "print entries to the vector yi" << std::endl;
          for (int i = 0; i < yi.size(); ++i)
            std::cout << "[ " << yi[i] << " ]" << std::endl;
        }

      //while loop to carry out sdc iterations till the given tolerance is reached failing which it loops till maxSDCiterations is reached.
      auto stepNorm = sdcExEulerIterationStep(grid, yi, dyi, rightHandSide, normFunc, verbose);
      while (stepNorm > tolerance  && loop < maxSDCiterations)
        {
          //std::cout << "step norm = " << stepNorm << "\n";
          if(verbose) std::cout << "Iteration: " << loop << "\n";
          stepNorm = sdcExEulerIterationStep(grid, yi, dyi, rightHandSide, normFunc, verbose);
          if(verbose) std::cout << "norm = " << stepNorm << std::endl;
          ++loop;
        }
      return yi.back();
    }

  private:
    double tbegin = 0, tend = 1;
    size_t maxSDCiterations = 100;
    double tolerance = 0.1;
    bool verbose = false;
    TimeGrid grid;
  };

}

#endif // SDC_EULER_HH
