/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2002-2012 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*
 * preconditioner_factory.hh
 *
 * @TODO
 * export Dune's linker symbols (at least for Dune::Preconditioner, Dune::InverseOperator)
 * in order to allow the usage of dynamic_cast
 *
 *  Created on: 09.02.2012
 *      Author: lars lubkoll
 */

#ifndef PRECONDITIONER_FACTORY_HH_
#define PRECONDITIONER_FACTORY_HH_

#include <iostream>
#include <stdexcept>

#include <boost/preprocessor/cat.hpp>

#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>

#ifndef KASKADE_NO_TRIVIAL_PRECONDITIONER
#include "linalg/trivialpreconditioner.hh"
#endif

#ifndef KASKADE_NO_ILU_PRECONDITIONER
#include "linalg/iluprecond.hh"
#endif

#ifndef KASKADE_NO_ADDITIVESCHWARZ_PRECONDITIONER
#include "linalg/additiveschwarz.hh"
#endif

#ifndef KASKADE_NO_HYPRE_PRECONDITIONER
#include "linalg/hyprecond.hh"
#endif

#ifndef KASKADE_NO_DIRECT_PRECONDITIONER
#include "linalg/umfpack_solve.hh"
#ifndef KASKADE_NO_SUPERLU_SOLVER
#include "linalg/superlu_solve.hh"
#endif
#ifndef KASKADE_NO_MUMPS_SOLVER
#include "linalg/mumps_solve.hh"
#endif
#include "linalg/directPreconditioner.hh"
#endif

#include "utilities/detailed_exception.hh"
#include "utilities/enums.hh"
#include "utilities/factory2.hh"

/////////// MACROS ///////////////////////////////////////////////////////////////////////////////////////

#define KASKADE_CREATE_ITERATIVE_SOLVER(Name,Constructor) \
    if( parameter.preconditioner == Name) return new Solver(A, *static_cast<Constructor<Operator>* >(preconditioner), parameter.iteEps, parameter.iteSteps, parameter.verbose);

#define KASKADE_REGISTER_PRECONDITIONER(Name) \
    registered = factory.add(Name, ObjectCreation::BOOST_PP_CAT(Name,Creator)) && parameterFactory.add(Name, ObjectCreation::BOOST_PP_CAT(Name,ParameterCreator)); \
    print(#Name, registered, verbose);

//////////////////////////////////////////////////////////////////////////////////////////////////////////

namespace Kaskade
{
  struct PreconditionerParameter { virtual ~PreconditionerParameter(){} };

#ifndef KASKADE_NO_DIRECT_PRECONDITIONER
  struct DirectParameter : public PreconditionerParameter
  {

    explicit DirectParameter(DirectType solverName_=DirectType::DirectType::MUMPS, MatrixProperties property_=MatrixProperties::MatrixProperties::GENERAL) : solverName(solverName_), property(property_){}

    DirectType solverName;
    MatrixProperties property;
  };
#endif /* PrecondType::DIRECT */

#ifndef KASKADE_NO_TRIVIAL_PRECONDITIONER
  struct TrivialParameter : public PreconditionerParameter{};
#endif /* PrecondType::NONE */

#ifndef KASKADE_NO_ILU_PRECONDITIONER
  struct ILUTParameter : public PreconditionerParameter
  {
    explicit ILUTParameter(int lfil_=140, double dropTol_=0.01) : lfil(lfil_), dropTol(dropTol_){}

    int lfil;
    double dropTol;
  };

  struct ILUKParameter : public PreconditionerParameter
  {
    explicit ILUKParameter(int fill_lev_=3) : fill_lev(fill_lev_){}

    int fill_lev;
  };

  struct ARMSParameter : public PreconditionerParameter
  {
    explicit ARMSParameter(int lfil_=140, int lev_reord_=1, double dropTol_=0.01, double tolind_=0.2)
    : lfil(lfil_), lev_reord(lev_reord_), dropTol(dropTol_), tolind(tolind_){}

    int lfil, lev_reord;
    double dropTol, tolind;
  };
#endif /* ILU */

#ifndef KASKADE_NO_ADDITIVESCHWARZ_PRECONDITIONER
  struct AdditiveSchwarzParameter : public PreconditionerParameter
  {
    AdditiveSchwarzParameter(size_t firstIdx, size_t secondIdx) : firstIndex(firstIdx), secondIndex(secondIdx){}

    size_t firstIndex, secondIndex;
  };
#endif /* PrecondType::ADDITIVESCHWARZ */

#ifndef KASKADE_NO_HYPRE_PRECONDITIONER
  struct BoomerAMGParameter : public PreconditionerParameter
  {
    explicit BoomerAMGParameter(int steps_=50, int coarsentype_=21, int interpoltype_=0, int cycletype_=1, int relaxtype_=3,
        int variant_=0, int overlap_=1, double tol_=1e-9, double strongThreshold_=0.6)
    : steps(steps_), coarsentype(coarsentype_), interpoltype(interpoltype_), cycletype(cycletype_),
      relaxtype(relaxtype_), variant(variant_), overlap(overlap_), tol(tol_), strongThreshold(strongThreshold_)
    {}

    int steps, coarsentype, interpoltype, cycletype, relaxtype, variant, overlap;
    double tol, strongThreshold;
  };

  struct EuclidParameter : public PreconditionerParameter
  {
    explicit EuclidParameter(int level_=1, double droptol_=0.01, int printlevel_=0, int bj_=0)
    : level(level_), droptol(droptol_), printlevel(printlevel_), bj(bj_){}

    int level;
    double droptol;
    int printlevel, bj;
  };
#endif /* HYPRE */

  /**
   * \cond internals
   */
  namespace ObjectCreation{
    template <class Solver, class Operator, class Parameter>
    Solver* CreateIterativeLinearSolver(Operator &A, Dune::Preconditioner<typename Operator::domain_type, typename Operator::range_type>* preconditioner, Parameter const& parameter)
    {
#ifndef KASKADE_NO_TRIVIAL_PRECONDITIONER
      KASKADE_CREATE_ITERATIVE_SOLVER(PrecondType::NONE, TrivialPreconditioner)
#endif /* TRIVIAL */
#ifndef KASKADE_NO_DIRECT_PRECONDITIONER
      KASKADE_CREATE_ITERATIVE_SOLVER(PrecondType::DIRECT, DirectPreconditioner)
#endif

#ifndef KASKADE_NO_ILU_PRECONDITIONER
      KASKADE_CREATE_ITERATIVE_SOLVER(PrecondType::ILUT, ILUTPreconditioner)
      KASKADE_CREATE_ITERATIVE_SOLVER(PrecondType::ILUK, ILUKPreconditioner)
      KASKADE_CREATE_ITERATIVE_SOLVER(PrecondType::ARMS, ARMSPreconditioner)
#endif /* ILU */
#ifndef KASKADE_NO_ADDITIVESCHWARZ_PRECONDITIONER
      KASKADE_CREATE_ITERATIVE_SOLVER(PrecondType::ADDITIVESCHWARZ, AdditiveSchwarzPreconditioner)
#endif /* PrecondType::ADDITIVESCHWARZ */
#ifndef KASKADE_NO_HYPRE_PRECONDITIONER
      KASKADE_CREATE_ITERATIVE_SOLVER(PrecondType::BOOMERAMG, BoomerAMG)
      KASKADE_CREATE_ITERATIVE_SOLVER(PrecondType::EUCLID, Euclid)
#endif /* HYPRE */
      throw std::runtime_error("Could not create iterative solver. No preconditioner found.\n");
    }


    /******************************** PRECONDITIONER ********************************/

    /********************************************************************************/
    /*********************************** TRIVIAL ************************************/
    /********************************************************************************/
#ifndef KASKADE_NO_TRIVIAL_PRECONDITIONER
    struct CreateTrivialPreconditionerParameter{
      TrivialParameter* operator()() const { return new TrivialParameter; }
    } NONEParameterCreator;

    struct CreateTrivialPreconditioner
    {
      template <class Operator>
      TrivialPreconditioner<Operator>* operator()(Operator&, PreconditionerParameter const&) const {
        return new TrivialPreconditioner<Operator>;
      }
    } NONECreator;
#endif /* TRIVIAL */
    /********************************************************************************/
    /************************************ PrecondType::DIRECT ************************************/
    /********************************************************************************/
#ifndef KASKADE_NO_DIRECT_PRECONDITIONER
    struct CreateDIRECTPreconditionerParameter{
      DirectParameter* operator()() const { return new DirectParameter; }
    } DIRECTParameterCreator;

    struct CreateDIRECTPreconditioner
    {
      template <class Operator>
      DirectPreconditioner<Operator>* operator()(Operator& A, PreconditionerParameter const& param) const {
        DirectParameter const& parameter = dynamic_cast<DirectParameter const&>(param);
        return new DirectPreconditioner<Operator>(A, parameter.solverName, parameter.property);
      }
    } DIRECTCreator;
#endif /* PrecondType::DIRECT */

    /********************************************************************************/
    /************************************* PrecondType::ILUT *************************************/
    /********************************************************************************/
#ifndef KASKADE_NO_ILU_PRECONDITIONER
    struct CreateILUTPreconditionerParameter{
      ILUTParameter* operator()() const { return new ILUTParameter; }
    } ILUTParameterCreator;

    struct CreateILUTPreconditioner
    {
      template <class Operator>
      ILUTPreconditioner<Operator>* operator()(Operator& A, PreconditionerParameter const& param) const {
        ILUTParameter const& parameter = dynamic_cast<ILUTParameter const&>(param);
        return new ILUTPreconditioner<Operator>(A, parameter.lfil, parameter.dropTol);
      }
    } ILUTCreator;
    /********************************************************************************/
    /************************************* PrecondType::ILUK *************************************/
    /********************************************************************************/
    struct CreateILUKPreconditionerParameter{
      ILUKParameter* operator()() const { return new ILUKParameter; }
    } ILUKParameterCreator;

    struct CreateILUKPreconditioner
    {
      template <class Operator>
      ILUKPreconditioner<Operator>* operator()(Operator& A, PreconditionerParameter const& param) const {
        ILUKParameter const& parameter = dynamic_cast<ILUKParameter const&>(param);
        return new ILUKPreconditioner<Operator>(A, parameter.fill_lev);
      }
    } ILUKCreator;

    /********************************************************************************/
    /************************************* PrecondType::ARMS *************************************/
    /********************************************************************************/
    struct CreateARMSPreconditionerParameter{
      ARMSParameter* operator()() const { return new ARMSParameter; }
    } ARMSParameterCreator;

    struct CreateARMSPreconditioner{
      template <class Operator>
      ARMSPreconditioner<Operator>* operator()(Operator& A, PreconditionerParameter const& param) const
      {
        ARMSParameter const& parameter = dynamic_cast<ARMSParameter const&>(param);
        return new ARMSPreconditioner<Operator>(A, parameter.lfil, parameter.lev_reord, parameter.dropTol, parameter.tolind);
      }
    } ARMSCreator;
#endif /* ILU */
    /********************************************************************************/
    /******************************* ADDITIVE SCHWARZ *******************************/
    /********************************************************************************/
#ifndef KASKADE_NO_ADDITIVESCHWARZ_PRECONDITIONER
    struct CreateADDITIVESCHWARZPreconditionerParameter{
      AdditiveSchwarzParameter* operator()() const { return new AdditiveSchwarzParameter(0,1); }
    } ADDITIVESCHWARZParameterCreator;

    struct CreateADDITIVESCHWARZPreconditioner{
      template <class Operator>
      AdditiveSchwarzPreconditioner<Operator>* operator()(Operator& A, PreconditionerParameter const& param) const
      {
        AdditiveSchwarzParameter const& parameter = dynamic_cast<AdditiveSchwarzParameter const&>(param);
        return new AdditiveSchwarzPreconditioner<Operator>(A, parameter.firstIndex, parameter.secondIndex);
      }
    } ADDITIVESCHWARZCreator;
#endif /* PrecondType::ADDITIVESCHWARZ */

    /********************************************************************************/
    /********************************** Boomer AMG **********************************/
    /********************************************************************************/
#ifndef KASKADE_NO_HYPRE_PRECONDITIONER
    struct CreateBoomerAMGPreconditionerParameter{
      BoomerAMGParameter* operator()() const { return new BoomerAMGParameter; }
    } BOOMERAMGParameterCreator;

    struct CreateBoomerAMGPreconditioner{

      template <class Operator>
      BoomerAMG<Operator>* operator()(Operator& A, PreconditionerParameter const& param) const {
        BoomerAMGParameter const& parameter = dynamic_cast<BoomerAMGParameter const&>(param);
        return new BoomerAMG<Operator>(A, parameter.steps, parameter.coarsentype, parameter.interpoltype,
            parameter.tol, parameter.cycletype, parameter.relaxtype, parameter.strongThreshold,
            parameter.variant, parameter.overlap);
      }
    } BOOMERAMGCreator;

    /********************************************************************************/
    /************************************ Euclid ************************************/
    /********************************************************************************/
    struct CreateEuclidPreconditionerParameter{
      EuclidParameter* operator()() const { return new EuclidParameter; }
    } EUCLIDParameterCreator;

    struct CreateEuclidPreconditioner{
      template <class Operator>
      Euclid<Operator>* operator()(Operator& A, PreconditionerParameter const& param) const
      {
        EuclidParameter const& parameter = dynamic_cast<EuclidParameter const&>(param);
        return new Euclid<Operator>(A, parameter.level, parameter.droptol, parameter.printlevel, parameter.bj);
      }
    } EUCLIDCreator;
#endif /* HYPRE */
  }
  /**
   * \endcond
   */

  /// Preconditioner factory
  /**
   * Preconditioner identifier:
   *   PrecondType::NONE
   *   PrecondType::DIRECT
   *   PrecondType::ILUT
   *   PrecondType::ILUK
   *   PrecondType::ADDITIVESCHWARZ
   *   PrecondType::ARMS
   *   PrecondType::BOOMERAMG
   *   PrecondType::EUCLID
   *
   *  If you want to add a new preconditioner you have to:
   *  1. Add a unique identifier to the enum PrecondType.
   *  2. Create a default constructible container that stores the preconditioner's parameters. This container must
   *     inherit from PreconditionerParameter.
   *  3. For both the parameter object as well as the preconditioner create a functor that is
   *     responsible for the creation of the objects and returns them through a unique pointer.
   *  4. Use the macro KASKADE_REGISTER_PRECONDITIONER to add the triple (identifier,preconditioner,parameter) to
   *     the factory (or do it by hand).
   *
   *   In order to save compile time it is possible to disable preconditioners using the following preprocessor macros:
   *   Macro:                                          Preconditioner:            Headers (in folder linalg):
   *   KASKADE_NO_DIRECT_PRECONDITIONER                PrecondType::DIRECT                     mumps_solve.hh, umfpack_solve.hh, superlu_solve.hh
   *   KASKADE_NO_TRIVIAL_PRECONDITIONER               PrecondType::NONE                       trivialpreconditioner.hh
   *   KASKADE_NO_ILU_PRECONDITIONER                   PrecondType::ILUT, PrecondType::ILUK, PrecondType::ARMS           iluprecond.hh
   *   KASKADE_NO_ADDITIVESCHWARZ_PRECONDITIONER       PrecondType::ADDITIVESCHWARZ            additiveschwarz.hh
   *   KASKADE_NO_HYPRE_PRECONDITIONER                 PrecondType::BOOMERAMG, PrecondType::EUCLID          hyprecond.hh
   */
  template <class Assembler, class Domain, class Range, class Scalar=typename Domain::field_type>
  class PreconditionerFactory{
  public:
    typedef ApplyOperator<Scalar,Domain,Range,Assembler> Operator;
    typedef Factory2<PrecondType, Dune::Preconditioner<Domain,Range>, Operator&, PreconditionerParameter const&> FactoryImpl;
    typedef Factory2<PrecondType, PreconditionerParameter> ParameterFactoryImpl;

  public:
    /** \typedef base type for preconditioners */
    typedef Dune::Preconditioner<Domain, Range> PreconditionerType;

    explicit PreconditionerFactory(bool verbose_=false)
    : verbose(verbose_), factory(), parameterFactory(), parameter(nullptr), identifier(PrecondType::NONE)
    {
      bool registered;
      if(verbose) std::cout << "Registering preconditioners:" << std::endl;

#ifndef KASKADE_NO_TRIVIAL_PRECONDITIONER
      KASKADE_REGISTER_PRECONDITIONER(PrecondType::NONE)
#endif /* TRIVIAL */

#ifndef KASKADE_NO_DIRECT_PRECONDITIONER
      KASKADE_REGISTER_PRECONDITIONER(PrecondType::DIRECT)
#endif /* PrecondType::DIRECT */

#ifndef KASKADE_NO_ILU_PRECONDITIONER
      KASKADE_REGISTER_PRECONDITIONER(PrecondType::ILUT)
      KASKADE_REGISTER_PRECONDITIONER(PrecondType::ILUK)
      KASKADE_REGISTER_PRECONDITIONER(PrecondType::ARMS)
#endif /* ILU */

#ifndef KASKADE_NO_ADDITIVESCHWARZ_PRECONDITIONER
      KASKADE_REGISTER_PRECONDITIONER(PrecondType::ADDITIVESCHWARZ)
#endif /* PrecondType::ADDITIVESCHWARZ */

#ifndef KASKADE_NO_HYPRE_PRECONDITIONER
      KASKADE_REGISTER_PRECONDITIONER(PrecondType::BOOMERAMG)
      KASKADE_REGISTER_PRECONDITIONER(PrecondType::EUCLID)
#endif /* HYPRE */

      if(verbose) std::cout << std::endl;
    }

    explicit PreconditionerFactory(PrecondType identifier_, bool verbose_ = false)
     : PreconditionerFactory(verbose)
    {
      setPreconditioner(identifier_);
    }

    template <class Parameter>
    PreconditionerFactory(PrecondType identifier_, Parameter const& parameter_ , bool verbose_=false)
     : PreconditionerFactory(verbose_)
    {
      setParameter(identifier_,parameter_);
    }

    template <class Parameter>
    void setParameter(PrecondType identifier_, Parameter const& parameter_)
    {
      parameter.reset(new Parameter(parameter_));
      identifier = identifier_;
    }

    /// Also resets the parameter object to default values.
    void setPreconditioner(PrecondType identifier_)
    {
      identifier = identifier_;
      setParameter();
    }

    /// Create Preconditioner.
    /*
     * Throws if identifier has not been registered in the constructor.
     *
     * \param ReturnType type of the desired preconditioner, actual return type is std::unique_ptr<ReturnType> (optional, default is std::unique_ptr<Dune::Preconditioner<Domain,Range> >)
     *
     * \param identifier one of above mentioned identifiers
     * \param A operator (must inherit from MatrixRepresentedOperator
     * \param parameter (must inherit from PreconditionerParameter
     * \param verbose
     * \return preconditioner
     */
    template <class Preconditioner=Dune::Preconditioner<Domain,Range> >
    std::unique_ptr<Preconditioner> create(PrecondType identifier_, Operator A, PreconditionerParameter const& parameters) const
    {
      return std::unique_ptr<Preconditioner>( static_cast<Preconditioner*>( factory.create( identifier, A, parameters ).release() ) );
    }

    /// Create Preconditioner using default parameters.
    /*
     * Throws if identifier has not been registered in the constructor.
     *
     * \param ReturnType type of the desired preconditioner, actual return type is std::unique_ptr<ReturnType> (optional, default is std::unique_ptr<Dune::Preconditioner<Domain,Range> >)
     *
     * \param identifier one of above mentioned identifiers
     * \param A operator (must inherit from MatrixRepresentedOperator
     * \param parameter (must inherit from PreconditionerParameter
     * \param verbose
     * \return preconditioner
     */
    template <class Preconditioner=Dune::Preconditioner<Domain,Range> >
    std::unique_ptr<Preconditioner> create(PrecondType identifier_, Operator A)
    {
      setPreconditioner(identifier_);
      return create<Preconditioner>(A);
    }

    /// Create Preconditioner using stored parameters and identifier (resp. default parameters if no parameters are stored).
    /*
     * Throws if identifier has not been registered in the constructor or no identifier has been set.
     *
     * \param ReturnType type of the desired preconditioner, actual return type is std::unique_ptr<ReturnType> (optional, default is std::unique_ptr<Dune::Preconditioner<Domain,Range> >)
     *
     * \param identifier one of above mentioned identifiers
     * \param A operator (must inherit from MatrixRepresentedOperator
     * \param parameter (must inherit from PreconditionerParameter
     * \param verbose
     * \return preconditioner
     */
    template <class Preconditioner=Dune::Preconditioner<Domain,Range> >
    std::unique_ptr<Preconditioner> create(Operator A) const
    {
      return std::unique_ptr<Preconditioner>( static_cast<Preconditioner*>( factory.create( identifier, A, *parameter ).release() ) );
    }

    /// Get preconditioner's parameter object.
    template <class Parameter=PreconditionerParameter>
    Parameter get(PrecondType identifier_=PrecondType::NONE)
    {
      if(identifier_ != identifier && identifier_ != PrecondType::NONE) setPreconditioner(identifier_);
      else if(parameter == nullptr) setParameter();
      return Parameter( *( static_cast<Parameter*>(parameter.get()) ) );
    }

  private:
    void setParameter()
    {
      parameter = parameterFactory.create(identifier);
    }

    /// Print information
    void print(std::string name, int registered, bool verbose) const
    {
      if(!verbose) return;
      std::cout << name << ": ";
      if(registered) std::cout << "successful" << std::endl;
      else std::cout << "failed" << std::endl;
    }

    bool verbose;
    mutable FactoryImpl factory;
    mutable ParameterFactoryImpl parameterFactory;
    std::unique_ptr<PreconditionerParameter> parameter;
    PrecondType identifier;
  };
}

#undef KASKADE_CREATE_ITERATIVE_SOLVER
#undef KASKADE_REGISTER_PRECONDITIONER // undef macros used for code generation

#endif /* PRECONDITIONER_FACTORY_HH_ */
