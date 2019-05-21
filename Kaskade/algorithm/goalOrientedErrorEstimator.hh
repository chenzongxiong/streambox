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

#ifndef GOAL_ORIENTED_ERROR_ESTIMATOR_HH
#define GOAL_ORIENTED_ERROR_ESTIMATOR_HH

#include <cmath>
#include <iostream>
#include <type_traits>
#include <utility>

#include <boost/fusion/include/at_c.hpp>
#include <boost/mpl/at.hpp>

#include "algorithm/newton_bridge.hh"
#include "algorithm/errorDistribution.hh"
#include "linalg/jacobiPreconditioner.hh"
#include "fem/hierarchicErrorEstimator.hh"
#include "fem/variables.hh"
#include "linalg/direct.hh"
#include "utilities/enums.hh"

namespace Kaskade
{
  bool biggerThanAbs(const std::pair<double,int>& x1, const std::pair<double,int>&  x2)
  {
    return fabs(x1.first)>fabs(x2.first);
  }

  double add(const double& x1, const std::pair<double,int>&  x2)
  {
    return x1+x2.first;
  }

  template <template <class,class> class TemplateFunctional, class OriginalVariableSetDescription, class ExtensionVariableSetDescription, class ExtensionSpace>
  class GoalOrientedErrorEstimator : public AbstractHierarchicalErrorEstimator
  {
    static constexpr int noOfVariables = ExtensionVariableSetDescription::noOfVariables;

    // define functional types
    typedef TemplateFunctional<OriginalVariableSetDescription,OriginalVariableSetDescription> F_hh;
    typedef TemplateFunctional<OriginalVariableSetDescription,ExtensionVariableSetDescription> F_he;
    typedef TemplateFunctional<ExtensionVariableSetDescription,OriginalVariableSetDescription> F_eh;
    typedef TemplateFunctional<ExtensionVariableSetDescription,ExtensionVariableSetDescription> F_ee;

    static constexpr int yIdx = F_hh::yIdx;
    static constexpr int uIdx = F_hh::uIdx;
    static constexpr int pIdx = F_hh::lIdx;

    // define an assembler for each functional
    typedef VariationalFunctionalAssembler<LinearizationAt<F_hh> > Assembler_hh;
    typedef VariationalFunctionalAssembler<LinearizationAt<F_eh> > Assembler_eh;
    typedef VariationalFunctionalAssembler<LinearizationAt<F_he> > Assembler_he;
    typedef VariationalFunctionalAssembler<LinearizationAt<F_ee> > Assembler_ee;

    typedef typename OriginalVariableSetDescription::VariableSet VariableSet;
    typedef typename ExtensionVariableSetDescription::VariableSet ExtensionVariableSet;

    typedef typename ExtensionVariableSetDescription::GridView::template Codim<0>::Iterator CellIterator;

  public:
    typedef typename F_hh::Scalar Scalar;
    static constexpr int dim = VariableSet::Descriptions::Grid::dimension;

    // error estimator and corresponding assembler
    typedef HierarchicErrorEstimator<LinearizationAt<F_hh>,ExtensionVariableSetDescription,ExtensionVariableSetDescription,HierarchicErrorEstimatorDetail::TakeAllD2<LinearizationAt<F_hh> > > ErrorEstimator;
    typedef VariationalFunctionalAssembler<ErrorEstimator> ExtensionAssembler;

    // coefficient vectors
    typedef typename OriginalVariableSetDescription::template CoefficientVectorRepresentation<>::type Vector;
    typedef typename OriginalVariableSetDescription::template CoefficientVectorRepresentation<> VectorCreator;
    typedef typename OriginalVariableSetDescription::template CoefficientVectorRepresentation<uIdx,uIdx+1>::type UVector;
    typedef typename OriginalVariableSetDescription::template CoefficientVectorRepresentation<yIdx,yIdx+1>::type YVector;
    typedef typename OriginalVariableSetDescription::template CoefficientVectorRepresentation<uIdx,uIdx+1> UVectorCreator;
    typedef typename OriginalVariableSetDescription::template CoefficientVectorRepresentation<yIdx,yIdx+1> YVectorCreator;

    // coefficient vectors in extended space
    typedef typename ExtensionVariableSetDescription::template CoefficientVectorRepresentation<>::type ExtensionVector;
    typedef typename ExtensionVariableSetDescription::template CoefficientVectorRepresentation<uIdx,uIdx+1>::type UExtensionVector;
    typedef typename ExtensionVariableSetDescription::template CoefficientVectorRepresentation<yIdx,yIdx+1>::type YExtensionVector;
    typedef typename ExtensionVariableSetDescription::template CoefficientVectorRepresentation<> ExtensionVectorCreator;
    typedef typename ExtensionVariableSetDescription::template CoefficientVectorRepresentation<uIdx,uIdx+1> UExtensionVectorCreator;
    typedef typename ExtensionVariableSetDescription::template CoefficientVectorRepresentation<yIdx,yIdx+1> YExtensionVectorCreator;

    GoalOrientedErrorEstimator(ExtensionVariableSetDescription& extensionVariableSetDescription_, ExtensionSpace& extensionSpace_, Scalar fraction_=0.8, size_t maxNoOfCells_=5000, bool onlyLowerTriangle_=false)
    : extensionVariableSetDescription(extensionVariableSetDescription_), extensionSpace(extensionSpace_), extensionAssembler(extensionVariableSetDescription.spaces),
      fraction(fraction_), maxNoOfCells(maxNoOfCells_), onlyLowerTriangle(onlyLowerTriangle_)
    {}

    ~GoalOrientedErrorEstimator(){}

    void operator()(AbstractVector const& x_, AbstractVector const& dx_, int step)
    {
      if(extensionSpace.gridManager().grid().size(0) > maxNoOfCells)
      {
        std::cout << "Skipping error estimation." << std::endl;
        return;
      }

      std::cout << "Starting error estimation." << std::endl;

      Bridge::Vector<VariableSet> const& x = dynamic_cast<const Bridge::Vector<VariableSet>&>(x_);
      Bridge::Vector<VariableSet> const& dx = dynamic_cast<const Bridge::Vector<VariableSet>&>(dx_);

      //OriginalVariableSetDescription const& variableSetDescription = x.get().descriptions;

      using namespace boost::fusion;

      MatrixProperties property = MatrixProperties::GENERAL;
      DirectType directType = DirectType::UMFPACK;

      std::cout << "computing residual in extension space." << std::endl;
      //ExtensionAssembler extensionAssembler(extensionVariableSetDescription.spaces);
      extensionAssembler.assemble(ErrorEstimator(LinearizationAt<F_hh>(*f_hh,x.get()),dx.get()));
      AssembledGalerkinOperator<ExtensionAssembler,0,noOfVariables,0,noOfVariables> H_ee(extensionAssembler); // in space extension
      ExtensionVector rhs_e(extensionAssembler.rhs());
      rhs_e *= -1;
      ExtensionVector dx_e(ExtensionVectorCreator::init(extensionVariableSetDescription));
      // compute correction in extension space
      directInverseOperator(H_ee,directType,property).apply(rhs_e,dx_e);
      // create all necessary assemblers
      Assembler_hh assembler_hh(x.get().descriptions.spaces);
      Assembler_ee assembler_ee(extensionVariableSetDescription.spaces);
      Assembler_eh assembler_eh(extensionVariableSetDescription.spaces);
      Assembler_he assembler_he(extensionVariableSetDescription.spaces);

      assembler_hh.assemble(linearization(*f_hh,x.get()));
      assembler_he.assemble(linearization(*f_he,x.get()));
      ExtensionVariableSet x_e(extensionVariableSetDescription);
      x_e *= 0;
      assembler_eh.assemble(linearization(*f_eh,x_e));
      assembler_ee.assemble(linearization(*f_ee,x_e));

      /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
      std::cout << "computing weight functions 1" << std::endl;
      // compute weight function for state variable
      // 1. get rhs
      YVector tmpRhs_h(YVectorCreator::init(x.get().descriptions));
      YExtensionVector tmpRhs_e(YExtensionVectorCreator::init(extensionVariableSetDescription));

      UVector du_h( at_c<uIdx>(dx.get().data).coefficients() );
      UExtensionVector du_e( at_c<uIdx>(dx_e.data) );

      {
        AssembledGalerkinOperator<Assembler_hh,pIdx,pIdx+1,uIdx,uIdx+1> B_hh(assembler_hh, false);
        AssembledGalerkinOperator<Assembler_he,pIdx,pIdx+1,uIdx,uIdx+1> B_he(assembler_he, false);
        AssembledGalerkinOperator<Assembler_eh,pIdx,pIdx+1,uIdx,uIdx+1> B_eh(assembler_eh, false);
        AssembledGalerkinOperator<Assembler_ee,pIdx,pIdx+1,uIdx,uIdx+1> B_ee(assembler_ee, false);

        B_hh.apply(du_h, tmpRhs_h);
        B_eh.applyscaleadd(1.0, du_e, tmpRhs_h);
        tmpRhs_h *= -1;
        B_he.apply(du_h, tmpRhs_e);
        B_ee.applyscaleadd(1.0, du_e, tmpRhs_e);
        tmpRhs_e *= -1;
      }

      // 2. compute weight
      AssembledGalerkinOperator<Assembler_hh,pIdx,pIdx+1,yIdx,yIdx+1> A_hh(assembler_hh, false);
      AssembledGalerkinOperator<Assembler_eh,pIdx,pIdx+1,yIdx,yIdx+1> A_eh(assembler_eh, false);
      AssembledGalerkinOperator<Assembler_ee,pIdx,pIdx+1,yIdx,yIdx+1> A_ee(assembler_ee, false);
      AssembledGalerkinOperator<Assembler_he,pIdx,pIdx+1,yIdx,yIdx+1> A_he(assembler_he, false);
      YVector wy_h(YVectorCreator::init(x.get().descriptions));
      wy_h *= 0;
      YExtensionVector wy_e(YExtensionVectorCreator::init(extensionVariableSetDescription));

      //JacobiPreconditioner<AssembledGalerkinOperator<Assembler_hh,pIdx,pIdx+1,yIdx,yIdx+1> >(A_hh).apply(tmpRhs_h,wy_h);
      //A_hh.apply(tmpRhs_h, wy_h);
      InverseLinearOperator<DirectSolver<YVector,YVector> > ds_hh(DirectSolver<YVector,YVector>(A_hh,DirectType::MUMPS,MatrixProperties::SYMMETRIC));
      ds_hh.apply(tmpRhs_h, wy_h);
      A_he.applyscaleadd(-1.0,wy_h,tmpRhs_e);
      //JacobiPreconditioner<AssembledGalerkinOperator<Assembler_ee,pIdx,pIdx+1,yIdx,yIdx+1> >(A_ee).apply(tmpRhs_e,wy_e);
      directInverseOperator(A_ee,directType,property).apply(tmpRhs_e, wy_e);

      /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
      std::cout << "computing weight functions 2" << std::endl;
      // compute weight function for adjoint variable
      // 1. get rhs
      {
        AssembledGalerkinOperator<Assembler_hh,yIdx,yIdx+1,yIdx,yIdx+1> H_hh(assembler_hh, false);
        AssembledGalerkinOperator<Assembler_he,yIdx,yIdx+1,yIdx,yIdx+1> H_he(assembler_he, false);
        AssembledGalerkinOperator<Assembler_eh,yIdx,yIdx+1,yIdx,yIdx+1> H_eh(assembler_eh, false);
        AssembledGalerkinOperator<Assembler_ee,yIdx,yIdx+1,yIdx,yIdx+1> H_ee(assembler_ee, false);
        H_hh.apply(wy_h,tmpRhs_h);
        H_eh.applyscaleadd(1.0,wy_e,tmpRhs_h);
        tmpRhs_h *= -1;
        H_he.apply(wy_h,tmpRhs_e);
        H_ee.applyscaleadd(1.0,wy_e,tmpRhs_e);
        tmpRhs_e *= -1;
      }

      // 2. compute weight
      YVector wp_h(YVectorCreator::init(x.get().descriptions));
      YExtensionVector wp_e(YExtensionVectorCreator::init(extensionVariableSetDescription));
      {
        //JacobiPreconditioner<AssembledGalerkinOperator<Assembler_hh,pIdx,pIdx+1,yIdx,yIdx+1> >(A_hh).apply(tmpRhs_h,wp_h);
        directInverseOperator(A_hh,directType,property).apply(tmpRhs_h, wp_h);
        A_he.applyscaleadd(-1.0,wp_h, tmpRhs_e);
        //JacobiPreconditioner<AssembledGalerkinOperator<Assembler_ee,pIdx,pIdx+1,yIdx,yIdx+1> >(A_ee).apply(tmpRhs_e,wp_e);
        directInverseOperator(A_ee,directType,property).apply(tmpRhs_e, wp_e);
      }

      Scalar estimatedError = 0.5 * ( at_c<0>(wy_e.data)*at_c<yIdx>(rhs_e.data)
          + at_c<0>(du_e.data)*at_c<uIdx>(rhs_e.data)
          + at_c<0>(wp_e.data)*at_c<pIdx>(rhs_e.data) );

      std::cout << "estimated error: " << estimatedError << std::endl;

      // get number of cells corresponding to edges resp. edge shape functions used in the error estimator
      CellIterator cend = extensionVariableSetDescription.gridView.template end<0>();
      std::vector<size_t> neighBoringCells(extensionVariableSetDescription.gridView.size(dim-1),0);
      for(CellIterator ci = extensionVariableSetDescription.gridView.template begin<0>(); ci!=cend; ++ci)
      {
        // iterate over cells edges
        for(int edgeIdInEntity=0; edgeIdInEntity<3; ++edgeIdInEntity)
        {
          int id = extensionVariableSetDescription.gridView.indexSet().subIndex(*ci, edgeIdInEntity, dim-1);
          if(id < neighBoringCells.size())
            neighBoringCells[id]++;
        }
      }
//      for(size_t i=0; i<neighBoringCells.size(); ++i) if(neighBoringCells[i]==0) std::cout << "missing cells at edge id: " << i << std::endl;
//
      // estimate cell-wise error distribution
      typename ExtensionVariableSetDescription::VariableSet mySol(extensionVariableSetDescription);
      mySol *= 0;
      std::vector<Scalar> edgeErrorIndicators(at_c<yIdx>(dx_e.data).N());
      std::cout << "edgeErrorIndicators: " << edgeErrorIndicators.size() << std::endl;
      std::cout << "mySol: " << boost::fusion::at_c<1>(mySol.data).coefficients().N() << std::endl;
      for(size_t i=0; i<edgeErrorIndicators.size(); ++i)
      {
        edgeErrorIndicators[i] = 0.5 * ( at_c<0>(wy_e.data)[i]*at_c<yIdx>(rhs_e.data)[i] + at_c<0>(du_e.data)[0]*at_c<uIdx>(rhs_e.data)[0] + at_c<0>(wp_e.data)[i]*at_c<pIdx>(rhs_e.data)[i] );
        boost::fusion::at_c<1>(mySol.data).coefficients()[i][0] = edgeErrorIndicators[i];
        edgeErrorIndicators[i] /= neighBoringCells[i];
      }

      errorDistribution.clear();
      errorDistribution.resize(extensionVariableSetDescription.gridView.size(0),std::make_pair(0.0,0));
      auto const& is = extensionSpace.gridManager().grid().leafIndexSet();
//
//      typedef ErrorDistribution<F_hh,ExtensionVariableSetDescription> EnergyError;
//      EnergyError energyError(*f_hh,x.get(),mySol);
//      typedef VariationalFunctionalAssembler<LinearizationAt<EnergyError> > EnergyErrorAssembler;
//      EnergyErrorAssembler eeAssembler(extensionSpace.gridManager(), energyError.getSpaces());
//
//      eeAssembler.assemble(linearization(energyError,x.get()), EnergyErrorAssembler::RHS);
//      typename EnergyError::ErrorVector distError( eeAssembler.rhs() );
//      //if(!nearMinimizer_) distError *= 0;
//      typename EnergyError::AnsatzVars::VariableSet mde( energyError.getVariableSetDescription() );
//      mde = distError;

      for (CellIterator ci=extensionVariableSetDescription.gridView.template begin<0>(); ci!=cend; ++ci)
//        errorDistribution[is.index(*ci)] = std::make_pair( fabs(boost::fusion::at_c<0>(mde.data).value(*ci,Dune::FieldVector<Scalar,dim>(0.3))) , is.index(*ci));
      {
        size_t index = is.index(*ci);
        errorDistribution[index] = std::make_pair(0,index);
        size_t nEdges = (dim==2) ? 3 : 6;
        for(int edgeIdInEntity=0; edgeIdInEntity<nEdges; ++edgeIdInEntity)
        {
          typedef typename ExtensionVariableSetDescription::GridView::Grid Grid;
//          std::cout << "edgeInInEntity: " << edgeIdInEntity << std::endl;
          size_t globalEdgeId = Dune::GenericReferenceElements<typename Grid::ctype,Grid::dimension>::simplex().subEntity(edgeIdInEntity,Grid::dimension-1,0,Grid::dimension);
          //size_t globalEdgeId = is.subIndex(*ci, edgeIdInEntity, dim-1);
          errorDistribution[index].first += edgeErrorIndicators[globalEdgeId];
        }
     }

      totalError = std::accumulate(errorDistribution.begin(), errorDistribution.end(), 0.0, add);
      //      typedef FEFunctionSpace<DiscontinuousLagrangeMapper<Scalar,typename ExtensionSpace::GridView> > AnsatzSpace;
      //      typedef boost::fusion::vector<AnsatzSpace const*> AnsatzSpaces;
      //      typedef Variable<SpaceIndex<0>,Components<1>,VariableId<0> > AnsatzVariableInformation;
      //      typedef boost::fusion::vector<AnsatzVariableInformation> VariableDescriptions;
      //      typedef VariableSetDescription<AnsatzSpaces,VariableDescriptions> AnsatzVars;
      //      AnsatzSpace discontinuousSpace(extensionSpace.gridManager(),extensionSpace.gridManager().grid().leafView(),0);
      //      AnsatzSpaces discontinuousSpaces(&discontinuousSpace);
      //      std::string dname[] = { "error" };
      //      AnsatzVars errorRepVarSetDesc(discontinuousSpaces,dname);
      //      typename AnsatzVars::VariableSet discontinuousErrorRepresentation(errorRepVarSetDesc);
      //
      //      std::cout << "before transfer" << std::endl;
      //      spaceTransfer(boost::fusion::at_c<0>(discontinuousErrorRepresentation.data),
      //                    boost::fusion::at_c<0>(edgeErrorElement.data));
      //      std::cout << "after" << std::endl;
      //
      ////     typename AnsatzSpace::template Element<1> errorRep(discontinuousSpace);
      //
      //      for(CellIterator ci = extensionVariableSetDescription.gridView.template begin<0>(); ci!=cend; ++ci)
      //      {
      ////        typedef typename std::remove_pointer<typename boost::mpl::at_c<typename OriginalVariableSetDescription::Spaces,0>::type>::type Space;
      ////        typedef typename Space::Mapper::GlobalIndexRange IndexMapper;
      //        typedef typename ExtensionSpace::Mapper::GlobalIndexRange IndexMapper;
      //        Scalar err = 0;
      //        IndexMapper mapper = extensionSpace.mapper().globalIndices(*ci);
      //        typename IndexMapper::iterator iend = mapper.end();
      //        size_t index = extensionVariableSetDescription.gridView.indexSet().index(*ci);
      //        std::cout << "index: " << index << std::endl;
      //        if(index==2)
      //        {
      //          for(typename IndexMapper::iterator iter=mapper.begin(); iter!=iend; ++iter)
      //          {
      //            std::cout << "it: " << *iter << std::endl;
      //            std::cout << "err: " << err << std::endl;
      //            std::cout << "eei: " << edgeErrorIndicators[*iter] << std::endl;
      ////            std::cout << "nbc: " << neighBoringCells[*iter] << std::endl;
      //
      //            err += edgeErrorIndicators[*iter];//neighBoringCells[*iter];
      //          }
      //        }
      //        else
      //        for(typename IndexMapper::iterator iter=mapper.begin(); iter!=iend; ++iter) err += edgeErrorIndicators[*iter];///neighBoringCells[*iter];
      //
      ////        errorDistribution[extensionVariableSetDescription.gridView.indexSet().index(*ci)] = err;
      //
      //        errorDistribution[index] = std::make_pair(err,index);
      //        if(std::isinf(err))
      //        {
      //          std::cout << "err: inf at index " << index << std::endl;
      //          exit(-1);
      //        }
      //        if(std::isnan(err))
      //        {
      //          std::cout << "err: nan at index " << index << std::endl;
      //          exit(-1);
      //        }
      //      }
      //      std::cout << "total error " << std::endl;
      //      totalError = std::accumulate(errorDistribution.begin(), errorDistribution.end(), 0.0, add);
      //      std::cout << "total error done" << std::endl;
      ////            std::vector<Scalar> errorDistribution(variableSetDescription.gridView.size(0),0);
      ////            for(CellIterator ci = variableSetDescription.gridView.template begin<0>(); ci!=cend; ++ci)
      ////            {
      ////              for(int edgeIdInEntity=0; edgeIdInEntity<(dim==2) ? 3 : 6; ++edgeIdInEntity)
      ////              {
      ////                size_t globalEdgeId = variableSetDescription.gridView.indexSet().subIndex(*ci, edgeIdInEntity, dim-1);
      ////
      ////              }
      ////            }

    }

    template <typename... Args>
    void initializeFunctionals(Args... args)
    {
      f_hh.reset(new F_hh(args...));
      f_eh.reset(new F_eh(args...));
      f_he.reset(new F_he(args...));
      f_ee.reset(new F_ee(args...));
    }

    bool nearMinimizer() const { return true; }

    void refineGrid()
    {
      if(extensionSpace.gridManager().grid().size(0) > maxNoOfCells)
      {
        std::cout << "ERROR ESTIMATOR: Grid has more than " << maxNoOfCells << " cells -> No more refinement." << std::endl;
        return;
      }

      std::sort(errorDistribution.begin(), errorDistribution.end(), biggerThanAbs);

      double errLevel = errorDistribution[3*errorDistribution.size()/4].first;
      auto const& is = extensionSpace.gridManager().grid().leafIndexSet();

      CellIterator cend = extensionVariableSetDescription.gridView.template end<0>();
      for (CellIterator ci=extensionVariableSetDescription.gridView.template begin<0>(); ci!=cend; ++ci) // iterate over cells
        if(errorDistribution[is.index(*ci)].first > errLevel) extensionSpace.gridManager().mark(1,*ci);

//      Scalar bulkCriterionTolerance = fraction*totalError;
//      std::cout << "totalError: " << totalError << std::endl;
//      std::cout << "bulkCriterionTolerance: " << bulkCriterionTolerance << std::endl;
//      Scalar bulkErrorSquared = 0;
//
//      std::vector<std::pair<Scalar,size_t> > bulkErrorDistribution;
//      auto const& is = extensionSpace.gridManager().grid().leafIndexSet();
//      CellIterator cend = extensionVariableSetDescription.gridView.template end<0>();
//
//      while(bulkErrorSquared < bulkCriterionTolerance)
//      {
//        bulkErrorDistribution.push_back(errorDistribution[bulkErrorDistribution.size()]);
//        bulkErrorSquared += bulkErrorDistribution.back().first;
//      }
//
//      std::cout << "number of candidates for refinement: " << bulkErrorDistribution.size() << std::endl;
//
//
//
//      // Refine mesh.
//      for (CellIterator ci=extensionVariableSetDescription.gridView.template begin<0>(); ci!=cend; ++ci) // iterate over cells
//      {
//        for(int i=0; i<bulkErrorDistribution.size(); ++i) // iterate over chosen part of the error distribution
//          if(is.index(*ci) == bulkErrorDistribution[i].second)
//            extensionSpace.gridManager().mark(1,*ci);
//      }

      extensionSpace.gridManager().adaptAtOnce();
    }

    double estimatedAbsoluteError() const
    {
      return fabs(totalError);
    }

    size_t gridSize() const final
    {
      return extensionSpace.gridManager().grid().size(0);
    }

  private:
    //OriginalVariableSetDescription const& variableSetDescription;
    ExtensionVariableSetDescription& extensionVariableSetDescription;
    ExtensionSpace& extensionSpace;
    ExtensionAssembler extensionAssembler;
    Scalar fraction;
    Scalar totalError;
    size_t maxNoOfCells;
    bool onlyLowerTriangle;
    std::vector<std::pair<double,size_t> > errorDistribution;
    //    std::vector<Scalar> errorDistribution;
    std::unique_ptr<F_hh> f_hh;
    std::unique_ptr<F_eh> f_eh;
    std::unique_ptr<F_he> f_he;
    std::unique_ptr<F_ee> f_ee;
  };
}

#endif
