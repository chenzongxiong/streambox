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

#ifndef BOUNDARY_NORMAL_COLLECTOR_HH
#define BOUNDARY_NORMAL_COLLECTOR_HH

#include <algorithm>
#include "utilities/interpolation/tools.hh"

namespace Kaskade
{

  namespace BoundaryNormalCollectorDetail
  {

    template <class GridView, class NormalCollection>
    void checkForBoundaryVertices(GridView const& gridView, NormalCollection& collection)
    {
      std::for_each(gridView.template begin<0>(), gridView.template end<0>(), [&](typename GridView::template Codim<0>::Iterator::Entity const& entity)
      {
        Dune::GeometryType const gt(entity.type().id(), GridView::dimension);
        std::for_each(gridView.ibegin(entity), gridView.iend(entity), [&](typename GridView::IntersectionIterator::Intersection const& face)
        {
          if(face.boundary()) for(int vertexId=0; vertexId<face.geometry().corners(); ++vertexId) collection[gridView.indexSet().subIndex(entity, Dune::GenericReferenceElements<typename GridView::ctype,GridView::dimension>::general(gt).subEntity(face.indexInInside(),1,vertexId,GridView::dimension), GridView::dimension)].onBoundary = true;
        });
      });
    }

    template <class GridView, class NormalCollection>
    void checkForBoundaryEdges(GridView const& gridView, NormalCollection& collection)
    {
      std::for_each(gridView.template begin<0>(), gridView.template end<0>(), [&](typename GridView::template Codim<0>::Iterator::Entity const& entity)
      {
        Dune::GeometryType const gt(entity.type().id(), GridView::dimension);
        std::for_each(gridView.ibegin(entity), gridView.iend(entity), [&](typename GridView::IntersectionIterator::Intersection const& face)
        {
          if(face.boundary()) for(int edgeId = 0; edgeId < 3; ++edgeId)
            collection[gridView.indexSet().subIndex(entity, Dune::GenericReferenceElements<typename GridView::ctype,GridView::dimension>::general(gt).subEntity(face.indexInInside(),1,edgeId,GridView::dimension-1), GridView::dimension-1)].onBoundary = true;
        });
      });
    }

    // helper function 2D
    template <class GridView, class NormalContainer>
    inline void outerBoundaryCheckImpl(GridView const& gridView, NormalContainer& container, std::integral_constant<int,2>)
    {
      checkForBoundaryVertices(gridView, container.vertexNormals);
    }

    // helper function 3D
    template <class GridView, class NormalContainer>
    inline void outerBoundaryCheckImpl(GridView const& gridView, NormalContainer& container, std::integral_constant<int,3>)
    {
      checkForBoundaryVertices(gridView, container.vertexNormals);
      checkForBoundaryEdges(gridView, container.edgeNormals);
    }

    // emulate template specification via overloading and std::integral_constant
    template <class GridView, class NormalContainer>
    inline void outerBoundaryCheck(GridView const& gridView, NormalContainer& container)
    {
      outerBoundaryCheckImpl(gridView, container, std::integral_constant<int,GridView::dimension>());
    }

    // store normals and associate them with vertices
    template <class GridView, class Entity, class Face, class NormalCollection, class Scalar, class MergedPolicy>
    void getEntitiesVertexNormals(GridView const& gridView, Entity const& entity, Face const& face, NormalCollection &collection, Scalar tissue, Dune::GeometryType const& gt, MergedPolicy const& policy)
    {
      const int dim = GridView::dimension;
      tissue = round(tissue);
      int local_fidx = face.indexInInside();
      for(int i=0; i<face.geometry().corners(); ++i){
        // get local index in cell from local index in face
        int local_vidx = Dune::GenericReferenceElements<Scalar,dim>::general(gt).subEntity(local_fidx,1,i,dim);
        // get global index from local index in cell
        int vidx = gridView.indexSet().subIndex(entity, local_vidx, dim);

        Scalar tissueNeighbour = -1;

        // boundary faces -> outer boundary policy
        if( face.boundary() ){
          if( policy.ignoreOuterBoundary(face) ){
            collection[vidx].isRelevant = false;
            continue;
          }
        } // if not on boundary:
        else{
          // face not on boundary, but contains boundary vertex
          // ignore contribution from the inside of the domain
          if(collection[vidx].onBoundary) continue;

          // get tissue of neighbouring cell sharing the intersection
          tissueNeighbour = round( policy.phaseId(*face.outside()) );

          if( policy.ignoreInnerBoundary(tissue, tissueNeighbour) ){
            collection[vidx].isRelevant = false;
            continue;
          }
        }

        // store normal with corresponding phase index and mark the boundary vertex
        // on intersections the smaller phase index determines the orientation of the normal
        // i.e. the smaller phase index is always the associated to the inside phase
        if(tissue < tissueNeighbour-1.e-9){
          collection[vidx].normals.push_back(face.centerUnitOuterNormal());
          collection[vidx].weights.push_back(1.0/face.geometry().volume());
          collection[vidx].phaseIds.push_back(tissue);
          collection[vidx].ids.push_back(gridView.indexSet().index(entity));
        }
        // store the inverted normal of the "outside" phase
        if(tissue > tissueNeighbour+1.e-9){
          collection[vidx].normals.push_back(-1.*face.centerUnitOuterNormal());
          collection[vidx].weights.push_back(1.0/face.geometry().volume());
          collection[vidx].phaseIds.push_back(tissue);
          collection[vidx].ids.push_back(gridView.indexSet().index(entity));
        }
      }
    }

    /// store normals and associate them with edges (in 3D)
    template <class GridView, class Entity, class Face, class NormalCollection, class Scalar, class MergedPolicy>
    static void getEntitiesEdgeNormals(GridView const& gridView, Entity const& entity, Face const& face, NormalCollection& collection,
        Scalar tissue, Dune::GeometryType const& gt, MergedPolicy const& policy)
    {
      constexpr int dim = GridView::dimension;
      tissue = round(tissue);

      // index of face in cell
      int local_fidx = face.indexInInside();
      // iterate over edges
      for(int i=0; i<3; ++i){
        // get local edge id
        int local_eidx = Dune::GenericReferenceElements<Scalar,dim>::general(gt).subEntity(local_fidx,1,i,dim-1);
        // get global edge id
        int eidx = gridView.indexSet().subIndex(entity, local_eidx, dim-1);

        Scalar tissueNeighbour = -1;
        // get tissue of neighbouring cell sharing the intersection
        if(!face.boundary()){
          tissueNeighbour = round(policy.phaseId(*face.outside()));
          if(policy.ignoreInnerBoundary(tissue, tissueNeighbour)){
            collection[eidx].isRelevant = false;
            return;
          }
        }
        else
          tissueNeighbour = tissue + 1;

        // store normal with corresponding phase index and mark the boundary vertex
        // on intersections the smaller phase index determines the orientation of the normal
        // i.e. the smaller phase index is always the associated to the inside phase
        if(tissue < tissueNeighbour-1.e-9){
          collection[eidx].normals.push_back(face.centerUnitOuterNormal());
          collection[eidx].weights.push_back(1.0);
          collection[eidx].isRelevant = (!policy.ignoreOuterBoundary(face)) && collection[eidx].isRelevant;
          collection[eidx].phaseIds.push_back(tissue);
          collection[eidx].ids.push_back(gridView.indexSet().index(entity));
        }
        // store the inverted normal of the "outside" phase
        if(tissue > tissueNeighbour+1.e-9){
          collection[eidx].normals.push_back(-1.*face.centerUnitOuterNormal());
          collection[eidx].weights.push_back(1.0);
          collection[eidx].phaseIds.push_back(tissue);
          collection[eidx].ids.push_back(gridView.indexSet().index(entity));
        }
      }
    }

    // helper function 2D
    template <class GridView, class Entity, class Face, class NormalContainer, class Scalar, class Policy>
    inline void collectEntityNormals(GridView const& gridView, Entity const& entity, Face const& face, NormalContainer& container, Scalar tissue,
        Dune::GeometryType const gt, Policy const& policy, std::integral_constant<int,2>)
    {
      getEntitiesVertexNormals(gridView, entity, face, container.vertexNormals, tissue, gt, policy);
    }

    // helper function 3D
    template <class GridView, class Entity, class Face, class NormalContainer, class Scalar, class Policy>
    inline void collectEntityNormals(GridView const& gridView, Entity const& entity, Face const& face, NormalContainer& container, Scalar tissue,
        Dune::GeometryType const gt, Policy const& policy, std::integral_constant<int,3>)
    {
      getEntitiesVertexNormals(gridView, entity, face, container.vertexNormals, tissue, gt, policy);
      getEntitiesEdgeNormals(gridView, entity, face, container.edgeNormals, tissue, gt, policy);
    }


    // collect normals
    template <class GridView, class NormalContainer, class MergedPolicy>
    void getEntityNormals(GridView const& gridView, NormalContainer &container, MergedPolicy const& policy)
    {
      typedef typename GridView::template Codim<0>::Iterator CellIterator;
      typedef typename CellIterator::Entity Entity;
      typedef typename GridView::IntersectionIterator FaceIterator;
      typedef typename FaceIterator::Intersection Face;
      typedef typename NormalContainer::EntryType::Scalar Scalar;
      //iterate over cells
      std::for_each(gridView.template begin<0>(), gridView.template end<0>(),[&](Entity const& entity)
      {
        Scalar tissue = policy.phaseId(entity);
        Dune::GeometryType const gt(entity.type().id(), GridView::dimension);
        // iterate over faces
        FaceIterator fi = gridView.ibegin(entity),
                     fend = gridView.iend(entity);
        for(;fi!=fend;++fi) collectEntityNormals(gridView, entity, *fi, container, tissue, gt, policy, std::integral_constant<int,GridView::dimension>());
//        std::for_each(gridView.ibegin(entity), gridView.iend(entity),[&](Face const& face)
//        {
//          // emulate template specification via overloading and std::integral_constant
//          collectEntityNormals(gridView, entity, face, container, tissue, gt, policy, std::integral_constant<int,GridView::dimension>());
//        });
      });
    }

    template <class NormalContainer>
    inline void computeMeanNormalsImpl(NormalContainer& container, std::integral_constant<int,2>)
    {
      for(auto& collection : container.vertexNormals) InterpolationTools::computeMeanNormal<true>(collection);
    }

    template <class NormalContainer>
    inline void computeMeanNormalsImpl(NormalContainer& container, std::integral_constant<int,3>)
    {
      for(auto& collection : container.vertexNormals) InterpolationTools::computeMeanNormal<true>(collection);
      for(auto& collection : container.edgeNormals) InterpolationTools::computeMeanNormal<false>(collection);
    }
    // mean normals
    template <class NormalContainer>
    inline void computeMeanNormals(NormalContainer& container)
    {
      computeMeanNormalsImpl(container, std::integral_constant<int,NormalContainer::EntryType::dim>());
    }
  }

  /// Collects face normals associated to corners (2D,3D) and edges(3D only).
  /**
   * Users can specify the collectors behaviour via 3 policies (see namespace Policy):
   * a) PhasePolicy: specifies phases
   * b) OuterBoundaryPolicy: specifies the behaviour on the domain boundary
   * c) InnerBoundaryPolicy: specifies the behaviour on inner boundaries (intersections of Codim<0>-entities where different phases meet). A specification
   *    of this policy only makes sense if different phases have been specified (i.e. in a (discontinuous) FunctionSpaceElement).
   */
  template<class GridView, class OuterBoundaryPolicy = Policy::ConsiderOuterBoundary, class InnerBoundaryPolicy = Policy::IgnoreInnerBoundary,
           class Phase = typename GridView::ctype, template <class> class PhasePolicy = Policy::NoPhaseInfo>
  class BoundaryNormalCollector
  {
  public:
    static int const dim = GridView::dimension;
    typedef typename GridView::ctype Scalar;
    typedef InterpolationTools::NormalContainer<Scalar,dim> NormalContainer;

    /// Constructor for the case that no phase ids are provided
    BoundaryNormalCollector(GridView const& gridView_, OuterBoundaryPolicy const& outerBoundaryPolicy = OuterBoundaryPolicy(),
        InnerBoundaryPolicy const& innerBoundaryPolicy = InnerBoundaryPolicy())
    : gridView(gridView_), mergedPolicy(outerBoundaryPolicy,innerBoundaryPolicy, Policy::NoPhaseInfo<Scalar>())
    {}

    /// Constructor for the case that phase ids are provided.
    /**
     * Providing an phase-element together with Policy::NoPhaseInfo will result in a compile-time error.
     */
    BoundaryNormalCollector(GridView const& gridView_, Phase const& phase, OuterBoundaryPolicy const& outerBoundaryPolicy = OuterBoundaryPolicy(),
        InnerBoundaryPolicy const& innerBoundaryPolicy = InnerBoundaryPolicy())
    : gridView(gridView_), mergedPolicy(outerBoundaryPolicy,innerBoundaryPolicy,PhasePolicy<Phase>(phase))
    {}

    /// Get normals of all faces incident to each vertex.
    NormalContainer normals() const
    {
      NormalContainer container(gridView);
      // identify codim-2- or codim-3-entities lying on the boundary without
      // being part of a codim-1-entity on the boundary
      BoundaryNormalCollectorDetail::outerBoundaryCheck(gridView, container);
      BoundaryNormalCollectorDetail::getEntityNormals(gridView, container, mergedPolicy);
      return container;
    }

    /// Get average of normals incident to each vertex.
    NormalContainer meanNormals() const
    {
      NormalContainer container(gridView);
      // identify codim-2- or codim-3-entities lying on the boundary without
      // being part of a codim-1-entity on the boundary
      BoundaryNormalCollectorDetail::outerBoundaryCheck(gridView, container);
      BoundaryNormalCollectorDetail::getEntityNormals(gridView, container, mergedPolicy);
      BoundaryNormalCollectorDetail::computeMeanNormals(container);
      return container;
    }

    // no default and copy constructor, no copy assignment
    BoundaryNormalCollector() = delete;
    BoundaryNormalCollector(BoundaryNormalCollector const&) = delete;
    BoundaryNormalCollector& operator=(BoundaryNormalCollector const&) = delete;

    ~BoundaryNormalCollector() = default;

  private:
    GridView gridView;
    Policy::MergedPolicy<OuterBoundaryPolicy, InnerBoundaryPolicy, PhasePolicy<Phase> > mergedPolicy;
  };
}
#endif /* BOUNDARY_NORMAL_COLLECTOR_HH */
