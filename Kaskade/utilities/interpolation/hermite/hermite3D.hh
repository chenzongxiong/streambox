/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2002-2011 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef HERMITE_INTERPOLATION_3D_HH
#define HERMITE_INTERPOLATION_3D_HH

#include "utilities/geometry/geomtools.hh"
#include "utilities/interpolation/hermite/base.hh"

/**
 * @file
 * @brief  3D hermite interpolation.
 * @author Lars Lubkoll
 */
/**
 * \cond internals
 */
namespace Kaskade
{
  /// Implementation for 3D
  /*
   * Hermite interpolation on a simplicial grid. It is assumed that the function values
   * of the vertices are zero. Moreover the tangent planes at each vertex and each edge
   * are given through its outer normals. Note that the deformation can only act in
   * normal direction. On the faces this is clear. On the edges the normal direction is
   * determined by the normals of the two adjacent faces.
   * The interpolation is split into two independant parts:
   * - First use cubic hermite interpolation at the edges in order to achieve tangent
   *   plane continuity at the vertices. This approach also minimizes the second derivative
   *   of the curves parametrization over [0,1] in L2.
   * - Second adjust the remaining shape functions. Note that with the chosen "3rd" order
   *   hierarchic shape functions it is not possible to guarantee tangent plane continuity
   *   over the edges. Thus tangent plane continuity is only required at the edge's center.
   *
   * \param RT floating point type
   * \param GridView GridView
   * \param ShapeFunctionSet type of the shape function set of the codim 0 reference element
   * \param InterpolationPolicy class specifying the interpolation details
   * \param Container container class storing vertex and edge normals
   */
  template<class Scalar_, class GridView, class InterpolationPolicy>
  class HIPImpl<Scalar_, GridView, InterpolationPolicy, 3> : public HIPBase<Scalar_,3>
  {
  public:
    typedef HIPBase<Scalar_,3> Base;
    static int const dim = 3;
    typedef typename Base::Scalar Scalar;
    typedef typename Base::Vector Vector;
    typedef typename GridView::template Codim<0>::Iterator CellIterator;
    typedef typename CellIterator::Entity Entity;
    typedef InterpolationTools::NormalContainer<Scalar,dim> Container;
    typedef typename Container::CollectionType NormalCollection;
    typedef typename GridView::IntersectionIterator LocalFaceIterator;
    typedef typename LocalFaceIterator::Intersection Face;
    typedef typename GridView::Grid Grid;
    typedef Dune::GenericReferenceElements<Scalar,dim> ReferenceTetrahedron;
    typedef Dune::GenericReferenceElements<Scalar,dim-1> ReferenceTriangle;


    /// Constructor for the use with interpolation polynomial collection.
    /**
     * Constructs the representation of the zero function.
     *
     * \param numberOfShapeFunctions size of shape function set
     */
    HIPImpl(int numberOfShapeFunctions) : Base(numberOfShapeFunctions), origin(0) {}

    /// Constructor
    /*
     * \param entity entity on which the polynomial is defined
     * \param gridView grid view
     * \param container container holding vertex normals and edge normals in container.vertexNormals/container.edgeNormals
     * \param shapefunctionset shape function set of the reference codim 0 element
     * \param policy policy object specifiying interpolation details
     */
    template <class ShapeFunctionSet>
    HIPImpl(Entity const& entity, GridView const& gridView,
        ShapeFunctionSet const& shapeFunctionSet, Container const& container,
        InterpolationPolicy const& policy) :
        Base(shapeFunctionSet.size()), origin(0)
        {
      init(entity, gridView, shapeFunctionSet, container, policy);
        }

    /// Initialize interpolation polynomial.
    /**
     * \param entity entity on which the polynomial is defined
     * \param gridView grid view
     * \param container container holding vertex normals and edge normals in container.vertexNormals/container.edgeNormals
     * \param shapefunctionset shape function set of the reference codim 0 element
     * \param policy policy object specifiying interpolation details
     */
    template <class ShapeFunctionSet>
    void init(Entity const& entity, GridView const& gridView, ShapeFunctionSet const& shapeFunctionSet, Container const& container, InterpolationPolicy const& policy)
    {
      initSamplingPointValuesOnEdges(entity, gridView, shapeFunctionSet, container, policy);
      initSamplingPointValuesOnFaces(entity, gridView, shapeFunctionSet, container, policy);
    }

  private:
    /// Interpolation on edges.
    template <class ShapeFunctionSet>
    void initSamplingPointValuesOnEdges(Entity const& entity, GridView const  &gridView, ShapeFunctionSet const& shapeFunctionSet, Container const& container, InterpolationPolicy const& policy){
      // get global vertex ids
      std::vector<int> globalVertexIds = InterpolationTools::getGlobalVertexIds(entity, gridView);
      // get tissue id of the current cell
      Scalar tissue = round(policy.phaseId( entity ));

      // get "normal" vectors associated with the vertices
      std::vector<Vector> entityVertexNormals(dim+1, Vector());
      InterpolationTools::getVertexNormals(globalVertexIds, container.vertexNormals, entityVertexNormals, tissue);

      LocalFaceIterator fend = gridView.iend(entity);
      for(LocalFaceIterator face = gridView.ibegin(entity); face!=fend; ++face){

        Scalar tissueNeighbour = -1;
        // get neighbouring tissue id if neighbour exists and material information is provided
        if( face->neighbor() )
          tissueNeighbour = round(policy.phaseId( *face->outside() ));

        // criteria for skipping calculation
        // if no inner boundary exists (no intersection with cells associated with other materials) continue
        if(fabs(tissue-tissueNeighbour) < 1e-9) continue;

        if( policy.ignoreOuterBoundary(*face) ) continue;

        // id of face in cell
        int faceId = face->indexInInside();

        // iterate over edges
        for(int edgeIdInFace = 0; edgeIdInFace < dim; ++edgeIdInFace){

          // get corner and edge ids
          int edgeIdInEntity = ReferenceTetrahedron::simplex().subEntity(faceId,1,edgeIdInFace,dim-1);
          std::vector<int> edgeVertexIdsInEntity = InterpolationTools::getEdgeCornerIdsInEntity(entity, edgeIdInEntity);

          // local id of the edge in the cell
          // global edge id
          int globalEdgeId = gridView.indexSet().subIndex(entity, edgeIdInEntity, dim-1);

          // if edge has been marked as not relevant skip it
          if( !container.edgeNormals[globalEdgeId].isRelevant ) continue;

          // create normalized vector in edge direction
          Vector edgeDirection = entity.geometry().corner(edgeVertexIdsInEntity[1]) - entity.geometry().corner(edgeVertexIdsInEntity[0]);
          GeomTools::normalize(edgeDirection);

          // get normal vector of the plane spanned the edge direction and the edge "normal"
          Vector planeNormal;

          Vector edgeNormal = getEdgeNormal(globalEdgeId, tissue, container.edgeNormals);

          GeomTools::crossProduct(edgeNormal, edgeDirection, planeNormal);

          std::vector<Vector> projectedVertexNormals(2), desiredGradientDirections(2,Vector(1e-301));
          Dune::FieldVector<Scalar,2> desiredGradients;
          std::vector<int> shapeFunctionIds(2);
          for(int vertexIdInEdge=0; vertexIdInEdge<2; ++vertexIdInEdge){
            // ids of 2nd and 3rd order shape functions associated with the considered edge
            shapeFunctionIds[vertexIdInEdge] = 4 + edgeIdInEntity + vertexIdInEdge*6;

            // skip calculation of rhs values
            if( !container.vertexNormals[ globalVertexIds[ edgeVertexIdsInEntity[vertexIdInEdge] ] ].isRelevant ){
              desiredGradients[vertexIdInEdge] = 0;
              continue;
            }

            // get convenient coordinate system
            projectedVertexNormals[vertexIdInEdge] = entityVertexNormals[edgeVertexIdsInEntity[vertexIdInEdge]];
            GeomTools::projectOnPlane(projectedVertexNormals[vertexIdInEdge], planeNormal);
            GeomTools::normalize(projectedVertexNormals[vertexIdInEdge]);
            GeomTools::crossProduct(planeNormal, projectedVertexNormals[vertexIdInEdge], desiredGradientDirections[vertexIdInEdge]);

            // first determine the desired gradients at the (deformed) edges and apply threshold
            // set the gradient = 0 if on the corresponding edge no smoothing is desired
            if(container.vertexNormals[globalVertexIds[ edgeVertexIdsInEntity[vertexIdInEdge] ] ].isRelevant){
              desiredGradients[vertexIdInEdge] = (edgeNormal * desiredGradientDirections[vertexIdInEdge])
                  /(edgeDirection * desiredGradientDirections[vertexIdInEdge]);
              policy.applyGradientThreshold(desiredGradients[vertexIdInEdge]);
            }
            else
            {
              desiredGradients[vertexIdInEdge] = 0;
            }
          }

          // evaluate derivatives of shape functions at the edge's corners
          Dune::FieldMatrix<Scalar,2,2> A;
          for(int vertexIdInEdge=0; vertexIdInEdge<2; ++vertexIdInEdge){
            for(int sfId = 0; sfId<2; ++sfId){
              Dune::FieldMatrix<Scalar,dim,dim> jacobian = entity.geometry().jacobianInverseTransposed(entity.geometry().local(entity.geometry().corner(edgeVertexIdsInEntity[vertexIdInEdge])));
              Vector sf_grad = shapeFunctionSet[shapeFunctionIds[sfId]].evaluateDerivative( entity.geometry().local( entity.geometry().corner(edgeVertexIdsInEntity[vertexIdInEdge]) ) )[0];
              Vector sf_grad_chain;
              jacobian.mv(sf_grad, sf_grad_chain);
              A[vertexIdInEdge][sfId] = sf_grad_chain * edgeDirection;
            }
          }

          // determine coefficients of the hierarchic ansatz functions associated with this edge
          Dune::FieldVector<Scalar,2> coefficients;
          A.solve(coefficients, desiredGradients);

          // store calculated coefficients
          for(int sfId=0; sfId<2; ++sfId)
            this->insertEntry(coefficients[sfId], edgeNormal, shapeFunctionIds[sfId]);
        }
      }
    }


    /// Interpolation on faces.
    template <class ShapeFunctionSet>
    void initSamplingPointValuesOnFaces(Entity const& entity,
        GridView const  &gridView,
        ShapeFunctionSet const& shapeFunctionSet,
        Container const& container,
        InterpolationPolicy const& policy)
    {
      // get tissue id of the current cell
      Scalar tissue = round(policy.phaseId(entity));

      // iterate over faces
      LocalFaceIterator fend = gridView.iend(entity);
      for(LocalFaceIterator face = gridView.ibegin(entity); face!=fend; ++face){

        // id of face in cell
        int faceId = face->indexInInside();

        Scalar tissueNeighbour = -1;
        // get neighbouring tissue id if neighbour exists and material information is provided
        if( face->neighbor() )
          tissueNeighbour = round(policy.phaseId( *face->outside() ));

        // if no inner boundary exists (no intersection with cells associated with other materials) continue
        if( fabs(tissue-tissueNeighbour) < 1e-9 ) continue;
        // continue if policy decides to skip face
        if( policy.ignoreOuterBoundary(*face) ) continue;

        // interpolation points in local coordinates
        std::vector<Vector> iNodes(dim);

        // affine coordinate mapping -> jacobian is constant on each entity
        Dune::FieldMatrix<Scalar,dim,dim> jacobian = entity.geometry().jacobianInverseTransposed(origin);

        // initialize variables
        Vector desiredGradients;	// right hand side vector
        std::vector<int> edgeIdInEntity(dim,0), globalEdgeIds(dim,0), shapeFunctionIds(6,0); // store used indices
        std::vector<std::vector<int> > edgeVertexIdsInEntity(dim);	// more indices
        std::vector<Vector> edgeNormal(dim,origin),
            planeNormal(dim,origin), // note that plane denotes the plane spanned by the edge and the "real" surface normal
            projectedPlaneNormal(dim,origin),
            deformedProjectedPlaneNormal(dim,origin),
            deformedFaceNormal(dim,origin); // store used vectors

        Vector faceNormal = face->centerUnitOuterNormal();
        // adjust orientation such that a unique normal direction is associated with every face
        if(tissueNeighbour+1e-9 < tissue && tissueNeighbour > 1e-9)
          faceNormal *= -1;

        // get ids
        for(int edgeIdInFace = 0; edgeIdInFace<dim; ++edgeIdInFace){
          edgeIdInEntity[edgeIdInFace] = ReferenceTetrahedron::simplex().subEntity(faceId,1,edgeIdInFace,dim-1);
          globalEdgeIds[edgeIdInFace] = gridView.indexSet().subIndex(entity, edgeIdInEntity[edgeIdInFace], dim-1);
          edgeVertexIdsInEntity[edgeIdInFace] = getExtendedEdgeVertexIds(edgeIdInFace, faceId);
          shapeFunctionIds[edgeIdInFace] = 4 + edgeIdInEntity[edgeIdInFace];
          shapeFunctionIds[edgeIdInFace+dim] = 10 + edgeIdInEntity[edgeIdInFace];
        }

        // iterate over edges and collect directions and interpolation nodes
        for(int edgeIdInFace = 0; edgeIdInFace<dim; ++edgeIdInFace){

          // get interpolation points
          iNodes[edgeIdInFace] = entity.geometry().local( 0.5 * ( entity.geometry().corner(edgeVertexIdsInEntity[edgeIdInFace][0]) + entity.geometry().corner(edgeVertexIdsInEntity[edgeIdInFace][1]) ) );

          // get normalized edge vector
          Vector edgeDirection = entity.geometry().corner(edgeVertexIdsInEntity[edgeIdInFace][0]) - entity.geometry().corner(edgeVertexIdsInEntity[edgeIdInFace][1]);
          GeomTools::normalize(edgeDirection);

          // check if edge is relevant and if this is the case get associated edge normal
          // and compute the associated tangential directions in the necessary coordinate systems
          if(container.edgeNormals[ globalEdgeIds[edgeIdInFace ] ].isRelevant)
          {
            // get edge normal
            edgeNormal[edgeIdInFace] = getEdgeNormal(globalEdgeIds[edgeIdInFace], tissue, container.edgeNormals);
            // get tangential direction
            GeomTools::crossProduct(edgeNormal[edgeIdInFace], edgeDirection, planeNormal[edgeIdInFace]);
            // and project it on triangle
            projectedPlaneNormal[edgeIdInFace] = planeNormal[edgeIdInFace];
            GeomTools::projectOnPlane(projectedPlaneNormal[edgeIdInFace], faceNormal);
            GeomTools::normalize(projectedPlaneNormal[edgeIdInFace]);
            // correct orientation such that projectedPlaneNormal points into the face
            int remainingVertexId = edgeVertexIdsInEntity[edgeIdInFace][2];
            Vector anotherEdgeDirection = entity.geometry().corner(remainingVertexId) - entity.geometry().corner(edgeVertexIdsInEntity[edgeIdInFace][0]);
            if(anotherEdgeDirection * projectedPlaneNormal[edgeIdInFace] < 0)
              projectedPlaneNormal[edgeIdInFace] *= -1;

            // project projected tangential direction on deformed surface and normalize
            Dune::FieldMatrix<Scalar,dim,dim> deformationGradient = getPartialDeformationGradient(iNodes[edgeIdInFace], entity, shapeFunctionIds, shapeFunctionSet);
            // use derivative, not gradient for the transformation of directions
            deformationGradient.mtv(projectedPlaneNormal[edgeIdInFace], deformedProjectedPlaneNormal[edgeIdInFace]);
            GeomTools::normalize(deformedProjectedPlaneNormal[edgeIdInFace]);
            deformationGradient.mtv(faceNormal, deformedFaceNormal[edgeIdInFace]);
            // ??? GeomTools::normalize(deformedFaceNormal[edgeIdInFace]);
          }
          else
          {
            int remainingVertexId = edgeVertexIdsInEntity[edgeIdInFace][2];
            deformedProjectedPlaneNormal[edgeIdInFace] = entity.geometry().corner(remainingVertexId);
            deformedProjectedPlaneNormal[edgeIdInFace]-= entity.geometry().global(iNodes[edgeIdInFace]);
            GeomTools::normalize(deformedProjectedPlaneNormal[edgeIdInFace]);
          }
          // else deformedProjectedPlaneNormal[edgeIdInFace] = origin (initial value)
        }

        /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
        /* fill Matrix */
        /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
        // interpolate coefficients for the hierarchic shape functions associated with the face
        // shape functions associated with the vertices vanish -> consider 2nd and 3rd order terms only.
        Dune::FieldMatrix<Scalar,dim,dim> A(0);

        for(int edgeIdInFace=0; edgeIdInFace<dim; ++edgeIdInFace)
          for(int sfID=0; sfID<dim; ++sfID){
            // get gradients on reference tetrahedron
            Vector grad = shapeFunctionSet[16 + faceId*dim + sfID].evaluateDerivative(iNodes[edgeIdInFace])[0];

            Dune::FieldMatrix<Scalar,dim,dim> inverseDeformationGradient = getPartialInverseDeformationGradient(iNodes[edgeIdInFace], entity, shapeFunctionIds, shapeFunctionSet);
            Vector chainRuleTmp(0);
            // apply chain rule -> gradient on the actual tetrahedron
            jacobian.mv(grad, chainRuleTmp);
            // apply chain rule -> gradient on the deformed tetrahedron
            inverseDeformationGradient.mv(chainRuleTmp,grad);

            // evaluate directional derivative in deformed surface
            A[edgeIdInFace][sfID] = grad * deformedProjectedPlaneNormal[edgeIdInFace];
          }


        /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
        /* fill right hand side */
        /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
        // first determine the desired gradients at the (deformed) edges and apply threshold
        // set the gradient = 0 if on the corresponding edge no smoothing is desired
        for(int iNode = 0; iNode<dim; ++iNode){
          // desired normal
          if(container.edgeNormals[ globalEdgeIds[iNode] ].isRelevant){
            desiredGradients[iNode] = (planeNormal[iNode] * deformedFaceNormal[iNode]) / (planeNormal[iNode] * deformedProjectedPlaneNormal[iNode]);
            policy.applyGradientThreshold(desiredGradients[iNode]);
          }
          else
            desiredGradients[iNode] = 0;
        }

        /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
        /* then subtract contributions from other shape functions */
        /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
        for(int iNode = 0; iNode<dim; ++iNode){
          // subtract contributions from other shape functions
          for(int edgeIdInFace = 0; edgeIdInFace < dim; ++edgeIdInFace){
            Dune::FieldMatrix<Scalar,dim,dim> inverseDeformationGradient = getPartialInverseDeformationGradient(iNodes[iNode], entity, shapeFunctionIds, shapeFunctionSet);
            Vector chainRuleTmp(0);

            // second order shape function
            Vector grad = shapeFunctionSet[ shapeFunctionIds[edgeIdInFace] ].evaluateDerivative(iNodes[iNode])[0];

            jacobian.mv(grad,chainRuleTmp);
            inverseDeformationGradient.mv(chainRuleTmp,grad);

            desiredGradients[iNode] -= grad * deformedProjectedPlaneNormal[iNode] * this->values[ shapeFunctionIds[edgeIdInFace] ];

            // third order shape function
            grad = shapeFunctionSet[ shapeFunctionIds[edgeIdInFace+dim]].evaluateDerivative(iNodes[iNode])[0];

            jacobian.mv(grad,chainRuleTmp);
            inverseDeformationGradient.mv(chainRuleTmp,grad);

            desiredGradients[iNode] -= grad * deformedProjectedPlaneNormal[iNode] * this->values[ shapeFunctionIds[edgeIdInFace+dim] ];
          }
        }

        // solve LSE
        Vector coefficients;
        A.solve(coefficients, desiredGradients);
        if(!(coefficients==coefficients))
        {
          for(int edgeIdInFace=0; edgeIdInFace<dim; ++edgeIdInFace)
            for(int sfID=0; sfID<dim; ++sfID){
              // get gradients on reference tetrahedron
              Vector grad = shapeFunctionSet[16 + faceId*dim + sfID].evaluateDerivative(iNodes[edgeIdInFace])[0];

              Dune::FieldMatrix<Scalar,dim,dim> inverseDeformationGradient = getPartialInverseDeformationGradient(iNodes[edgeIdInFace], entity, shapeFunctionIds, shapeFunctionSet);
              Vector chainRuleTmp(0);
              // apply chain rule -> gradient on the actual tetrahedron
              jacobian.mv(grad, chainRuleTmp);


              // apply chain rule -> gradient on the deformed tetrahedron
              inverseDeformationGradient.mv(chainRuleTmp,grad);

              std::cout << grad << " * " << deformedProjectedPlaneNormal[edgeIdInFace] << std::endl << std::endl;

            }
          std::cout << A;
          std::cout << desiredGradients << " -> " << coefficients << std::endl << std::endl;
          coefficients = 0.0;

        }


        for(int sfId=0; sfId<dim; ++sfId)
          this->insertEntry(coefficients[sfId], faceNormal, 16 + dim*faceId + sfId);
      }
    }

    /// Get a vector with 3 ids associated with a faces vertices.
    /* The first two ids are the edge ids wrt. the entity. The last is the remaining vertex id wrt. the entity.
     *
     * \param edgeIdInFace id of the edge wrt. the face.
     * \return vector containing the local vertex ids of the edges vertices and additionally in the last entry the id of the remaining vertex of the face.
     */
    std::vector<int> getExtendedEdgeVertexIds(int edgeIdInFace, int faceId) const {
      std::vector<int> result(dim);
      result[0] = ReferenceTriangle::simplex().subEntity(edgeIdInFace,1,0,dim-1);
      result[1] = ReferenceTriangle::simplex().subEntity(edgeIdInFace,1,1,dim-1);
      result[2] = dim - result[1] - result[0];
      result[0] = ReferenceTetrahedron::simplex().subEntity(faceId,1, result[0], dim);
      result[1] = ReferenceTetrahedron::simplex().subEntity(faceId,1, result[1], dim);
      result[2] = ReferenceTetrahedron::simplex().subEntity(faceId,1, result[2], dim);
      return result;
    }

    /// Get edge normal corresponding to the material with id tissueId.
    /*
     * \param eid global id of the edge
     * \param tissueId material id
     * \result edge normal of edge eid associated with tissueId
     */
    Vector getEdgeNormal(int eid, Scalar tissueId, NormalCollection const& edgeNormals){
      for(int i=0; i<edgeNormals[eid].phaseIds.size(); ++i)
        if(fabs(edgeNormals[eid].phaseIds[i] - round(tissueId)) < 1e-9)
          return edgeNormals[eid].normals[i];
      std::cout << "Warning: No edge normal found! " << std::endl;
      return origin;
    }

    /// Get deformation gradient. The relevant shape function ids must be specified.
    /*
     * \param iNode evaluation point
     * \param entity entity containing the evaluation point
     * \param shapeFunctionIds ids of the used shape functions
     * \param shapeFunctionSet shape function set
     *
     * \return deformation gradient
     */
    template <class ShapeFunctionSet>
    Dune::FieldMatrix<Scalar,dim,dim> getPartialDeformationGradient(Vector const& iNode, Entity const& entity, std::vector<int> const& shapeFunctionIds, ShapeFunctionSet const& shapeFunctionSet) const {
      Dune::FieldMatrix<Scalar,dim,dim> deformationGradient(0);
      Dune::FieldMatrix<Scalar,dim,dim> jacobian = entity.geometry().jacobianInverseTransposed(Vector(0));
      for(int i=0;i<dim; ++i){
        for(int j=0; j<dim; ++j){
          for(int k=0; k<shapeFunctionIds.size(); ++k){
            // evaluation on reference element
            Vector grad = shapeFunctionSet[ shapeFunctionIds[k] ].evaluateDerivative( iNode )[0];
            Vector chain(0);
            // chain rule -> derivative on actual element
            jacobian.mv(grad,chain);
            deformationGradient[j][i] += this->values[ shapeFunctionIds[k] ] * chain[j] * this->directions[ shapeFunctionIds[k] ][i];
          }
        }
        // add identity
        deformationGradient[i][i] += 1;
      }
      return deformationGradient;
    }

    /// Get inverse of deformation gradient. The relevant shape function ids must be specified.
    /*
     * \param iNode evaluation point
     * \param entity entity containing the evaluation point
     * \param shapeFunctionIds ids of the used shape functions
     * \param shapeFunctionSet shape function set
     *
     * \return inverse of deformation gradient
     */
    template <class ShapeFunctionSet>
    Dune::FieldMatrix<Scalar,dim,dim> getPartialInverseDeformationGradient(Vector const& iNode, Entity const& entity, std::vector<int> const& shapeFunctionIds, ShapeFunctionSet const& shapeFunctionSet) const {
      Dune::FieldMatrix<Scalar,dim,dim> inverseDeformationGradient = getPartialDeformationGradient(iNode, entity, shapeFunctionIds, shapeFunctionSet);
      inverseDeformationGradient.invert();
      return inverseDeformationGradient;
    }

    Vector const origin;
  };
} // namespace Kaskade
/**
 * \endcond
 */
#endif // HERMITE_INTERPOLATION_3D_HH
