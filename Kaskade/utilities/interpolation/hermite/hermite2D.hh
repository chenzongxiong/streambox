/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/Numerik/numsoft/kaskade7/                        */
/*                                                                           */
/*  Copyright (C) 2002-2013 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**
 * \file
 * \author Lars Lubkoll
 */

#ifndef HERMITE_INTERPOLATION_2D
#define HERMITE_INTERPOLATION_2D

#include "utilities/geometry/geomtools.hh"
#include "utilities/interpolation/hermite/base.hh"

/**
 * \cond internals
 */
namespace Kaskade
{
  /// Implementation 2D
  template<class Scalar_, class GridView, class InterpolationPolicy>
  class HIPImpl<Scalar_, GridView, InterpolationPolicy, 2> : public HIPBase<Scalar_,2>
  {
  public:
    typedef HIPBase<Scalar_,2> Base;
    static int const dim = 2;
    typedef typename Base::Scalar Scalar;
    typedef typename Base::Vector Vector;
    typedef typename GridView::template Codim<dim>::Iterator VertexIterator;
    typedef typename GridView::template Codim<0>::Iterator::Entity Entity;
    typedef InterpolationTools::NormalContainer<Scalar,dim> Container;
    typedef typename Container::CollectionType NormalCollection;
    typedef typename GridView::IntersectionIterator LocalFaceIterator;
    typedef typename LocalFaceIterator::Intersection Face;

    /// Constructor for the use with interpolation polynomial collection.
    /**
     * Constructs the representation of the zero function.
     *
     * \param numberOfShapeFunctions size of shape function set
     */
    HIPImpl(int numberOfShapeFunctions) : Base(numberOfShapeFunctions){}

    /// Constructor
    /*
     * \param entity entity on which the polynomial is defined
     * \param gridView grid view
     * \param container container holding vertex normals in container.vertexNormals
     * \param shapefunctionset shape function set of the reference codim 0 element
     * \param policy policy object specifiying interpolation details
     */
    template <class ShapeFunctionSet>
    HIPImpl(Entity const& entity, GridView const& gridView,
        ShapeFunctionSet const& shapeFunctionSet, Container const& container,
        InterpolationPolicy const& policy) :
        Base(shapeFunctionSet.size())
    {
      init(entity, gridView, shapeFunctionSet, container, policy);
    }

    /// Initialize interpolation polynomial.
    /**
     * \param entity entity on which the polynomial is defined
     * \param gridView grid view
     * \param container container holding vertex normals in container.vertexNormals
     * \param shapefunctionset shape function set of the reference codim 0 element
     * \param policy policy object specifiying interpolation details
     * \param tol tolerance for the decision if two floating point numbers coincide   */
    template <class ShapeFunctionSet>
    void init(Entity const& entity, GridView const& gridView, ShapeFunctionSet const& shapeFunctionSet, Container const& container, InterpolationPolicy const& policy)
    {
      // get global vertex ids of entity
      std::vector<int> vid = InterpolationTools::getGlobalVertexIds(entity,gridView);
      // get directions of tangents associated with the entities vertices
      std::vector<Vector> tangents(dim+1, Vector(0));
      getVertexTangents(vid, container.vertexNormals, tangents);

      // get material id
      Scalar tissue = round(policy.phaseId(entity));
      // affine coordinate mapping -> jacobian is constant on each entity
      Dune::FieldMatrix<Scalar,dim,dim> jacobian = entity.geometry().jacobianInverseTransposed(Vector(0));

      // iterate over edges
      LocalFaceIterator fend = gridView.iend(entity);
      for(LocalFaceIterator face = gridView.ibegin(entity); face!=fend; ++face){
        // use policy to exclude interpolation on specific faces
        if( face->boundary() ) if( policy.ignoreOuterBoundary(*face) ) continue;

        Scalar tissueNeighbour = -1;
        if( face->neighbor() ) tissueNeighbour = policy.phaseId(*face->outside());

        // if no inner boundary exists (no intersection with cells associated with other materials) continue
        if(fabs(tissue-tissueNeighbour) < 1e-9) continue;

        // id of face in entity
        int edgeId = face->indexInInside();
        // get faces vertex ids wrt to entity
        std::vector<int> localEdgeVertexIds = InterpolationTools::getEdgeCornerIdsInEntity(entity, edgeId),
            globalEdgeVertexIds = InterpolationTools::getGlobalEdgeCornerIds(entity, edgeId, gridView),
            shapeFunctionIds = getFacesShapeFunctionIds(edgeId);

        VertexIterator iter = gridView.template begin<dim>(), iter2 = gridView.template begin<dim>();

        for(int i=0; i<globalEdgeVertexIds[0]; ++i) ++iter;
        for(int i=0; i<globalEdgeVertexIds[1]; ++i) ++iter2;
        // get outer normal of face
        Vector faceNormal = face->centerUnitOuterNormal();

        // get face direction
        Vector edgeDirection = face->geometry().corner(1) - face->geometry().corner(0);
        GeomTools::normalize(edgeDirection);
        /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
        /* fill right hand side */
        /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
        // first determine the desired gradients and apply threshold
        Vector desiredGradients;
        for(int vertexIdInEdge = 0; vertexIdInEdge<dim; ++vertexIdInEdge){
          if(container.vertexNormals[ globalEdgeVertexIds[vertexIdInEdge] ].isRelevant)
          {
            desiredGradients[vertexIdInEdge] = (tangents[localEdgeVertexIds[vertexIdInEdge]] * faceNormal) / (tangents[localEdgeVertexIds[vertexIdInEdge]] * edgeDirection);
            policy.applyGradientThreshold(desiredGradients[vertexIdInEdge]);
          }
          else
            desiredGradients[vertexIdInEdge] = 0;
        }

        /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
        /* fill Matrix */
        /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
        // interpolate coefficients for the lagrange shape functions associated with the face
        Dune::FieldMatrix<Scalar,dim,dim> A(0);

        for(int vertexIdInFace=0; vertexIdInFace<dim; ++vertexIdInFace){
          //      	std::cout << "corners: " << face->geometry().corner(vertexIdInFace) << std::endl;
          for(int sfId=0; sfId<dim; ++sfId){
            // get gradients on reference tetrahedron
            Vector grad = ( (dim==2 && vertexIdInFace==0) ? shapeFunctionSet[ shapeFunctionIds[sfId] ].evaluateDerivative( entity.geometry().local(iter->geometry().corner(0)) )[0]
                                                                                                                                                                                 : shapeFunctionSet[ shapeFunctionIds[sfId] ].evaluateDerivative( entity.geometry().local(iter2->geometry().corner(0)) )[0] );
            Vector chainRule(0);
            // apply chain rule -> gradient on the actual tetrahedron
            jacobian.mv(grad, chainRule);
            // evaluate directional derivative in deformed surface
            A[vertexIdInFace][sfId] = chainRule * edgeDirection;
          }
        }

        Vector coefficients(0);
        A.solve(coefficients, desiredGradients);

        for(int sfId=0; sfId<dim; ++sfId){
          Scalar orientation = 1;
          if(tissueNeighbour+1e-9 < tissue && tissueNeighbour > -1e-9) orientation = -1;
          this->insertEntry(orientation*coefficients[sfId], orientation*faceNormal, shapeFunctionIds[sfId]);
        }
      }
    }

  private:
    /// Get the vertices' "tangents"
    /**
     * Each entry in tangents will be an orthogonal (wrt euclidean scalar product) vector to the vertices
     * normal direction.
     *
     * \param vid vertex indices
     * \param vertexNormals collection of directions associated with vertices and materials
     * \param tangents storage for the output vectors
     */
    void getVertexTangents(std::vector<int> const &vid, NormalCollection const& vertexNormals, std::vector<Vector> &tangents) const {
      std::vector<Vector> normals(3,Vector());
      normals[0] = vertexNormals[vid[0]].normals[0];
      normals[1] = vertexNormals[vid[1]].normals[0];
      normals[2] = vertexNormals[vid[2]].normals[0];
      tangents[0][0] = normals[0][1];
      tangents[0][1] = -normals[0][0];
      tangents[1][0] = normals[1][1];
      tangents[1][1] = -normals[1][0];
      tangents[2][0] = normals[2][1];
      tangents[2][1] = -normals[2][0];
    }

    /// Get ids of shape function associated with face.
    /**
     * \param edgeId index of the face
     * \return shape functions ids associated with the face
     */
    std::vector<int> getFacesShapeFunctionIds(int edgeId){
      std::vector<int> shapeFunctionIds(dim,-1);
      switch(edgeId){
      case 0:
        shapeFunctionIds[0] = 4;
        shapeFunctionIds[1] = 7;
        break;
      case 1:
        shapeFunctionIds[0] = 2;
        shapeFunctionIds[1] = 1;
        break;
      case 2:
        shapeFunctionIds[0] = 8;
        shapeFunctionIds[1] = 6;
      default: break;
      }
      return shapeFunctionIds;
    }
  };
} /* end namespace Kaskade */
/**
 * \endcond
 */
#endif
