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

#ifndef INTERPOLATION_TOOLS_HH
#define INTERPOLATION_TOOLS_HH

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>
#include "fem/barycentric.hh"
#include "fem/lagrangespace.hh"
#include "fem/hierarchicspace.hh"
#include "fem/fetransfer.hh"
#include "utilities/member_variable_macro.hh"
#include "utilities/geometry/geomtools.hh"
#include "utilities/detailed_exception.hh"
#include "utilities/linalg/scalarproducts.hh"

/// Necessary information for Dune::GeometryGrid
namespace Dune
{

  namespace Capabilities
  {

    template< class Grid >
    struct hasHierarchicIndexSet
    {
      static const bool v = false;
    };

    template< class Grid >
    struct hasHierarchicIndexSet< const Grid >
    {
      static const bool v = hasHierarchicIndexSet< Grid >::v;
    };

  }

}

namespace Kaskade
{

  namespace InterpolationTools
  {
    /// Get maximum value in a vector.
    /*
     * Throws an exception if an empty vector is passed.
     * Vector entries must be less-than comparable, provide a copy-constructor and the assignment-operator.
     *
     * \param vec input vector
     * \result maximal entry in vec
     */
    template<class Type>//, class ReturnType=Type>
    Type getMaximumValue(std::vector<Type> const& vec)
    {
      if(vec.size() == 0) throw DetailedException("Cannot detect maximum value in empty vector." , __FILE__, __LINE__);
      Type result = vec[0];
      std::for_each(++(vec.begin()), vec.end(), [&](Type const& entry) { if(result < entry) result = entry; });

      return result;
    }

    /// Check if entry is contained in vec.
    template<class Type>
    bool contains(Type const entry, std::vector< Type > const& vec)
    {
      for(size_t i=0; i<vec.size(); ++i)
        if(vec[i] == entry) return true;
      return false;
    }
    /// Container storing normals together with some additional information.
    /*
     * This container may be used to store normals together with additional information that are associated
     * to i.e. a vertex or an edge.
     */
    template<class Scalar_, int dim_>
    struct NormalCollection{
      static int const dim = dim_;
      typedef Scalar_ Scalar;
      typedef Dune::FieldVector<Scalar,dim> Vector;

      NormalCollection() : onBoundary(false), normals(0), weights(0), phaseIds(0), ids(0), isRelevant(true) {}

      NormalCollection(NormalCollection&& other) : onBoundary(false), normals(0), weights(0), phaseIds(0), ids(0), isRelevant(true)
      {
        *this = std::move(other);
      }

      std::ostream& print(std::ostream& os = std::cout) const
      {
        os << std::boolalpha << "Relevant: " << isRelevant << std::endl;
        for(int i=0; i<normals.size(); ++i) os << "normal: " << normals[i] << ", material id: " << phaseIds[i] << std::endl;
        return os;
      }

      void clear() {
        normals.clear();
        weights.clear();
        phaseIds.clear();
        ids.clear();
        isRelevant = false;
      }

      void irrelevant()
      {
        clear();
        normals.push_back(Vector(0));
        weights.push_back(0);
        phaseIds.push_back(-1);
        ids.push_back(-1);
      }

      bool onBoundary;
      std::vector<Vector> normals;
      std::vector<Scalar> weights;
      std::vector<int> phaseIds, ids;
      bool isRelevant;

      // move assignment operator
      NormalCollection& operator=(NormalCollection&& other)
      {
        // steal resources
        onBoundary = other.onBoundary;
        normals = std::move(other.normals);
        weights = std::move(other.normals);
        phaseIds = std::move(other.phaseIds);
        ids = std::move(other.ids);
        isRelevant = other.isRelevant;

        // invalidate other
        other.onBoundary = false;
        other.isRelevant = true;
      }

      /// Print some information on a collection of normals
      void print(std::ostream& os){
        os << std::boolalpha << "is relevant: " << isRelevant << std::endl;
        os << "normals:" << std::endl;
        for(size_t j=0; j<normals.size(); ++j) os << normals[j] << ", phase id: " << phaseIds[j] << std::endl;
      }

      // friends
      friend std::ostream& operator <<(std::ostream &os, NormalCollection<Scalar,dim> const& collection)
      {
        return collection.print(os);
      }
    };

    // compute mean normals
    /**
     * Using the data of NormalCollection this function calculates a mean normal for each material
     * Ignores entities where more than 2 materials meet and which have been marked as irrelevant.
     */
    template <bool inVertex, class NormalCollection>
    void computeMeanNormal(NormalCollection &collection){

      // if vertex does not lie on the boundary: normals.size() == 0
      // if edge has more than 4 normals (more than 2 materials meet at this edge) -> do nothing
      if(collection.normals.size() == 0 || (collection.normals.size() > 4 && !inVertex) ){
        collection.irrelevant();
        return;
      }
      if(!collection.isRelevant)
      {
        collection.irrelevant();
        return;
      }

      // if vertex/edge belongs to two faces of the same entity -> no smoothing
      // TODO delegate to policy
      if(!collection.onBoundary)
      {
        std::vector<int> idCount;
        for(size_t i=0; i<collection.ids.size(); ++i){
          if( InterpolationTools::contains(collection.ids[i], idCount) ){
            collection.irrelevant();
            return;
          }
          else idCount.push_back( collection.ids[i] );
        }
      }
      // else
      // collect material types
      std::vector<bool> materials;
      if(collection.onBoundary){
        materials = std::vector<bool>(1, true);
      }
      else{
        materials = std::vector<bool>( (InterpolationTools::getMaximumValue(collection.phaseIds) + 1), false );
        for(size_t i=0; i<collection.normals.size(); ++i)
          materials[collection.phaseIds[i]] = true;
      }

      // compute mean normals
      std::vector<int> reducedMaterialList;
      std::vector<typename NormalCollection::Vector> reducedNormalList;
      for(size_t i=0; i<materials.size(); ++i){
        // if material does not exist in adjacent cells do nothing
        if(!materials[i]) continue;

        // if material exists store its id
        reducedMaterialList.push_back(i);
        std::vector<typename NormalCollection::Vector> tmp;
        // and collect all corresponding normals
        for(size_t j=0; j<collection.normals.size(); ++j)
          if(collection.phaseIds[j] == i || collection.onBoundary)
            tmp.push_back(collection.weights[j]*collection.normals[j]);

        // compute mean normal
        typename NormalCollection::Vector meanNormal(0);
        for(size_t j=0; j<tmp.size(); ++j)
          meanNormal += tmp[j];

        GeomTools::normalize(meanNormal);

        reducedNormalList.push_back(meanNormal);
      }

      collection.normals = reducedNormalList;
      collection.phaseIds = reducedMaterialList;

      // if intersection with other inner boundary layers
      if(collection.phaseIds.size() > 2) collection.irrelevant();

    }

    // -------------------------------------------------------------------------------------------------
    // Choose function space / shape function set
    // -------------------------------------------------------------------------------------------------
    /**
     * \cond internal
     */
    /// Choose a shape function set according to the grid's dimension.
    template<class Scalar, class Grid, int dim>
    struct ChoiceOfShapeFunctionSet
    {
      typedef typename std::conditional<
          dim==3,
          HierarchicSimplexShapeFunctionSet<Grid,Scalar>,
          LagrangeSimplexShapeFunctionSet<Grid,Scalar>
      >::type type;

      ChoiceOfShapeFunctionSet() : shapeFunctionSet(3){}

      type const shapeFunctionSet;
    };

    // -------------------------------------------------------------------------------------------------
    // Data Container.
    // -------------------------------------------------------------------------------------------------

    struct Empty
    {
      Empty() = default;

      Empty(const Empty&) = delete;
      Empty& operator=(const Empty&) = delete;
    };
    /**
     * \endcond
     */


    /// Simple data container holding only vertex/edge normals and a threshold value for corner detection.
    template<class Scalar, int dim>
    struct NormalContainer;

    /// 2D specialization.
    template<class Scalar>
    struct NormalContainer<Scalar, 2> : public Empty {
      typedef InterpolationTools::NormalCollection<Scalar,2> EntryType;
      typedef std::vector<EntryType> CollectionType;

      explicit NormalContainer(size_t numberOfVertices=0) : vertexNormals(numberOfVertices) {}
      
      template <class GridView>
      explicit NormalContainer(GridView const& gridView) : vertexNormals(gridView.size(GridView::dimension))
      {}

      NormalContainer(NormalContainer&& other) : vertexNormals(std::move(other.vertexNormals)) {}

      NormalContainer& operator=(NormalContainer&& other)
      {
        vertexNormals = std::move(other.vertexNormals);
      }

      std::ostream& print(std::ostream& os) const
      {
        int id = 0;
        for(EntryType const& entry : vertexNormals)
        {
          std::cout << "entry nr. " << id++ << std::endl;
          entry.print(os);
        }
        return os;
      }

      friend std::ostream& operator<<(std::ostream& os, NormalContainer const& container) { return container.print(os); }

      CollectionType vertexNormals;
    };

    /// 3D specilization
    template<class Scalar>
    struct NormalContainer<Scalar, 3> : public Empty{
      typedef InterpolationTools::NormalCollection<Scalar,3> EntryType;
      typedef std::vector< EntryType > CollectionType;

      explicit NormalContainer(size_t numberOfVertices=0, size_t numberOfEdges=0) : vertexNormals(numberOfVertices), edgeNormals(numberOfEdges) {}

      template <class GridView>
      explicit NormalContainer(GridView const& gridView) : vertexNormals(gridView.size(GridView::dimension)), edgeNormals(gridView.size(GridView::dimension-1))
      {}
      
      NormalContainer(NormalContainer&& other) : vertexNormals(std::move(other.vertexNormals)), edgeNormals(std::move(other.edgeNormals)) {}

      NormalContainer& operator=(NormalContainer&& other)
      {
        vertexNormals = std::move(other.vertexNormals);
        edgeNormals = std::move(other.edgeNormals);
      }

      std::ostream& print(std::ostream& os) const
      {
        int id = 0;
        for(EntryType const& entry : vertexNormals)
        {
          std::cout << "entry nr. " << id++ << std::endl;
          entry.print(os);
        }
        id = 0;
        for(EntryType const& entry : edgeNormals)
        {
          std::cout << "entry nr. " << id++ << std::endl;
          entry.print(os);
        }
        return os;
      }

      friend std::ostream& operator<<(std::ostream& os, NormalContainer const& container) { return container.print(os); }

      CollectionType vertexNormals, edgeNormals;
    };
  } /* end of namespace InterpolationTools 1 */

  namespace Policy{

    /**************************************************************************************************/
    /******************************* Boundary Policies - Inner Boundary *******************************/
    /**************************************************************************************************/

    /// Default policy for inner boundaries. Ignores them.
    struct IgnoreInnerBoundary{
      bool ignoreInnerBoundary(int tissue, int tissueNeighbour) const { return true; }
    };

    /// Policy for inner boundaries. Considers all of them.
    struct ConsiderInnerBoundary{
      bool ignoreInnerBoundary(int tissue, int tissueNeighbour) const { return false; }
    };

    /// Allows to specify parts of the inner boundary. Inner boundaries are characterized by a change of the material id.
    /*
     * \param IdVector vector containing the boundary ids (i.e. std::vector<int>)
     * \param skipSpecifiedIds bool value indicating default treatment (i.e. true: consider only specified parts of the boundary, false: ignore specified parts)
     */
    template<class IdVector, bool skipSpecifiedIds = false>
    struct SpecifyInnerBoundary
    {
      /*
       * \param boundaryIds vector of boundary ids
       * \param specifiedIds vector of ids that will be ignored or considered only
       */
      SpecifyInnerBoundary(std::vector<std::pair<int,int> > specifiedInterfaces)
      {
        for(size_t i=0; i<specifiedInterfaces.size(); ++i) addInterface( specifiedInterfaces[i] );
      }

      /// Add interface specified by two material ids
      void addInterface(std::pair<int,int> ids){
        if(ids.first == ids.second) std::cout << "Ignoring interface between same material ids (id: " << ids.first << ")" << std::endl;
        else specifiedInterfaces_.push_back(std::make_pair(std::min(ids.first,ids.second), std::max(ids.first,ids.second)));
      }

      /// Add interface specified by two material ids
      void addInterface(int first, int second){
        addInterface(std::make_pair(first,second));
      }

      /// Check if interface specified by two material ids should be ignored
      bool ignoreInnerBoundary(int tissue, int tissueNeighbour) const
      {
        if(tissue == tissueNeighbour) return true;
        int const tmp1 = std::min(tissue, tissueNeighbour),
            tmp2 = std::max(tissue, tissueNeighbour);
        bool contains = std::find( specifiedInterfaces_.begin(), specifiedInterfaces_.end(), std::make_pair(tmp1, tmp2) ) != specifiedInterfaces_.end();
        return (contains == skipSpecifiedIds);
      }

    private:
      std::vector<std::pair<int,int> > specifiedInterfaces_;
    };

    /**************************************************************************************************/
    /******************************* Boundary Policies - Outer Boundary *******************************/
    /**************************************************************************************************/

    /// Default policy for outer boundary. Ignores them.
    struct IgnoreOuterBoundary{
      template <class Face>
      bool ignoreOuterBoundary(Face const&) const { return true; }
    };

    /// Policy for inner boundaries. Considers all of them.
    struct ConsiderOuterBoundary{
      template <class Face>
      bool ignoreOuterBoundary(Face const&) const { return false; }
    };

    /// Allows to specify parts of the outer boundary.
    /*
     * \param IdVector vector containing the boundary ids (i.e. std::vector<int>)
     * \param skipSpecifiedIds bool value indicating default treatment (i.e. true: consider only specified parts of the boundary, false: ignore specified parts)
     */
    template<class IdVector, bool skipSpecifiedIds = false>
    struct SpecifyOuterBoundary
    {
      /*
       * \param boundaryIds vector of boundary ids
       * \param specifiedIds vector of ids that will be ignored or considered only
       */
      SpecifyOuterBoundary(IdVector const& boundaryIds, IdVector const& specifiedIds) :
        boundaryIds_(boundaryIds), considerFace_(InterpolationTools::getMaximumValue(specifiedIds) + 1, skipSpecifiedIds)
      {
        for(size_t i=0; i<specifiedIds.size(); ++i)
          considerFace_[specifiedIds[i]] = !skipSpecifiedIds;
      }

      /// Check if outer bounday face should be ignored
      template<class Face>
      bool ignoreOuterBoundary(Face const& face) const
      {
        unsigned int index = face.boundarySegmentIndex();
        if(index < boundaryIds_.size()) return considerFace_[boundaryIds_[index]];
        else return !skipSpecifiedIds;

        return false;
      }

    private:
      IdVector const& boundaryIds_;
      std::vector<bool> considerFace_;
    };

    /**************************************************************************************************/
    /******************************* Interpolation Policies - Threshold *******************************/
    /**************************************************************************************************/

    template <class Scalar>
    struct NoGradientThreshold
    {
      void applyGradientThreshold(Scalar &gradient) const {}
    };

    /// Threshold for feature detection.
    /**
     * This is a simple threshold strategy which works on 1-dimensional gradients, i.e. deformations defined along edges, edge normals in the neighbouring surface triangles,...
     * If the threshold value is negative it will be ignored!!!
     */
    template <class Scalar>
    struct ApplyGradientThreshold
    {
      explicit ApplyGradientThreshold(Scalar const threshold=1.0) : threshold_(threshold){}

      void applyGradientThreshold(Scalar &gradient) const
      {
        if(threshold_ < 0) return;
        if(fabs(gradient) > threshold_)
          gradient = 0;
      }

    private:
      Scalar threshold_;
    };

    /**************************************************************************************************/
    /**************************************** Material Policy *****************************************/
    /**************************************************************************************************/
    /**
     * \cond internals
     */
    /// Automatically chosen policy in BoundaryNormalCollector and HermiteInterpolation.
    template <class Scalar>
    struct NoPhaseInfo
    {
      NoPhaseInfo(){}

      template <class PhaseElement>
      explicit NoPhaseInfo(PhaseElement const&){}

      template <class Entity>
      Scalar phaseId(Entity const& entity) const { return 0; }
    };

    /// Automatically chosen policy in BoundaryNormalCollector and HermiteInterpolation.
    /**
     * The phases are identified via ids stored in a FSElement.
     */
    template<class PhaseElement>
    struct PhaseAsFSElement
    {
      explicit PhaseAsFSElement(PhaseElement const& phase_) : phase(phase_){}

      PhaseAsFSElement(PhaseAsFSElement const& other) : phase(other.phase){}

      template<class Entity>
      typename PhaseElement::Scalar phaseId(Entity const& entity) const
      {
        return phase.value(entity, Dune::FieldVector<typename PhaseElement::Scalar,PhaseElement::Space::dim>(0));
      }

      PhaseElement const& phase;
    };

    /**************************************************************************************************/
    /************************************* ShapeFunctionSet Policy ************************************/
    /**************************************************************************************************/

    /// Automatically chosen policy in InterpolationPolynomialCollection.
    template <class InterpolationPolynomial>
    struct NoShapeFunctionSet{
      typedef typename InterpolationPolynomial::GridView::Grid::template Codim<0>::Entity Entity;

      template <typename... Parameters>
      InterpolationPolynomial init(Entity const& entity, typename InterpolationPolynomial::GridView const& gridView, Parameters const&... parameters)
      {
        return InterpolationPolynomial(entity, gridView, parameters...);
      }

      typename InterpolationPolynomial::range_type
      evaluate(InterpolationPolynomial& polynomial, typename InterpolationPolynomial::domain_type const& x) const
      {
        return polynomial.evaluate(x);
      }
    };

    /// Automatically chosen policy in InterpolationPolynomialCollection.
    template <class InterpolationPolynomial>
    struct UseShapeFunctionSet
    : private InterpolationTools::ChoiceOfShapeFunctionSet<typename InterpolationPolynomial::Scalar, typename InterpolationPolynomial::GridView::Grid, InterpolationPolynomial::GridView::dimension>
    {
      typedef typename InterpolationPolynomial::GridView::Grid::template Codim<0>::Entity Entity;

      template <typename... Parameters>
      InterpolationPolynomial init(Entity const& entity, typename InterpolationPolynomial::GridView const& gridView, Parameters const&... parameters)
      {
        return InterpolationPolynomial(entity, gridView, shapeFunctionSet, parameters...);
      }

      typename InterpolationPolynomial::range_type
      evaluate(InterpolationPolynomial const& polynomial, typename InterpolationPolynomial::domain_type const& x) const
      {
        return polynomial.evaluate(x, shapeFunctionSet);
      }

    private:
      using InterpolationTools::ChoiceOfShapeFunctionSet<typename InterpolationPolynomial::Scalar,typename InterpolationPolynomial::GridView::Grid,InterpolationPolynomial::GridView::dimension>::shapeFunctionSet;
    };

    // Create struct that checks if a type (i.e. an interpolation polynomial) has a const bool variable
    // called 'needsShapeFunctionSet'
    namespace Detail{
      KASKADE_CREATE_MEMBER_VARIABLE_CHECK(bool const, needsShapeFunctionSet, Has_needsShapeFunctionSet)

            template <class InterpolationPolynomial, bool has_needsShapeFunctionSet>
      struct ChooseShapeFunctionSetPolicy
      {
        typedef NoShapeFunctionSet<InterpolationPolynomial> type;
      };

      template <class InterpolationPolynomial>
      struct ChooseShapeFunctionSetPolicy<InterpolationPolynomial,true>
      {
        // use UseShapeFunctionSet if InterpolationPolynomial::needsShapeFunctionSet=true
        typedef typename std::conditional<InterpolationPolynomial::needsShapeFunctionSet,UseShapeFunctionSet<InterpolationPolynomial>,NoShapeFunctionSet<InterpolationPolynomial> >::type type;
      };
    }

    template <class InterpolationPolynomial>
    struct ShapeFunctionSetPolicy : public Detail::ChooseShapeFunctionSetPolicy<InterpolationPolynomial,Detail::Has_needsShapeFunctionSet<InterpolationPolynomial>::value>::type
    {};
    /**
     * \endcond
     */

    template <class Scalar>
    struct RelativeDeformation
    {
      explicit RelativeDeformation(Scalar scale_=0.25) : scale(scale_)
      {}

      template <class Entity>
      void init(Entity const& entity)
      {
        initialVolume = entity.geometry().volume();
      }

      template <class Entity>
      bool operator()(Entity const& entity, Scalar referenceVolume = -1.) const
      {
        return (entity.geometry().volume() < ((referenceVolume<0) ? scale*initialVolume : referenceVolume ));
      }

    private:
      Scalar scale, initialVolume;
    };

    template <class Scalar>
    struct ConstantDamping
    {
      explicit ConstantDamping(Scalar const dampingFactor_ = 0.9) : dampingFactor(dampingFactor_), exponent(0)
      {
        assert(dampingFactor > 0);
      }

      void increase() { (dampingFactor < 1) ? --exponent : ++exponent; }

      void decrease() { (dampingFactor < 1) ? ++exponent : --exponent; }

      void reset() { exponent = 0; }

      Scalar operator()() const
      {
        return 1;//pow(dampingFactor,exponent);
      }

    private:
      Scalar const dampingFactor;
      int exponent;
    };

    /**************************************************************************************************/
    /**************************************** Merging Policies ****************************************/
    /**************************************************************************************************/
    template <class... Policies> struct MergedPolicy : public Policies...
    {
      explicit MergedPolicy(const Policies&... policies) : Policies(policies)...{}
    };

  } /* end of namespace Policy */


  // -------------------------------------------------------------------------------------------------
  // Namespace hiding tools.
  // -------------------------------------------------------------------------------------------------
  namespace InterpolationTools
  {
    /// Get global vertex ids of the entities corners
    /*
     * \param entity entity
     * \param gridView GridView containing the entity
     * \result global indices of the entity's vertices
     */
    template <class Entity, class GridView>
    std::vector<int> getGlobalVertexIds(Entity const &entity, GridView const &gridView)
    {
      std::vector<int> result(GridView::dimension+1,-1);
      for(size_t i=0; i<result.size(); ++i) result[i] = gridView.indexSet().subIndex(entity,i,GridView::dimension);

      return result;
    }

    /// Get local coordinates of the edges vertices in the codim<0> entity.
    /*
     * \param entity entity containing the edge
     * \param edgeId local index of the edge in the entity
     * \result local indices of the edge's vertices wrt entity
     */
    template <class Entity>
    std::vector<int> getEdgeCornerIdsInEntity(Entity const& entity, int edgeId)
    {
      Dune::GeometryType const gt(entity.type().id(), Entity::dimension);
      std::vector<int> result(2,0);
      // get local indices in cell
      result[0] = Dune::GenericReferenceElements<typename Entity::ctype,Entity::dimension>::general(gt).subEntity(edgeId,Entity::dimension-1,0,Entity::dimension);
      result[1] = Dune::GenericReferenceElements<typename Entity::ctype,Entity::dimension>::general(gt).subEntity(edgeId,Entity::dimension-1,1,Entity::dimension);
      return result;
    }

    /// Get global indices of the edges vertices.
    /*
     * \param entity entity containing the edge
     * \param edgeId local id in entity
     * \param gridView gridView containing the entity and edge
     * \return global indices of the edges vertices
     */
    template <class Entity, class GridView>
    std::vector<int> getGlobalEdgeCornerIds(Entity const& entity, int edgeId, GridView const& gridView)
    {
      std::vector<int> localIds = getEdgeCornerIdsInEntity(entity, edgeId);
      std::vector<int> result(2,0);
      result[0] = gridView.indexSet().subIndex(entity,localIds[0],GridView::dimension);
      result[1] = gridView.indexSet().subIndex(entity,localIds[1],GridView::dimension);
      return result;
    }

    /// Get the normals assoicated to the global vertex indices in vid.
    /*
     * \param vid global vertex indices
     * \param vertexNormals global collection of directions associated to the vertices
     * \param localNormals directions corresponding to the vertices indexed through vid
     */
    template <class NormalCollection, class Vector>
    void getVertexNormals(std::vector<int> const &vid, NormalCollection const& vertexNormals, std::vector<Vector> &localNormals, typename Vector::field_type const tissueId)
    {
      for(int i=0; i<vid.size(); ++i){
        for(int j=0; j<vertexNormals[vid[i]].normals.size(); ++j)
          if(fabs(vertexNormals[vid[i]].phaseIds[j] - tissueId ) < 1e-9){
            localNormals[i] = vertexNormals[vid[i]].normals[0];
            break;
          }
      }
    }

  } /* end namespace InterpolationTools 2 */


  template <class FaceIterator, class Scalar, int dim>
  Dune::FieldVector<Scalar,dim+1> getWeights(FaceIterator& iter, FaceIterator& iend, Dune::FieldVector<Scalar,dim> const& x, std::vector<Dune::FieldVector<Scalar,dim> >& normals)
  {
    Dune::FieldVector<Scalar,dim+1> result;
    size_t counter = 0;
    LinAlg::EuclideanScalarProduct scalarProduct;
    for(; iter!=iend; ++iter, ++counter)
    {
      normals[counter] = -1.*iter->centerUnitOuterNormal();
      result[counter] = -1.*scalarProduct(x - iter->geometry().corner(0), normals[counter]);
      assert(result[counter] <= 0);
    }
    return result;
  }

  template <class Vector>
  inline typename Vector::field_type dist(Vector const& point, Vector const& start, Vector const& end)
  {
    typedef typename Vector::field_type Scalar;
    Vector edge =  end - start;
    Vector p = point - start;
    GeomTools::normalize(edge);

    Scalar alpha = p*edge;
    edge *= alpha;
    p -= edge;

    return sqrt(p*p);
  }

  template <class Entity, class Scalar, int dim>
  inline Scalar dist(Entity const& entity, Dune::FieldVector<Scalar,dim> const& x)
  {
    Scalar d = dist(x, entity.geometry().corner(0), entity.geometry().corner(1));
    Scalar tmp = dist(x, entity.geometry().corner(1), entity.geometry().corner(2));
    if(tmp < d) d = tmp;
    tmp = dist(x, entity.geometry().corner(2), entity.geometry().corner(3));
    if(tmp < d) d = tmp;
    return d;
  }

  template <class Scalar>
  inline Scalar invert(Scalar const entry){ return -1./entry; }

  template <class Scalar>
  inline Scalar invert2(Scalar const entry){ return 1./(entry*entry); }

  template <class Scalar>
  inline Scalar invert3(Scalar const entry){
    Scalar div = sqrt(pow(-1.*entry,1.8));
    return (div > 1e-9) ? 1./div : 1e9;
    //  1./sqrt(-1.*pow(entry,1.3));
  }

  template <class Vector>
  inline Vector getLinearWeights(Vector const& weights, typename Vector::field_type (*fun)(typename Vector::field_type const))
  {
    Vector result(0);
    for(int i=0; i<weights.size(); ++i) result[i] = 1.0 - weights[i];
    return GeomTools::normalize(result);
  }

  template <class Vector>
  inline Vector getNormalizedWeights(Vector const& weights, typename Vector::field_type (*fun)(typename Vector::field_type const))
  {
    Vector result(0);
    for(int i=0; i<weights.size(); ++i) result[i] = fun(weights[i]);
    return GeomTools::normalize(result);
  }

  ///
  template <class GridView, class Entity, class Vector, class Deformation>
  Vector interpolateDeformation(GridView const& gridView, Entity const& entity, Vector const& x_local, Deformation const& deformation)
  {
    typedef typename Vector::field_type Scalar;
    typedef typename GridView::IntersectionIterator FaceIterator;
    int constexpr dim = GridView::dimension;
    std::vector<Vector> normals(dim+1);
    FaceIterator iter = gridView.ibegin(entity),
        iend = gridView.iend(entity);
    Dune::FieldVector<Scalar, dim+1> weights= getWeights<FaceIterator>(iter, iend, entity.geometry().global(x_local), normals);
    std::vector<Vector> localEvaluationPoints(dim+1);
    for(int i=0; i<dim+1; ++i) localEvaluationPoints[i] = entity.geometry().local(entity.geometry().global(x_local) + (weights[i])*normals[i]);

    GeomTools::normalize(weights);

    Dune::FieldVector<Scalar, dim+1> normalizedWeights = getNormalizedWeights(weights, invert);//getLinearWeights(weights, invert3);//GeomTools::getNormalized(weights);
    Vector result(0);
    for(size_t i=0; i<localEvaluationPoints.size(); ++i) result += normalizedWeights[i] * deformation.value(entity, localEvaluationPoints[i]);
    return result;
  }

  template <class Interpolation>
  struct Distance
  {
    explicit Distance(Interpolation const& interpolation_)
    : interpolation(interpolation_)
    {}

    template <class Cell, class Vector>
    typename Vector::field_type operator()(Cell const& cell, Vector const& x) const
    {
      Vector const& c0 = cell.geometry().corner(0),
          c1 = cell.geometry().corner(1),
          c2 = cell.geometry().corner(2);

      Vector pos = 0.5*(c0 + c1);
      Vector def = interpolation.value(cell,cell.geometry().local(pos));
      if(def*def > 1e-9) return dist(x,c0,c1);
      pos = 0.5*(c1 + c2);
      def = interpolation.value(cell,cell.geometry().local(pos));
      if(def*def > 1e-9) return dist(x,c1,c2);
      pos = 0.5*(c2 + c0);
      def = interpolation.value(cell,cell.geometry().local(pos));
      if(def*def > 1e-9) return dist(x,c2,c0);
    }

    template <class Cell, class Vector=Dune::FieldVector<double,Cell::dimension> >
    Vector normal(Cell const& cell) const
    {
      Vector const& c0 = cell.geometry().corner(0),
          c1 = cell.geometry().corner(1),
          c2 = cell.geometry().corner(2);

      Vector pos = 0.5*(c0 + c1);
      Vector def = interpolation.value(cell,cell.geometry().local(pos));
      if(def*def > 1e-9) return getNormal(c0,c1);
      pos = 0.5*(c1 + c2);
      def = interpolation.value(cell,cell.geometry().local(pos));
      if(def*def > 1e-9) return getNormal(c1,c2);
      pos = 0.5*(c2 + c0);
      def = interpolation.value(cell,cell.geometry().local(pos));
      if(def*def > 1e-9) return getNormal(c2,c0);
    }

  private:
    template <class Vector>
    Vector getNormal(Vector const& start, Vector const& end) const
    {
      Vector result(end);
      result -= start;
      auto tmp = result[0];
      result[0] = -result[1];
      result[1] = tmp;
      return GeomTools::normalize(result);
    }

    Interpolation const& interpolation;
  };

  template <class Distance, class Cell, class Vector, class Scalar=typename Vector::field_type>
  void relaxVertex(Distance const& distance, Cell const& cell, Vector& x, Scalar const& h, Scalar eta = 0.3, size_t rad=3)
  {
    Scalar r = rad*sqrt(h*h);

    Scalar alpha = eta*h*(1+distance(cell,x)/r);

    Vector dir = distance.normal(cell);
    dir *= alpha;
    x -= dir;
  }

  template <class Grid, class Deformation, class Vector, class DeformationPolicy, class DampingPolicy>
  bool adjustDampingFactor(Grid& grid, Deformation const& deformation, Vector const& initialVolume, DeformationPolicy const& deformationPolicy, DampingPolicy& dampingPolicy)
  {
    typedef typename Grid::ctype ctype;
    typedef typename Grid::template Codim<0>::EntityPointer CellPointer;
    //  constexpr int dim = Grid::dimension;

    // iterate over cells
    std::for_each(grid.leafView().template begin<0>(), grid.leafView().template end<0>(), [&](typename Grid::LeafGridView::template Codim<0>::Iterator::Entity const& entity)
        {
      LocalGeometryInCoarseGridAncestor<Grid> lgcga(grid, CellPointer(entity));
      typedef decltype(entity) EntityRef;
      typedef typename std::remove_reference<EntityRef>::type Entity;
      constexpr int dim = Entity::Geometry::dimension;

      for(int i=0; i<entity.geometry().corners(); ++i)
      {
        typename Dune::template UGGridFamily<dim, dim>::Traits::template Codim<dim>::EntityPointer ep(entity.template subEntity<dim>(i));

        Dune::FieldVector<ctype,dim> p(Dune::template GenericReferenceElements<ctype,dim>::general(entity.type()).position(i,dim));
        Dune::FieldVector<ctype,dim> v(lgcga.global(p));
        Dune::FieldVector<ctype,dim> w(lgcga.getFather()->geometry().global(v));

        Dune::FieldVector<ctype,dim+1> baryInFatherCell = barycentric(v);
        bool inCell = true;
        for(int i=0; i<dim+1; ++i) if(std::abs(baryInFatherCell[i]) < 1e-9) inCell = false;

        w = lgcga.getFather()->geometry().global(v);
        Distance<Deformation> distance(deformation);
        if(!inCell)
          w += dampingPolicy() * deformation.value(*(lgcga.getFather()),v);
        else w += dampingPolicy() * interpolateDeformation(grid.levelView(0), *(lgcga.getFather()), v, deformation);
        //else
          //{
          //  relaxVertex(distance,*(lgcga.getFather()),w,1);
        //}
        grid.setPosition(ep,w);
      }
        });

    ctype minvol = 1000;
    ctype minrel = 1000;
    auto cend = grid.leafView().template end<0>();
    for(auto ci = grid.leafView().template begin<0>(); ci!=cend; ++ci)
    {
      if(ci->geometry().volume()/initialVolume[grid.leafView().indexSet().index(*ci)] < minrel) minrel = ci->geometry().volume()/initialVolume[grid.leafView().indexSet().index(*ci)];
      if(ci->geometry().volume() < minvol) minvol = ci->geometry().volume();
      if(deformationPolicy(*ci, 0.3*initialVolume[grid.leafView().indexSet().index(*ci)])) return false;
    }

    std::cout << "Minvol: " << minvol << std::endl;
    std::cout << "Minrel: " << minrel << std::endl;

    return true;
  }

  template<class Grid, class DefFunction, template <class> class AcceptDeformationPolicy=Policy::RelativeDeformation, template <class> class DampingPolicy = Policy::ConstantDamping >
  void deformGrid(Grid& grid, DefFunction const& deformation, bool allCells, bool doreset,
      AcceptDeformationPolicy<typename Grid::ctype> deformationPolicy = AcceptDeformationPolicy<typename Grid::ctype>(),
      DampingPolicy<typename Grid::ctype> dampingPolicy = DampingPolicy<typename Grid::ctype>())
  {
    typedef typename Grid::LeafGridView::template Codim<0>::Iterator CellIterator;

    // store initial cell volumes
    std::vector<typename Grid::ctype> initialVolume(grid.leafView().indexSet().size(0),0);
    std::for_each(grid.leafView().template begin<0>(), grid.leafView().template end<0>(), [&initialVolume,&grid](typename CellIterator::Entity const& entity)
        {
      for(int i=0; i<entity.geometry().corners(); ++i) initialVolume[grid.leafView().indexSet().index(entity)] = entity.geometry().volume();
        });

    std::cout << "Initial damping factor: " << dampingPolicy() << std::endl;
    while(!adjustDampingFactor(grid,deformation,initialVolume,deformationPolicy,dampingPolicy))
    {
      break;
      dampingPolicy.decrease();
      std::cout << "Damping factor: " << dampingPolicy() << std::endl;
    }
  }

} /* end of namespace Kaskade */

#endif
