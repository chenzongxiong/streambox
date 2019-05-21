/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2002-2013 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef IO_TOOLS_HH
#define IO_TOOLS_HH

#include <dune/common/fvector.hh>

#if HAVE_UG
#include <dune/grid/uggrid/uggridfactory.hh>
#endif

#if HAVE_ALBERTA
#include <dune/grid/albertagrid/gridfactory.hh>
#endif

#if ENABLE_ALUGRID
#include <dune/grid/alugrid/3d/alu3dgridfactory.hh>
#endif

namespace Kaskade
{
  namespace IOTools
  {
    /**
     * \cond internals
     */

    // guess what happens here
    inline int copyToInt(Dune::FieldVector<unsigned char,1> t){ return (int)t[0]; }

    // compare floating point vector with double vector
    template <class FloatVec, class DoubleVec>
    bool coincide(FloatVec const& fv, DoubleVec const& dv)
    {
      for(size_t i=0; i<FloatVec::size; ++i)
        if(fabs(fv[i]-dv[i]) > 1.0e-9) return false;

      return true;
    }

    // check if vertex is in intersection
    template <class GlobalCoordinate, class Intersection>
    bool vertexInIntersection(GlobalCoordinate const& vertex, Intersection const& intersection){
      for(int cornerId=0; cornerId<intersection.geometry().corners(); ++cornerId)
        if(coincide(vertex,intersection.geometry().corner(cornerId))) return true;

      return false;
    }

    /// check if boundary segment given through the entries in corners and through an intersection coincide
    template <class GlobalCoordinate, class Intersection>
    bool segmentInIntersection(std::vector<GlobalCoordinate> const& corners, Intersection const& intersection)
    {
      for(int cornerId=0; cornerId<intersection.geometry().corners(); ++cornerId)
        if( !vertexInIntersection(corners[cornerId], intersection) ) return false;

      return true;
    }

    /// Wrapper hiding the differences in the grids constructors;
    template <int dim, typename Grid>
    struct FactoryGenerator;

#ifdef HAVE_UG
    template <int dim>
    struct FactoryGenerator<dim, Dune::UGGrid<dim> >
    {
      typedef Dune::UGGrid<dim> Grid;
      typedef Dune::GridFactory<Grid> Factory;

      static Factory createFactory(unsigned int initialHeapSize=0)
      {
        if(initialHeapSize > 0) Grid::setDefaultHeapSize(initialHeapSize);

        return Factory();
      }
    };
#endif /* HAVE_UG */

#ifdef ENABLE_ALUGRID
    template <int dim>
    struct FactoryGenerator<dim, Dune::ALUSimplexGrid<dim,dim> >
    {
      typedef Dune::ALUSimplexGrid<dim,dim> Grid;
      typedef Dune::GridFactory<Grid> Factory;

      static Factory createFactory(unsigned int){ return Factory(); }
    };

    template <int dim>
    struct FactoryGenerator<dim, Dune::ALUCubeGrid<dim,dim> >
    {
      typedef Dune::ALUCubeGrid<dim,dim> Grid;
      typedef Dune::GridFactory<Grid> Factory;

      static Factory createFactory(unsigned int){ return Factory(); }
    };
#endif /* ENABLE_ALUGRID */

    /// Wrapper for choosing the correct way of getting the boundary segment indices
    template <int dim, class Grid> struct BoundarySegmentIndexWrapper;

#ifdef HAVE_UG
    template <int dim>
    struct BoundarySegmentIndexWrapper<dim, Dune::UGGrid<dim> >
    {
      template <class Factory, class VertexVector, class BoundaryVector, class TmpIdVector, class IdVector>
      static void readBoundarySegmentIndices(Dune::UGGrid<dim> const& grid, Factory const& factory, VertexVector const& vertices, BoundaryVector const& boundary, TmpIdVector const& boundaryIdsFromFile, IdVector& boundaryIds)
      {
        if(boundaryIds.size() < boundaryIdsFromFile.size()) boundaryIds.resize(boundaryIdsFromFile.size());
        std::transform(boundaryIdsFromFile.begin(), boundaryIdsFromFile.end(), boundaryIds.begin(), copyToInt);
      }
    };
#endif /* HAVE_UG */

#ifdef ENABLE_ALUGRID
    template <int dim>
    struct BoundarySegmentIndexWrapper<dim, Dune::ALUSimplexGrid<dim,dim> >
    {
      typedef Dune::ALUSimplexGrid<dim,dim> Grid;

      template <class Factory, class VertexVector, class BoundaryVector, class TmpIdVector, class IdVector>
      static void readBoundarySegmentIndices(Grid const& grid, Factory const& factory, VertexVector const& vertices, BoundaryVector const& boundary, TmpIdVector const& boundaryIdsFromFile, IdVector& boundaryIds)
      {
        typedef typename Grid::LeafGridView::template Codim<0>::Iterator CellIterator;
        typedef typename Grid::LeafGridView::IntersectionIterator FaceIterator;
        CellIterator const cend = grid.leafView().template end<0>();
        for(CellIterator ci=grid.leafView().template begin<0>(); ci != cend; ++ci)
        {
          FaceIterator const fend = grid.leafView().iend(*ci);
          for(FaceIterator fi=grid.leafView().ibegin(*ci); fi!=fend; ++fi)
          {

            if(fi->boundary())
            {
              for(size_t boundarySegment=0; boundarySegment<boundary.size(); ++boundarySegment)
              {
                typedef typename FaceIterator::Intersection::Geometry::GlobalCoordinate GlobalCoordinate;
                std::vector<GlobalCoordinate> corners(dim);
                for(int cornerId=0; cornerId<dim; ++cornerId)
                  for(int entry=0; entry<dim; ++entry)
                    corners[cornerId][entry]= vertices[boundary[boundarySegment][cornerId]][entry];

                if(segmentInIntersection(corners, *fi))
                {
                  size_t index = factory.insertionIndex(*fi);
                  if(boundaryIds.size() < index+1)
                    boundaryIds.resize(index+1);
                  boundaryIds[index] = boundaryIdsFromFile[boundarySegment][0];
                }
              } // end boundarySegment
            } // end if
          } // end FaceIterator
        } // end CellIterator

      }
    };
#endif /* ENABLE_ALUGRID */

    /**
     * \endcond
     */
  } // namespace IOTools
} // namespace Kaskade

#endif // IO_TOOLS_HH
