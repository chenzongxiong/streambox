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

#ifndef GRIDSCANNER_HH
#define GRIDSCANNER_HH

#include <vector>
#include <string>
#include <fstream>

#include "dune/grid/common/grid.hh"
#include "dune/common/fmatrix.hh"

#include "fem/shapefunctioncache.hh"

namespace Kaskade
{
  /**
   * @file
   * @brief  Two classes for the visualization of higher order grid functions
   * @author Anton Schiela
   *
   */
  /// Class that can sample a Function(View) on a uniform rectangular grid
  /** This is useful, when using higher order FunctionSpaceElement s or highly
   * nonlinear FunctionViews on a relatively coarse grid. The output to
   * PARAVIEW can then be performed in high resolution. The use together with EdgeScanner
   * yields a nice visualization of these functions.
   *
   * Moreover, it is possible to zoom into interesting regions by defining the
   * area of sampling appropriately.
   *
   * Here a sample of how an output of UniformSampler and EdgeScanner may look like:
   *
   * \image html gridscanner.jpg
   *
   */
  template<class Function>
  class UniformSampler {
  public:

    typedef typename Function::Space Space;
    typedef typename Space::Grid Grid;
    static int const spaceDim = Grid::dimension;


    /// Performs a uniform sampling of fu on the cuboid marked by lb, ub with stepsize st
    /**

    @param lb defines the lower bounds of the cuboid
    @param ub defines the upper bounds of the cuboid
    @param st determines the stepsize that is taken for sampling in each direction
    @param fu is a FunctionSpaceElement or a FunctionViews object: the function to be scanned
    @param dataIsScalar determines whether a vectorial or a multicomponent scalar VTK file is written
     */
    UniformSampler(Dune::FieldVector<typename Grid::ctype, spaceDim> const& lb,
                   Dune::FieldVector<typename Grid::ctype, spaceDim> const& ub,
                   Dune::FieldVector<typename Grid::ctype, spaceDim> const& st,
                   Function const& fu,
                   bool dataIsScalar = Function::components==1):
                   lowerbound(lb), upperbound(ub), step(st), scalar(dataIsScalar)
    {
      size=getSize();
      for(int i=0; i < spaceDim; i++)
        step[i]=(upperbound[i]-lowerbound[i])/(size[i]-1);
      scanData(fu);
    }

    /// Print data in .vtk (Legacy) format. The extension ".vtk" is automatically appended
    void writeVTK(std::string name)
    {
      name = name + ".vtr";
      std::ofstream out(name.c_str());
      if(!out)
        throw FileIOException("Failed to open file for writing.",name,__FILE__,__LINE__);
      
      out << "<VTKFile type=\"RectilinearGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n"
          << "  <RectilinearGrid WholeExtent=\"";
      for (int i=0; i<spaceDim; ++i) out << "0 " << size[i]-1 << ' ';
      for (int i=spaceDim; i<3; ++i) out << "0 0 ";
      out << "\">\n"
          << "    <Piece Extent=\"";
      for (int i=0; i<spaceDim; ++i) out << "0 " << size[i]-1 << ' ';
      for (int i=spaceDim; i<3; ++i) out << "0 0 ";
      
      int vtkComponents = Function::components==1? 1: 3;
      
      out << "\">\n"
          << "      <PointData>\n"
          << "        <DataArray type=\"Float32\" format=\"ascii\" Name=\"u\" NumberOfComponents=\"" << vtkComponents << "\">\n";

      // Note that x is innermost loop according to Paraview manual
      for(int i=0; i<data.size(); ++i) {
        for (int j=0; j<std::min(Function::components,vtkComponents); ++j)
          out << static_cast<float>(data[i][j]) << ' ';
        for (int j=std::min(Function::components,vtkComponents); j<vtkComponents; ++j)
          out << "0 ";
        out << std::endl;
      }
          
          
      out << "        </DataArray>\n"
          << "      </PointData>\n"
          << "      <Coordinates>\n";
          
      for (int i=0; i<3; ++i)
      {
        out << "        <DataArray type=\"Float32\" format=\"ascii\" NumberOfComponents=\"1\" Name=\"" << (i==0? 'x': i==1? 'y': 'z') << "\">\n          ";
        if (i<spaceDim)
          for (int j=0; j<size[i]; ++j)
            out << lowerbound[i] + j*step[i] << ' ';
        else
          out << 0;
        out << "\n"
            << "        </DataArray>\n";
      }
      
      out << "      </Coordinates>\n"
          << "    </Piece>\n"
          << "  </RectilinearGrid>\n"
          << "</VTKFile>\n";
    }

  private:

    Dune::FieldVector<int, spaceDim> getSize()
    {
      Dune::FieldVector<int, spaceDim> size;
      for(int i=0; i<spaceDim; i++)
        size[i]=(int)floor((upperbound[i]-lowerbound[i])/step[i])+1;
      return size;
    }


    void scanData(Function const& fu)
    {

      ShapeFunctionCache<Grid,typename Function::Scalar> sfCache;

      Dune::FieldVector<int, spaceDim> ind, lrange, urange;
      int datasize=1;
      for(int i=0; i< spaceDim; i++) datasize *= size[i];
      data.resize(datasize);
      for(int i=0; i<datasize; ++i) data[i] = 0.0;
      typename Space::Evaluator fkt(fu.space(),&sfCache);
      typedef typename Grid::template Codim<0>::Entity Cell;

      for (auto ci=fu.space().gridView().template begin<0>(); ci!=fu.space().gridView().template end<0>(); ++ci)
      {
        fkt.moveTo(*ci);
        Dune::FieldVector<int, spaceDim> at;
        getIndexRange<Cell, Grid::dimension>(lrange,urange,*ci);
        scanLocally<Cell, Grid::dimension>(fu, fkt, data,*ci,lrange,urange,at,spaceDim);
      }
    }

    // Compute the bounding index box of the given cell
    template<class Entity, int spaceDim>
    void getIndexRange(Dune::FieldVector<int, spaceDim> &lower,
                       Dune::FieldVector<int, spaceDim> &upper,
                       Entity const& ci)
    {
      Dune::FieldVector<typename Grid::ctype, spaceDim> ur(lowerbound),lr(upperbound);
      for(int i=0; i<ci.geometry().corners(); i++)
      {
        auto x = ci.geometry().corner(i);
        for(int k=0; k<spaceDim; k++)
        {
          lr[k]= std::min(lr[k],x[k]);
          ur[k]= std::max(ur[k],x[k]);
        }
      }
      for(int i=0; i< spaceDim; i++)
      {
        lower[i] = std::max((int)ceil((lr[i]-lowerbound[i])/step[i]),  0);
        upper[i] = std::min((int)ceil((ur[i]-lowerbound[i])/step[i])+1,size[i]);
      }
    }

    template<int spaceDim>
    Dune::FieldVector<typename Grid::ctype,spaceDim> getCoord(const Dune::FieldVector<int, spaceDim> &ind)
    {
      Dune::FieldVector<typename Grid::ctype,spaceDim> x;
      for(int i=0; i<spaceDim; i++)
        x[i] = lowerbound[i]+step[i]*ind[i];
      return x;
    }


    // Convert integer index coordinate tuple to linear array address.
    // According to the Paraview user's guide, x is the innermost loop.
    int getDataIndex(Dune::FieldVector<int, spaceDim> const& ind)
    {
      int index = 0;
      int mult = 1;
      for(int i=0; i<spaceDim; i++)
      {
        index += ind[i]*mult;
        mult *= size[i];
      }
      return index;
    }

    template<class Entity, int spaceDim>
    void scanLocally(Function const& fu,
                     typename Function::Space::Evaluator& fkt,
                     std::vector<typename Function::ValueType>& data,
                     Entity const& ci,
                     Dune::FieldVector<int, spaceDim> const& lrange,
                     Dune::FieldVector<int, spaceDim> const& urange,
                     Dune::FieldVector<int, spaceDim>& at,
                     int dimnumber)
    {
      if(dimnumber>0)
      {
        for(at[dimnumber-1]=lrange[dimnumber-1]; at[dimnumber-1]<urange[dimnumber-1];  ++at[dimnumber-1])
          scanLocally<Entity,spaceDim>(fu, fkt,data,ci,lrange,urange,at,dimnumber-1);
      }
      else
      {
        auto x = getCoord<spaceDim>(at);
        auto xi = ci.geometry().local(x);
        if (checkInside(ci.geometry().type(),xi) <= 0)
        {
          fkt.evaluateAt(xi);
          int dataInd = getDataIndex(at);
          data[dataInd] = fu.value(fkt);
        }
      }
    }

    Dune::FieldVector<typename Grid::ctype, spaceDim> lowerbound, upperbound, step;
    bool scalar;
    Dune::FieldVector<int, spaceDim> lrange, urange, size;
    std::vector<typename Function::ValueType> data;
  };

  //---------------------------------------------------------------------

  /// Class to scan functions along the edges of a computational grid
  /** This is useful for the visalization of higer order, or nonlinear
    functions on a relatively coarse grid. Together with the
    UniformSampler on can create two PARAVIEW files, which help to
    visualize the obtained results. The UniformSampler provides the
    surface (or volume) plot. The EdgeScanner marks, where the grid
    lines are. 

   * Here a sample of how an output of UniformSampler and EdgeScanner may look like:
   *
   * \image html gridscanner.jpg


    This class is meant to be used on coarse grids. On fine grids this
    kind of visualization is not needed and creates large datasets.
   */
  template<class Function>
  class EdgeScanner
  {
  public:
    /// Scans Function fu along the edges of the computational grid
    /**
    - fu is a FunctionSpaceElement or a FunctionViews object: the function to be scanned

    - nEdgePoints is the number of points used for each edge
     */
    EdgeScanner(Function const& fu, int nEdgePoints)
    {
      scanEdges(fu,nEdgePoints);
    }

    /// Print data in .vtk (Legacy) format. The extension ".vtk" is automatically appended
    void writeVTK(std::string name)
    {
      name = name +".vtk";
      std::ofstream out(name.c_str());
      if(!out)
      {
        std::cout << "Error when opening file!" << std::endl;
        return;
      }
      out << "# vtk DataFile Version 2.0" <<std::endl;
      out << "Scanned Data" <<std::endl;
      out << "ASCII" <<std::endl;
      out << "DATASET POLYDATA" << std::endl;

      int npoints=0;
      for(int i=0; i<data.size(); ++i)
        npoints+=data[i].size();

      out << "POINTS " << npoints << " float" << std::endl;
      for(int i=0; i<data.size(); ++i)
        for(int k=0; k<data[i].size(); ++k)
          out << static_cast<float>(data[i][k][0]) << " "
          << static_cast<float>(data[i][k][1]) << " "
          << static_cast<float>(data[i][k][2]) << std::endl;
      out << "LINES " << data.size() << " " << data.size()+npoints << std::endl;
      int dind=0;
      for(int i=0; i<data.size(); ++i)
      {
        out << data[i].size();
        for(int k=0; k<data[i].size();++k)
        {
          out << " " << dind;
          dind++;
        }
        out << std::endl;
      }
      out << "POINT_DATA " << npoints << std::endl;

      out << "SCALARS data float 1" << std::endl;
      out << "LOOKUP_TABLE default" << std::endl;

      for(int i=0; i<data.size(); ++i)
      {
        for(int k=0; k<data[i].size();++k)
          out << static_cast<float>(data[i][k][3]) << " ";
        out << std::endl;
      }

    }

  private:

    void scanEdges(Function const& fu, int nEdgePoints)
    {
      typedef typename Function::Space Space;
      typedef typename Space::Grid Grid;

      ShapeFunctionCache<Grid,typename Function::Scalar> sfCache;

      typename Space::Evaluator fkt(fu.space(),&sfCache);
      typedef typename Space::IndexSet::template Codim<0>::template Partition<Dune::All_Partition>::Iterator CellIterator;

      Dune::FieldVector<typename Grid::ctype,1> onEdge;

      for (CellIterator ci=fu.space().indexSet().template begin<0,Dune::All_Partition>();
          ci!=fu.space().indexSet().template end<0,Dune::All_Partition>(); ++ci)
      {
        fkt.moveTo(*ci);
        Dune::GenericReferenceElement<typename Grid::ctype,Grid::dimension> const& 
          r(Dune::GenericReferenceElements<typename Grid::ctype,Grid::dimension>::general(ci->type()));
        for(int i=0; i < r.size(Grid::dimension-1); ++i)
        {
          std::vector<Dune::FieldVector<typename Grid::ctype,4> > dcell;
          for(int k=0; k<nEdgePoints; k++)
          {
            onEdge[0]=k*1.0/(nEdgePoints-1);
            Dune::FieldVector<typename Grid::ctype,Grid::dimension> v(r.template global<Grid::dimension-1>(onEdge,i,Grid::dimension-1));
            fkt.evaluateAt(v);
            v=ci->geometry().global(v);
            Dune::FieldVector<typename Grid::ctype,4> dt(0.0);

            for(int j=0; j<Grid::dimension;++j) dt[j]=v[j];
            dt[3]=fu.value(fkt);
            dcell.push_back(dt);
          }
          data.push_back(dcell);
        }
      }
    }

    std::vector< std::vector<Dune::FieldVector<typename Function::Space::Grid::ctype,4> > > data;
  };
} /* end of namespace Kaskade */
#endif
