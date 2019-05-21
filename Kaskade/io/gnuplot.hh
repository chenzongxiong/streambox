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

#ifndef GNUPLOT_HH
#define GNUPLOT_HH

#include "gnuplotwriter.hh"

#include "io/iobase.hh"

/**
 * @file
 * @brief  Output of mesh and solution for visualization software
 * @author Martin Weiser, Bodo Erdmann, Lutz Weimann
 *
 * This file contains functionality needed to export mesh data and
 * solution data.
 *
 * - \a writeGnuplotFile: write data to Gnuplot Data file, e.g. for visualization with Gnuplot software
 * - \a TimeSeriesGnuplotWriter: A class for writing a (time-) sequence of variable sets into Gnuplot
 * XML files and creating a Paraview description file for playing the
 * animation. The Paraview description file is completed and closed on
 * destruction of the writer object.

 *
 */

//---------------------------------------------------------------------

namespace Kaskade
{
/**
 * \brief Default manipulator for Gnuplot output.
 * 
 * This manipulator simply puts out the function value for functions with 1 to 3 components. As Paraview cannot handle
 * higher number of components, we write just the euclidean norm for higher-dimensional functional values.
 */
template <int components>
struct GnuplotIdManipulator 
{
  int ncomps() const { return components>3? 3: components; }

  template <class Scalar>
  Scalar operator()(int comp, Dune::FieldVector<Scalar,components> const& val) const { 
    if (components <= 3)
      return val[comp]; 
    else
      return val.two_norm();
  }
};

/**
 * \ingroup IO
 * 
 * A wrapper class that implements the Dune::GnuplotWriter::GnuplotFunction
 * interface for finite element functions. 
 */
template <class F, class GridView, class Manipulator=GnuplotIdManipulator<F::components> >
class GnuplotFunction: public GnuplotWriter<GridView>::GnuplotFunction
{
  typedef typename F::Space::Grid Grid;
  
public:
  /**
   * The specified function name \a name_ is modified by erasing all
   * non-alphanumerical characters.
   */
  GnuplotFunction(F const& f_, std::string const& name_, Manipulator const& manipulator=Manipulator()):
    f(f_), manipulator_(manipulator)
  {
    // At least for paraview, the name of the function should not
    // contain operator symbols such as =,-,+,/ etc.
    for (std::string::const_iterator i=name_.begin(); i!=name_.end(); ++i)
      if (std::isalnum(*i))
        myname += *i;
  }
  
  virtual int ncomps () const { return manipulator_.ncomps(); }


  virtual double evaluate (int comp, const typename Grid::Traits::template Codim<0>::Entity& e,
                           const Dune::FieldVector<typename Grid::ctype,Grid::dimension>& xi) const {
    // with -O2/-O3 this seems to be miscompiled sometimes. Adding this brute force output helps for some weird reason.
    // std::cout << myname << " at " << xi << ": " << f.value(e,xi) << " -> " << manipulator_(comp,f.value(e,xi)) << '\n';
    
    return manipulator_(comp,f.value(e,xi));
  }
  
 
  // get name
  virtual std::string name () const { return myname; }
    

private:
  F const& f;
  std::string myname;
  Manipulator manipulator_;
};

//---------------------------------------------------------------------

/**
 * \cond internals
 */
namespace IoDetail {

template <class GridView, class Filter, class RAIter>
struct AddDataTo<GridView,GnuplotWriter<GridView>,Filter,RAIter>
{
  AddDataTo(GridView const& gridView_,GnuplotWriter<GridView>& writer_, Filter const& filter_, RAIter const& names_):
    gridView(gridView_),
    writer(writer_),
    filter(filter_),
    names(names_)
  {}
  
  template <class Pair>
  void operator()(Pair const& pair) const 
  {
    using namespace boost::fusion;

    int const id = boost::remove_reference<typename result_of::value_at_c<Pair,0>::type>::type::id;

    if (!filter(names[id]))
      return;
    
    typedef typename boost::remove_reference<typename result_of::value_at_c<Pair,1>::type>::type Function;
    
    GnuplotFunction<Function,GridView>* Gnuplotf = new GnuplotFunction<Function,GridView>(at_c<1>(pair),names[id]);
    if (Function::Space::continuity >= 0)
      writer.addVertexData(Gnuplotf);
    else
      writer.addCellData(Gnuplotf);
  }

private:
  GridView const& gridView;
  GnuplotWriter<GridView>&    writer;
  Filter const& filter;
  RAIter        names;
};


} // End of namespace IoDetail

/**
 * \endcond
 */

/**
 *   This procedure writes data to Gnuplot file in ascii format, e.g. for visualization in Gnuplot software.
 *   @param gridView gridView, i.e. LeafGridView 
 *   @param description variable set description 
 *   @param vars variables in a variable set
 *   @param filename name of output file 
 *   @param options.info defines the amount of info-output, either none or summary or detail is possible
 *   @param order linear (order=1) or quadratic (order=2) output, default value is linear
 */
template <class GridView,class VariableSet>
void writeGnuplotFile(GridView const& gridView,
                  VariableSet const& vars,
                  std::string const& filename,
                  IoOptions options=ioOptions_default,
		  int order = 1) 
{
  GnuplotOptions::Order od;
  if (order==1) {
    od = GnuplotOptions::linear; }
  else if (order==2) {
    od = GnuplotOptions::quadratic; }
  else {
    std::cout << "order > 2 not supported" << std::endl;
    std::cout << "order=1 is chosen instead" << std::endl;
    od = GnuplotOptions::linear;
  }
  GnuplotOptions::Info inf;
  if (options.info==IoOptions::none) {
    inf = GnuplotOptions::none;}
  else if (options.info==IoOptions::summary) {
    inf = GnuplotOptions::summary; }
  else if (options.info==IoOptions::detail) {
    inf = GnuplotOptions::detail; }
  else {
    inf = GnuplotOptions::none;
  };
    
  typedef GnuplotWriter<GridView> GnuplotWriter;
  GnuplotWriter writer(gridView, od, inf);

  writePartialFile(gridView,writer,vars,filename,UnaryTrue(),inf);
}

template <class GridView, class F>
void writeGnuplot(std::string const& filename, GridView gridView, F const& f, int order=1) 
{
  GnuplotOptions::Order od;
  if (order==1) 
    od = GnuplotOptions::linear; 
  else if (order==2) 
    od = GnuplotOptions::quadratic; 
  else {
    std::cerr << "order > 2 not supported" << std::endl;
    std::cerr << "order=1 is chosen instead" << std::endl;
    od = GnuplotOptions::linear;
  }

  GnuplotWriter<GridView> writer(gridView,od);
  GnuplotFunction<F,GridView>* Gnuplotf = new GnuplotFunction<F,GridView>(f,filename.c_str());
  if (F::Space::continuity >= 0)
    writer.addVertexData(Gnuplotf);
  else
    writer.addCellData(Gnuplotf);
  writer.write(filename.c_str(),GnuplotOptions::ascii);
}

  /**
   * \ingroup IO
   * \brief Writes a set of finite element functions into a single VTK file.
   * @param vars variables in a variable set
   * @param filename name of output file
   * @param options defines the format for the output, either ascii or binary is possible
   */
  template <class VariableSet>
  void writeGnuplotFile(VariableSet const& vars, std::string const& filename, IoOptions options=IoOptions())
  {
    writeGnuplotFile(vars.descriptions.gridView,vars,filename,options);
  }

} // end namespace Kaskade
//---------------------------------------------------------------------

#endif
