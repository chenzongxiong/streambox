/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2002-2016 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef AMIRA_HH
#define AMIRA_HH

#define HAVE_AMIRAMESH 1
#undef max
#include <dune/grid/io/file/amirameshwriter.hh>
#include "io/iobase.hh"


/**
 * @file
 * @brief  Output of mesh and solution for visualization software
 * @author Martin Weiser, Bodo Erdmann
 *
 * This file contains functionality needed to write mesh data and
 * solution data for Amira.
 */

//---------------------------------------------------------------------

namespace Kaskade
{
  /**
   * \cond internals
   */
  namespace IoDetail {

    template <class GridView, class Filter, class RAIter>
    struct AddDataTo<GridView,Dune::LeafAmiraMeshWriter<typename GridView::Grid>,Filter,RAIter>
    {
      typedef typename GridView::Grid Grid;
      typedef Dune::LeafAmiraMeshWriter<Grid> AmiraWriter;

      AddDataTo(GridView const& gridView_,AmiraWriter& writer_, Filter const& filter_, RAIter const& names_):
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

        Function const& f = at_c<1>(pair);
        //    std::string const& name = names[id];

        // Get the index set. We hope this is actually the leaf index set
        // (probably yes, but not guaranteed). Unfortunately there seems
        // to be no possibility to check index sets for equality.
        typedef typename Function::Space::IndexSet IndexSet;
        IndexSet const& indexSet = f.space().indexSet();


        //    typedef typename Grid::template Codim<0>::template Partition<Dune::All_Partition>::LeafIterator CellIterator;
        typedef typename GridView::template Codim<0>::template Partition<Dune::All_Partition>::Iterator CellIterator;
        //--> typedef typename Grid::Codim<0>::template Partition<Dune::All_Partition>::LeafIterator CellIterator;
        CellIterator cend = f.space().grid().template leafend<0>();
        // Create evaluator. Todo: Create and use a vertex-based
        // quadrature rule, create and use a shape function cache based on
        // this quadrature rule.
        typename Function::Space::Evaluator eval(f.space());


        if (Function::Space::continuity >= 0) {
          // Continuous function. Written as values at grid vertices.
          Dune::BlockVector<typename Function::ValueType> data(indexSet.size(GridView::dimension));

          // Step through all cells.

          for (CellIterator ci=f.space().grid().template leafbegin<0>(); ci!=cend; ++ci) {
            // Obtain reference element.
            Dune::ReferenceElement<typename Grid::ctype, GridView::dimension> const& refElem = Dune::ReferenceElements<typename Grid::ctype, GridView::dimension>::general(ci->type());
            // Move evaluator to current cell.
            eval.moveTo(*ci);

            // Step through all vertices of the current cell.
            for (int i=0; i<ci->template count<GridView::dimension>(); ++i) {
              // Move evaluator to current vertex.
              eval.evaluateAt(refElem.position(i,GridView::dimension));
              // Store value of function in data vector (inefficient, occurs multiple times...).
              data[indexSet.subIndex(*ci,i,GridView::dimension)] = f.value(eval);
            }
          }

          writer.addVertexData(data,gridView);
        } else {
          // Discontinuous function. Written as values at cells.
          Dune::BlockVector<typename Function::ValueType> data(indexSet.size(0));

          assert("implement me!"==0);

          writer.addVertexData(data,gridView);
        }
      }

    private:
      GridView const& gridView;
      AmiraWriter&  writer;
      Filter const& filter;
      RAIter        names;

    };

  } // End of namespace IoDetail

  /**
   * \endcond
   */

  //---------------------------------------------------------------------

  /**
   *   This procedure writes data in AmiraMesh ascii/binary format, e.g. for visualization in AMIRA software.
   *   @param gridView gridView, i.e. LeafGridView
   *   @param description variable set description
   *   @param vars data of variables
   *   @param filename name of output file
   *   @param options defines the format for the output,e.g., ascii/binary
   */
  template <class GridView, class VariableSet>
  void writeAMIRAFile(GridView const& gridView,
      VariableSet const& vars,
      std::string const& filename,
      IoOptions options=ioOptions_default)
  {
    typedef typename GridView::Grid Grid;
    typedef Dune::LeafAmiraMeshWriter<Grid> AmiraWriter;
    AmiraWriter writer(gridView.grid());

    char fullname[128];
    sprintf(fullname,"%s.am",filename.c_str());

    if (options.outputType == IoOptions::ascii)
      writePartialFile(gridView,writer,vars,fullname,UnaryTrue(),1);
    else if (options.outputType == IoOptions::binary)
      writePartialFile(gridView,writer,vars,fullname,UnaryTrue(),0);
    else std::cout << "Amira:  graphic output for nonvalid driver" << std::endl;
  }

  /*
template <class VariableSetDescription>
void writeAMIRAFile(VariableSetDescription const& description,
                    typename VariableSetDescription::VariableSet const& vars,
                    std::string const& filename) 
{
    typedef Dune::LeafAmiraMeshWriter<typename VariableSetDescription::Grid> AmiraWriter;
    AmiraWriter writer(description.gridManager.grid());

  	char fullname[128];
  	sprintf(fullname,"%s.am",filename.c_str());

	IoOptions options;
  	if (options.outputType == IoOptions::ascii)
		writePartialFile(writer,description,vars,fullname,UnaryTrue(),1);
  	else if (options.outputType == IoOptions::binary)
		writePartialFile(writer,description,vars,fullname,UnaryTrue(),0);
  	else std::cout << "Amira:  graphic output for nonvalid driver" << std::endl;
}
   */
} // end namespace Kaskade

#endif
