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

#ifndef PARTIALINDEXSET_HH
#define PARTIALINDEXSET_HH

#include <cassert>
#include <vector>

#include <boost/mpl/range_c.hpp>

#include <boost/fusion/algorithm.hpp>

#include "dune/grid/common/indexidset.hh"
#include "dune/grid/common/referenceelements.hh"

#include "fem/fetransfer.hh"

/**
 *\todo recompute internal representation whenever grid changes.
 *
 * Part has to provide a method bool contains(typename Grid::Traits::template Codim<0>::Entity const&) const;
 */
template<class Grid, class IndexSet, class Part>
class PartialIndexSet {
  typedef PartialIndexSet<Grid,IndexSet,Part> Self;
  typedef typename Grid::Traits::template Codim<0>::Entity Cell;
  
public:
  template <int cd>
  struct Codim
  {
    template <Dune::PartitionIteratorType pitype>
    struct Partition
    {
      class Iterator
      {
        typedef typename IndexSet::template Codim<cd>::template Partition<pitype>::Iterator Iter;
        
      public:
        typedef typename Grid::template Codim<0>::Entity Entity;
        
        Iterator(Iter const& cur_, Iter const& end_, Part const& part_):
          cur(cur_), end(end_), part(part_) {
          ahead();
        }

        Iterator& operator++() 
        {
          if (cur!=end)
            ++cur;
          ahead();
          return *this;
        }

        Entity& operator*() const { return *cur; }
        Entity* operator->() const { return &*cur; }

        bool operator==(Iterator const& i) const { return cur==i.cur; }
        bool operator!=(Iterator const& i) const { return !(*this == i); }
        
        operator typename Grid::Traits::template Codim<0>::EntityPointer() const 
        {
          return cur;
        }
        
      private:
        void ahead() 
        {
          while (cur!=end && !part.contains(*cur))
            ++cur;
        }
            
        Iter cur, end;
        Part const& part;
      };
    };
  };

  PartialIndexSet(GridSignals& signals, IndexSet const& indexSet_, Part const& part_):
    indexSet(indexSet_), part(part_)
  {
    typedef typename IndexSet::template Codim<0>::
        template Partition<Dune::All_Partition>::Iterator CellIterator;
    CellIterator end = indexSet.template end<0,Dune::All_Partition>();
    for (CellIterator ci=indexSet.template begin<0,Dune::All_Partition>(); ci!=end; ++ci) 
      if (part.contains(*ci)) {
        insert(ci->type(),indexSet.index(*ci));
        boost::fusion::for_each(boost::mpl::range_c<int,1,Grid::dimension+1>(),
                                Update(*this,*ci));
      }
  }
  
      
  template <int cc>
  int index (typename Grid::Traits::template Codim<cc>::Entity const& e) const 
  {
    typename Map::const_iterator i = idx.find(e.type());
    if (i==idx.end()) return -1;
    else              return i->second.first[indexSet.index(e)];
  }

  template<class EntityType>
  int index (const EntityType& e) const
  {
    return index<EntityType::codimension>(e);
  }

  template<int cc>
  int subIndex (Cell const& e, int i) const
  {
    // unfortunately the following is not implemented in every grid interface...
    // return index(*e.template entity<cc>(i));

    Dune::GeometryType subentityType =
      Dune::ReferenceElements<typename Grid::ctype,Grid::dimension>::general(e.type()).type(i,cc);
    typename Map::const_iterator it = idx.find(subentityType);
    if (it==idx.end()) return -1;
    else               return it->second.first[indexSet.template subIndex<cc>(e,i)];
  }
  
  const std::vector<Dune::GeometryType>& geomTypes (int codim) const
  {
    return indexSet.geomTypes(codim);
  }

  int size (Dune::GeometryType type) const
  {
    typename Map::const_iterator i= idx.find(type);
    if (i==idx.end()) return 0;
    else              return i->second.second;
  }

  int size (int codim) const
  {
    int count = 0;
    for (typename Map::const_iterator i=idx.begin(); i!=idx.end(); ++i)
      if (i->first.dim() == Grid::dimension-codim)
        count += i->second.second;
    return count;
  }

  template<class EntityType>
  bool contains (const EntityType& e) const
  {
    return indexSet.contains(e) && index(e)>=0;
  }

  template<int cd, Dune::PartitionIteratorType pitype>
  typename Codim<cd>::template Partition<pitype>::Iterator begin () const
  {
    return typename Codim<cd>::template Partition<pitype>::Iterator(indexSet.template begin<cd,pitype>(),
                                                                    indexSet.template end<cd,pitype>(),
                                                                    part);
  }
  
  template<int cd, Dune::PartitionIteratorType pitype>
  typename Codim<cd>::template Partition<pitype>::Iterator end () const
  {
    return typename Codim<cd>::template Partition<pitype>::Iterator(indexSet.template end<cd,pitype>(),
                                                                    indexSet.template end<cd,pitype>(),
                                                                    part);
  }

private:
  IndexSet const& indexSet;
  Part const& part;

  // For each geometry type, separate contiguous indices are
  // maintained, which are stored according to the underlying index
  // set. A negative value denotes entities which are not contained in
  // the partial index set. The second part in the data type is a
  // running counter used for generating contiguous indices.
  typedef std::map<Dune::GeometryType,std::pair<std::vector<int>,int> > Map;
  Map idx;

  void insert(Dune::GeometryType gt, int index) 
  {
    std::pair<std::vector<int>,int>& gtIdx = idx[gt];

    if (gtIdx.first.empty()) {
      gtIdx.first.resize(indexSet.size(gt),-1);
      gtIdx.second = 0;
    }

    if (gtIdx.first[index]<0) 
      gtIdx.first[index] = gtIdx.second++;
  }
  
  // A MPL/fusion functor which inserts all the cell's subentities of
  // given codimension into the index map.
  struct Update 
  {
    Update(Self& pis_, Cell const& cell_): pis(pis_), cell(cell_) {}
    
    template <class Integer>
    void operator()(Integer const /* codim */) const 
    {
      int const codim = Integer::value;
      int const count = cell.template count<codim>();
      

      for (int i=0; i<count; ++i)
        pis.insert(Dune::ReferenceElements<typename Grid::ctype,Grid::dimension>::general(cell.type()).type(i,codim),
                   pis.indexSet.template subIndex<codim>(cell,i));
    }

  private:
    Self& pis;
    Cell const& cell;
  };
};


#endif
