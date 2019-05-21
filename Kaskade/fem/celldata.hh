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

#ifndef CELLDATA_HH
#define CELLDATA_HH

#include "fem/lagrangespace.hh"

namespace Kaskade
{
  /** \ingroup adapt
   * \brief Class that stores information for each cell of a grid
   *
   * In many cases it is useful to store information on a cell in a vector.
   * An example is the storage of error indicators, or some debugging information,
   * such as which entity has been marked for refinement, or has been refined.
   * The class CellData provides infrastructure for these tasks with an emphasis
   * on error estimation.
   *
   * Since a CellData object contains many EntityPointers, it may become very large.
   * If memory space is an urgent issue, it may be a good idea to replace the current
   * data structure by a different one. In any case it might be a good idea to restrict
   * the livetime of these objects to when they are actually needed.
   */

  template<class Grd, class T=Dune::FieldVector<double,1> >
  class CellData
  {
  public:
    typedef Grd Grid;
    typedef T ValueType;
    typedef typename Grid::template Codim<0>::Entity Cell;
    /// A pair, where the data and an EntityPointer is stored
    typedef std::pair<T, typename Cell::EntityPointer>  CellDataPair;
    typedef std::vector<CellDataPair> CellDataVector;
    static int const dim = Grid::dimension;

  private:

    struct GreaterFirst
    {
      bool operator() (CellDataPair const & p1, CellDataPair const & p2) { return fabs(p1.first) > fabs(p2.first);}
    };

  public:
    /// Default constructor
    CellData() : data(), sorted(false), maxed(false), corrEst(false) {};

    /// Construct Cell Data from a given CellDataVector
    explicit CellData(CellDataVector const& data_) : data(data_), sorted(false),  maxed(false) {}

    ///Construct Cell Data from a given CellDataVector of different scalar type
    template<class S>
    explicit CellData(std::vector< std::pair<S, typename Cell::EntityPointer> > const& data_) : sorted(false),  maxed(false)
    {
      for(int i=0; i<data_.size(); ++i){
        data.push_back(std::make_pair(data_[i].first,data_[i].second));
      }
    }

    template<class Space, class Vector>
    explicit CellData(Vector const& data_,Space const& ds) : data(), sorted(false), maxed(false)
    {
      Space dspace(ds.gridManager(), ds.gridView(), 0);
      typename Space::Evaluator isfs(dspace);
      typedef typename Space::GridView::template Codim<0>::Iterator CellIterator;
      T fv;
      for (CellIterator ci=dspace.gridView().template begin<0>();
          ci!=dspace.gridView(). template end<0>(); ++ci)
      {
        isfs.moveTo(*ci);
        fv[0]=T(data_[isfs.globalIndices()[0]]);
        data.push_back(std::make_pair(fv,ci));
      }
    }

    /// Assignment operator from a CellData object
    CellData<Grid,T>& operator=(const CellData<Grid,T>& ei)
    {
      data = ei.data;
      sorted=ei.sorted;
      maxed=false;
      return *this;
    }

    /// Assignment operator between two CellDataVector
    CellData<Grid,T>& operator=(const CellDataVector& ei)
    {
      data=ei;
      sorted=false;
      maxed=false;
      return *this;
    }

    /// Compute the sum over all entries
    double sum() const
    {
      double totalErrorSq=0.0;
      for(int i=0; i<data.size(); ++i)
        totalErrorSq+=data[i].first;
      return totalErrorSq;
    }

    /// Compute the sum over the absolute values of all entries
    double abssum() const
    {
      double totalErrorSq=0.0;
      for(int i=0; i<data.size(); ++i)
        totalErrorSq+=fabs(data[i].first);
      return totalErrorSq;
    }

    /// Sort all entries from large to small
    void sort()
    {
      if(!sorted) std::sort(data.begin(),data.end(),GreaterFirst());
      sorted=true;
    }

    /// Compute the maximal value of all entries
    double maxElement() const
    {
      if(maxed) return maxEle;
      maxed=true;
      if(sorted)
      {
        maxEle=fabs(data[0].first);
        return maxEle;
      }
      maxEle=-1e300;
      for(int i=0; i<data.size(); i++)
        maxEle = std::max(maxEle,fabs(data[i].first));
      return maxEle;
    }

    /// Scale all entries, such that their maximal value is 1
    void normalize()
    {
      double m=maxElement();
      for(int i=0; i<data.size(); i++) data[i].first *= 1.0/m;
      maxEle=1.0;
    }

    // Routines for reading the data. Designed to allow flexibility in the choice of the container

    /// Return a const_iterator on data.begin()
    typename CellDataVector::const_iterator begin() const {return data.begin();}

    /// Return a const_iterator on data.end()
    typename CellDataVector::const_iterator end() const {return data.end();}

    /// Return a const_iterator on data.begin()
    typename CellDataVector::iterator begin() {return data.begin();}

    /// Return a const_iterator on data.end()
    typename CellDataVector::iterator end() {return data.end();}

    /// Create a function space element from *this
    template<typename FSElement>
    void toFunction(FSElement& fse) const
    {
      typedef typename FSElement::Space ImageSpace;

      Dune::Matrix<typename FSElement::StorageValueType> localCoefficients;
      std::vector<typename Grid::ctype> volume(fse.space().degreesOfFreedom(),0.0);
      *fse = typename ImageSpace::Mapper::ShapeFunctionSet::Scalar(0.0);

      typename ImageSpace::Evaluator isfs(fse.space());
      for(int i=0; i<data.size(); i++)
      {
        typename CellDataPair::second_type ci(data[i].second);
        isfs.moveTo(*ci);
        typename Grid::ctype myvolume(ci->geometry().volume());
        for (int j=0; j<isfs.size(); ++j) volume[isfs.globalIndices()[j]]+=myvolume;
        localCoefficients.setSize(isfs.size(),1);
        for (int j=0; j<isfs.size(); ++j)
        {
          localCoefficients[j][0]=data[i].first;
          localCoefficients[j][0] *= myvolume;
          (*fse)[isfs.globalIndices()[j]]+= localCoefficients[j][0];
        }
      }
      for(int i=0; i<volume.size(); ++i) (*fse)[i]/=volume[i];
    }

    void setCorrectedEstimate(double est)
    {
      corrEst=true;
      correctedestimate=est;
    }

    double getCorrectedEstimate() const
    {
      return correctedestimate;
    }

    bool estIsCorrected() const
    {
      return corrEst;
    }


  private:
    CellDataVector data;
    mutable bool sorted,maxed;
    mutable T maxEle;
    bool corrEst;
    double correctedestimate;
  };


  /// Create a CellDataVector, in which the largest entities of *this are marked by 1, using the "bulk criterion"
  /** To be used together with GridManager::mark(...)*/
  template<class Grid, class T>
  typename CellData<Grid,int>::CellDataVector markByBulkCriterion(CellData<Grid,T>& ic, double Theta)
  {
    double partialErrorSq(0.0), pSum(0.0);
    ic.sort();
    double totalErrorSq=ic.abssum();
    Theta*=Theta;
    typename CellData<Grid,int>::CellDataVector iv;
    typename CellData<Grid,T>::CellDataVector::const_iterator iind;
    for(iind=ic.begin(); iind != ic.end(); iind++)
    {
      typename CellData<Grid,int>::CellDataPair markpair(0,iind->second);
      if(partialErrorSq <= Theta*totalErrorSq)
      {
        markpair.first=1;
        partialErrorSq += fabs(iind->first);
        pSum += iind->first;
      }
      iv.insert(iv.end(),markpair);
    }
    //std::cout << "Expected change:" << sqrt(fabs(pSum)) << " relative:" << sqrt(fabs(1-(ic.sum()-(1-0.25*0.25)*pSum)/ic.sum())) << std::endl;
    return iv;
  }

  /// Delete all error entries of entities which are above a certain level
  /** This may be sensible, if the geometry resolution is coarse, such that too deep refinement is not useful*/
  template<class Grid, class T>
  void deleteAllAboveLevel(CellData<Grid,T>& ic, int maxlevel)
  {
    typename CellData<Grid,T>::CellDataVector::iterator iind;
    for(iind=ic.begin(); iind != ic.end(); iind++)
    {
      if(iind->second->level() > maxlevel) iind->first = 0.0;
    }
  }

  /// Create a CellDataVector, in which the largest entities of *this are marked by 1, using the "max criterion"
  // /** To be used together with GridManager::mark(...)*/
  template<class Grid, class T>
  typename CellData<Grid,int>::CellDataVector markByMaxCriterion(CellData<Grid,T>& ic, double Ratio)
  {
    double m=ic.maxElement();
    Ratio*=Ratio;
    typename CellData<Grid,int>::CellDataVector iv;
    typename CellData<Grid,T>::CellDataVector::const_iterator iind;
    for(iind=ic.begin(); iind!=ic.end();++iind)
    {
      typename CellData<Grid,int>::CellDataPair markpair(0,iind->second);
      if(fabs(iind->first) >= Ratio*m) markpair.first=1;
      iv.insert(iv.end(),markpair);
    }
    return iv;
  }

  /// Create a CellDataVector, in which the largest entities of *this are marked by 1, using the "max criterion"
  // /** To be used together with GridManager::mark(...)*/
  template<class Grid, class T>
  typename CellData<Grid,int>::CellDataVector markWithinInterval(CellData<Grid,T>& ic, double lower, double upper, int times)
  {
    typename CellData<Grid,int>::CellDataVector iv;
    typename CellData<Grid,T>::CellDataVector::const_iterator iind;
    for(iind=ic.begin(); iind!=ic.end();++iind)
    {
      typename CellData<Grid,int>::CellDataPair markpair(0,iind->second);
      if(iind->first>=lower && iind->first <=upper) markpair.first=times;
      iv.insert(iv.end(),markpair);
    }
    return iv;
  }

  /// Create a CellDataVector, in which the largest entities of *this are marked by 1, using the "max criterion"
  // /** To be used together with GridManager::mark(...)*/
  template<class Grid>
  typename CellData<Grid,int>::CellDataVector unmarkOutOf(CellData<Grid,int> const& ic, Dune::FieldVector<typename Grid::ctype,Grid::dimension> center, double radius)
  {
    typename CellData<Grid,int>::CellDataVector iv;
    typename CellData<Grid,int>::CellDataVector::const_iterator iind;
    for(iind=ic.begin(); iind!=ic.end();++iind)
    {
      int isWithin=0;
      typename CellData<Grid,int>::CellDataPair markpair(*iind);
      Dune::FieldVector<typename Grid::ctype,Grid::dimension> difference=markpair.second->geometry().center()-center;
      if(difference.two_norm() <= radius) isWithin = 1;
      markpair.first *= isWithin;
      iv.insert(iv.end(),markpair);
    }
    return iv;
  }

  /// Scale error indicators for H^1, such that L_2 indicators are the result
  /** This is done by the heuristics that the L_2 error is h * H^1-error.
   */
  template<class Indicator>
  Indicator errorL2(Indicator const& indi)
  {
    const int dim(Indicator::dim);
    typename Indicator::CellDataVector iv;

    typename CellData<typename Indicator::Grid>::CellDataVector::const_iterator iind;
    for(iind=indi.begin(); iind != indi.end(); iind++)
    {
      typename CellData<typename Indicator::Grid>::CellDataPair
      markpair(iind->first[0]*std::pow(iind->second->geometry().volume(),2.0/dim),iind->second);
      iv.insert(iv.end(),markpair);
    }
    return Indicator(iv);
  }


  template<typename FSElement>
  void extendMarks(FSElement& fse,  FSElement const& fu, int neighboursForMarking)
  {
    typedef typename FSElement::Space ImageSpace;
    typedef typename ImageSpace::Grid Grid;
    typedef typename Grid::template Codim<0>::Entity Entity;

    Dune::Matrix<typename FSElement::StorageValueType > localCoefficients;
    *fse = typename ImageSpace::Mapper::ShapeFunctionSet::Scalar(0.0);

    typedef typename ImageSpace::IndexSet::template Codim<0>::template Partition<Dune::All_Partition>::Iterator CellIterator;
    typename ImageSpace::Evaluator isfs(fse.space());
    typename FSElement::Space::Evaluator dsfs(fu.space());
    for (CellIterator ci=fse.space().indexSet().template begin<0,Dune::All_Partition>();
        ci!=fse.space().indexSet().template end<0,Dune::All_Partition>(); ++ci)
    {
      isfs.moveTo(*ci);
      dsfs.moveTo(*ci);

      std::vector<Dune::FieldVector<typename FSElement::Space::Grid::ctype, FSElement::Space::dim> >const &
      iNodes(isfs.shapeFunctions().interpolationNodes());

      localCoefficients.setSize(iNodes.size(),1);
      dsfs.evaluateAt(iNodes[0]);
      if(fu.value(dsfs) > 0.5) 
        (*fse)[isfs.globalIndices()[0]]= neighboursForMarking;
      typename Entity::LeafIntersectionIterator faceEnd = ci->ileafend();
      for (typename Entity::LeafIntersectionIterator face=ci->ileafbegin(); face!=faceEnd; ++face) {
        if (!face.neighbor())
          continue; // interior face
        dsfs.moveTo(*face.outside());
        dsfs.evaluateAt(iNodes[0]);
        if(fu.value(dsfs) > 0.5)
          (*fse)[isfs.globalIndices()[0]]+=1.0;
      }
    }
    for(int i=0; i<(*fse).size(); ++i)
      if((*fse)[i] >= neighboursForMarking) (*fse)[i]=1.0;

  }

  namespace Bridge
  {

    template<class ErrorEst, class GridMan>
    ErrorEst extendMarkings(ErrorEst const& orig, GridMan & gridMan)
    {
      typedef FEFunctionSpace<DiscontinuousLagrangeMapper<double,typename GridMan::Grid> > PWCSpace;
      typedef typename PWCSpace::template Element<1>::type FSE;

      std::cout << "Before Extension:" << orig.getData().sum() << std::endl;

      PWCSpace pwc(gridMan,gridMan.grid().leafIndexSet(),0);
      FSE marks(pwc),smark(pwc);
      std::vector<double> vmark;
      {
        orig.getData().toFunction(marks);
        extendMarks(smark,marks,2);
        for(int i=0; i<(*smark).size();++i)
          vmark.push_back((*smark)[i][0]);
      }
      ErrorEst extmarks(CellData<typename GridMan::Grid, Dune::FieldVector<double, 1> >(vmark,pwc));
      std::cout << "After Extension:" << extmarks.getData().sum() << std::endl;
      return extmarks;
    }
  }

/*  class AbstractErrorEstimate;

  template<class Grid, class T=Dune::FieldVector<double,1> >
  class CellDataErrorEstimate : public AbstractErrorEstimate
  {
  public:
    static int const dim = CellData<Grid,T>::dim;

    CellDataErrorEstimate(CellData<Grid,T> cd_) : cd(cd_) {};

    virtual double absoluteError() const {
      if(cd.estIsCorrected())
        return cd.getCorrectedEstimate();
      else
        return std::sqrt(fabs(cd.abssum())); }

    CellData<Grid,T>& getData() { return cd; }
    CellData<Grid,T> const& getData() const { return cd; }
  private:
    CellData<Grid,T> cd;
  };*/
} /* end of namespace Kaskade */

#endif
