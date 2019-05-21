#ifndef GET_FROM_TYPE_OR_FUNCTOR_HH_
#define GET_FROM_TYPE_OR_FUNCTOR_HH_

#include "utilities/detailed_exception.hh"

// forward declaration
namespace Dune
{
  template <class,int,int> class FieldMatrix;
}

namespace Kaskade
{
  /// Get object from reference or functor. Is itself a functor.
  /**
   * Source should not be inherited from Type to guarantee proper behaviour!!
   *
   * \param Type object type
   * \param Source source, i.e. functor type if object is obtained via a functor.
   */
  template <class Type, class Source=Type> struct Get;

  template <class Type, class Source>
  struct Get{
    Get(){ source_ = nullptr; }
    explicit Get(Source const& source) : source_(&source){}
    Get(Get const& get) : source_(get.source_){}

    Get& operator=(Get const& get){ source_ = get.source_; return *this; }

    Type operator()() const
    {
      if(source_ != nullptr) return (*source_)();
      else throw DetailedException(std::string("Source has not been initialized!"),__FILE__,__LINE__);
    }

  private:
    Source const* source_;
  };

  template <class Type>
  struct Get<Type,Type>{
    Get(){ type_ = nullptr; }
    explicit Get(Type const& type) : type_(&type){}
    Get(Get const& get) : type_(get.type_){}

    Get& operator=(Get const& get){ type_ = get.type_; return *this; }

    Type const& operator()() const
    {
      if(type_ != nullptr) return *type_;
      else throw DetailedException(std::string("Type has not been initialized!"),__FILE__,__LINE__);
    }

  private:
    Type const* type_;
  };

  template <class Scalar, int dim, class Type>
  struct GetSubType{
    typedef typename Type::SubType type;
  };

  template <class Scalar, int dim>
  struct GetSubType<Scalar,dim,Dune::FieldMatrix<Scalar,dim,dim> >{
    typedef Dune::FieldMatrix<Scalar,dim-1,dim-1> type;
  };

} // end of namespace Kaskade

#endif
