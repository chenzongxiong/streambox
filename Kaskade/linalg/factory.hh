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

#ifndef FACTORY_HH
#define FACTORY_HH

#include <cassert>
#include <map>
#include <memory>
#include <string>

#include "utilities/detailed_exception.hh"

namespace Kaskade
{
  /**
   * \ingroup utilities
   * \brief Abstract base class of creators that can be registered at the pluggable factory.
   */
  template <class Product>
  struct Creator 
  {
    typedef int Argument;
    
    std::unique_ptr<Product> operator()(Argument /* dummy */) const
    {
      return new Product;
    }
    
    template <class Key>
    static std::string const& lookupFailureHint(Key const& /* ignored */)
    {
      static std::string msg = "Key not registered.";
      return msg;
    }
  };
  
  
  /**
   * \ingroup utilities
   * \brief A pluggable factory.
   * 
   * \tparam K the key type the values of which specify which product is to be created
   * \tparam P the product type
   * 
   * \see Creator
   */
  template <class K, class P>
  class Factory 
  {
    typedef std::map<K,Creator<P> const*> Map;
    
  public:
    typedef K Key;
    typedef P Product;
    
    
    static std::unique_ptr<Product> create(Key key, typename Creator<Product>::Argument arg)
    {
      return makeProduct(map().find(key),Creator<Product>::lookupFailureHint(key),arg);
    }
    
    static std::unique_ptr<Product> createAny(typename Creator<Product>::Argument arg)
    {
      return makeProduct(begin(map()),"no key registered",arg);
    }
    
    
    
    /**
     * \brief Registers a creator for creation of Product objects with the
     * given \arg key. 
     * 
     * Note that the creator object has to exist for the
     * whole lifetime of the program.
     */
    static void plugin(Key key, Creator<Product> const* creator) 
    {
      map()[key] = creator;
    }
    
    
  private:
    
    static std::unique_ptr<Product> makeProduct(typename Map::const_iterator i, std::string message, typename Creator<Product>::Argument arg)
    {
      if (i==map().end())
        throw(LookupException(message,__FILE__,__LINE__));
      
      return i->second->create(arg);
    }
    
    static Map& map()  
    {
      static Map m;     // this realizes a singleton
      return m;
    }
  };
  
}  // namespace Kaskade
#endif
