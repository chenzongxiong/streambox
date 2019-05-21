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

#ifndef FACTORY_2_HH
#define FACTORY_2_HH

#include <map>
#include <memory>
#include <string>

#include <boost/lexical_cast.hpp>

#include "utilities/functor.hh"

namespace Kaskade
{
  /**
   * \brief A pluggable factory.
   *
   * Allows to create objects of type AbstractProduct. The concrete object is determined by
   * an identifier that is associated with a creator functor or function pointer. If the functor
   * or function pointer take arguments specify their types with the variadic template parameter
   * Arguments.
   *
   */
  template
  <
    class Identifier,
    class AbstractProduct,
    typename... Arguments
  >
  class Factory2
  {
    /// Exception for the case that an entry is not found.
    class Exception : public std::exception
    {
    public:
      explicit Exception(Identifier const& id_) : id(id_){}

      const char* what()
      {
        std::string message("No entry with id \"");
        message += boost::lexical_cast<std::string>(id);
        message += "\" found in factory.\n";
        return message.c_str();
      }

    private:
      Identifier const& id;
    };

  public:
    // Store creator in functor (creation via functor and function pointer possible)
    typedef Functor<AbstractProduct*,Arguments...> ProductCreator;

  private:
    // Associative container
    typedef std::map<Identifier,ProductCreator> Map;

  public:
    /// Add entry to factory.
    /**
     * \param id unique identifier associated with creator
     * \param creator functor or function pointer realizing the creation of AbstractProduct
     *
     * \param false if id already exists (in this case creator is not added to the factory)
     */
    bool add(Identifier const& id, ProductCreator creator)
    {
      return map.insert(typename Map::value_type(id,creator)).second;
    }

    /// Remove entry from factory.
    /**
     * \param id unique identifier associated with the entry to be deleted.
     * \return true if entry was deleted, false if no entry corresponding to id has been found
     */
    bool remove(Identifier const& id)
    {
      return map.erase(id)==1;
    }


    std::unique_ptr<AbstractProduct> create(Identifier const& id, Arguments... args)
    {
      typename Map::iterator i = map.find(id);
      if(i != map.end()) return std::unique_ptr<AbstractProduct>(i->second(args...));
      else throw Exception(id);

      return std::unique_ptr<AbstractProduct>(nullptr);
    }

  private:
    Map map;
  };

}  // namespace Kaskade
#endif
