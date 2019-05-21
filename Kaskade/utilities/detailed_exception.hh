/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2012-2013 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef DETAILED_EXCEPTION_HH
#define DETAILED_EXCEPTION_HH

#include <string>
#include <stdexcept>
#include <boost/lexical_cast.hpp>

namespace Kaskade{
  /**
   * \ingroup exceptions
   * \brief A wrapper class for conveniently providing exceptions with context information.
   * 
   * This is intended to be the base class for all runtime error exceptions in Kaskade, where runtime errors 
   * are understood as "unexpected situations that prevent the normal program flow and cannot be reliably
   * prevented by careful programming". The prototypical case are I/O errors. Further examples may be the occurence of 
   * inconsistent data, or program failure due to unsuited user input.
   * 
   * Example:
   * \code
   * throw DetailedException("Division by zero not defined",__FILE__,__LINE__);
   * \endcode
   */
  class DetailedException : public std::runtime_error{
  public:
    DetailedException(std::string const& message, std::string const& file, int const line) :
      std::runtime_error(std::string("Exception in ") + file + " at line " + boost::lexical_cast<std::string>(line) + ".\n" + message)
    {
    }
  };
  
  /**
   * \ingroup exceptions
   * \brief An exception class for IO errors.
   */
  class IOException: public DetailedException 
  {
  public:
    IOException(std::string const& message, std::string const& file, int const line) :
      DetailedException(message,file,line) 
    {
    }
  };
  
  /**
   * \ingroup exceptions
   * \brief An exception class for file IO errors.
   */
  class FileIOException: public IOException
  {
  public:
    FileIOException(std::string const& message, std::string const& filename_, std::string const& file, int const line) :
      IOException(std::string("File I/O exception for file ") + filename_ + ".\n" + message,file,line), filename(filename_)
    {
    }
    
    std::string filename;
  };
  
  /**
   * \ingroup exceptions
   * \brief An exception that can be thrown whenever a key lookup fails.
   */
  class LookupException: public DetailedException
  {
  public:
    LookupException (std::string const& message, std::string const& file, int const line):  
      DetailedException(std::string("Key lookup failed:\n")+message,file,line)
    {
    }
  };
  
  /**
   * \ingroup exceptions
   * \brief A base class for linear algebra exceptions.
   */
  class LinearAlgebraException: public DetailedException
  {
  public:
    LinearAlgebraException(std::string const& message, std::string const& file, int const line) :
      DetailedException(std::string("Linear algebra failure:\n")+message,file,line) 
    {
    }
  };
  
  /**
   * \ingroup exceptions
   * \brief To be raised if the matrix is not positive definite.
   */
  class NonpositiveMatrixException: public LinearAlgebraException
  {
  public:
    NonpositiveMatrixException(std::string const& message, std::string const& file, int const line) :
      LinearAlgebraException(std::string("Matrix is not positive definite:\n")+message,file,line) 
    {
    }
  };
  
  /**
   * \ingroup exceptions
   * \brief To be raised if the matrix is singular.
   */
  class SingularMatrixException: public LinearAlgebraException
  {
  public:
    SingularMatrixException(std::string const& message, std::string const& file, int const line):
      LinearAlgebraException(std::string("Matrix is singular.\n")+message,file,line)
    {
    }
  };

  /**
   * \ingroup exceptions
   * \brief To be raised if the matrix is nonsquare.
   */
  class NonsquareMatrixException: public LinearAlgebraException
  {
  public:
    NonsquareMatrixException(std::string const& message, std::string const& file, int const line):
      LinearAlgebraException(std::string("Matrix is nonsquare.\n")+message,file,line)
    {
    }
  };

  /**
   * \ingroup exceptions
   * \brief To be raised if a direct solver fails.
   * 
   * If the solver fails due to a singular matrix to be factorized, use SingularMatrixException instead.
   */
  class DirectSolverException: public LinearAlgebraException
  {
  public:
    DirectSolverException(std::string const& message, std::string const& file, int const line):
      LinearAlgebraException(std::string("Direct solver failed.\n")+message,file,line)
    {
    }
  };
  
}

#endif
