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

#ifndef KASKADE_SAVE_FILE_HH_
#define KASKADE_SAVE_FILE_HH_

#include <string>

namespace Kaskade
{
  /// create file name
  /**
   * Checks if a file called <fileName><ending> exists. If this is the case tries fileName0...0ending - fileName9...9ending,
   *  where the number of added digits can be controlled via variable length
   * \param desiredName desired file name
   * \param ending (i.e. '.doc', '.cpp', '.vtu', ...)
   * \param useEnding true: create file name with ending, false: use ending only for checking file name candidates
   * \param length length of the extension if file name exists
   */
  std::string createFileName(std::string const& desiredName, std::string const& ending = std::string(), bool useEnding = true, int length = 3);
  
  bool createFolder(std::string folder_name);
  
  std::string createFileNameInFolder(std::string folderName, std::string const& desiredName, std::string const& ending = std::string(), bool const useEnding = true, int const length = 3);
}

#endif /* KASKADE_SAVE_FILE_HH_ */
