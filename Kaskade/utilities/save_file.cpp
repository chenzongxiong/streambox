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

#include <fstream>
#include <iostream>

#include <sys/stat.h>

#include <boost/lexical_cast.hpp>

#include "io/iobase.hh"
#include "utilities/date_and_time.hh"
#include "utilities/save_file.hh"

namespace Kaskade
{
  std::string createFileName(std::string const& desiredName, std::string const& ending, bool useEnding, int length)
  {
    std::ifstream infile;
    std::string fileName = desiredName + "_" + paddedString(0,length) + ending; // file name does not yet exist

    // open and close file -> set failbit
    infile.open(fileName.c_str());
    infile.close();

    size_t index = 1;
    size_t const maxIndex = (int)pow(10,length);
    while(!infile.fail())
      {
        // reset streams internal failbit state
        infile.clear(std::ios::failbit);

        // change filename
        fileName = desiredName + "_" + paddedString(index, length) + ending;

        // open and close file -> set failbit
        infile.open(fileName.c_str());
        infile.close();

        ++index;
        if(index==maxIndex)
          {
            std::cerr << __FILE__ << ": " << __LINE__ << ": Tried " << maxIndex << " names, all of them are currently used." << std::endl;
          }
      }
    if(!useEnding) fileName = desiredName + "_" + paddedString(index-1, length);

    return fileName;
  }

  bool createFolder(std::string folder_name)
  {
    return ( mkdir(folder_name.c_str(), 0777) == 0 );
  }

  std::string createFileNameInFolder(std::string folderName, std::string const& desiredName, std::string const& ending, bool useEnding, int length)
  {
    if(folderName.empty()) folderName = desiredName;

    createFolder( appendDate(folderName) );

    folderName += "/" + desiredName;

    return createFileName(folderName, ending, useEnding, length);
  }
}
