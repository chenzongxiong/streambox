/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2015 Zuse Institute Berlin                                 */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <iostream>
#include <fstream>
#include <cstring>
#include <algorithm>

bool filter(char* buffer)
{
  if ( std::strstr(buffer,"computing time") || std::strstr(buffer,"#tasks") )
    return true;
  return false;
}

int main(int argc, char *argv[])
{
  if ( argc != 3 )
  {
    std::cout << "usage: filterOutput sourcefile resultfile" << std::endl;
  }
  else
  {
    std::fstream sourceStream, resultStream;
    int const bufferLength=2048;
    char buffer[bufferLength];
    bool inDataToCheck = false;
    int line=0;
    sourceStream.open(argv[1],std::fstream::in);
    if ( ! sourceStream.is_open() )
    {
      std::cout << "*** Fatal error: open of input-stream failed ***" << std::endl;
      exit(10);
    };
    resultStream.open(argv[2],std::fstream::out);
    if ( ! resultStream.is_open() )
    {
      std::cout << "*** Fatal error: open of output-stream failed ***" << std::endl;
      exit(11);
    };
    while ( ! sourceStream.eof() )
    {
      sourceStream.getline(buffer,bufferLength-1);
      line++;
      if ( ! sourceStream.good() && ! sourceStream.eof() )
      {
        std::cout << "*** error: line " << line << " from sourcefile could't be read\n" <<
                     "*** check line for length and possibly increase the bufferLength from\n" <<
                     "*** the current value " << bufferLength << " to a larger value.\n";
        exit(1);
      };
      if ( ! std::strncmp(buffer,"Start ",6) )
        inDataToCheck=true;
      if ( inDataToCheck ) 
        if ( ! filter(buffer) )
          resultStream << buffer << std::endl;
      if ( ! std::strncmp(buffer,"End ",4) || ! std::strncmp(buffer,"Test ",5) )
        inDataToCheck=false;
    }
  }
}
