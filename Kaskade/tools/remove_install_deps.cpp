/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2016 Zuse Institute Berlin                                 */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <cstring>
#include <string>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <unistd.h>

//#define KASKADE7 "Kaskade7.3"

int main(int argc, char *argv[])
{
  if (argc!=2)
  {
    std::cout << "usage: remove_install_deps Makefile" << std::endl;
  };
  char tempFileNameTemplate[]="/tmp/remove_usr_dep_out_XXXXXX";
  char tempFileName[80];
  strcpy(tempFileName,tempFileNameTemplate);
  mktemp(tempFileName);
  char cmd[80];
  std::string pwd = (std::string) std::getenv("PWD");
#ifdef DEBUG
  std::cout << "PWD=" << pwd << std::endl;
#endif
  std::string buffer="", compressedBuffer="", dependency, kaskade7Dir;
  const std::string kaskade7=KASKADE7, beginDeps="#DONOTDELETE";
  const int kaskade7Len=kaskade7.length();
  std::streampos inFileBegin, inFileEnd, outFileBegin, outFileEnd;
  std::fstream makefileStream, outStream;
  bool  returnOut=false;
  int lengthKaskDir=0, pos, colonFound, removeCount=0;

  std::ifstream inFile( argv[1], std::ios::binary | std::fstream::in );
  inFileBegin = inFile.tellg();
  inFile.seekg( 0, std::ios::end );
  inFileEnd = inFile.tellg();
  inFile.close();

  makefileStream.open(argv[1],std::fstream::in);
  if ( makefileStream.rdstate() & std::ifstream::failbit )
  {
    std::cout <<  argv[0] << ": couldn't open file " << argv[1] << "." << std::endl;
    exit (10);
  };
  outStream.open(tempFileName,std::fstream::out);
  outFileBegin = outStream.tellg();
  while ( ! makefileStream.eof() && compressedBuffer != beginDeps )
  {
    std::getline(makefileStream,buffer);
    if ( makefileStream.eof() ) 
      break;
    if ( buffer.find('#') != std::string::npos )
    {
      compressedBuffer="";
      for (int i=0;i<buffer.length();i++)
      {
        const char ch=buffer[i];
        if ( (ch!=' ') && (ch!='\t') ) 
         compressedBuffer.append(1,ch);
      };
    };
    outStream << buffer << std::endl;
  };
  
  if ( compressedBuffer != beginDeps )
  {
    std::cout << argv[0] << ": didn't change " << argv[1] << "." << std::endl;
    exit(0);
  };
  
// find last occurrence of "Kaskade7.x" string in current path
  pos = 0;
  while ( true )
  {
    pos = pwd.find(kaskade7,pos);
    if ( pos == std::string::npos )
      break;
    pos = pos+kaskade7Len;
    lengthKaskDir = pos;
  }
  if ( lengthKaskDir == 0 )
  {
    std::cout << kaskade7 << " not found in current directory!" << std::endl <<
      "Check whether the program remove_install_deps was compiled with the correct -DKASKADE7=... flag." <<
      std::endl;
    exit(9);
  }
  kaskade7Dir=pwd.substr(0,lengthKaskDir);

#ifdef DEBUG
  std::cout << "kaskade7Dir=" << kaskade7Dir << std::endl;
#endif
  
  while ( true )
  {
    makefileStream >> dependency;
    if ( makefileStream.eof() )
      break;
    colonFound = dependency.find(':');
    if ( colonFound != std::string::npos )
      outStream << std::endl;
    switch ( (const int) dependency[0] )
    {
      case ':':
      {
        outStream << dependency;
        returnOut = true;
        break;
      };
      case '\\':
      {
        if ( returnOut )
          outStream << " \\" << std::endl;
        returnOut = false;
        break;
      };
      case '/':
      {
        if ( dependency.substr(0,lengthKaskDir) == kaskade7Dir )
        {
          outStream << " $(KASKADE7)" << dependency.substr(lengthKaskDir);
          returnOut = true;
        }
        else
          removeCount++;
        break;
      };
      default:
      {
        if ( colonFound == std::string::npos )
          outStream << " ";
        outStream << dependency;
        returnOut = true;
        break;
      }
    }
  };
  outStream.flush();
  outFileEnd = outStream.tellg();
  makefileStream.close();
  outStream.close();
  strcpy(cmd,"cat ");
  strcat(cmd,tempFileName);
  strcat(cmd," > ");
#ifdef DEBUG
  strcat(cmd,"TestMakefile_removed");
  std::cout << "inFileBeg=" << inFileBegin << ", outFileBeg=" << outFileBegin << std::endl;
  std::cout << "inFileEnd=" << inFileEnd << ", outFileEnd=" << outFileEnd << std::endl;
  std::cout << "inFileSiz=" << inFileEnd-inFileBegin << ", outFileSiz=" <<
    outFileEnd-outFileBegin << std::endl;
#else
  strcat(cmd,argv[1]);
#endif
  system(cmd);
  std::cout << argv[0] << ": number of dependencies removed are " << removeCount <<
    ", filesize reduced by " << (inFileEnd-inFileBegin)-(outFileEnd-outFileBegin) <<
    " bytes." << std::endl;
  std::cout << argv[0] << ": done with " << argv[1] << "." << std::endl;
  unlink(tempFileName);
}
