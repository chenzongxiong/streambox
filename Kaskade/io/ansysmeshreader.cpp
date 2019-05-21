/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*  This file is part of the library KASKADE 7                               */
/*    see http://www.zib.de/projects/kaskade7-finite-element-toolbox         */
/*                                                                           */
/*  Copyright (C) 2013-2013 Zuse Institute Berlin                            */
/*                                                                           */
/*  KASKADE 7 is distributed under the terms of the ZIB Academic License.    */
/*    see $KASKADE/academic.txt                                              */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include "io/ansysmeshreader.hh"

namespace Kaskade
{
  void AnsysMeshReader::sortElements(std::vector<std::vector<unsigned int> >& elements) const
  {
    for(std::vector<unsigned int>& element : elements)
    {
      // do not adjust local vertex order for pyramids and simplices
      if( (element[4]==element[5]) && (element[5]==element[6]) && (element[6]==element[7]) ) return;
      
      unsigned int tmp = element[0];
      
      element[0] = element[7];  // 0 <- 7
      element[7] = element[1];  // 7 <- 1
      element[1] = element[3];  // 1 <- 3
      element[3] = element[2];  // 3 <- 2
      element[2] = element[6];  // 2 <- 6
      element[6] = element[5];  // 6 <- 5
      element[5] = tmp;         // 5 <- 0
      // 4 <- 4
    }
  }
  
  void AnsysMeshReader::adjustElement(std::vector<unsigned int>& element, std::vector<std::pair<unsigned int,unsigned int> > const& offsets)
  {
    for(int corner=0; corner<element.size(); ++corner)
    {
      size_t id = offsets.size()-1;
      
      while(offsets[id].first > element[corner] && id > 0) 
        --id;
      element[corner] -= offsets[id].second;
    }
  }
  
  void AnsysMeshReader::readElements(std::ifstream& file, std::string& buffer, std::vector<std::vector<unsigned int> >& elements, std::vector<std::pair<unsigned int,unsigned int> > const& offsets)
  {
    std::vector<unsigned int> element(8);
    char tmp;
    
    getline(file,buffer,delimiter);
    removeLeadingWhiteSpaces(buffer);
    
    while(file.good())
    {
      for(size_t i=0; i<element.size(); ++i) element[i] = readEntry<unsigned int>(file,buffer) - 1;
      
      getline(file,buffer);
      
      adjustElement(element,offsets);
      elements.push_back(element);
      file >> tmp;
      file.unget();
      
      if(tmp=='*')
      {
        getline(file,buffer);
        return;
      }
      
      getline(file,buffer,delimiter);
      removeLeadingWhiteSpaces(buffer);
    }
  }
  
  void AnsysMeshReader::readTetrahedra(std::ifstream& file, std::string& buffer, std::vector<std::vector<unsigned int> >& elements, 
                                       const std::vector<std::pair<unsigned int,unsigned int> >& offsets, short elementSetId)
  {
    std::vector<unsigned int> tetra(4);
    char tmp;
    getline(file,buffer,delimiter);
    removeLeadingWhiteSpaces(buffer);
    
    while(file.good())
    {
      for(size_t i=0; i<4; ++i) tetra[i] = readEntry<unsigned int>(file,buffer) - 1;
      
      getline(file,buffer);
      
      adjustElement(tetra,offsets);
      elements.push_back(tetra);
      elementSetIds.push_back(elementSetId);
      // read char into tmp
      file >> tmp;
      // restore last char of file
      file.unget();
      
      if(tmp == '*')
      {
        getline(file,buffer);
        return;
      }
      
      getline(file,buffer,delimiter);
      removeLeadingWhiteSpaces(buffer);
    }
  }
  
  void AnsysMeshReader::getline(std::ifstream& file, std::string& line, char const delim)
  {
    std::getline(file,line,delim);
    if (delim=='\n') 
      ++currentLine;
    else 
      currentLine += std::count(std::begin(line),std::end(line),'\n');
  }

}