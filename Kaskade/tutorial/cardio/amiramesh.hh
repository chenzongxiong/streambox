#ifndef AMIRAMESH_HH
#define AMIRAMESH_HH

#include <stdio.h>

#include <vector>
#include <map>

#include <amiramesh/AmiraMesh.h>

extern int ComputeBoundary(int nPoints, int nTet, int *tets, int *bounds);
extern void SetMaterialsIds(int nPoints, int *newIds, int nUsedPoints,
  int nBndTr, int *bndTrI, unsigned char *bndIds,
  int nUsedBndTr, int *usedBndTrI, unsigned char *usedBndIds);

template <int dim, int worldDim, class Grid>
class DuneAmiraMesh {

public:
  typedef double DoubleCoord[worldDim];
  typedef float FloatCoord[worldDim];
  typedef int TetrahedraIndices[4];
  typedef int TriangleIndices[3];

private:
  AmiraMesh* mesh;
  HxParamBundle* materials, *boundaryIds;

  int nPoints, nTet, nBndTr, parSize, nMaterials, nBoundaryIds;
  FloatCoord *floatCoords;
  DoubleCoord *doubleCoords;
  TetrahedraIndices *tetI;
  TriangleIndices *bndTrI;
  unsigned char *matIds;
  unsigned char *bndIds;
  unsigned char *tetBndIds;
  int *newIds;
  int *renumberedUGIds, *revRenumberedUGIds;

  bool selected[256];
  unsigned int cnts[256];
  int nUsedPoints, nUsedTet, nUsedBndTr;
  FloatCoord *usedFloatCoords;
  TetrahedraIndices *usedTetI;
  TetrahedraIndices *usedBounds;
  TriangleIndices *usedBndTrI;
  unsigned char *usedMatIds;
  unsigned char *usedBndIds;
  unsigned char *usedBndCharId;
  std::map<std::string,int> matsMap;
  std::map<std::string,int> bIdsMap;

public:
  DuneAmiraMesh(const char *fileName)
    : floatCoords(0), doubleCoords(0), tetI(0), bndTrI(0), matIds(0), bndIds(0), tetBndIds(0),
      newIds(0), renumberedUGIds(0), revRenumberedUGIds(0),
      nUsedPoints(0), nUsedTet(0), nUsedBndTr(0),
      usedFloatCoords(0), usedTetI(0), usedBounds(0), usedBndTrI(0)
    {
      mesh = AmiraMesh::read(fileName);
      if (mesh==0)
        {
          std::cerr << "Failed to read " << fileName << '\n';
          throw -34101;
        }
      nPoints = mesh->nElements("Vertices");
	  if (nPoints <= 0 )
	    nPoints = mesh->nElements("Nodes");
	  nTet = mesh->nElements("Tetrahedra");
	  nBndTr = mesh->nElements("BoundaryTriangles");
	  for (int k = 0; k<256; k++)
	    {
	      selected[k] = false;
	      cnts[k] = 0;
	    }
      parSize = mesh->parameters.size();
      materials = mesh->parameters.materials();
      nMaterials = materials->nBundles();
      boundaryIds = mesh->parameters.boundaryIds();
      nBoundaryIds = boundaryIds->nBundles();

      int k, id;
      const char *str;

      for (k=0; k<nMaterials; k++)
        {
          materials->bundle(k)->findNum("Id", id);
          str = materials->bundle(k)->name();
          matsMap[str] = id;
        }
      for (k=0; k<nBoundaryIds; k++)
        {
          boundaryIds->bundle(k)->findNum("Id", id);
          boundaryIds->bundle(k)->findString("Info", str);
          bIdsMap[str] = id;
        }
      newIds = new int[nPoints];
    };

  ~DuneAmiraMesh()
    {
      if (usedFloatCoords)
        delete[] usedFloatCoords;
      if (usedTetI)
        delete[] usedTetI;
      if (usedBounds)
        delete[] usedBounds;
      if (newIds)
        delete[] newIds;
    };

  void Info(FILE *f)
    {
      mesh->info(f);
    };

  int GetNMatrials() { return nMaterials; };
  int GetMaterialId(const char *name) { return matsMap[name]; };
  const char *GetMaterialName(int id)
    {
      std::map<std::string,int>::iterator iter;
      for (iter=matsMap.begin(); iter!=matsMap.end(); iter++)
        if (iter->second==id)
          return iter->first.c_str();
      return 0;
    };
  int GetNBoundaryIds() { return nBoundaryIds; };
  int GetBoundaryId(const char *name) { return bIdsMap[name]; };
  const char *GetBoundaryIdName(int id)
    {
      std::map<std::string,int>::iterator iter;
      for (iter=bIdsMap.begin(); iter!=bIdsMap.end(); iter++)
        if (iter->second==id)
          return iter->first.c_str();
      return 0;
    };

  int GetNPoints() { return nPoints; };
  int GetNTet() { return nTet; };
  int GetNBdrTr() { return nBndTr; };
  FloatCoord *GetFloatCoord()
    {
      if (floatCoords==0)
        {
    	  AmiraMesh::Data* coordinateData = mesh->findData("Nodes", HxFLOAT, 3, "Coordinates");
          if (coordinateData==0)
            {
              throw -34104;
            }
          floatCoords = (FloatCoord*) coordinateData->dataPtr();
        }
      return floatCoords;
    };
  TetrahedraIndices *GetTetrahedraIndices()
    {
      if (tetI==0)
        {
          AmiraMesh::Data* tetrahedronData =  mesh->findData("Tetrahedra", HxINT32, 4, "Nodes");
          if (tetrahedronData==0)
            {
              throw -34105;
            }
	      tetI = (TetrahedraIndices*) tetrahedronData->dataPtr();
	      for (int k=0; k<nTet; k++)
	        for (int i=0; i<4; i++)
	            tetI[k][i]--;
        }
      return tetI;
    };
  TriangleIndices *GetBoundaryTriangleIndices()
    {
      if (bndTrI==0)
        {
          AmiraMesh::Data* triangleData =  mesh->findData("BoundaryTriangles", HxINT32, 3, "Nodes");
          if (triangleData==0)
            {
              throw -34106;
            }
	      bndTrI = (TriangleIndices*) triangleData->dataPtr();
	      for (int k=0; k<nBndTr; k++)
	        for (int i=0; i<3; i++)
	          bndTrI[k][i]--;
        }
      return bndTrI;
    };
  unsigned char *GetMatIds()
    {
	  if (matIds==0)
	    {
		  AmiraMesh::Data *matIdsData = mesh->findData("Tetrahedra", HxBYTE, 1, "Materials");
		  if (matIdsData)
			matIds = (unsigned char *)matIdsData->dataPtr();
		}
	  return matIds;
    }
  unsigned char *GetBndIds()
    {
	  if (bndIds==0)
	    {
		  AmiraMesh::Data *bndIdsData = mesh->findData("BoundaryTriangles", HxBYTE, 1, "Id");
		  if (bndIdsData)
			bndIds = (unsigned char *)bndIdsData->dataPtr();
		}
	  return bndIds;
    }
  unsigned char *GetMatBndIds()
    {
	  if (tetBndIds==0)
	    {
		  AmiraMesh::Data *tetrahedronBndIds = mesh->findData("Tetrahedra", HxBYTE, 5, "BndIds");
		  if (tetrahedronBndIds)
		    tetBndIds = (unsigned char *)tetrahedronBndIds->dataPtr();
		}
	  return tetBndIds;
    }
  void SelectAll() { for (int k=0; k<256; k++) selected[k] = true; };
  void UnSelectAll() { for (int k=0; k<256; k++) selected[k] = false; };
  bool Select(unsigned int k)
    {
      if (k>255)
        throw -34102;
      bool rc = selected[k];
      selected[k] = true;
      return rc;
    };
  void CountTet(FILE *f = 0)
    {
	  if (matIds==0)
	    GetMatIds();

	  int k;
	  for (k=0; k<256; k++)
	    cnts[k] = 0;
	  for (k=0; k<nTet; k++)
	    cnts[matIds[k]]++;
	  if (f!=0)
	    for (k=0; k<256; k++)
	      if (cnts[k]>0)
	        fprintf(f, "%3d %6d/%d %c\n", k, cnts[k], nTet, selected[k]?'*':'-');
    };
  void CountUsedTet(FILE *f = 0)
    {
	  if (usedMatIds==0)
	    throw -34103;

	  int k;
	  for (k=0; k<256; k++)
	    cnts[k] = 0;
	  for (k=0; k<nUsedTet; k++)
	    cnts[usedMatIds[k]]++;
	  if (f!=0)
	    for (k=0; k<256; k++)
	      if (cnts[k]>0)
	        fprintf(f, "%3d %6d/%d\n", k, cnts[k], nUsedTet);
    };
  void CountUsedTr(FILE *f = 0)
    {
	  if (usedBndIds==0)
	    throw -34104;

	  int i, k;
      int magic[] = {1,2,3,0};
	  for (k=0; k<256; k++)
	    cnts[k] = 0;
// 	  for (k=0; k<nUsedBndTr; k++)
// 	    cnts[usedBndIds[k]]++;
	  for (k=0; k<nUsedTet; k++)
	    for (i=0; i<4; i++)
	      if (usedBounds[k][i]==-1)
	        {
	          cnts[usedBndCharId[k*5+magic[i]+1]]++;
	        }
	  
	  if (f!=0)
	    for (k=0; k<256; k++)
	      if (cnts[k]>0)
	        fprintf(f, "%3d %6d/%d\n", k, cnts[k], nUsedBndTr);
    };
  void Restrict()
    {
      if (floatCoords==0)
        GetFloatCoord();
      if (tetI==0)
        GetTetrahedraIndices();
	  if (matIds==0)
	    GetMatIds();
// 	  if (bndTrI==0)
// 	    GetBoundaryTriangleIndices();
	  if (bndIds==0)
	    GetBndIds();

	  int i, j, k, count;
	  std::vector<unsigned int> newId(nPoints);
	  std::vector<bool> used(nPoints);
	  for (k=0; k<nPoints; k++)
	    {
	      used[k] = false;
	      newId[k] = 9999999;
	    }
	  nUsedTet = 0;
	  for (k=0; k<nTet; k++)
	    {
	      if (selected[matIds[k]])
	        {
	          for (i=0; i<4; i++)
	            {
	              if ((tetI[k][i]<0)||(tetI[k][i]>=nPoints))
	                {
	                  printf("Ojeh: k=%d, i=%d, index=%d\n", k, i, tetI[k][i]);
	                  throw -34103;
	                }
	              used[tetI[k][i]] = true;
	            }
	          nUsedTet++;
	        }
	    }
	  nUsedPoints = 0;
	  for (k=0; k<nPoints; k++)
	    {
	      if (used[k])
	        {
	          newId[k] = nUsedPoints++;
	        }
	    }

      if (usedFloatCoords)
        delete[] usedFloatCoords;
      usedFloatCoords = new FloatCoord[nUsedPoints];

      count = 0;
	  for (k=0; k<nPoints; k++)
	    {
		  newIds[k] = -1;
		  if (used[k])
		    {
			  for (i=0; i<3; i++)
			    usedFloatCoords[count][i] = floatCoords[k][i];
			  newIds[k] = count;
			  count++;
		    }
	    }

      if (usedTetI)
        delete[] usedTetI;
      usedTetI = new TetrahedraIndices[nUsedTet];
      usedMatIds = new unsigned char[nUsedTet];
      if (usedBounds)
        delete[] usedBounds;
      usedBounds = new TetrahedraIndices[nUsedTet];

      count = 0;
	  for (k=0; k<nTet; k++)
	    {
		  if (selected[matIds[k]])
		    {
			  for (i=0; i<4; i++)
			    {
			      usedTetI[count][i] =  newId[tetI[k][i]];
			      usedBounds[count][i] =  0;
			    }
			  usedMatIds[count] = matIds[k];
			  count++;
		    }
	    }
	  nUsedBndTr = ComputeBoundary(nUsedPoints, nUsedTet, &usedTetI[0][0], &usedBounds[0][0]);

	  usedBndTrI = new TriangleIndices[nUsedBndTr];
	  usedBndIds = new unsigned char[nUsedBndTr];

      count = 0;
	  for (k=0; k<nUsedTet; k++)
	    {
		  for (i=0; i<4; i++)
		    {
			  if (usedBounds[k][i]==-1)
			    {
			      int pos = 0;
			      for (j=0; j<4; j++)
			        {
			          if (j!=i)
			            usedBndTrI[count][pos++] = usedTetI[k][j];
			        }
			      count++;
			    }
		    }
	    }

	  SetMaterialsIds(nPoints, newIds, nUsedPoints, nBndTr, &bndTrI[0][0], bndIds,
	                  nUsedBndTr, &usedBndTrI[0][0], usedBndIds);

      usedBndCharId = new unsigned char[5*nUsedTet];
      int ii, magic[] = {1,2,3,0};
      count = 0;
	  for (k=0; k<nUsedTet; k++)
	    {
		  usedBndCharId[k*5] = usedMatIds[k];
		  for (i=0; i<4; i++)
		    {
			  ii = magic[i];
			  if (usedBounds[k][i]==-1)
			    {
			      usedBndCharId[k*5+ii+1] = usedBndIds[count];
			      count++;
			    }
			  else
			    {
			      usedBndCharId[k*5+ii+1] = -1;
			    }
		    }
	     }
	  return;
    };
  int GetUsedNPoints() { return nUsedPoints; };
  int GetUsedNTet() { return nUsedTet; };
  int GetUsedNBdrTr() { return nUsedBndTr; };
  void InsertUGGrid(Dune::GridFactory<Grid> &factory)
    {
	  int k;
	  Dune::FieldVector<double,dim> v;
	  std::vector<unsigned int> vid(4);

      if (floatCoords==0)
        GetFloatCoord();
      if (tetI==0)
        GetTetrahedraIndices();

//      grid.createBegin();
	  for (k=0; k<nPoints; k++)
		{
		  v[0] = floatCoords[k][0];
		  v[1] = floatCoords[k][1];
		  v[2] = floatCoords[k][2];
		  factory.insertVertex(v);
		}

	  for (k=0; k<nTet; k++)
		{
		  vid[0] = tetI[k][0];
		  vid[1] = tetI[k][1];
		  vid[2] = tetI[k][2];
		  vid[3] = tetI[k][3];
		  factory.insertElement(Dune::GeometryType(Dune::GeometryType::simplex,3),vid);
		}
//	  grid.createEnd();
	  
      return;
    }
  template <class Functional>
  void CheckRenumbering(Grid &grid, Functional &F)
    {
	  int i, k;

	  if (renumberedUGIds==0)
	    renumberedUGIds = new int[nPoints];
	  if (revRenumberedUGIds==0)
	    revRenumberedUGIds = new int[nPoints];

      typedef typename Functional::AnsatzVars::IndexSet  IndexSet;
      typedef typename IndexSet::template Codim<0>::template Partition<Dune::All_Partition>::Iterator ElementLeafIterator;
      ElementLeafIterator tetIt(grid.template leafbegin<0>());

      for (k=0; k<nPoints; k++)
        {
          renumberedUGIds[k] = -1;
          revRenumberedUGIds[k] = -1;
        }
      for (k=0;  tetIt!=grid.template leafend<0>(); ++tetIt, k++)
        {
          for (i=0; i<tetIt->geometry().corners(); i++)
            {
              int newPos = grid.leafIndexSet().template subIndex<dim>(*tetIt,i);
              revRenumberedUGIds[tetI[k][i]] = newPos;
              renumberedUGIds[newPos] = tetI[k][i];
            }
        }
      for (k=0; k<nPoints; k++)
        {
          if (renumberedUGIds[k] == -1)
            {
              printf ("Failed to find new number for %d\n", k);
              abort();
            }
          if (revRenumberedUGIds[k] == -1)
            {
              printf ("Failed to find new number for %d\n", k);
              abort();
            }
        }
      printf("CheckRenumbering succeeded\n");
    }
  void InsertRestrictedUGGrid(Grid *grid)
    {
      return;
    }
  void WriteRestrictedAlberta(const char *fileName)
    {
      FILE *f = fopen(fileName, "w");
      int i, k, lineCount = 0;

      fprintf(f, "DIM: %d\n", dim);
      fprintf(f, "DIM_OF_WORLD: %d\n", dim);

      fprintf(f, "number of vertices: %d\n", nUsedPoints);
      fprintf(f, "number of elements: %d\n", nUsedTet);

      fprintf(f, "vertex coordinates:\n");
      lineCount += 5;

	  for (k=0; k<nUsedPoints; k++)
	    {
		  for (i=0; i<3; i++)
			fprintf(f, " %f", usedFloatCoords[k][i]);
		  fprintf(f, "\n");
		  lineCount++;
		}

      fprintf(f, "element vertices:\n");
      lineCount++;
	  for (k=0; k<nUsedTet; k++)
	    {
		  for (i=0; i<4; i++)
			fprintf(f, " %d", usedTetI[k][i]);
		  fprintf(f, "\n");
		  lineCount++;
		}

      fprintf(f, "element boundaries:\n");
      lineCount++;
	  for (k=0; k<nUsedTet; k++)
	    {
		  for (i=0; i<4; i++)
		    fprintf(f, " %d", (usedBounds[k][i]==-1)?1:0);
		  fprintf(f, "\n");
		  lineCount++;
	    }

      fclose(f);
      std::cout << "Alberta file " << fileName << ", " << lineCount << " lines written\n";
    }
  void WriteRestrictedAmira(const char *fileName)
    {
	  int i, k;

	  AmiraMesh *outMesh = new AmiraMesh;
	  outMesh->parameters = mesh->parameters;

	  McPrimType floatType(&usedFloatCoords[0][0]);
	  AmiraMesh::Location outVert("Nodes", nUsedPoints);
	  AmiraMesh::Data outCoord("Coordinates", &outVert, floatType, 3, &usedFloatCoords[0][0]);
	  outMesh->insert(&outVert);
	  outMesh->insert(&outCoord);

	  for (k=0; k<nUsedTet; k++)
	    {
		  for (i=0; i<4; i++)
			usedTetI[k][i]++;
		}

	  for (k=0; k<nUsedBndTr; k++)
	    {
		  for (i=0; i<3; i++)
			usedBndTrI[k][i]++;
		}

	  McPrimType intType(&usedTetI[0][0]);
	  AmiraMesh::Location outTet("Tetrahedra", nUsedTet);
	  AmiraMesh::Data outNodes("Nodes", &outTet, intType, 4, &usedTetI[0][0]);
	  outMesh->insert(&outTet);
	  outMesh->insert(&outNodes);

	  McPrimType byteType(&usedMatIds[0]);
	  AmiraMesh::Location outTetData("TetrahedronData", nUsedTet);
	  AmiraMesh::Data outBndIds("BndIds", &outTetData, byteType, 5, &usedBndCharId[0]);
	  outMesh->insert(&outBndIds);

	  AmiraMesh::Location outBtr("BoundaryTriangles", nUsedBndTr);
	  AmiraMesh::Data outBtrNodes("Nodes", &outBtr, intType, 3, &usedBndTrI[0][0]);
	  outMesh->insert(&outBtr);
	  outMesh->insert(&outBtrNodes);

	  AmiraMesh::Location outBtrData("BoundaryTriangleData", nUsedBndTr);
	  AmiraMesh::Data outBtrId("Id", &outBtrData, byteType, 1, &usedBndIds[0]);
	  outMesh->insert(&outBtrId);

	  outMesh->write(fileName, 1);

	  for (k=0; k<nUsedBndTr; k++)
	    {
		  for (i=0; i<3; i++)
			usedBndTrI[k][i]--;
		}

	  for (k=0; k<nUsedTet; k++)
	    {
		  for (i=0; i<4; i++)
			usedTetI[k][i]--;
		}
	  delete outMesh;
    }
  template <class VertexIterator, class ElementIterator>
  void WriteAmiraSolution(Grid *grid, const char *fileName, int size, int nSol, int solSize, double *x)
    {

      if (matIds==0) GetMatIds();
      if (bndTrI==0) GetBoundaryTriangleIndices();
      if (bndIds==0) GetBndIds();
      if (tetBndIds==0) GetMatBndIds();

	  int i, k;

	  AmiraMesh *outMesh = new AmiraMesh;
	  outMesh->parameters = mesh->parameters;

	  McPrimType floatType(&floatCoords[0][0]);
	  AmiraMesh::Location outVert("Nodes", nPoints);
	  AmiraMesh::Data outCoord("Coordinates", &outVert, floatType, 3);

	  outMesh->insert(&outVert);
	  outMesh->insert(&outCoord);

      VertexIterator vertex    = grid->template leafbegin<dim>();
      VertexIterator endvertex = grid->template leafend<dim>();
    
      for (; vertex!=endvertex; ++vertex)
        {
        
          int index = grid->leafIndexSet().template index<dim>(*vertex);
          const Dune::FieldVector<double, dim>& coords = vertex->geometry()[0];
        
        // Copy coordinates
          for (int i=0; i<dim; i++)
            ((float*)outCoord.dataPtr())[dim*index+i] = coords[i];
        
        }

	  McPrimType intType(&tetI[0][0]);
	  AmiraMesh::Location outTet("Tetrahedra", nTet);
	  AmiraMesh::Data outNodes("Nodes", &outTet, intType, 4);
	  outMesh->insert(&outTet);
	  outMesh->insert(&outNodes);

      int *dPtr = (int*)outNodes.dataPtr();
      ElementIterator eIt    = grid->template leafbegin<0>();
      ElementIterator eEndIt = grid->template leafend<0>();
   
      for (i=0; eIt!=eEndIt; ++eIt)
        {
          if (eIt->type().isTetrahedron())
            {
           
              for (int j=0; j<4; j++) 
                 dPtr[i++] = grid->leafIndexSet().template subIndex<dim>(*eIt,j)+1;
           
            }
          else
            {
              std::cout << "Can't write GeometryType " << eIt->type() << "\n";
              abort();
            }
            
        }

	  McPrimType byteType(&bndIds[0]);
	  AmiraMesh::Location outMatLoc("Tetrahedra", nTet);
	  AmiraMesh::Data outMatData("Materials", &outMatLoc, byteType, 1);
	  outMesh->insert(&outMatData);
	  unsigned char *mIds = (unsigned char*)outMatData.dataPtr();
	  for (k=0; k<nTet; k++)
	    mIds[k] = tetBndIds[5*k];

// 
// 	  AmiraMesh::Location outBtr("BoundaryTriangles", nBndTr);
// 	  AmiraMesh::Data outBtrNodes("Nodes", &outBtr, intType, 3, &bndTrI[0][0]);
// 	  outMesh->insert(&outBtr);
// 	  outMesh->insert(&outBtrNodes);
// 
// 	  AmiraMesh::Location outBtrLoc("BoundaryTriangles", nBndTr);
// 	  AmiraMesh::Data outBtrId("Id", &outBtrLoc, byteType, 1, &bndIds[0]);
// 	  outMesh->insert(&outBtrId);

	  AmiraMesh::Location outU("Nodes", nPoints);
	  AmiraMesh::Data outUData("UValues", &outU, floatType, 3);
	  outMesh->insert(&outUData);
	  float *u = (float*)outUData.dataPtr();
	  for (k=0; k<3*nPoints; k++)
	    u[k] = x[k];
	  
	  AmiraMesh::Location outV("Nodes", nPoints);
	  AmiraMesh::Data outVData("VValues", &outV, floatType, 3);
	  outMesh->insert(&outVData);
	  float *v = (float*)outVData.dataPtr();
	  for (k=0; k<3*nPoints; k++)
	    v[k] = x[3*nPoints+k];

	  AmiraMesh::Field outUField("u", 3, floatType, AmiraMesh::t_linear, &outUData);
	  outMesh->insert(&outUField);
	  AmiraMesh::Field outVField("v", 3, floatType, AmiraMesh::t_linear, &outVData);
	  outMesh->insert(&outVField);

	  outMesh->write(fileName, 1);
//	  delete outMesh;
      printf("WriteAmiraSolution: written to %s\n", fileName);
    }
};

#endif
