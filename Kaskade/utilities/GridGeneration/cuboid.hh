#include "cube.hh"

class cuboid{

public:
		std::vector<Dune::FieldVector<double, 3> > vertices;
        std::vector<std::vector< unsigned int> > tetraeder;
  
	//creates a standart cube with weidth 1, depth 1 and heigth 1
	cuboid(){
		cube cube1;
		this->vertices=cube1.getVertices();
		this->tetraeder=cube1.getTetraeder();
	}

	//creates a cuboid with weidth x , depth y and heigth 1
	cuboid(int x, int y, int z){
		std::vector<cube> cubes;
		for(int i=0; i<x; i++){
			for(int j=0; j<y; j++){
				for(int k=0; k<z; k++){
				  cube cube1(i,j,k);  cubes.push_back(cube1);
				}
			}
		}
		createCuboid(cubes);
		
	}

	//creates a cuboid with an given vector of qubes
	//vertices that are twice or more will be deleted and the tedrader will be updated
	void createCuboid(std::vector<cube> cubes){
		std::vector<Dune::FieldVector<double,3> > verts;
		std::vector<std::vector<unsigned int> > tet;
		std::vector<unsigned int>  new_ind;
		Dune::FieldVector<double,3> totest;
		Dune::FieldVector<double,3> sub;
		std::vector<int> toerase;
		std::vector<bool> erased(verts.size(),false);


		for(int k=0;k<cubes.size();k++){
			for(int i=0;i<9;i++){
				verts.push_back(cubes[k].getVertices()[i]);
				new_ind.push_back(k*9+i);
			}
			for(int i=0;i<12;i++){
				for(int j=0;j<4;j++){
					cubes[k].getTetraeder()[i][j]+=k*9;
				}
				tet.push_back(cubes[k].getTetraeder()[i]);
			}
		}

		for(unsigned int i=0;i<verts.size();i++){
			totest=verts[i];
			for(unsigned int j=i+1;j<verts.size();j++){
				sub=totest-verts[j];
				if(sub.two_norm()<pow(10,-12) && !erased[j]){
					erased[j]=true;
					new_ind[j]=new_ind[i];
					for(int l=j+1;l<new_ind.size();l++){
						if(new_ind[l]>new_ind[j]) new_ind[l]-=1;
					}
					toerase.push_back(j);
				}
			}
		}

		for(int k=0;k<tet.size();k++){
			for(int l=0;l<4;l++){
				tet[k][l]=new_ind[tet[k][l]];
			}
		}
		sort(toerase.begin(),toerase.end());
		for(int k=toerase.size()-1;k>=0;k--){
			verts.erase(verts.begin()+toerase[k]);
		}

		this->vertices=verts;
		this->tetraeder=tet;  
	}

	
		std::vector<Dune::FieldVector<double,3> >& getVerts(){return this->vertices;}
        std::vector<std::vector<unsigned int> >& getTets(){return this->tetraeder;}


	/*std::vector<Dune::FieldVector<double,3> > deldouble(std::vector<Dune::FieldVector<double,3> > tosort){
		std::vector<int> toerase;
		std::vector<bool> erased(tosort.size(),false);
		for(int i=0; i<tosort.size();i++){
			for(int j=i+1;j<tosort.size();j++){
				if(tosort[i]==tosort[j] && !erased[j]){
					erased[j]=true;
					toerase.push_back(j);
				}
			}
		}	
		for(int k=toerase.size()-1;k>=0;k--){
			tosort.erase(tosort.begin()+toerase[k]);
		}
		return tosort;
	}*/

};
