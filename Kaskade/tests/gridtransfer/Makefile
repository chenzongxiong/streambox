FLAGSO=-Wall -Wno-sign-compare -Wno-deprecated -fPIC  -O3 -funroll-loops -DNDEBUG
FLAGSG=-Wall -Wno-sign-compare -Wno-deprecated -fPIC  -g
FLAGS=$(FLAGSG)
# CXX=/local/bzfweise/gcc-4.0/installed/bin/g++ $(FLAGS)
CXX = g++ $(FLAGS)
DLIB =

DUNE=/work/bzfweise/dune/installed
FUSION=/work/bzfweise/fusion/fusion_v2
BOOST=/work/bzfweise/boost-1.33
zibDUNE=$(HOME)/dune/zibDUNE

OBJ=umfpack_solve.o


INCLUDE=-I. -I$(zibDUNE)/fem -I$(DUNE)/include -I$(DUNE)/include/dune -I$(FUSION) -I$(BOOST)

DUNELIBS=$(DUNE)/lib/*.a
UGLIBS=-lug3 -lug2 -lgg3 -lgg2 -ldevS -ldevX -ldomS3 -ldomS2

BLAS=-lblas -lgfortran -lg2c -lm
UMFPACK=-lumfpack -lamd -llapack $(BLAS)


LIB=-L$(DUNE)/lib -L$(zibDUNE)/fem -L.


default: gt

gt: gt.cpp ../../fem/fetransfer.hh
	$(CXX) $(INCLUDE) $< $(DUNELIBS) $(LIB) $(zibDUNE)/fem/libzibfem.so $(UGLIBS) $(UMFPACK) $(DLIB) -o $@ 
