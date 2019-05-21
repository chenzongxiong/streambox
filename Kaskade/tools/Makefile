include ../Makefile.Local
include ../Makefile.Rules

FLAGS =  -DKASKADE7=\"$(KASKADE7NAME)\" $(WFLAGS) $(OPTFLAGS) $(CLANGINC)
KASKADE7NAME = $(shell basename $(KASKADE7))

OBJ = filterOutput.o sumDiffFile.o
SRC = filterOutput.cpp sumDiffFile.cpp

all: filterOutput sumDiffFile

filterOutput: filterOutput.o
	$(CXX) $< $(FTNLIB) -o $@

sumDiffFile: sumDiffFile.o
	$(CXX) $< $(FTNLIB) -o $@
	
remove_install_deps: remove_install_deps.o
	$(CXX) $< $(FTNLIB) -o $@

install: filterOutput sumDiffFile
	if [ ! -d "${HOME}/bin" ]; then mkdir -p ${HOME}/bin; fi
	cp -p filterOutput sumDiffFile ${HOME}/bin

clean:
	rm -f filterOutput sumDiffFile remove_install_deps gccerr.txt *.o *.bak *.tmp
	
# DO NOT DELETE

