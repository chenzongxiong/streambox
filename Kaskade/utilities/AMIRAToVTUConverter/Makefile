include ../../Makefile.Local
include ../../Makefile.Rules

INCLUDE = $(DUNEINC) $(UGINC) $(BOOSTINC) $(AMIRAINC) $(KASKADEINC) $(UMFPACKINC)
FLAGS = $(WFLAGS) $(OPTFLAGS)# $(DEBUGFLAGS)


convert: amira_to_vtu.o Makefile  $(KASKADE7)/libs/libkaskade.a 
	$(CXX) -g $< $(KASKADE7)/libs/umfpack_solve.o \
	$(KASKADELIB)  \
 $(DUNELIB) \
 $(UGLIB) \
 $(BOOSTLIB) \
 $(AMIRALIB) \
 $(UMFPACKLIB) \
 $(BLASLIB) $(FTNLIB) -o $@ 
 
clean:
	rm -f convert amira_to_vtu.o *.bak

depend:
	cp Makefile.gen Makefile; $(KASKADE7)/gccmakedep -f Makefile $(INCLUDE) amira_to_vtu.cpp; ../../remove_install_deps Makefile

# DO NOT DELETE
