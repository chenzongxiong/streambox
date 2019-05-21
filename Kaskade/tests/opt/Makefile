

include ../../Makefile.Local
include ../../Makefile.Rules

ALGPATH=$(KASKADE7)/algorithm/


TARGET = mopt
OPTIMIZER = hyp_ip
OBJ= $(KASKADE7)/libs/umfpack_solve.o $(KASKADE7)/libs/mumps_solve.o $(KASKADE7)/libs/pardiso_solve.o 
INCLUDE = $(DUNEINC) $(UGINC) $(FUSIONINC) $(BOOSTINC) $(UMFPACKINC) $(IPOPTINC) $(AMIRAINC) $(KASKADEINC) $(MUMPSINC) $(SUPERLUINC) $(ALGPATH) $(ALGPATH)/opt_aux/src/include/
FLAGS = $(WFLAGS) $(DEBUGFLAGS)

$(TARGET): $(TARGET).o  Makefile $(KASKADE7)/libs/libkaskade.a 
	$(CXX) $< $(OBJ) \
$(KASKADELIB) \
$(PARDISOLIB)\
 $(BOOSTLIB)\
$(MUMPSLIB)\
 $(UGLIB)\
 $(UMFPACKLIB)\
 $(SUPERLULIB)\
$(PARDISOLIB)\
$(BLASLIB) $(FTNLIB) -o $@

$(TARGET).o: $(TARGET).cpp $(ALGPATH)searchspace.hh $(ALGPATH)$(OPTIMIZER).hh testfunctions.hh  Makefile

clean:
	rm $(TARGET) $(TARGET).o
