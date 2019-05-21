Makefile.Local:
	false "Please install Makefile.Local, use Makefile-htc.Local or Makefile-mac.Local as template"

libs:
	mkdir libs

kasklib: Makefile.Local libs
	cd algorithm; make
	cd fem; make
	cd fem/diffops; make
	cd io; make
	cd linalg; make
	cd mg; make
	cd timestepping; make
	cd utilities; make

depend: Makefile.Local remove_install_deps
	cd algorithm; make -f Makefile.gen depend
	cd fem; make -f Makefile.gen depend
	cd fem/diffops; make -f Makefile.gen depend
	cd io; make -f Makefile.gen depend
	cd linalg; make -f Makefile.gen depend
	cd mg; make -f Makefile.gen depend
	cd timestepping; make -f Makefile.gen depend
	cd utilities; make -f Makefile.gen depend
	cd tutorial/laplacian; make -f Makefile.gen depend
	cd tutorial/stationary_heattransfer; make -f Makefile.gen depend
	cd tutorial/instationary_heattransfer; make -f Makefile.gen depend
	cd tutorial/stokes; make -f Makefile.gen depend
	cd tutorial/instationary_NavierStokes; make -f Makefile.gen depend
	cd tutorial/Embedded_errorEstimation; make -f Makefile.gen depend
	cd tutorial/HB_errorEstimation; make -f Makefile.gen depend
	cd tutorial/sst_pollution; make -f Makefile.gen depend
	cd tutorial/elastomechanics; make -f Makefile.gen depend
	cd tutorial/artificial_1d_testProblem; make -f Makefile.gen depend
	cd tests; make -f Makefile.gen depend
	cd tests/poisson; make -f Makefile.gen depend
	cd tests/statHeatTransfer; make -f Makefile.gen depend
	cd tests/embErrorEst; make -f Makefile.gen depend
	cd tests/reentrantCorner; make -f Makefile.gen depend
	cd tests/hbPrec; make -f Makefile.gen depend

clean: Makefile.Local
	cd algorithm; make clean
	cd fem; make clean
	cd fem/diffops; make clean
	cd io; make clean
	cd linalg; make clean
	cd mg; make clean
	cd timestepping; make clean
	cd utilities; make clean
	cd tutorial/laplacian; make clean
	cd tutorial/stationary_heattransfer; make clean
	cd tutorial/instationary_heattransfer; make clean
	cd tutorial/stokes; make clean
	cd tutorial/instationary_NavierStokes; make clean
	cd tutorial/Embedded_errorEstimation; make clean
	cd tutorial/HB_errorEstimation; make clean
	cd tutorial/sst_pollution; make clean
	cd tutorial/elastomechanics; make clean
	cd tutorial/artificial_1d_testProblem; make clean
	cd tests; make clean
	cd tests/poisson; make clean
	cd tests/statHeatTransfer; make clean
	cd tests/embErrorEst; make clean
	cd tests/reentrantCorner; make clean
	cd tests/hbPrec; make clean
	cd tools; make clean
	rm -f remove_install_deps *.bak libs/*.o libs/*.a

cleantutorial: Makefile.Local
	cd tutorial/laplacian; make clean
	cd tutorial/stationary_heattransfer; make clean
	cd tutorial/instationary_heattransfer; make clean
	cd tutorial/stokes; make clean
	cd tutorial/instationary_NavierStokes; make clean
	cd tutorial/Embedded_errorEstimation; make clean
	cd tutorial/HB_errorEstimation; make clean
	cd tutorial/sst_pollution; make clean
	cd tutorial/elastomechanics; make clean
	cd tutorial/artificial_1d_testProblem; make clean

remove_install_deps: Makefile.Local tools/remove_install_deps.o
	cd tools; make remove_install_deps; cp -p remove_install_deps ..

tools/remove_install_deps.o: Makefile.Local
	cd tools; make remove_install_deps.o;

install: Makefile.Local
	make depend; make kasklib

include Makefile.Local
KASKADE7NAME = $(shell basename $(KASKADE7))

tutorial: Makefile.Local
	cd tutorial/laplacian; make; export DYLD_LIBRARY_PATH=$(BOOST)/lib; ./laplace
	cd tutorial/stationary_heattransfer; make; export DYLD_LIBRARY_PATH=$(BOOST)/lib; ./heat; ./heat3d
	cd tutorial/artificial_1d_testProblem; make; export DYLD_LIBRARY_PATH=$(BOOST)/lib; ./atp --graphicalOutput 0; 
	cd tutorial/Embedded_errorEstimation; make; export DYLD_LIBRARY_PATH=$(BOOST)/lib; ./heat
	cd tutorial/HB_errorEstimation; make; export DYLD_LIBRARY_PATH=$(BOOST)/lib; ./heat; ./heat3d --file cube.am
	cd tutorial/sst_pollution; make sst sst_nleqErr sst_giant; \
	export DYLD_LIBRARY_PATH=$(BOOST)/lib; \
	./sst; ./sst_nleqErr; ./sst_giant 
	cd tutorial/stokes; make; export DYLD_LIBRARY_PATH=$(BOOST)/lib; \
	./stokes --solver.direct MUMPS; \
	./stokes --solver.direct UMFPACK; \
	./stokes --solver.direct UMFPACK3264; \
	./stokes --solver.direct SUPERLU;
	cd tutorial/elastomechanics; make; export DYLD_LIBRARY_PATH=$(BOOST)/lib; ./elastomechanics; 
	cd tutorial/instationary_heattransfer; make; export DYLD_LIBRARY_PATH=$(BOOST)/lib; ./ms;
	cd tutorial/instationary_NavierStokes; make; export DYLD_LIBRARY_PATH=$(BOOST)/lib; ./navierStokes

manual: Makefile.Local
	cd tutorial/ProgrammersManual; pdflatex manual.tex; pdflatex manual.tex

test: Makefile.Local
	cd tests; make tests

benchmarks: tests

distribution: Makefile.Local
	make clean
	cd algorithm; rm Makefile
	cd fem; rm Makefile
	cd fem/diffops; rm Makefile
	cd io; rm Makefile
	cd linalg; rm Makefile
	cd mg; rm Makefile
	cd timestepping; rm Makefile
	cd utilities; rm Makefile
	cd tutorial/laplacian; rm Makefile
	cd tutorial/stationary_heattransfer; rm Makefile
	cd tutorial/instationary_heattransfer; rm Makefile
	cd tutorial/stokes; rm Makefile
	cd tutorial/instationary_NavierStokes; rm Makefile
	cd tutorial/Embedded_errorEstimation; rm Makefile
	cd tutorial/HB_errorEstimation; rm Makefile
	cd tutorial/sst_pollution; rm Makefile 
	cd tutorial/elastomechanics; rm Makefile
	cd tutorial/artificial_1d_testProblem; rm Makefile
	cd tests; rm Makefile
	cd tests/poisson; rm Makefile
	cd tests/statHeatTransfer; rm Makefile
	cd tests/embErrorEst; rm Makefile
	cd tests/reentrantCorner; rm Makefile
	cd tests/hbPrec; rm Makefile
	rm Makefile.Local
	cd ..; \
	tar --exclude-vcs --exclude $(KASKADE7NAME)/doc/html --exclude $(KASKADE7NAME)/doc/latex \
	--exclude $(KASKADE7NAME)/problems --exclude $(KASKADE7NAME)/experiments \
	--exclude $(KASKADE7NAME)/tutorial/ZR-16-17 \
	--exclude $(KASKADE7NAME)/tests/sequential --exclude $(KASKADE7NAME)/tests/parallel --exclude $(KASKADE7NAME)/tests/thirdparty --exclude $(KASKADE7NAME)/tests/clang \
	--exclude $(KASKADE7NAME)/Makefile-sequential.Local --exclude $(KASKADE7NAME)/Makefile-parallel.Local --exclude $(KASKADE7NAME)/Makefile-thirdparty.Local $(KASKADE7NAME)/Makefile-clang.Local \
	--exclude $(KASKADE7NAME)/setenv.sh --exclude $(KASKADE7NAME)/setenv-parallel.sh --exclude $(KASKADE7NAME)/setenv-sequential.sh \
	--exclude $(KASKADE7NAME)/setenv-thirdparty.sh \
	--exclude $(KASKADE7NAME)/setenv-clang.sh \
	-cf $(KASKADE7NAME).tar $(KASKADE7NAME); \
	gzip $(KASKADE7NAME).tar; \
	tar --exclude-vcs --exclude $(KASKADE7NAME)/doc \
	 --exclude $(KASKADE7NAME)/problems --exclude $(KASKADE7NAME)/experiments \
	--exclude $(KASKADE7NAME)/tutorial/ZR-16-17 \
	--exclude $(KASKADE7NAME)/tests/sequential --exclude $(KASKADE7NAME)/tests/parallel --exclude $(KASKADE7NAME)/tests/thirdparty --exclude $(KASKADE7NAME)/tests/clang \
	--exclude $(KASKADE7NAME)/Makefile-sequential.Local --exclude $(KASKADE7NAME)/Makefile-parallel.Local --exclude $(KASKADE7NAME)/Makefile-thirdparty.Local $(KASKADE7NAME)/Makefile-clang.Local \
	--exclude $(KASKADE7NAME)/setenv.sh --exclude $(KASKADE7NAME)/setenv-parallel.sh --exclude $(KASKADE7NAME)/setenv-sequential.sh \
	--exclude $(KASKADE7NAME)/setenv-thirdparty.sh \
	--exclude $(KASKADE7NAME)/setenv-clang.sh \
	-cf $(KASKADE7NAME)_No_Doc.tar $(KASKADE7NAME); \
	gzip $(KASKADE7NAME)_No_Doc.tar
