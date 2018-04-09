# Top level makefile
#
include Makefile.inc

all : lsdriver openmp lbstime

lsdriver: lsdriver.o openmp lbstime
	$(F90) lsdriver.o -o lsdriver $(LDLIBS) $(SYSLIBS)

lsdriver.o: lsdriver.f90
	$(F90) $(FFLAGS) $(FUNC) lsdriver.f90 -c  

openmp: 
	cd openmp && $(MAKE)

lbstime: 
	cd lbstime && $(MAKE)

clean:
	cd openmp && $(MAKE) clean
	cd lbstime && $(MAKE) clean
	rm *.o
	touch *.f90

pristine:
	cd openmp && $(MAKE) pristine 
	cd lbstime && $(MAKE) pristine
	rm *.o	
	rm driver
	touch *.f90

#This next target get "made" every time
.PHONY: openmp lbstime
