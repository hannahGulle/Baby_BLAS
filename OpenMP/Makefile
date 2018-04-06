# Top level makefile
#
include Makefile.inc

all : driver openmp lbstime 

driver: driver.o openmp lbstime
	$(F90) driver.o -o driver $(MYLIBS) $(SYSLIBS)  

driver.o: driver.f90
	$(F90) $(FFLAGS) $(FUNC) driver.f90 -c  

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
