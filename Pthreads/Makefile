# Top level makefile
#
include Makefile.inc

all : driver pthreads lbstime

driver: driver.o lbstime pthreads
	$(F90) driver.o -o driver $(LDLIBS)  

driver.o: driver.f90
	$(F90) $(FFLAGS) $(FUNC) driver.f90 -c  

pthreads: 
	cd pthreads && $(MAKE)

lbstime: 
	cd lbstime && $(MAKE)

clean:
	cd pthreads && $(MAKE) clean
	cd lbstime && $(MAKE) clean
	rm *.o
	touch *.f90

pristine:
	cd pthreads && $(MAKE) pristine 
	cd lbstime && $(MAKE) pristine
	rm *.o	
	rm driver
	touch *.f90

#This next target get "made" every time
.PHONY: pthreads lbstime
