# Top level makefile
#
include Makefile.inc

all : driver serial pthreads lbstime

driver: driver.o serial lbstime pthreads
	$(F90) driver.o -o driver $(LDLIBS)  

driver.o: driver.f90
	$(F90) $(FFLAGS) driver.f90 -c  

serial: 
	cd serial && $(MAKE)

pthreads: 
	cd pthreads && $(MAKE)

lbstime: 
	cd lbstime && $(MAKE)

clean:
	cd serial && $(MAKE) clean
	cd pthreads && $(MAKE) clean
	cd lbstime && $(MAKE) clean
	rm *.o
	touch *.f90

pristine:
	cd serial && $(MAKE) pristine 
	cd pthreads && $(MAKE) pristine 
	cd lbstime && $(MAKE) pristine
	rm *.o	
	rm driver
	touch *.f90

#This next target get "made" every time
.PHONY: serial pthreads lbstime
