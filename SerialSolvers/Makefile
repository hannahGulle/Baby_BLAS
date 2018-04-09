# Top level makefile
#
include Makefile.inc

all : lsdriver serial lbstime

lsdriver: lsdriver.o serial lbstime
	$(F90) lsdriver.o -o lsdriver $(LDLIBS) $(SYSLIBS)

lsdriver.o: lsdriver.f90
	$(F90) $(FFLAGS) $(FUNC) lsdriver.f90 -c  

serial: 
	cd serial && $(MAKE)

lbstime: 
	cd lbstime && $(MAKE)

clean:
	cd serial && $(MAKE) clean
	cd lbstime && $(MAKE) clean
	rm *.o
	touch *.f90

pristine:
	cd serial && $(MAKE) pristine 
	cd lbstime && $(MAKE) pristine
	rm *.o	
	rm driver
	touch *.f90

#This next target get "made" every time
.PHONY: serial lbstime
