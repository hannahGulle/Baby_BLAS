# Define the makefile variables that will be used in building system

F90 = gfortran
CC  = gcc
FFLAGS = -cpp -O2 #-DACCURACY_TEST
CFLAGS = -O2
LDLIBS = -L./pthreads -lpthreads_bblas -L./serial -lserial_bblas -L./lbstime -llbstime -lpthread
