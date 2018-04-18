# Define the makefile variables that will be used in building system

F90 = gfortran
CC  = gcc
FFLAGS = -cpp -O3
CFLAGS = -std=c99 -O3
FUNC = -DILS
LDLIBS = -L./pthreads -lpthreads_bblas -L./lbstime -llbstime -lpthread
SYSLIBS = -lpthread
