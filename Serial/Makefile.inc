# Define the makefile variables that will be used in building system

F90 = gfortran
CC  = gcc
FFLAGS = -cpp -O3
CFLAGS = -O3
FUNC = -DDOT
LDLIBS = -L./serial -lserial_bblas -L./lbstime -llbstime
SYSLIBS = -lpapi
