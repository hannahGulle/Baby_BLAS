# Define the makefile variables that will be used in building system

F90 = gfortran
CC  = gcc
FFLAGS = -cpp -O0
CFLAGS = -O0
FUNC = -DMMM
LDLIBS = -L./serial -lserial_bblas -L./lbstime -llbstime
SYSLIBS = -lpapi
