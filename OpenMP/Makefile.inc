# Define the makefile variables that will be used in building system

F90 = gfortran 
CC  = gcc
FFLAGS = -cpp -O0 
CFLAGS = -O0 -fopenmp
FUNC = -DMMM
MYLIBS = -L./openmp -lopenmp_bblas -L./lbstime -llbstime 
SYSLIBS = -lgomp -lpapi
