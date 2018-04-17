# Define the makefile variables that will be used in building system

F90 = gfortran 
CC  = gcc
FFLAGS = -cpp -O3
CFLAGS = -O3 -fopenmp
FUNC = -DMVV
MYLIBS = -L./openmp -lopenmp_bblas -L./lbstime -llbstime 
SYSLIBS = -lgomp -lpapi
