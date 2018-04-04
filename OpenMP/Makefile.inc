# Define the makefile variables that will be used in building system

F90 = gfortran 
CC  = gcc
FFLAGS = -cpp -O2 
CFLAGS = -O2 -fopenmp
MYLIBS = -L./openmp -lopenmp_bblas -L./pthreads -lpthreads_bblas\
         -L./serial -lserial_bblas -L./lbstime -llbstime 
SYSLIBS = -lgomp -lpthread
