# Define the makefile variables that will be used in building system

F90 = gfortran
CC  = gcc
FFLAGS = -cpp -O0
CFLAGS = -std=c99 -O0
FUNC = -DDOT
LDLIBS = -L./pthreads -lpthreads_bblas -L./lbstime -llbstime -lpthread
SYSLIBS = -lpapi -lpthread
