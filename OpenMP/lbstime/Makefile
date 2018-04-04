include ../Makefile.inc

# Timing Library targets 

all: lib

walltime.o : walltime.c
	$(CC)  -c walltime.c

cputime.o : cputime.c
	$(CC)  -c cputime.c

lib: cputime.o walltime.o
	ar -rc liblbstime.a cputime.o walltime.o
	ranlib liblbstime.a

# Default Targets for Cleaning up the Environment

clean :
	rm *.o
	touch *.c

pristine:
	make clean
	rm *.a
