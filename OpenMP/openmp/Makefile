include ../Makefile.inc

# OpenMP library targets 

all: lib

mmm.o : mmm.c
	$(CC) $(CFLAGS) -fopenmp -c mmm.c

vvm.o : vvm.c
	$(CC) $(CFLAGS) -fopenmp -c vvm.c

mvv.o : mvv.c
	$(CC) $(CFLAGS) -fopenmp -c mvv.c

dot.o : dot.c
	$(CC) $(CFLAGS) -fopenmp -c dot.c

lib: mmm.o vvm.o mvv.o dot.o
	ar -rc libopenmp_bblas.a  dot.o mmm.o vvm.o mvv.o
	ranlib libopenmp_bblas.a

# Default Targets for Cleaning up the Environment

clean :
	rm *.o
	touch *.c

pristine:
	make clean
	rm *.a
