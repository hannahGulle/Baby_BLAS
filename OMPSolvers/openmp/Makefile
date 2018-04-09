include ../Makefile.inc

# Timing Library targets 

all: lib

# $(CFLAGS) from previous directory Makefile.inc
ils.o : ils.c
	$(CC) $(CFLAGS) -fopenmp -c ils.c

dls.o : dls.c
	$(CC) $(CFLAGS) -fopenmp -c dls.c

lib: dls.o ils.o
	ar -rc libopenmp_bblas.a dls.o ils.o
	ranlib libopenmp_bblas.a

# Default Targets for Cleaning up the Environment

clean :
	rm *.o
	touch *.c

pristine:
	make clean
	rm *.a
