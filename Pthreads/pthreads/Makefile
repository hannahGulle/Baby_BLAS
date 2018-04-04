include ../Makefile.inc

# Pthread library targets 

all: lib

mmm.o : mmm.c
	$(CC) $(CFLAGS)  -c mmm.c

lib: mmm.o 
	ar -rc libpthreads_bblas.a  mmm.o 
	ranlib libpthreads_bblas.a

# Default Targets for Cleaning up the Environment

clean :
	rm *.o
	touch *.c

pristine:
	make clean
	rm *.a
