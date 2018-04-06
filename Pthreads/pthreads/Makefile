include ../Makefile.inc

# Pthread library targets 

all: lib

mmm.o : mmm.c
	$(CC) $(CFLAGS) -c mmm.c

vvm.o : vvm.c
	$(CC) $(CFLAGS) -c vvm.c

mvv.o : mvv.c
	$(CC) $(CFLAGS) -c mvv.c

dot.o : dot.c
	$(CC) $(CFLAGS) -c dot.c

lib: mmm.o vvm.o mvv.o dot.o
	ar -rc libpthreads_bblas.a dot.o mvv.o mmm.o vvm.o
	ranlib libpthreads_bblas.a

# Default Targets for Cleaning up the Environment

clean :
	rm *.o
	touch *.c

pristine:
	make clean
	rm *.a
