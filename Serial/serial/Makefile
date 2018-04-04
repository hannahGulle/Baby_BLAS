include ../Makefile.inc

# Timing Library targets 

all: lib

# $(CFLAGS) from previous directory Makefile.inc
vvm.o : vvm.c
	$(CC) $(CFLAGS) -c vvm.c

mmm.o : mmm.c
	$(CC) $(CFLAGS) -DSTRIP8  -c mmm.c

dot.o : dot.c
	$(CC) $(CFLAGS) -c dot.c

mvv.o : mvv.c
	$(CC) $(CFLAGS) -c mvv.c

lib: vvm.o mmm.o dot.o mvv.o
	ar -rc libserial_bblas.a vvm.o mmm.o dot.o mvv.o
	ranlib libserial_bblas.a

# Default Targets for Cleaning up the Environment

clean :
	rm *.o
	touch *.c

pristine:
	make clean
	rm *.a
