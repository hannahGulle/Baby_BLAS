include ../Makefile.inc

# Timing Library targets 

all: lib

# $(CFLAGS) from previous directory Makefile.inc
vvm.o : vvm.c
	$(CC) $(CFLAGS) -c vvm.c

mmm.o : mmm.c
	$(CC) $(CFLAGS) -c -DSTRIP8 mmm.c

dot.o : dot.c
	$(CC) $(CFLAGS) -c dot.c

mvv.o : mvv.c
	$(CC) $(CFLAGS) -c mvv.c

dls.o : dls.c
	$(CC) $(CFLAGS) -c dls.c

lib: vvm.o mmm.o dot.o mvv.o dls.c
	ar -rc libserial_bblas.a vvm.o mmm.o dot.o mvv.o dls.c
	ranlib libserial_bblas.a

# Default Targets for Cleaning up the Environment

clean :
	rm *.o
	touch *.c

pristine:
	make clean
	rm *.a
