include ../Makefile.inc

# Timing Library targets 

all: lib

# $(CFLAGS) from previous directory Makefile.inc
ils.o : ils.c
	$(CC) $(CFLAGS) -c ils.c

dls.o : dls.c
	$(CC) $(CFLAGS) -c dls.c

lib: dls.o ils.o
	ar -rc libserial_bblas.a dls.o ils.o
	ranlib libserial_bblas.a

# Default Targets for Cleaning up the Environment

clean :
	rm *.o
	touch *.c

pristine:
	make clean
	rm *.a
