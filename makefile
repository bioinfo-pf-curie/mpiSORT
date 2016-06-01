CC=			mpicc
CXX=		mpic++
CFLAGS=		-DBUG -g -Wall -O2 -o 
WRAP_MALLOC=-DUSE_MALLOC_WRAPPERS
AR=			ar
DFLAGS=		-DMPI-DHAVE_PTHREAD -D_FILE_OFFSET_BITS=64 -D_USE_FILE_OFFSET64 -D_LARGEFILE64_SOURCE -D_USE_MISC $(WRAP_MALLOC)
LOBJS=		
AOBJS=		parabitonicsort.o diffuse.o format.o merge_utils.o merge.o mergeSort.o mpi_globals.o parser.o preWrite.o write2.o tokenizer.o bufferized_read.o
			
PROG=		psort
INCLUDES=	
LIBS=		-lm -lz
SUBDIRS=	.

.SUFFIXES:.c .o .cc

.c.o:
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

psort:$(AOBJS) mpiSort.o
		$(CC) $(CFLAGS) $(DFLAGS) $(AOBJS) mpiSort.o -o $@ -L. $(LIBS)

clean:
		rm -f gmon.out *.o a.out $(PROG) *~ *.a

depend:
	( LC_ALL=C ; export LC_ALL; makedepend -Y -- $(CFLAGS) $(DFLAGS) -- *.c )

# DO NOT DELETE THIS LINE -- make depend depends on it.

parabitonicsort.o: parabitonicsort.h
bufferized_read.o: bufferized_read.h
parser.o: parser.h tokenizer.h
format.o: merge.h parser.h diffuse.h
diffuse.o: diffuse.h merge.h
merge.o: merge.h format.h diffuse.h merge_utils.h 
main_parallel_version.o: parser.h 
