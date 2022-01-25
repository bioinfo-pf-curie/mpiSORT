
bin_PROGRAMS = mpiSORT
noinst_HEADERS = $(HDR)

mpiSORT_SOURCES =


if HAVE_HTSLIB
AM_CPPFLAGS = -I@LIBHTSDIR@/htslib -I@LIBHTSDIR@ -DHAVE_HTSLIB
AM_LDFLAGS = -lhts -lrt -lpthread
mpiSORT_SOURCES += $(SRC)
mpiSORT_LDADD = @LIBHTSDIR@/libhts.a
mpiSORT_LDFLAGS = -I@LIBHTSDIR@/htslib -I@LIBHTSDIR@
else
mpiSORT_SOURCES += $(SRC)
endif

SRC = mpiSort.c mpiSortUtils.c mergeSort.c parallelBitonicSort.c \
        parallelBitonicSort2.c parallelBitonicSort3.c qkSort.c parser.c \
        tokenizer.c preWrite.c write.c  writeUtils.c sortAnyDim.c

HDR = mpiSortUtils.h mergeSort.h parallelBitonicSort.h \
        parallelBitonicSort2.h parallelBitonicSort3.h qkSort.h \
        parser.h tokenizer.h preWrite.h write.h writeUtils.h \
        bgzf2.c bgzf2.h khash.h sortAnyDim.h \
        compat.h




