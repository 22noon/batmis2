CC ?= gcc
CXX ?= g++
CFLAGS ?= -w -O3 -funroll-loops -mpopcnt
CXXFLAGS ?= -O3 -Wall -mpopcnt
AR ?= ar
RANLIB ?= ranlib

FMINDEX_CORE_OBJS = BWT_fmindex.o BWTConstruct.o MiscUtilities.o MemManager.o TextConverter.o QSufSort.o DNACount.o

.PHONY: all clean fmindex-only-lib

all: libfmindex_core.a bwt_build_index

fmindex-only-lib: libfmindex_core.a

libfmindex_core.a: $(FMINDEX_CORE_OBJS)
	$(AR) cru $@ $(FMINDEX_CORE_OBJS)
	$(RANLIB) $@

BWT_fmindex.o: BWT.c
	$(CC) $(CFLAGS) -I. -c $< -o $@

%.o: %.c
	$(CC) $(CFLAGS) -I. -c $< -o $@

bwt_build_index: bwt_build_index.cpp libfmindex_core.a
	$(CXX) -std=c++17 $(CXXFLAGS) -I. -o $@ bwt_build_index.cpp libfmindex_core.a -lm

clean:
	rm -f *.o libfmindex_core.a bwt_build_index
