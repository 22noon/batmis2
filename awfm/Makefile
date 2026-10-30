# Builds the AvxWindowFmIndex (AWFM) based BatMis tools.
# Requires a built AvxWindowFmIndex checkout (see README.batmisAWFM.md).

CXX ?= g++
CXXFLAGS ?= -O3 -funroll-loops
BINDIR ?= ../bin

AWFM_ROOT ?= ../../../AvxWindowFmIndex
AWFM_SRC ?= $(AWFM_ROOT)/src
AWFM_BUILD ?= $(AWFM_ROOT)/build/build
AWFM_FASTAVECTOR_SRC ?= $(AWFM_ROOT)/lib/FastaVector/src
AWFM_INCLUDES = -I$(AWFM_SRC) -I$(AWFM_FASTAVECTOR_SRC)
AWFM_LIBS ?= -L$(AWFM_BUILD) -lawfmindex_static -ldivsufsort64 -lfastavector_static -fopenmp
AWFM_STATIC = $(AWFM_BUILD)/libawfmindex_static.a

PROGRAMS = batmisAWFM batmisAWFM_index_builder

.PHONY: all batman_zero_mismatch_indexX install clean

all: $(PROGRAMS)

%.o: %.cpp $(AWFM_STATIC)
	$(CXX) -std=c++17 -fpermissive $(AWFM_INCLUDES) $(CXXFLAGS) -c $< -o $@

batmisAWFM: batman_zero_mismatch_indexX.o
	$(CXX) -std=c++17 $(CXXFLAGS) -o $@ $< $(AWFM_LIBS) -lm

batmisAWFM_index_builder: batmisAWFM_index_builder.o
	$(CXX) -std=c++17 $(CXXFLAGS) -o $@ $< $(AWFM_LIBS) -lm

# Old name for batmisAWFM
batman_zero_mismatch_indexX: batmisAWFM
	cp batmisAWFM $@

$(AWFM_STATIC):
	@echo "Missing $(AWFM_STATIC)." >&2
	@echo "Build AvxWindowFmIndex first:" >&2
	@echo "  cmake -S $(AWFM_ROOT) -B $(AWFM_ROOT)/build && cmake --build $(AWFM_ROOT)/build -j" >&2
	@echo "or set AWFM_ROOT=<path to AvxWindowFmIndex>." >&2
	@exit 1

install: $(PROGRAMS)
	mkdir -p $(BINDIR)
	cp $(PROGRAMS) $(BINDIR)

clean:
	rm -f *.o $(PROGRAMS) batman_zero_mismatch_indexX
