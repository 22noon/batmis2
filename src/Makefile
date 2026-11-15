# BatMis mappers (README.md):
#   batmis1   the search on a pair of one-direction FM indexes (batmis1.cpp)
#   batmis2   the same search on bidirectional SA ranges (batmis2.cpp)
#
#   make                both
#   make test           both test suites (test-batmis1, test-batmis2)
#   make install        to $(BINDIR), with the index builder
# Build variants: DEFS=-DBATMAN_BATCH_STATS (step counts, printed with --timing);
# batmis2_verify checks every bidirectional range against the other index.

CXX ?= g++
CXXFLAGS ?= -O3 -funroll-loops -mpopcnt -pthread
DEFS ?=
BINDIR ?= ../bin

INDEX_DIR ?= ../index
FMINDEX_LIB = $(INDEX_DIR)/libfmindex_core.a
INDEX_BUILDER = $(INDEX_DIR)/bwt_build_index
# htslib (SAM/BAM output), built from source; its configure picked these libraries.
HTSLIB_DIR ?= ../htslib-1.24
HTSLIB = $(HTSLIB_DIR)/libhts.a
HTSLIB_LIBS ?= -llzma -lbz2 -lz -lm -lpthread
TESTS = ../tests

HEADERS = fm_occ.hpp genome_map.hpp sam_writer.hpp
COMPILE = $(CXX) -std=c++17 -fpermissive -I. -I$(INDEX_DIR) -I$(HTSLIB_DIR) $(CXXFLAGS) $(DEFS)
LIBS = $(FMINDEX_LIB) $(HTSLIB) $(HTSLIB_LIBS)

.PHONY: all test test-batmis1 test-batmis2 install clean FORCE

all: batmis1 batmis2

batmis1: batmis1.cpp $(HEADERS) $(FMINDEX_LIB) $(HTSLIB)
	$(COMPILE) -o $@ $< $(LIBS)

batmis2: batmis2.cpp $(HEADERS) $(FMINDEX_LIB) $(HTSLIB)
	$(COMPILE) -o $@ $< $(LIBS)

# batmis2 checking each range's other-index start and each other-row locate (slow).
batmis2_verify: batmis2.cpp $(HEADERS) $(FMINDEX_LIB) $(HTSLIB)
	$(COMPILE) -DBATMIS2_VERIFY -o $@ $< $(LIBS)

# 16-read chunks, so that the tests' few reads span several blocks, chunks and threads.
batmis1_chunk16: batmis1.cpp $(HEADERS) $(FMINDEX_LIB) $(HTSLIB)
	$(COMPILE) -DBATMAN_BATCH_BLOCK=16 -o $@ $< $(LIBS)

batmis2_chunk16: batmis2.cpp $(HEADERS) $(FMINDEX_LIB) $(HTSLIB)
	$(COMPILE) -DBATMAN_BATCH_BLOCK=16 -o $@ $< $(LIBS)

# Test programs, per mapper (they include its source).
test_located_hits_batmis1: $(TESTS)/test_located_hits.cpp batmis1.cpp $(HEADERS) $(FMINDEX_LIB) $(HTSLIB)
	$(COMPILE) -DMAPPER_SOURCE='"batmis1.cpp"' -o $@ $< $(LIBS)

test_located_hits_batmis2: $(TESTS)/test_located_hits.cpp batmis2.cpp $(HEADERS) $(FMINDEX_LIB) $(HTSLIB)
	$(COMPILE) -DMAPPER_SOURCE='"batmis2.cpp"' -o $@ $< $(LIBS)

test_sa_decode_batmis1: $(TESTS)/test_sa_decode.cpp batmis1.cpp $(HEADERS) $(FMINDEX_LIB) $(HTSLIB)
	$(COMPILE) -DMAPPER_SOURCE='"batmis1.cpp"' -o $@ $< $(LIBS)

test_sa_decode_batmis2: $(TESTS)/test_sa_decode.cpp batmis2.cpp $(HEADERS) $(FMINDEX_LIB) $(HTSLIB)
	$(COMPILE) -DMAPPER_SOURCE='"batmis2.cpp"' -o $@ $< $(LIBS)

test_genome_map: $(TESTS)/test_genome_map.cpp genome_map.hpp
	$(CXX) -std=c++17 -I. $(CXXFLAGS) -o $@ $<

TEST_VIEW = $(HTSLIB_DIR)/test/test_view

$(TEST_VIEW): FORCE
	$(MAKE) -C $(HTSLIB_DIR) test/test_view

# Each mapper's suite: located hits (lanes, walked against located, brute force for
# -n 1 to 7) and SAM/BAM output (brute force on both strands, Sa_Decode).
test-batmis1: batmis1 batmis1_chunk16 test_located_hits_batmis1 test_sa_decode_batmis1 test_genome_map $(INDEX_BUILDER) $(TEST_VIEW)
	MAPPER=batmis1 $(TESTS)/test_located_hits.sh . $(INDEX_BUILDER) ./test_located_hits_batmis1
	MAPPER=batmis1 $(TESTS)/test_sam_output.sh . $(INDEX_BUILDER) ./test_genome_map $(TEST_VIEW) ./test_sa_decode_batmis1

# batmis2's also compares it with batmis1 and runs batmis2_verify (test_batmis2.sh).
test-batmis2: batmis1 batmis2 batmis2_chunk16 batmis2_verify test_located_hits_batmis2 test_sa_decode_batmis2 test_genome_map \
              $(INDEX_BUILDER) $(TEST_VIEW)
	MAPPER=batmis2 $(TESTS)/test_located_hits.sh . $(INDEX_BUILDER) ./test_located_hits_batmis2
	MAPPER=batmis2 $(TESTS)/test_sam_output.sh . $(INDEX_BUILDER) ./test_genome_map $(TEST_VIEW) ./test_sa_decode_batmis2
	$(TESTS)/test_batmis2.sh . $(INDEX_BUILDER)

test: test-batmis1 test-batmis2

# index/ owns the library and the index builder; let its Makefile decide.
$(FMINDEX_LIB): FORCE
	$(MAKE) -C $(INDEX_DIR) libfmindex_core.a

$(INDEX_BUILDER): FORCE
	$(MAKE) -C $(INDEX_DIR) bwt_build_index

$(HTSLIB): FORCE
	$(MAKE) -C $(HTSLIB_DIR) lib-static

install: batmis1 batmis2 $(INDEX_BUILDER)
	mkdir -p $(BINDIR)
	cp batmis1 batmis2 $(INDEX_BUILDER) $(BINDIR)

clean:
	rm -f batmis1 batmis2 batmis2_verify batmis1_chunk16 batmis2_chunk16 \
	      test_located_hits_batmis1 test_located_hits_batmis2 test_sa_decode_batmis1 test_sa_decode_batmis2 \
	      test_genome_map

FORCE:
