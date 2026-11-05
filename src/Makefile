# Builds the batman_zero_min exact-match searchers.
#   make                     batman_zero_min.cpp, native BatMis FM index (../index)
#   make batman_zero_min_bwa batman_zero_min_bwa.cpp, BWA index (needs a built bwa-fmd)
# Native build variant via DEFS: DEFS=-DBATMAN_ZERO_MIN_USE_PAIRED_RANK
#   make batman_zero_min_inlineocc  fm_occ.hpp's inline Occ instead of the library's
#                                   (BWTOccValue, BWTAllOccValue), for comparison
#   make test-occ                fm_occ.hpp's Occ against the library's

CXX ?= g++
CXXFLAGS ?= -O3 -funroll-loops -mpopcnt -pthread
DEFS ?=
BINDIR ?= ../bin

INDEX_DIR ?= ../index
FMINDEX_LIB = $(INDEX_DIR)/libfmindex_core.a
# htslib (SAM/BAM output), built from source; its configure picked these libraries.
HTSLIB_DIR ?= ../htslib-1.24
HTSLIB = $(HTSLIB_DIR)/libhts.a
HTSLIB_LIBS ?= -llzma -lbz2 -lz -lm -lpthread
BWA_ROOT ?= ../../../bwa-fmd

.PHONY: all test-occ batman_zero_min_fmindex batman_zero_structured test_exact_search_intervals test-structured test-located test-sam bench_sa_decode install clean FORCE

all: batman_zero_min

batman_zero_min.o: batman_zero_min.cpp fm_occ.hpp genome_map.hpp sam_writer.hpp
	$(CXX) -std=c++17 -fpermissive -I$(INDEX_DIR) -I$(HTSLIB_DIR) $(CXXFLAGS) $(DEFS) -c $< -o $@

batman_zero_min: batman_zero_min.o $(FMINDEX_LIB) $(HTSLIB)
	$(CXX) -std=c++17 $(CXXFLAGS) -o $@ batman_zero_min.o $(FMINDEX_LIB) $(HTSLIB) $(HTSLIB_LIBS)

# The same with 16-read chunks, so that the tests' few reads span several blocks,
# chunks and threads (-t).
batman_zero_min_chunk16: batman_zero_min.cpp fm_occ.hpp genome_map.hpp sam_writer.hpp $(FMINDEX_LIB) $(HTSLIB)
	$(CXX) -std=c++17 -fpermissive -I$(INDEX_DIR) -I$(HTSLIB_DIR) $(CXXFLAGS) $(DEFS) -DBATMAN_BATCH_BLOCK=16 -o $@ $< \
	    $(FMINDEX_LIB) $(HTSLIB) $(HTSLIB_LIBS)

batman_zero_min_inlineocc: batman_zero_min.cpp fm_occ.hpp genome_map.hpp sam_writer.hpp $(FMINDEX_LIB) $(HTSLIB)
	$(CXX) -std=c++17 -fpermissive -I$(INDEX_DIR) -I$(HTSLIB_DIR) $(CXXFLAGS) $(DEFS) -DBATMAN_INLINE_OCC -o $@ $< \
	    $(FMINDEX_LIB) $(HTSLIB) $(HTSLIB_LIBS)

test_occ_inline: ../tests/test_occ_inline.cpp fm_occ.hpp $(FMINDEX_LIB)
	$(CXX) -std=c++17 -I. -I$(INDEX_DIR) $(CXXFLAGS) -DBATMAN_INLINE_OCC -o $@ $< $(FMINDEX_LIB) -lm

test-occ: test_occ_inline
	$(MAKE) -C $(INDEX_DIR) bwt_build_index
	../tests/test_occ_inline.sh $(INDEX_DIR)/bwt_build_index ./test_occ_inline

# Old name for batman_zero_min
batman_zero_min_fmindex: batman_zero_min

batman_zero_structured.o: batman_zero_structured.cpp exact_search.hpp
	$(CXX) -std=c++17 -fpermissive -I$(INDEX_DIR) $(CXXFLAGS) -c $< -o $@

exact_search.o: exact_search.cpp exact_search.hpp
	$(CXX) -std=c++17 -fpermissive -I$(INDEX_DIR) $(CXXFLAGS) -c $< -o $@

batman_zero_structured: batman_zero_structured.o exact_search.o $(FMINDEX_LIB)
	$(CXX) -std=c++17 $(CXXFLAGS) -o $@ batman_zero_structured.o exact_search.o $(FMINDEX_LIB) -lm

test_exact_search_intervals.o: ../tests/test_exact_search_intervals.cpp exact_search.hpp
	$(CXX) -std=c++17 -fpermissive -I. -I$(INDEX_DIR) $(CXXFLAGS) -c $< -o $@

test_exact_search_intervals: test_exact_search_intervals.o exact_search.o $(FMINDEX_LIB)
	$(CXX) -std=c++17 $(CXXFLAGS) -o $@ test_exact_search_intervals.o exact_search.o $(FMINDEX_LIB) -lm

test-structured: batman_zero_min batman_zero_structured test_exact_search_intervals
	$(MAKE) -C $(INDEX_DIR) bwt_build_index
	../tests/test_structured_exact_search.sh . $(INDEX_DIR)/bwt_build_index ./test_exact_search_intervals

test_located_hits: ../tests/test_located_hits.cpp batman_zero_min.cpp fm_occ.hpp genome_map.hpp sam_writer.hpp $(FMINDEX_LIB) $(HTSLIB)
	$(CXX) -std=c++17 -fpermissive -I. -I$(INDEX_DIR) -I$(HTSLIB_DIR) $(CXXFLAGS) $(DEFS) -o $@ $< $(FMINDEX_LIB) $(HTSLIB) $(HTSLIB_LIBS)

test-located: batman_zero_min batman_zero_min_chunk16 test_located_hits
	$(MAKE) -C $(INDEX_DIR) bwt_build_index
	../tests/test_located_hits.sh . $(INDEX_DIR)/bwt_build_index ./test_located_hits

test_genome_map: ../tests/test_genome_map.cpp genome_map.hpp
	$(CXX) -std=c++17 -I. $(CXXFLAGS) -o $@ $<

test_sa_decode: ../tests/test_sa_decode.cpp batman_zero_min.cpp fm_occ.hpp genome_map.hpp sam_writer.hpp $(FMINDEX_LIB) $(HTSLIB)
	$(CXX) -std=c++17 -fpermissive -I. -I$(INDEX_DIR) -I$(HTSLIB_DIR) $(CXXFLAGS) $(DEFS) -o $@ $< $(FMINDEX_LIB) $(HTSLIB) $(HTSLIB_LIBS)

# Replays recorded decode requests (tests/bench_sa_decode.cpp); not part of the tests.
bench_sa_decode: ../tests/bench_sa_decode.cpp batman_zero_min.cpp fm_occ.hpp genome_map.hpp sam_writer.hpp $(FMINDEX_LIB) $(HTSLIB)
	$(CXX) -std=c++17 -fpermissive -I. -I$(INDEX_DIR) -I$(HTSLIB_DIR) $(CXXFLAGS) $(DEFS) -DBATMAN_BATCH_STATS -o $@ $< $(FMINDEX_LIB) $(HTSLIB) $(HTSLIB_LIBS)

test_chain_decode: ../tests/test_chain_decode.cpp batman_zero_min.cpp fm_occ.hpp genome_map.hpp sam_writer.hpp $(FMINDEX_LIB) $(HTSLIB)
	$(CXX) -std=c++17 -fpermissive -I. -I$(INDEX_DIR) -I$(HTSLIB_DIR) $(CXXFLAGS) $(DEFS) -DBATMAN_BATCH_STATS -o $@ $< $(FMINDEX_LIB) $(HTSLIB) $(HTSLIB_LIBS)

test-sam: batman_zero_min batman_zero_min_chunk16 test_genome_map test_sa_decode test_chain_decode
	$(MAKE) -C $(INDEX_DIR) bwt_build_index
	$(MAKE) -C $(HTSLIB_DIR) test/test_view
	../tests/test_sam_output.sh . $(INDEX_DIR)/bwt_build_index ./test_genome_map $(HTSLIB_DIR)/test/test_view \
	    ./test_sa_decode ./test_chain_decode

# index/ owns the library; let its Makefile decide whether it is up to date.
$(FMINDEX_LIB): FORCE
	$(MAKE) -C $(INDEX_DIR) libfmindex_core.a

$(HTSLIB): FORCE
	$(MAKE) -C $(HTSLIB_DIR) lib-static

batman_zero_min_bwa.o: batman_zero_min_bwa.cpp
	$(CXX) -std=c++17 -fpermissive -I$(BWA_ROOT) $(CXXFLAGS) -c $< -o $@

batman_zero_min_bwa: batman_zero_min_bwa.o
	$(CXX) -std=c++17 $(CXXFLAGS) -o $@ $< $(BWA_ROOT)/libbwa.a -lz -lpthread -lm

install: batman_zero_min
	mkdir -p $(BINDIR)
	cp batman_zero_min $(BINDIR)
	if [ -f batman_zero_min_bwa ]; then cp batman_zero_min_bwa $(BINDIR); fi
	if [ -f $(INDEX_DIR)/bwt_build_index ]; then cp $(INDEX_DIR)/bwt_build_index $(BINDIR); fi

clean:
	rm -f batman_zero_min.o batman_zero_min_bwa.o batman_zero_min batman_zero_min_bwa \
	      batman_zero_min_spot.o batman_zero_min_spot batman_zero_structured.o \
	      exact_search.o batman_zero_structured test_exact_search_intervals.o \
	      test_exact_search_intervals test_located_hits test_genome_map test_sa_decode test_chain_decode bench_sa_decode \
	      batman_zero_min_chunk16 batman_zero_min_inlineocc test_occ_inline

FORCE:
