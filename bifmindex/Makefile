# Bidirectional FM index over the BatMis index pair (see README.md).
#   make         bi_check, bi_exact, bi_mm and bi_batman
#   bi_batman: batman_zero_min ported to bidirectional SA ranges (bi_batman.cpp);
#   bi_batman_verify checks every range's other-index start (-DBI_VERIFY_OTHER)
#   make test    builds a toy index and checks the step and exact search (test.sh)

CXX ?= g++
CXXFLAGS ?= -O3 -funroll-loops -mpopcnt

INDEX_DIR ?= ../index
FMINDEX_LIB = $(INDEX_DIR)/libfmindex_core.a
SRC_DIR ?= ../src
HTSLIB_DIR ?= ../htslib-1.24
HTSLIB = $(HTSLIB_DIR)/libhts.a
HTSLIB_LIBS ?= -llzma -lbz2 -lz -lm -lpthread
BATMAN_DEPS = bi_batman.cpp $(SRC_DIR)/fm_occ.hpp $(SRC_DIR)/genome_map.hpp $(SRC_DIR)/sam_writer.hpp $(FMINDEX_LIB) $(HTSLIB)

.PHONY: all test clean FORCE

all: bi_check bi_exact bi_mm bi_batman

bi_check: bi_check.cpp bi_index.hpp $(FMINDEX_LIB)
	$(CXX) -std=c++17 -I$(INDEX_DIR) $(CXXFLAGS) -o $@ $< $(FMINDEX_LIB) -lm

bi_exact: bi_exact.cpp bi_index.hpp bi_ref.hpp bi_batch.hpp bi_locate.hpp $(FMINDEX_LIB)
	$(CXX) -std=c++17 -I$(INDEX_DIR) $(CXXFLAGS) -o $@ $< $(FMINDEX_LIB) -lm

bi_mm: bi_mm.cpp bi_index.hpp bi_ref.hpp bi_mismatch.hpp bi_two_mismatch.hpp bi_locate.hpp $(FMINDEX_LIB)
	$(CXX) -std=c++17 -I$(INDEX_DIR) $(CXXFLAGS) -o $@ $< $(FMINDEX_LIB) -lm

bi_batman: $(BATMAN_DEPS)
	$(CXX) -std=c++17 -fpermissive -I$(SRC_DIR) -I$(INDEX_DIR) -I$(HTSLIB_DIR) $(CXXFLAGS) -pthread -o $@ $< $(FMINDEX_LIB) $(HTSLIB) $(HTSLIB_LIBS)

bi_batman_verify: $(BATMAN_DEPS)
	$(CXX) -std=c++17 -fpermissive -I$(SRC_DIR) -I$(INDEX_DIR) -I$(HTSLIB_DIR) $(CXXFLAGS) -pthread -DBI_VERIFY_OTHER -o $@ $< $(FMINDEX_LIB) $(HTSLIB) $(HTSLIB_LIBS)

test: all
	$(MAKE) -C $(INDEX_DIR) bwt_build_index
	./test.sh . $(INDEX_DIR)/bwt_build_index

# index/ owns the library; let its Makefile decide whether it is up to date.
$(FMINDEX_LIB): FORCE
	$(MAKE) -C $(INDEX_DIR) libfmindex_core.a

clean:
	rm -f bi_check bi_exact bi_mm bi_batman bi_batman_verify

FORCE:
