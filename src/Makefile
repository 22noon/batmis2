# Builds the batman_zero_min exact-match searchers.
#   make                     batman_zero_min.cpp, native BatMis FM index (../index)
#   make batman_zero_min_bwa batman_zero_min_bwa.cpp, BWA index (needs a built bwa-fmd)
# Native build variant via DEFS: DEFS=-DBATMAN_ZERO_MIN_USE_PAIRED_RANK

CXX ?= g++
CXXFLAGS ?= -O3 -funroll-loops -mpopcnt
DEFS ?=
BINDIR ?= ../bin

INDEX_DIR ?= ../index
FMINDEX_LIB = $(INDEX_DIR)/libfmindex_core.a
BWA_ROOT ?= ../../../bwa-fmd

.PHONY: all batman_zero_min_fmindex batman_zero_structured test_exact_search_intervals test-structured install clean FORCE

all: batman_zero_min

batman_zero_min.o: batman_zero_min.cpp
	$(CXX) -std=c++17 -fpermissive -I$(INDEX_DIR) $(CXXFLAGS) $(DEFS) -c $< -o $@

batman_zero_min: batman_zero_min.o $(FMINDEX_LIB)
	$(CXX) -std=c++17 $(CXXFLAGS) -o $@ batman_zero_min.o $(FMINDEX_LIB) -lm

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

# index/ owns the library; let its Makefile decide whether it is up to date.
$(FMINDEX_LIB): FORCE
	$(MAKE) -C $(INDEX_DIR) libfmindex_core.a

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
	      test_exact_search_intervals

FORCE:
