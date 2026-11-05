#ifndef BATMIS_EXACT_SEARCH_HPP
#define BATMIS_EXACT_SEARCH_HPP

#include <array>
#include <cstddef>
#include <cstdio>
#include <cstdint>

extern "C" {
#include "BWT.h"
}

namespace batmis {

// Search results stay in the BatMIS BWT's inclusive SA-row convention.
// Skip preserves the minimized mapper's sampled-SA singleton semantics.
struct ExactResult {
    unsigned long Start = 0;
    unsigned long End = 0;
    unsigned char Skip = 0;
};

// Per-index/per-read-length exact-search state. The table is built once after
// the first read establishes the fixed seed length; search() allocates nothing.
class ExactSearch {
public:
    explicit ExactSearch(const BWT* index);
    bool prepare(int readLength);
    bool search(const unsigned char* encodedRead, ExactResult& result) const;

private:
    struct SeedRange {
        unsigned Start;
        unsigned End;
        unsigned Label;
    };

    struct SearchRange {
        unsigned long Start;
        unsigned long End;
        int Level;
        unsigned char Skip;
    };

    static constexpr unsigned kLookupCapacity = 4096;
    const BWT* index_;
    int readLength_ = 0;
    int lookupSize_ = 0;
    std::array<unsigned, kLookupCapacity> seedStart_{};
    std::array<unsigned, kLookupCapacity> seedEnd_{};

    void extend(unsigned long& start, unsigned long& end, unsigned char base) const;
    SeedRange initialBaseRange(int base) const;
    void buildSeed(SeedRange range, int level, int base);
    bool searchExact(const unsigned char* read, SearchRange& range) const;
};

}  // namespace batmis

#endif
