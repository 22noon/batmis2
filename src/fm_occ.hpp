// fm_occ.hpp: occurrence counts (Occ) of a BatMis FM index, for batmis1 and batmis2.
//
// Occ(Row, c) counts c among the BWT characters of rows [0, Row), the $ row skipped:
// Occ_All / Occ_One are the library's (BWT.c's BWTAllOccValue / BWTOccValue).
// Count_Codes counts each base among stored characters directly, for ranges narrower
// than a second Occ lookup is worth (Count_Range; benchmark/branch_threshold.md).
// An inline Occ (counting from the nearest checkpoint with popcounts instead of decode
// tables) was tried and saved instructions but no time here (benchmark/inline_occ.md);
// it is in attic/src/fm_occ.hpp.
#pragma once

#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstring>

extern "C" {
#include "MemManager.h"
#include "BWT.h"
}

// For a function that only prefetches. GCC knows of no side effect in __builtin_prefetch,
// so a call to such a function that it does not inline is removed as dead code, silently
// (GCC 13.3 -O3 dropped every Prefetch_Row and Mm_Prefetch that way; see
// benchmark/interleaved_index.md). Forced inline, the prefetches land in a caller that
// has other effects and stay. tests/check_prefetch.sh checks the binaries.
#define PREFETCH_ONLY __attribute__((always_inline)) inline

namespace fm_occ {

// Counts each base among the stored characters [First, Last) of Code (character k in
// bits 31-2(k%16)..30-2(k%16) of word k/16). Three bit-planes instead of four base
// comparisons: the low bit of a pair set means C or T, the high bit G or T, both T; A is
// what remains. Only the first and last words are masked; whole words in between are
// counted two at a time.
inline void Count_Codes(const unsigned* Code, unsigned long First, unsigned long Last, unsigned Count[4]) {
    if (First >= Last) {
        Count[0] = Count[1] = Count[2] = Count[3] = 0;
        return;
    }
    const unsigned long W0 = First / 16, W1 = (Last - 1) / 16;
    const unsigned Head = ~0u >> (2 * (First % 16));           // characters k >= First % 16
    const unsigned Tail = ~0u << (30 - 2 * ((Last - 1) % 16));  // characters k <= (Last - 1) % 16
    unsigned long CT = 0, GT = 0, T = 0;
    auto Word32 = [&](unsigned W, unsigned M) {
        M &= 0x55555555u;
        const unsigned L = W & M, H = (W >> 1) & M;
        CT += __builtin_popcount(L);
        GT += __builtin_popcount(H);
        T += __builtin_popcount(L & H);
    };
    if (W0 == W1) {
        Word32(Code[W0], Head & Tail);
    } else {
        Word32(Code[W0], Head);
        unsigned long W = W0 + 1;
        for (; W + 2 <= W1; W += 2) {
            uint64_t X;
            std::memcpy(&X, Code + W, 8);
            // (X >> 1) moves bit 32 into bit 31, an odd bit that the mask drops.
            const uint64_t L = X & 0x5555555555555555ULL, H = (X >> 1) & 0x5555555555555555ULL;
            CT += __builtin_popcountll(L);
            GT += __builtin_popcountll(H);
            T += __builtin_popcountll(L & H);
        }
        if (W < W1) Word32(Code[W], ~0u);
        Word32(Code[W1], Tail);
    }
    const unsigned long N = Last - First;
    Count[3] = T;
    Count[1] = CT - T;
    Count[2] = GT - T;
    Count[0] = N - CT - GT + T;
}

// Occ at Row for every base.
inline void Occ_All(const BWT* B, unsigned long Row, unsigned Out[4]) {
    BWTAllOccValue(B, static_cast<unsigned>(Row), Out);
}

// Occ at Row for base c.
inline unsigned Occ_One(const BWT* B, unsigned long Row, unsigned c) {
    return BWTOccValue(B, static_cast<unsigned>(Row), c);
}

}  // namespace fm_occ
