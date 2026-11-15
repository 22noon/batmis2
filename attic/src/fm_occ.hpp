// fm_occ.hpp: occurrence counts (Occ) of a BatMis FM index, inline and without decode
// tables (benchmark/inline_occ.md).
//
// Occ(Row, c) counts c among the BWT characters of rows [0, Row), the $ row skipped. The
// index stores a checkpoint every OCC_INTERVAL characters (a major count every
// OCC_INTERVAL_MAJOR plus a 16-bit minor one), and the library (BWT.c's BWTOccValue /
// BWTAllOccValue) adds or subtracts the characters between the nearest checkpoint and Row
// through decode tables. With -DBATMAN_INLINE_OCC, those at most OCC_INTERVAL / 2
// characters are counted directly from bwtCode, with popcounts; the values are the
// library's (tests/test_occ_inline.cpp). By default Occ_All / Occ_One are the library
// calls: here the inline count saves 6-9% of instructions but no cycles, and -n 1 runs
// 10-15% slower (benchmark/inline_occ.md). bifmindex/ uses it, where it does pay.
#pragma once

#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstring>

extern "C" {
#include "MemManager.h"
#include "BWT.h"
}

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

// Counts base c among the stored characters [First, Last) of Code: a pair equals c where
// both its bits equal c's, i.e. where Y = ~(word ^ c's pattern) has both set.
inline unsigned Count_Code(const unsigned* Code, unsigned long First, unsigned long Last, unsigned c) {
    if (First >= Last) return 0;
    const unsigned long W0 = First / 16, W1 = (Last - 1) / 16;
    const unsigned Head = ~0u >> (2 * (First % 16));
    const unsigned Tail = ~0u << (30 - 2 * ((Last - 1) % 16));
    const unsigned Pattern = c * 0x55555555u;
    const uint64_t Pattern64 = c * 0x5555555555555555ULL;
    unsigned long N = 0;
    auto Word32 = [&](unsigned W, unsigned M) {
        const unsigned Y = ~(W ^ Pattern);
        N += __builtin_popcount(Y & (Y >> 1) & M & 0x55555555u);
    };
    if (W0 == W1) {
        Word32(Code[W0], Head & Tail);
    } else {
        Word32(Code[W0], Head);
        unsigned long W = W0 + 1;
        for (; W + 2 <= W1; W += 2) {
            uint64_t X;
            std::memcpy(&X, Code + W, 8);
            const uint64_t Y = ~(X ^ Pattern64);
            N += __builtin_popcountll(Y & (Y >> 1) & 0x5555555555555555ULL);
        }
        if (W < W1) Word32(Code[W], ~0u);
        Word32(Code[W1], Tail);
    }
    return static_cast<unsigned>(N);
}

// Occ at Row for every base.
inline void Occ_All(const BWT* B, unsigned long Row, unsigned Out[4]) {
#ifndef BATMAN_INLINE_OCC
    BWTAllOccValue(B, static_cast<unsigned>(Row), Out);
#else
    const unsigned long I = Row - (Row > B->inverseSa0);
    const unsigned long E = (I + OCC_INTERVAL / 2 - 1) / OCC_INTERVAL, At = E * OCC_INTERVAL;
    const unsigned* Major = B->occValueMajor + E * OCC_INTERVAL / OCC_INTERVAL_MAJOR * ALPHABET_SIZE;
    const unsigned* Minor = B->occValue + E / OCC_VALUE_PER_WORD * ALPHABET_SIZE;
    const unsigned Shift = E % OCC_VALUE_PER_WORD ? 0 : 16;  // two minor counts a word, the even one high
    unsigned Count[4];
    if (At <= I) Count_Codes(B->bwtCode, At, I, Count);
    else Count_Codes(B->bwtCode, I, At, Count);
    for (int b = 0; b < 4; b++) {
        const unsigned Check = Major[b] + ((Minor[b] >> Shift) & 0xFFFFu);
        Out[b] = At <= I ? Check + Count[b] : Check - Count[b];
    }
#endif
}

// Occ at Row for base c.
inline unsigned Occ_One(const BWT* B, unsigned long Row, unsigned c) {
#ifndef BATMAN_INLINE_OCC
    return BWTOccValue(B, static_cast<unsigned>(Row), c);
#else
    const unsigned long I = Row - (Row > B->inverseSa0);
    const unsigned long E = (I + OCC_INTERVAL / 2 - 1) / OCC_INTERVAL, At = E * OCC_INTERVAL;
    const unsigned Shift = E % OCC_VALUE_PER_WORD ? 0 : 16;
    const unsigned Check = B->occValueMajor[E * OCC_INTERVAL / OCC_INTERVAL_MAJOR * ALPHABET_SIZE + c] +
                           ((B->occValue[E / OCC_VALUE_PER_WORD * ALPHABET_SIZE + c] >> Shift) & 0xFFFFu);
    if (At <= I) return Check + Count_Code(B->bwtCode, At, I, c);
    return Check - Count_Code(B->bwtCode, I, At, c);
#endif
}

}  // namespace fm_occ
