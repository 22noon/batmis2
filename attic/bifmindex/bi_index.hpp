// bi_index.hpp: a bidirectional FM index over the BatMis index pair.
//
// bwt_build_index writes <g>.bwt/.fmv/.sa, the FM index of the text T (here Fw), and
// rev<g>.bwt/.fmv/.sa, the FM index of T reversed (Rev). A pattern P has a range in each:
// Fw's rows are T's suffixes that start with P, Rev's are reversed T's suffixes that start
// with P reversed. Both have one row per occurrence of P, so the same size, and a
// Bi_Range keeps them in step as P grows at either end:
//   - Extend_Left (P -> cP) is a backward step with c in Fw. Rev's rows of P are ordered
//     by the base after P reversed, which is the base before P in T: the Fw BWT
//     character of the matching row. So Rev's new range starts after the rows whose
//     character is smaller than c: Fw's $ row (P at T's start), if it is in the range,
//     then the A, C, ... rows of Fw's range.
//   - Extend_Right (P -> Pc) is the mirror: a backward step with c in Rev, Fw shifted by
//     Rev's counts.
//   - A single row (one occurrence) is two independent rows: a step moves only its own
//     index's row, by LF if the row's BWT character is the base added, and the other row
//     stays, as the rows before it in the general step are none (the singleton step).
//   - A narrow range counts its BWT characters directly (Count_Range): the counts give
//     the new range's size and, by the bases smaller than c, the other range's offset,
//     so one lookup at Start is enough (the narrow step).
// Conventions as index/BWT.h: rows 0..n (n = textLength), ranges inclusive, the $ row is
// inverseSa0 (its character is not stored), Occ(i, c) counts c in rows [0, i).
// See README.md.
#pragma once

#include <cstddef>  // size_t and FILE, which the index headers use without including
#include <cstdio>

extern "C" {
#include "MemManager.h"
#include "BWT.h"
}

#include <cstdint>
#include <cstring>
#include <emmintrin.h>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace bifm {

// -DBI_PROFILE keeps the three step kinds out of line, so a profile splits by kind.
#ifdef BI_PROFILE
#define BI_STEP __attribute__((noinline))
#else
#define BI_STEP
#endif

// Rows [Start, End] of one index; empty when End < Start.
struct Range {
    unsigned long Start = 0, End = 0;
    bool Empty() const { return End < Start; }
    unsigned long Size() const { return Empty() ? 0 : End - Start + 1; }
};

// The ranges of one pattern in Fw and in Rev (same size).
struct Bi_Range {
    unsigned long Fw_Start = 0, Fw_End = 0, Rev_Start = 0;
    bool Empty() const { return Fw_End < Fw_Start; }
    unsigned long Size() const { return Empty() ? 0 : Fw_End - Fw_Start + 1; }
    Range Fw() const { return {Fw_Start, Fw_End}; }
    Range Rev() const { return {Rev_Start, Rev_Start + Size() - 1}; }
};

// A hit's range chain in one index (Rev: Rev ranges of P'[0, k), else Fw ranges of
// P'[L - k, L)): levels Low .. L, level k the range after k bases (S[k], E[k]), as a
// search records them (bi_locate.hpp). P' is the read with base Sub_Base[s] at Sub_At[s]
// for s < Subs.
struct Range_Chain {
    bool Rev = false;
    unsigned Low = 1;
    std::vector<unsigned long> S, E;  // indexed by level; [0, Low) unset
    unsigned Subs = 0;
    unsigned Sub_At[2] = {0, 0};
    unsigned char Sub_Base[2] = {0, 0};
    void Substitute(unsigned At, unsigned char Base) { Sub_At[Subs] = At, Sub_Base[Subs++] = Base; }

    // Starts a chain of a read of L bases.
    void Reset(std::size_t L, bool In_Rev, unsigned First_Level) {
        Rev = In_Rev, Low = First_Level, Subs = 0;
        if (S.size() < L + 1) S.resize(L + 1), E.resize(L + 1);
    }
    // Records R as level Level, in the chain's index.
    void Note(std::size_t Level, const Bi_Range& R) {
        if (Rev) S[Level] = R.Rev_Start, E[Level] = R.Rev_Start + R.Size() - 1;
        else S[Level] = R.Fw_Start, E[Level] = R.Fw_End;
    }
};

// Base code (A C G T = 0 1 2 3, any case); -1 for anything else.
inline int Base_Code(char c) {
    switch (c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default: return -1;
    }
}

// The classic backward step in one index: P -> cP.
inline Range Backward_Step(const BWT* B, Range R, int c) {
    return {B->cumulativeFreq[c] + BWTOccValue(B, static_cast<unsigned>(R.Start), c) + 1,
            B->cumulativeFreq[c] + BWTOccValue(B, static_cast<unsigned>(R.End + 1), c)};
}

class Bi_Index {
public:
    // Loads <Prefix>.bwt/.fmv/.sa and rev<Prefix>.bwt/.fmv/.sa (the reverse index sits
    // next to the forward one, its file name prefixed with "rev").
    explicit Bi_Index(const std::string& Prefix) {
        const std::size_t Slash = Prefix.find_last_of('/');
        const std::string Dir = Slash == std::string::npos ? "" : Prefix.substr(0, Slash + 1);
        const std::string Base = Slash == std::string::npos ? Prefix : Prefix.substr(Slash + 1);
        const std::string Rev_Prefix = Dir + "rev" + Base;
        MMMasterInitialize(3, 0, FALSE, nullptr);
        Pool = MMPoolCreate(524288);
        Fw_Index = Load(Prefix);
        Rev_Index = Load(Rev_Prefix);
        if (Fw_Index->textLength != Rev_Index->textLength)
            throw std::runtime_error("the forward and reverse indexes are of different texts");
    }
    Bi_Index(const Bi_Index&) = delete;
    Bi_Index& operator=(const Bi_Index&) = delete;

    const BWT* Fw() const { return Fw_Index; }
    const BWT* Rev() const { return Rev_Index; }
    unsigned long Text_Length() const { return Fw_Index->textLength; }

    // The empty pattern: every row of both indexes.
    Bi_Range Empty_Pattern() const { return {0, Text_Length(), 0}; }

    bool Singleton = true;  // the singleton step for single rows (off: always the general step)
    // Ranges with End - Start below this take the narrow step (0: never), as
    // batman_zero_min's BRANCHTHRESHOLD.
    unsigned long Threshold = 128;

    // P -> cP.
    Bi_Range Extend_Left(const Bi_Range& R, int c) const { return Extend(Fw_Index, R.Fw(), R.Rev_Start, c, true); }
    // P -> Pc.
    Bi_Range Extend_Right(const Bi_Range& R, int c) const { return Extend(Rev_Index, R.Rev(), R.Fw_Start, c, false); }

    // All four extensions at once, Out[c] for base c (empty ones {1, 0, 0}), from one set
    // of counts: the four children of a node, for branching searches.
    void Children_Left(const Bi_Range& R, Bi_Range Out[4]) const { Children(Fw_Index, R.Fw(), R.Rev_Start, true, Out); }
    void Children_Right(const Bi_Range& R, Bi_Range Out[4]) const { Children(Rev_Index, R.Rev(), R.Fw_Start, false, Out); }

    // The start-of-search table: the Bi_Range of every k-mer (Build_Lookup), so a search
    // starts k bases in. One table serves both directions.
    int Lookup_Depth() const { return Depth; }
    void Build_Lookup(int k) {
        Depth = k;
        Lookup.assign(k ? 1ul << (2 * k) : 0, Lookup_Entry{1, 0, 0, 0});
        if (k) Build_Node(Empty_Pattern(), 0, 0);
    }
    // The Bi_Range of P[0, k), k = Lookup_Depth().
    Bi_Range Look_Up(const char* P) const {
        unsigned long Key = 0;
        for (int i = 0; i < Depth; i++) Key = Key << 2 | static_cast<unsigned char>(P[i]);
        const Lookup_Entry& E = Lookup[Key];
        return {E.Fw_Start, E.Fw_End, E.Rev_Start};
    }

    // The narrow step: as the general step, with Own's counts taken directly from its BWT
    // characters (Count_Range) instead of a second lookup at End + 1.
    BI_STEP static Bi_Range Extend_Narrow(const BWT* B, Range Own, unsigned long Other_Start, int c, bool Left) {
        unsigned Count[4];
        Count_Range(B, Own.Start, Own.End, Count);
        if (!Count[c]) return {1, 0, 0};  // empty
        unsigned long Before = B->inverseSa0 >= Own.Start && B->inverseSa0 <= Own.End;
        for (int b = 0; b < c; b++) Before += Count[b];
        const unsigned long Start = B->cumulativeFreq[c] + Occ_Base(B, Own.Start, c) + 1;
        if (Left) return {Start, Start + Count[c] - 1, Other_Start + Before};
        return {Other_Start + Before, Other_Start + Before + Count[c] - 1, Start};
    }

    // Counts each base among the BWT characters of rows Start..End of B, the $ row
    // skipped (batman_zero_min's Count_Range): three bit-planes over bwtCode, where a
    // character's low bit set means C or T, its high bit G or T, both T.
    static void Count_Range(const BWT* B, unsigned long Start, unsigned long End, unsigned Count[4]) {
        Count_Codes(B->bwtCode, Start - (Start > B->inverseSa0), End + 1 - (End + 1 > B->inverseSa0), Count);
    }

    // Counts each base among the stored characters [First, Last) of Code (character k in
    // bits 31-2(k%16)..30-2(k%16) of word k/16): the low bit of a pair set means C or T,
    // the high bit G or T, both T. Only the first and last words are masked; whole words in
    // between are counted two at a time.
    static void Count_Codes(const unsigned* Code, unsigned long First, unsigned long Last, unsigned Count[4]) {
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

    // Counts base c among the stored characters [First, Last) of Code: a pair equals c
    // where both its bits equal c's, i.e. where Y = ~(word ^ c's pattern) has both set.
    static unsigned Count_Code(const unsigned* Code, unsigned long First, unsigned long Last, unsigned c) {
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

    // Occ in B at Row (rows [0, Row), $ skipped) for every base, or for base c: the nearest
    // occurrence checkpoint (every OCC_INTERVAL characters, its value split into a major
    // count every OCC_INTERVAL_MAJOR and a 16-bit minor one, as BWT.c's
    // BWTAllOccValueExplicit), plus or minus a direct count of the at most
    // OCC_INTERVAL / 2 characters between it and Row. Inline, with no decode tables; equal
    // to BWTAllOccValue / BWTOccValue. Fused off: those library calls.
    static inline bool Fused = true;
    static void Occ_All(const BWT* B, unsigned long Row, unsigned Out[4]) {
        if (!Fused) return BWTAllOccValue(B, static_cast<unsigned>(Row), Out);
        const unsigned long I = Row - (Row > B->inverseSa0);
        const unsigned long E = (I + OCC_INTERVAL / 2 - 1) / OCC_INTERVAL, At = E * OCC_INTERVAL;
        const unsigned* Major = B->occValueMajor + E * OCC_INTERVAL / OCC_INTERVAL_MAJOR * ALPHABET_SIZE;
        const unsigned* Minor = B->occValue + E / OCC_VALUE_PER_WORD * ALPHABET_SIZE;
        const unsigned Shift = E % OCC_VALUE_PER_WORD ? 0 : 16;
        unsigned Count[4];
        if (At <= I) Count_Codes(B->bwtCode, At, I, Count);
        else Count_Codes(B->bwtCode, I, At, Count);
        for (int b = 0; b < 4; b++) {
            const unsigned Check = Major[b] + ((Minor[b] >> Shift) & 0xFFFFu);
            Out[b] = At <= I ? Check + Count[b] : Check - Count[b];
        }
    }
    // One base's Occ for the steps: Occ_One, or without Fused, BWTAllOccValue's count (in
    // a chain of LF steps about 40% faster than BWTOccValue; README.md, "Singleton step").
    static unsigned Occ_Base(const BWT* B, unsigned long Row, int c) {
        if (Fused) return Occ_One(B, Row, c);
        unsigned At[4];
        BWTAllOccValue(B, static_cast<unsigned>(Row), At);
        return At[c];
    }
    static unsigned Occ_One(const BWT* B, unsigned long Row, int c) {
        if (!Fused) return BWTOccValue(B, static_cast<unsigned>(Row), static_cast<unsigned>(c));
        const unsigned long I = Row - (Row > B->inverseSa0);
        const unsigned long E = (I + OCC_INTERVAL / 2 - 1) / OCC_INTERVAL, At = E * OCC_INTERVAL;
        const unsigned Shift = E % OCC_VALUE_PER_WORD ? 0 : 16;
        const unsigned Check = B->occValueMajor[E * OCC_INTERVAL / OCC_INTERVAL_MAJOR * ALPHABET_SIZE + c] +
                               ((B->occValue[E / OCC_VALUE_PER_WORD * ALPHABET_SIZE + c] >> Shift) & 0xFFFFu);
        if (At <= I) return Check + Count_Code(B->bwtCode, At, I, static_cast<unsigned>(c));
        return Check - Count_Code(B->bwtCode, I, At, static_cast<unsigned>(c));
    }

    // The general step: a backward step with c in index B (range Own), and the other
    // index's range Other_Start narrowed to its rows followed by c, the rows of Own whose
    // BWT character is smaller ($ first) coming before them. Left: B is Fw.
    BI_STEP static Bi_Range Extend_General(const BWT* B, Range Own, unsigned long Other_Start, int c, bool Left) {
        unsigned At_Start[4], At_End[4];
        Occ_All(B, Own.Start, At_Start);
        Occ_All(B, Own.End + 1, At_End);
        unsigned long Before = B->inverseSa0 >= Own.Start && B->inverseSa0 <= Own.End;
        for (int b = 0; b < c; b++) Before += At_End[b] - At_Start[b];
        const unsigned long Start = B->cumulativeFreq[c] + At_Start[c] + 1, End = B->cumulativeFreq[c] + At_End[c];
        if (End < Start) return {1, 0, 0};  // empty
        if (Left) return {Start, End, Other_Start + Before};
        return {Other_Start + Before, Other_Start + Before + (End + 1 - Start) - 1, Start};
    }

    // Prefetches the occurrence and BWT words an Occ lookup at Row reads
    // (batman_zero_min's Prefetch_Row).
    static void Prefetch_Row(const BWT* B, unsigned long Row) {
        Row -= (Row > B->inverseSa0);
        const unsigned long Explicit = (Row + OCC_INTERVAL / 2 - 1) / OCC_INTERVAL;
        __builtin_prefetch(B->occValue + Explicit / OCC_VALUE_PER_WORD * ALPHABET_SIZE);
        __builtin_prefetch(B->occValueMajor + Explicit * OCC_INTERVAL / OCC_INTERVAL_MAJOR * ALPHABET_SIZE);
        __builtin_prefetch(B->bwtCode + Row / CHAR_PER_WORD);
        __builtin_prefetch(B->bwtCode + Explicit * OCC_INTERVAL / CHAR_PER_WORD);
    }
    // Prefetches what the next step from R reads: Extend_Left (Left) or Extend_Right.
    void Prefetch_Step(const Bi_Range& R, bool Left) const {
        const BWT* B = Left ? Fw_Index : Rev_Index;
        const Range Own = Left ? R.Fw() : R.Rev();
        Prefetch_Row(B, Own.Start);
        if (Own.Start == Own.End) return;
        if (Own.End - Own.Start < Threshold) __builtin_prefetch(B->bwtCode + (Own.End - (Own.End > B->inverseSa0)) / CHAR_PER_WORD);
        else Prefetch_Row(B, Own.End + 1);
    }
    // Prefetches P[0, k)'s lookup entry.
    void Prefetch_Look_Up(const char* P) const {
        unsigned long Key = 0;
        for (int i = 0; i < Depth; i++) Key = Key << 2 | static_cast<unsigned char>(P[i]);
        __builtin_prefetch(Lookup.data() + Key);
    }

    // The BWT character of Row in B (not the $ row).
    static unsigned Bwt_Char(const BWT* B, unsigned long Row) {
        const unsigned long I = Row - (Row > B->inverseSa0);
        return B->bwtCode[I / 16] << (2 * (I % 16)) >> 30;
    }

    // Exact search of P (base codes) by Extend_Left from its last base, or by
    // Extend_Right from its first; stops once the range is empty.
    // With a lookup table, from the range of P's last (Left) or first k bases.
    // With Ch, each step's range is recorded (Fw going left, Rev going right).
    Bi_Range Search_Left(std::string_view P, Range_Chain* Ch = nullptr) const {
        Bi_Range R = Empty_Pattern();
        const std::size_t L = P.size();
        std::size_t i = L;
        if (Depth && L >= static_cast<std::size_t>(Depth)) {
            i -= Depth;
            R = Look_Up(P.data() + i);
        }
        if (Ch) Ch->Reset(L, false, i == L ? 1 : static_cast<unsigned>(L - i));
        if (Ch && i < L && !R.Empty()) Ch->Note(L - i, R);
        while (i-- > 0 && !R.Empty()) {
            R = Extend_Left(R, P[i]);
            if (Ch && !R.Empty()) Ch->Note(L - i, R);
        }
        return R;
    }
    Bi_Range Search_Right(std::string_view P, Range_Chain* Ch = nullptr) const {
        Bi_Range R = Empty_Pattern();
        std::size_t i = 0;
        if (Depth && P.size() >= static_cast<std::size_t>(Depth)) {
            R = Look_Up(P.data());
            i = Depth;
        }
        if (Ch) Ch->Reset(P.size(), true, i ? static_cast<unsigned>(i) : 1);
        if (Ch && i && !R.Empty()) Ch->Note(i, R);
        for (; i < P.size() && !R.Empty(); i++) {
            R = Extend_Right(R, P[i]);
            if (Ch && !R.Empty()) Ch->Note(i + 1, R);
        }
        return R;
    }

    // Classic searches in one index: P's Fw range by backward steps from its last base;
    // its Rev range by backward steps in Rev from its first base (the "forward search").
    Range Backward_Search(std::string_view P) const {
        Range R{0, Text_Length()};
        for (std::size_t i = P.size(); i-- > 0 && !R.Empty();) R = Backward_Step(Fw_Index, R, P[i]);
        return R;
    }
    Range Forward_Search(std::string_view P) const {
        Range R{0, Text_Length()};
        for (std::size_t i = 0; i < P.size() && !R.Empty(); i++) R = Backward_Step(Rev_Index, R, P[i]);
        return R;
    }

    // Where in T the occurrence at Fw row Row starts; and at Rev row Row, for a pattern of
    // Length bases (reversed T's position q is T's n - q - Length).
    unsigned long Locate_Fw(unsigned long Row) const { return BWTSaValue(Fw_Index, static_cast<unsigned>(Row)); }
    unsigned long Locate_Rev(unsigned long Row, unsigned long Length) const {
        return Text_Length() - BWTSaValue(Rev_Index, static_cast<unsigned>(Row)) - Length;
    }

private:
    MMPool* Pool = nullptr;
    BWT* Fw_Index = nullptr;
    BWT* Rev_Index = nullptr;

    // 16 bytes, so a lookup reads one cache line. Rows fit 32 bits (BWT.h's indexes are
    // unsigned).
    struct Lookup_Entry {
        uint32_t Fw_Start, Fw_End, Rev_Start, Unused;
    };
    std::vector<Lookup_Entry> Lookup;
    int Depth = 0;

    // Fills the entries under R, the range of the k-mers' last Level bases, whose key bits
    // so far are Key (first base in the high bits): the children by one count each.
    void Build_Node(const Bi_Range& R, int Level, unsigned long Key) {
        if (Level == Depth) {
            Lookup[Key] = {static_cast<uint32_t>(R.Fw_Start), static_cast<uint32_t>(R.Fw_End),
                           static_cast<uint32_t>(R.Rev_Start), 0};
            return;
        }
        Bi_Range Child[4];
        Children_Left(R, Child);
        for (int c = 0; c < 4; c++)
            if (!Child[c].Empty()) Build_Node(Child[c], Level + 1, Key | static_cast<unsigned long>(c) << (2 * Level));
    }

    // Children_Left / Children_Right in B (Fw for Left): the four counts of Own's BWT
    // characters (from the singleton row, Count_Range, or lookups at both ends), each
    // base's range starting after the smaller ones ($ first) in the other index.
    void Children(const BWT* B, Range Own, unsigned long Other_Start, bool Left, Bi_Range Out[4]) const {
        for (int c = 0; c < 4; c++) Out[c] = {1, 0, 0};
        unsigned At_Start[4], Count[4];
        if (Singleton && Own.Start == Own.End) {
            if (Own.Start == B->inverseSa0) return;
            const unsigned c = Bwt_Char(B, Own.Start);
            Out[c] = Singleton_Step(B, Own.Start, Other_Start, static_cast<int>(c), Left);
            return;
        }
        Occ_All(B, Own.Start, At_Start);
        __m128i N;  // the four counts
        const __m128i S = _mm_loadu_si128(reinterpret_cast<const __m128i*>(At_Start));
        if (Own.End - Own.Start < Threshold) {
            Count_Range(B, Own.Start, Own.End, Count);
            N = _mm_loadu_si128(reinterpret_cast<const __m128i*>(Count));
        } else {
            unsigned At_End[4];
            Occ_All(B, Own.End + 1, At_End);
            N = _mm_sub_epi32(_mm_loadu_si128(reinterpret_cast<const __m128i*>(At_End)), S);
            _mm_storeu_si128(reinterpret_cast<__m128i*>(Count), N);
        }
        // Each base's start in the other index: Other_Start, [$ in range], and the counts of
        // the smaller bases, an exclusive prefix sum: [0, n0, n0+n1, n0+n1+n2].
        __m128i P = _mm_slli_si128(N, 4);
        P = _mm_add_epi32(P, _mm_slli_si128(P, 4));
        P = _mm_add_epi32(P, _mm_slli_si128(P, 8));
        const unsigned Base = static_cast<unsigned>(Other_Start + (B->inverseSa0 >= Own.Start && B->inverseSa0 <= Own.End));
        P = _mm_add_epi32(P, _mm_set1_epi32(static_cast<int>(Base)));
        // Each base's start in this index: C[c] + Occ(Start, c) + 1.
        const __m128i Q = _mm_add_epi32(
            _mm_add_epi32(_mm_loadu_si128(reinterpret_cast<const __m128i*>(B->cumulativeFreq)), S), _mm_set1_epi32(1));
        unsigned Other[4], Own_Start[4];
        _mm_storeu_si128(reinterpret_cast<__m128i*>(Other), P);
        _mm_storeu_si128(reinterpret_cast<__m128i*>(Own_Start), Q);
        for (int c = 0; c < 4; c++) {
            if (!Count[c]) continue;
            Out[c] = Left ? Bi_Range{Own_Start[c], Own_Start[c] + Count[c] - 1, Other[c]}
                          : Bi_Range{Other[c], Other[c] + Count[c] - 1, Own_Start[c]};
        }
    }

    BWT* Load(const std::string& Prefix) const {
        const std::string Bwt = Prefix + ".bwt", Occ = Prefix + ".fmv", Sa = Prefix + ".sa";
        BWT* B = BWTLoad(Pool, Bwt.c_str(), Occ.c_str(), Sa.c_str(), nullptr, nullptr, nullptr);
        if (!B) throw std::runtime_error("cannot load the index " + Prefix);
        return B;
    }

    // A step with c in B (Fw for Left), from Own (B's range) and the other index's
    // Other_Start: the singleton, narrow or general step by Own's width.
    Bi_Range Extend(const BWT* B, Range Own, unsigned long Other_Start, int c, bool Left) const {
        if (Singleton && Own.Start == Own.End) return Singleton_Step(B, Own.Start, Other_Start, c, Left);
        if (Own.End - Own.Start < Threshold) return Extend_Narrow(B, Own, Other_Start, c, Left);
        return Extend_General(B, Own, Other_Start, c, Left);
    }

    // The singleton step: Own row Row of B, the other index's row Other. The occurrence
    // extends with c only if Row's BWT character is c (the $ row extends with nothing);
    // then Row moves by LF and Other stays: in the general step, no row of [Row, Row] is
    // $ or smaller than c. Left: B is Fw.
    BI_STEP static Bi_Range Singleton_Step(const BWT* B, unsigned long Row, unsigned long Other, int c, bool Left) {
        if (Row == B->inverseSa0 || Bwt_Char(B, Row) != static_cast<unsigned>(c)) return {1, 0, 0};  // empty
        const unsigned long Next = B->cumulativeFreq[c] + Occ_Base(B, Row, c) + 1;
        if (Left) return {Next, Next, Other};
        return {Other, Other, Next};
    }
};

}  // namespace bifm
