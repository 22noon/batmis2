// batmis2: batmis1 (batmis1.cpp) with bidirectional SA ranges.
// Every SARange also keeps Other_Start, the start of the same string's range in the other
// index (revfmi <-> fwfmi), so a range is one in both: a step with c moves it past the
// rows of the range whose BWT character is $ or below c (Get_SARange, Branch_Detect); a
// single-row step leaves it; a 1-base range starts at the same row in both indexes; an
// index switch ends at the seed's own Start. Two uses (Switch_Walk, Both_Rows): an index
// switch reads the switched range off the seed instead of re-walking the seed in the
// other index, and a single row is located from either index's sampled row. The search
// is otherwise batmis1's, with the same hits (byte-identical output with --switch-walk
// --own-row). -DBATMIS2_VERIFY checks Other_Start and every other-row locate against
// the other index. What follows is batmis1's description.
//
// batmis1: the BatMis short-read mapper on a pair of one-direction FM indexes, the
// reverse index rev<genome>.bwt/.fmv and the forward index <genome>.bwt/.fmv (built by
// index/bwt_build_index): exact search, and with -n 1 to 7 hits with up to that many
// mismatches.
// By default it also loads the SA samples (.sa) and the packed genome (.pac):
// a single-row path that reaches a sampled row is finished by comparing the read
// with the genome there (see Locate). --no-ref skips them to save memory.
#include <algorithm>
#include <array>
#include <cassert>
#include <cctype>
#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <string_view>
#include <vector>
#include <atomic>
#include <chrono>
#include <condition_variable>
#include <functional>
#include <mutex>
#include <thread>

extern "C" {
#include "MemManager.h"
#include "BWT.h"
}

#include "fm_occ.hpp"
#include "genome_map.hpp"
#include "sam_writer.hpp"

#include <zlib.h>

namespace {

constexpr int MAXTAG = 4096;
constexpr int MAX_MISMATCHES_BOUND = 16;
constexpr int MAX_MISMATCH_SEARCH = 7;  // -n up to this (Mm_Builder)
constexpr int SAINTERVAL = 8;
// Ranges with End-Start below this count their BWT characters directly (Count_Range)
// instead of a second stored-count lookup. Counting is faster up to ~256 rows (cold) and
// ~400 rows (cached), and much slower beyond; see benchmark/branch_threshold.md.
#ifdef BATMAN_BRANCHTHRESHOLD
constexpr int BRANCHTHRESHOLD = BATMAN_BRANCHTHRESHOLD;
#else
constexpr int BRANCHTHRESHOLD = 128;
#endif
constexpr unsigned char NO_BASE = 4;

// An SA range during the search. Level-1 bases of the current segment (lh, rh or
// the whole read) are accounted for, with Mismatches mismatches at Mismatch_Pos
// (read positions) against the genome bases in Mismatch_Char (2 bits each).
// Skip marks a single-row path (Start == End). With the reference loaded (Use_Ref)
// it also means located: Start = End stay on the sampled row (never 0) and Pos is
// the read's 0-based start in the forward genome, the same for revfmi and fwfmi.
// A located tag needs no rows any more: Search_* hands it to Extend_* (main
// resets Level to 1 for the next segment), and a deferred one is compared again
// at Pos when resumed. Only a finished hit has the rest of its mismatches recorded.
struct SARange {
    unsigned long Start;
    unsigned long End;
    int Level;
    char Mismatches;
    unsigned char Skip;          // nonzero: single row; with the reference, located (Pos valid)
    unsigned char Forbidden_Base = NO_BASE;  // threshold reject: base already tried at Level
    int Stage;                   // a reported mismatch hit: its stage (Result_List, Hit_Stage)
    unsigned char Mismatch_Pos[MAX_MISMATCHES_BOUND];
    unsigned Mismatch_Char;
    unsigned long Pos;           // located: forward-genome start of the read
    unsigned long Other_Start = 0;  // the same string's range start in the other index
    int Anchor = 0;  // the string's start in the read (revfmi ranges) or its end (fwfmi): Locate_Other
};

struct Range {
    unsigned Start;
    unsigned End;
    int Label;
    unsigned Other;  // Other_Start
};

// A result list of the mismatch search (Result_List). The caller clears it before each search.
struct ForwardResults {
    std::vector<SARange> hits;        // exactly Count mismatches
    std::vector<SARange> mismatches;  // deferred parents at Count mismatches: resume at Level, skipping Forbidden_Base
    unsigned long hitCount = 0;       // locations in hits
};

// Each read's result lists for -n 1 to 7 (Res[RL_LISTS]). The stages run in output
// order and add their hits to RL_HITS; End_Stage tags each with its stage, and keeps
// the stage's rejects (collected in RL_HITS.mismatches) as its chain's latest round,
// so each chain of resumed rejects holds one round, whatever the mismatch count.
enum Result_List : unsigned {
    RL_EXACT,              // 0|0 hits; its rejects are the forward chain's first round
    RL_HITS,               // every mismatch hit, in stage order; the running stage's rejects
    RL_FWD,                // the chains' latest rejects, one list each (Mm_Builder): 0|k,
    RL_BWD,                // k|0, then the seed families' in order (1|k, k|1, 2|k, k|2, 3|k for -n 7)
    RL_CHAINS_END = RL_FWD + 7,
    RL_SEED_1 = RL_CHAINS_END,  // the running stage's seeds: in the crossing's index,
    RL_SEED_2,                  // switched to it, and (2- and 3-mismatch seeds extended in
    RL_SEED_3,                  // two steps) before the switch
    RL_SEED_4,
    RL_LISTS,
    RL_NONE = 0xFF         // End_Stage: drop the rejects
};

// A mismatch stage's output: the index its ranges are in (F revfmi, B fwfmi) and the
// mismatch count. Its hits carry its position in the search's table as SARange::Stage.
struct Hit_Stage {
    const char* Dir;
    int Mismatches;
};
// Tags the hits added since Tagged with Stage.
inline void Tag_Hits(ForwardResults* Res, std::size_t& Tagged, int Stage) {
    std::vector<SARange>& Hits = Res[RL_HITS].hits;
    for (; Tagged < Hits.size(); Tagged++) Hits[Tagged].Stage = Stage;
}

// Ends a stage: tags its hits, keeps its rejects as list Keep's (replacing that chain's
// previous round; RL_NONE drops them) and drops its seeds.
inline void End_Stage(ForwardResults* Res, std::size_t& Tagged, int Stage, unsigned Keep) {
    Tag_Hits(Res, Tagged, Stage);
    std::vector<SARange>& Rejects = Res[RL_HITS].mismatches;
    if (Keep != RL_NONE) Res[Keep].mismatches.swap(Rejects);
    Rejects.clear();
    for (unsigned k : {RL_SEED_1, RL_SEED_2, RL_SEED_3, RL_SEED_4}) {
        Res[k].hits.clear();
        Res[k].mismatches.clear();
        Res[k].hitCount = 0;
    }
}

MMPool* mmPool = nullptr;
using ActiveBWT = BWT;
ActiveBWT* revfmi = nullptr;
ActiveBWT* fwfmi = nullptr;  // loaded only for -n 1

int LOOKUP_DEPTH = 6;  // --lookup: requested table depth
int LOOKUPSIZE = 6;    // depth in use, set by Build_Tables
int STRINGLENGTH = 36;

// Start-of-search tables, 4^LOOKUPSIZE entries each, filled by Build_Tables: the SA
// range of every LOOKUPSIZE-base prefix in revfmi (Forward, first base in the low
// bits) and suffix in fwfmi (Backward, last base in the low bits). Start and End
// sit together, so a lookup is one cache line. Depth: benchmark/lookup_table.md.
struct Lookup_Entry {
    unsigned Start;
    unsigned End;
    unsigned Other;  // Other_Start
};
std::vector<Lookup_Entry> Forward_Lookup;
std::vector<Lookup_Entry> Backward_Lookup;
thread_local unsigned Branch_Characters[4] = {0, 0, 0, 0};  // Branch_Detect: per thread
thread_local SARange Branch_Ranges[4];

// Reference for located singletons (see Locate). Off with --no-ref.
bool Use_Ref = false;
// --switch-walk: index switches re-walk the seed's string in the other index, as
// batmis1 (Mm_Start_Switch); off (the default), they read it off the seed's
// Other_Start (Switched_Range).
bool Switch_Walk = false;
// A single row is located from either index's row (default; --own-row: its own index's
// only, as batmis1): the other row stays put on a single-row path, so it is
// tested where the path becomes single (Mm_Search_Enter, Mm_Single_Enter), the moving own
// row after each step. Needs both indexes' SA samples (fwfmi loaded for -n 0 too).
bool Both_Rows = true;
// Whether the backward table (fwfmi) is needed: -n 1 and up.
bool Need_Backward = true;
// Without the reference, but with the SA samples loaded (SAM/BAM output): a single-row
// path records its genome position in Pos at the first sampled row (Mark_Sampled).
bool Track_Pos = false;
constexpr unsigned long NO_POS = ~0UL;  // Pos unknown (Track_Pos: the read runs off the genome)
std::vector<uint8_t> Pac;             // forward genome (.pac): 2 bits a base, first base in bits 7-6, 16 bytes of padding
unsigned long Text_Length = 0;
uint64_t Read_Packed[MAXTAG / 32 + 1];  // current read: 32 bases a word, first base in bits 63-62

int EncodeBase(char c) {
    switch (std::toupper(static_cast<unsigned char>(c))) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default: return -1;
    }
}

// EncodeBase as a table (Run_Batched encodes every base of the input).
const std::array<signed char, 256> Base_Codes = [] {
    std::array<signed char, 256> T;
    for (int c = 0; c < 256; c++) T[c] = static_cast<signed char>(EncodeBase(static_cast<char>(c)));
    return T;
}();

ActiveBWT* initFMI(const char* bwtCodeFileName, const char* bwtOccValueFileName, const char* saValueFileName = nullptr) {
    const int poolSize = 524288;
    MMMasterInitialize(3, 0, FALSE, nullptr);
    mmPool = MMPoolCreate(poolSize);
    return BWTLoad(mmPool, bwtCodeFileName, bwtOccValueFileName, saValueFileName, nullptr, nullptr, nullptr);
}

// Packed defaults to the current read's Read_Packed; batched searches pass their own.
void Pack_Read(std::string_view Read, uint64_t* Packed = Read_Packed) {
    std::memset(Packed, 0, (Read.size() + 31) / 32 * sizeof(uint64_t));
    for (std::size_t i = 0; i < Read.size(); i++) {
        Packed[i / 32] |= static_cast<uint64_t>(Read[i]) << (62 - 2 * (i % 32));
    }
}

// The 32 genome bases from G0, packed as Read_Packed (first base in bits 63-62).
inline uint64_t Genome_Word(unsigned long G0) {
    const uint8_t* P = Pac.data() + G0 / 4;
    uint64_t G;
    std::memcpy(&G, P, 8);
    G = __builtin_bswap64(G);
    const unsigned Shift = 2 * (G0 % 4);
    if (Shift) G = (G << Shift) | (P[8] >> (8 - Shift));
    return G;
}

// For read word W against genome word G: the low bit of each base's pair (bit
// 62-2k for base k) is set where the bases differ, for bases k in [Lo, Hi) only.
inline uint64_t Mismatch_Bits(uint64_t G, int W, int Lo, int Hi, const uint64_t* Packed = Read_Packed) {
    uint64_t X = G ^ Packed[W];
    X = (X | (X >> 1)) & 0x5555555555555555ULL;
    return X & (~0ULL >> (2 * Lo)) & (Hi == 32 ? ~0ULL : ~(~0ULL >> (2 * Hi)));
}

// Mismatches between read[A..B) and genome[S+A..S+B), 32 bases at a time:
// XOR with the genome, fold each base's two bits into one, popcount.
unsigned Genome_Mismatches(unsigned long S, int A, int B, const uint64_t* Packed = Read_Packed) {
    unsigned D = 0;
    for (int W = A / 32; W * 32 < B; W++) {
        D += __builtin_popcountll(Mismatch_Bits(Genome_Word(S + W * 32), W, std::max(A - W * 32, 0), std::min(B - W * 32, 32), Packed));
    }
    return D;
}

// Appends the mismatches of read[A..B) against the genome at Tag.Pos to Tag's
// Mismatch_Pos/Mismatch_Char, in the order a walk would meet them: ascending
// read position for Right, descending otherwise.
void Record_Mismatches(SARange& Tag, bool Right, int A, int B, const uint64_t* Packed = Read_Packed) {
    auto Add = [&](int ReadPos, unsigned GenomeBase) {
        if (Tag.Mismatches >= MAX_MISMATCHES_BOUND) return;
        Tag.Mismatch_Char |= GenomeBase << (Tag.Mismatches * 2);
        Tag.Mismatch_Pos[Tag.Mismatches] = ReadPos;
        Tag.Mismatches++;
    };
    if (A >= B) return;
    if (Right) {
        for (int W = A / 32; W * 32 < B; W++) {
            const uint64_t G = Genome_Word(Tag.Pos + W * 32);
            for (uint64_t Bits = Mismatch_Bits(G, W, std::max(A - W * 32, 0), std::min(B - W * 32, 32), Packed); Bits; ) {
                const int K = (__builtin_clzll(Bits) - 1) / 2;	// lowest read position first
                Add(W * 32 + K, (G >> (62 - 2 * K)) & 3);
                Bits &= ~(1ULL << (62 - 2 * K));
            }
        }
    } else {
        for (int W = (B - 1) / 32; W >= A / 32; W--) {
            const uint64_t G = Genome_Word(Tag.Pos + W * 32);
            for (uint64_t Bits = Mismatch_Bits(G, W, std::max(A - W * 32, 0), std::min(B - W * 32, 32), Packed); Bits; Bits &= Bits - 1) {
                const int K = (62 - __builtin_ctzll(Bits)) / 2;	// highest read position first
                Add(W * 32 + K, (G >> (62 - 2 * K)) & 3);
            }
        }
    }
}

// A singleton path at a sampled row (Start % SAINTERVAL == 0) can stop walking: the
// stored SA value gives the read's genome position, and the rest of the read is
// compared with the packed genome. Done bases of the segment read[Offset..Offset+Len)
// are matched. Right: revfmi, which indexes the reversed genome; the row stands for
// read[0..Offset+Done). Otherwise fwfmi; the row stands for read[Offset+Len-Done..).
// Returns false if the read would run off the genome (where the walk meets $).
#ifdef BATMIS2_VERIFY
extern unsigned long Other_Checked, Other_Failed;
#endif

// Other: from the other index's row (Other_Start) instead: fwfmi's gives the string's
// start in T, read position Anchor for a revfmi range; revfmi's gives n - its end, read
// position Anchor for an fwfmi range (Both_Rows).
bool Locate(const ActiveBWT* fmi, SARange& Tag, bool Right, int Offset, int Len, int Done, bool Other = false) {
    long S;
    if (Other) {
        const ActiveBWT* O = Right ? fwfmi : revfmi;
        const long P = O->saValue[Tag.Other_Start / SAINTERVAL];
        S = Right ? P - Tag.Anchor : static_cast<long>(Text_Length) - P - Tag.Anchor;
#ifdef BATMIS2_VERIFY
        // Against the own row's walk (BWTSaValue).
        const long Q = BWTSaValue(fmi, static_cast<unsigned>(Tag.Start));
        const long Own = Right ? static_cast<long>(Text_Length) - (Offset + Done) - Q : Q - (Offset + Len - Done);
        Other_Checked++;
        if (Own != S && Other_Failed++ < 5)
            std::fprintf(stderr, "Verify_Other: other row locates at %ld, own row at %ld (Right %d Offset %d Len %d Done %d Anchor %d Mismatches %d)\n",
                         S, Own, Right, Offset, Len, Done, Tag.Anchor, Tag.Mismatches);
#endif
    } else {
        const long P = fmi->saValue[Tag.Start / SAINTERVAL];
        S = Right ? static_cast<long>(Text_Length) - (Offset + Done) - P : P - (Offset + Len - Done);
    }
    Tag.Skip = 1;
    Tag.End = Tag.Start;
    if (S < 0 || S + STRINGLENGTH > static_cast<long>(Text_Length)) return false;
    Tag.Pos = S;
    return true;
}

// Whether a single row's other-index row is sampled, to locate from (Both_Rows).
inline bool Other_Sampled(const SARange& T) {
    return Use_Ref && Both_Rows && T.Other_Start && T.Other_Start % SAINTERVAL == 0;
}

// A single-row path without the reference reaching a sampled row: from now on Skip (a
// single row; the walk goes on). With Track_Pos, also the read's forward-genome start in
// Pos, read from the SA sample as Locate does (NO_POS if the read would run off the
// genome): later steps follow the same occurrence, so it stays valid, and the hit needs
// no SA walk at output. Done bases of the segment are matched at this row.
inline void Mark_Sampled(const ActiveBWT* fmi, SARange& Tag, bool Right, int Offset, int Len, int Done) {
    Tag.Skip++;
    Tag.End = Tag.Start;
    if (Track_Pos && !Locate(fmi, Tag, Right, Offset, Len, Done)) Tag.Pos = NO_POS;
}

// Mismatches in the part of the segment a located Tag has not matched yet.
unsigned Rest_Mismatches(const SARange& Tag, bool Right, int Offset, int Len, int Done, const uint64_t* Packed = Read_Packed) {
    return Right ? Genome_Mismatches(Tag.Pos, Offset + Done, Offset + Len, Packed)
                 : Genome_Mismatches(Tag.Pos, Offset, Offset + Len - Done, Packed);
}

// The first read position in [A, B) that differs from the genome at S, in the order a
// walk meets them (Right: lowest, else highest), or -1.
int First_Mismatch(unsigned long S, int A, int B, bool Right, const uint64_t* Packed = Read_Packed) {
    if (A >= B) return -1;
    if (Right) {
        for (int W = A / 32; W * 32 < B; W++)
            if (const uint64_t Bits = Mismatch_Bits(Genome_Word(S + W * 32), W, std::max(A - W * 32, 0), std::min(B - W * 32, 32), Packed))
                return W * 32 + (__builtin_clzll(Bits) - 1) / 2;
    } else {
        for (int W = (B - 1) / 32; W >= A / 32; W--)
            if (const uint64_t Bits = Mismatch_Bits(Genome_Word(S + W * 32), W, std::max(A - W * 32, 0), std::min(B - W * 32, 32), Packed))
                return W * 32 + (62 - __builtin_ctzll(Bits)) / 2;
    }
    return -1;
}

// Whether the rest of the segment matches at a located Tag (as Rest_Mismatches == 0).
// If not, Tag.Level is set to the failure depth: the segment's first Level bases (in
// the search direction) do not occur, since the located row was their only prefix.
bool Located_Matches(SARange& Tag, bool Right, int Offset, int Len, int Done, const uint64_t* Packed = Read_Packed) {
    const int X = Right ? First_Mismatch(Tag.Pos, Offset + Done, Offset + Len, true, Packed)
                        : First_Mismatch(Tag.Pos, Offset, Offset + Len - Done, false, Packed);
    if (X < 0) return true;
    Tag.Level = Right ? X - Offset + 1 : Offset + Len - X;
    return false;
}

// An exact search at a sampled row: matches if the rest of the segment does.
// On failure Tag.Level is the failure depth (Located_Matches; the whole segment when
// the read would run off the genome).
bool Search_Located(const ActiveBWT* fmi, SARange& Tag, bool Right, int Offset, int Len, int Done,
                    const uint64_t* Packed = Read_Packed) {
    if (!Locate(fmi, Tag, Right, Offset, Len, Done)) {
        Tag.Level = Len;
        return false;
    }
    return Located_Matches(Tag, Right, Offset, Len, Done, Packed);
}

// An extension of a located Tag, whose first Level-1 segment bases are matched with
// Mismatches mismatches. As the walk would: a hit at exactly Count
// mismatches (with the rest of its mismatches recorded); above Count, deferred
// unchanged, to be compared again when resumed with a larger Count.
void Extend_Located(const SARange& Tag, bool Right, int Offset, int Len, int Count, ForwardResults& Out,
                    const uint64_t* Packed = Read_Packed) {
    const int Done = Tag.Level - 1;
    const unsigned D = Tag.Mismatches + Rest_Mismatches(Tag, Right, Offset, Len, Done, Packed);
    if (D == static_cast<unsigned>(Count)) {
        SARange Hit = Tag;
        if (Right) Record_Mismatches(Hit, true, Offset + Done, Offset + Len, Packed);
        else Record_Mismatches(Hit, false, Offset, Offset + Len - Done, Packed);
        Out.hits.push_back(Hit);
        Out.hitCount++;
    } else if (D > static_cast<unsigned>(Count)) {
        Out.mismatches.push_back(Tag);
    }
}

// Prefetches the occurrence and BWT words an Occ lookup reads for row Index
// (as BWTPrefetchOccValueExplicit/BWTPrefetchBWT in index/BWT.c).
PREFETCH_ONLY void Prefetch_Row(const ActiveBWT* fmi, unsigned long Index) {
    Index -= (Index > fmi->inverseSa0);
    const unsigned long Explicit = (Index + OCC_INTERVAL / 2 - 1) / OCC_INTERVAL;
    __builtin_prefetch(fmi->occValue + Explicit / OCC_VALUE_PER_WORD * ALPHABET_SIZE);
    __builtin_prefetch(fmi->occValueMajor + Explicit * OCC_INTERVAL / OCC_INTERVAL_MAJOR * ALPHABET_SIZE);
    __builtin_prefetch(fmi->bwtCode + Index / CHAR_PER_WORD);
    __builtin_prefetch(fmi->bwtCode + Explicit * OCC_INTERVAL / CHAR_PER_WORD);
}

// Counts each base among the BWT characters of rows Start..End. Row inverseSa0
// holds $, which is not stored in bwtCode, so it is skipped (as in BWTOccValue).
inline void Count_Range(const ActiveBWT* fmi, unsigned long Start, unsigned long End, unsigned Count[4]) {
    fm_occ::Count_Codes(fmi->bwtCode, Start - (Start > fmi->inverseSa0), End + 1 - (End + 1 > fmi->inverseSa0), Count);
}

#ifdef BATMIS2_VERIFY
// Checks Other_Start of a range of up to 32 rows in fmi: the other index's rows from
// Other give the same occurrences, a string length apart (fwfmi gives starts, revfmi
// n - ends). Counts the checks and failures (printed at exit).
unsigned long Other_Checked = 0, Other_Failed = 0;
void Verify_Other(const ActiveBWT* fmi, unsigned long Start, unsigned long End, unsigned long Other) {
    if (!fwfmi || !Start || End < Start || End - Start >= 32) return;
    const ActiveBWT* O = fmi == revfmi ? fwfmi : revfmi;
    std::vector<long> A, B;
    for (unsigned long r = Start; r <= End; r++) {
        const long v = BWTSaValue(fmi, static_cast<unsigned>(r));
        A.push_back(fmi == fwfmi ? v : static_cast<long>(Text_Length) - v);
    }
    for (unsigned long r = Other; r <= Other + (End - Start); r++) {
        const long v = BWTSaValue(O, static_cast<unsigned>(r));
        B.push_back(O == fwfmi ? v : static_cast<long>(Text_Length) - v);
    }
    std::sort(A.begin(), A.end());
    std::sort(B.begin(), B.end());
    const std::vector<long>& Starts = fmi == fwfmi ? A : B;
    const std::vector<long>& Ends = fmi == fwfmi ? B : A;
    Other_Checked++;
    for (std::size_t k = 0; k < A.size(); k++)
        if (Ends[k] - Starts[k] != Ends[0] - Starts[0] || Ends[0] <= Starts[0]) {
            if (Other_Failed++ < 5) std::fprintf(stderr, "Verify_Other: rows %lu..%lu other %lu differ\n", Start, End, Other);
            return;
        }
}
struct Verify_Report {
    ~Verify_Report() { std::fprintf(stderr, "Verify_Other: checked=%lu failed=%lu\n", Other_Checked, Other_Failed); }
} Verify_Report_At_Exit;
#endif

// Get_SARange with Other_Start: it moves past the range's rows whose BWT character is $
// or smaller than newChar, which come first in the other index. Ranges below
// BRANCHTHRESHOLD have the counts already; wider ones take all four bases' Occ at both
// ends (batmis1 takes newChar's alone), the bidirectional step's extra cost.
void Get_SARange(long newChar, unsigned long& start, unsigned long& end, unsigned long& other, ActiveBWT* fmi) {
    const unsigned long Dollar = fmi->inverseSa0 >= start && fmi->inverseSa0 <= end;
    unsigned count[4];
    if (end - start < BRANCHTHRESHOLD) {
        Count_Range(fmi, start, end, count);
        if (!count[newChar]) {
            start = 0;
            return;
        }
        start = fmi->cumulativeFreq[newChar] + fm_occ::Occ_One(fmi, start, newChar) + 1;
        end = start + count[newChar] - 1;
    } else {
        unsigned At_Start[4], At_End[4];
        fm_occ::Occ_All(fmi, start, At_Start);
        fm_occ::Occ_All(fmi, end + 1, At_End);
        for (int b = 0; b < 4; b++) count[b] = At_End[b] - At_Start[b];
        start = fmi->cumulativeFreq[newChar] + At_Start[newChar] + 1;
        end = fmi->cumulativeFreq[newChar] + At_End[newChar];
        if (end < start) {
            start = 0;
            return;
        }
    }
    other += Dollar;
    for (int b = 0; b < newChar; b++) other += count[b];
#ifdef BATMIS2_VERIFY
    Verify_Other(fmi, start, end, other);
#endif
}

inline void Get_SARange(long newChar, SARange& range, ActiveBWT* fmi) {
    Get_SARange(newChar, range.Start, range.End, range.Other_Start, fmi);
}

// Returns Base unchanged, as a value the compiler cannot relate to anything
// else. Take it from the read's base before comparing with the BWT character:
// the LF step on a match then uses the read's register, which is ready before
// the BWT load completes, so the occurrence lookup can start in parallel
// instead of waiting for that load (~2x faster singleton walks). Without it the
// compiler may pass the (equal) BWT character instead.
inline unsigned Read_Base(unsigned Base) { asm("" : "+r"(Base)); return Base; }

Range InitialBaseRange(ActiveBWT* fmi, int base) {
    Range range{};
    range.Start = fmi->cumulativeFreq[base] + 1;
    range.End = fmi->cumulativeFreq[base + 1];
    if (range.Start > range.End) range.Start = 0;
    range.Other = range.Start;  // a 1-base string: the same rows in both indexes
    return range;
}

void Build_Preindex_Forward(Range range, int level, int bit) {
    if (LOOKUPSIZE == level) {
        assert(range.Start <= range.End);
        range.Label = range.Label | (bit << (2 * (level - 1)));
        if (range.Start) {
            unsigned long start = range.Start, end = range.End, other = range.Other;
            Get_SARange(bit, start, end, other, revfmi);
            range.Start = static_cast<unsigned>(start);
            range.End = static_cast<unsigned>(end);
            range.Other = static_cast<unsigned>(other);
        }
        Forward_Lookup[range.Label] = {range.Start, range.End, range.Other};
        assert(range.Start <= range.End);
    } else {
        assert(range.Start <= range.End);
        range.Label = range.Label | (bit << (2 * (level - 1)));
        if (range.Start) {
            unsigned long start = range.Start, end = range.End, other = range.Other;
            Get_SARange(bit, start, end, other, revfmi);
            range.Start = static_cast<unsigned>(start);
            range.End = static_cast<unsigned>(end);
            range.Other = static_cast<unsigned>(other);
        }
        level++;
        for (int i = 0; i < 4; i++) {
            Build_Preindex_Forward(range, level, i);
        }
    }
}

// Same as Build_Preindex_Forward, on fwfmi: the lookup key is the last
// LOOKUPSIZE bases of the read, last base first.
void Build_Preindex_Backward(Range range, int level, int bit) {
    assert(range.Start <= range.End);
    range.Label = range.Label | (bit << (2 * (level - 1)));
    if (range.Start) {
        unsigned long start = range.Start, end = range.End, other = range.Other;
        Get_SARange(bit, start, end, other, fwfmi);
        range.Start = static_cast<unsigned>(start);
        range.End = static_cast<unsigned>(end);
        range.Other = static_cast<unsigned>(other);
    }
    if (LOOKUPSIZE == level) {
        Backward_Lookup[range.Label] = {range.Start, range.End, range.Other};
        assert(range.Start <= range.End);
    } else {
        level++;
        for (int i = 0; i < 4; i++) {
            Build_Preindex_Backward(range, level, i);
        }
    }
}

void Build_Tables(ActiveBWT* fmi) {
    // The -n 1 halves search on from Level LOOKUPSIZE+1, so the key must be shorter
    // than the left half. Reads too short for that (under 8 bp) keep the minimum 3.
    LOOKUPSIZE = std::max(3, std::min(LOOKUP_DEPTH, STRINGLENGTH / 2 - 1));
    Forward_Lookup.assign(1ul << (2 * LOOKUPSIZE), Lookup_Entry{});
    if (fwfmi && Need_Backward) Backward_Lookup.assign(1ul << (2 * LOOKUPSIZE), Lookup_Entry{});

    Range lRange{};
    for (int i = 0; i < 4; i++) {
        lRange = InitialBaseRange(fmi, i);
        lRange.Label = i;
        for (int j = 0; j < 4; j++) {
            Build_Preindex_Forward(lRange, 2, j);
        }
    }

    if (fwfmi && Need_Backward) {
        for (int i = 0; i < 4; i++) {
            lRange = InitialBaseRange(fwfmi, i);
            lRange.Label = i;
            for (int j = 0; j < 4; j++) {
                Build_Preindex_Backward(lRange, 2, j);
            }
        }
    }
}

// Search_Right along a single SA path (tag.End == tag.Start or tag.Skip).
// Kept out of line: inlining it into the read loop costs ~10%.
__attribute__((noinline)) bool Search_Right_Singleton(std::string_view currentTag, int offset, SARange& tag, ActiveBWT* fmi) {
    const int stringLength = static_cast<int>(currentTag.size());
    unsigned long index;
    unsigned int now;

    if (tag.Start % SAINTERVAL == 0 && !tag.Skip) {
        if (Use_Ref) {
            if (!Search_Located(fmi, tag, true, offset, stringLength, tag.Level - 1)) return false;
            tag.Level = stringLength;
            return true;
        }
        Mark_Sampled(fmi, tag, true, offset, stringLength, tag.Level - 1);
    }

    for (;;) {
        // Row inverseSa0 holds $ (not stored in bwtCode): the path ends there.
        if (tag.Start == fmi->inverseSa0) {
            return false;
        }
        index = tag.Start - (tag.Start > fmi->inverseSa0);
        now = fmi->bwtCode[index / 16] << ((index % 16) * 2) >> (BITS_IN_WORD - 2);
        const unsigned int expected = static_cast<unsigned char>(currentTag[tag.Level - 1]);
        const unsigned base = Read_Base(expected);
        if (expected == now) {
            tag.Start = fmi->cumulativeFreq[base] + fm_occ::Occ_One(fmi, tag.Start, base) + 1;
            if (!tag.Skip && tag.Start % SAINTERVAL == 0) {
                if (Use_Ref) {
                    if (!Search_Located(fmi, tag, true, offset, stringLength, tag.Level)) return false;
                    tag.Level = stringLength;
                    return true;
                }
                Mark_Sampled(fmi, tag, true, offset, stringLength, tag.Level);
            }

            if (tag.Level == stringLength) {
                //if (!tag.Skip) {
                    tag.End = tag.Start;
                //}
                return true;
            }
            tag.Level++;
            continue;
        }
        return false;
    }
}

// Extends tag to the right over currentTag[Level-1] with no mismatches.
bool Search_Right(std::string_view currentTag, int offset, SARange& tag, ActiveBWT* fmi) {
    const int stringLength = static_cast<int>(currentTag.size());

    if (!tag.Start) {
        return false;
    }

    for (;;) {
        if (tag.End == tag.Start || tag.Skip) {
            return Search_Right_Singleton(currentTag, offset, tag, fmi);
        } else {
            Get_SARange(currentTag[tag.Level - 1], tag, fmi);

            if (tag.Start != 0) {
                if (tag.Level == stringLength) {
                    return true;
                }
                tag.Level++;
                continue;
            }
            return false;
        }
    }
}

void Branch_Detect (const struct SARange Tag,BWT *fmi,int Start)
{

	Branch_Characters[0]=0;Branch_Characters[1]=0;Branch_Characters[2]=0;Branch_Characters[3]=0;
	if(Tag.End-Tag.Start<BRANCHTHRESHOLD)//only small number of branches
	{
		// Occurrences at Start, plus the range's own characters (no second lookup).
		unsigned int ALIGN_16 Occ_Start[ALPHABET_SIZE];
		unsigned Count[4];
		fm_occ::Occ_All(fmi, Tag.Start, Occ_Start);
		Count_Range(fmi, Tag.Start, Tag.End, Count);
		// Other_Start: past the range's $ row and the smaller bases' rows.
		unsigned long Other = Tag.Other_Start + (fmi->inverseSa0 >= Tag.Start && fmi->inverseSa0 <= Tag.End);
		for (int Branch=0;Branch<4;Branch++)
		{
			Branch_Characters[Branch]=Count[Branch];
			if (Count[Branch])
			{
				Branch_Ranges[Branch].Start = fmi->cumulativeFreq[Branch] + Occ_Start[Branch] + 1;
				Branch_Ranges[Branch].End = Branch_Ranges[Branch].Start + Count[Branch]-1;// Calculate SAranges
				Branch_Ranges[Branch].Other_Start = Other;
#ifdef BATMIS2_VERIFY
				Verify_Other(fmi, Branch_Ranges[Branch].Start, Branch_Ranges[Branch].End, Other);
#endif
			}
			Other += Count[Branch];
		}
	}
	else
	{
		// Occurrences of all four bases at both ends.
		unsigned int ALIGN_16 Occ_Start[ALPHABET_SIZE], Occ_End[ALPHABET_SIZE];
		fm_occ::Occ_All(fmi, Tag.Start, Occ_Start);
		fm_occ::Occ_All(fmi, Tag.End + 1, Occ_End);
		unsigned long Other = Tag.Other_Start + (fmi->inverseSa0 >= Tag.Start && fmi->inverseSa0 <= Tag.End);
		for (int Branch=0;Branch<4;Branch++)
		{
				Branch_Ranges[Branch].Start = fmi->cumulativeFreq[Branch] + Occ_Start[Branch] + 1;
				Branch_Ranges[Branch].End = fmi->cumulativeFreq[Branch] + Occ_End[Branch];
				Branch_Ranges[Branch].Other_Start = Other;
				if(!(Branch_Ranges[Branch].End<Branch_Ranges[Branch].Start)) Branch_Characters[Branch]=1;
#ifdef BATMIS2_VERIFY
				if (Branch_Characters[Branch]) Verify_Other(fmi, Branch_Ranges[Branch].Start, Branch_Ranges[Branch].End, Other);
#endif
				Other += Occ_End[Branch] - Occ_Start[Branch];
		}

	}
}
// Range of the read's first LOOKUPSIZE bases in revfmi, from the lookup table,
// with Level at the next base to match. A read no longer than the lookup key is
// matched whole from its first base instead (the table is keyed on LOOKUPSIZE
// bases); its Level is then past the end of the read, and there is nothing left
// to search.
SARange Start_Range(std::string_view read) {
    SARange range{};
    const int length = static_cast<int>(read.size());
    if (length <= LOOKUPSIZE) {
        const Range first = InitialBaseRange(revfmi, read[0]);
        range.Start = first.Start;
        range.End = first.End;
        range.Other_Start = first.Other;
        for (int i = 1; i < length && range.Start; i++) Get_SARange(read[i], range, revfmi);
        range.Level = length + 1;
        return range;
    }
    int c = 0;
    for (int i = 0; i < LOOKUPSIZE; i++) c |= read[i] << (2 * i);
    range.Start = Forward_Lookup[c].Start;
    range.End = Forward_Lookup[c].End;
    range.Other_Start = Forward_Lookup[c].Other;
    range.Level = LOOKUPSIZE + 1;
    return range;
}

// Range of Seg in revfmi (Right: from its first bases) or fwfmi (from its last bases),
// from the lookup tables, with Level at the next base to match. A segment no longer
// than the lookup key is matched whole from its first (Right) or last base, and its
// Level is then past its end.
SARange Segment_Start(std::string_view Seg, bool Right) {
    if (Right) return Start_Range(Seg);
    SARange R{};
    const int Len = static_cast<int>(Seg.size());
    if (Len <= LOOKUPSIZE) {
        const Range Last = InitialBaseRange(fwfmi, Seg[Len - 1]);
        R.Start = Last.Start;
        R.End = Last.End;
        R.Other_Start = Last.Other;
        for (int i = 1; i < Len && R.Start; i++) Get_SARange(Seg[Len - 1 - i], R, fwfmi);
        R.Level = Len + 1;
        return R;
    }
    int bc = 0;
    for (int i = 0; i < LOOKUPSIZE; i++) bc |= Seg[Len - 1 - i] << (2 * i);
    R.Start = Backward_Lookup[bc].Start;
    R.End = Backward_Lookup[bc].End;
    R.Other_Start = Backward_Lookup[bc].Other;
    R.Level = LOOKUPSIZE + 1;
    return R;
}

// Read segments, as fractions u of a part of the read (Mm_Part). In a half, u counts from
// its outer end (LH: the read's start; RH: its end): u is at LH·u in LH and at
// LH + R·(1 − u) in RH (R = L − LH), floor in both, so RH's mirror of an LH segment has the
// same fractions. Across the whole read, u counts from its start. A segment is [Lo, Hi) in u.
// Mm_Step::Seg and Seg2 index the segment table, which the program builder fills.
enum Mm_Part : unsigned char { PART_READ, PART_LH, PART_RH };
struct Mm_Fraction {
    unsigned char Num, Den;
};
struct Mm_Segment_Spec {
    Mm_Part Part;
    Mm_Fraction Lo, Hi;
};
constexpr unsigned char SEG_ALL = 0, SEG_LH = 1, SEG_RH = 2;  // the table's first entries
std::vector<Mm_Segment_Spec> Segment_Table = {
    {PART_READ, {0, 1}, {1, 1}}, {PART_LH, {0, 1}, {1, 1}}, {PART_RH, {0, 1}, {1, 1}}};

inline int Part_Position(Mm_Part Part, Mm_Fraction U, int L) {
    const int LH = L / 2, R = L - LH;
    if (Part == PART_READ) return L * U.Num / U.Den;
    if (Part == PART_LH) return LH * U.Num / U.Den;
    return LH + R * (U.Den - U.Num) / U.Den;
}

// Segment Seg of a read of L bases: [Start, End).
inline void Segment_Bounds(unsigned Seg, int L, int& Start, int& End) {
    const Mm_Segment_Spec& S = Segment_Table[Seg];
    Start = Part_Position(S.Part, S.Part == PART_RH ? S.Hi : S.Lo, L);
    End = Part_Position(S.Part, S.Part == PART_RH ? S.Lo : S.Hi, L);
}

inline std::string_view Read_Segment_Of(std::string_view All, unsigned Seg, int& Offset) {
    int End;
    Segment_Bounds(Seg, static_cast<int>(All.size()), Offset, End);
    return All.substr(Offset, End - Offset);
}

// Mismatch windows for the windowed extensions: read positions [Lo, Hi) that must hold
// exactly Need mismatches. Positions outside every window must match.
struct Mm_Windows {
    int N = 0;
    int Lo[3], Hi[3], Need[3];
};

// A window set (Mm_Step::Win indexes the window table; 0 is none): up to three segments,
// each with its count and how the exact half searches' failure depths bound it
// (One_Mismatch passes: read[0, LH_Fail) and read[RH_Fail, L) do not occur, so a hit's
// first LH mismatch is before LH_Fail and its last RH mismatch at or after RH_Fail):
//   WIN_FREE        not bounded;
//   WIN_NARROW      it holds that first (LH) or last (RH) mismatch alone: Hi <= LH_Fail,
//                   or Lo >= RH_Fail;
//   WIN_AFTER_SEED  its case's seed is the half's exact outer part before it, so the case
//                   cannot hold once LH_Fail <= Lo (or RH_Fail >= Hi).
enum Mm_Window_Kind : unsigned char { WIN_FREE, WIN_NARROW, WIN_AFTER_SEED };
struct Mm_Window_Spec {
    unsigned char N;
    unsigned char Seg[3], Need[3];
    Mm_Window_Kind Kind[3];
};
std::vector<Mm_Window_Spec> Window_Table(1);  // 0: no window

// The windows of Set for a read of L bases, bounded by the failure depths.
inline Mm_Windows Window_Set(unsigned Set, int L, int LH_Fail, int RH_Fail) {
    const Mm_Window_Spec& S = Window_Table[Set];
    const int LH = L / 2;
    Mm_Windows W;
    for (int i = 0; i < S.N; i++) {
        int Lo, Hi;
        Segment_Bounds(S.Seg[i], L, Lo, Hi);
        const bool Left = Hi <= LH;
        if (S.Kind[i] == WIN_NARROW) {
            if (Left) Hi = std::min(Hi, LH_Fail);
            else Lo = std::max(Lo, RH_Fail);
        } else if (S.Kind[i] == WIN_AFTER_SEED) {
            if (Left && LH_Fail <= Lo) Hi = Lo;
            if (!Left && RH_Fail >= Hi) Lo = Hi;
        }
        W.Lo[W.N] = Lo;
        W.Hi[W.N] = Hi;
        W.Need[W.N] = S.Need[i];
        W.N++;
    }
    return W;
}

// Whether every window can hold its count (else the case has no hits).
inline bool Windows_Possible(const Mm_Windows& W) {
    for (int i = 0; i < W.N; i++)
        if (W.Hi[i] - W.Lo[i] < W.Need[i]) return false;
    return true;
}

inline int Window_Of(const Mm_Windows& W, int x) {
    for (int i = 0; i < W.N; i++)
        if (x >= W.Lo[i] && x < W.Hi[i]) return i;
    return -1;
}

inline int Window_Count(const SARange& T, const Mm_Windows& W, int w) {
    int c = 0;
    for (int i = 0; i < T.Mismatches && i < MAX_MISMATCHES_BOUND; i++)
        c += T.Mismatch_Pos[i] >= W.Lo[w] && T.Mismatch_Pos[i] < W.Hi[w];
    return c;
}

// Whether read position x may take a mismatch: it is in a window that has room.
inline bool Window_Can_Change(const SARange& T, const Mm_Windows& W, int x) {
    const int w = Window_Of(W, x);
    return w >= 0 && Window_Count(T, W, w) < W.Need[w];
}

// After position x is matched (Right: positions ascending, else descending): whether
// x's window can still reach its count with the positions it has left.
inline bool Window_Keep(const SARange& T, const Mm_Windows& W, int x, bool Right) {
    const int w = Window_Of(W, x);
    if (w < 0) return true;
    const int Left = Right ? W.Hi[w] - 1 - x : x - W.Lo[w];
    return Window_Count(T, W, w) + Left >= W.Need[w];
}

// A located hit: every window has its count and nothing else in [A, B) changed.
inline bool Windows_Hold(const SARange& T, const Mm_Windows& W, int A, int B) {
    for (int w = 0; w < W.N; w++)
        if (Window_Count(T, W, w) != W.Need[w]) return false;
    for (int i = 0; i < T.Mismatches && i < MAX_MISMATCHES_BOUND; i++)
        if (T.Mismatch_Pos[i] >= A && T.Mismatch_Pos[i] < B && Window_Of(W, T.Mismatch_Pos[i]) < 0) return false;
    return true;
}


// The search engine (--batch B, or auto: Batch_Chooser). Each step of
// a search waits for one cache miss, and one read's steps depend on each other, so B
// reads ("lanes") are searched at once: each round prefetches what every lane's next
// step reads, then steps every lane. A lane runs its read's pass program (Mm_Program,
// built by Mm_Builder) one task at a time: an exact search of a segment (with its
// singleton walk), a depth-first extension (deferring at its Count; its singleton walk
// and the located path), a windowed extension, or an index switch. Each task is a loop
// cut at each memory access, so a read's hits do not depend on B. Exact search (-n 0)
// gives the same ranges as Search_Right.
// Branch_Detect's globals are written and read within one step.
// See benchmark/batched_reads.md.
constexpr unsigned MAX_BATCH = 64;

enum Mm_Sub : unsigned char { MM_POP, MM_WIDE, MM_SINGLE, MM_LOCATE, MM_LOCATED, MM_SWITCH };
// A lane's task: an exact search, a depth-first extension (Count mismatches, deferring
// at the threshold), a windowed extension (Mm_Windows, no deferral) or an index switch.
enum Mm_Task : unsigned char { TASK_SEARCH, TASK_EXTEND, TASK_WINDOW, TASK_SWITCH };

// A read's pass program, on the read's segments (Segment_Table).
enum Mm_Op : unsigned char {
    OP_SEARCH,          // exact search of Seg (Right: forward table + revfmi; else backward table + fwfmi);
                        // with Win, skipped when those windows cannot hold (Windows_Possible); with
                        // Stored, the range kept in that slot by OP_SEARCH_KEEP, if there is one
    OP_SEARCH_KEEP,     // OP_SEARCH of Seg that keeps, in slot Stored, the range of its prefix (Right)
                        // or suffix Seg2 (Mm_Keep_Rest)
    OP_EXTEND,          // extend the last search's range over Seg with Count mismatches into list Out
    OP_COUNT_EXACT,     // the read's exact hit count from list From's hits (each capped at MaxHits)
    OP_RESUME,          // extend each range list From deferred, in order, as OP_EXTEND
    OP_CROSS,           // extend each hit of list From (a seed) over Seg, keeping its mismatches; with
                        // Switch, a seed not located first gets the other index's range of Seg2
                        // (with its mismatch applied) by a Start-only walk (Mm_Start_Switch)
    OP_END_STAGE,       // End_Stage: tag the stage's hits with stage Count, keep its rejects in list Out
    OP_STOP_IF_MAXHITS, // stop once the hits so far reach MaxHits
    OP_END
};
struct Mm_Step {
    Mm_Op Op;
    bool Right;
    unsigned char Seg, Count, Out, From;
    bool Limited;  // extensions stop at MaxHits (less the hits already counted)
    bool Switch;
    unsigned char Win, Seg2;  // Win: windowed extension, or a search's case (Window_Table); Seg2: OP_CROSS's switch
                              // segment, or OP_SEARCH_KEEP's kept segment
    unsigned char Stored;     // OP_SEARCH / OP_SEARCH_KEEP: 1 + the slot of the kept range (0: none)
    unsigned char Strand;     // two-strand programs (Zig_Zag): 0 the read, 1 its reverse complement
};
constexpr unsigned MM_PROGRAM_MAX = 320;  // longest program (-n 7 zig-zag: 307 steps), for the stats
constexpr unsigned MM_LISTS = RL_LISTS;  // each read's lists (Result_List)

// The mismatch programs, built by one rule (benchmark/five_mismatch.md). A hit splits as
// i | j, i mismatches in LH and j in RH:
//   - 0|j and j|0: the exact half, then its rejects resumed round after round (a chain);
//   - i|j with 1 <= i <= j: seeds of LH with exactly i mismatches, crossing RH in round
//     2i and resumed as a chain after; i|j with i > j >= 1: the mirror, seeds of RH with
//     exactly j, crossing LH in round 2j + 1;
//   - seeds of a half with exactly s mismatches are the -n s search on that half, keeping
//     its exactly-s cases (Mm_Builder::Seeds);
//   - -n 2 alone finds 1|1 by the read's thirds instead.
// The stages run, and report, by mismatch count; with -m the run stops between rounds
// once the hits reach MaxHits.
struct Mm_Program {
    std::vector<Mm_Step> Steps;
    std::vector<Hit_Stage> Stages;
};

class Mm_Builder {
  public:
    explicit Mm_Builder(int K) : K(K) {}

    Mm_Program Build() {
        if (K == 0) {
            Add({OP_SEARCH, true, SEG_ALL});
        } else {
            One_Mismatch();
            if (K == 2) Two_Mismatch_Thirds();
            for (int C = 2; K >= 3 && C <= K; C++) {
                Add({OP_STOP_IF_MAXHITS});
                for (const Chain& Ch : Chains) {  // the chains' round C
                    Add({OP_RESUME, Ch.Right, Ch.Right ? SEG_RH : SEG_LH, U(C), RL_HITS, Ch.List, true});
                    End_Stage(Ch.Right, C, C < K ? Ch.List : RL_NONE);
                }
                Add({OP_STOP_IF_MAXHITS});
                Family(C % 2 ? PART_RH : PART_LH, C / 2, C);  // the round's new seeds
            }
        }
        Add({OP_END});
        return P;
    }

  private:
    struct Chain {
        unsigned char List;
        bool Right;  // resumed right over RH (else left over LH)
    };
    int K;
    Mm_Program P;
    std::vector<Chain> Chains;

    static unsigned char U(int X) { return static_cast<unsigned char>(X); }
    void Add(const Mm_Step& S) { P.Steps.push_back(S); }
    void End_Stage(bool Right, int Mismatches, unsigned Keep) {
        Add({OP_END_STAGE, false, 0, U(P.Stages.size()), U(Keep)});
        P.Stages.push_back({Right ? "F" : "B", Mismatches});
    }

    static unsigned char Segment(Mm_Part Part, Mm_Fraction Lo, Mm_Fraction Hi) {
        for (std::size_t i = 0; i < Segment_Table.size(); i++) {
            const Mm_Segment_Spec& S = Segment_Table[i];
            if (S.Part == Part && S.Lo.Num * Lo.Den == Lo.Num * S.Lo.Den && S.Hi.Num * Hi.Den == Hi.Num * S.Hi.Den)
                return U(i);
        }
        Segment_Table.push_back({Part, Lo, Hi});
        return U(Segment_Table.size() - 1);
    }
    struct Window {
        unsigned char Seg, Need;
        Mm_Window_Kind Kind;
    };
    static unsigned char Windows(std::initializer_list<Window> Ws) {
        Mm_Window_Spec W{};
        for (const Window& X : Ws) {
            W.Seg[W.N] = X.Seg;
            W.Need[W.N] = X.Need;
            W.Kind[W.N] = X.Kind;
            W.N++;
        }
        for (std::size_t i = 1; i < Window_Table.size(); i++) {
            const Mm_Window_Spec& S = Window_Table[i];
            if (S.N == W.N && std::equal(S.Seg, S.Seg + S.N, W.Seg) && std::equal(S.Need, S.Need + S.N, W.Need) &&
                std::equal(S.Kind, S.Kind + S.N, W.Kind))
                return U(i);
        }
        Window_Table.push_back(W);
        return U(Window_Table.size() - 1);
    }

    // 0|0, 0|1 from the deferred 0|0 ranges, then 1|0. The exact half scans keep a part's
    // range for later seeds: the read's outer thirds (-n 2), or the outer quarters.
    void One_Mismatch() {
        if (K == 1) {
            Add({OP_SEARCH, true, SEG_LH});
        } else {
            const unsigned char Kept = K == 2 ? Segment(PART_READ, {0, 1}, {1, 3}) : Segment(PART_LH, {0, 1}, {1, 2});
            Add({OP_SEARCH_KEEP, true, SEG_LH, 0, 0, 0, false, false, 0, Kept, 1});
        }
        Add({OP_EXTEND, true, SEG_RH, 0, RL_EXACT});
        Add({OP_COUNT_EXACT, false, 0, 0, 0, RL_EXACT});
        Add({OP_RESUME, true, SEG_RH, 1, RL_HITS, RL_EXACT, true});
        End_Stage(true, 1, RL_FWD);
        Add({OP_STOP_IF_MAXHITS});
        if (K == 1) {
            Add({OP_SEARCH, false, SEG_RH});
        } else {
            const unsigned char Kept = K == 2 ? Segment(PART_READ, {2, 3}, {1, 1}) : Segment(PART_RH, {0, 1}, {1, 2});
            Add({OP_SEARCH_KEEP, false, SEG_RH, 0, 0, 0, false, false, 0, Kept, 2});
        }
        Add({OP_EXTEND, false, SEG_LH, 1, RL_HITS, 0, true});
        End_Stage(false, 1, RL_BWD);
        Chains = {{RL_FWD, true}, {RL_BWD, false}};
    }

    // -n 2: 0|2 and 2|0 from the chains, then 1|1 by the read's thirds T1 = [0, a),
    // T2 = [a, b), T3 = [b, L) (a < LH < b), p in LH and q in RH:
    //   (a) p in [a, LH): T1 exact, extended right over [a, L);
    //   (b) p in T1, q in [LH, b): T3 exact, extended left over [0, b);
    //   (c) p in T1, q in T3: T2 exact, extended right over T3, switched to fwfmi and
    //       extended left over T1.
    void Two_Mismatch_Thirds() {
        Add({OP_STOP_IF_MAXHITS});
        Add({OP_RESUME, true, SEG_RH, 2, RL_HITS, RL_FWD, true});
        End_Stage(true, 2, RL_NONE);
        Add({OP_RESUME, false, SEG_LH, 2, RL_HITS, RL_BWD, true});
        End_Stage(false, 2, RL_NONE);
        const unsigned char T1 = Segment(PART_READ, {0, 1}, {1, 3}), T2 = Segment(PART_READ, {1, 3}, {2, 3}),
                            T3 = Segment(PART_READ, {2, 3}, {1, 1}), A_End = Segment(PART_READ, {1, 3}, {1, 1}),
                            Start_B = Segment(PART_READ, {0, 1}, {2, 3}), A_LH = Segment(PART_READ, {1, 3}, {1, 2}),
                            LH_B = Segment(PART_READ, {1, 2}, {2, 3});
        const unsigned char Win_A = Windows({{A_LH, 1, WIN_NARROW}, {SEG_RH, 1, WIN_NARROW}}),
                            Win_B = Windows({{T1, 1, WIN_NARROW}, {LH_B, 1, WIN_NARROW}}),
                            Win_C3 = Windows({{T3, 1, WIN_NARROW}}), Win_C1 = Windows({{T1, 1, WIN_NARROW}});
        Add({OP_STOP_IF_MAXHITS});
        Add({OP_SEARCH, true, T1, 0, 0, 0, false, false, Win_A, 0, 1});
        Add({OP_EXTEND, true, A_End, 0, RL_HITS, 0, true, false, Win_A});
        End_Stage(true, 2, RL_NONE);
        Add({OP_STOP_IF_MAXHITS});
        Add({OP_SEARCH, false, T3, 0, 0, 0, false, false, Win_B, 0, 2});
        Add({OP_EXTEND, false, Start_B, 0, RL_HITS, 0, true, false, Win_B});
        End_Stage(false, 2, RL_NONE);
        Add({OP_STOP_IF_MAXHITS});
        Add({OP_SEARCH, true, T2});
        Add({OP_EXTEND, true, T3, 0, RL_SEED_1, 0, false, false, Win_C3});
        Add({OP_CROSS, false, T1, 0, RL_HITS, RL_SEED_1, true, true, Win_C1, A_End});
        End_Stage(false, 2, RL_NONE);
    }

    // Round C's seed family: seeds of half Side with exactly S mismatches, crossing the other
    // half with Count C. RL_SEED_1 holds seeds in the crossing's index (revfmi for LH, which
    // crosses right; fwfmi for RH), RL_SEED_2 those switched to it first. Its rejects start a
    // new chain if a later round needs them.
    void Family(Mm_Part Side, int S, int C) {
        const bool D = Side == PART_LH;  // the crossing direction: right from LH, left from RH
        Seeds(Side, S);
        Add({OP_CROSS, D, D ? SEG_RH : SEG_LH, U(C), RL_HITS, RL_SEED_1, true, false});
        Add({OP_CROSS, D, D ? SEG_RH : SEG_LH, U(C), RL_HITS, RL_SEED_2, true, true, 0, D ? SEG_LH : SEG_RH});
        const unsigned char List = U(RL_FWD + Chains.size());
        if (C < K) {
            if (List >= RL_CHAINS_END) {
                std::fprintf(stderr, "-n %d needs more chain lists than RL_CHAINS_END allows\n", K);
                std::exit(1);
            }
            Chains.push_back({List, D});
        }
        End_Stage(D, C, C < K ? List : RL_NONE);
    }

    // Seeds of half Side with exactly S mismatches: the -n S search on the half, keeping its
    // exactly-S cases. The half's parts are fractions from its outer end (Mm_Segment_Spec),
    // so the RH cases are the LH ones mirrored. Quarters Q_out = [0, 1/2), Q_in = [1/2, 1);
    // thirds H_out = [0, 1/3), H_mid = [1/3, 2/3), H_in = [2/3, 1); p is the outer mismatch
    // (the first in LH, the last in RH), which alone the failure bound limits.
    //   S = 1, the -n 1 search on the half: Q_out exact, extended inward over Q_in; Q_in
    //          exact, extended outward over Q_out;
    //   S = 2, the -n 2 search on the half, as BatMis when both mismatches share a quarter:
    //          A, both in Q_out: Q_in exact, extended outward over Q_out;
    //          B, both in Q_in: Q_out exact (the half scan's), extended inward over Q_in;
    //          one in each, by the thirds (longer than BatMis's eighths):
    //          C1, p in [1/3, 1/2): H_out exact, extended inward over [1/3, 1);
    //          C2, p in H_out, the other in [1/2, 2/3): H_in exact, extended outward over [0, 2/3);
    //          C3, p in H_out, the other in H_in: H_mid exact, extended inward over H_in,
    //          switched and extended outward over H_out.
    void Seeds(Mm_Part Side, int S) {
        const bool D = Side == PART_LH;  // inward: right in LH, left in RH
        const unsigned char Slot = D ? 1 : 2;  // the half scan's kept outer quarter
        auto Seg = [&](int LoNum, int LoDen, int HiNum, int HiDen) {
            return Segment(Side, {U(LoNum), U(LoDen)}, {U(HiNum), U(HiDen)});
        };
        const unsigned char Q_Out = Seg(0, 1, 1, 2), Q_In = Seg(1, 2, 1, 1);
        if (S == 1) {
            const unsigned char W_In = Windows({{Q_In, 1, WIN_NARROW}}), W_Out = Windows({{Q_Out, 1, WIN_NARROW}});
            Add({OP_SEARCH, D, Q_Out, 0, 0, 0, false, false, W_In, 0, Slot});
            Add({OP_EXTEND, D, Q_In, 0, RL_SEED_1, 0, false, false, W_In});
            Add({OP_SEARCH, !D, Q_In});
            Add({OP_EXTEND, !D, Q_Out, 0, RL_SEED_2, 0, false, false, W_Out});
            return;
        }
        if (S == 3) {
            Three_Seeds(D, Slot, Seg, Q_Out, Q_In);
            return;
        }
        if (S != 2) {
            std::fprintf(stderr, "-n %d needs seeds with %d mismatches in a half (not built)\n", K, S);
            std::exit(1);
        }
        const unsigned char H_Out = Seg(0, 1, 1, 3), H_Mid = Seg(1, 3, 2, 3), H_In = Seg(2, 3, 1, 1);
        const unsigned char W_A = Windows({{Q_Out, 2, WIN_FREE}}), W_B = Windows({{Q_In, 2, WIN_AFTER_SEED}}),
                            W_C1 = Windows({{Seg(1, 3, 1, 2), 1, WIN_NARROW}, {Q_In, 1, WIN_FREE}}),
                            W_C2 = Windows({{H_Out, 1, WIN_NARROW}, {Seg(1, 2, 2, 3), 1, WIN_FREE}}),
                            W_H_In = Windows({{H_In, 1, WIN_FREE}}), W_H_Out = Windows({{H_Out, 1, WIN_NARROW}});
        Add({OP_SEARCH, !D, Q_In});  // A
        Add({OP_EXTEND, !D, Q_Out, 0, RL_SEED_2, 0, false, false, W_A});
        Add({OP_SEARCH, D, Q_Out, 0, 0, 0, false, false, W_B, 0, Slot});  // B
        Add({OP_EXTEND, D, Q_In, 0, RL_SEED_1, 0, false, false, W_B});
        Add({OP_SEARCH, D, H_Out, 0, 0, 0, false, false, W_C1});  // C1
        Add({OP_EXTEND, D, Seg(1, 3, 1, 1), 0, RL_SEED_1, 0, false, false, W_C1});
        Add({OP_SEARCH, !D, H_In, 0, 0, 0, false, false, W_C2});  // C2
        Add({OP_EXTEND, !D, Seg(0, 1, 2, 3), 0, RL_SEED_2, 0, false, false, W_C2});
        Add({OP_SEARCH, D, H_Mid, 0, 0, 0, false, false, W_H_Out});  // C3
        Add({OP_EXTEND, D, H_In, 0, RL_SEED_3, 0, false, false, W_H_In});
        Add({OP_CROSS, !D, H_Out, 0, RL_SEED_2, RL_SEED_3, false, true, W_H_Out, Seg(1, 3, 1, 1)});
    }

    // Seeds with exactly 3 mismatches in a half: the -n 3 search on it, by its quarters
    // (Q_out, Q_in) and eighths E1 = [0, 1/4), E2 = [1/4, 1/2), E3 = [1/2, 3/4), E4 = [3/4, 1),
    // from its outer end; p is the outer mismatch, r the inner one:
    //   (0, 3) in (Q_out, Q_in): Q_out exact (the half scan's), extended inward over Q_in;
    //   (3, 0): Q_in exact, extended outward over Q_out;
    //   (1, 2), p in E2: E1 exact, extended inward over [1/4, 1);
    //   (1, 2), p in E1: E2 exact, extended outward over E1, switched and extended inward
    //           over Q_in;
    //   (2, 1), r in E3: E4 exact, extended outward over [0, 3/4);
    //   (2, 1), r in E4: E3 exact, extended inward over E4, switched and extended outward
    //           over Q_out.
    // The smallest seed is an eighth of the half: three mismatches at its quarter points
    // leave no longer exact part (benchmark/seven_mismatch.md).
    template <class Seg_Of>
    void Three_Seeds(bool D, unsigned char Slot, Seg_Of& Seg, unsigned char Q_Out, unsigned char Q_In) {
        const unsigned char E1 = Seg(0, 1, 1, 4), E2 = Seg(1, 4, 1, 2), E3 = Seg(1, 2, 3, 4), E4 = Seg(3, 4, 1, 1);
        const unsigned char W_03 = Windows({{Q_In, 3, WIN_AFTER_SEED}}), W_30 = Windows({{Q_Out, 3, WIN_FREE}}),
                            W_12a = Windows({{E2, 1, WIN_NARROW}, {Q_In, 2, WIN_FREE}}),
                            W_E1 = Windows({{E1, 1, WIN_NARROW}}), W_Q_In2 = Windows({{Q_In, 2, WIN_FREE}}),
                            W_21a = Windows({{Q_Out, 2, WIN_FREE}, {E3, 1, WIN_FREE}}),
                            W_E4 = Windows({{E4, 1, WIN_FREE}}), W_Q_Out2 = Windows({{Q_Out, 2, WIN_FREE}});
        Add({OP_SEARCH, D, Q_Out, 0, 0, 0, false, false, W_03, 0, Slot});  // (0, 3)
        Add({OP_EXTEND, D, Q_In, 0, RL_SEED_1, 0, false, false, W_03});
        Add({OP_SEARCH, !D, Q_In});  // (3, 0)
        Add({OP_EXTEND, !D, Q_Out, 0, RL_SEED_2, 0, false, false, W_30});
        Add({OP_SEARCH, D, E1, 0, 0, 0, false, false, W_12a});  // (1, 2), p in E2
        Add({OP_EXTEND, D, Seg(1, 4, 1, 1), 0, RL_SEED_1, 0, false, false, W_12a});
        Add({OP_SEARCH, !D, E2, 0, 0, 0, false, false, W_E1});  // (1, 2), p in E1
        Add({OP_EXTEND, !D, E1, 0, RL_SEED_3, 0, false, false, W_E1});
        Add({OP_CROSS, D, Q_In, 0, RL_SEED_1, RL_SEED_3, false, true, W_Q_In2, Q_Out});
        Add({OP_SEARCH, !D, E4, 0, 0, 0, false, false, W_21a});  // (2, 1), r in E3
        Add({OP_EXTEND, !D, Seg(0, 1, 3, 4), 0, RL_SEED_2, 0, false, false, W_21a});
        Add({OP_SEARCH, D, E3});  // (2, 1), r in E4
        Add({OP_EXTEND, D, E4, 0, RL_SEED_4, 0, false, false, W_E4});
        Add({OP_CROSS, !D, Q_Out, 0, RL_SEED_2, RL_SEED_4, false, true, W_Q_Out2, Q_In});
    }
};

// The program and stage table of -n K (0 to MAX_MISMATCH_SEARCH), built on first use (thread-safe).
const Mm_Program& Program_For(int K) {
    static const std::vector<Mm_Program> Programs = [] {  // built once, before any thread reads them
        std::vector<Mm_Program> P;
        for (int k = 0; k <= MAX_MISMATCH_SEARCH; k++) P.push_back(Mm_Builder(k).Build());
        return P;
    }();
    return Programs[K];
}
inline const Mm_Step* Program_Of(int K) { return Program_For(K).Steps.data(); }
inline const Hit_Stage* Stages_Of(int K) { return Program_For(K).Stages.data(); }

#ifdef BATMAN_BATCH_STATS
// Lane occupancy and step counts (build with -DBATMAN_BATCH_STATS; printed with --timing).
// Per thread; Run_Batched adds the workers' to the main thread's (Add).
struct Batch_Stats_Counts {
    unsigned long Rounds = 0, Lane_Rounds = 0, Drain_Rounds = 0, Lanes = 0;
    unsigned long Steps[MM_PROGRAM_MAX][4][6] = {};  // [program step][task][sub]
    unsigned long Skipped[MM_PROGRAM_MAX] = {};      // searches skipped by their windows
    void Add(const Batch_Stats_Counts& O) {
        Rounds += O.Rounds;
        Lane_Rounds += O.Lane_Rounds;
        Drain_Rounds += O.Drain_Rounds;
        Lanes = std::max(Lanes, O.Lanes);
        for (unsigned p = 0; p < MM_PROGRAM_MAX; p++) {
            for (unsigned t = 0; t < 4; t++)
                for (unsigned u = 0; u < 6; u++) Steps[p][t][u] += O.Steps[p][t][u];
            Skipped[p] += O.Skipped[p];
        }
    }
};
thread_local Batch_Stats_Counts Batch_Stats;
#endif

// A read strand's state in a lane (two-strand programs): the lane's fields of the same
// names hold the active strand's, Mm_Set_Strand swaps them with these.
struct Mm_Strand {
    std::string_view All;
    unsigned Index;
    ForwardResults* Res;
    unsigned long ReadHits;
    int LH_Fail, RH_Fail;
    SARange Stored[2];
    bool Stored_Valid[2];
    std::size_t Tagged;
    uint64_t Packed[MAXTAG / 32 + 1];
};

struct Mm_Lane {
    std::string_view All;
    unsigned Index;
    ForwardResults* Res;  // the read's result lists (-n 1), or null
    unsigned long ReadHits;
    const Mm_Step* Program;
    unsigned PC;
    bool InTask, Done_Read;
    std::size_t Resume;   // OP_RESUME / OP_CROSS: next deferred range or seed
    bool Switching;       // the task is an index switch (Tag points into Patched)
    bool Locate_Other;    // MM_LOCATE: from the other index's row (Both_Rows)
    int LH_Fail, RH_Fail; // the exact LH / RH searches' failure bounds (Window_Set)
    SARange Stored[2];    // OP_SEARCH_KEEP's kept ranges, if Stored_Valid
    bool Stored_Valid[2];
    bool Keeping;         // OP_SEARCH_KEEP is searching its kept segment (else the rest)
    std::size_t Tagged;   // hits of Res[RL_HITS] tagged with their stage (End_Stage)
    // The current task.
    Mm_Task Task;
    bool Right, Found;
    ActiveBWT* Fmi;
    std::string_view Tag;
    int Offset, Count, Done;  // Done: bases matched when the sampled row was reached
    ForwardResults* Out;
    unsigned TaskMaxHits;
    std::vector<SARange> Stack;
    SARange Cur;
    Mm_Sub Sub;
    Mm_Windows W;            // TASK_WINDOW
    unsigned long Width;     // TASK_SWITCH: the range's width, the same in both indexes
    unsigned long Seed_Start;  // TASK_SWITCH: the seed's Start, the result's Other_Start
    unsigned long Seed_Other;  // TASK_SWITCH: the seed's Other_Start, where the walk must end
    uint64_t Packed[MAXTAG / 32 + 1];
    char Patched[MAXTAG];    // TASK_SWITCH: the read with the seed's mismatch applied
    // Two-strand programs: the active strand (the fields above) and the other one.
    bool Two_Strands;
    unsigned char Strand;
    Mm_Strand Other;
};

// Makes strand S the lane's active one (between steps only): swaps the per-strand fields.
inline void Mm_Set_Strand(Mm_Lane& L, unsigned S) {
    if (S == L.Strand || !L.Two_Strands) return;
    Mm_Strand& O = L.Other;
    const std::size_t Words = (L.All.size() + 31) / 32;  // both strands have the read's length
    std::swap(L.All, O.All);
    std::swap(L.Index, O.Index);
    std::swap(L.Res, O.Res);
    std::swap(L.ReadHits, O.ReadHits);
    std::swap(L.LH_Fail, O.LH_Fail);
    std::swap(L.RH_Fail, O.RH_Fail);
    std::swap(L.Stored, O.Stored);
    std::swap(L.Stored_Valid, O.Stored_Valid);
    std::swap(L.Tagged, O.Tagged);
    if (Use_Ref) std::swap_ranges(L.Packed, L.Packed + Words, O.Packed);
    L.Strand = static_cast<unsigned char>(S);
}

// The read base at Level and its read position, for the task's direction.
inline unsigned Mm_Base(const Mm_Lane& L, int Level) {
    return static_cast<unsigned char>(L.Right ? L.Tag[Level - 1] : L.Tag[L.Tag.size() - Level]);
}
inline int Mm_Read_Pos(const Mm_Lane& L, int Level) {
    return L.Right ? L.Offset + Level - 1 : L.Offset + static_cast<int>(L.Tag.size()) - Level;
}

// The top of an exact search's loop: a wide range steps both ends, a single row walks.
inline void Mm_Search_Enter(Mm_Lane& L) {
    SARange& T = L.Cur;
    if (T.End != T.Start && !T.Skip) {
        L.Sub = MM_WIDE;
        return;
    }
    L.Sub = MM_SINGLE;
    if (!T.Skip && (T.Start % SAINTERVAL == 0 || Other_Sampled(T))) {
        if (Use_Ref) {
            L.Done = T.Level - 1;
            L.Locate_Other = T.Start % SAINTERVAL != 0;
            L.Sub = MM_LOCATE;
            return;
        }
        Mark_Sampled(L.Fmi, T, L.Right, L.Offset, static_cast<int>(L.Tag.size()), T.Level - 1);
    }
}

// Ends a search task: Found, and (left searches) a zero Start on failure.
inline bool Mm_Search_End(Mm_Lane& L, bool Found) {
    L.Found = Found;
    if (!Found && !L.Right) L.Cur.Start = 0;
    return true;
}

// One step of a search task. Returns true when the task has finished.
inline bool Mm_Search_Step(Mm_Lane& L) {
    SARange& T = L.Cur;
    const int Len = static_cast<int>(L.Tag.size());
    switch (L.Sub) {
        case MM_WIDE:
            Get_SARange(Mm_Base(L, T.Level), T, L.Fmi);
            if (!T.Start) return Mm_Search_End(L, false);
            if (T.Level == Len) return Mm_Search_End(L, true);
            T.Level++;
            Mm_Search_Enter(L);
            return false;
        case MM_SINGLE: {
            if (T.Start == L.Fmi->inverseSa0) return Mm_Search_End(L, false);
            const unsigned long Index = T.Start - (T.Start > L.Fmi->inverseSa0);
            const unsigned Now = L.Fmi->bwtCode[Index / 16] << ((Index % 16) * 2) >> (BITS_IN_WORD - 2);
            const unsigned Expected = Mm_Base(L, T.Level);
            const unsigned Base = Read_Base(Expected);
            if (Expected != Now) return Mm_Search_End(L, false);
            T.Start = L.Fmi->cumulativeFreq[Base] + fm_occ::Occ_One(L.Fmi, T.Start, Base) + 1;
            if (!T.Skip && T.Start % SAINTERVAL == 0) {
                if (Use_Ref) {
                    L.Done = T.Level;
                    L.Locate_Other = false;
                    L.Sub = MM_LOCATE;
                    return false;
                }
                Mark_Sampled(L.Fmi, T, L.Right, L.Offset, Len, T.Level);
            }
            if (T.Level == Len) {
                if (L.Right || !T.Skip) T.End = T.Start;  // a left walk keeps a marked End
                return Mm_Search_End(L, true);
            }
            T.Level++;
            return false;
        }
        case MM_LOCATE:
            if (!Locate(L.Fmi, T, L.Right, L.Offset, Len, L.Done, L.Locate_Other)) {
                T.Level = Len;  // as Search_Located
                return Mm_Search_End(L, false);
            }
            L.Sub = MM_LOCATED;
            return false;
        case MM_LOCATED:
            if (!Located_Matches(T, L.Right, L.Offset, Len, L.Done, L.Packed)) return Mm_Search_End(L, false);
            if (L.Right) T.Level = Len;
            return Mm_Search_End(L, true);
        default:
            return true;
    }
}

// End of a singleton walk inside an extension, with its task's MaxHits check.
inline bool Mm_Single_Done(Mm_Lane& L) {
    L.Sub = MM_POP;
    return L.TaskMaxHits && L.Out->hitCount >= L.TaskMaxHits;
}

// A single-row node of an extension (L.Cur popped), up to its walk. Returns true when the
// task has finished; leaves Sub at MM_SINGLE to walk, or MM_LOCATE / MM_LOCATED.
inline bool Mm_Single_Enter(Mm_Lane& L) {
    SARange& T = L.Cur;
    if (T.Start == 0) return Mm_Single_Done(L);
    if (Use_Ref && T.Skip) {
        L.Sub = MM_LOCATED;
        return false;
    }
    L.Sub = MM_SINGLE;
    if (!T.Skip && (T.Start % SAINTERVAL == 0 || Other_Sampled(T))) {
        if (Use_Ref) {
            L.Locate_Other = T.Start % SAINTERVAL != 0;
            L.Sub = MM_LOCATE;
            return false;
        }
        Mark_Sampled(L.Fmi, T, L.Right, L.Offset, static_cast<int>(L.Tag.size()), T.Level - 1);
    }
    return false;
}

// One step of an extension task. Returns true when the task has finished.
inline bool Mm_Extend_Step(Mm_Lane& L) {
    const int Len = static_cast<int>(L.Tag.size());
    ForwardResults& Out = *L.Out;
    switch (L.Sub) {
        case MM_POP: {
            if (L.Stack.empty()) return true;
            SARange Range = L.Stack.back();
            L.Stack.pop_back();
            if (Range.End == Range.Start || Range.Skip) {
                L.Cur = Range;
                if (Mm_Single_Enter(L)) return true;
                if (L.Sub != MM_SINGLE) return false;
                return Mm_Extend_Step(L);  // its row was prefetched for this step
            }
            if (Range.Mismatches == L.Count) {  // at threshold: exact step only
                const unsigned char Expected = Mm_Base(L, Range.Level);
                SARange Temp = Range;
                Get_SARange(Expected, Temp, L.Fmi);
                if (!Temp.Start || Temp.End - Temp.Start != Range.End - Range.Start) {
                    Range.Forbidden_Base = Expected;
                    Out.mismatches.push_back(Range);
                }
                if (!Temp.Start) return false;
                Temp.Forbidden_Base = NO_BASE;
                if (Temp.Level == Len) {
                    Out.hits.push_back(Temp);
                    Out.hitCount += Temp.Skip ? 1 : Temp.End - Temp.Start + 1;
                    return L.TaskMaxHits && Out.hitCount >= L.TaskMaxHits;
                }
                Temp.Level++;
                L.Stack.push_back(Temp);
                return false;
            }
            Branch_Detect(Range, L.Fmi, L.Offset);
            for (int Branch = 0; Branch < 4; Branch++) {
                if (!Branch_Characters[Branch] || Branch == Range.Forbidden_Base) continue;
                SARange Temp = Range;
                Temp.Forbidden_Base = NO_BASE;
                Temp.Start = Branch_Ranges[Branch].Start;
                Temp.End = Branch_Ranges[Branch].End;
                Temp.Other_Start = Branch_Ranges[Branch].Other_Start;
                if (Mm_Base(L, Temp.Level) != static_cast<unsigned>(Branch)) {
                    Temp.Mismatch_Char = Temp.Mismatch_Char | (Branch << Temp.Mismatches * 2);
                    Temp.Mismatch_Pos[Temp.Mismatches] = Mm_Read_Pos(L, Temp.Level);
                    Temp.Mismatches++;
                }
                if (Temp.Mismatches > L.Count) continue;
                if (Temp.Level == Len) {
                    if (Temp.Mismatches == L.Count) {
                        Out.hits.push_back(Temp);
                        Out.hitCount += Temp.Skip ? 1 : Temp.End - Temp.Start + 1;
                        if (L.TaskMaxHits && Out.hitCount >= L.TaskMaxHits) return true;
                    }
                    continue;
                }
                Temp.Level++;
                L.Stack.push_back(Temp);
            }
            return false;
        }
        case MM_SINGLE: {
            SARange& T = L.Cur;
            if (T.Start == L.Fmi->inverseSa0) return Mm_Single_Done(L);
            const unsigned long Index = T.Start - (T.Start > L.Fmi->inverseSa0);
            const unsigned Now = L.Fmi->bwtCode[Index / 16] << ((Index % 16) * 2) >> (BITS_IN_WORD - 2);
            if (Now == T.Forbidden_Base) return Mm_Single_Done(L);
            const unsigned Expected = Mm_Base(L, T.Level);
            const unsigned Base = Read_Base(Expected);
            if (Now != Expected && T.Mismatches == L.Count) {  // at threshold: defer the parent
                if (!T.Skip) T.End = T.Start;
                T.Forbidden_Base = Expected;
                Out.mismatches.push_back(T);
                return Mm_Single_Done(L);
            }
            T.Forbidden_Base = NO_BASE;
            if (Now == Expected) T.Start = L.Fmi->cumulativeFreq[Base] + fm_occ::Occ_One(L.Fmi, T.Start, Base) + 1;
            else T.Start = L.Fmi->cumulativeFreq[Now] + fm_occ::Occ_One(L.Fmi, T.Start, Now) + 1;
            if (!T.Skip && T.Start % SAINTERVAL == 0) {
                Mark_Sampled(L.Fmi, T, L.Right, L.Offset, Len, T.Level);
            }
            if (Expected != Now) {
                T.Mismatch_Char = T.Mismatch_Char | (Now << T.Mismatches * 2);
                T.Mismatch_Pos[T.Mismatches] = Mm_Read_Pos(L, T.Level);
                T.Mismatches++;
            }
            if (Use_Ref && T.Skip) {  // just reached a sampled row
                T.Level++;
                L.Locate_Other = false;
                L.Sub = MM_LOCATE;
                return false;
            }
            if (T.Mismatches <= L.Count) {
                if (T.Level == Len) {
                    if (T.Mismatches == L.Count) {
                        if (!T.Skip) T.End = T.Start;
                        Out.hits.push_back(T);
                        Out.hitCount += T.Skip ? 1 : T.End - T.Start + 1;
                    }
                    return Mm_Single_Done(L);
                }
                T.Level++;
            }
            return false;
        }
        case MM_LOCATE:
            if (!Locate(L.Fmi, L.Cur, L.Right, L.Offset, Len, L.Cur.Level - 1, L.Locate_Other)) return Mm_Single_Done(L);
            L.Sub = MM_LOCATED;
            return false;
        case MM_LOCATED:
            Extend_Located(L.Cur, L.Right, L.Offset, Len, L.Count, Out, L.Packed);
            return Mm_Single_Done(L);
        default:
            return true;
    }
}

inline void Mm_Start_Search(Mm_Lane& L, bool Right, std::string_view Tag, int Offset, const SARange& R, ActiveBWT* Fmi) {
    L.Task = TASK_SEARCH;
    L.Right = Right;
    L.Tag = Tag;
    L.Offset = Offset;
    L.Fmi = Fmi;
    L.Cur = R;
    L.Cur.Anchor = Right ? Offset : Offset + static_cast<int>(Tag.size());  // the search's string starts (ends) with Tag
    L.Found = false;
    Mm_Search_Enter(L);
}

// Returns false if there is nothing to extend (no range).
inline bool Mm_Start_Extend(Mm_Lane& L, bool Right, std::string_view Tag, int Offset, const SARange& R, int Count,
                            ActiveBWT* Fmi, ForwardResults* Out, unsigned MaxHits) {
    if (!R.Start) return false;
    L.Task = TASK_EXTEND;
    L.Right = Right;
    L.Tag = Tag;
    L.Offset = Offset;
    L.Count = Count;
    L.Fmi = Fmi;
    L.Out = Out;
    L.TaskMaxHits = MaxHits;
    L.Stack.assign(1, R);
    L.Sub = MM_POP;
    return true;
}

// A windowed extension (1|1): as Mm_Start_Extend, with windows W instead of Count.
inline bool Mm_Start_Window(Mm_Lane& L, bool Right, std::string_view Tag, int Offset, const SARange& R,
                            const Mm_Windows& W, ForwardResults* Out, unsigned MaxHits) {
    if (!Mm_Start_Extend(L, Right, Tag, Offset, R, 0, Right ? revfmi : fwfmi, Out, MaxHits)) return false;
    L.Task = TASK_WINDOW;
    L.W = W;
    return true;
}

// A windowed node T that has matched its last base, or a child to extend further.
// Returns true when MaxHits stops the task.
inline bool Mm_Window_Push(Mm_Lane& L, SARange& T) {
    if (T.Level == static_cast<int>(L.Tag.size())) {
        L.Out->hits.push_back(T);
        L.Out->hitCount += T.Skip ? 1 : T.End - T.Start + 1;
        return L.TaskMaxHits && L.Out->hitCount >= L.TaskMaxHits;
    }
    T.Level++;
    L.Stack.push_back(T);
    return false;
}

inline void Mm_Record(SARange& T, unsigned Base, int Pos) {
    T.Mismatch_Char |= Base << T.Mismatches * 2;
    T.Mismatch_Pos[T.Mismatches] = Pos;
    T.Mismatches++;
}

// One step of a windowed extension: an extension's loop and singleton walk, where a
// mismatch is allowed only at a position whose window has room (Window_Can_Change), a
// node is dropped once a window can no longer reach its count (Window_Keep), and nothing
// is deferred. Returns true when the task has finished.
inline bool Mm_Window_Step(Mm_Lane& L) {
    const int Len = static_cast<int>(L.Tag.size());
    switch (L.Sub) {
        case MM_POP: {
            if (L.Stack.empty()) return true;
            SARange Range = L.Stack.back();
            L.Stack.pop_back();
            if (Range.End == Range.Start || Range.Skip) {
                L.Cur = Range;
                if (Mm_Single_Enter(L)) return true;
                if (L.Sub != MM_SINGLE) return false;
                return Mm_Window_Step(L);  // its row was prefetched for this step
            }
            const int x = Mm_Read_Pos(L, Range.Level);
            const unsigned Expected = Mm_Base(L, Range.Level);
            if (!Window_Can_Change(Range, L.W, x)) {  // exact step only
                SARange T = Range;
                Get_SARange(Expected, T, L.Fmi);
                if (!T.Start || !Window_Keep(T, L.W, x, L.Right)) return false;
                return Mm_Window_Push(L, T);
            }
            Branch_Detect(Range, L.Fmi, L.Offset);
            for (int Branch = 0; Branch < 4; Branch++) {
                if (!Branch_Characters[Branch]) continue;
                SARange T = Range;
                T.Start = Branch_Ranges[Branch].Start;
                T.End = Branch_Ranges[Branch].End;
                T.Other_Start = Branch_Ranges[Branch].Other_Start;
                if (static_cast<unsigned>(Branch) != Expected) Mm_Record(T, Branch, x);
                if (!Window_Keep(T, L.W, x, L.Right)) continue;
                if (Mm_Window_Push(L, T)) return true;
            }
            return false;
        }
        case MM_SINGLE: {
            SARange& T = L.Cur;
            if (T.Start == L.Fmi->inverseSa0) return Mm_Single_Done(L);
            const unsigned long Index = T.Start - (T.Start > L.Fmi->inverseSa0);
            const unsigned Now = L.Fmi->bwtCode[Index / 16] << ((Index % 16) * 2) >> (BITS_IN_WORD - 2);
            const unsigned Expected = Mm_Base(L, T.Level);
            const unsigned Base = Read_Base(Expected);
            const int x = Mm_Read_Pos(L, T.Level);
            if (Now != Expected && !Window_Can_Change(T, L.W, x)) return Mm_Single_Done(L);
            if (Now == Expected) T.Start = L.Fmi->cumulativeFreq[Base] + fm_occ::Occ_One(L.Fmi, T.Start, Base) + 1;
            else T.Start = L.Fmi->cumulativeFreq[Now] + fm_occ::Occ_One(L.Fmi, T.Start, Now) + 1;
            if (!T.Skip && T.Start % SAINTERVAL == 0) {
                Mark_Sampled(L.Fmi, T, L.Right, L.Offset, Len, T.Level);
            }
            if (Now != Expected) Mm_Record(T, Now, x);
            if (!Window_Keep(T, L.W, x, L.Right)) return Mm_Single_Done(L);
            if (Use_Ref && T.Skip) {  // just reached a sampled row
                T.Level++;
                L.Locate_Other = false;
                L.Sub = MM_LOCATE;
                return false;
            }
            if (T.Level == Len) {
                if (!T.Skip) T.End = T.Start;
                L.Out->hits.push_back(T);
                L.Out->hitCount++;
                return Mm_Single_Done(L);
            }
            T.Level++;
            return false;
        }
        case MM_LOCATE:
            if (!Locate(L.Fmi, L.Cur, L.Right, L.Offset, Len, L.Cur.Level - 1, L.Locate_Other)) return Mm_Single_Done(L);
            L.Sub = MM_LOCATED;
            return false;
        case MM_LOCATED: {
            SARange Hit = L.Cur;
            const int Done = Hit.Level - 1;
            if (L.Right) Record_Mismatches(Hit, true, L.Offset + Done, L.Offset + Len, L.Packed);
            else Record_Mismatches(Hit, false, L.Offset, L.Offset + Len - Done, L.Packed);
            if (Windows_Hold(Hit, L.W, L.Offset, L.Offset + Len)) {
                L.Out->hits.push_back(Hit);
                L.Out->hitCount++;
            }
            return Mm_Single_Done(L);
        }
        default:
            return true;
    }
}

// An index switch: the range of Seg of the read, with Seed's mismatches applied, in
// the other index (Right: revfmi, else fwfmi), for a seed found in one index and to be
// extended in the other. The string occurs (it is the seed's), so every range on the way
// is non-empty: Start alone is stepped from the lookup table (S' = C[c] + Occ(c, S) + 1,
// exact in that case) and End = Start + width - 1, the width being the same in both
// indexes. There is no located shortcut on the way: until the last base the range is
// wider than the seed's (shorter strings occur more often), so Start is not yet the
// seed's row, and without End it is not known when it becomes so.
inline void Mm_Start_Switch(Mm_Lane& L, const SARange& Seed, unsigned Seg, bool Right) {
    std::memcpy(L.Patched, L.All.data(), L.All.size());
    for (int i = 0; i < Seed.Mismatches; i++)  // one, or two for -n 4's 2|2 seeds
        L.Patched[Seed.Mismatch_Pos[i]] = static_cast<char>(Seed.Mismatch_Char >> 2 * i & 3);
    int Offset;
    const std::string_view S = Read_Segment_Of(L.All, Seg, Offset);
    L.Task = TASK_SWITCH;
    L.Right = Right;
    L.Fmi = Right ? revfmi : fwfmi;
    L.Tag = std::string_view(L.Patched + Offset, S.size());
    L.Offset = Offset;
    L.Width = Seed.Skip ? 1 : Seed.End - Seed.Start + 1;
    L.Seed_Start = Seed.Start;
    L.Seed_Other = Seed.Other_Start;
    L.Cur = Segment_Start(L.Tag, Right);
    L.Cur.Anchor = Right ? Offset : Offset + static_cast<int>(S.size());
    L.Cur.Mismatches = Seed.Mismatches;
    L.Cur.Mismatch_Char = Seed.Mismatch_Char;
    std::memcpy(L.Cur.Mismatch_Pos, Seed.Mismatch_Pos, sizeof L.Cur.Mismatch_Pos);
    L.Switching = true;
    L.Found = false;
    L.Sub = MM_SWITCH;
}

// The range Mm_Start_Switch's walk ends at, at once: the seed's string in the other index
// starts at the seed's Other_Start, with the seed's width, and its own Other_Start is the
// seed's Start (bidirectional ranges, checked by -DBATMIS2_VERIFY against the walk).
// Anchor: where the seed's string starts in the read (Right: the new range grows right)
// or ends.
inline SARange Switched_Range(const SARange& Seed, int Anchor) {
    SARange R{};
    R.Anchor = Anchor;
    R.Start = Seed.Other_Start;
    R.End = R.Start + (Seed.Skip ? 1 : Seed.End - Seed.Start + 1) - 1;
    R.Other_Start = Seed.Start;
    R.Mismatches = Seed.Mismatches;
    R.Mismatch_Char = Seed.Mismatch_Char;
    std::memcpy(R.Mismatch_Pos, Seed.Mismatch_Pos, sizeof R.Mismatch_Pos);
    return R;
}

// One step of an index switch. Returns true when it has finished (Found, range in Cur).
inline bool Mm_Switch_Step(Mm_Lane& L) {
    SARange& T = L.Cur;
    const int Len = static_cast<int>(L.Tag.size());
    if (!T.Start) return true;
    if (T.Level > Len) {
        T.End = T.Start + L.Width - 1;
        T.Other_Start = L.Seed_Start;  // the same string, in the seed's index
#ifdef BATMIS2_VERIFY
        Verify_Other(L.Fmi, T.Start, T.End, T.Other_Start);
        if (T.Start != L.Seed_Other) {
            if (Other_Failed++ < 5) std::fprintf(stderr, "Verify_Other: switch walk ends at %lu, seed's Other_Start %lu\n", T.Start, L.Seed_Other);
        }
        Other_Checked++;
#endif
        L.Found = true;
        return true;
    }
    const unsigned c = Mm_Base(L, T.Level);
    T.Start = L.Fmi->cumulativeFreq[c] + fm_occ::Occ_One(L.Fmi, T.Start, c) + 1;
    T.Level++;
    return false;
}

// A segment of the lane's read and its offset in the read.
inline std::string_view Mm_Segment(const Mm_Lane& L, unsigned Seg, int& Offset) {
    return Read_Segment_Of(L.All, Seg, Offset);
}

// Hits the read has so far: its exact count plus its mismatch hits (RL_HITS), unless
// Skip is RL_HITS (an extension into RL_HITS counts that list itself); with two strands,
// the other strand's too.
inline unsigned long Mm_Used(const Mm_Lane& L, int Skip = -1) {
    unsigned long Used = L.ReadHits + (L.Res && Skip != static_cast<int>(RL_HITS) ? L.Res[RL_HITS].hitCount : 0);
    if (L.Two_Strands) Used += L.Other.ReadHits + (L.Other.Res ? L.Other.Res[RL_HITS].hitCount : 0);
    return Used;
}

// OP_CROSS: extends seed R (its mismatches kept) over the step's segment. Returns false
// if there is nothing to extend.
inline bool Mm_Cross_Extend(Mm_Lane& L, const Mm_Step& P, SARange R, unsigned MaxHits) {
    int Offset;
    const std::string_view Seg = Mm_Segment(L, P.Seg, Offset);
    R.Level = 1;
    R.Forbidden_Base = NO_BASE;
    const unsigned Budget = P.Limited && MaxHits ? MaxHits - Mm_Used(L, P.Out) : 0;
    if (P.Win)
        return Mm_Start_Window(L, P.Right, Seg, Offset, R, Window_Set(P.Win, static_cast<int>(L.All.size()), L.LH_Fail, L.RH_Fail),
                               &L.Res[P.Out], Budget);
    return Mm_Start_Extend(L, P.Right, Seg, Offset, R, P.Count, P.Right ? revfmi : fwfmi, &L.Res[P.Out], Budget);
}

// A failed exact search of LH or RH: the first Depth bases of the segment (in the search
// direction) do not occur, which bounds the 1|1 windows (Window_Set).
inline void Mm_Search_Failed(Mm_Lane& L, const Mm_Step& P, int Depth) {
    if (P.Seg == SEG_LH) L.LH_Fail = Depth;
    else if (P.Seg == SEG_RH) L.RH_Fail = static_cast<int>(L.All.size()) - Depth;
}

// OP_SEARCH_KEEP: L.Cur has matched the kept segment, its first Stop bases. Keeps it in
// the step's slot and starts the search of the rest of Seg.
inline void Mm_Keep_Rest(Mm_Lane& L, const Mm_Step& P, std::string_view Seg, int Offset, int Stop) {
    L.Stored[P.Stored - 1] = L.Cur;
    L.Stored_Valid[P.Stored - 1] = true;
    L.Keeping = false;
    SARange R = L.Cur;
    const bool Located = Use_Ref && R.Skip;
    if (!Located) R.Level = Stop + 1;
    Mm_Start_Search(L, P.Right, Seg, Offset, R, P.Right ? revfmi : fwfmi);
    if (Located) {  // compare the rest with the genome
        L.Done = Stop;
        L.Sub = MM_LOCATED;
    }
}

// Runs the lane's program until it starts a task or ends (Done_Read). Called when a
// read starts and whenever its task finishes.
void Mm_Run(Mm_Lane& L, unsigned MaxHits) {
    for (;;) {
        const Mm_Step& P = L.Program[L.PC];
        Mm_Set_Strand(L, P.Strand);
        int Offset;
        switch (P.Op) {
            case OP_SEARCH: {
                if (L.InTask) {  // the search has finished; Level is the failure depth if it failed
                    L.InTask = false;
                    if (!L.Found) Mm_Search_Failed(L, P, L.Cur.Level);
                    L.PC++;
                    continue;
                }
                if (P.Win && !Windows_Possible(Window_Set(P.Win, static_cast<int>(L.All.size()), L.LH_Fail, L.RH_Fail))) {
#ifdef BATMAN_BATCH_STATS
                    Batch_Stats.Skipped[L.PC]++;
#endif
                    L.Found = false;  // the case has no hits: skip the search and its extension
                    L.PC++;
                    continue;
                }
                if (P.Stored && L.Stored_Valid[P.Stored - 1]) {
                    L.Cur = L.Stored[P.Stored - 1];
                    L.Found = true;
                    L.PC++;
                    continue;
                }
                const std::string_view Seg = Mm_Segment(L, P.Seg, Offset);
                SARange R = Segment_Start(Seg, P.Right);
                R.Anchor = P.Right ? Offset : Offset + static_cast<int>(Seg.size());  // also for a range taken whole
                if (R.Level > static_cast<int>(Seg.size()) || !R.Start) {  // no longer than the key, or no match
                    L.Cur = R;
                    L.Found = R.Start != 0;
                    if (!L.Found) Mm_Search_Failed(L, P, R.Level - 1);
                    L.PC++;
                    continue;
                }
                Mm_Start_Search(L, P.Right, Seg, Offset, R, P.Right ? revfmi : fwfmi);
                L.InTask = true;
                return;
            }
            case OP_SEARCH_KEEP: {
                const std::string_view Seg = Mm_Segment(L, P.Seg, Offset);
                int Part_Offset;
                const std::string_view Part = Mm_Segment(L, P.Seg2, Part_Offset);
                const int Stop = static_cast<int>(Part.size());
                ActiveBWT* Fmi = P.Right ? revfmi : fwfmi;
                if (L.InTask) {
                    L.InTask = false;
                    if (L.Keeping && L.Found) {  // the kept segment matched: keep it, search the rest
                        Mm_Keep_Rest(L, P, Seg, Offset, Stop);
                        L.InTask = true;
                        return;
                    }
                    if (!L.Found) Mm_Search_Failed(L, P, L.Cur.Level);
                    L.PC++;
                    continue;
                }
                L.Stored_Valid[P.Stored - 1] = false;
                SARange R = Segment_Start(Seg, P.Right);
                R.Anchor = P.Right ? Offset : Offset + static_cast<int>(Seg.size());  // also kept as is (Mm_Keep_Rest)
                if (R.Level > static_cast<int>(Seg.size()) || !R.Start) {  // as OP_SEARCH
                    L.Cur = R;
                    L.Found = R.Start != 0;
                    if (!L.Found) Mm_Search_Failed(L, P, R.Level - 1);
                    L.PC++;
                    continue;
                }
                L.InTask = true;
                if (R.Level - 1 > Stop) {  // the lookup table covers the kept segment: a plain search
                    L.Keeping = false;
                    Mm_Start_Search(L, P.Right, Seg, Offset, R, Fmi);
                } else if (R.Level - 1 == Stop) {
                    L.Cur = R;
                    Mm_Keep_Rest(L, P, Seg, Offset, Stop);
                } else {
                    L.Keeping = true;
                    Mm_Start_Search(L, P.Right, Part, Part_Offset, R, Fmi);
                }
                return;
            }
            case OP_EXTEND: {
                if (L.InTask) {
                    L.InTask = false;
                    L.PC++;
                    continue;
                }
                if (!L.Found || !L.Cur.Start || (P.Limited && MaxHits && Mm_Used(L) >= MaxHits)) {
                    L.PC++;  // nothing to extend, or no hits left to find (a budget of 0 is no limit)
                    continue;
                }
                SARange R = L.Cur;
                R.Level = 1;
                R.Mismatches = 0;
                R.Mismatch_Char = 0;
                const std::string_view Seg = Mm_Segment(L, P.Seg, Offset);
                const unsigned Budget = P.Limited && MaxHits ? MaxHits - Mm_Used(L, P.Out) : 0;
                const bool Started =
                    P.Win ? Mm_Start_Window(L, P.Right, Seg, Offset, R, Window_Set(P.Win, static_cast<int>(L.All.size()), L.LH_Fail, L.RH_Fail),
                                            &L.Res[P.Out], Budget)
                          : Mm_Start_Extend(L, P.Right, Seg, Offset, R, P.Count, P.Right ? revfmi : fwfmi, &L.Res[P.Out], Budget);
                if (Started) {
                    L.InTask = true;
                    return;
                }
                L.PC++;
                continue;
            }
            case OP_COUNT_EXACT:
                L.ReadHits = 0;
                for (const SARange& hit : L.Res[P.From].hits) {
                    unsigned long hitCount = hit.Skip ? 1UL : (hit.End >= hit.Start ? (hit.End - hit.Start + 1UL) : 0UL);
                    if (MaxHits > 0 && hitCount > MaxHits) hitCount = MaxHits;
                    L.ReadHits += hitCount;
                }
                L.PC++;
                continue;
            case OP_RESUME: {
                L.InTask = false;
                const std::vector<SARange>& Deferred = L.Res[P.From].mismatches;
                const std::string_view Seg = Mm_Segment(L, P.Seg, Offset);
                while (L.Resume < Deferred.size()) {
                    if (MaxHits && Mm_Used(L) >= MaxHits) break;
                    const SARange Parent = Deferred[L.Resume++];
                    const unsigned Budget = P.Limited && MaxHits ? MaxHits - Mm_Used(L, P.Out) : 0;
                    if (Mm_Start_Extend(L, P.Right, Seg, Offset, Parent, P.Count, P.Right ? revfmi : fwfmi,
                                        &L.Res[P.Out], Budget)) {
                        L.InTask = true;
                        return;
                    }
                }
                L.Resume = 0;
                L.PC++;
                continue;
            }
            case OP_CROSS: {
                const std::vector<SARange>& Seeds = L.Res[P.From].hits;
                if (L.InTask) {
                    L.InTask = false;
                    if (L.Task == TASK_SWITCH) {  // the seed's switch has finished
                        L.Switching = false;
                        if (L.Found && Mm_Cross_Extend(L, P, L.Cur, MaxHits)) {
                            L.InTask = true;
                            return;
                        }
                    }
                }
                while (L.Resume < Seeds.size()) {
                    if (MaxHits && Mm_Used(L) >= MaxHits) {
                        L.Done_Read = true;
                        return;
                    }
                    const SARange& Seed = Seeds[L.Resume++];
                    if (P.Switch && !(Use_Ref && Seed.Skip)) {
                        if (!Switch_Walk) {  // the switched range at once, then the crossing
                            int Off;
                            const std::string_view S2 = Read_Segment_Of(L.All, P.Seg2, Off);
                            const int Anchor = P.Right ? Off : Off + static_cast<int>(S2.size());
                            if (Mm_Cross_Extend(L, P, Switched_Range(Seed, Anchor), MaxHits)) {
                                L.InTask = true;
                                return;
                            }
                            continue;
                        }
                        Mm_Start_Switch(L, Seed, P.Seg2, P.Right);
                        L.InTask = true;
                        return;
                    }
                    if (Mm_Cross_Extend(L, P, Seed, MaxHits)) {
                        L.InTask = true;
                        return;
                    }
                }
                L.Resume = 0;
                L.PC++;
                continue;
            }
            case OP_END_STAGE:
                End_Stage(L.Res, L.Tagged, P.Count, P.Out);
                L.PC++;
                continue;
            case OP_STOP_IF_MAXHITS:
                if (MaxHits && Mm_Used(L) >= MaxHits) {
                    L.Done_Read = true;
                    return;
                }
                L.PC++;
                continue;
            case OP_END:
                L.Done_Read = true;
                return;
        }
    }
}

PREFETCH_ONLY void Mm_Prefetch(const Mm_Lane& L) {
    const SARange* T = &L.Cur;
    switch (L.Sub) {
        case MM_POP:
            if (L.Stack.empty()) return;
            T = &L.Stack.back();
            if (Use_Ref && T->Skip) {
                __builtin_prefetch(Pac.data() + (T->Pos + L.Offset) / 4);
                return;
            }
            [[fallthrough]];
        case MM_WIDE:
            Prefetch_Row(L.Fmi, T->Start);
            if (!T->Skip && T->End != T->Start) Prefetch_Row(L.Fmi, T->End + 1);
            return;
        case MM_SINGLE:
            Prefetch_Row(L.Fmi, T->Start);
            return;
        case MM_LOCATE:
            if (L.Locate_Other) __builtin_prefetch((L.Right ? fwfmi : revfmi)->saValue + T->Other_Start / SAINTERVAL);
            else __builtin_prefetch(L.Fmi->saValue + T->Start / SAINTERVAL);
            return;
        case MM_LOCATED:
            __builtin_prefetch(Pac.data() + (T->Pos + L.Offset) / 4);
            __builtin_prefetch(Pac.data() + (T->Pos + L.Offset + L.Tag.size()) / 4);
            return;
        case MM_SWITCH:
            Prefetch_Row(L.Fmi, T->Start);
            return;
    }
}

// One step of the lane's task. Returns true when the task has finished.
inline bool Mm_Step_Task(Mm_Lane& L) {
    switch (L.Task) {
        case TASK_SEARCH: return Mm_Search_Step(L);
        case TASK_EXTEND: return Mm_Extend_Step(L);
        case TASK_WINDOW: return Mm_Window_Step(L);
        default: return Mm_Switch_Step(L);
    }
}

// Runs Program on Count reads of Len bases (codes, back to back), B at a time.
// -n 0 (Program_Of(0)): Ranges[i] / Found[i] as Search_Right from Start_Range.
// -n 1 to 5 (Program_Of(k)): Res[MM_LISTS * i ...] get read i's result lists
// (Result_List), ReadHits[i] its exact hit count.
// With Two_Strands (a Zig_Zag program), Reads holds 2 * Count reads, read Count + i the
// reverse complement of read i, and one lane searches both (Mm_Strand); Res and ReadHits
// are indexed by read as before. Only reads [First, Last) are searched (Batch_Chooser).
void Batch_Search(const char* Reads, unsigned Count, int Len, unsigned B, unsigned MaxHits, const Mm_Step* Program,
                  std::vector<ForwardResults>* Res, std::vector<unsigned long>* ReadHits, std::vector<SARange>* Ranges,
                  std::vector<unsigned char>* Found, bool Two_Strands = false, unsigned First = 0,
                  unsigned Last = ~0U) {
    static thread_local Mm_Lane Lanes[MAX_BATCH];
    unsigned Active = 0, Next = First;
    Last = std::min(Last, Count);
    auto Finish_Strand = [&](Mm_Lane& L) {
        if (L.Res) {  // a read stopped inside a stage (OP_CROSS at MaxHits): tag its last hits
            for (const Mm_Step* P = Program + L.PC; P->Op != OP_END; P++)
                if (P->Op == OP_END_STAGE && P->Strand == L.Strand) {
                    Tag_Hits(L.Res, L.Tagged, P->Count);
                    break;
                }
        }
        if (ReadHits) (*ReadHits)[L.Index] = L.ReadHits;
        if (Ranges) {
            (*Ranges)[L.Index] = L.Cur;
            (*Found)[L.Index] = L.Found;
        }
    };
    auto Finish = [&](Mm_Lane& L) {
        Finish_Strand(L);
        if (!L.Two_Strands) return;
        Mm_Set_Strand(L, 1 - L.Strand);
        Finish_Strand(L);
    };
    // Starts a strand: read I's lists, packed read and search state (the active fields).
    auto Start_Strand = [&](Mm_Lane& L, unsigned I) {
        L.All = std::string_view(Reads + static_cast<std::size_t>(I) * Len, Len);
        L.Index = I;
        L.Res = Res ? &(*Res)[MM_LISTS * static_cast<std::size_t>(I)] : nullptr;
        L.Tagged = 0;
        if (L.Res)
            for (unsigned k = 0; k < MM_LISTS; k++) {
                L.Res[k].hits.clear();
                L.Res[k].mismatches.clear();
                L.Res[k].hitCount = 0;
            }
        L.ReadHits = 0;
        L.LH_Fail = Len;
        L.RH_Fail = 0;
        L.Stored_Valid[0] = L.Stored_Valid[1] = false;
        if (Use_Ref) Pack_Read(L.All, L.Packed);
    };
    auto Fill = [&](Mm_Lane& L) {
        while (Next < Last) {
            const unsigned I = Next++;
            L.Two_Strands = false;
            L.Strand = 0;
            if (Two_Strands) {  // the reverse complement first, then swapped out
                Start_Strand(L, Count + I);
                L.Two_Strands = true;
                L.Strand = 1;
                Mm_Set_Strand(L, 0);
            }
            Start_Strand(L, I);
            L.Program = Program;
            L.PC = 0;
            L.InTask = L.Done_Read = false;
            L.Resume = 0;
            L.Switching = false;
            L.Found = false;
            Mm_Run(L, MaxHits);
            if (!L.Done_Read) return true;
            Finish(L);
        }
        return false;
    };
    for (unsigned i = 0; i < B && Fill(Lanes[Active]); i++) Active++;
    while (Active) {
#ifdef BATMAN_BATCH_STATS
        Batch_Stats.Rounds++;
        Batch_Stats.Lane_Rounds += Active;
        Batch_Stats.Lanes = B;
        Batch_Stats.Drain_Rounds += Active < B;
#endif
        for (unsigned i = 0; i < Active; i++) Mm_Prefetch(Lanes[i]);
        for (unsigned i = 0; i < Active;) {
            Mm_Lane& L = Lanes[i];
#ifdef BATMAN_BATCH_STATS
            Batch_Stats.Steps[L.PC][L.Task][L.Sub]++;
#endif
            if (Mm_Step_Task(L)) Mm_Run(L, MaxHits);
            if (!L.Done_Read) {
                i++;
                continue;
            }
            Finish(L);
            if (Fill(L)) {
                i++;
                continue;
            }
            std::swap(L, Lanes[--Active]);  // no reads left: close the gap (keeps the stacks' memory)
            if (L.Switching) L.Tag = std::string_view(L.Patched + L.Offset, L.Tag.size());  // the moved lane's own copy
        }
    }
}

void Batch_Exact(const char* Reads, unsigned Count, int Len, unsigned B, std::vector<SARange>& Ranges,
                 std::vector<unsigned char>& Found) {
    Batch_Search(Reads, Count, Len, B, 0, Program_Of(0), nullptr, nullptr, &Ranges, &Found);
}

void Batch_One_Mismatch(const char* Reads, unsigned Count, int Len, unsigned B, unsigned MaxHits,
                        std::vector<ForwardResults>& Res, std::vector<unsigned long>& ReadHits) {
    Batch_Search(Reads, Count, Len, B, MaxHits, Program_Of(1), &Res, &ReadHits, nullptr, nullptr);
}

void Batch_Two_Mismatch(const char* Reads, unsigned Count, int Len, unsigned B, unsigned MaxHits,
                        std::vector<ForwardResults>& Res, std::vector<unsigned long>& ReadHits) {
    Batch_Search(Reads, Count, Len, B, MaxHits, Program_Of(2), &Res, &ReadHits, nullptr, nullptr);
}

void Batch_Three_Mismatch(const char* Reads, unsigned Count, int Len, unsigned B, unsigned MaxHits,
                          std::vector<ForwardResults>& Res, std::vector<unsigned long>& ReadHits) {
    Batch_Search(Reads, Count, Len, B, MaxHits, Program_Of(3), &Res, &ReadHits, nullptr, nullptr);
}

void Batch_Four_Mismatch(const char* Reads, unsigned Count, int Len, unsigned B, unsigned MaxHits,
                         std::vector<ForwardResults>& Res, std::vector<unsigned long>& ReadHits) {
    Batch_Search(Reads, Count, Len, B, MaxHits, Program_Of(4), &Res, &ReadHits, nullptr, nullptr);
}

void Batch_Five_Mismatch(const char* Reads, unsigned Count, int Len, unsigned B, unsigned MaxHits,
                         std::vector<ForwardResults>& Res, std::vector<unsigned long>& ReadHits) {
    Batch_Search(Reads, Count, Len, B, MaxHits, Program_Of(5), &Res, &ReadHits, nullptr, nullptr);
}

// The two-strand (zig-zag) form of a one-strand program, for Batch_Search with
// Two_Strands. P is cut into groups at each OP_STOP_IF_MAXHITS, and after OP_COUNT_EXACT
// so that the exact hits are a group of their own. Each group runs on the read, then on
// its reverse complement, with a MaxHits stop after each. Each strand's stages keep their
// order, and a group holds one mismatch count (the first: exact hits, never budget-
// limited), so the hits found are still the best: a strand stops only once the hits
// reach MaxHits with nothing better left in either strand, and the reverse complement's
// half of a group could only tie (ties: the read's strand first).
std::vector<Mm_Step> Zig_Zag(const Mm_Step* P) {
    std::vector<Mm_Step> Out, Group;
    auto Flush = [&] {
        if (Group.empty()) return;
        for (unsigned char Strand : {0, 1}) {
            for (Mm_Step Step : Group) {
                Step.Strand = Strand;
                Out.push_back(Step);
            }
            Out.push_back({OP_STOP_IF_MAXHITS, false, 0, 0, 0, 0, false, false});
        }
        Group.clear();
    };
    for (; P->Op != OP_END; P++) {
        if (P->Op == OP_STOP_IF_MAXHITS) {
            Flush();
            continue;
        }
        Group.push_back(*P);
        if (P->Op == OP_COUNT_EXACT) Flush();
    }
    Flush();
    Out.push_back({OP_END, false, 0, 0, 0, 0, false, false});
    if (Out.size() > MM_PROGRAM_MAX) {
        std::fprintf(stderr, "zig-zag program too long (%zu steps)\n", Out.size());
        std::exit(1);
    }
    return Out;
}

// -n 1 to 7 on both strands of each read (Zig_Zag): Reads holds 2 * Count reads, read
// Count + i the reverse complement of read i; Res and ReadHits are indexed the same way.
// The zig-zag program of -n K (1 to MAX_MISMATCH_SEARCH).
const Mm_Step* Zig_Zag_Program(int K) {
    static const std::vector<std::vector<Mm_Step>> Programs = [] {
        std::vector<std::vector<Mm_Step>> P;
        for (int k = 1; k <= MAX_MISMATCH_SEARCH; k++) P.push_back(Zig_Zag(Program_Of(k)));
        return P;
    }();
    return Programs[K - 1].data();
}

void Batch_Both_Strands(const char* Reads, unsigned Count, int Len, unsigned B, unsigned MaxHits, int MaxMismatches,
                        std::vector<ForwardResults>& Res, std::vector<unsigned long>& ReadHits) {
    Batch_Search(Reads, Count, Len, B, MaxHits, Zig_Zag_Program(MaxMismatches), &Res, &ReadHits, nullptr, nullptr, true);
}

// --batch auto: the lane count chosen on the fly, and a profile of the candidates printed
// at the end (Auto_Profile). The best count depends on the machine (how many cache misses
// a core keeps in flight), the index and the reads, while the output does not depend on
// it. With working prefetches the time is flat from about 8 lanes up
// (benchmark/prefetch_fix.md), so a fixed count is the default (DEFAULT_BATCH, or
// $BATMIS_BATCH) and this is for checking that default on a new machine or workload.
// The first chunks try the candidates, one chunk each,
// in rounds of alternating order (cancelling a drift such as caches warming up); each
// chunk's first reads only warm up, so that its time is the candidate's steady state
// (switching the lane count every few hundred reads favours many lanes). The candidate
// with the lowest median search time per read is kept for the rest of the run, or,
// within 3% of it, the one with the fewest lanes. With threads (-t), each thread has its
// own chooser, timed while the others run: they share the cores' and the memory's
// capacity for misses in flight, and a thread may run on a P- or an E-core.
struct Batch_Chooser {
    static constexpr unsigned Candidates[] = {4, 8, 16};
    static constexpr unsigned N = sizeof Candidates / sizeof Candidates[0];
    static constexpr unsigned Warm_Up = 256, Rounds = 3;
    unsigned Chosen = 0;  // 0: still trying
    unsigned Trials = 0;
    std::vector<double> Seconds[N];  // per read, one per trial block

    // Searches a block's reads [First, Last) (one thread's chunk) with
    // Search_Range(First, Last, Lanes).
    template <class Search>
    void Run(unsigned First, unsigned Last, Search&& Search_Range) {
        if (Chosen) {
            Search_Range(First, Last, Chosen);
            return;
        }
        const unsigned Round = Trials / N, C = Round % 2 ? N - 1 - Trials % N : Trials % N;
        if (Last - First < 2 * Warm_Up) {  // a short (last) chunk: too noisy to count
            Search_Range(First, Last, Candidates[C]);
            return;
        }
        Search_Range(First, First + Warm_Up, Candidates[C]);
        const auto Start = std::chrono::steady_clock::now();
        Search_Range(First + Warm_Up, Last, Candidates[C]);
        Seconds[C].push_back(std::chrono::duration<double>(std::chrono::steady_clock::now() - Start).count() /
                             (Last - First - Warm_Up));
        if (++Trials == N * Rounds) Choose();
    }

    // Each candidate's median seconds per read over its trial chunks so far; false if
    // one has none yet.
    bool Medians(double Median[N]) const {
        for (unsigned c = 0; c < N; c++) {
            std::vector<double> S = Seconds[c];
            if (S.empty()) return false;
            std::nth_element(S.begin(), S.begin() + S.size() / 2, S.end());
            Median[c] = S[S.size() / 2];
        }
        return true;
    }

    void Choose() {
        double Median[N], Best = 0;
        Medians(Median);
        for (unsigned c = 0; c < N; c++)
            if (!c || Median[c] < Best) Best = Median[c];
        for (unsigned c = 0; c < N && !Chosen; c++)
            if (Median[c] <= Best * 1.03) Chosen = Candidates[c];
        Report = "batch=auto:" + std::to_string(Chosen) + " (ns/read";
        for (unsigned c = 0; c < N; c++)
            Report += " " + std::to_string(Candidates[c]) + ":" + std::to_string(static_cast<long>(Median[c] * 1e9));
        Report += ")";
    }

    // For --timing; "batch=auto:-" if the input ended before the trial did.
    std::string Report = "batch=auto:-";
};

constexpr unsigned DEFAULT_BATCH = 16;  // --batch, when neither it nor $BATMIS_BATCH is given

// --batch auto's profile, one line for stderr: each candidate's search time per read (with
// threads, the median over the threads that timed it), what was chosen, and how far the
// default lane count is from the fastest candidate. A few runs of it on a machine show
// whether the default suits it; if not, set --batch or $BATMIS_BATCH.
std::string Auto_Profile(const std::vector<Batch_Chooser>& Choosers) {
    constexpr unsigned N = Batch_Chooser::N;
    std::vector<double> Times[N];
    unsigned Chose[N] = {}, Finished = 0;
    for (const Batch_Chooser& C : Choosers) {
        double Median[N];
        if (C.Medians(Median))
            for (unsigned c = 0; c < N; c++) Times[c].push_back(Median[c]);
        for (unsigned c = 0; c < N; c++) Chose[c] += C.Chosen == Batch_Chooser::Candidates[c];
        Finished += C.Chosen != 0;
    }
    const unsigned Threads = static_cast<unsigned>(Choosers.size());
    const unsigned long Needed = static_cast<unsigned long>(N) * Batch_Chooser::Rounds;
    std::string Out = "batch auto: ";
    if (Times[0].empty())
        return Out + "too few reads to time the candidates (the trial takes " + std::to_string(Needed) +
               " chunks per thread)";
    double Typical[N];
    unsigned Fastest = 0, Default = N;
    Out += "ns/read";
    for (unsigned c = 0; c < N; c++) {
        std::nth_element(Times[c].begin(), Times[c].begin() + Times[c].size() / 2, Times[c].end());
        Typical[c] = Times[c][Times[c].size() / 2];
        if (Typical[c] < Typical[Fastest]) Fastest = c;
        if (Batch_Chooser::Candidates[c] == DEFAULT_BATCH) Default = c;
        Out += " " + std::to_string(Batch_Chooser::Candidates[c]) + ":" + std::to_string(static_cast<long>(Typical[c] * 1e9));
    }
    if (Threads == 1) {
        Out += Finished ? "; chosen " + std::to_string(Choosers[0].Chosen) : std::string("; trial unfinished");
    } else {
        Out += "; chosen";
        for (unsigned c = 0; c < N; c++) Out += " " + std::to_string(Batch_Chooser::Candidates[c]) + ":" + std::to_string(Chose[c]);
        Out += " threads";
        if (Finished < Threads) Out += " (" + std::to_string(Threads - Finished) + " did not finish the trial)";
    }
    if (Default < N) {
        char Text[96];
        std::snprintf(Text, sizeof Text, "; default %u is %+.1f%% against the fastest (%u)", DEFAULT_BATCH,
                      (Typical[Default] / Typical[Fastest] - 1) * 100, Batch_Chooser::Candidates[Fastest]);
        Out += Text;
    }
    return Out;
}

// "\t<dir>\t<start>\t<end>" for an SA range, or "\tP\t<pos>\t<pos>" for a located hit.
std::string Hit_Text(const SARange& hit, const char* dir) {
    if (Use_Ref && hit.Skip) {
        const std::string pos = std::to_string(hit.Pos);
        return "\tP\t" + pos + "\t" + pos;
    }
    return std::string("\t") + dir + "\t" + std::to_string(hit.Start) + "\t" + std::to_string(hit.Skip ? hit.Start : hit.End);
}

// The reads file, plain or gzip-compressed (zlib reads either; "-" is stdin), line by
// line. With Inflate_Ahead, a thread of its own reads and inflates the file into a ring
// of buffers while Get_Line splits the lines of the one before; otherwise Get_Line reads
// each buffer itself.
class Read_File {
public:
    explicit Read_File(const std::string& Path)
        : File(Path == "-" ? gzdopen(fileno(stdin), "rb") : gzopen(Path.c_str(), "rb")) {
        if (File) gzbuffer(File, 1 << 20);
        for (Chunk& C : Ring) C.Data.resize(1 << 20);
    }
    ~Read_File() {
        {
            std::lock_guard<std::mutex> Lock(Mutex);
            Stop = true;
        }
        Changed.notify_all();
        if (Inflater.joinable()) Inflater.join();
        if (File) gzclose(File);
    }
    Read_File(const Read_File&) = delete;
    Read_File& operator=(const Read_File&) = delete;

    bool Is_Open() const { return File != nullptr; }
    // A read error or corrupt compressed data ended the input early.
    bool Failed() const { return Error; }

    // Starts the inflate thread; call before the first Get_Line.
    void Inflate_Ahead() {
        if (File && !Inflater.joinable()) Inflater = std::thread([this] { Inflate(); });
    }

    // The next line, without its '\n' (as std::getline); false at the end.
    bool Get_Line(std::string& Line) {
        if (Has_Pending) {
            Line.swap(Pending);
            Has_Pending = false;
            return true;
        }
        Line.clear();
        for (;;) {
            if (Pos == End && !Fill()) return !Line.empty();
            const char* Start = Data + Pos;
            const char* Newline = static_cast<const char*>(std::memchr(Start, '\n', End - Pos));
            if (Newline) {
                Line.append(Start, Newline);
                Pos += Newline - Start + 1;
                return true;
            }
            Line.append(Start, End - Pos);
            Pos = End;
        }
    }
    // Line is returned again by the next Get_Line.
    void Unget_Line(std::string& Line) {
        Pending.swap(Line);
        Has_Pending = true;
    }

private:
    struct Chunk {
        std::vector<char> Data;
        std::size_t Size = 0;  // 0: the end of the input
        bool Error = false;
    };
    static constexpr unsigned RING = 4;

    // Reads the next buffer into C; false at the end or on an error (C.Error).
    bool Read(Chunk& C) {
        const int N = gzread(File, C.Data.data(), static_cast<unsigned>(C.Data.size()));
        int Status = Z_OK;
        if (N <= 0) gzerror(File, &Status);  // Z_BUF_ERROR: the compressed data ended early
        C.Error = N < 0 || Status != Z_OK;
        C.Size = N > 0 ? static_cast<std::size_t>(N) : 0;
        return C.Size > 0;
    }

    // The inflate thread: chunk I goes to Ring[I % RING] once chunk I - RING is released.
    void Inflate() {
        for (unsigned long I = 0;; I++) {
            {
                std::unique_lock<std::mutex> Lock(Mutex);
                Changed.wait(Lock, [&] { return Stop || I < Released + RING; });
                if (Stop) return;
            }
            const bool More = Read(Ring[I % RING]);
            std::lock_guard<std::mutex> Lock(Mutex);
            Produced++;
            Changed.notify_all();
            if (!More) return;
        }
    }

    // The next buffer of the file into Data[0, End); false at the end.
    bool Fill() {
        Pos = End = 0;
        if (!File || Ended) return false;
        Chunk* C = &Ring[0];
        if (Inflater.joinable()) {
            std::unique_lock<std::mutex> Lock(Mutex);
            if (Holding) Released++;  // the chunk used up
            Holding = true;
            Changed.notify_all();
            Changed.wait(Lock, [&] { return Produced > Released; });
            C = &Ring[Released % RING];
        } else {
            Read(*C);
        }
        Data = C->Data.data();
        End = C->Size;
        if (!End) {
            Ended = true;
            Error = C->Error;
        }
        return End > 0;
    }

    gzFile File;
    Chunk Ring[RING];
    const char* Data = nullptr;
    std::size_t Pos = 0, End = 0;
    std::string Pending;
    bool Has_Pending = false, Error = false, Ended = false;
    // The inflate thread's hand-over: chunks inflated, and released by Fill.
    std::thread Inflater;
    std::mutex Mutex;
    std::condition_variable Changed;
    unsigned long Produced = 0, Released = 0;
    bool Holding = false, Stop = false;
};

// The next FASTA or FASTQ record; with qual, also its quality line (empty for FASTA).
bool ReadNextRecord(Read_File& in, std::string& id, std::string& seq, std::string* qual = nullptr) {
    std::string line;
    while (in.Get_Line(line)) {
        if (line.empty()) {
            continue;
        }

        if (line[0] == '@') {
            id = line.substr(1);
            if (!in.Get_Line(seq)) {
                return false;
            }
            if (!in.Get_Line(line)) {
                return false;
            }
            if (!in.Get_Line(qual ? *qual : line)) {
                return false;
            }
            return true;
        }

        if (line[0] == '>') {
            id = line.substr(1);
            seq.clear();
            if (qual) qual->clear();
            while (in.Get_Line(line)) {
                if (!line.empty() && (line[0] == '>' || line[0] == '@')) {
                    in.Unget_Line(line);
                    break;
                }
                seq += line;
            }
            return !seq.empty();
        }
    }
    return false;
}

void BuildIndexPaths(const std::string& genomePrefix, std::string& revBwt, std::string& revOcc,
                     std::string& fwBwt, std::string& fwOcc, std::string* revSa = nullptr, std::string* fwSa = nullptr) {
    std::size_t slash = genomePrefix.find_last_of('/');
    std::string dir;
    std::string base;

    if (slash == std::string::npos) {
        dir = "";
        base = genomePrefix;
    } else {
        dir = genomePrefix.substr(0, slash + 1);
        base = genomePrefix.substr(slash + 1);
    }

    revBwt = dir + "rev" + base + ".bwt";
    revOcc = dir + "rev" + base + ".fmv";
    fwBwt = genomePrefix + ".bwt";
    fwOcc = genomePrefix + ".fmv";
    if (revSa) *revSa = dir + "rev" + base + ".sa";
    if (fwSa) *fwSa = genomePrefix + ".sa";
}

bool FileExists(const std::string& name) {
    return std::ifstream(name.c_str()).good();
}

// Reads <genome>.pac, padded so that Genome_Mismatches can load 9 bytes anywhere.
bool LoadPac(const std::string& name) {
    std::ifstream in(name.c_str(), std::ios::binary | std::ios::ate);
    if (!in) return false;
    const std::streamsize size = in.tellg();
    in.seekg(0);
    Pac.assign(static_cast<std::size_t>(size) + 16, 0);
    return static_cast<bool>(in.read(reinterpret_cast<char*>(Pac.data()), size));
}

// The output lines of one -n 1 to 3 read (main and the batched runner), from its lists
// (Result_List): the exact hits (each capped at maxHits), then the mismatch hits in
// stage order until maxHits, each with its stage's index and mismatch count (Stages).
// Returns the read's reported hit count; with none, the caller writes "F 0 0 0".
unsigned long Format_Hits(const std::string& id, const ForwardResults* Res, const Hit_Stage* Stages,
                          unsigned long readHits, unsigned maxHits, std::string& lines) {
    for (const SARange& hit : Res[RL_EXACT].hits) {
        unsigned long hitCount = hit.Skip ? 1UL : (hit.End >= hit.Start ? (hit.End - hit.Start + 1UL) : 0UL);
        if (maxHits > 0 && hitCount > maxHits) {
            hitCount = maxHits;
        }
        lines += id + Hit_Text(hit, "F") + "\t" + std::to_string(hitCount) + "\n";
    }
    for (const SARange& hit : Res[RL_HITS].hits) {
        if (maxHits > 0 && readHits >= maxHits) break;
        unsigned long hitCount = hit.Skip ? 1UL : hit.End - hit.Start + 1UL;
        if (maxHits > 0 && readHits + hitCount > maxHits) {
            hitCount = maxHits - readHits;
        }
        readHits += hitCount;
        const Hit_Stage& Stage = Stages[hit.Stage];
        lines += id + Hit_Text(hit, Stage.Dir) + "\t" + std::to_string(hitCount) + "\t" + std::to_string(Stage.Mismatches) + "\n";
    }
    return readHits;
}

// SAM/BAM output (--sam / --bam): the genome's sequences and N runs, the writer, and
// whether the reverse complement of each read is searched too.
struct Sam_Output {
    Genome_Map Map;
    Sam_Writer Writer;
    bool Reverse = true;
    bool Write = true;  // false with --quiet
    std::vector<uint64_t> Positions;
};

// ---- SA rows to SA values for SAM/BAM output (Sa_Decode) ----
//
// BWTSaValue walks a row back with LF (BWTPsiMinusValue) until it reaches a sampled row
// (every saInterval-th row), about 7 dependent cache misses, far more where LF strides
// regularly (a homopolymer run, e.g. an N run stored as C: a block of rows shifts by a
// multiple of saInterval and the same rows miss the samples step after step).
// Sa_Decode does many rows at once, each round stepping every unit of work, all of
// whose rows are Skipped steps from their requests:
//   - a span: rows Lo..Hi, of which the ones whose (row - Lo) % saInterval is a class
//     in Mask are resolved. Its sampled rows form one class, read off from saValue
//     together. If all its rows have one BWT character (Occ_All at both ends),
//     LF shifts it whole: the same rows, classes and outputs, at O(1) a step however
//     long it stalls. Otherwise it splits into its runs of one character: runs of at
//     least SA_SPAN_MIN rows stay spans, the others' unresolved rows become groups.
//   - a group: unresolved rows in row order. Its sampled rows (and the $ row) are read
//     off; the rest map by LF to at most 4 child groups (the rows with each BWT
//     character; LF keeps their order). Occ comes from one Occ_All at the first
//     row, then by counting the BWT characters between neighbouring rows, or a fresh
//     Occ_All across a gap over SA_GAP_LOOKUP rows.
//   - a walk: a single row (SA_SINGLE_WALK), stepped like BWTSaValue; the walks of a
//     round go side by side with prefetching, so their cache misses overlap.
// SA_DECODE 0 walks every row; SA_SPAN_MIN 0 makes no spans.
#ifndef SA_DECODE
#define SA_DECODE 1
#endif
#ifndef SA_GAP_LOOKUP
#define SA_GAP_LOOKUP 64
#endif
#ifndef SA_SINGLE_WALK
#define SA_SINGLE_WALK 1
#endif
#ifndef SA_SPAN_MIN
#define SA_SPAN_MIN 16
#endif

// Rows Start..End of Fmi, whose SA values go to Values[Out ...].
struct Sa_Request {
    const ActiveBWT* Fmi;
    unsigned long Start, End;
    std::size_t Out;
};

#ifdef BATMAN_BATCH_STATS
struct Decode_Stats_Counts {  // per thread, as Batch_Stats
    unsigned long Rows = 0, Direct = 0, Walk_Steps = 0, Range_Steps = 0, Span_Steps = 0;
    double Seconds = 0;  // SAM output up to writing: candidates, chains, Sa_Decode, the records built
    void Add(const Decode_Stats_Counts& O) {
        Rows += O.Rows;
        Direct += O.Direct;
        Walk_Steps += O.Walk_Steps;
        Range_Steps += O.Range_Steps;
        Span_Steps += O.Span_Steps;
    }
};
thread_local Decode_Stats_Counts Decode_Stats;
#endif

void Sa_Decode(const std::vector<Sa_Request>& Requests, std::vector<unsigned>& Values) {
    struct Walk {
        const ActiveBWT* Fmi;
        unsigned long Row;
        unsigned Skipped;
        std::size_t Out;
    };
    struct Item {
        unsigned long Row;
        std::size_t Out;
    };
    struct Group {  // Items[Begin .. Begin + Count)
        const ActiveBWT* Fmi;
        unsigned Skipped;
        std::size_t Begin, Count;
    };
    struct Span {  // rows Lo..Hi, row r's output Out + (r - Lo)
        const ActiveBWT* Fmi;
        unsigned Skipped;
        unsigned long Lo, Hi;
        std::size_t Out;
        uint64_t Mask;  // resolved classes of (r - Lo) % saInterval
    };
    static thread_local std::vector<Walk> Walks, Next_Walks;
    static thread_local std::vector<Item> Items, Next_Items, Kids[4];
    static thread_local std::vector<Group> Groups, Next_Groups;
    static thread_local std::vector<Span> Spans, Next_Spans;
    Walks.clear();
    Items.clear();
    Groups.clear();
    Spans.clear();
#ifdef BATMAN_BATCH_STATS
    static FILE* Dump = std::getenv("BATMAN_DECODE_DUMP") ? std::fopen(std::getenv("BATMAN_DECODE_DUMP"), "w") : nullptr;
#endif
    for (const Sa_Request& R : Requests) {
        const unsigned long Count = R.End - R.Start + 1;
#ifdef BATMAN_BATCH_STATS
        Decode_Stats.Rows += Count;
        if (Dump) std::fprintf(Dump, "%d %lu %lu\n", R.Fmi == fwfmi, R.Start, R.End);
#endif
        if (SA_DECODE == 0 || (SA_SINGLE_WALK && Count == 1)) {
            for (unsigned long Row = R.Start; Row <= R.End; Row++) Walks.push_back({R.Fmi, Row, 0, R.Out + (Row - R.Start)});
        } else if (SA_SPAN_MIN && Count >= SA_SPAN_MIN && R.Fmi->saInterval <= 64) {
            Spans.push_back({R.Fmi, 0, R.Start, R.End, R.Out, 0});
        } else {
            Groups.push_back({R.Fmi, 0, Items.size(), Count});
            for (unsigned long Row = R.Start; Row <= R.End; Row++) Items.push_back({Row, R.Out + (Row - R.Start)});
        }
    }
    auto Char_At = [](const ActiveBWT* F, unsigned long Row) -> unsigned {
        const unsigned long Index = Row - (Row > F->inverseSa0);
        return F->bwtCode[Index / 16] << ((Index % 16) * 2) >> (BITS_IN_WORD - 2);
    };
    // Kids[c] as the next round's groups.
    auto Emit_Kids = [&](const ActiveBWT* F, unsigned Skipped) {
        for (const std::vector<Item>& K : Kids) {
            if (K.empty()) continue;
            Next_Groups.push_back({F, Skipped, Next_Items.size(), K.size()});
            Next_Items.insert(Next_Items.end(), K.begin(), K.end());
        }
    };
    // One step of a group: read off its sampled rows, then its children by BWT character.
    auto Group_Step = [&](const Group& G) {
        const ActiveBWT* F = G.Fmi;
        const unsigned Interval = F->saInterval;
        Item* It = Items.data() + G.Begin;
        std::size_t Left = 0;
        for (std::size_t i = 0; i < G.Count; i++) {
            const Item X = It[i];
            if (X.Row % Interval == 0) {
                Values[X.Out] = F->saValue[X.Row / Interval] + G.Skipped;
#ifdef BATMAN_BATCH_STATS
                Decode_Stats.Direct++;
#endif
            } else if (X.Row == F->inverseSa0) {  // the $ row: LF goes to row 0 (saValue[0] is -1)
                Values[X.Out] = F->saValue[0] + G.Skipped + 1;
            } else {
                It[Left++] = X;
            }
        }
        if (!Left) return;
        if (SA_SINGLE_WALK && Left == 1) {
#ifdef BATMAN_BATCH_STATS
            Decode_Stats.Walk_Steps++;
#endif
            Next_Walks.push_back({F, BWTPsiMinusValue(F, static_cast<unsigned>(It[0].Row)), G.Skipped + 1, It[0].Out});
            return;
        }
        for (std::vector<Item>& K : Kids) K.clear();
        unsigned ALIGN_16 Occ[ALPHABET_SIZE];
        unsigned long Prev = 0;
        for (std::size_t i = 0; i < Left; i++) {
            const unsigned long Row = It[i].Row;
            if (i == 0 || Row - Prev > SA_GAP_LOOKUP) {
                fm_occ::Occ_All(F, Row, Occ);
#ifdef BATMAN_BATCH_STATS
                Decode_Stats.Range_Steps++;
#endif
            } else {
                for (unsigned long r = Prev + 1; r < Row; r++)
                    if (r != F->inverseSa0) Occ[Char_At(F, r)]++;
            }
            const unsigned c = Char_At(F, Row);
            Kids[c].push_back({F->cumulativeFreq[c] + Occ[c] + 1, It[i].Out});
            Occ[c]++;
            Prev = Row;
        }
        Emit_Kids(F, G.Skipped + 1);
    };
    // One step of a span: read off its sampled class, then shift it whole or split it.
    auto Span_Step = [&](const Span& P) {
        const ActiveBWT* F = P.Fmi;
        const unsigned Interval = F->saInterval;
        const uint64_t Full = Interval == 64 ? ~uint64_t{0} : (uint64_t{1} << Interval) - 1;
        auto Unresolved = [&](uint64_t Mask, unsigned long r) { return !(Mask >> ((r - P.Lo) % Interval) & 1); };
        for (std::vector<Item>& K : Kids) K.clear();
        if (F->inverseSa0 >= P.Lo && F->inverseSa0 <= P.Hi) {  // the $ row has no character: as a group, same step
            for (unsigned long r = P.Lo; r <= P.Hi; r++)
                if (Unresolved(P.Mask, r)) Kids[0].push_back({r, P.Out + (r - P.Lo)});
            Emit_Kids(F, P.Skipped);
            return;
        }
        uint64_t Mask = P.Mask;
        const unsigned k = (Interval - P.Lo % Interval) % Interval;  // the class of the sampled rows
        if (!(Mask >> k & 1)) {
            for (unsigned long r = P.Lo + k; r <= P.Hi; r += Interval) {
                Values[P.Out + (r - P.Lo)] = F->saValue[r / Interval] + P.Skipped;
#ifdef BATMAN_BATCH_STATS
                Decode_Stats.Direct++;
#endif
            }
            Mask |= uint64_t{1} << k;
        }
        if (Mask == Full) return;
#ifdef BATMAN_BATCH_STATS
        Decode_Stats.Span_Steps++;
#endif
        const unsigned long Width = P.Hi - P.Lo + 1;
        unsigned ALIGN_16 Occ[ALPHABET_SIZE], Occ_Hi[ALPHABET_SIZE];
        fm_occ::Occ_All(F, P.Lo, Occ);
        if (P.Hi < F->textLength) {
            fm_occ::Occ_All(F, P.Hi + 1, Occ_Hi);
            for (unsigned c = 0; c < 4; c++)
                if (Occ_Hi[c] - Occ[c] == Width) {
                    const unsigned long Lo = F->cumulativeFreq[c] + Occ[c] + 1;
                    Next_Spans.push_back({F, P.Skipped + 1, Lo, Lo + Width - 1, P.Out, Mask});
                    return;
                }
        }
        for (unsigned long s = P.Lo; s <= P.Hi;) {
            const unsigned c = Char_At(F, s);
            unsigned long e = s;
            while (e < P.Hi && Char_At(F, e + 1) == c) e++;
            const unsigned long Lo = F->cumulativeFreq[c] + Occ[c] + 1;
            if (e - s + 1 >= SA_SPAN_MIN) {
                const unsigned d = (s - P.Lo) % Interval;  // class j of the run is class j + d of the span
                const uint64_t Run_Mask = d ? ((Mask >> d) | (Mask << (Interval - d))) & Full : Mask;
                Next_Spans.push_back({F, P.Skipped + 1, Lo, Lo + (e - s), P.Out + (s - P.Lo), Run_Mask});
            } else {
                for (unsigned long r = s; r <= e; r++)
                    if (Unresolved(Mask, r)) Kids[c].push_back({Lo + (r - s), P.Out + (r - P.Lo)});
            }
            Occ[c] += e - s + 1;
            s = e + 1;
        }
        Emit_Kids(F, P.Skipped + 1);
    };
    while (!Walks.empty() || !Groups.empty() || !Spans.empty()) {
        for (const Walk& W : Walks) {
            if (W.Row % W.Fmi->saInterval == 0) __builtin_prefetch(W.Fmi->saValue + W.Row / W.Fmi->saInterval);
            else Prefetch_Row(W.Fmi, W.Row + 1);
        }
        for (const Group& G : Groups) Prefetch_Row(G.Fmi, Items[G.Begin].Row);
        for (const Span& P : Spans) {
            Prefetch_Row(P.Fmi, P.Lo);
            Prefetch_Row(P.Fmi, P.Hi + 1);
        }
        Next_Walks.clear();
        Next_Groups.clear();
        Next_Items.clear();
        Next_Spans.clear();
        for (const Walk& W : Walks) {
            if (W.Row % W.Fmi->saInterval == 0) {
                Values[W.Out] = W.Fmi->saValue[W.Row / W.Fmi->saInterval] + W.Skipped;
                continue;
            }
#ifdef BATMAN_BATCH_STATS
            Decode_Stats.Walk_Steps++;
#endif
            Next_Walks.push_back({W.Fmi, BWTPsiMinusValue(W.Fmi, static_cast<unsigned>(W.Row)), W.Skipped + 1, W.Out});
        }
        for (const Group& G : Groups) Group_Step(G);
        for (const Span& P : Spans) Span_Step(P);
        std::swap(Walks, Next_Walks);
        std::swap(Groups, Next_Groups);
        std::swap(Items, Next_Items);
        std::swap(Spans, Next_Spans);
    }
}

// A read's hits for SAM/BAM output, in output order: by mismatches, the read's strand
// first on a tie (stable).
struct Sam_Candidate {
    int Mismatches, Strand;
    const SARange* Hit;
    bool Forward_Index;
};

void Sam_Candidates(const ForwardResults* const Res[2], const Hit_Stage* const Stages[2], std::vector<Sam_Candidate>& Out) {
    Out.clear();
    for (int Strand = 0; Strand < 2; Strand++) {
        if (!Res[Strand]) continue;
        for (const SARange& h : Res[Strand][RL_EXACT].hits) Out.push_back({0, Strand, &h, false});
        for (const SARange& h : Res[Strand][RL_HITS].hits) {
            const Hit_Stage& Stage = Stages[Strand][h.Stage];
            Out.push_back({Stage.Mismatches, Strand, &h, Stage.Dir[0] == 'B'});
        }
    }
    std::stable_sort(Out.begin(), Out.end(),
                     [](const Sam_Candidate& a, const Sam_Candidate& b) { return a.Mismatches < b.Mismatches; });
}

// Whether a hit's position is known without decoding: located, or recorded at a sampled
// row (Track_Pos).
inline bool Hit_Has_Pos(const SARange& hit) { return hit.Skip && (Use_Ref || Track_Pos) && hit.Pos != NO_POS; }

// The rows Encode_Sam_Read will decode for a read, as decode requests: its candidates'
// rows in output order, up to MaxHits locations (0: all) counting the ones with a known
// Pos. Returns the number of values requested (from Values_At on). Origins gets each
// request's candidate.
std::size_t Sam_Requests(const std::vector<Sam_Candidate>& Candidates, unsigned MaxHits, std::size_t Values_At,
                         std::vector<Sa_Request>& Requests, std::vector<const Sam_Candidate*>& Origins) {
    std::size_t Locations = 0, Values = 0;
    for (const Sam_Candidate& C : Candidates) {
        if (MaxHits && Locations >= MaxHits) break;
        const SARange& h = *C.Hit;
        if (Hit_Has_Pos(h)) {
            Locations++;
            continue;
        }
        unsigned long Last = h.Skip ? h.Start : h.End;
        if (MaxHits) Last = std::min<unsigned long>(Last, h.Start + (MaxHits - Locations) - 1);
        Requests.push_back({C.Forward_Index ? fwfmi : revfmi, h.Start, Last, Values_At + Values});
        Origins.push_back(&C);
        Locations += Last - h.Start + 1;
        Values += Last - h.Start + 1;
    }
    return Values;
}

// Builds one read's SAM/BAM records (into E) from its strands' lists (Res[1]: the reverse
// complement's, or null): the hits in order of mismatches (forward strand first on a
// tie), each location placed on its sequence (Genome_Map::Place; locations leaving a
// sequence or covering an N are dropped, but count toward MaxHits as if written), at
// most MaxHits locations (0: all). The first record is primary. No record: one
// unmapped record. Returns the mapped records. The SA values of the rows Sam_Requests
// listed, which are all the rows that count, are Decoded[0 .. Count).
unsigned long Encode_Sam_Read(const Sam_Output& S, Sam_Encoder& E, const std::string& Id, const std::string& Seq_Text, const std::string& Qual,
                             std::string_view Codes, const std::vector<Sam_Candidate>& Candidates, unsigned MaxHits,
                             const unsigned* Decoded, std::size_t Count) {
    struct Record {
        int Tid, Strand;
        int64_t Pos;
        const SARange* Hit;
    };
    static thread_local std::vector<Record> Records;
    Records.clear();
    const int L = static_cast<int>(Codes.size());
    std::size_t Next = 0, Locations = 0;
    auto Add = [&](const Sam_Candidate& C, uint64_t G) {
        int Tid;
        int64_t Pos;
        if (S.Map.Place(G, L, Tid, Pos)) Records.push_back({Tid, C.Strand, Pos, C.Hit});
        return MaxHits && ++Locations >= MaxHits;
    };
    for (const Sam_Candidate& C : Candidates) {
        const SARange& h = *C.Hit;
        if (Hit_Has_Pos(h)) {
            if (Add(C, h.Pos)) break;
            continue;
        }
        bool Full = false;
        const unsigned long Last = h.Skip ? h.Start : h.End;
        for (unsigned long Row = h.Start; Row <= Last && !Full && Next < Count; Row++) {
            const unsigned Sa = Decoded[Next++];
            Full = Add(C, C.Forward_Index ? Sa : Text_Length - L - Sa);
        }
        if (Full) break;
    }
    if (!S.Write) return Records.size();
    const std::string Name = Id.substr(0, Id.find_first_of(" \t"));
    if (Records.empty()) {
        E.Unmapped(Name, Seq_Text, Qual);
        return 0;
    }
    std::string Seq[2], Qual_Of[2];
    Seq[0].resize(L);
    Seq[1].resize(L);
    for (int i = 0; i < L; i++) {
        Seq[0][i] = "ACGT"[static_cast<unsigned char>(Codes[i])];
        Seq[1][L - 1 - i] = "TGCA"[static_cast<unsigned char>(Codes[i])];
    }
    Qual_Of[0] = Qual.size() == static_cast<std::size_t>(L) ? Qual : std::string();
    Qual_Of[1] = std::string(Qual_Of[0].rbegin(), Qual_Of[0].rend());
    std::vector<std::pair<int, int>> Mismatches;
    for (std::size_t k = 0; k < Records.size(); k++) {
        const Record& R = Records[k];
        Mismatches.clear();
        for (int i = 0; i < R.Hit->Mismatches && i < MAX_MISMATCHES_BOUND; i++)
            Mismatches.push_back({R.Hit->Mismatch_Pos[i], static_cast<int>(R.Hit->Mismatch_Char >> (2 * i) & 3)});
        const uint16_t Flag = (R.Strand ? BAM_FREVERSE : 0) | (k ? BAM_FSECONDARY : 0);
        E.Mapped(Name, Flag, R.Tid, R.Pos, Seq[R.Strand], Qual_Of[R.Strand], Mismatches, static_cast<int>(Records.size()));
    }
    return Records.size();
}

// The threads of Run_Batched (-t): Run(Job) calls Job(t) on each of the Size() threads
// (the caller is thread 0) and returns when all have finished. At_Exit runs on each
// worker as it ends (the stats builds add its counts to the main thread's).
class Thread_Pool {
public:
    Thread_Pool(unsigned Count, std::function<void()> At_Exit) {
        for (unsigned t = 1; t < Count; t++)
            Threads.emplace_back([this, t, At_Exit] {
                Work(t);
                At_Exit();
            });
    }
    ~Thread_Pool() {
        {
            std::lock_guard<std::mutex> Lock(Mutex);
            Stop = true;
        }
        Wake.notify_all();
        for (std::thread& T : Threads) T.join();
    }
    unsigned Size() const { return static_cast<unsigned>(Threads.size()) + 1; }

    void Run(const std::function<void(unsigned)>& Job) {
        if (Threads.empty()) {
            Job(0);
            return;
        }
        {
            std::lock_guard<std::mutex> Lock(Mutex);
            Current = &Job;
            Pending = static_cast<unsigned>(Threads.size());
            Generation++;
        }
        Wake.notify_all();
        Job(0);
        std::unique_lock<std::mutex> Lock(Mutex);
        Done.wait(Lock, [&] { return Pending == 0; });
    }

private:
    void Work(unsigned t) {
        unsigned long Seen = 0;
        for (;;) {
            const std::function<void(unsigned)>* Job;
            {
                std::unique_lock<std::mutex> Lock(Mutex);
                Wake.wait(Lock, [&] { return Stop || Generation != Seen; });
                if (Stop) return;
                Seen = Generation;
                Job = Current;
            }
            (*Job)(t);
            std::lock_guard<std::mutex> Lock(Mutex);
            if (--Pending == 0) Done.notify_one();
        }
    }

    std::vector<std::thread> Threads;
    std::mutex Mutex;
    std::condition_variable Wake, Done;
    const std::function<void(unsigned)>* Current = nullptr;
    unsigned Pending = 0;
    unsigned long Generation = 0;
    bool Stop = false;
};

// A block of records as read, not yet checked or encoded (Run_Batched).
struct Input_Block {
    std::vector<std::string> Ids, Seqs, Quals;
    unsigned N = 0;
    bool More = true;     // false: the input ended in this block
    double Seconds = 0;   // reading it
};

// Reads the blocks in turn. With Ahead, a thread of its own reads the next block into
// a second buffer while the caller works on the one Next returned (double buffering);
// otherwise Next reads it. Next returns null once the input has ended; the block it
// returned before stays valid until it is called again.
class Block_Reader {
public:
    Block_Reader(Read_File& In, unsigned Size, bool Qual, bool Ahead) : In(In), Qual(Qual) {
        if (Ahead) In.Inflate_Ahead();
        for (Input_Block& B : Blocks) {
            B.Ids.resize(Size);
            B.Seqs.resize(Size);
            B.Quals.resize(Qual ? Size : 0);
        }
        if (Ahead) Thread = std::thread([this] { Read_Ahead(); });
    }
    ~Block_Reader() {
        {
            std::lock_guard<std::mutex> Lock(Mutex);
            Stop = true;
        }
        Changed.notify_all();
        if (Thread.joinable()) Thread.join();
    }

    const Input_Block* Next() {
        if (!Thread.joinable()) {
            if (Ended) return nullptr;
            Read(Blocks[0]);
            Ended = !Blocks[0].More;
            return &Blocks[0];
        }
        std::unique_lock<std::mutex> Lock(Mutex);
        Changed.wait(Lock, [&] { return Filled > Taken || Ended; });
        if (Filled == Taken) return nullptr;
        const Input_Block* B = &Blocks[Taken % 2];
        Taken++;  // the block before it is free
        Changed.notify_all();
        return B;
    }

private:
    void Read(Input_Block& B) {
        const auto Start = std::chrono::steady_clock::now();
        B.N = 0;
        while (B.N < B.Ids.size() && (B.More = ReadNextRecord(In, B.Ids[B.N], B.Seqs[B.N], Qual ? &B.Quals[B.N] : nullptr)))
            B.N++;
        B.Seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - Start).count();
    }

    void Read_Ahead() {
        for (unsigned long I = 0;; I++) {
            {  // block I's buffer is free once block I - 1 has been taken
                std::unique_lock<std::mutex> Lock(Mutex);
                Changed.wait(Lock, [&] { return Stop || I < 2 || Taken >= I; });
                if (Stop) return;
            }
            Read(Blocks[I % 2]);
            std::lock_guard<std::mutex> Lock(Mutex);
            Filled++;
            Ended = !Blocks[I % 2].More;
            Changed.notify_all();
            if (Ended) return;
        }
    }

    Read_File& In;
    const bool Qual;
    Input_Block Blocks[2];
    std::thread Thread;
    std::mutex Mutex;
    std::condition_variable Changed;
    unsigned long Filled = 0, Taken = 0;  // blocks read, and returned by Next
    bool Ended = false, Stop = false;
};

// Reads the records in blocks, searches each
// block with Batch_Exact (-n 0, or reads too short to split into halves) or the -n k
// program (Batch_One_Mismatch ... Batch_Five_Mismatch, Batch_Both_Strands) and writes
// its lines in input order. Encoding the
// reads counts as input here. Returns false on an unsupported read length.
// With sam, writes SAM/BAM records instead (Encode_Sam_Read, Sam_Writer::Write). The reverse complements of
// the block's reads, if searched, follow the reads: one lane searches both strands of a
// read in zig-zag order (Batch_Both_Strands); the exact search takes them as more reads.
// With Threads > 1, a block holds 2 * Threads chunks of reads, which the threads take in
// turn to search (each with its own lanes and Batch_Chooser); then each thread formats
// (or, for SAM/BAM, decodes) a contiguous part of the block, and the main thread writes
// the parts in order. Before the search, each thread checks and encodes its part of the
// block's reads, and a reader thread (Block_Reader) reads the next block meanwhile. The
// output does not depend on Threads. inputSeconds: waiting for the reader, and encoding;
// readSeconds: the reading itself (overlapped with the rest when Threads > 1, and part of
// inputSeconds otherwise).
std::string Batch_Report;  // the lanes used, for --timing (Run_Batched)

bool Run_Batched(Read_File& in, std::ostream& out, unsigned B, unsigned Threads, int maxMismatches,
                 unsigned maxHits, bool quiet, uint64_t& records, uint64_t& reportedHits, double& tableBuildSeconds,
                 double& inputSeconds, double& readSeconds, double& searchSeconds, double& outputSeconds,
                 Sam_Output* sam = nullptr) {
    using Clock = std::chrono::steady_clock;
#ifdef BATMAN_BATCH_BLOCK
    constexpr unsigned CHUNK = BATMAN_BATCH_BLOCK;
#else
    constexpr unsigned CHUNK = 4096;
#endif
    const unsigned BLOCK = Threads > 1 ? 2 * Threads * CHUNK : CHUNK;
    const unsigned strands = sam && sam->Reverse ? 2 : 1;  // searched copies of each read
#ifdef BATMAN_BATCH_STATS
    Batch_Stats_Counts* const mainBatch = &Batch_Stats;
    Decode_Stats_Counts* const mainDecode = &Decode_Stats;
    std::mutex statsMutex;
    auto atExit = [&] {
        std::lock_guard<std::mutex> Lock(statsMutex);
        mainBatch->Add(Batch_Stats);
        mainDecode->Add(Decode_Stats);
    };
#else
    auto atExit = [] {};
#endif
    Thread_Pool pool(Threads, atExit);
    std::vector<Batch_Chooser> choosers(Threads);  // B = 0: --batch auto
    Block_Reader reader(in, BLOCK, sam != nullptr, Threads > 1);
    std::vector<unsigned char> searched(BLOCK);  // valid read of the right length
    std::vector<unsigned> slot(BLOCK);            // its place among the searched reads
    std::vector<SARange> ranges(strands * BLOCK);
    std::vector<unsigned char> found(strands * BLOCK);
    std::vector<ForwardResults> results(MM_LISTS * strands * BLOCK);  // -n 1 to 5: each read's result lists
    std::vector<unsigned long> readHits(strands * BLOCK);
    // SAM/BAM: each read's hits in output order, and the SA rows decoded for them (in its
    // thread's part).
    std::vector<std::vector<Sam_Candidate>> candidates(sam ? BLOCK : 0);
    std::vector<std::size_t> decodeBegin(sam ? BLOCK : 0), decodeCount(sam ? BLOCK : 0);
    struct Part {  // a thread's part of the block's output
        std::vector<Sa_Request> decodeRequests;
        std::vector<const Sam_Candidate*> origins;  // the current read's requests' candidates
        std::vector<unsigned> decoded;
        std::string lines;
        Sam_Encoder records;  // SAM/BAM: the part's records
        uint64_t hits = 0;
    };
    std::vector<Part> parts(Threads);
    std::vector<unsigned> partReads(Threads);  // searched reads in each part
    std::string reads;  // codes of the searched reads, back to back
    bool initialized = false;
    for (;;) {
        const auto waitStart = Clock::now();
        const Input_Block* block = reader.Next();
        if (!block) break;
        const unsigned n = block->N;
        const std::vector<std::string>& ids = block->Ids;
        const std::vector<std::string>& seqs = block->Seqs;  // the reads as given
        const std::vector<std::string>& quals = block->Quals;
        records += n;
        readSeconds += block->Seconds;
        double tableSeconds = 0;
        if (!initialized && n) {
            STRINGLENGTH = static_cast<int>(seqs[0].size());
            if (STRINGLENGTH <= 0 || STRINGLENGTH >= MAXTAG) {
                std::fprintf(stderr, "Unsupported read length: %d\n", STRINGLENGTH);
                return false;
            }
            const auto tableStart = Clock::now();
            Build_Tables(revfmi);
            tableSeconds = std::chrono::duration<double>(Clock::now() - tableStart).count();
            tableBuildSeconds += tableSeconds;
            initialized = true;
        }
        // Thread t's part of the block: records [Part_Begin(t), Part_Begin(t + 1)).
        auto Part_Begin = [&](unsigned t) { return static_cast<unsigned>(uint64_t(n) * t / Threads); };
        // A read is searched if it has the first read's length and only ACGT (any case);
        // its codes go to reads at its place among the searched reads (slot), and with
        // strands == 2 its reverse complement's at count + slot.
        const std::size_t L = static_cast<std::size_t>(STRINGLENGTH);
        pool.Run([&](unsigned t) {
            unsigned c = 0;
            for (unsigned j = Part_Begin(t); j < Part_Begin(t + 1); j++) {
                const std::string& seq = seqs[j];
                bool valid = seq.size() == L;
                for (std::size_t i = 0; i < L && valid; i++) valid = Base_Codes[static_cast<unsigned char>(seq[i])] >= 0;
                searched[j] = valid;
                c += valid;
            }
            partReads[t] = c;
        });
        unsigned count = 0;
        for (unsigned t = 0; t < Threads; t++) {
            const unsigned c = partReads[t];
            partReads[t] = count;  // the part's first slot
            count += c;
        }
        reads.resize(strands * count * L);
        pool.Run([&](unsigned t) {
            unsigned r = partReads[t];
            for (unsigned j = Part_Begin(t); j < Part_Begin(t + 1); j++) {
                slot[j] = r;
                if (!searched[j]) continue;
                char* code = &reads[r * L];
                const char* seq = seqs[j].data();
                for (std::size_t i = 0; i < L; i++) code[i] = Base_Codes[static_cast<unsigned char>(seq[i])];
                if (strands == 2) {
                    char* rc = &reads[(count + r) * L];
                    for (std::size_t i = 0; i < L; i++) rc[i] = static_cast<char>(3 - code[L - 1 - i]);
                }
                r++;
            }
        });
        inputSeconds += std::chrono::duration<double>(Clock::now() - waitStart).count() - tableSeconds;
        if (!initialized) continue;  // no reads yet

        const bool oneMismatch = maxMismatches >= 1 && STRINGLENGTH / 2 > LOOKUPSIZE;
        const auto searchStart = Clock::now();
        // The block's search: the -n k program (both strands of a read in one lane, zig-zag)
        // or the exact search (each strand a read), in chunks (one with a single thread).
        const bool zigZag = oneMismatch && strands == 2;
        const unsigned units = zigZag ? count : strands * count;
        const unsigned chunk = Threads > 1 ? CHUNK : std::max(units, 1U);
        auto Search_Range = [&](unsigned First, unsigned Last, unsigned Lanes) {
            if (oneMismatch)
                Batch_Search(reads.data(), units, STRINGLENGTH, Lanes, maxHits,
                             zigZag ? Zig_Zag_Program(maxMismatches) : Program_Of(maxMismatches), &results, &readHits,
                             nullptr, nullptr, zigZag, First, Last);
            else
                Batch_Search(reads.data(), units, STRINGLENGTH, Lanes, 0, Program_Of(0), nullptr, nullptr, &ranges, &found,
                             false, First, Last);
        };
        std::atomic<unsigned> nextChunk{0};
        pool.Run([&](unsigned t) {
            for (unsigned first; (first = nextChunk.fetch_add(chunk)) < units;) {
                const unsigned last = std::min(units, first + chunk);
                if (B) Search_Range(first, last, B);
                else choosers[t].Run(first, last, Search_Range);
            }
        });
        searchSeconds += std::chrono::duration<double>(Clock::now() - searchStart).count();

        const auto outputStart = Clock::now();
        if (sam) {
            // Each read's hits in output order, the rows to decode for them (Sam_Requests),
            // all of the part's decoded together (Sa_Decode), then its records, which the
            // main thread writes part by part.
            pool.Run([&](unsigned t) {
                Part& p = parts[t];
                p.decodeRequests.clear();
                std::size_t values = 0;
                for (unsigned j = Part_Begin(t); j < Part_Begin(t + 1); j++) {
                    candidates[j].clear();
                    if (!searched[j]) continue;
                    if (oneMismatch) {
                        const ForwardResults* res[2] = {nullptr, nullptr};
                        const Hit_Stage* stages[2] = {Stages_Of(maxMismatches), Stages_Of(maxMismatches)};
                        for (unsigned s = 0; s < strands; s++) res[s] = &results[MM_LISTS * (s * count + slot[j])];
                        Sam_Candidates(res, stages, candidates[j]);
                    } else {
                        for (unsigned s = 0; s < strands; s++) {
                            const std::size_t at = s * count + slot[j];
                            if (found[at] && ranges[at].Start) candidates[j].push_back({0, static_cast<int>(s), &ranges[at], false});
                        }
                    }
                    decodeBegin[j] = values;
                    p.origins.clear();
                    decodeCount[j] = Sam_Requests(candidates[j], maxHits, values, p.decodeRequests, p.origins);
                    values += decodeCount[j];
                    p.decoded.resize(values);
                }
                p.decoded.resize(values);
                Sa_Decode(p.decodeRequests, p.decoded);
                p.hits = 0;
                if (sam->Write) sam->Writer.Start(p.records);
                for (unsigned j = Part_Begin(t); j < Part_Begin(t + 1); j++) {
                    if (!searched[j]) {
                        if (sam->Write) p.records.Unmapped(ids[j].substr(0, ids[j].find_first_of(" \t")), seqs[j], quals[j]);
                        continue;
                    }
                    p.hits += Encode_Sam_Read(*sam, p.records, ids[j], seqs[j], quals[j],
                                              std::string_view(reads.data() + static_cast<std::size_t>(slot[j]) * STRINGLENGTH, STRINGLENGTH),
                                              candidates[j], maxHits, p.decoded.data() + decodeBegin[j], decodeCount[j]);
                }
            });
#ifdef BATMAN_BATCH_STATS
            Decode_Stats.Seconds += std::chrono::duration<double>(Clock::now() - outputStart).count();
#endif
            for (const Part& p : parts) {
                reportedHits += p.hits;
                if (sam->Write && !sam->Writer.Write(p.records)) {
                    std::fprintf(stderr, "Error writing the SAM/BAM output\n");
                    return false;
                }
            }
            outputSeconds += std::chrono::duration<double>(Clock::now() - outputStart).count();
            continue;
        }

        pool.Run([&](unsigned t) {
            std::string& lines = parts[t].lines;
            uint64_t& partHits = parts[t].hits;
            lines.clear();
            partHits = 0;
            for (unsigned j = Part_Begin(t); j < Part_Begin(t + 1) && oneMismatch; j++) {
                unsigned long hits = 0;
                if (searched[j]) {
                    const ForwardResults* r = &results[MM_LISTS * static_cast<std::size_t>(slot[j])];
                    const std::size_t at = lines.size();
                    hits = Format_Hits(ids[j], r, Stages_Of(maxMismatches), readHits[slot[j]], maxHits, lines);
                    if (quiet) lines.resize(at);
                }
                partHits += hits;
                if (!hits && !quiet) lines += ids[j] + "\tF\t0\t0\t0\n";
            }
            for (unsigned j = Part_Begin(t); j < Part_Begin(t + 1) && !oneMismatch; j++) {
                const SARange& range = ranges[slot[j]];
                if (!searched[j] || !found[slot[j]] || range.Start == 0) {
                    if (!quiet) lines += ids[j] + "\tF\t0\t0\t0\n";
                    continue;
                }
                unsigned long hitCount = range.Skip ? 1UL : (range.End >= range.Start ? (range.End - range.Start + 1UL) : 0UL);
                if (maxHits > 0 && hitCount > maxHits) hitCount = maxHits;
                partHits += hitCount;
                if (!quiet) lines += ids[j] + Hit_Text(range, "F") + "\t" + std::to_string(hitCount) + "\n";
            }
        });
        for (const Part& p : parts) {
            reportedHits += p.hits;
            if (!quiet) out << p.lines;
        }
        outputSeconds += std::chrono::duration<double>(Clock::now() - outputStart).count();
    }
    if (B) {
        Batch_Report = "batch=" + std::to_string(B);
    } else if (Threads == 1) {
        Batch_Report = choosers[0].Report;
    } else {  // each thread's choice
        Batch_Report = "batch=auto:[";
        for (unsigned t = 0; t < Threads; t++)
            Batch_Report += (t ? "," : "") + (choosers[t].Chosen ? std::to_string(choosers[t].Chosen) : std::string("-"));
        Batch_Report += "]";
    }
    if (!B) std::fprintf(stderr, "%s\n", Auto_Profile(choosers).c_str());
    return true;
}

void Usage(const char* prog) {
    std::fprintf(stderr,
    "Usage: %s -q <reads.fa/fq[.gz]> -g <genome_prefix> [-n <0-7>] [-o <hits.txt>] [-m <maxhits>] [--no-ref] [--lookup <k>] [--batch <B|auto>] [-t <threads>] [--sam|--bam] [--forward-only] [--quiet|-Q] [--timing]\n"
        "  -q  the reads: FASTA or FASTQ, plain or gzip-compressed (- for stdin)\n"
        "  -n  maximum mismatches, 0-7 (default 0); 1 to 7 also need the forward index\n"
        "      <genome_prefix>.bwt/.fmv\n"
        "  -m  maximum hits per read (default 1, 0 = all)\n"
        "  --no-ref  do not load the SA samples (.sa) and packed genome (.pac); single-row\n"
        "            paths then walk the index to the end (less memory, slower)\n"
        "  --sam, --bam  write SAM / BAM (htslib; also chosen by an -o ending .sam / .bam;\n"
        "            stdout without -o): both strands, NM/MD/NH tags, hits on N runs or across\n"
        "            sequences dropped. Needs the SA samples (.sa) even with --no-ref, which\n"
        "            then records positions at sampled rows of single-row paths\n"
        "  --forward-only  SAM/BAM: search the reads as given, not their reverse complements\n"
        "  --switch-walk  index switches re-walk the seed in the other index (batmis1's\n"
        "            way) instead of reading it off the bidirectional range\n"
        "  --own-row  locate single rows from their own index's row only (batmis1's\n"
        "            way), not from either index's\n"
        "  --lookup  depth k of the start-of-search tables, 3-14 (default 6): 8*4^k bytes\n"
        "            per index; deeper is faster once there are more than about 4^k/3 reads\n"
        "            (see benchmark/lookup_table.md)\n"
        "  --batch   search B reads at a time, 1-64, so that their memory accesses\n"
        "            overlap (default 16, or $BATMIS_BATCH). auto times B = 4, 8, 16\n"
        "            on the first reads, keeps the fastest and prints the times to\n"
        "            stderr: a check of the default on a new machine or workload\n"
        "            (see benchmark/prefetch_fix.md)\n"
        "  -t, --threads  search with N threads (default 1, 0 = all CPUs); the output\n"
        "            does not depend on N (see benchmark/threads.md)\n"
        "Output, one line per SA range or located hit:\n"
        "  exact:      query_id\tF\tSA_start\tSA_end\thit_count\n"
        "  1 or 2 mismatches: query_id\tF|B\tSA_start\tSA_end\thit_count\t1|2\n"
        "  F = range in the reverse index (rev<genome_prefix>), B = range in the forward index\n"
        "  P = located hit: query_id\tP\tpos\tpos\t1[\t1], pos = 0-based start in the genome\n"
        "  query_id\tF\t0\t0\t0 when a read has no hits\n",
        prog);
}

}  // namespace

int main(int argc, char** argv) {
    std::string queryFile;
    std::string genomePrefix;
    std::string outputFile;
    unsigned maxHits = 1;
    int maxMismatches = 0;
    bool quiet = false;
    bool timing = false;
    bool useRef = true;
    int batchSize = DEFAULT_BATCH;  // --batch; 0: auto
    // "auto" or 1-MAX_BATCH, from --batch or $BATMIS_BATCH.
    auto Parse_Batch = [&](const char* Text, const char* From) {
        batchSize = std::strcmp(Text, "auto") == 0 ? 0 : std::atoi(Text);
        if (std::strcmp(Text, "auto") == 0 || (batchSize >= 1 && batchSize <= static_cast<int>(MAX_BATCH))) return true;
        std::fprintf(stderr, "%s must be auto or 1-%u\n", From, MAX_BATCH);
        return false;
    };
    if (const char* Env = std::getenv("BATMIS_BATCH"))
        if (!Parse_Batch(Env, "BATMIS_BATCH")) return 1;
    int threads = 1;    // -t; 0: all CPUs
    int samFormat = 0;  // 1 SAM, 2 BAM (--sam / --bam, or the -o extension)
    bool forwardOnly = false;

    for (int i = 1; i < argc; i++) {
        if (std::strcmp(argv[i], "-q") == 0 && i + 1 < argc) {
            queryFile = argv[++i];
        } else if (std::strcmp(argv[i], "-g") == 0 && i + 1 < argc) {
            genomePrefix = argv[++i];
        } else if (std::strcmp(argv[i], "-o") == 0 && i + 1 < argc) {
            outputFile = argv[++i];
        } else if (std::strcmp(argv[i], "-m") == 0 && i + 1 < argc) {
            maxHits = static_cast<unsigned>(std::strtoul(argv[++i], nullptr, 10));
        } else if (std::strcmp(argv[i], "-n") == 0 && i + 1 < argc) {
            ++i;
            if (std::strlen(argv[i]) != 1 || argv[i][0] < '0' || argv[i][0] > '0' + MAX_MISMATCH_SEARCH) {
                Usage(argv[0]);
                return 1;
            }
            maxMismatches = argv[i][0] - '0';
        } else if (std::strcmp(argv[i], "--quiet") == 0 || std::strcmp(argv[i], "-Q") == 0) {
            quiet = true;
        } else if (std::strcmp(argv[i], "--timing") == 0) {
            timing = true;
        } else if (std::strcmp(argv[i], "--batch") == 0 && i + 1 < argc) {
            if (!Parse_Batch(argv[++i], "--batch")) return 1;
        } else if ((std::strcmp(argv[i], "-t") == 0 || std::strcmp(argv[i], "--threads") == 0) && i + 1 < argc) {
            threads = std::atoi(argv[++i]);
            if (threads < 0 || threads > 1024) {
                std::fprintf(stderr, "-t must be 0-1024\n");
                return 1;
            }
            if (threads == 0) threads = static_cast<int>(std::max(1U, std::thread::hardware_concurrency()));
        } else if (std::strcmp(argv[i], "--lookup") == 0 && i + 1 < argc) {
            LOOKUP_DEPTH = std::atoi(argv[++i]);
            if (LOOKUP_DEPTH < 3 || LOOKUP_DEPTH > 14) {
                std::fprintf(stderr, "--lookup must be 3-14\n");
                return 1;
            }
        } else if (std::strcmp(argv[i], "--no-ref") == 0) {
            useRef = false;
        } else if (std::strcmp(argv[i], "--switch-walk") == 0) {
            Switch_Walk = true;
        } else if (std::strcmp(argv[i], "--own-row") == 0) {
            Both_Rows = false;
        } else if (std::strcmp(argv[i], "--sam") == 0) {
            samFormat = 1;
        } else if (std::strcmp(argv[i], "--bam") == 0) {
            samFormat = 2;
        } else if (std::strcmp(argv[i], "--forward-only") == 0) {
            forwardOnly = true;
        } else if (std::strcmp(argv[i], "-h") == 0 || std::strcmp(argv[i], "--help") == 0) {
            Usage(argv[0]);
            return 0;
        } else {
            Usage(argv[0]);
            return 1;
        }
    }

    if (queryFile.empty() || genomePrefix.empty()) {
        Usage(argv[0]);
        return 1;
    }
    auto endsWith = [&](const char* suffix) {
        const std::size_t n = std::strlen(suffix);
        return outputFile.size() >= n && outputFile.compare(outputFile.size() - n, n, suffix) == 0;
    };
    if (!samFormat && endsWith(".bam")) samFormat = 2;
    if (!samFormat && endsWith(".sam")) samFormat = 1;

    using Clock = std::chrono::steady_clock;
    const auto totalStart = Clock::now();
    double indexLoadSeconds = 0.0;
    double tableBuildSeconds = 0.0;
    double inputParseSeconds = 0.0;
    double inputReadSeconds = 0.0;
    double searchSeconds = 0.0;
    double outputSeconds = 0.0;

    std::string revBwt;
    std::string revOcc;
    std::string fwBwt;
    std::string fwOcc;
    std::string revSa;
    std::string fwSa;
    BuildIndexPaths(genomePrefix, revBwt, revOcc, fwBwt, fwOcc, &revSa, &fwSa);
    const std::string pacName = genomePrefix + ".pac";
    if (useRef) {
        for (const std::string* name : std::initializer_list<const std::string*>{&revSa, &fwSa, &pacName}) {
            if ((name != &fwSa || maxMismatches > 0) && !FileExists(*name)) {
                std::fprintf(stderr, "Warning: %s not found; walking single-row paths instead (as --no-ref)\n", name->c_str());
                useRef = false;
                break;
            }
        }
    }

    const bool loadSa = useRef || samFormat;  // SAM/BAM: positions need the SA samples
    if (samFormat) {
        for (const std::string* name : {&revSa, &fwSa}) {
            if ((name != &fwSa || maxMismatches > 0) && !FileExists(*name)) {
                std::fprintf(stderr, "SAM/BAM output needs %s\n", name->c_str());
                return 1;
            }
        }
    }

    const auto loadStart = Clock::now();
    revfmi = initFMI(revBwt.c_str(), revOcc.c_str(), loadSa ? revSa.c_str() : nullptr);
    if (!revfmi) {
        std::fprintf(stderr, "Failed to load reverse index: %s / %s\n", revBwt.c_str(), revOcc.c_str());
        return 1;
    }
    Need_Backward = maxMismatches > 0;
    if (!useRef) Both_Rows = false;  // the other row is used to locate only
    if (maxMismatches > 0 || Both_Rows) {
        fwfmi = BWTLoad(mmPool, fwBwt.c_str(), fwOcc.c_str(), loadSa ? fwSa.c_str() : nullptr, nullptr, nullptr, nullptr);
        if (!fwfmi) {
            std::fprintf(stderr, "Failed to load forward index: %s / %s\n", fwBwt.c_str(), fwOcc.c_str());
            return 1;
        }
    }
    if (useRef) {
        if (revfmi->saInterval != SAINTERVAL || (fwfmi && fwfmi->saInterval != SAINTERVAL)) {
            std::fprintf(stderr, "Warning: SA sample interval is not %d; walking single-row paths instead\n", SAINTERVAL);
            useRef = false;
        } else if (!LoadPac(pacName)) {
            std::fprintf(stderr, "Warning: cannot read %s; walking single-row paths instead\n", pacName.c_str());
            useRef = false;
        }
    }
    Use_Ref = useRef;
    Track_Pos = samFormat && !useRef && revfmi->saInterval == SAINTERVAL && (!fwfmi || fwfmi->saInterval == SAINTERVAL);
    Text_Length = revfmi->textLength;

    Sam_Output sam;
    if (samFormat) {
        std::string error, commandLine = argv[0];
        for (int i = 1; i < argc; i++) commandLine += std::string(" ") + argv[i];
        if (!sam.Map.Load(genomePrefix, error) ||
            (!quiet && (sam.Writer.Program = "batmis2", !sam.Writer.Open(outputFile.empty() ? "-" : outputFile, samFormat == 2, sam.Map, commandLine, error, threads)))) {
            std::fprintf(stderr, "%s\n", error.c_str());
            return 1;
        }
        if (sam.Map.Text_Length() != Text_Length) {
            std::fprintf(stderr, "%s.ann.location does not match the index (%llu bases, index %lu)\n", genomePrefix.c_str(),
                         static_cast<unsigned long long>(sam.Map.Text_Length()), Text_Length);
            return 1;
        }
        sam.Reverse = !forwardOnly;
        sam.Write = !quiet;
    }
    indexLoadSeconds = std::chrono::duration<double>(Clock::now() - loadStart).count();

    Read_File in(queryFile);
    if (!in.Is_Open()) {
        std::fprintf(stderr, "Cannot open query file: %s\n", queryFile.c_str());
        return 1;
    }

    std::ostream* out = &std::cout;
    std::ofstream outFile;
    if (!outputFile.empty() && !samFormat) {
        outFile.open(outputFile.c_str());
        if (!outFile) {
            std::fprintf(stderr, "Cannot open output file: %s\n", outputFile.c_str());
            return 1;
        }
        out = &outFile;
    }

    uint64_t records = 0;
    uint64_t reportedHits = 0;

    if (!Run_Batched(in, *out, batchSize, threads, maxMismatches, maxHits, quiet, records, reportedHits, tableBuildSeconds,
                     inputParseSeconds, inputReadSeconds, searchSeconds, outputSeconds, samFormat ? &sam : nullptr)) {
        return 1;
    }
    if (in.Failed()) {
        std::fprintf(stderr, "Error reading %s (corrupt or truncated compressed data?)\n", queryFile.c_str());
        return 1;
    }

    if (timing) {
        const double totalSeconds = std::chrono::duration<double>(Clock::now() - totalStart).count();
        std::cerr << std::fixed << std::setprecision(6)
                  << "batmis2 timing: index_load=" << indexLoadSeconds
                  << "s table_build=" << tableBuildSeconds
                  << "s input_parse=" << inputParseSeconds
                  << "s input_read=" << inputReadSeconds
                  << "s search=" << searchSeconds
                  << "s output_format=" << outputSeconds
                  << "s total=" << totalSeconds
                  << "s records=" << records
                  << " reported_hits=" << reportedHits << ' ' << Batch_Report << '\n';
#ifdef BATMAN_BATCH_STATS
        if (Batch_Stats.Rounds) {
            const char* subs[] = {"pop", "wide", "single", "locate", "located", "switch"};
            const char* tasks[] = {"srch.", "ext.", "win.", "sw."};
            std::cerr << "batch: rounds=" << Batch_Stats.Rounds
                      << " occupancy=" << double(Batch_Stats.Lane_Rounds) / Batch_Stats.Rounds << "/" << Batch_Stats.Lanes
                      << " drain_rounds=" << Batch_Stats.Drain_Rounds << " steps_per_record="
                      << double(Batch_Stats.Lane_Rounds) / records << '\n';
            for (unsigned p = 0; p < MM_PROGRAM_MAX; p++) {
                std::cerr << "  program step " << p << ":";
                for (int e = 0; e < 4; e++)
                    for (int u = 0; u < 6; u++)
                        if (Batch_Stats.Steps[p][e][u])
                            std::cerr << ' ' << tasks[e] << subs[u] << '='
                                      << double(Batch_Stats.Steps[p][e][u]) / records;
                if (Batch_Stats.Skipped[p]) std::cerr << " skipped=" << double(Batch_Stats.Skipped[p]) / records;
                std::cerr << '\n';
            }
        }
        if (Decode_Stats.Rows)
            std::cerr << "decode: seconds=" << Decode_Stats.Seconds << " rows=" << Decode_Stats.Rows << " read_off=" << Decode_Stats.Direct
                      << " group_lookups=" << Decode_Stats.Range_Steps << " span_steps=" << Decode_Stats.Span_Steps
                      << " row_steps=" << Decode_Stats.Walk_Steps << " lookups_per_row="
                      << double(Decode_Stats.Range_Steps + 2 * Decode_Stats.Span_Steps + Decode_Stats.Walk_Steps) / Decode_Stats.Rows
                      << '\n';
#endif
    }


    return 0;
}
