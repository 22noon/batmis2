// batman_zero_min: exact-match search using the native BatMis FM index
// (rev<genome>.bwt/.fmv built by index/bwt_build_index). With -n 1 it also
// finds one-mismatch hits, using the forward index <genome>.bwt/.fmv too.
// By default it also loads the SA samples (.sa) and the packed genome (.pac):
// a single-row path that reaches a sampled row is finished by comparing the read
// with the genome there (see Locate). --no-ref skips them to save memory.
// The BWA-index variant lives in batman_zero_min_bwa.cpp.
// Optional build variant: -DBATMAN_ZERO_MIN_USE_PAIRED_RANK.
#include <algorithm>
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
#include <chrono>

extern "C" {
#include "MemManager.h"
#include "BWT.h"
}

namespace {

constexpr int MAXTAG = 4096;
constexpr int MAX_MISMATCHES_BOUND = 16;
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
    int Tag;
    unsigned char Mismatch_Pos[MAX_MISMATCHES_BOUND];
    unsigned Mismatch_Char;
    unsigned long Pos;           // located: forward-genome start of the read
};

struct Range {
    unsigned Start;
    unsigned End;
    int Label;
};

// Output of Extend_Right/Extend_Left. The caller clears it before each search.
struct ForwardResults {
    std::vector<SARange> hits;        // exactly Count mismatches
    std::vector<SARange> mismatches;  // deferred parents at Count mismatches: resume at Level, skipping Forbidden_Base
    unsigned long hitCount = 0;       // locations in hits
};

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
};
std::vector<Lookup_Entry> Forward_Lookup;
std::vector<Lookup_Entry> Backward_Lookup;
unsigned Branch_Characters[4] = {0, 0, 0, 0};
SARange Branch_Ranges[4];

// Reference for located singletons (see Locate). Off with --no-ref.
bool Use_Ref = false;
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
bool Locate(const ActiveBWT* fmi, SARange& Tag, bool Right, int Offset, int Len, int Done) {
    const long P = fmi->saValue[Tag.Start / SAINTERVAL];
    const long S = Right ? static_cast<long>(Text_Length) - (Offset + Done) - P : P - (Offset + Len - Done);
    Tag.Skip = 1;
    Tag.End = Tag.Start;
    if (S < 0 || S + STRINGLENGTH > static_cast<long>(Text_Length)) return false;
    Tag.Pos = S;
    return true;
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

// Search_Right/Search_Left at a sampled row: matches if the rest of the segment does.
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

// Extend_Right/Extend_Left for a located Tag, whose first Level-1 segment bases are
// matched with Mismatches mismatches. As the walk would: a hit at exactly Count
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

// Prefetches the occurrence and BWT words BWTOccValue reads for row Index
// (as BWTPrefetchOccValueExplicit/BWTPrefetchBWT in index/BWT.c).
inline void Prefetch_Row(const ActiveBWT* fmi, unsigned long Index) {
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
    const unsigned long First = Start - (Start > fmi->inverseSa0);
    const unsigned long Last = End + 1 - (End + 1 > fmi->inverseSa0);	// exclusive
    const unsigned* Code = fmi->bwtCode;
    if (First >= Last) { Count[0] = Count[1] = Count[2] = Count[3] = 0; return; }

    // Characters are packed high bits first: character k occupies bits 31-2k..30-2k.
    // Three bit-planes instead of four base comparisons: low bit of the pair set
    // means C or T, high bit set means G or T, both means T; A is what remains.
    // Only the first and last words are masked; whole words in between are counted
    // two at a time.
    const unsigned long W0 = First / 16, W1 = (Last - 1) / 16;
    const unsigned Head = ~0u >> (2 * (First % 16));			// characters k >= First%16
    const unsigned Tail = ~0u << (30 - 2 * ((Last - 1) % 16));	// characters k <= (Last-1)%16
    unsigned long CT = 0, GT = 0, T = 0;
    auto Word32 = [&](unsigned W, unsigned M) {
        M &= 0x55555555u;
        const unsigned L = W & M, H = (W >> 1) & M;
        CT += __builtin_popcount(L); GT += __builtin_popcount(H); T += __builtin_popcount(L & H);
    };
    if (W0 == W1) Word32(Code[W0], Head & Tail);
    else {
        Word32(Code[W0], Head);
        unsigned long W = W0 + 1;
        for (; W + 2 <= W1; W += 2) {
            uint64_t X; std::memcpy(&X, Code + W, 8);
            // (X >> 1) moves bit 32 into bit 31, an odd bit that the mask drops.
            const uint64_t L = X & 0x5555555555555555ULL, H = (X >> 1) & 0x5555555555555555ULL;
            CT += __builtin_popcountll(L); GT += __builtin_popcountll(H); T += __builtin_popcountll(L & H);
        }
        if (W < W1) Word32(Code[W], ~0u);
        Word32(Code[W1], Tail);
    }
    const unsigned long N = Last - First;
    Count[3] = T; Count[1] = CT - T; Count[2] = GT - T; Count[0] = N - CT - GT + T;
}

// Narrows Start..End to the rows preceded by newChar (Start=0 if none). Ranges
// with End-Start below BRANCHTHRESHOLD count their characters directly (Count_Range)
// instead of a second occurrence lookup; see benchmark/branch_threshold.md.
void Get_SARange(long newChar, unsigned long& start, unsigned long& end, ActiveBWT* fmi) {
    if (end - start < BRANCHTHRESHOLD) {
        unsigned count[4];
        Count_Range(fmi, start, end, count);
        if (!count[newChar]) {
            start = 0;
            return;
        }
        start = fmi->cumulativeFreq[newChar] + BWTOccValue(fmi, start, newChar) + 1;
        end = start + count[newChar] - 1;
        return;
    }
#ifdef BATMAN_ZERO_MIN_USE_PAIRED_RANK
    unsigned int occValue[2];
    BWTOccValueTwoIndex(fmi, static_cast<unsigned int>(start),
                        static_cast<unsigned int>(end + 1),
                        static_cast<unsigned int>(newChar), occValue);
    start = fmi->cumulativeFreq[newChar] + occValue[0] + 1;
    end = fmi->cumulativeFreq[newChar] + occValue[1];
#else
    start = fmi->cumulativeFreq[newChar] + BWTOccValue(fmi, start, newChar) + 1;
    end = fmi->cumulativeFreq[newChar] + BWTOccValue(fmi, end + 1, newChar);
#endif
    if (end < start) {
        start = 0;
    }
}

inline void Get_SARange(long newChar, SARange& range, ActiveBWT* fmi) {
    Get_SARange(newChar, range.Start, range.End, fmi);
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
    return range;
}

void Build_Preindex_Forward(Range range, int level, int bit) {
    if (LOOKUPSIZE == level) {
        assert(range.Start <= range.End);
        range.Label = range.Label | (bit << (2 * (level - 1)));
        if (range.Start) {
            unsigned long start = range.Start;
            unsigned long end = range.End;
            Get_SARange(bit, start, end, revfmi);
            range.Start = static_cast<unsigned>(start);
            range.End = static_cast<unsigned>(end);
        }
        Forward_Lookup[range.Label] = {range.Start, range.End};
        assert(range.Start <= range.End);
    } else {
        assert(range.Start <= range.End);
        range.Label = range.Label | (bit << (2 * (level - 1)));
        if (range.Start) {
            unsigned long start = range.Start;
            unsigned long end = range.End;
            Get_SARange(bit, start, end, revfmi);
            range.Start = static_cast<unsigned>(start);
            range.End = static_cast<unsigned>(end);
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
        unsigned long start = range.Start;
        unsigned long end = range.End;
        Get_SARange(bit, start, end, fwfmi);
        range.Start = static_cast<unsigned>(start);
        range.End = static_cast<unsigned>(end);
    }
    if (LOOKUPSIZE == level) {
        Backward_Lookup[range.Label] = {range.Start, range.End};
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
    if (fwfmi) Backward_Lookup.assign(1ul << (2 * LOOKUPSIZE), Lookup_Entry{});

    Range lRange{};
    for (int i = 0; i < 4; i++) {
        lRange = InitialBaseRange(fmi, i);
        lRange.Label = i;
        for (int j = 0; j < 4; j++) {
            Build_Preindex_Forward(lRange, 2, j);
        }
    }

    if (fwfmi) {
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
        tag.Skip++;
        tag.End = tag.Start;
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
            tag.Start = fmi->cumulativeFreq[base] + BWTOccValue(fmi, tag.Start, base) + 1;
            if (!tag.Skip && tag.Start % SAINTERVAL == 0) {
                if (Use_Ref) {
                    if (!Search_Located(fmi, tag, true, offset, stringLength, tag.Level)) return false;
                    tag.Level = stringLength;
                    return true;
                }
                tag.Skip++;
                tag.End = tag.Start;
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
		BWTAllOccValue(fmi, Tag.Start, Occ_Start);
		Count_Range(fmi, Tag.Start, Tag.End, Count);
		for (int Branch=0;Branch<4;Branch++)
		{
			Branch_Characters[Branch]=Count[Branch];
			if (Count[Branch])
			{
				Branch_Ranges[Branch].Start = fmi->cumulativeFreq[Branch] + Occ_Start[Branch] + 1;
				Branch_Ranges[Branch].End = Branch_Ranges[Branch].Start + Count[Branch]-1;// Calculate SAranges
			}
		}
	}
	else
	{
		// Occurrences of all four bases at both ends in one call (one decode per end).
		unsigned int ALIGN_16 Occ_Start[ALPHABET_SIZE], Occ_End[ALPHABET_SIZE];
		BWTAllOccValueTwoIndex(fmi, Tag.Start, Tag.End+1, Occ_Start, Occ_End);
		for (int Branch=0;Branch<4;Branch++)
		{
				Branch_Ranges[Branch].Start = fmi->cumulativeFreq[Branch] + Occ_Start[Branch] + 1;
				Branch_Ranges[Branch].End = fmi->cumulativeFreq[Branch] + Occ_End[Branch];
				if(!(Branch_Ranges[Branch].End<Branch_Ranges[Branch].Start)) Branch_Characters[Branch]=1;
		}

	}
}
void Extend_Right_Singleton(std::string_view Current_Tag,int Offset,struct SARange & Tag,int Count,BWT *fmi,ForwardResults& Out)
{
	const int StringLength=static_cast<int>(Current_Tag.size());
	unsigned long Index,Now;
	if (Tag.Start==0) return;
	if (Use_Ref && Tag.Skip) {Extend_Located(Tag,true,Offset,StringLength,Count,Out);return;}//located earlier
	if(Tag.Start % SAINTERVAL == 0 && !Tag.Skip)
	{
		if (Use_Ref) {if (Locate(fmi,Tag,true,Offset,StringLength,Tag.Level-1)) Extend_Located(Tag,true,Offset,StringLength,Count,Out);return;}
		Tag.Skip++;
		Tag.End=Tag.Start;
	}

	for(;;)
	{
		
		if (Tag.Start==fmi->inverseSa0) return;//$ row: no base to extend with
		Index=Tag.Start-(Tag.Start>fmi->inverseSa0);//adjust for missing $
		Now=fmi->bwtCode[(Index) / 16] << (((Index) % 16) * 2)>> (BITS_IN_WORD - 2);
		if (Now==Tag.Forbidden_Base) return;//already tried by the pass that deferred this range
		const unsigned char Expected=Current_Tag[Tag.Level-1];
		const unsigned Base=Read_Base(Expected);
		if (Now!=Expected && Tag.Mismatches==Count)//at threshold: defer the parent
		{
			if (!Tag.Skip) Tag.End=Tag.Start;
			Tag.Forbidden_Base=Expected;
			Out.mismatches.push_back(Tag);
			return;
		}
		Tag.Forbidden_Base=NO_BASE;
		// On a match step with the read's base (see Read_Base).
		if (Now==Expected) Tag.Start = fmi->cumulativeFreq[Base] + BWTOccValue(fmi, Tag.Start, Base) + 1;
		else Tag.Start = fmi->cumulativeFreq[Now] + BWTOccValue(fmi, Tag.Start, Now) + 1;

		if(!Tag.Skip && Tag.Start % SAINTERVAL == 0)
		{
			Tag.Skip++;
			Tag.End=Tag.Start;
		}

		if (Expected!=Now)
		{

			Tag.Mismatch_Char=Tag.Mismatch_Char | (Now<<Tag.Mismatches*2);
			Tag.Mismatch_Pos[Tag.Mismatches]=(Offset+Tag.Level-1);
			Tag.Mismatches++;
			
		}
		if (Use_Ref && Tag.Skip)//just reached a sampled row
		{
			Tag.Level++;
			if (Locate(fmi,Tag,true,Offset,StringLength,Tag.Level-1)) Extend_Located(Tag,true,Offset,StringLength,Count,Out);
			return;
		}

		if (Tag.Mismatches<=Count)
		{
			if(Tag.Level== StringLength)
			{
				if(Tag.Mismatches==Count)//avoid printing exact matches...
				{
					//if (Tag.Skip) Tag.Start= Tag.End; else Tag.End=Tag.Start;
					if (!Tag.Skip) Tag.End=Tag.Start;
					Out.hits.push_back(Tag);
					Out.hitCount += Tag.Skip ? 1 : Tag.End-Tag.Start+1;
				}
				return;
			}
			else {Tag.Level++;continue;}
		}
	}
}

// Extends Tag to the right over Current_Tag[Level-1] allowing up to Count
// mismatches. Offset is Current_Tag's position in the read (for Mismatch_Pos).
// Ranges with exactly Count mismatches go to Out.hits. Ranges already at Count
// mismatches only extend with the read's base; where that loses rows, the
// parent goes to Out.mismatches with Forbidden_Base set, to be resumed with a
// larger Count (as minimal_bwt/src/mismatch_search.cpp).
// Stops once Out.hitCount reaches MaxHits (0 = no limit). Kept out of line, as Extend_Left.
__attribute__((noinline)) void Extend_Right(std::string_view Current_Tag,int Offset,const struct SARange & Tag,int Count,BWT *fmi,ForwardResults& Out,unsigned MaxHits)
{
	if (!Tag.Start) return;
	const int StringLength=static_cast<int>(Current_Tag.size());
	std::vector<SARange> Stack;
	Stack.push_back(Tag);
	struct SARange Range,Temp_Range;
	while(!Stack.empty())//While Stack non-empty....
	{
		Range=Stack.back();
		Stack.pop_back();		//Pop the range
		if (Range.End==Range.Start || Range.Skip)//does this SArange have only one branch?
		{
			Extend_Right_Singleton(Current_Tag,Offset,Range,Count,fmi,Out);
			if(MaxHits && Out.hitCount>=MaxHits) return;
		}
		else if (Range.Mismatches==Count)//at threshold: exact step only
		{
			const unsigned char Expected=Current_Tag[Range.Level-1];
			Temp_Range=Range;
			Get_SARange(Expected,Temp_Range,fmi);
			// A parent holding the $ row also differs in size; deferring it is harmless (its resume finds nothing).
			if (!Temp_Range.Start || Temp_Range.End-Temp_Range.Start!=Range.End-Range.Start)
			{
				Range.Forbidden_Base=Expected;//other bases remain: defer the parent
				Out.mismatches.push_back(Range);
			}
			if (!Temp_Range.Start) continue;
			Temp_Range.Forbidden_Base=NO_BASE;
			if(Temp_Range.Level== StringLength)
			{
				Out.hits.push_back(Temp_Range);
				Out.hitCount += Temp_Range.Skip ? 1 : Temp_Range.End-Temp_Range.Start+1;
				if (MaxHits && Out.hitCount>=MaxHits) return;
			}
			else
			{
				Temp_Range.Level++;
				Stack.push_back(Temp_Range);
			}
		}
		else
		{
			Branch_Detect(Range,fmi,Offset);//One_Branch(Range,revfmi);
			for(int Branch=0;Branch<4;Branch++)
			{
				if (Branch_Characters[Branch] && Branch!=Range.Forbidden_Base)//This character actually branches
				{
					Temp_Range=Range;
					Temp_Range.Forbidden_Base=NO_BASE;
					Temp_Range.Start = Branch_Ranges[Branch].Start;//fmi->cumulativeFreq[Branch] + BWTOccValue(fmi, Temp_Range.Start, Branch) + 1;
					Temp_Range.End = Branch_Ranges[Branch].End;//Temp_Range.Start + Branch_Characters[Branch]-1;// Calculate SAranges

					if (Current_Tag[Temp_Range.Level-1]!=Branch)
					{

						Temp_Range.Mismatch_Char=Temp_Range.Mismatch_Char | (Branch<<Temp_Range.Mismatches*2);
						Temp_Range.Mismatch_Pos[Temp_Range.Mismatches]=(Offset+Temp_Range.Level-1);
						Temp_Range.Mismatches++;

					}


					if (Temp_Range.Mismatches<=Count)//we are guaranteed a valid SA range, check only for mismatches
					{
						if(Temp_Range.Level== StringLength)
						{
							if (Temp_Range.Mismatches == Count) //dont print exact matches
							{
								Out.hits.push_back(Temp_Range);
								Out.hitCount += Temp_Range.Skip ? 1 : Temp_Range.End-Temp_Range.Start+1;
								if (MaxHits && Out.hitCount>=MaxHits) return;
							}
							else continue;
						}
						else
						{

							Temp_Range.Level++;
							Stack.push_back(Temp_Range);//Push range
						}
					}
				} 
			}
		}
	}
	return;
}

// Search_Left along a single SA path (Tag.End==Tag.Start or Tag.Skip).
// Kept out of line, as Search_Right_Singleton.
__attribute__((noinline)) void Search_Left_Singleton(std::string_view Current_Tag,int Offset,struct SARange & Tag,BWT *fmi)
{
	const int StringLength=static_cast<int>(Current_Tag.size());
	int Level=Tag.Level;
	unsigned long Index,Now;

	if(Tag.Start % SAINTERVAL == 0 && !Tag.Skip)
	{
		if (Use_Ref)
		{
			if (!Search_Located(fmi,Tag,false,Offset,StringLength,Level-1)) Tag.Start=0;//Level: the failure depth
			return;
		}
		Tag.Skip++;
		Tag.End=Tag.Start;
	}
	for(;;)
	{
		if (Tag.Start==fmi->inverseSa0) {Tag.Start=0;Tag.Level=Level;return;}//$ row: no base to extend with
		Index=Tag.Start-(Tag.Start>fmi->inverseSa0);//adjust for missing $
		Now=fmi->bwtCode[(Index) / 16] << (((Index) % 16) * 2)>> (BITS_IN_WORD - 2);//FMIBwtValue(fmi,Index);
		const unsigned long Expected=static_cast<unsigned char>(Current_Tag[StringLength-Level]);
		const unsigned Base=Read_Base(Expected);
		if (Expected == Now)
		{
			Tag.Start = fmi->cumulativeFreq[Base] + BWTOccValue(fmi, Tag.Start, Base) + 1;

			if(!Tag.Skip && Tag.Start % SAINTERVAL == 0)
			{
				if (Use_Ref)
				{
					if (!Search_Located(fmi,Tag,false,Offset,StringLength,Level)) Tag.Start=0;
					return;
				}
				Tag.Skip++;Tag.End=Tag.Start;
			}

			if(Level== StringLength)//no need to print as we are still halfway..
			{
				if (!Tag.Skip) Tag.End=Tag.Start;
				return;
			}
			else {Level++;continue;}
		}
		else
		{
			Tag.Start=0;Tag.Level=Level;//the failure depth, as Search_Left
			return;
		}
	}
}

// Extends Tag to the left over Current_Tag[size-Level] with no mismatches.
// On return Tag.Start==0 if there is no match, otherwise Tag holds the SA range.
void Search_Left(std::string_view Current_Tag,int Offset,struct SARange & Tag,BWT *fmi)
{
	const int StringLength=static_cast<int>(Current_Tag.size());

	if (!Tag.Start) return;

	for(;;)
	{
		if(Tag.End==Tag.Start || Tag.Skip)//Only one branch?
		{
			Search_Left_Singleton(Current_Tag,Offset,Tag,fmi);
			return;
		}
		else
		{
			Get_SARange(Current_Tag[StringLength-Tag.Level],Tag,fmi);

			if (Tag.Start!=0)
			{
				if(Tag.Level== StringLength)
				{
					return;
				}
				else {Tag.Level++;continue;}
			}
			else
			{
				return;
			}

		}
	}
}

void Extend_Left_Singleton(std::string_view Current_Tag,int Offset,struct SARange & Tag,int Count,BWT *fmi,ForwardResults& Out)
{
	const int StringLength=static_cast<int>(Current_Tag.size());
	unsigned long Index,Now;
	if (Tag.Start==0) return;
	if (Use_Ref && Tag.Skip) {Extend_Located(Tag,false,Offset,StringLength,Count,Out);return;}//located earlier
	if(Tag.Start % SAINTERVAL == 0 && !Tag.Skip)
	{
		if (Use_Ref) {if (Locate(fmi,Tag,false,Offset,StringLength,Tag.Level-1)) Extend_Located(Tag,false,Offset,StringLength,Count,Out);return;}
		Tag.Skip++;
		Tag.End=Tag.Start;
	}

	for(;;)
	{
		if (Tag.Start==fmi->inverseSa0) return;//$ row: no base to extend with
		Index=Tag.Start-(Tag.Start>fmi->inverseSa0);//adjust for missing $
		Now=fmi->bwtCode[(Index) / 16] << (((Index) % 16) * 2)>> (BITS_IN_WORD - 2);//FMIBwtValue(fmi,Index);
		if (Now==Tag.Forbidden_Base) return;//already tried by the pass that deferred this range
		const unsigned char Expected=Current_Tag[StringLength-Tag.Level];
		const unsigned Base=Read_Base(Expected);
		if (Now!=Expected && Tag.Mismatches==Count)//at threshold: defer the parent
		{
			if (!Tag.Skip) Tag.End=Tag.Start;
			Tag.Forbidden_Base=Expected;
			Out.mismatches.push_back(Tag);
			return;
		}
		Tag.Forbidden_Base=NO_BASE;
		// On a match step with the read's base (see Read_Base).
		if (Now==Expected) Tag.Start = fmi->cumulativeFreq[Base] + BWTOccValue(fmi, Tag.Start, Base) + 1;
		else Tag.Start = fmi->cumulativeFreq[Now] + BWTOccValue(fmi, Tag.Start, Now) + 1;

		if(!Tag.Skip && Tag.Start % SAINTERVAL == 0)
		{
			Tag.Skip++;
			Tag.End=Tag.Start;
		}

		if (Expected!=Now)
		{
			Tag.Mismatch_Char=Tag.Mismatch_Char | (Now<<Tag.Mismatches*2);
			Tag.Mismatch_Pos[Tag.Mismatches]=(Offset+StringLength-Tag.Level);
			Tag.Mismatches++;

		}
		if (Use_Ref && Tag.Skip)//just reached a sampled row
		{
			Tag.Level++;
			if (Locate(fmi,Tag,false,Offset,StringLength,Tag.Level-1)) Extend_Located(Tag,false,Offset,StringLength,Count,Out);
			return;
		}

		if (Tag.Mismatches<=Count)
		{
			if(Tag.Level== StringLength)
			{
				if(Tag.Mismatches==Count)
				{
					//if (Tag.Skip) Tag.Start= Tag.End; else Tag.End=Tag.Start;
					if (!Tag.Skip) Tag.End=Tag.Start;
					Out.hits.push_back(Tag);
					Out.hitCount += Tag.Skip ? 1 : Tag.End-Tag.Start+1;
				}
				return;
			}
			else {Tag.Level++;continue;}
		}
	}
}

// Extends Tag to the left over Current_Tag[size-Level] allowing up to Count
// mismatches. Offset and results as for Extend_Right.
// Kept out of line: inlined into the read loop, r24 searches ran ~15% slower.
__attribute__((noinline)) void Extend_Left(std::string_view Current_Tag,int Offset,const struct SARange & Tag,int Count,BWT *fmi,ForwardResults& Out,unsigned MaxHits)
{
	if (!Tag.Start) return;
	const int StringLength=static_cast<int>(Current_Tag.size());
	std::vector<SARange> Stack;
	Stack.push_back(Tag);
	struct SARange Range,Temp_Range;
	while(!Stack.empty())//While Stack non-empty....
	{
		Range=Stack.back();
		Stack.pop_back();	//Pop the range
		if (Range.End==Range.Start || Range.Skip)//does this SArange have only one branch?
		{
			Extend_Left_Singleton(Current_Tag,Offset,Range,Count,fmi,Out);
			if(MaxHits && Out.hitCount>=MaxHits) return;
		}
		else if (Range.Mismatches==Count)//at threshold: exact step only
		{
			const unsigned char Expected=Current_Tag[StringLength-Range.Level];
			Temp_Range=Range;
			Get_SARange(Expected,Temp_Range,fmi);
			// A parent holding the $ row also differs in size; deferring it is harmless (its resume finds nothing).
			if (!Temp_Range.Start || Temp_Range.End-Temp_Range.Start!=Range.End-Range.Start)
			{
				Range.Forbidden_Base=Expected;//other bases remain: defer the parent
				Out.mismatches.push_back(Range);
			}
			if (!Temp_Range.Start) continue;
			Temp_Range.Forbidden_Base=NO_BASE;
			if(Temp_Range.Level== StringLength)
			{
				Out.hits.push_back(Temp_Range);
				Out.hitCount += Temp_Range.Skip ? 1 : Temp_Range.End-Temp_Range.Start+1;
				if (MaxHits && Out.hitCount>=MaxHits) return;
			}
			else
			{
				Temp_Range.Level++;
				Stack.push_back(Temp_Range);
			}
		}
		else
		{
			Branch_Detect(Range,fmi,Offset);
			for(int Branch=0;Branch<4;Branch++)
			{
				if (Branch_Characters[Branch] && Branch!=Range.Forbidden_Base)//This character actually branches
				{
					Temp_Range=Range;//adjust
					Temp_Range.Forbidden_Base=NO_BASE;
					Temp_Range.Start = Branch_Ranges[Branch].Start;//fmi->cumulativeFreq[Branch] + BWTOccValue(fmi, Temp_Range.Start, Branch) + 1;
					Temp_Range.End = Branch_Ranges[Branch].End;//Temp_Range.Start + Branch_Characters[Branch]-1;// Calculate SAranges

					if (Current_Tag[StringLength-Temp_Range.Level] != Branch)
					{
						Temp_Range.Mismatch_Char=Temp_Range.Mismatch_Char | (Branch<<Temp_Range.Mismatches*2);
						Temp_Range.Mismatch_Pos[Temp_Range.Mismatches]=Offset+StringLength-Temp_Range.Level;
						Temp_Range.Mismatches++;
					}

					if (Temp_Range.Mismatches<=Count)//we are guaranteed a valid SA range, check only for mismatches
					{
						if(Temp_Range.Level== StringLength)
						{
							if(Temp_Range.Mismatches==Count)
							{
								Out.hits.push_back(Temp_Range);
								Out.hitCount += Temp_Range.Skip ? 1 : Temp_Range.End-Temp_Range.Start+1;
								if(MaxHits && Out.hitCount>=MaxHits) return;
							}
							else continue;
						}
						else
						{
							Temp_Range.Level++;
							Stack.push_back(Temp_Range);//Push range
						}
					}
				}
			}
		}
	}
	return;
}

// BatMis one-mismatch scan of one read split into lh|rh, starting from range
// (lh's lookup-table range): lh exact, then rh extended with no mismatches (0|0,
// into exactRight). The ranges that extension defers are resumed with one
// mismatch (0|1, into rightMismatch). Finally rh exact with a mismatch in lh
// (1|0, into leftMismatch). Returns the locations of the 0|0 hits (each range
// capped at MaxHits, 0 = no limit), which bound the later passes.
// Bounds (-n 2), when given: if the exact lh search fails, *LH_Fail = d with read[0, d)
// absent; if the exact rh search fails, *RH_Fail = e with read[e, L) absent. Left as
// they are otherwise.
__attribute__((noinline)) unsigned long One_Mismatch_Search(std::string_view lh, std::string_view rh, SARange range, unsigned MaxHits,
                                                            ForwardResults& exactRight, ForwardResults& rightMismatch, ForwardResults& leftMismatch,
                                                            int* LH_Fail = nullptr, int* RH_Fail = nullptr) {
    const int LH = static_cast<int>(lh.size());
    const int RH = static_cast<int>(rh.size());
    unsigned long readHits = 0;
    for (ForwardResults* r : {&exactRight, &rightMismatch, &leftMismatch}) {
        r->hits.clear();
        r->mismatches.clear();
        r->hitCount = 0;
    }

    // 0|0
    const bool Table_Hit = range.Start != 0;  // else the table's first Level - 1 bases do not occur
    if (Search_Right(lh, 0, range, revfmi) && range.Start != 0) {
        range.Level = 1;
        range.Mismatches = 0;
        range.Mismatch_Char = 0;
        Extend_Right(rh, LH, range, 0, revfmi, exactRight, 0);
    } else if (LH_Fail) {
        *LH_Fail = Table_Hit ? range.Level : range.Level - 1;
    }
    for (const SARange& hit : exactRight.hits) {
        unsigned long hitCount = hit.Skip ? 1UL : (hit.End >= hit.Start ? (hit.End - hit.Start + 1UL) : 0UL);
        if (MaxHits > 0 && hitCount > MaxHits) {
            hitCount = MaxHits;
        }
        readHits += hitCount;
    }

    // 0|1, resuming the ranges the exact extension deferred. Their rows
    // have usually left the cache, so fetch them all before resuming.
    for (const SARange& parent : exactRight.mismatches) {
        if (Use_Ref && parent.Skip) {
            __builtin_prefetch(Pac.data() + (parent.Pos + LH) / 4);  // located: only the genome is read
            continue;
        }
        Prefetch_Row(revfmi, parent.Start);
        if (!parent.Skip && parent.End != parent.Start) Prefetch_Row(revfmi, parent.End + 1);
    }
    for (const SARange& parent : exactRight.mismatches) {
        if (MaxHits && readHits + rightMismatch.hitCount >= MaxHits) break;
        Extend_Right(rh, LH, parent, 1, revfmi, rightMismatch,
                     MaxHits ? MaxHits - readHits : 0);
    }

    // 1|0
    if (MaxHits == 0 || readHits + rightMismatch.hitCount < MaxHits) {
        int bc = 0;
        for (int i = 0; i < LOOKUPSIZE; i++) {
            bc |= rh[RH - 1 - i] << (2 * i);
        }
        SARange rightRange{};
        rightRange.Start = Backward_Lookup[bc].Start;
        rightRange.End = Backward_Lookup[bc].End;
        rightRange.Level = LOOKUPSIZE + 1;
        const bool Right_Table_Hit = rightRange.Start != 0;
        Search_Left(rh, LH, rightRange, fwfmi);
        if (!rightRange.Start && RH_Fail) *RH_Fail = LH + RH - (Right_Table_Hit ? rightRange.Level : rightRange.Level - 1);
        if (rightRange.Start) {
            rightRange.Level = 1;
            Extend_Left(lh, 0, rightRange, 1, fwfmi, leftMismatch,
                        MaxHits ? MaxHits - readHits - rightMismatch.hitCount : 0);
        }
    }
    return readHits;
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
        for (int i = 1; i < length && range.Start; i++) Get_SARange(read[i], range, revfmi);
        range.Level = length + 1;
        return range;
    }
    int c = 0;
    for (int i = 0; i < LOOKUPSIZE; i++) c |= read[i] << (2 * i);
    range.Start = Forward_Lookup[c].Start;
    range.End = Forward_Lookup[c].End;
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
        for (int i = 1; i < Len && R.Start; i++) Get_SARange(Seg[Len - 1 - i], R, fwfmi);
        R.Level = Len + 1;
        return R;
    }
    int bc = 0;
    for (int i = 0; i < LOOKUPSIZE; i++) bc |= Seg[Len - 1 - i] << (2 * i);
    R.Start = Backward_Lookup[bc].Start;
    R.End = Backward_Lookup[bc].End;
    R.Level = LOOKUPSIZE + 1;
    return R;
}

// Exact search of Seg (at Offset in the read): Search_Right in revfmi or Search_Left in
// fwfmi from Segment_Start. Returns whether it occurs; R is its range.
bool Search_Segment(std::string_view Seg, int Offset, bool Right, SARange& R) {
    R = Segment_Start(Seg, Right);
    if (R.Level > static_cast<int>(Seg.size())) return R.Start != 0;
    if (!R.Start) return false;
    if (Right) return Search_Right(Seg, Offset, R, revfmi) && R.Start != 0;
    Search_Left(Seg, Offset, R, fwfmi);
    return R.Start != 0;
}

// Result lists of a two-mismatch search, in output order after the first three
// (One_Mismatch_Search's exactRight, rightMismatch, leftMismatch). 1|1 is split by the
// read's thirds T1 = [0, a), T2 = [a, b), T3 = [b, L) (Thirds); p and q are its
// mismatches in LH and RH.
enum Two_Mismatch_List : unsigned {
    TM_EXACT, TM_RIGHT1, TM_LEFT1,  // 0|0, 0|1, 1|0
    TM_RIGHT2, TM_LEFT2,            // 0|2 (resumed 0|1 rejects), 2|0 (resumed 1|0 rejects)
    TM_CROSS_A,                     // 1|1, p in [a, LH): T1 exact, extended right over [a, L) (revfmi)
    TM_CROSS_B,                     // 1|1, p in T1, q in [LH, b): T3 exact, extended left over [0, b) (fwfmi)
    TM_SEED_C,                      // 1|1, p in T1, q in T3: T2 exact, extended right over T3 (seeds)
    TM_CROSS_C,                     // ... each seed in fwfmi, extended left over T1
    TM_LISTS
};

// The thirds' boundaries a < LH < b of a read of L bases.
inline void Thirds(int L, int& a, int& b) {
    a = L / 3;
    b = 2 * L / 3;
}

// Read segments used by the searches: the read, its halves, its thirds, [a, L) and [0, b).
enum Read_Segment : unsigned char { SEG_ALL, SEG_LH, SEG_RH, SEG_T1, SEG_T2, SEG_T3, SEG_A_END, SEG_START_B };

inline std::string_view Read_Segment_Of(std::string_view All, unsigned Seg, int& Offset) {
    const int L = static_cast<int>(All.size()), LH = L / 2;
    int a, b;
    Thirds(L, a, b);
    static const int Lo[8] = {0, 0, -1, 0, -2, -3, -2, 0};  // -1 LH, -2 a, -3 b
    const int Start = Lo[Seg] == -1 ? LH : Lo[Seg] == -2 ? a : Lo[Seg] == -3 ? b : 0;
    const int End = Seg == SEG_LH ? LH : Seg == SEG_T1 ? a : Seg == SEG_T2 || Seg == SEG_START_B ? b : L;
    Offset = Start;
    return All.substr(Start, End - Start);
}

// Mismatch windows for the 1|1 extensions: read positions [Lo, Hi) that must hold
// exactly Need mismatches. Positions outside every window must match.
struct Mm_Windows {
    int N = 0;
    int Lo[3], Hi[3], Need[3];
};
enum Mm_Window_Set : unsigned char { WIN_NONE, WIN_A, WIN_B, WIN_C3, WIN_C1 };

// The windows of Set for a read of L bases, narrowed by the exact half searches' failure
// depths (One_Mismatch_Search): read[0, LH_Fail) and read[RH_Fail, L) do not occur, so
// a 1|1 hit has p < LH_Fail and q >= RH_Fail (L and 0: no bound).
inline Mm_Windows Window_Set(unsigned Set, int L, int LH_Fail, int RH_Fail) {
    int a, b;
    Thirds(L, a, b);
    const int LH = L / 2;
    Mm_Windows W;
    auto Add = [&](int Lo, int Hi) {
        if (Hi <= LH) Hi = std::min(Hi, LH_Fail);
        else Lo = std::max(Lo, RH_Fail);
        W.Lo[W.N] = Lo;
        W.Hi[W.N] = Hi;
        W.Need[W.N] = 1;
        W.N++;
    };
    switch (Set) {
        case WIN_A: Add(a, LH); Add(LH, L); break;  // p in [a, LH), q in RH
        case WIN_B: Add(0, a); Add(LH, b); break;   // p in T1, q in [LH, b)
        case WIN_C3: Add(b, L); break;              // q in T3
        case WIN_C1: Add(0, a); break;              // p in T1
        default: break;
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

// Hits the read has so far: its exact count plus the mismatch lists, leaving out Skip.
inline unsigned long Two_Mismatch_Used(const ForwardResults* Res, unsigned long ReadHits, int Skip = -1) {
    unsigned long Used = ReadHits;
    for (unsigned k : {TM_RIGHT1, TM_LEFT1, TM_RIGHT2, TM_LEFT2, TM_CROSS_A, TM_CROSS_B, TM_CROSS_C})
        if (static_cast<int>(k) != Skip) Used += Res[k].hitCount;
    return Used;
}

// Batched search (--batch B > 1). Each step of a search waits for one cache miss,
// and one read's steps depend on each other, so B reads ("lanes") are searched at
// once: each round prefetches what every lane's next step reads, then steps every
// lane. A lane runs its read's pass program (Mm_Program: exact search, then extend),
// one task at a time: a search of a segment (Search_Right / Search_Left and their
// singleton walks) or one depth-first extension (Extend_Right / Extend_Left, their
// singleton walks and the located path). The step functions are those loops cut at
// each memory access, in the same order, so each read's hits and deferrals come out
// exactly as Search_Right (-n 0) and One_Mismatch_Search (-n 1) make them.
// Branch_Detect's globals are written and read within one step.
// See benchmark/batched_reads.md.
constexpr unsigned MAX_BATCH = 64;

enum Mm_Sub : unsigned char { MM_POP, MM_WIDE, MM_SINGLE, MM_LOCATE, MM_LOCATED, MM_SWITCH };
// A lane's task: an exact search, a depth-first extension (Count mismatches, deferring
// at the threshold), a windowed extension (Mm_Windows, no deferral) or an index switch.
enum Mm_Task : unsigned char { TASK_SEARCH, TASK_EXTEND, TASK_WINDOW, TASK_SWITCH };

// A read's pass program, on the read's segments (Read_Segment).
enum Mm_Op : unsigned char {
    OP_SEARCH,          // exact search of Seg (Right: forward table + revfmi; else backward table + fwfmi);
                        // with Win, skipped when those windows cannot hold (Windows_Possible)
    OP_EXTEND,          // extend the last search's range over Seg with Count mismatches into list Out
    OP_COUNT_EXACT,     // the read's exact hit count from list From's hits (each capped at MaxHits)
    OP_RESUME,          // extend each range list From deferred, in order, as OP_EXTEND
    OP_CROSS,           // extend each hit of list From (a seed) over Seg, keeping its mismatches; with
                        // Switch, a seed not located first gets the other index's range of Seg2
                        // (with its mismatch applied) by a Start-only walk (Mm_Start_Switch)
    OP_STOP_IF_MAXHITS, // stop once the hits so far reach MaxHits
    OP_END
};
struct Mm_Step {
    Mm_Op Op;
    bool Right;
    unsigned char Seg, Count, Out, From;
    bool Limited;  // extensions stop at MaxHits (less the hits already counted)
    bool Switch;
    unsigned char Win, Seg2;  // Win: windowed extension, or a search's case (Mm_Window_Set); Seg2: OP_CROSS's switch segment
};
constexpr unsigned MM_PROGRAM_MAX = 24;  // longest program (-n 2: 21 steps), for the stats
constexpr unsigned MM_LISTS = TM_LISTS;  // the lists of Two_Mismatch_Search (-n 1 uses the first three)
// The lists whose hits are reported, and count towards MaxHits (not the 1|1 seeds).
constexpr unsigned MM_HIT_LISTS = 1u << TM_RIGHT1 | 1u << TM_LEFT1 | 1u << TM_RIGHT2 | 1u << TM_LEFT2 |
                                   1u << TM_CROSS_A | 1u << TM_CROSS_B | 1u << TM_CROSS_C;

// -n 0: Search_Right over the read.
constexpr Mm_Step Program_Exact[] = {{OP_SEARCH, true, 0, 0, 0, 0, false, false},
                                     {OP_END, false, 0, 0, 0, 0, false, false}};
// -n 1, as One_Mismatch_Search: 0|0, 0|1 from the deferred 0|0 ranges, then 1|0.
// Lists: 0 exactRight, 1 rightMismatch, 2 leftMismatch.
#define MM_ONE_MISMATCH_PASSES                                   \
    {OP_SEARCH, true, 1, 0, 0, 0, false, false},                 \
    {OP_EXTEND, true, 2, 0, TM_EXACT, 0, false, false},          \
    {OP_COUNT_EXACT, false, 0, 0, 0, TM_EXACT, false, false},    \
    {OP_RESUME, true, 2, 1, TM_RIGHT1, TM_EXACT, true, false},   \
    {OP_STOP_IF_MAXHITS, false, 0, 0, 0, 0, false, false},       \
    {OP_SEARCH, false, 2, 0, 0, 0, false, false},                \
    {OP_EXTEND, false, 1, 1, TM_LEFT1, 0, true, false}
constexpr Mm_Step Program_One_Mismatch[] = {MM_ONE_MISMATCH_PASSES, {OP_END, false, 0, 0, 0, 0, false, false}};
// -n 2, as Two_Mismatch_Search: the -n 1 passes, 0|2 and 2|0 from their rejects, then
// 1|1 by the thirds: (a) T1 exact, (b) T3 exact, (c) T2 exact, each with its windows,
// narrowed by where the exact LH / RH searches failed; (a) and (b) are skipped when
// their windows cannot hold.
constexpr Mm_Step Program_Two_Mismatch[] = {
    MM_ONE_MISMATCH_PASSES,
    {OP_RESUME, true, SEG_RH, 2, TM_RIGHT2, TM_RIGHT1, true, false},
    {OP_RESUME, false, SEG_LH, 2, TM_LEFT2, TM_LEFT1, true, false},
    {OP_STOP_IF_MAXHITS, false, 0, 0, 0, 0, false, false},
    {OP_SEARCH, true, SEG_T1, 0, 0, 0, false, false, WIN_A},
    {OP_EXTEND, true, SEG_A_END, 0, TM_CROSS_A, 0, true, false, WIN_A},
    {OP_STOP_IF_MAXHITS, false, 0, 0, 0, 0, false, false},
    {OP_SEARCH, false, SEG_T3, 0, 0, 0, false, false, WIN_B},
    {OP_EXTEND, false, SEG_START_B, 0, TM_CROSS_B, 0, true, false, WIN_B},
    {OP_STOP_IF_MAXHITS, false, 0, 0, 0, 0, false, false},
    {OP_SEARCH, true, SEG_T2, 0, 0, 0, false, false},
    {OP_EXTEND, true, SEG_T3, 0, TM_SEED_C, 0, false, false, WIN_C3},
    {OP_CROSS, false, SEG_T1, 0, TM_CROSS_C, TM_SEED_C, true, true, WIN_C1, SEG_A_END},
    {OP_END, false, 0, 0, 0, 0, false, false},
};

#ifdef BATMAN_BATCH_STATS
// Lane occupancy and step counts (build with -DBATMAN_BATCH_STATS; printed with --timing).
struct {
    unsigned long Rounds = 0, Lane_Rounds = 0, Drain_Rounds = 0, Lanes = 0;
    unsigned long Steps[MM_PROGRAM_MAX][4][6] = {};  // [program step][task][sub]
    unsigned long Skipped[MM_PROGRAM_MAX] = {};      // searches skipped by their windows
} Batch_Stats;
#endif

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
    int LH_Fail, RH_Fail; // the exact LH / RH searches' failure bounds (Window_Set)
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
    uint64_t Packed[MAXTAG / 32 + 1];
    char Patched[MAXTAG];    // TASK_SWITCH: the read with the seed's mismatch applied
};

// The read base at Level and its read position, for the task's direction.
inline unsigned Mm_Base(const Mm_Lane& L, int Level) {
    return static_cast<unsigned char>(L.Right ? L.Tag[Level - 1] : L.Tag[L.Tag.size() - Level]);
}
inline int Mm_Read_Pos(const Mm_Lane& L, int Level) {
    return L.Right ? L.Offset + Level - 1 : L.Offset + static_cast<int>(L.Tag.size()) - Level;
}

// Search_Right / Search_Left: the top of their loop.
inline void Mm_Search_Enter(Mm_Lane& L) {
    SARange& T = L.Cur;
    if (T.End != T.Start && !T.Skip) {
        L.Sub = MM_WIDE;
        return;
    }
    L.Sub = MM_SINGLE;
    if (T.Start % SAINTERVAL == 0 && !T.Skip) {
        if (Use_Ref) {
            L.Done = T.Level - 1;
            L.Sub = MM_LOCATE;
            return;
        }
        T.Skip++;
        T.End = T.Start;
    }
}

// Ends a search task; Search_Right reports Found, Search_Left a zero Start on failure.
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
            T.Start = L.Fmi->cumulativeFreq[Base] + BWTOccValue(L.Fmi, T.Start, Base) + 1;
            if (!T.Skip && T.Start % SAINTERVAL == 0) {
                if (Use_Ref) {
                    L.Done = T.Level;
                    L.Sub = MM_LOCATE;
                    return false;
                }
                T.Skip++;
                T.End = T.Start;
            }
            if (T.Level == Len) {
                if (L.Right || !T.Skip) T.End = T.Start;  // as Search_Right / Search_Left_Singleton
                return Mm_Search_End(L, true);
            }
            T.Level++;
            return false;
        }
        case MM_LOCATE:
            if (!Locate(L.Fmi, T, L.Right, L.Offset, Len, L.Done)) {
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

// End of a singleton walk inside an extension: Extend_Right/Left's MaxHits check.
inline bool Mm_Single_Done(Mm_Lane& L) {
    L.Sub = MM_POP;
    return L.TaskMaxHits && L.Out->hitCount >= L.TaskMaxHits;
}

// Extend_Right/Left_Singleton up to their loop (L.Cur popped). Returns true when the
// task has finished; leaves Sub at MM_SINGLE to walk, or MM_LOCATE / MM_LOCATED.
inline bool Mm_Single_Enter(Mm_Lane& L) {
    SARange& T = L.Cur;
    if (T.Start == 0) return Mm_Single_Done(L);
    if (Use_Ref && T.Skip) {
        L.Sub = MM_LOCATED;
        return false;
    }
    L.Sub = MM_SINGLE;
    if (T.Start % SAINTERVAL == 0 && !T.Skip) {
        if (Use_Ref) {
            L.Sub = MM_LOCATE;
            return false;
        }
        T.Skip++;
        T.End = T.Start;
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
            if (Now == Expected) T.Start = L.Fmi->cumulativeFreq[Base] + BWTOccValue(L.Fmi, T.Start, Base) + 1;
            else T.Start = L.Fmi->cumulativeFreq[Now] + BWTOccValue(L.Fmi, T.Start, Now) + 1;
            if (!T.Skip && T.Start % SAINTERVAL == 0) {
                T.Skip++;
                T.End = T.Start;
            }
            if (Expected != Now) {
                T.Mismatch_Char = T.Mismatch_Char | (Now << T.Mismatches * 2);
                T.Mismatch_Pos[T.Mismatches] = Mm_Read_Pos(L, T.Level);
                T.Mismatches++;
            }
            if (Use_Ref && T.Skip) {  // just reached a sampled row
                T.Level++;
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
            if (!Locate(L.Fmi, L.Cur, L.Right, L.Offset, Len, L.Cur.Level - 1)) return Mm_Single_Done(L);
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
    L.Found = false;
    Mm_Search_Enter(L);
}

// Returns false if Extend_Right/Left would return at once (no range).
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

// One step of a windowed extension: Extend_Right/Left's loop and singleton walk, where a
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
            if (Now == Expected) T.Start = L.Fmi->cumulativeFreq[Base] + BWTOccValue(L.Fmi, T.Start, Base) + 1;
            else T.Start = L.Fmi->cumulativeFreq[Now] + BWTOccValue(L.Fmi, T.Start, Now) + 1;
            if (!T.Skip && T.Start % SAINTERVAL == 0) {
                T.Skip++;
                T.End = T.Start;
            }
            if (Now != Expected) Mm_Record(T, Now, x);
            if (!Window_Keep(T, L.W, x, L.Right)) return Mm_Single_Done(L);
            if (Use_Ref && T.Skip) {  // just reached a sampled row
                T.Level++;
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
            if (!Locate(L.Fmi, L.Cur, L.Right, L.Offset, Len, L.Cur.Level - 1)) return Mm_Single_Done(L);
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

// An index switch: the range of Seg of the read, with Seed's (first) mismatch applied, in
// the other index (Right: revfmi, else fwfmi), for a seed found in one index and to be
// extended in the other. The string occurs (it is the seed's), so every range on the way
// is non-empty: Start alone is stepped from the lookup table (S' = C[c] + Occ(c, S) + 1,
// exact in that case) and End = Start + width - 1, the width being the same in both
// indexes. There is no located shortcut on the way: until the last base the range is
// wider than the seed's (shorter strings occur more often), so Start is not yet the
// seed's row, and without End it is not known when it becomes so.
inline void Mm_Start_Switch(Mm_Lane& L, const SARange& Seed, unsigned Seg, bool Right) {
    std::memcpy(L.Patched, L.All.data(), L.All.size());
    L.Patched[Seed.Mismatch_Pos[0]] = static_cast<char>(Seed.Mismatch_Char & 3);
    int Offset;
    const std::string_view S = Read_Segment_Of(L.All, Seg, Offset);
    L.Task = TASK_SWITCH;
    L.Right = Right;
    L.Fmi = Right ? revfmi : fwfmi;
    L.Tag = std::string_view(L.Patched + Offset, S.size());
    L.Offset = Offset;
    L.Width = Seed.Skip ? 1 : Seed.End - Seed.Start + 1;
    L.Cur = Segment_Start(L.Tag, Right);
    L.Cur.Mismatches = Seed.Mismatches;
    L.Cur.Mismatch_Char = Seed.Mismatch_Char;
    std::memcpy(L.Cur.Mismatch_Pos, Seed.Mismatch_Pos, sizeof L.Cur.Mismatch_Pos);
    L.Switching = true;
    L.Found = false;
    L.Sub = MM_SWITCH;
}

// One step of an index switch. Returns true when it has finished (Found, range in Cur).
inline bool Mm_Switch_Step(Mm_Lane& L) {
    SARange& T = L.Cur;
    const int Len = static_cast<int>(L.Tag.size());
    if (!T.Start) return true;
    if (T.Level > Len) {
        T.End = T.Start + L.Width - 1;
        L.Found = true;
        return true;
    }
    const unsigned c = Mm_Base(L, T.Level);
    T.Start = L.Fmi->cumulativeFreq[c] + BWTOccValue(L.Fmi, T.Start, c) + 1;
    T.Level++;
    return false;
}

// A segment of the lane's read and its offset in the read.
inline std::string_view Mm_Segment(const Mm_Lane& L, unsigned Seg, int& Offset) {
    return Read_Segment_Of(L.All, Seg, Offset);
}

// Hits the read has so far: its exact count plus the mismatch lists, leaving out Skip.
inline unsigned long Mm_Used(const Mm_Lane& L, int Skip = -1) {
    unsigned long Used = L.ReadHits;
    if (L.Res)
        for (unsigned k = 1; k < MM_LISTS; k++)
            if ((MM_HIT_LISTS >> k & 1) && static_cast<int>(k) != Skip) Used += L.Res[k].hitCount;
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

// Runs the lane's program until it starts a task or ends (Done_Read). Called when a
// read starts and whenever its task finishes.
void Mm_Run(Mm_Lane& L, unsigned MaxHits) {
    for (;;) {
        const Mm_Step& P = L.Program[L.PC];
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
                const std::string_view Seg = Mm_Segment(L, P.Seg, Offset);
                const SARange R = Segment_Start(Seg, P.Right);
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
            case OP_EXTEND: {
                if (L.InTask) {
                    L.InTask = false;
                    L.PC++;
                    continue;
                }
                if (!L.Found || !L.Cur.Start) {
                    L.PC++;
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

inline void Mm_Prefetch(const Mm_Lane& L) {
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
            __builtin_prefetch(L.Fmi->saValue + T->Start / SAINTERVAL);
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
// -n 0 (Program_Exact): Ranges[i] / Found[i] as Search_Right from Start_Range.
// -n 1 (Program_One_Mismatch): Res[3 * i ...] get read i's exactRight, rightMismatch
// and leftMismatch, ReadHits[i] its exact hit count (One_Mismatch_Search's value).
void Batch_Search(const char* Reads, unsigned Count, int Len, unsigned B, unsigned MaxHits, const Mm_Step* Program,
                  std::vector<ForwardResults>* Res, std::vector<unsigned long>* ReadHits, std::vector<SARange>* Ranges,
                  std::vector<unsigned char>* Found) {
    static Mm_Lane Lanes[MAX_BATCH];
    unsigned Active = 0, Next = 0;
    auto Finish = [&](const Mm_Lane& L) {
        if (ReadHits) (*ReadHits)[L.Index] = L.ReadHits;
        if (Ranges) {
            (*Ranges)[L.Index] = L.Cur;
            (*Found)[L.Index] = L.Found;
        }
    };
    auto Fill = [&](Mm_Lane& L) {
        while (Next < Count) {
            const unsigned I = Next++;
            L.All = std::string_view(Reads + static_cast<std::size_t>(I) * Len, Len);
            L.Index = I;
            L.Res = Res ? &(*Res)[MM_LISTS * static_cast<std::size_t>(I)] : nullptr;
            if (L.Res)
                for (unsigned k = 0; k < MM_LISTS; k++) {
                    L.Res[k].hits.clear();
                    L.Res[k].mismatches.clear();
                    L.Res[k].hitCount = 0;
                }
            L.ReadHits = 0;
            L.Program = Program;
            L.PC = 0;
            L.InTask = L.Done_Read = false;
            L.Resume = 0;
            L.Switching = false;
            L.Found = false;
            L.LH_Fail = Len;
            L.RH_Fail = 0;
            if (Use_Ref) Pack_Read(L.All, L.Packed);
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
    Batch_Search(Reads, Count, Len, B, 0, Program_Exact, nullptr, nullptr, &Ranges, &Found);
}

void Batch_One_Mismatch(const char* Reads, unsigned Count, int Len, unsigned B, unsigned MaxHits,
                        std::vector<ForwardResults>& Res, std::vector<unsigned long>& ReadHits) {
    Batch_Search(Reads, Count, Len, B, MaxHits, Program_One_Mismatch, &Res, &ReadHits, nullptr, nullptr);
}

void Batch_Two_Mismatch(const char* Reads, unsigned Count, int Len, unsigned B, unsigned MaxHits,
                        std::vector<ForwardResults>& Res, std::vector<unsigned long>& ReadHits) {
    Batch_Search(Reads, Count, Len, B, MaxHits, Program_Two_Mismatch, &Res, &ReadHits, nullptr, nullptr);
}

// The per-read path's windowed extensions and index switches: the lane steps above, run
// to the end on a private lane, so both paths share one implementation.
Mm_Lane& Seq_Lane(std::string_view All) {
    static Mm_Lane L;
    L.All = All;
    std::memcpy(L.Packed, Read_Packed, (All.size() + 31) / 32 * sizeof(uint64_t));
    L.Switching = false;
    return L;
}

// Extends Seed (its Mismatches kept) over segment Seg with windows W into Out.
void Extend_Windows(std::string_view All, unsigned Seg, SARange Seed, bool Right, const Mm_Windows& W, ForwardResults& Out,
                    unsigned MaxHits) {
    Mm_Lane& L = Seq_Lane(All);
    int Offset;
    const std::string_view Tag = Read_Segment_Of(All, Seg, Offset);
    Seed.Level = 1;
    Seed.Forbidden_Base = NO_BASE;
    if (!Mm_Start_Window(L, Right, Tag, Offset, Seed, W, &Out, MaxHits)) return;
    while (!Mm_Window_Step(L)) {
    }
}

// Seed's range of segment Seg (its mismatch applied) in the other index (Mm_Start_Switch).
bool Switch_Index(std::string_view All, const SARange& Seed, unsigned Seg, bool Right, SARange& R) {
    Mm_Lane& L = Seq_Lane(All);
    Mm_Start_Switch(L, Seed, Seg, Right);
    while (!Mm_Switch_Step(L)) {
    }
    R = L.Cur;
    return L.Found;
}

// Up to two mismatches, the halves LH | RH as pigeonholes:
//   0|0, 0|1, 1|0  One_Mismatch_Search;
//   0|2            the 0|1 rejects (deferred ranges) extended over RH with Count 2;
//   2|0            the 1|0 rejects extended over LH with Count 2;
//   1|1            by the thirds T1 = [0, a), T2 = [a, b), T3 = [b, L) (a < LH < b), with
//                  mismatch windows so that each hit is found once (p in LH, q in RH),
//                  narrowed by where the exact LH / RH searches failed (Window_Set); a
//                  case whose windows cannot hold is skipped:
//     (a) p in [a, LH):          T1 exact, extended right over [a, L) (revfmi);
//     (b) p in T1, q in [LH, b): T3 exact, extended left over [0, b) (fwfmi);
//     (c) p in T1, q in T3:      T2 exact, extended right over T3 (revfmi), each such seed
//                                switched to fwfmi (unless located) and extended left over T1.
// Stops between stages at MaxHits. Res has TM_LISTS lists; returns the exact hit count.
// With the reference, the read must be packed in Read_Packed (Pack_Read).
__attribute__((noinline)) unsigned long Two_Mismatch_Search(std::string_view All, const SARange& Range, unsigned MaxHits,
                                                            ForwardResults* Res) {
    const int LH = static_cast<int>(All.size()) / 2;
    const std::string_view lh = All.substr(0, LH), rh = All.substr(LH);
    for (unsigned k = TM_RIGHT2; k < TM_LISTS; k++) {
        Res[k].hits.clear();
        Res[k].mismatches.clear();
        Res[k].hitCount = 0;
    }
    int LH_Fail = static_cast<int>(All.size()), RH_Fail = 0;
    const unsigned long ReadHits =
        One_Mismatch_Search(lh, rh, Range, MaxHits, Res[TM_EXACT], Res[TM_RIGHT1], Res[TM_LEFT1], &LH_Fail, &RH_Fail);
    auto Done = [&] { return MaxHits && Two_Mismatch_Used(Res, ReadHits) >= MaxHits; };
    auto Budget = [&](unsigned Out) { return MaxHits ? MaxHits - Two_Mismatch_Used(Res, ReadHits, Out) : 0; };

    // 0|2
    for (const SARange& Parent : Res[TM_RIGHT1].mismatches) {
        if (Done()) return ReadHits;
        Extend_Right(rh, LH, Parent, 2, revfmi, Res[TM_RIGHT2], Budget(TM_RIGHT2));
    }
    // 2|0
    for (const SARange& Parent : Res[TM_LEFT1].mismatches) {
        if (Done()) return ReadHits;
        Extend_Left(lh, 0, Parent, 2, fwfmi, Res[TM_LEFT2], Budget(TM_LEFT2));
    }
    // 1|1
    SARange Seed;
    int Offset;
    auto Exact = [&](unsigned Seg, bool Right) {
        const std::string_view Part = Read_Segment_Of(All, Seg, Offset);  // sets Offset: not in the same call
        if (!Search_Segment(Part, Offset, Right, Seed)) return false;
        Seed.Mismatches = 0;
        Seed.Mismatch_Char = 0;
        return true;
    };
    const int L = static_cast<int>(All.size());
    const Mm_Windows A = Window_Set(WIN_A, L, LH_Fail, RH_Fail), B = Window_Set(WIN_B, L, LH_Fail, RH_Fail);
    if (Done()) return ReadHits;
    if (Windows_Possible(A) && Exact(SEG_T1, true)) Extend_Windows(All, SEG_A_END, Seed, true, A, Res[TM_CROSS_A], Budget(TM_CROSS_A));
    if (Done()) return ReadHits;
    if (Windows_Possible(B) && Exact(SEG_T3, false)) Extend_Windows(All, SEG_START_B, Seed, false, B, Res[TM_CROSS_B], Budget(TM_CROSS_B));
    if (Done()) return ReadHits;
    if (Exact(SEG_T2, true)) Extend_Windows(All, SEG_T3, Seed, true, Window_Set(WIN_C3, L, LH_Fail, RH_Fail), Res[TM_SEED_C], 0);
    const Mm_Windows C1 = Window_Set(WIN_C1, L, LH_Fail, RH_Fail);
    for (const SARange& S : Res[TM_SEED_C].hits) {
        if (Done()) return ReadHits;
        SARange R = S;
        if (!(Use_Ref && S.Skip) && !Switch_Index(All, S, SEG_A_END, false, R)) continue;
        Extend_Windows(All, SEG_T1, R, false, C1, Res[TM_CROSS_C], Budget(TM_CROSS_C));
    }
    return ReadHits;
}

// "\t<dir>\t<start>\t<end>" for an SA range, or "\tP\t<pos>\t<pos>" for a located hit.
std::string Hit_Text(const SARange& hit, const char* dir) {
    if (Use_Ref && hit.Skip) {
        const std::string pos = std::to_string(hit.Pos);
        return "\tP\t" + pos + "\t" + pos;
    }
    return std::string("\t") + dir + "\t" + std::to_string(hit.Start) + "\t" + std::to_string(hit.Skip ? hit.Start : hit.End);
}

bool ReadNextRecord(std::ifstream& in, std::string& id, std::string& seq) {
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }

        if (line[0] == '@') {
            id = line.substr(1);
            if (!std::getline(in, seq)) {
                return false;
            }
            if (!std::getline(in, line)) {
                return false;
            }
            if (!std::getline(in, line)) {
                return false;
            }
            return true;
        }

        if (line[0] == '>') {
            id = line.substr(1);
            seq.clear();
            std::streampos pos;
            while (true) {
                pos = in.tellg();
                if (!std::getline(in, line)) {
                    break;
                }
                if (!line.empty() && (line[0] == '>' || line[0] == '@')) {
                    in.seekg(pos);
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

// One stage's hits in the output: the list, its index (F revfmi, B fwfmi) and mismatches.
struct Hit_Stage {
    const ForwardResults* List;
    const char* Dir;
    int Mismatches;
};

// The output lines of one -n 1 / -n 2 read (main and the batched runner): the exact
// hits (each capped at maxHits), then the mismatch stages' hits in order until maxHits.
// Returns the read's reported hit count; with none, the caller writes "F 0 0 0".
unsigned long Format_Hits(const std::string& id, const ForwardResults& exact, const Hit_Stage* stages, int count,
                          unsigned long readHits, unsigned maxHits, std::string& lines) {
    for (const SARange& hit : exact.hits) {
        unsigned long hitCount = hit.Skip ? 1UL : (hit.End >= hit.Start ? (hit.End - hit.Start + 1UL) : 0UL);
        if (maxHits > 0 && hitCount > maxHits) {
            hitCount = maxHits;
        }
        lines += id + Hit_Text(hit, "F") + "\t" + std::to_string(hitCount) + "\n";
    }
    for (int s = 0; s < count; s++) {
        for (const SARange& hit : stages[s].List->hits) {
            if (maxHits > 0 && readHits >= maxHits) break;
            unsigned long hitCount = hit.Skip ? 1UL : hit.End - hit.Start + 1UL;
            if (maxHits > 0 && readHits + hitCount > maxHits) {
                hitCount = maxHits - readHits;
            }
            readHits += hitCount;
            lines += id + Hit_Text(hit, stages[s].Dir) + "\t" + std::to_string(hitCount) + "\t" +
                     std::to_string(stages[s].Mismatches) + "\n";
        }
    }
    return readHits;
}

unsigned long Format_One_Mismatch(const std::string& id, const ForwardResults& exactRight,
                                  const ForwardResults& rightMismatch, const ForwardResults& leftMismatch,
                                  unsigned long readHits, unsigned maxHits, std::string& lines) {
    const Hit_Stage stages[] = {{&rightMismatch, "F", 1}, {&leftMismatch, "B", 1}};
    return Format_Hits(id, exactRight, stages, 2, readHits, maxHits, lines);
}

// Two_Mismatch_Search's lists in output order: 0|1, 1|0, 0|2, 2|0, then 1|1 (a), (b), (c).
unsigned long Format_Two_Mismatch(const std::string& id, const ForwardResults* Res, unsigned long readHits,
                                  unsigned maxHits, std::string& lines) {
    const Hit_Stage stages[] = {{&Res[TM_RIGHT1], "F", 1},  {&Res[TM_LEFT1], "B", 1},   {&Res[TM_RIGHT2], "F", 2},
                                {&Res[TM_LEFT2], "B", 2},   {&Res[TM_CROSS_A], "F", 2}, {&Res[TM_CROSS_B], "B", 2},
                                {&Res[TM_CROSS_C], "B", 2}};
    return Format_Hits(id, Res[TM_EXACT], stages, 7, readHits, maxHits, lines);
}

// --batch B > 1: reads the records in blocks, searches each block with Batch_Exact
// (-n 0, or reads too short to split into halves) or Batch_One_Mismatch (-n 1) and
// writes its lines in input order, as the per-read loop in main does. Encoding the
// reads counts as input here. Returns false on an unsupported read length.
bool Run_Batched(std::ifstream& in, std::ostream& out, unsigned B, int maxMismatches, unsigned maxHits, bool quiet,
                 uint64_t& records, uint64_t& reportedHits, double& tableBuildSeconds,
                 double& inputParseSeconds, double& searchSeconds, double& outputSeconds) {
    using Clock = std::chrono::steady_clock;
#ifdef BATMAN_BATCH_BLOCK
    constexpr unsigned BLOCK = BATMAN_BATCH_BLOCK;
#else
    constexpr unsigned BLOCK = 4096;
#endif
    std::vector<std::string> ids(BLOCK);
    std::vector<unsigned char> searched(BLOCK);  // valid read of the right length
    std::vector<unsigned> slot(BLOCK);            // its place among the searched reads
    std::vector<SARange> ranges(BLOCK);
    std::vector<unsigned char> found(BLOCK);
    std::vector<ForwardResults> results(MM_LISTS * BLOCK);  // -n 1 / -n 2: each read's result lists
    std::vector<unsigned long> readHits(BLOCK);
    std::string reads;  // codes of the searched reads, back to back
    std::string seq, lines;
    bool initialized = false, more = true;
    while (more) {
        unsigned n = 0, count = 0;
        double tableSeconds = 0.0;
        reads.clear();
        const auto parseStart = Clock::now();
        while (n < BLOCK && (more = ReadNextRecord(in, ids[n], seq))) {
            ++records;
            if (!initialized) {
                STRINGLENGTH = static_cast<int>(seq.size());
                if (STRINGLENGTH <= 0 || STRINGLENGTH >= MAXTAG) {
                    std::fprintf(stderr, "Unsupported read length: %d\n", STRINGLENGTH);
                    return false;
                }
                const auto tableStart = Clock::now();
                Build_Tables(revfmi);
                tableSeconds = std::chrono::duration<double>(Clock::now() - tableStart).count();
                initialized = true;
            }
            bool valid = static_cast<int>(seq.size()) == STRINGLENGTH;
            const std::size_t at = reads.size();
            if (valid) {
                reads.resize(at + STRINGLENGTH);
                for (int i = 0; i < STRINGLENGTH && valid; i++) {
                    const int code = EncodeBase(seq[i]);
                    valid = code >= 0;
                    reads[at + i] = static_cast<char>(code);
                }
                if (!valid) reads.resize(at);
            }
            searched[n] = valid;
            slot[n] = count;
            count += valid;
            n++;
        }
        tableBuildSeconds += tableSeconds;
        inputParseSeconds += std::chrono::duration<double>(Clock::now() - parseStart).count() - tableSeconds;

        const bool oneMismatch = maxMismatches >= 1 && STRINGLENGTH / 2 > LOOKUPSIZE;
        const auto searchStart = Clock::now();
        if (oneMismatch && maxMismatches == 2) Batch_Two_Mismatch(reads.data(), count, STRINGLENGTH, B, maxHits, results, readHits);
        else if (oneMismatch) Batch_One_Mismatch(reads.data(), count, STRINGLENGTH, B, maxHits, results, readHits);
        else Batch_Exact(reads.data(), count, STRINGLENGTH, B, ranges, found);
        searchSeconds += std::chrono::duration<double>(Clock::now() - searchStart).count();

        const auto outputStart = Clock::now();
        lines.clear();
        for (unsigned j = 0; j < n && oneMismatch; j++) {
            unsigned long hits = 0;
            if (searched[j]) {
                const ForwardResults* r = &results[MM_LISTS * static_cast<std::size_t>(slot[j])];
                const std::size_t at = lines.size();
                hits = maxMismatches == 2 ? Format_Two_Mismatch(ids[j], r, readHits[slot[j]], maxHits, lines)
                                          : Format_One_Mismatch(ids[j], r[0], r[1], r[2], readHits[slot[j]], maxHits, lines);
                if (quiet) lines.resize(at);
            }
            reportedHits += hits;
            if (!hits && !quiet) lines += ids[j] + "\tF\t0\t0\t0\n";
        }
        for (unsigned j = 0; j < n && !oneMismatch; j++) {
            const SARange& range = ranges[slot[j]];
            if (!searched[j] || !found[slot[j]] || range.Start == 0) {
                if (!quiet) lines += ids[j] + "\tF\t0\t0\t0\n";
                continue;
            }
            unsigned long hitCount = range.Skip ? 1UL : (range.End >= range.Start ? (range.End - range.Start + 1UL) : 0UL);
            if (maxHits > 0 && hitCount > maxHits) hitCount = maxHits;
            reportedHits += hitCount;
            if (!quiet) lines += ids[j] + Hit_Text(range, "F") + "\t" + std::to_string(hitCount) + "\n";
        }
        if (!quiet) out << lines;
        outputSeconds += std::chrono::duration<double>(Clock::now() - outputStart).count();
    }
    return true;
}

void Usage(const char* prog) {
    std::fprintf(stderr,
    "Usage: %s -q <reads.fa/fq> -g <genome_prefix> [-n <0|1|2>] [-o <hits.txt>] [-m <maxhits>] [--no-ref] [--lookup <k>] [--batch <B>] [--quiet|-Q] [--timing]\n"
        "  -n  maximum mismatches, 0-2 (default 0); 1 and 2 also need the forward index\n"
        "      <genome_prefix>.bwt/.fmv\n"
        "  -m  maximum hits per read (default 1, 0 = all)\n"
        "  --no-ref  do not load the SA samples (.sa) and packed genome (.pac); single-row\n"
        "            paths then walk the index to the end (less memory, slower)\n"
        "  --lookup  depth k of the start-of-search tables, 3-14 (default 6): 8*4^k bytes\n"
        "            per index; deeper is faster once there are more than about 4^k/3 reads\n"
        "            (see benchmark/lookup_table.md)\n"
        "  --batch   search B reads at a time, 1-64 (default 1), so that their memory\n"
        "            accesses overlap (see benchmark/batched_reads.md)\n"
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
    int batchSize = 1;

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
            if (std::strcmp(argv[i], "0") == 0) {
                maxMismatches = 0;
            } else if (std::strcmp(argv[i], "1") == 0) {
                maxMismatches = 1;
            } else if (std::strcmp(argv[i], "2") == 0) {
                maxMismatches = 2;
            } else {
                Usage(argv[0]);
                return 1;
            }
        } else if (std::strcmp(argv[i], "--quiet") == 0 || std::strcmp(argv[i], "-Q") == 0) {
            quiet = true;
        } else if (std::strcmp(argv[i], "--timing") == 0) {
            timing = true;
        } else if (std::strcmp(argv[i], "--batch") == 0 && i + 1 < argc) {
            batchSize = std::atoi(argv[++i]);
            if (batchSize < 1 || batchSize > static_cast<int>(MAX_BATCH)) {
                std::fprintf(stderr, "--batch must be 1-%u\n", MAX_BATCH);
                return 1;
            }
        } else if (std::strcmp(argv[i], "--lookup") == 0 && i + 1 < argc) {
            LOOKUP_DEPTH = std::atoi(argv[++i]);
            if (LOOKUP_DEPTH < 3 || LOOKUP_DEPTH > 14) {
                std::fprintf(stderr, "--lookup must be 3-14\n");
                return 1;
            }
        } else if (std::strcmp(argv[i], "--no-ref") == 0) {
            useRef = false;
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

    using Clock = std::chrono::steady_clock;
    const auto totalStart = Clock::now();
    double indexLoadSeconds = 0.0;
    double tableBuildSeconds = 0.0;
    double inputParseSeconds = 0.0;
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

    const auto loadStart = Clock::now();
    revfmi = initFMI(revBwt.c_str(), revOcc.c_str(), useRef ? revSa.c_str() : nullptr);
    if (!revfmi) {
        std::fprintf(stderr, "Failed to load reverse index: %s / %s\n", revBwt.c_str(), revOcc.c_str());
        return 1;
    }
    if (maxMismatches > 0) {
        fwfmi = BWTLoad(mmPool, fwBwt.c_str(), fwOcc.c_str(), useRef ? fwSa.c_str() : nullptr, nullptr, nullptr, nullptr);
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
    Text_Length = revfmi->textLength;
    indexLoadSeconds = std::chrono::duration<double>(Clock::now() - loadStart).count();

    std::ifstream in(queryFile.c_str());
    if (!in) {
        std::fprintf(stderr, "Cannot open query file: %s\n", queryFile.c_str());
        return 1;
    }

    std::ostream* out = &std::cout;
    std::ofstream outFile;
    if (!outputFile.empty()) {
        outFile.open(outputFile.c_str());
        if (!outFile) {
            std::fprintf(stderr, "Cannot open output file: %s\n", outputFile.c_str());
            return 1;
        }
        out = &outFile;
    }

    std::string id;
    std::string seq;
    std::string read;  // base codes 0-3
    bool initialized = false;
    uint64_t records = 0;
    uint64_t reportedHits = 0;
    ForwardResults exactRight;     // 0|0 range and deferred parents for 0|1, reverse index
    ForwardResults rightMismatch;  // 0|1 ranges, reverse index
    ForwardResults leftMismatch;   // 1|0 ranges, forward index
    ForwardResults twoMismatch[TM_LISTS];  // -n 2: Two_Mismatch_Search's lists

    const bool batched = batchSize > 1;
    if (batched && !Run_Batched(in, *out, batchSize, maxMismatches, maxHits, quiet, records, reportedHits,
                                tableBuildSeconds, inputParseSeconds, searchSeconds, outputSeconds)) {
        return 1;
    }
    while (!batched) {
        const auto parseStart = Clock::now();
        const bool haveRecord = ReadNextRecord(in, id, seq);
        inputParseSeconds += std::chrono::duration<double>(Clock::now() - parseStart).count();
        if (!haveRecord) break;
        ++records;

        if (!initialized) {
            STRINGLENGTH = static_cast<int>(seq.size());
            if (STRINGLENGTH <= 0 || STRINGLENGTH >= MAXTAG) {
                std::fprintf(stderr, "Unsupported read length: %d\n", STRINGLENGTH);
                return 1;
            }
            const auto tableStart = Clock::now();
            read.resize(STRINGLENGTH);
            Build_Tables(revfmi);
            tableBuildSeconds += std::chrono::duration<double>(Clock::now() - tableStart).count();
            initialized = true;
        }

        if (static_cast<int>(seq.size()) != STRINGLENGTH) {
            if (!quiet) {
                const auto outputStart = Clock::now();
                *out << id << "\tF\t0\t0\t0\n";
                outputSeconds += std::chrono::duration<double>(Clock::now() - outputStart).count();
            }
            continue;
        }

        const auto searchStart = Clock::now();
        bool valid = true;
        for (int i = 0; i < STRINGLENGTH; i++) {
            int code = EncodeBase(seq[i]);
            if (code < 0) {
                valid = false;
                break;
            }
            read[i] = static_cast<char>(code);
        }

        if (!valid) {
            searchSeconds += std::chrono::duration<double>(Clock::now() - searchStart).count();
            if (!quiet) {
                const auto outputStart = Clock::now();
                *out << id << "\tF\t0\t0\t0\n";
                outputSeconds += std::chrono::duration<double>(Clock::now() - outputStart).count();
            }
            continue;
        }

        const std::string_view all(read);
        if (Use_Ref) Pack_Read(all);
        SARange range = Start_Range(all);

        const int LH = STRINGLENGTH / 2;
        const std::string_view lh = all.substr(0, LH);
        const std::string_view rh = all.substr(LH);
        if (maxMismatches >= 1 && LH > LOOKUPSIZE) {
            std::string lines;
            const unsigned long readHits =
                maxMismatches == 2
                    ? Format_Two_Mismatch(id, twoMismatch, Two_Mismatch_Search(all, range, maxHits, twoMismatch), maxHits,
                                          lines)
                    : Format_One_Mismatch(id, exactRight, rightMismatch, leftMismatch,
                                          One_Mismatch_Search(lh, rh, range, maxHits, exactRight, rightMismatch, leftMismatch),
                                          maxHits, lines);
            reportedHits += readHits;
            searchSeconds += std::chrono::duration<double>(Clock::now() - searchStart).count();

            if (!quiet) {
                const auto outputStart = Clock::now();
                if (readHits) {
                    *out << lines;
                } else {
                    *out << id << "\tF\t0\t0\t0\n";
                }
                outputSeconds += std::chrono::duration<double>(Clock::now() - outputStart).count();
            }
            continue;
        }

        bool found = range.Level > STRINGLENGTH ? range.Start != 0 : Search_Right(all, 0, range, revfmi);
        if (!found || range.Start == 0) {
            searchSeconds += std::chrono::duration<double>(Clock::now() - searchStart).count();
            if (!quiet) {
                const auto outputStart = Clock::now();
                *out << id << "\tF\t0\t0\t0\n";
                outputSeconds += std::chrono::duration<double>(Clock::now() - outputStart).count();
            }
            continue;
        }

        unsigned long hitCount = range.Skip ? 1UL : (range.End >= range.Start ? (range.End - range.Start + 1UL) : 0UL);
        if (maxHits > 0 && hitCount > maxHits) {
            hitCount = maxHits;
        }
        reportedHits += hitCount;
        searchSeconds += std::chrono::duration<double>(Clock::now() - searchStart).count();

        if (!quiet) {
            const auto outputStart = Clock::now();
            *out << id << Hit_Text(range, "F") << "\t" << hitCount << "\n";
            outputSeconds += std::chrono::duration<double>(Clock::now() - outputStart).count();
        }
    }

    if (timing) {
        const double totalSeconds = std::chrono::duration<double>(Clock::now() - totalStart).count();
        std::cerr << std::fixed << std::setprecision(6)
                  << "batman_zero_min timing: index_load=" << indexLoadSeconds
                  << "s table_build=" << tableBuildSeconds
                  << "s input_parse=" << inputParseSeconds
                  << "s search=" << searchSeconds
                  << "s output_format=" << outputSeconds
                  << "s total=" << totalSeconds
                  << "s records=" << records
                  << " reported_hits=" << reportedHits << '\n';
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
#endif
    }


    return 0;
}
