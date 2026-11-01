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
void Record_Mismatches(SARange& Tag, bool Right, int A, int B) {
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
            for (uint64_t Bits = Mismatch_Bits(G, W, std::max(A - W * 32, 0), std::min(B - W * 32, 32)); Bits; ) {
                const int K = (__builtin_clzll(Bits) - 1) / 2;	// lowest read position first
                Add(W * 32 + K, (G >> (62 - 2 * K)) & 3);
                Bits &= ~(1ULL << (62 - 2 * K));
            }
        }
    } else {
        for (int W = (B - 1) / 32; W >= A / 32; W--) {
            const uint64_t G = Genome_Word(Tag.Pos + W * 32);
            for (uint64_t Bits = Mismatch_Bits(G, W, std::max(A - W * 32, 0), std::min(B - W * 32, 32)); Bits; Bits &= Bits - 1) {
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

// Search_Right/Search_Left at a sampled row: matches if the rest of the segment does.
bool Search_Located(const ActiveBWT* fmi, SARange& Tag, bool Right, int Offset, int Len, int Done) {
    return Locate(fmi, Tag, Right, Offset, Len, Done) && Rest_Mismatches(Tag, Right, Offset, Len, Done) == 0;
}

// Extend_Right/Extend_Left for a located Tag, whose first Level-1 segment bases are
// matched with Mismatches mismatches. As the walk would: a hit at exactly Count
// mismatches (with the rest of its mismatches recorded); above Count, deferred
// unchanged, to be compared again when resumed with a larger Count.
void Extend_Located(const SARange& Tag, bool Right, int Offset, int Len, int Count, ForwardResults& Out) {
    const int Done = Tag.Level - 1;
    const unsigned D = Tag.Mismatches + Rest_Mismatches(Tag, Right, Offset, Len, Done);
    if (D == static_cast<unsigned>(Count)) {
        SARange Hit = Tag;
        if (Right) Record_Mismatches(Hit, true, Offset + Done, Offset + Len);
        else Record_Mismatches(Hit, false, Offset, Offset + Len - Done);
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
			if (!Search_Located(fmi,Tag,false,Offset,StringLength,Level-1)) Tag.Start=0;
			return;
		}
		Tag.Skip++;
		Tag.End=Tag.Start;
	}
	for(;;)
	{
		if (Tag.Start==fmi->inverseSa0) {Tag.Start=0;return;}//$ row: no base to extend with
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
			Tag.Start=0;//Tag.End=0;
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
__attribute__((noinline)) unsigned long One_Mismatch_Search(std::string_view lh, std::string_view rh, SARange range, unsigned MaxHits,
                                                            ForwardResults& exactRight, ForwardResults& rightMismatch, ForwardResults& leftMismatch) {
    const int LH = static_cast<int>(lh.size());
    const int RH = static_cast<int>(rh.size());
    unsigned long readHits = 0;
    for (ForwardResults* r : {&exactRight, &rightMismatch, &leftMismatch}) {
        r->hits.clear();
        r->mismatches.clear();
        r->hitCount = 0;
    }

    // 0|0
    if (Search_Right(lh, 0, range, revfmi) && range.Start != 0) {
        range.Level = 1;
        range.Mismatches = 0;
        range.Mismatch_Char = 0;
        Extend_Right(rh, LH, range, 0, revfmi, exactRight, 0);
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
        Search_Left(rh, LH, rightRange, fwfmi);
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

// Batched exact search (-n 0 with --batch B > 1). Each step of a search waits for
// one cache miss, and one read's steps depend on each other, so B reads ("lanes")
// are searched at once: each round prefetches what every lane's next step reads,
// then steps every lane. Each lane follows Search_Right and Search_Right_Singleton
// exactly (same ranges, Skip and located positions), so the output is unchanged.
// See benchmark/batched_reads.md.
constexpr unsigned MAX_BATCH = 64;

enum Lane_State : unsigned char { LANE_WIDE, LANE_SINGLE, LANE_LOCATE, LANE_COMPARE, LANE_DONE };

struct Exact_Lane {
    SARange R;
    std::string_view Read;
    unsigned Index;        // the read's place in the block
    int Done;              // LOCATE/COMPARE: bases matched when the sampled row was reached
    Lane_State State;
    uint64_t Packed[MAXTAG / 32 + 1];
};

// Search_Right_Singleton's check on entry and after each step: at a sampled row,
// locate (with the reference) or mark the range as a single row.
inline void Lane_Sampled(Exact_Lane& L, int Done) {
    if (L.R.Skip || L.R.Start % SAINTERVAL != 0) return;
    if (Use_Ref) {
        L.Done = Done;
        L.State = LANE_LOCATE;
        return;
    }
    L.R.Skip++;
    L.R.End = L.R.Start;
}

inline void Lane_Prefetch(const Exact_Lane& L) {
    switch (L.State) {
        case LANE_WIDE:
            Prefetch_Row(revfmi, L.R.Start);
            Prefetch_Row(revfmi, L.R.End + 1);
            break;
        case LANE_SINGLE:
            Prefetch_Row(revfmi, L.R.Start);
            break;
        case LANE_LOCATE:
            __builtin_prefetch(revfmi->saValue + L.R.Start / SAINTERVAL);
            break;
        case LANE_COMPARE:
            __builtin_prefetch(Pac.data() + (L.R.Pos + L.Done) / 4);
            __builtin_prefetch(Pac.data() + (L.R.Pos + L.Read.size()) / 4);
            break;
        default:
            break;
    }
}

// One step of a lane. Returns true when the lane has finished (LANE_DONE), with
// Found as Search_Right would return.
inline bool Lane_Step(Exact_Lane& L, bool& Found) {
    const int Len = static_cast<int>(L.Read.size());
    switch (L.State) {
        case LANE_WIDE:
            Get_SARange(L.Read[L.R.Level - 1], L.R, revfmi);
            if (!L.R.Start) break;
            if (L.R.Level == Len) {
                Found = true;
                return true;
            }
            L.R.Level++;
            if (L.R.End == L.R.Start) {
                L.State = LANE_SINGLE;
                Lane_Sampled(L, L.R.Level - 1);
            }
            return false;
        case LANE_SINGLE: {
            if (L.R.Start == revfmi->inverseSa0) break;  // $ row
            const unsigned long Index = L.R.Start - (L.R.Start > revfmi->inverseSa0);
            const unsigned Now = revfmi->bwtCode[Index / 16] << ((Index % 16) * 2) >> (BITS_IN_WORD - 2);
            const unsigned Expected = static_cast<unsigned char>(L.Read[L.R.Level - 1]);
            const unsigned Base = Read_Base(Expected);
            if (Expected != Now) break;
            L.R.Start = revfmi->cumulativeFreq[Base] + BWTOccValue(revfmi, L.R.Start, Base) + 1;
            Lane_Sampled(L, L.R.Level);
            if (L.State == LANE_LOCATE) return false;
            if (L.R.Level == Len) {
                L.R.End = L.R.Start;
                Found = true;
                return true;
            }
            L.R.Level++;
            return false;
        }
        case LANE_LOCATE:
            if (!Locate(revfmi, L.R, true, 0, Len, L.Done)) break;
            L.State = LANE_COMPARE;
            return false;
        case LANE_COMPARE:
            if (Rest_Mismatches(L.R, true, 0, Len, L.Done, L.Packed)) break;
            L.R.Level = Len;
            Found = true;
            return true;
        default:
            break;
    }
    Found = false;
    return true;
}

// Searches Count reads of Len bases (codes, back to back) and sets Ranges[i] and
// Found[i] as Search_Right from Start_Range would, B reads at a time.
void Batch_Exact(const char* Reads, unsigned Count, int Len, unsigned B, std::vector<SARange>& Ranges,
                 std::vector<unsigned char>& Found) {
    static Exact_Lane Lanes[MAX_BATCH];
    unsigned Active = 0, Next = 0;
    // Starts the next read in lane L; reads that need no search finish at once.
    auto Fill = [&](Exact_Lane& L) {
        while (Next < Count) {
            const unsigned I = Next++;
            L.Read = std::string_view(Reads + static_cast<std::size_t>(I) * Len, Len);
            L.Index = I;
            L.R = Start_Range(L.Read);
            if (L.R.Level > Len || !L.R.Start) {  // read no longer than the lookup key, or no match
                Ranges[I] = L.R;
                Found[I] = L.R.Start != 0;
                continue;
            }
            if (Use_Ref) Pack_Read(L.Read, L.Packed);
            L.State = LANE_WIDE;
            if (L.R.End == L.R.Start || L.R.Skip) {
                L.State = LANE_SINGLE;
                Lane_Sampled(L, L.R.Level - 1);
            }
            return true;
        }
        return false;
    };
    for (unsigned i = 0; i < B && Fill(Lanes[Active]); i++) Active++;
    while (Active) {
        for (unsigned i = 0; i < Active; i++) Lane_Prefetch(Lanes[i]);
        for (unsigned i = 0; i < Active;) {
            Exact_Lane& L = Lanes[i];
            bool F;
            if (!Lane_Step(L, F)) {
                i++;
                continue;
            }
            Ranges[L.Index] = L.R;
            Found[L.Index] = F;
            if (!Fill(L)) {
                L = Lanes[--Active];  // no reads left: close the gap
                continue;             // step the moved lane in this round too
            }
            i++;
        }
    }
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

// -n 0 with --batch B > 1: reads the records in blocks, searches each block with
// Batch_Exact and writes its lines in input order, as the per-read loop in main
// does. Encoding the reads counts as input here. Returns false on an unsupported
// read length.
bool Run_Batched_Exact(std::ifstream& in, std::ostream& out, unsigned B, unsigned maxHits, bool quiet,
                       uint64_t& records, uint64_t& reportedHits, double& tableBuildSeconds,
                       double& inputParseSeconds, double& searchSeconds, double& outputSeconds) {
    using Clock = std::chrono::steady_clock;
    constexpr unsigned BLOCK = 4096;
    std::vector<std::string> ids(BLOCK);
    std::vector<unsigned char> searched(BLOCK);  // valid read of the right length
    std::vector<unsigned> slot(BLOCK);            // its place among the searched reads
    std::vector<SARange> ranges(BLOCK);
    std::vector<unsigned char> found(BLOCK);
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

        const auto searchStart = Clock::now();
        Batch_Exact(reads.data(), count, STRINGLENGTH, B, ranges, found);
        searchSeconds += std::chrono::duration<double>(Clock::now() - searchStart).count();

        const auto outputStart = Clock::now();
        lines.clear();
        for (unsigned j = 0; j < n; j++) {
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
    "Usage: %s -q <reads.fa/fq> -g <genome_prefix> [-n <0|1>] [-o <hits.txt>] [-m <maxhits>] [--no-ref] [--lookup <k>] [--batch <B>] [--quiet|-Q] [--timing]\n"
        "  -n  maximum mismatches (default 0); 1 also needs the forward index <genome_prefix>.bwt/.fmv\n"
        "  -m  maximum hits per read (default 1, 0 = all)\n"
        "  --no-ref  do not load the SA samples (.sa) and packed genome (.pac); single-row\n"
        "            paths then walk the index to the end (less memory, slower)\n"
        "  --lookup  depth k of the start-of-search tables, 3-14 (default 6): 8*4^k bytes\n"
        "            per index; deeper is faster once there are more than about 4^k/3 reads\n"
        "            (see benchmark/lookup_table.md)\n"
        "  --batch   with -n 0, search B reads at a time, 1-64 (default 1), so that their\n"
        "            memory accesses overlap (see benchmark/batched_reads.md)\n"
        "Output, one line per SA range or located hit:\n"
        "  exact:      query_id\tF\tSA_start\tSA_end\thit_count\n"
        "  1 mismatch: query_id\tF|B\tSA_start\tSA_end\thit_count\t1\n"
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

    const bool batched = batchSize > 1 && maxMismatches == 0;
    if (batched && !Run_Batched_Exact(in, *out, batchSize, maxHits, quiet, records, reportedHits, tableBuildSeconds,
                                      inputParseSeconds, searchSeconds, outputSeconds)) {
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
        if (maxMismatches == 1 && LH > LOOKUPSIZE) {
            unsigned long readHits = One_Mismatch_Search(lh, rh, range, maxHits, exactRight, rightMismatch, leftMismatch);
            std::string lines;
            for (const SARange& hit : exactRight.hits) {
                unsigned long hitCount = hit.Skip ? 1UL : (hit.End >= hit.Start ? (hit.End - hit.Start + 1UL) : 0UL);
                if (maxHits > 0 && hitCount > maxHits) {
                    hitCount = maxHits;
                }
                lines += id + Hit_Text(hit, "F") + "\t" + std::to_string(hitCount) + "\n";
            }

            const ForwardResults* stages[2] = {&rightMismatch, &leftMismatch};
            const char* dirs[2] = {"F", "B"};
            for (int s = 0; s < 2; s++) {
                for (const SARange& hit : stages[s]->hits) {
                    if (maxHits > 0 && readHits >= maxHits) break;
                    unsigned long hitCount = hit.Skip ? 1UL : hit.End - hit.Start + 1UL;
                    if (maxHits > 0 && readHits + hitCount > maxHits) {
                        hitCount = maxHits - readHits;
                    }
                    readHits += hitCount;
                    lines += id + Hit_Text(hit, dirs[s]) + "\t" + std::to_string(hitCount) + "\t1\n";
                }
            }
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
    }


    return 0;
}
