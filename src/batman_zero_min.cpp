// batman_zero_min: exact-match search using the native BatMis FM index
// (rev<genome>.bwt/.fmv built by index/bwt_build_index). With -n 1 it also
// finds one-mismatch hits, using the forward index <genome>.bwt/.fmv too.
// The BWA-index variant lives in batman_zero_min_bwa.cpp.
// Optional build variant: -DBATMAN_ZERO_MIN_USE_PAIRED_RANK.
#include <cassert>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
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
constexpr int BRANCHTHRESHOLD = 0;

struct SARange {
    unsigned long Start;
    unsigned long End;
    int Level;
    char Mismatches;
    unsigned char Skip;
    int Tag;
    unsigned char Mismatch_Pos[MAX_MISMATCHES_BOUND];
    unsigned Mismatch_Char;
};

struct Range {
    unsigned Start;
    unsigned End;
    int Label;
};

// Output of Search_Forwards. The caller clears it before each search.
struct ForwardResults {
    std::vector<SARange> hits;        // exactly Count mismatches
    std::vector<SARange> mismatches;  // Count+1 mismatches; resume at Level+1 (complete if Level == StringLength)
    unsigned long hitCount = 0;       // locations in hits
};

MMPool* mmPool = nullptr;
using ActiveBWT = BWT;
ActiveBWT* revfmi = nullptr;
ActiveBWT* fwfmi = nullptr;  // loaded only for -n 1

int LOOKUPSIZE = 6;
int STRINGLENGTH = 36;

unsigned* Forward_Start_LookupX = nullptr;
unsigned* Forward_End_LookupX = nullptr;
unsigned* Backward_Start_LookupX = nullptr;
unsigned* Backward_End_LookupX = nullptr;
unsigned Branch_Characters[4] = {0, 0, 0, 0};
SARange Branch_Ranges[4];

int EncodeBase(char c) {
    switch (std::toupper(static_cast<unsigned char>(c))) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default: return -1;
    }
}

ActiveBWT* initFMI(const char* bwtCodeFileName, const char* bwtOccValueFileName) {
    const int poolSize = 524288;
    MMMasterInitialize(3, 0, FALSE, nullptr);
    mmPool = MMPoolCreate(poolSize);
    return BWTLoad(mmPool, bwtCodeFileName, bwtOccValueFileName, nullptr, nullptr, nullptr, nullptr);
}

void ExtendRange(unsigned long& start, unsigned long& end,
                 long newChar, ActiveBWT* fmi) {
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

void Get_SARange_Fast(long newChar, SARange& range, ActiveBWT* fmi) {
    ExtendRange(range.Start, range.End, newChar, fmi);
}

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
            ExtendRange(start, end, bit, revfmi);
            range.Start = static_cast<unsigned>(start);
            range.End = static_cast<unsigned>(end);
        }
        Forward_Start_LookupX[range.Label] = range.Start;
        Forward_End_LookupX[range.Label] = range.End;
        assert(Forward_Start_LookupX[range.Label] <= Forward_End_LookupX[range.Label]);
    } else {
        assert(range.Start <= range.End);
        range.Label = range.Label | (bit << (2 * (level - 1)));
        if (range.Start) {
            unsigned long start = range.Start;
            unsigned long end = range.End;
            ExtendRange(start, end, bit, revfmi);
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
        ExtendRange(start, end, bit, fwfmi);
        range.Start = static_cast<unsigned>(start);
        range.End = static_cast<unsigned>(end);
    }
    if (LOOKUPSIZE == level) {
        Backward_Start_LookupX[range.Label] = range.Start;
        Backward_End_LookupX[range.Label] = range.End;
        assert(Backward_Start_LookupX[range.Label] <= Backward_End_LookupX[range.Label]);
    } else {
        level++;
        for (int i = 0; i < 4; i++) {
            Build_Preindex_Backward(range, level, i);
        }
    }
}

void Build_Tables(ActiveBWT* fmi) {
    if (STRINGLENGTH < 28) {
        LOOKUPSIZE = 3;
    }

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

bool Search_Forwards_Exact(const unsigned char* currentTag, SARange& tag, int stringLength, ActiveBWT* fmi) {
    unsigned long index,first, last;
    unsigned int now;

    if (!tag.Start) {
        return false;
    }

    for (;;) {
        if (tag.End == tag.Start || tag.Skip) {
            if (tag.Start % SAINTERVAL == 0 && !tag.Skip) {
                tag.Skip++;
                tag.End = tag.Start;
            }

            for (;;) {
                index = tag.Start;
                if (index >= fmi->inverseSa0) {
                    index--;
                }
                now = fmi->bwtCode[index / 16] << ((index % 16) * 2) >> (BITS_IN_WORD - 2);
                if (currentTag[tag.Level - 1] == now) {
                    tag.Start = fmi->cumulativeFreq[now] + BWTOccValue(fmi, tag.Start, now) + 1;
                    if (tag.Skip) {
                        tag.Skip++;
                    } else if (tag.Start % SAINTERVAL == 0) {
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
        } else {
            if (tag.End - tag.Start < BRANCHTHRESHOLD) {
                Branch_Characters[0] = 0;
                Branch_Characters[1] = 0;
                Branch_Characters[2] = 0;
                Branch_Characters[3] = 0;

                if (tag.Start + 1 >= fmi->inverseSa0) {
                    first = tag.Start;
                    last = tag.End;
                } else {
                    first = tag.Start + 1;
                    last = tag.End + 1;
                }

                for (unsigned long pos = first; pos <= last; pos++) {
                    now = fmi->bwtCode[(pos - 1) / 16] << (((pos - 1) % 16) * 2) >> (BITS_IN_WORD - 2);
                    Branch_Characters[now]++;
                }

                now = currentTag[tag.Level - 1];
                if (Branch_Characters[now]) {
                    tag.Start = fmi->cumulativeFreq[now] + BWTOccValue(fmi, tag.Start, now) + 1;
                    tag.End = tag.Start + Branch_Characters[now] - 1;
                } else {
                    tag.Start = 0;
                }
            } else {
                Get_SARange_Fast(currentTag[tag.Level - 1], tag, fmi);
            }

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
		unsigned Last, First;
		char Now;

		if (Tag.Start+1 >= fmi->inverseSa0) {First=Tag.Start;Last=Tag.End;} else {First=Tag.Start+1;Last=Tag.End+1;} 

		for (unsigned long Pos=First;Pos<=Last;Pos++)
		{
			Now=fmi->bwtCode[(Pos-1) / 16] << (((Pos-1) % 16) * 2)>> (BITS_IN_WORD - 2);
			Branch_Characters[Now]++;	
		}

		for (int Branch=0;Branch<4;Branch++)
		{
			if (Branch_Characters[Branch])
			{
				Branch_Ranges[Branch].Start = fmi->cumulativeFreq[Branch] + BWTOccValue(fmi, Tag.Start, Branch) + 1;
				Branch_Ranges[Branch].End = Branch_Ranges[Branch].Start + Branch_Characters[Branch]-1;// Calculate SAranges
			}
		}
	}
	else
	{
		for (int Branch=0;Branch<4;Branch++)
		{
				Branch_Ranges[Branch].Start = fmi->cumulativeFreq[Branch] + BWTOccValue(fmi, Tag.Start, Branch) + 1;
				Branch_Ranges[Branch].End = fmi->cumulativeFreq[Branch] + BWTOccValue(fmi, Tag.End+1, Branch);
				if(!(Branch_Ranges[Branch].End<Branch_Ranges[Branch].Start)) Branch_Characters[Branch]=1;
		}

	}
}
void Search_Forwards_OneSA(const unsigned char* Current_Tag,struct SARange & Tag,int Count,int Start,int StringLength,BWT *fmi,ForwardResults& Out)
{
	unsigned long Index,Now;
	if (Tag.Start==0) return;
	if(Tag.Start % SAINTERVAL == 0 && !Tag.Skip) 
	{
		Tag.Skip++;
		Tag.End=Tag.Start;
	}

	for(;;)
	{
		
		Index=Tag.Start;
		if (Index >= fmi->inverseSa0) Index--;//adjust for missing $
		Now=fmi->bwtCode[(Index) / 16] << (((Index) % 16) * 2)>> (BITS_IN_WORD - 2);
		Tag.Start = fmi->cumulativeFreq[Now] + BWTOccValue(fmi, Tag.Start, Now) + 1;

		if (Tag.Skip) Tag.Skip++;
		else if(Tag.Start % SAINTERVAL == 0) 
		{
			Tag.Skip++;
			Tag.End=Tag.Start;
		}

		if (Current_Tag[Tag.Level+Start]!=Now)
		{

			Tag.Mismatch_Char=Tag.Mismatch_Char | (Now<<Tag.Mismatches*2);
			Tag.Mismatch_Pos[Tag.Mismatches]=(Start+Tag.Level);
			Tag.Mismatches++;
			
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
		else//log 2 mismatches
		{
			if (!Tag.Skip) Tag.End=Tag.Start;
			Out.mismatches.push_back(Tag);
			return;
		}
	}
}

// Extends Tag allowing up to Count mismatches. Ranges with exactly Count
// mismatches go to Out.hits, ranges needing one more to Out.mismatches.
// Stops once Out.hitCount reaches MaxHits (0 = no limit).
void Search_Forwards(const unsigned char* Current_Tag,const struct SARange & Tag,int Count,int Start,int StringLength,BWT *fmi,ForwardResults& Out,unsigned MaxHits)
{
	if (!Tag.Start) return;
	Start=Start-2;//Adjust for offset difference
	std::vector<SARange> Stack;
	Stack.push_back(Tag);
	struct SARange Range,Temp_Range;
	while(!Stack.empty())//While Stack non-empty....
	{
		Range=Stack.back();
		Stack.pop_back();		//Pop the range
		if (Range.End==Range.Start || Range.Skip)//does this SArange have only one branch?
		{
			Search_Forwards_OneSA(Current_Tag,Range,Count,Start,StringLength,fmi,Out);
			if(MaxHits && Out.hitCount>=MaxHits) return;
		}
		else
		{
			Branch_Detect(Range,fmi,Start);//One_Branch(Range,revfmi);
			for(int Branch=0;Branch<4;Branch++)
			{
				if (Branch_Characters[Branch])//This character actually branches
				{
					Temp_Range=Range;
					Temp_Range.Start = Branch_Ranges[Branch].Start;//fmi->cumulativeFreq[Branch] + BWTOccValue(fmi, Temp_Range.Start, Branch) + 1;
					Temp_Range.End = Branch_Ranges[Branch].End;//Temp_Range.Start + Branch_Characters[Branch]-1;// Calculate SAranges

					if (Current_Tag[Temp_Range.Level+Start]!=Branch)
					{

						Temp_Range.Mismatch_Char=Temp_Range.Mismatch_Char | (Branch<<Temp_Range.Mismatches*2);
						Temp_Range.Mismatch_Pos[Temp_Range.Mismatches]=(Start+Temp_Range.Level);
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
					else
					{
						Out.mismatches.push_back(Temp_Range);
						continue;
					}
				} 
			}
		}
	}
	return;
}

// Extends Tag backwards over Current_Tag[Start-Level] with no mismatches.
// On return Tag.Start==0 if there is no match, otherwise Tag holds the SA range.
void Search_Backwards_Exact(const unsigned char* Current_Tag,struct SARange & Tag,int Start,int StringLength,BWT *fmi)
{

	int Level;
	unsigned long Index,Now,First,Last;
	if (!Tag.Start) return;

	for(;;)
	{
		if(Tag.End==Tag.Start || Tag.Skip)//Only one branch?
		{
			Level=Tag.Level;
			if(Tag.Start % SAINTERVAL == 0 && !Tag.Skip)
			{
				Tag.Skip++;
				Tag.End=Tag.Start;
			}
			for(;;)
			{
				Index=Tag.Start;
				if (Index >= fmi->inverseSa0) Index--;//adjust for missing $
				Now=fmi->bwtCode[(Index) / 16] << (((Index) % 16) * 2)>> (BITS_IN_WORD - 2);//FMIBwtValue(fmi,Index);
				if (Current_Tag[Start-Level] == Now)
				{
					Tag.Start = fmi->cumulativeFreq[Now] + BWTOccValue(fmi, Tag.Start, Now) + 1;

					if (Tag.Skip) Tag.Skip++;
					else if(Tag.Start % SAINTERVAL == 0)
					{
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
		else
		{
			if(Tag.End-Tag.Start<BRANCHTHRESHOLD)//only small number of branches
			{
				Branch_Characters[0]=0;Branch_Characters[1]=0;Branch_Characters[2]=0;Branch_Characters[3]=0;

				if (Tag.Start+1 >= fmi->inverseSa0) {First=Tag.Start;Last=Tag.End;} else {First=Tag.Start+1;Last=Tag.End+1;}
				for (unsigned long Pos=First;Pos<=Last;Pos++)
				{
					Now=fmi->bwtCode[(Pos-1) / 16] << (((Pos-1) % 16) * 2)>> (BITS_IN_WORD - 2);
					Branch_Characters[Now]++;
				}

				Now=Current_Tag[Start-Tag.Level];
				if (Branch_Characters[Now])//we have a match...
				{

					Tag.Start = fmi->cumulativeFreq[Now] + BWTOccValue(fmi, Tag.Start, Now) + 1;
					Tag.End = Tag.Start + Branch_Characters[Now]-1;// Calculate SAranges

				}
				else
				{
					Tag.Start=0;//Tag.End=0;
				}
			}
			else Get_SARange_Fast(Current_Tag[Start-Tag.Level],Tag,fmi);

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

void Search_Backwards_OneSA(const unsigned char* Current_Tag,struct SARange & Tag,int Count,int Start,int StringLength,BWT *fmi,ForwardResults& Out)
{
	unsigned long Index,Now;
	if (Tag.Start==0) return;
	if(Tag.Start % SAINTERVAL == 0 && !Tag.Skip)
	{
		Tag.Skip++;
		Tag.End=Tag.Start;
	}

	for(;;)
	{
		Index=Tag.Start;
		if (Index >= fmi->inverseSa0) Index--;//adjust for missing $
		Now=fmi->bwtCode[(Index) / 16] << (((Index) % 16) * 2)>> (BITS_IN_WORD - 2);//FMIBwtValue(fmi,Index);
		Tag.Start = fmi->cumulativeFreq[Now] + BWTOccValue(fmi, Tag.Start, Now) + 1;

		if (Tag.Skip) Tag.Skip++;
		else if(Tag.Start % SAINTERVAL == 0)
		{
			Tag.Skip++;
			Tag.End=Tag.Start;
		}

		if (Current_Tag[Start-Tag.Level] != Now)
		{
			Tag.Mismatch_Char=Tag.Mismatch_Char | (Now<<Tag.Mismatches*2);
			Tag.Mismatch_Pos[Tag.Mismatches]=(Start-Tag.Level);
			Tag.Mismatches++;

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
		else
		{
			if (!Tag.Skip) Tag.End=Tag.Start;
			Out.mismatches.push_back(Tag);
			return;
		}
	}
}

// Extends Tag backwards over Current_Tag[Start-Level] allowing up to Count
// mismatches. Results as for Search_Forwards.
void Search_Backwards(const unsigned char* Current_Tag,const struct SARange & Tag,int Count,int Start,int StringLength,BWT *fmi,ForwardResults& Out,unsigned MaxHits)
{
	if (!Tag.Start) return;
	std::vector<SARange> Stack;
	Stack.push_back(Tag);
	struct SARange Range,Temp_Range;
	while(!Stack.empty())//While Stack non-empty....
	{
		Range=Stack.back();
		Stack.pop_back();	//Pop the range
		if (Range.End==Range.Start || Range.Skip)//does this SArange have only one branch?
		{
			Search_Backwards_OneSA(Current_Tag,Range,Count,Start,StringLength,fmi,Out);
			if(MaxHits && Out.hitCount>=MaxHits) return;
		}
		else
		{
			Branch_Detect(Range,fmi,Start);
			for(int Branch=0;Branch<4;Branch++)
			{
				if (Branch_Characters[Branch])//This character actually branches
				{
					Temp_Range=Range;//adjust
					Temp_Range.Start = Branch_Ranges[Branch].Start;//fmi->cumulativeFreq[Branch] + BWTOccValue(fmi, Temp_Range.Start, Branch) + 1;
					Temp_Range.End = Branch_Ranges[Branch].End;//Temp_Range.Start + Branch_Characters[Branch]-1;// Calculate SAranges

					if (Current_Tag[Start-Temp_Range.Level] != Branch)
					{
						Temp_Range.Mismatch_Char=Temp_Range.Mismatch_Char | (Branch<<Temp_Range.Mismatches*2);
						Temp_Range.Mismatch_Pos[Temp_Range.Mismatches]=Start-Temp_Range.Level;
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
					else
					{
						Out.mismatches.push_back(Temp_Range);
						continue;
					}
				}
			}
		}
	}
	return;
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
                     std::string& fwBwt, std::string& fwOcc) {
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
}

void Usage(const char* prog) {
    std::fprintf(stderr,
    "Usage: %s -q <reads.fa/fq> -g <genome_prefix> [-n <0|1>] [-o <hits.txt>] [-m <maxhits>] [--quiet|-Q] [--timing]\n"
        "  -n  maximum mismatches (default 0); 1 also needs the forward index <genome_prefix>.bwt/.fmv\n"
        "  -m  maximum hits per read (default 1, 0 = all)\n"
        "Output, one line per SA range:\n"
        "  exact:      query_id\tF\tSA_start\tSA_end\thit_count\n"
        "  1 mismatch: query_id\tF|B\tSA_start\tSA_end\thit_count\t1\n"
        "  F = range in the reverse index (rev<genome_prefix>), B = range in the forward index\n"
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
    BuildIndexPaths(genomePrefix, revBwt, revOcc, fwBwt, fwOcc);

    const auto loadStart = Clock::now();
    revfmi = initFMI(revBwt.c_str(), revOcc.c_str());
    if (!revfmi) {
        std::fprintf(stderr, "Failed to load reverse index: %s / %s\n", revBwt.c_str(), revOcc.c_str());
        return 1;
    }
    if (maxMismatches > 0) {
        fwfmi = BWTLoad(mmPool, fwBwt.c_str(), fwOcc.c_str(), nullptr, nullptr, nullptr, nullptr);
        if (!fwfmi) {
            std::fprintf(stderr, "Failed to load forward index: %s / %s\n", fwBwt.c_str(), fwOcc.c_str());
            return 1;
        }
    }
    indexLoadSeconds = std::chrono::duration<double>(Clock::now() - loadStart).count();

    const unsigned lookupArraySize = (2u << (2 * LOOKUPSIZE));
    Forward_Start_LookupX = static_cast<unsigned*>(std::malloc(sizeof(unsigned) * lookupArraySize));
    Forward_End_LookupX = static_cast<unsigned*>(std::malloc(sizeof(unsigned) * lookupArraySize));
    if (!Forward_Start_LookupX || !Forward_End_LookupX) {
        std::fprintf(stderr, "Out of memory for lookup tables\n");
        return 1;
    }
    if (fwfmi) {
        Backward_Start_LookupX = static_cast<unsigned*>(std::malloc(sizeof(unsigned) * lookupArraySize));
        Backward_End_LookupX = static_cast<unsigned*>(std::malloc(sizeof(unsigned) * lookupArraySize));
        if (!Backward_Start_LookupX || !Backward_End_LookupX) {
            std::fprintf(stderr, "Out of memory for lookup tables\n");
            return 1;
        }
    }

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
    unsigned char encoded[MAXTAG] = {0};
    bool initialized = false;
    uint64_t records = 0;
    uint64_t reportedHits = 0;
    ForwardResults rightMismatch;  // 0|1 ranges, reverse index
    ForwardResults leftMismatch;   // 1|0 ranges, forward index

    for (;;) {
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
            encoded[i] = static_cast<unsigned char>(code);
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

        int c;
        if (LOOKUPSIZE == 3) {
            c = encoded[0] | (encoded[1] << 2) | (encoded[2] << 4);
        } else {
            c = encoded[0] | (encoded[1] << 2) | (encoded[2] << 4) | (encoded[3] << 6) | (encoded[4] << 8) | (encoded[5] << 10);
        }

        SARange range{};
        range.Start = Forward_Start_LookupX[c];
        range.End = Forward_End_LookupX[c];
        range.Tag = 0;
        range.Mismatches = 0;
        range.Level = LOOKUPSIZE + 1;
        range.Skip = 0;

        const int LH = STRINGLENGTH / 2;
        const int RH = STRINGLENGTH - LH;
        if (maxMismatches == 1 && LH > LOOKUPSIZE) {
            // BatMis one-mismatch scan: exact 0|0, then left half exact with a
            // mismatch in the right half (0|1), then right half exact with a
            // mismatch in the left half (1|0).
            std::string lines;
            unsigned long readHits = 0;

            // 0|0, stopping at the left half to keep its exact range for 0|1.
            SARange leftRange{};
            bool found = Search_Forwards_Exact(encoded, range, LH, revfmi) && range.Start != 0;
            if (found) {
                leftRange = range;
                range.Level++;
                found = Search_Forwards_Exact(encoded, range, STRINGLENGTH, revfmi) && range.Start != 0;
            }
            if (found) {
                unsigned long hitCount = range.Skip ? 1UL : (range.End >= range.Start ? (range.End - range.Start + 1UL) : 0UL);
                if (maxHits > 0 && hitCount > maxHits) {
                    hitCount = maxHits;
                }
                readHits += hitCount;
                lines += id + "\tF\t" + std::to_string(range.Start) + "\t" + std::to_string(range.End) + "\t" + std::to_string(hitCount) + "\n";
            }

            rightMismatch.hits.clear();
            rightMismatch.mismatches.clear();
            rightMismatch.hitCount = 0;
            leftMismatch.hits.clear();
            leftMismatch.mismatches.clear();
            leftMismatch.hitCount = 0;

            // 0|1
            if (leftRange.Start && (maxHits == 0 || readHits < maxHits)) {
                leftRange.Level = 1;
                leftRange.Mismatches = 0;
                leftRange.Mismatch_Char = 0;
                Search_Forwards(encoded, leftRange, 1, LH + 1, RH, revfmi, rightMismatch,
                                maxHits ? maxHits - readHits : 0);
            }

            // 1|0
            if (maxHits == 0 || readHits + rightMismatch.hitCount < maxHits) {
                int bc = 0;
                for (int i = 0; i < LOOKUPSIZE; i++) {
                    bc |= encoded[STRINGLENGTH - 1 - i] << (2 * i);
                }
                SARange rightRange{};
                rightRange.Start = Backward_Start_LookupX[bc];
                rightRange.End = Backward_End_LookupX[bc];
                rightRange.Level = LOOKUPSIZE + 1;
                Search_Backwards_Exact(encoded, rightRange, STRINGLENGTH, RH, fwfmi);
                if (rightRange.Start) {
                    rightRange.Level = 1;
                    Search_Backwards(encoded, rightRange, 1, LH, LH, fwfmi, leftMismatch,
                                     maxHits ? maxHits - readHits - rightMismatch.hitCount : 0);
                }
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
                    lines += id + "\t" + dirs[s] + "\t" + std::to_string(hit.Start) + "\t" +
                             std::to_string(hit.Skip ? hit.Start : hit.End) + "\t" + std::to_string(hitCount) + "\t1\n";
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

        bool found = Search_Forwards_Exact(encoded, range, STRINGLENGTH, revfmi);
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
            *out << id << "\tF\t" << range.Start << "\t" << range.End << "\t" << hitCount << "\n";
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
