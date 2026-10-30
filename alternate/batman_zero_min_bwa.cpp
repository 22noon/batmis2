// batman_zero_min_bwa: exact-match search using a BWA (bwa-fmd) .bwt index.
// The native BatMis FM-index variant lives in batman_zero_min.cpp.
#include <cassert>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <chrono>

extern "C" {
#include "bwt.h"
}

namespace {

constexpr int MAXTAG = 4096;
constexpr int MAX_MISMATCHES_BOUND = 16;
constexpr int SAINTERVAL = 8;

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

using ActiveBWT = bwt_t;
ActiveBWT* revfmi = nullptr;

int LOOKUPSIZE = 6;
int STRINGLENGTH = 36;

unsigned* Forward_Start_LookupX = nullptr;
unsigned* Forward_End_LookupX = nullptr;

int EncodeBase(char c) {
    switch (std::toupper(static_cast<unsigned char>(c))) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default: return -1;
    }
}

ActiveBWT* initFMI(const char* bwtCodeFileName) {
    return bwt_restore_bwt(bwtCodeFileName);
}

void ExtendRange(unsigned long& start, unsigned long& end,
                 long newChar, ActiveBWT* fmi) {
    const bwtint_t left = start == 0
        ? static_cast<bwtint_t>(-1)
        : static_cast<bwtint_t>(start - 1);
    const bwtint_t right = static_cast<bwtint_t>(end);
    bwtint_t before = 0;
    bwtint_t through = 0;
    bwt_2occ(fmi, left, right, static_cast<ubyte_t>(newChar), &before, &through);
    const bwtint_t nextStart = fmi->L2[newChar] + 1 + before;
    const bwtint_t nextEnd = fmi->L2[newChar] + 1 + through;
    if (nextEnd <= nextStart) {
        start = 0;
        end = 0;
    } else {
        start = static_cast<unsigned long>(nextStart);
        end = static_cast<unsigned long>(nextEnd - 1);
    }
}

void Get_SARange_Fast(long newChar, SARange& range, ActiveBWT* fmi) {
    ExtendRange(range.Start, range.End, newChar, fmi);
}

Range InitialBaseRange(ActiveBWT* fmi, int base) {
    Range range{};
    const bwtint_t start = fmi->L2[base] + 1;
    const bwtint_t end = fmi->L2[base + 1] + 1;
    if (end > start) {
        range.Start = static_cast<unsigned>(start);
        range.End = static_cast<unsigned>(end - 1);
    }
    return range;
}

bool ExtendBwaSingleton(ActiveBWT* fmi, unsigned long row, unsigned char base,
                        unsigned long& nextStart) {
    const bwtint_t rankEnd = row == 0
        ? static_cast<bwtint_t>(-1)
        : static_cast<bwtint_t>(row - 1);
    const bwtint_t rank = bwt_occ(fmi, rankEnd, base);
    const bwtint_t start = fmi->L2[base] + 1 + rank;
    if (row == fmi->primary) return false;
    const bwtint_t packedRow = row - (row > fmi->primary);
    if (static_cast<unsigned char>(bwt_B0(fmi, packedRow)) != base) return false;
    nextStart = static_cast<unsigned long>(start);
    return true;
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
}

bool Search_Forwards_Exact(const unsigned char* currentTag, SARange& tag, int stringLength, ActiveBWT* fmi) {
    unsigned long index, now;

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
                unsigned long nextStart = 0;
                now = currentTag[tag.Level - 1];
                const bool matches = ExtendBwaSingleton(
                    fmi, index, static_cast<unsigned char>(now), nextStart);

                if (matches) {
                    tag.Start = nextStart;
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
            Get_SARange_Fast(currentTag[tag.Level - 1], tag, fmi);

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

void Usage(const char* prog) {
    std::fprintf(stderr,
    "Usage: %s -q <reads.fa/fq> -b <BWA .bwt> [-o <hits.txt>] [-m <maxhits>] [--quiet|-Q] [--timing]\n"
        "       (-g <genome_prefix> is accepted for compatibility and ignored)\n"
        "Output: query_id\tF\tSA_start\tSA_end\thit_count\n",
        prog);
}

}  // namespace

int main(int argc, char** argv) {
    std::string queryFile;
    std::string bwaBwtFile;
    std::string outputFile;
    unsigned maxHits = 1;
    bool quiet = false;
    bool timing = false;

    for (int i = 1; i < argc; i++) {
        if (std::strcmp(argv[i], "-q") == 0 && i + 1 < argc) {
            queryFile = argv[++i];
        } else if (std::strcmp(argv[i], "-g") == 0 && i + 1 < argc) {
            ++i;  // unused with a BWA index
        } else if (std::strcmp(argv[i], "-b") == 0 && i + 1 < argc) {
            bwaBwtFile = argv[++i];
        } else if (std::strcmp(argv[i], "-o") == 0 && i + 1 < argc) {
            outputFile = argv[++i];
        } else if (std::strcmp(argv[i], "-m") == 0 && i + 1 < argc) {
            maxHits = static_cast<unsigned>(std::strtoul(argv[++i], nullptr, 10));
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

    if (queryFile.empty() || bwaBwtFile.empty()) {
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

    const auto loadStart = Clock::now();
    revfmi = initFMI(bwaBwtFile.c_str());
    indexLoadSeconds = std::chrono::duration<double>(Clock::now() - loadStart).count();
    if (!revfmi) {
        std::fprintf(stderr, "Failed to load BWA BWT index: %s\n", bwaBwtFile.c_str());
        return 1;
    }

    const unsigned lookupArraySize = (2u << (2 * LOOKUPSIZE));
    Forward_Start_LookupX = static_cast<unsigned*>(std::malloc(sizeof(unsigned) * lookupArraySize));
    Forward_End_LookupX = static_cast<unsigned*>(std::malloc(sizeof(unsigned) * lookupArraySize));
    if (!Forward_Start_LookupX || !Forward_End_LookupX) {
        std::fprintf(stderr, "Out of memory for lookup tables\n");
        return 1;
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

    bwt_destroy(revfmi);

    return 0;
}
