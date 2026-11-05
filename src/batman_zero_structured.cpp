// Structured single-index exact mapper using the native BatMis FM index.
// (rev<genome>.bwt/.fmv built by index/bwt_build_index).
#include "exact_search.hpp"
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
#include "MemManager.h"
#include "BWT.h"
}

namespace {

constexpr int MAXTAG = 4096;
MMPool* mmPool = nullptr;
using ActiveBWT = BWT;
ActiveBWT* revfmi = nullptr;

int STRINGLENGTH = 36;


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

void BuildIndexPaths(const std::string& genomePrefix, std::string& revBwt, std::string& revOcc) {
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
}

void Usage(const char* prog) {
    std::fprintf(stderr,
    "Usage: %s -q <reads.fa/fq> -g <genome_prefix> [-o <hits.txt>] [-m <maxhits>] [--quiet|-Q] [--timing]\n"
        "Output: query_id\tF\tSA_start\tSA_end\thit_count\n",
        prog);
}

}  // namespace

int main(int argc, char** argv) {
    std::string queryFile;
    std::string genomePrefix;
    std::string outputFile;
    unsigned maxHits = 1;
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
    BuildIndexPaths(genomePrefix, revBwt, revOcc);

    const auto loadStart = Clock::now();
    revfmi = initFMI(revBwt.c_str(), revOcc.c_str());
    indexLoadSeconds = std::chrono::duration<double>(Clock::now() - loadStart).count();
    if (!revfmi) {
        std::fprintf(stderr, "Failed to load reverse index: %s / %s\n", revBwt.c_str(), revOcc.c_str());
        return 1;
    }

    batmis::ExactSearch exactSearch(revfmi);

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
            if (!exactSearch.prepare(STRINGLENGTH)) {
                std::fprintf(stderr, "Cannot prepare exact search for read length %d\n", STRINGLENGTH);
                return 1;
            }
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

        batmis::ExactResult range{};
        const bool found = exactSearch.search(encoded, range);
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
