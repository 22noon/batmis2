#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>

extern "C" {
#include "AwFmIndex.h"
}

namespace {

void Usage(const char* prog) {
    std::fprintf(
        stderr,
        "Usage: %s -f <reference.fa/fasta> -o <index.awfmi> "
        "[-k <kmer_seed_len>] [-s <sa_compression>] [--keep-sa-in-memory] [--store-seq]\n",
        prog);
}

bool ParseUint8(const char* text, uint8_t& out) {
    if (!text || !*text) {
        return false;
    }
    char* end = nullptr;
    const unsigned long val = std::strtoul(text, &end, 10);
    if (!end || *end != '\0' || val > 255UL) {
        return false;
    }
    out = static_cast<uint8_t>(val);
    return true;
}

}  // namespace

int main(int argc, char** argv) {
    std::string fastaPath;
    std::string awfmiPath;

    AwFmIndexConfiguration config{};
    config.suffixArrayCompressionRatio = 8;
    config.kmerLengthInSeedTable = 5;
    config.alphabetType = AwFmAlphabetDna;
    config.keepSuffixArrayInMemory = false;
    config.storeOriginalSequence = false;

    for (int i = 1; i < argc; ++i) {
        if ((std::strcmp(argv[i], "-f") == 0 || std::strcmp(argv[i], "--fasta") == 0) && i + 1 < argc) {
            fastaPath = argv[++i];
        } else if ((std::strcmp(argv[i], "-o") == 0 || std::strcmp(argv[i], "--output") == 0) && i + 1 < argc) {
            awfmiPath = argv[++i];
        } else if ((std::strcmp(argv[i], "-k") == 0 || std::strcmp(argv[i], "--kmer-seed") == 0) && i + 1 < argc) {
            if (!ParseUint8(argv[++i], config.kmerLengthInSeedTable)) {
                std::fprintf(stderr, "Invalid k-mer seed length: %s\n", argv[i]);
                return 1;
            }
        } else if ((std::strcmp(argv[i], "-s") == 0 || std::strcmp(argv[i], "--sa-compression") == 0) && i + 1 < argc) {
            if (!ParseUint8(argv[++i], config.suffixArrayCompressionRatio)) {
                std::fprintf(stderr, "Invalid SA compression ratio: %s\n", argv[i]);
                return 1;
            }
        } else if (std::strcmp(argv[i], "--keep-sa-in-memory") == 0) {
            config.keepSuffixArrayInMemory = true;
        } else if (std::strcmp(argv[i], "--store-seq") == 0) {
            config.storeOriginalSequence = true;
        } else if (std::strcmp(argv[i], "-h") == 0 || std::strcmp(argv[i], "--help") == 0) {
            Usage(argv[0]);
            return 0;
        } else {
            Usage(argv[0]);
            return 1;
        }
    }

    if (fastaPath.empty() || awfmiPath.empty()) {
        Usage(argv[0]);
        return 1;
    }

    AwFmIndex* index = nullptr;
    const AwFmReturnCode rc = awFmCreateIndexFromFasta(
        &index, &config, fastaPath.c_str(), awfmiPath.c_str());

    if (awFmReturnCodeIsFailure(rc)) {
        std::fprintf(
            stderr,
            "Failed to build index from %s to %s (AwFm return code %d)\n",
            fastaPath.c_str(), awfmiPath.c_str(), static_cast<int>(rc));
        return 1;
    }

    awFmDeallocIndex(index);

    std::fprintf(
        stdout,
        "Built AwFm index: %s\n"
        "  FASTA: %s\n"
        "  k-mer seed length: %u\n"
        "  SA compression ratio: %u\n"
        "  keep SA in memory: %s\n"
        "  store original sequence: %s\n",
        awfmiPath.c_str(),
        fastaPath.c_str(),
        static_cast<unsigned>(config.kmerLengthInSeedTable),
        static_cast<unsigned>(config.suffixArrayCompressionRatio),
        config.keepSuffixArrayInMemory ? "yes" : "no",
        config.storeOriginalSequence ? "yes" : "no");

    return 0;
}
