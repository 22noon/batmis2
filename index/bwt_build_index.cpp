/*
 * bwt_build_index.cpp
 *
 * Builds the BatMis forward and reverse FM indexes from a FASTA file in one
 * step, replacing scripts/build_index (filter + reverse + bwtformatdb x2 +
 * ann2loc.pl).
 *
 * Outputs for <genome> (written to <outdir>/<base> if an output directory is
 * given, otherwise next to the input):
 *   <genome>.N              non-ACGT runs per sequence
 *   <genome>.bin            concatenated filtered text (a/c/g/t)
 *   <genome>.ann.location   sequence names and lengths
 *   <genome>.pac            2-bit packed forward reference
 *   <genome>.bwt/.fmv/.sa   forward FM index
 *   rev<genome>.bwt/.fmv/.sa reverse FM index
 */
#include <cctype>
#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <string>
#include <vector>

#include <fcntl.h>
#include <getopt.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

extern "C" {
#include "TypeNLimit.h"
#include "MemManager.h"
#include "BWT.h"
#include "BWTConstruct.h"
#include "TextConverter.h"
}

namespace {

// bwtformatdb defaults
const unsigned int kPoolSize = 2097152;
const float kTargetNBit = 2.5f;
const unsigned int kInitialMaxBuildSize = 10000000;
const unsigned int kIncMaxBuildSize = 10000000;
const unsigned int kSaValueFreq = 8;

// References shorter than this are padded with a dummy sequence.
const unsigned int kMinTextLength = 16;
const char* const kPadSeqName = "__batmis_pad__";

const size_t kIoBufferSize = 16u << 20;
const size_t kReverseBlockSize = 16u << 20;

bool ShowProgress = false;

struct SeqInfo {
    std::string name;
    unsigned long long length;
};

// a/c/g/t -> 0..3; everything else is mapped to 'c' by the filter
unsigned char CharToCode[256];
unsigned char CharToFiltered[256];
bool IsAmbiguous[256];

void InitTables() {
    for (int i = 0; i < 256; i++) {
        CharToFiltered[i] = 'c';
        IsAmbiguous[i] = true;
        CharToCode[i] = 1;
    }
    const char* bases = "acgt";
    for (int k = 0; k < 4; k++) {
        unsigned char lc = bases[k];
        unsigned char uc = lc - 'a' + 'A';
        CharToFiltered[lc] = CharToFiltered[uc] = lc;
        IsAmbiguous[lc] = IsAmbiguous[uc] = false;
        CharToCode[lc] = CharToCode[uc] = (unsigned char)k;
    }
}

double Now() {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec / 1e9;
}

double StartTime;

void Log(const char* msg) {
    fprintf(stdout, "[%9.2fs] %s\n", Now() - StartTime, msg);
    fflush(stdout);
}

FILE* OpenFile(const std::string& name, const char* mode) {
    FILE* f = fopen64(name.c_str(), mode);
    if (!f) {
        fprintf(stderr, "Cannot open %s: %s\n", name.c_str(), strerror(errno));
        exit(1);
    }
    return f;
}

void MakeDirs(const std::string& path) {
    std::string cur;
    for (size_t i = 0; i <= path.size(); i++) {
        if (i == path.size() || path[i] == '/') {
            if (!cur.empty() && mkdir(cur.c_str(), 0777) != 0 && errno != EEXIST) {
                fprintf(stderr, "Cannot create directory %s: %s\n", cur.c_str(), strerror(errno));
                exit(1);
            }
        }
        if (i < path.size()) cur += path[i];
    }
}

// Appends 2-bit codes high bits first, the layout ConvertTextToBytePacked uses.
struct PackedText {
    std::vector<unsigned char> bytes;
    unsigned long long length = 0;

    void Append(unsigned char code) {
        unsigned shift = 6 - 2 * (length & 3);
        if ((length & 3) == 0) bytes.push_back(0);
        bytes.back() |= (unsigned char)(code << shift);
        length++;
    }
};

void CheckLength(unsigned long long length) {
    if (length >= 0xFFFFFFFFull) {
        fprintf(stderr, "Reference is too long (%llu bases); the FM index supports < 2^32.\n", length);
        exit(1);
    }
}

// Port of src/filter.cpp (default mode): writes .N and .bin, packs the
// filtered text and collects sequence names and lengths.
void FilterFasta(const std::string& fastaName, const std::string& prefix,
                 PackedText& packed, std::vector<SeqInfo>& seqs) {
    FILE* in = OpenFile(fastaName, "r");
    FILE* nFile = OpenFile(prefix + ".N", "w");
    FILE* binFile = OpenFile(prefix + ".bin", "w");
    setvbuf(binFile, NULL, _IOFBF, kIoBufferSize);

    char* line = NULL;
    size_t cap = 0;
    ssize_t len;
    unsigned long long offset = 0;
    bool inAmb = false;
    std::vector<char> out;

    while ((len = getline(&line, &cap, in)) != -1) {
        if (line[0] == '>') {
            if (inAmb) {
                fprintf(nFile, "%llu\n", offset);
                inAmb = false;
            }
            fputs(line, nFile);
            if (line[len - 1] != '\n') fputc('\n', nFile);
            if (!seqs.empty()) seqs.back().length = offset;

            size_t s = 1, e = 1;
            while (e < (size_t)len && !isspace((unsigned char)line[e])) e++;
            seqs.push_back(SeqInfo{std::string(line + s, e - s), 0});
            offset = 0;
            continue;
        }
        out.clear();
        for (ssize_t i = 0; i < len && line[i] != '\n' && line[i] != '\r'; i++) {
            unsigned char c = (unsigned char)line[i];
            if (IsAmbiguous[c]) {
                if (!inAmb) {
                    fprintf(nFile, "%llu:", offset);
                    inAmb = true;
                }
            } else if (inAmb) {
                fprintf(nFile, "%llu\n", offset);
                inAmb = false;
            }
            out.push_back((char)CharToFiltered[c]);
            packed.Append(CharToCode[c]);
            offset++;
        }
        if (!out.empty()) fwrite(out.data(), 1, out.size(), binFile);
    }
    if (inAmb) fprintf(nFile, "%llu\n", offset);
    if (!seqs.empty()) seqs.back().length = offset;

    // The BWT library crashes on texts of 2-3 characters, so very short
    // references get a dummy sequence appended last. Real coordinates are
    // unchanged; hits on kPadSeqName should be ignored downstream.
    if (packed.length > 0 && packed.length < kMinTextLength) {
        unsigned long long padLength = kMinTextLength - packed.length;
        fprintf(nFile, ">%s\n", kPadSeqName);
        seqs.push_back(SeqInfo{kPadSeqName, padLength});
        for (unsigned long long i = 0; i < padLength; i++) {
            fputc('a', binFile);
            packed.Append(CharToCode['a']);
        }
        fprintf(stdout, "Reference shorter than %u bases: appended %llu-base padding sequence %s.\n",
                kMinTextLength, padLength, kPadSeqName);
    }

    free(line);
    fclose(in);
    fclose(nFile);
    if (fclose(binFile) != 0) {
        fprintf(stderr, "Error writing %s.bin\n", prefix.c_str());
        exit(1);
    }
    CheckLength(packed.length);
    if (packed.length == 0) {
        fprintf(stderr, "No sequence data found in %s\n", fastaName.c_str());
        exit(1);
    }
}

// Replaces scripts/ann2loc.pl.
void WriteLocation(const std::string& prefix, const std::vector<SeqInfo>& seqs) {
    FILE* f = OpenFile(prefix + ".ann.location", "w");
    fprintf(f, "0\n");
    for (const SeqInfo& s : seqs) {
        fprintf(f, "%s\n%llu\n", s.name.c_str(), s.length);
    }
    fclose(f);
}

// Packs the reverse of the .bin text. Uses a memory map of the file when
// possible, otherwise reads it in blocks from the end.
void PackReverse(const std::string& binName, PackedText& packed) {
    int fd = open(binName.c_str(), O_RDONLY);
    if (fd < 0) {
        fprintf(stderr, "Cannot open %s: %s\n", binName.c_str(), strerror(errno));
        exit(1);
    }
    struct stat64 st;
    if (fstat64(fd, &st) != 0) {
        fprintf(stderr, "Cannot stat %s: %s\n", binName.c_str(), strerror(errno));
        exit(1);
    }
    unsigned long long n = (unsigned long long)st.st_size;
    packed.bytes.reserve((n + 3) / 4);

    void* map = getenv("BWT_BUILD_INDEX_NO_MMAP") ? MAP_FAILED
                : mmap(NULL, n, PROT_READ, MAP_PRIVATE, fd, 0);
    if (map != MAP_FAILED) {
        Log("Reversing text (memory mapped)..");
        const unsigned char* text = (const unsigned char*)map;
        // Walking backwards; the kernel still reads ahead within each region.
        madvise(map, n, MADV_SEQUENTIAL);
        for (unsigned long long i = n; i > 0; i--) {
            packed.Append(CharToCode[text[i - 1]]);
        }
        munmap(map, n);
    } else {
        Log("Reversing text (block reads)..");
        std::vector<unsigned char> buf(kReverseBlockSize);
        unsigned long long end = n;
        while (end > 0) {
            size_t chunk = (size_t)(end < kReverseBlockSize ? end : kReverseBlockSize);
            off64_t start = (off64_t)(end - chunk);
            size_t got = 0;
            while (got < chunk) {
                ssize_t r = pread64(fd, buf.data() + got, chunk - got, start + got);
                if (r <= 0) {
                    fprintf(stderr, "Read error on %s\n", binName.c_str());
                    exit(1);
                }
                got += (size_t)r;
            }
            for (size_t i = chunk; i > 0; i--) {
                packed.Append(CharToCode[buf[i - 1]]);
            }
            end -= chunk;
        }
    }
    close(fd);
}

// BWTIncCreate rejects small texts whose workspace estimate (targetNBit bits
// per char) cannot cover the fixed BWT/occ overhead; the actual allocation is
// clamped to BWTINC_MIN_MEMORY_IN_WORD anyway, so raise targetNBit for them.
float EffectiveTargetNBit(unsigned int n) {
    unsigned int rounded = (n + OCC_INTERVAL - 1) / OCC_INTERVAL * OCC_INTERVAL;
    unsigned int required = BWTResidentSizeInWord(n) + WORD_BETWEEN_OCC / 2
                          + BWTOccValueMinorSizeInWord(n) + 1;
    float needed = (float)required * BITS_IN_WORD / rounded;
    return needed > kTargetNBit ? needed + 0.01f : kTargetNBit;
}

// Same sequence of steps as bwtformatdb, minus FASTA parsing and the cached SA index.
void BuildFmIndex(MMPool* mmPool, PackedText& packed, const std::string& prefix, bool keepPac) {
    std::string pac = prefix + ".pac";
    std::string bwtName = prefix + ".bwt";
    std::string fmvName = prefix + ".fmv";
    std::string saName = prefix + ".sa";

    SaveBytePacked(pac.c_str(), packed.bytes.data(), (unsigned int)packed.length, ALPHABET_SIZE);
    float targetNBit = EffectiveTargetNBit((unsigned int)packed.length);
    std::vector<unsigned char>().swap(packed.bytes);
    if (targetNBit != kTargetNBit) {
        fprintf(stdout, "Small reference: targetNBit raised to %.2f\n", targetNBit);
    }

    Log(("Building BWT for " + prefix + "..").c_str());
    BWTInc* bwtInc = BWTIncConstructFromPacked(mmPool, pac.c_str(), ShowProgress,
                                               targetNBit, kInitialMaxBuildSize, kIncMaxBuildSize);
    BWTSaveBwtCodeAndOcc(bwtInc->bwt, bwtName.c_str(), fmvName.c_str());
    BWTIncFree(mmPool, bwtInc);

    Log("Building SA values..");
    BWT* bwt = BWTLoad(mmPool, bwtName.c_str(), fmvName.c_str(), NULL, NULL, NULL, NULL);
    BWTGenerateSaValue(bwt, kSaValueFreq, ShowProgress ? bwt->textLength / kSaValueFreq / 10 : 0);
    BWTSaveSaValue(bwt, saName.c_str());
    BWTFree(mmPool, bwt);

    if (!keepPac) remove(pac.c_str());
    Log(("Finished " + prefix).c_str());
}

void Usage(const char* prog) {
    fprintf(stderr,
            "Usage: %s [options] <genome.fa> [outdir]\n"
            "Builds the forward and reverse FM indexes used by BatMis.\n\n"
            "  -f, --forward-only   build only the forward index\n"
            "  -r, --reverse-only   build only the reverse index\n"
            "  -o, --outdir DIR     write index files to DIR\n"
            "  -p, --progress       show BWT construction progress\n"
            "  -h, --help           show this help\n",
            prog);
}

}  // namespace

int main(int argc, char** argv) {
    bool forward = true, reverse = true;
    bool forwardOnly = false, reverseOnly = false;
    std::string outDir;

    static const option longOptions[] = {
        {"forward-only", no_argument, NULL, 'f'},
        {"reverse-only", no_argument, NULL, 'r'},
        {"outdir", required_argument, NULL, 'o'},
        {"progress", no_argument, NULL, 'p'},
        {"help", no_argument, NULL, 'h'},
        {0, 0, 0, 0}};
    int opt;
    while ((opt = getopt_long(argc, argv, "fro:ph", longOptions, NULL)) != -1) {
        switch (opt) {
            case 'f': forwardOnly = true; break;
            case 'r': reverseOnly = true; break;
            case 'o': outDir = optarg; break;
            case 'p': ShowProgress = true; break;
            case 'h': Usage(argv[0]); return 0;
            default: Usage(argv[0]); return 1;
        }
    }
    if (forwardOnly && reverseOnly) {
        fprintf(stderr, "--forward-only and --reverse-only are mutually exclusive.\n");
        return 1;
    }
    if (forwardOnly) reverse = false;
    if (reverseOnly) forward = false;

    int positional = argc - optind;
    if (positional < 1 || positional > 2) {
        Usage(argv[0]);
        return 1;
    }
    std::string fasta = argv[optind];
    if (positional == 2) {
        if (!outDir.empty()) {
            fprintf(stderr, "Output directory given twice.\n");
            return 1;
        }
        outDir = argv[optind + 1];
    }

    std::string dir, base;
    size_t slash = fasta.find_last_of('/');
    if (slash == std::string::npos) {
        base = fasta;
    } else {
        dir = fasta.substr(0, slash + 1);
        base = fasta.substr(slash + 1);
    }
    if (!outDir.empty()) {
        MakeDirs(outDir);
        dir = outDir;
        if (dir.back() != '/') dir += '/';
    }
    std::string prefix = dir + base;
    std::string revPrefix = dir + "rev" + base;

    StartTime = Now();
    InitTables();
    MMMasterInitialize(1, 0, FALSE, NULL);
    MMPool* mmPool = MMPoolCreate(kPoolSize);

    std::vector<SeqInfo> seqs;
    {
        PackedText packed;
        Log(("Filtering " + fasta + "..").c_str());
        FilterFasta(fasta, prefix, packed, seqs);
        WriteLocation(prefix, seqs);
        fprintf(stdout, "%zu sequences, %llu bases.\n", seqs.size(), packed.length);
        // The forward .pac is the 2-bit reference and is always kept.
        if (forward) {
            BuildFmIndex(mmPool, packed, prefix, true);
        } else {
            SaveBytePacked((prefix + ".pac").c_str(), packed.bytes.data(), (unsigned int)packed.length, ALPHABET_SIZE);
        }
    }
    if (reverse) {
        PackedText packed;
        PackReverse(prefix + ".bin", packed);
        BuildFmIndex(mmPool, packed, revPrefix, false);
    }

    MMPoolFree(mmPool);
    Log("All done.");
    return 0;
}
