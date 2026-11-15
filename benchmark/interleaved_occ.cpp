// interleaved_occ: does storing the Occ checkpoints inside the BWT (the BWA layout) make
// a BatMis rank step faster? Stand-alone; see interleaved_index.md.
//
// Loads a BatMis index (.bwt + .fmv), builds three more layouts of it in memory, checks
// every layout's Occ against the library's, and times LF walks (read the row's BWT
// character, one Occ, go to that row: the mappers' single-row step).
//
//   LIB0 the library: bwtCode, occValue, occValueMajor; BWTOccValue (SSE2 decode); no
//        prefetch, which is what batmis1/2 run today (their Prefetch_Row is compiled away)
//   LIB  the library, with Prefetch_Row's four prefetches
//   SEP  the same three arrays (checkpoint every 256 rows, 16-bit counts, major table),
//        counted with the popcount code the interleaved layouts use
//   B    64-byte blocks: 4 x 32-bit counts + 192 bases
//   C    64-byte blocks: 4 x 16-bit counts + 224 bases, and a major table
//
// SEP against B and C is the layout comparison; LIB shows where the mappers are now.
//
//   interleaved_occ <prefix> [--steps N] [--thp] [--threads a,b,..] [--chains a,b,..]
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <thread>
#include <vector>
#include <sys/mman.h>

extern "C" {
#include "MemManager.h"
#include "BWT.h"
}

typedef uint64_t u64;

// A function that only prefetches has no side effect GCC knows of: when it is not inlined,
// GCC (13.3, -O3) deletes the calls to it. Lib::prefetch, a copy of batmis2's Prefetch_Row,
// lost all four prefetches that way, so every prefetch here is forced inline.
#define PF __attribute__((always_inline)) inline

static const u64 LOW = 0x5555555555555555ULL;

// Counts base c among characters [a, b) of W (32 per word, character j in bits 2j, 2j+1).
static inline unsigned Count(const u64* W, unsigned a, unsigned b, unsigned c) {
    unsigned r = 0;
    const u64 rep = c * LOW;
    for (unsigned w = a / 32; w * 32 < b; w++) {
        const unsigned lo = std::max(a, w * 32) - w * 32, hi = std::min(b, w * 32 + 32) - w * 32;
        const u64 y = ~(W[w] ^ rep);
        const u64 mhi = hi >= 32 ? ~0ULL : (1ULL << (2 * hi)) - 1;
        const u64 mlo = (1ULL << (2 * lo)) - 1;  // lo < 32 here
        r += __builtin_popcountll(y & (y >> 1) & LOW & mhi & ~mlo);
    }
    return r;
}

static void* Alloc(size_t bytes, bool thp) {
    const size_t A = 2ul << 20;
    bytes = (bytes + A - 1) / A * A;
    void* p = aligned_alloc(A, bytes);
    if (!p) { std::perror("alloc"); std::exit(1); }
    if (thp) madvise(p, bytes, MADV_HUGEPAGE);
    std::memset(p, 0, bytes);
    return p;
}

struct Lib {
    const BWT* b;
    unsigned inv;
    unsigned ch(unsigned row) const {
        const unsigned i = row - (row > inv);
        return b->bwtCode[i / 16] << ((i % 16) * 2) >> 30;
    }
    unsigned occ(unsigned row, unsigned c) const { return BWTOccValue(b, row, c); }
    PF void prefetch(unsigned row) const {  // batmis2's Prefetch_Row
        const unsigned i = row - (row > inv);
        const unsigned long e = (i + OCC_INTERVAL / 2 - 1) / OCC_INTERVAL;
        __builtin_prefetch(b->occValue + e / OCC_VALUE_PER_WORD * ALPHABET_SIZE);
        __builtin_prefetch(b->occValueMajor + e * OCC_INTERVAL / OCC_INTERVAL_MAJOR * ALPHABET_SIZE);
        __builtin_prefetch(b->bwtCode + i / CHAR_PER_WORD);
        __builtin_prefetch(b->bwtCode + e * OCC_INTERVAL / CHAR_PER_WORD);
    }
};

struct Lib0 : Lib {
    PF void prefetch(unsigned) const {}
};

struct Sep {
    const u64* W;          // 8 words of slack in front
    const uint16_t* minor; // [4] per checkpoint (every 256 rows)
    const unsigned* major; // [4] per 65536 rows
    unsigned inv;
    unsigned ch(unsigned row) const {
        const unsigned i = row - (row > inv);
        return W[i / 32] >> (2 * (i % 32)) & 3;
    }
    unsigned occ(unsigned row, unsigned c) const {
        const unsigned i = row - (row > inv);
        const unsigned k = (i + 128) >> 8, p = k << 8;
        // All four counts are loaded, so the address does not wait for c (the row's BWT
        // character in an LF walk): the BWT and count misses overlap.
        u64 m;
        std::memcpy(&m, minor + k * 4ul, 8);
        const unsigned at = major[(k >> 8) * 4 + c] + (m >> (16 * c) & 0xFFFF);
        if (i >= p) return at + Count(W + p / 32, 0, i - p, c);
        return at - Count(W + p / 32 - 8, i - (p - 256), 256, c);
    }
    PF void prefetch(unsigned row) const {
        const unsigned i = row - (row > inv);
        __builtin_prefetch(W + i / 32);
        __builtin_prefetch(minor + ((i + 128) >> 8) * 4ul);
    }
};

// Block: q[0], q[1] = Occ at the block's character 96 (4 x u32); q[2..7] = 192 bases.
struct LayB {
    const u64* Q;
    unsigned inv;
    unsigned ch(unsigned row) const {
        const unsigned i = row - (row > inv), o = i % 192;
        return Q[(i / 192) * 8ul + 2 + o / 32] >> (2 * (o % 32)) & 3;
    }
    unsigned occ(unsigned row, unsigned c) const {
        const unsigned i = row - (row > inv), o = i % 192;
        const u64* q = Q + (i / 192) * 8ul;
        const unsigned at = reinterpret_cast<const unsigned*>(q)[c];
        if (o >= 96) return at + Count(q + 5, 0, o - 96, c);
        return at - Count(q + 2, o, 96, c);
    }
    PF void prefetch(unsigned row) const {
        const unsigned i = row - (row > inv);
        __builtin_prefetch(Q + (i / 192) * 8ul);
    }
};

// Block: q[0..2] = bases 0..95; q[3] = Occ at character 96, less the major count
// (4 x u16); q[4..7] = bases 96..223. Major count every 256 blocks.
struct LayC {
    const u64* Q;
    const unsigned* major;
    unsigned inv;
    unsigned ch(unsigned row) const {
        const unsigned i = row - (row > inv), o = i % 224, w = o / 32;
        return Q[(i / 224) * 8ul + w + (w >= 3)] >> (2 * (o % 32)) & 3;
    }
    unsigned occ(unsigned row, unsigned c) const {
        const unsigned i = row - (row > inv), j = i / 224, o = i - j * 224;
        const u64* q = Q + j * 8ul;
        const unsigned at = major[(j >> 8) * 4 + c] + reinterpret_cast<const uint16_t*>(q + 3)[c];
        if (o >= 96) return at + Count(q + 4, 0, o - 96, c);
        return at - Count(q, o, 96, c);
    }
    PF void prefetch(unsigned row) const {
        const unsigned i = row - (row > inv);
        __builtin_prefetch(Q + (i / 224) * 8ul);
    }
};

static inline u64 Mix(u64& s) {
    u64 z = (s += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

// `chains` LF walks stepped round-robin, each row prefetched a round ahead.
template <class V>
__attribute__((noinline)) u64 Walk(const V& v, const unsigned* C, unsigned n, unsigned chains, unsigned long steps, u64 seed) {
    std::vector<unsigned> row(chains);
    u64 s = seed;
    for (unsigned& r : row) r = 1 + Mix(s) % n;
    for (unsigned long round = steps / chains; round--;) {
        for (unsigned k = 0; k < chains; k++) {
            unsigned r = row[k];
            if (r == v.inv) r = 1 + Mix(s) % n;  // the $ row: start another walk
            const unsigned c = v.ch(r);
            r = C[c] + v.occ(r, c) + 1;
            row[k] = r;
            if (chains > 1) v.prefetch(r);
        }
    }
    u64 sum = 0;
    for (unsigned r : row) sum += r;
    return sum;
}

template <class V>
void Time(const char* name, const V& v, const unsigned* C, unsigned n, unsigned chains, unsigned threads, unsigned long steps) {
    std::vector<u64> sums(threads);
    double best = 1e30;
    for (int rep = 0; rep < 3; rep++) {
        const auto t0 = std::chrono::steady_clock::now();
        std::vector<std::thread> T;
        for (unsigned t = 0; t < threads; t++)
            T.emplace_back([&, t] { sums[t] = Walk(v, C, n, chains, steps, 1234 + t); });
        for (auto& t : T) t.join();
        best = std::min(best, std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
    }
    u64 sum = 0;
    for (u64 x : sums) sum += x;
    const unsigned long done = steps / chains * chains;
    std::printf("%-4s threads=%-2u chains=%-2u  %7.1f ns/step/thread  %8.1f Msteps/s  sum=%llu\n", name, threads, chains,
                best * 1e9 / done, done * threads / best / 1e6, static_cast<unsigned long long>(sum));
    std::fflush(stdout);
}

static std::vector<unsigned> List(const char* s) {
    std::vector<unsigned> v;
    for (const char* p = s; *p;) {
        v.push_back(std::strtoul(p, const_cast<char**>(&p), 10));
        if (*p == ',') p++;
    }
    return v;
}

int main(int argc, char** argv) {
    if (argc < 2) { std::fprintf(stderr, "usage: interleaved_occ <prefix> [--steps N] [--thp] [--threads a,b] [--chains a,b]\n"); return 1; }
    const std::string prefix = argv[1];
    unsigned long steps = 20000000;
    bool thp = false;
    std::vector<unsigned> threads = {1}, chains = {1, 2, 4, 8, 16, 32};
    for (int i = 2; i < argc; i++) {
        const std::string a = argv[i];
        if (a == "--thp") thp = true;
        else if (a == "--steps") steps = std::strtoul(argv[++i], nullptr, 10);
        else if (a == "--threads") threads = List(argv[++i]);
        else if (a == "--chains") chains = List(argv[++i]);
    }

    MMMasterInitialize(3, 0, FALSE, nullptr);
    MMPool* pool = MMPoolCreate(524288);
    BWT* b = BWTLoad(pool, (prefix + ".bwt").c_str(), (prefix + ".fmv").c_str(), nullptr, nullptr, nullptr, nullptr);
    const unsigned n = b->textLength, inv = b->inverseSa0;
    // The library's arrays are moved to the allocator the other layouts use, so that
    // all of them get the same pages (--thp or not).
    auto Rehome = [&](unsigned*& p, size_t words) {
        unsigned* q = static_cast<unsigned*>(Alloc(words * 4 + 64, thp));
        std::memcpy(q, p, words * 4);
        p = q;
    };
    Rehome(b->bwtCode, b->bwtSizeInWord);
    Rehome(b->occValue, b->occSizeInWord);
    Rehome(b->occValueMajor, b->occMajorSizeInWord);
    std::printf("n=%u inverseSa0=%u thp=%d bwtCode%%64=%lu occValue%%64=%lu\n", n, inv, thp,
                reinterpret_cast<uintptr_t>(b->bwtCode) % 64, reinterpret_cast<uintptr_t>(b->occValue) % 64);

    // Build the layouts in one pass over the stored characters.
    const size_t wWords = n / 32 + 32, nCheck = n / 256 + 4, nMajS = n / 65536 + 2;
    const size_t nB = n / 192 + 4, nC = n / 224 + 4, nMajC = nC / 256 + 2;
    u64* Wbase = static_cast<u64*>(Alloc(wWords * 8 + 64, thp));
    u64* W = Wbase + 8;
    uint16_t* minor = static_cast<uint16_t*>(Alloc(nCheck * 8, thp));
    unsigned* majS = static_cast<unsigned*>(Alloc(nMajS * 16, thp));
    u64* QB = static_cast<u64*>(Alloc(nB * 64, thp));
    u64* QC = static_cast<u64*>(Alloc(nC * 64, thp));
    unsigned* majC = static_cast<unsigned*>(Alloc(nMajC * 16, thp));
    unsigned cnt[4] = {0, 0, 0, 0};
    for (unsigned long i = 0; i < n + 600ul; i++) {
        if (i % 65536 == 0 && i / 65536 < nMajS) std::memcpy(majS + i / 65536 * 4, cnt, 16);
        if (i % 256 == 0 && i / 256 < nCheck)
            for (int c = 0; c < 4; c++) minor[i / 256 * 4 + c] = cnt[c] - majS[(i >> 16) * 4 + c];
        if (i % 192 == 96 && i / 192 < nB) std::memcpy(QB + i / 192 * 8, cnt, 16);
        if (i % (224 * 256) == 0 && i / (224 * 256) < nMajC) std::memcpy(majC + i / (224 * 256) * 4, cnt, 16);
        if (i % 224 == 96 && i / 224 < nC) {
            uint16_t* m = reinterpret_cast<uint16_t*>(QC + i / 224 * 8 + 3);
            for (int c = 0; c < 4; c++) m[c] = cnt[c] - majC[(i / 224 >> 8) * 4 + c];
        }
        if (i >= n) continue;
        const u64 c = b->bwtCode[i / 16] << ((i % 16) * 2) >> 30;
        cnt[c]++;
        W[i / 32] |= c << (2 * (i % 32));
        unsigned o = i % 192;
        QB[i / 192 * 8 + 2 + o / 32] |= c << (2 * (o % 32));
        o = i % 224;
        QC[i / 224 * 8 + o / 32 + (o / 32 >= 3)] |= c << (2 * (o % 32));
    }
    std::printf("bytes: bwt %lu + minor %lu + major %lu | B %lu | C %lu + major %lu\n", n / 4ul, (n / 256ul + 1) * 8,
                (n / 65536ul + 1) * 16, (n / 192ul + 1) * 64, (n / 224ul + 1) * 64, (n / (224ul * 256) + 1) * 16);

    const Lib lib{b, inv};
    const Sep sep{W, minor, majS, inv};
    const LayB lb{QB, inv};
    const LayC lc{QC, majC, inv};

    // Occ and characters against the library: random rows, and the rows around $, the
    // ends and every kind of block edge.
    {
        std::vector<unsigned> rows;
        u64 s = 7;
        for (int i = 0; i < 20000000; i++) rows.push_back(Mix(s) % (n + 1ul));
        for (long d = -300; d <= 300; d++) {
            for (unsigned long base : {0ul, static_cast<unsigned long>(inv), static_cast<unsigned long>(n), 65536ul * 100, 224ul * 256 * 77, 192ul * 1000}) {
                const long r = static_cast<long>(base) + d;
                if (r >= 0 && r <= static_cast<long>(n)) rows.push_back(r);
            }
        }
        unsigned long bad = 0;
        for (unsigned r : rows) {
            for (unsigned c = 0; c < 4; c++) {
                const unsigned want = lib.occ(r, c);
                bad += sep.occ(r, c) != want;
                bad += lb.occ(r, c) != want;
                bad += lc.occ(r, c) != want;
            }
            if (r != inv && r >= 1 && r <= n) {
                const unsigned want = lib.ch(r);
                bad += (sep.ch(r) != want) + (lb.ch(r) != want) + (lc.ch(r) != want);
            }
        }
        std::printf("check: %zu rows x 4 bases x 3 layouts, %lu differences\n", rows.size(), bad);
        if (bad) return 1;
    }

    if (FILE* f = std::fopen("/proc/self/smaps_rollup", "r")) {
        char line[256];
        while (std::fgets(line, sizeof line, f))
            if (!std::strncmp(line, "AnonHugePages", 13) || !std::strncmp(line, "Rss", 3)) std::fputs(line, stdout);
        std::fclose(f);
    }
    const unsigned* C = b->cumulativeFreq;
    for (unsigned t : threads) {
        for (unsigned ch : chains) {
            Time("LIB0", Lib0{b, inv}, C, n, ch, t, steps);
            Time("LIB", lib, C, n, ch, t, steps);
            Time("SEP", sep, C, n, ch, t, steps);
            Time("B", lb, C, n, ch, t, steps);
            Time("C", lc, C, n, ch, t, steps);
        }
    }
    return 0;
}
