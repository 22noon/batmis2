// locate_steps: how many LF steps does it take to locate an SA row in a BatMis index?
// Stand-alone; see columba_comparison.md.
//
// BatMis samples SA *rows*: saValue[row / 8] for rows divisible by 8. Locating a row walks
// LF until it reaches such a row, so the walk is geometric (about 8 steps on average, with
// no upper bound). An index that samples *text positions* every s (Columba, s = 4) needs at
// most s - 1 steps, (s - 1) / 2 on average for a uniform position.
//
// For random rows of the index this measures the walk to a sampled row, its distribution,
// and its time (one walk at a time; and `chains` walks interleaved with prefetch, as the
// mappers' batches do).
//
//   locate_steps <prefix> [rows] [chains]
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

extern "C" {
#include "MemManager.h"
#include "BWT.h"
}

static inline uint64_t Mix(uint64_t& s) {
    uint64_t z = (s += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

struct Index {
    const BWT* b;
    unsigned ch(unsigned row) const {
        const unsigned i = row - (row > b->inverseSa0);
        return b->bwtCode[i / 16] << ((i % 16) * 2) >> 30;
    }
    // One LF step; the $ row goes to row 0 (as Sa_Decode does).
    unsigned lf(unsigned row) const {
        if (row == b->inverseSa0) return 0;
        const unsigned c = ch(row);
        return b->cumulativeFreq[c] + BWTOccValue(b, row, c) + 1;
    }
    void prefetch(unsigned row) const {
        const unsigned i = row - (row > b->inverseSa0);
        const unsigned long e = (i + OCC_INTERVAL / 2 - 1) / OCC_INTERVAL;
        __builtin_prefetch(b->occValue + e / OCC_VALUE_PER_WORD * ALPHABET_SIZE);
        __builtin_prefetch(b->bwtCode + i / CHAR_PER_WORD);
    }
};

int main(int argc, char** argv) {
    if (argc < 2) {
        std::fprintf(stderr, "usage: locate_steps <prefix> [rows] [chains]\n");
        return 1;
    }
    const std::string prefix = argv[1];
    const unsigned long Rows = argc > 2 ? std::strtoul(argv[2], nullptr, 10) : 1000000;
    const unsigned Chains = argc > 3 ? std::strtoul(argv[3], nullptr, 10) : 16;
    MMMasterInitialize(3, 0, FALSE, nullptr);
    MMPool* pool = MMPoolCreate(524288);
    const BWT* b = BWTLoad(pool, (prefix + ".bwt").c_str(), (prefix + ".fmv").c_str(), nullptr, nullptr, nullptr, nullptr);
    const Index X{b};
    const unsigned n = b->textLength, interval = 8;

    std::vector<unsigned> start(Rows);
    uint64_t s = 42;
    for (auto& r : start) r = 1 + Mix(s) % n;

    // Distribution, one walk at a time (also its time).
    std::vector<unsigned long> hist(64, 0);
    unsigned long total = 0, worst = 0;
    auto t0 = std::chrono::steady_clock::now();
    for (unsigned r : start) {
        unsigned steps = 0;
        while (r % interval) {
            r = X.lf(r);
            steps++;
        }
        hist[std::min(steps, 63u)]++;
        total += steps;
        worst = std::max<unsigned long>(worst, steps);
    }
    const double one = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();

    // The same walks, Chains at a time, round-robin with each next row prefetched.
    t0 = std::chrono::steady_clock::now();
    unsigned long check = 0;
    for (unsigned long base = 0; base < Rows; base += Chains) {
        const unsigned m = static_cast<unsigned>(std::min<unsigned long>(Chains, Rows - base));
        std::vector<unsigned> row(start.begin() + base, start.begin() + base + m);
        for (unsigned live = m; live;) {
            live = 0;
            for (unsigned k = 0; k < m; k++) {
                if (row[k] % interval == 0) continue;
                row[k] = X.lf(row[k]);
                X.prefetch(row[k]);
                live++;
            }
        }
        for (unsigned r : row) check += r;
    }
    const double many = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();

    std::printf("index %s: %u rows, SA row-sampled every %u\n", prefix.c_str(), n, interval);
    std::printf("LF steps to a sampled row, %lu random rows: mean %.2f, worst %lu\n", Rows, double(total) / Rows, worst);
    unsigned long cum = 0;
    for (unsigned k : {0u, 1u, 2u, 3u, 4u, 8u, 16u, 24u, 32u}) {
        cum = 0;
        for (unsigned i = 0; i <= k; i++) cum += hist[i];
        std::printf("  <= %2u steps: %6.2f%%\n", k, 100.0 * cum / Rows);
    }
    std::printf("time per locate walk: %.0f ns one at a time, %.0f ns with %u interleaved (check %lu)\n",
                one * 1e9 / Rows, many * 1e9 / Rows, Chains, check % 1000);
    std::printf("text-sampled every 4 (Columba): at most 3 steps, 1.5 on average\n");
    return 0;
}
