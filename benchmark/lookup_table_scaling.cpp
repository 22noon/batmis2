// Start-of-search lookup table: cost and benefit by depth k (4^k stored ranges),
// pre-built at startup versus filled on first use. Results: benchmark/lookup_table.md.
// Build from src/ (pulls in batman_zero_min.cpp for Get_SARange and the index loader):
//   g++ -std=c++17 -fpermissive -w -I. -I../index -O3 -funroll-loops -mpopcnt \
//       ../benchmark/lookup_table_scaling.cpp ../index/libfmindex_core.a -lm -o lookup_table_scaling
// Run:  taskset -c 2 ./lookup_table_scaling <genome prefix used with -g> <reads.fq> [max_k (13)] [reads (2000000)]
#define main orig_main
#include "batman_zero_min.cpp"
#undef main
#include <chrono>
#include <sys/resource.h>

// Per depth k, over the first bases of real reads (revfmi, as the forward table):
//   build  : eager recursion over all 4^k prefixes, one Get_SARange per node, into an
//            interleaved {Start, End} table; "lo" keys the first base in the low bits
//            (current layout), "hi" in the high bits (a subtree's leaves are contiguous).
//   lookup : one table read per read, independent / dependent (next key waits for the
//            last result, as a search does).
//   walk   : no table, k Get_SARange steps from the first base's range.
//   t6walk : what the search does today past k = 6: the k=6 table, then k-6 steps.
//   lazy   : table filled on first use (k=6 table + walk on a miss), zero pages from
//            calloc, so untouched parts cost nothing; per-read cost over N reads.

namespace {

using Clock = std::chrono::steady_clock;
struct Entry { unsigned Start, End; };  // {0,0} = unfilled (lazy); no match is stored as {0,1}

double Seconds(Clock::time_point t) { return std::chrono::duration<double>(Clock::now() - t).count(); }

std::vector<Entry> table;
int depth;
bool highKey;

void Build(unsigned long start, unsigned long end, int level, unsigned label) {
    if (level == depth) {
        table[label] = start ? Entry{static_cast<unsigned>(start), static_cast<unsigned>(end)} : Entry{0, 1};
        return;
    }
    for (int b = 0; b < 4; b++) {
        unsigned long s = start, e = end;
        if (s) Get_SARange(b, s, e, revfmi);
        Build(s, e, level + 1, highKey ? (label << 2) | b : label | (b << (2 * level)));
    }
}

void Build_Table(int k, bool high) {
    depth = k;
    highKey = high;
    table.assign(1ul << (2 * k), Entry{0, 0});
    for (unsigned b = 0; b < 4; b++) {
        const Range r = InitialBaseRange(revfmi, b);
        Build(r.Start, r.End, 1, b);
    }
}

Entry Walk(const char* read, int from, unsigned long start, unsigned long end, int k) {
    for (int i = from; i < k && start; i++) Get_SARange(read[i], start, end, revfmi);
    return start ? Entry{static_cast<unsigned>(start), static_cast<unsigned>(end)} : Entry{0, 1};
}

Entry Walk_From_Root(const char* read, int k) {
    const Range r = InitialBaseRange(revfmi, read[0]);
    return Walk(read, 1, r.Start, r.End, k);
}

unsigned Key(const char* read, int k) {
    unsigned c = 0;
    for (int i = 0; i < k; i++) c |= read[i] << (2 * i);
    return c;
}

long Peak_RSS_MB() {
    rusage u{};
    getrusage(RUSAGE_SELF, &u);
    return u.ru_maxrss / 1024;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 3) {
        std::fprintf(stderr, "Usage: %s <genome prefix> <reads.fq> [max_k] [reads]\n", argv[0]);
        return 2;
    }
    const int maxK = argc > 3 ? std::atoi(argv[3]) : 13;
    const unsigned long maxReads = argc > 4 ? std::strtoul(argv[4], nullptr, 10) : 2000000;
    constexpr int W = 16;  // bases kept per read
    std::string pre = argv[1], rb, ro, fb, fo;
    BuildIndexPaths(pre, rb, ro, fb, fo);
    revfmi = initFMI(rb.c_str(), ro.c_str());

    std::vector<char> reads;  // W codes per read
    {
        std::ifstream in(argv[2]);
        std::string id, seq;
        while (reads.size() / W < maxReads && ReadNextRecord(in, id, seq)) {
            if (static_cast<int>(seq.size()) < W) continue;
            bool ok = true;
            char r[W];
            for (int i = 0; i < W && ok; i++) ok = (r[i] = static_cast<char>(EncodeBase(seq[i]))) >= 0;
            if (ok) reads.insert(reads.end(), r, r + W);
        }
    }
    const unsigned long N = reads.size() / W;
    const char* R = reads.data();
    std::printf("text=%lu reads=%lu (%s)\n", revfmi->textLength, N, argv[2]);

    // The k=6 table that today's search starts from; the lazy table and t6walk use it.
    Build_Table(6, false);
    const std::vector<Entry> table6 = table;

    std::printf("%3s %9s %9s %9s %9s | %8s %8s %8s %8s | %8s %8s %8s | %9s %6s\n", "k", "entries", "MB",
                "build_lo", "build_hi", "lk_ind", "lk_dep", "walk", "t6walk", "lazy100k", "lazy1M", "lazyAll",
                "distinct", "errors");
    unsigned long sink = 0;
    for (int k = 3; k <= maxK; k++) {
        const unsigned long entries = 1ul << (2 * k);
        auto t = Clock::now();
        Build_Table(k, true);
        const double buildHi = Seconds(t);
        t = Clock::now();
        Build_Table(k, false);
        const double buildLo = Seconds(t);

        std::vector<unsigned> keys(N);
        for (unsigned long i = 0; i < N; i++) keys[i] = Key(R + i * W, k);

        // Correctness: table = walk from the root = k=6 table + walk.
        unsigned long errors = 0;
        for (unsigned long i = 0; i < N; i++) {
            const char* r = R + i * W;
            const Entry a = table[keys[i]], b = Walk_From_Root(r, k);
            Entry c;
            if (k >= 6) {
                const Entry s = table6[Key(r, 6)];
                c = s.Start ? Walk(r, 6, s.Start, s.End, k) : Entry{0, 1};
            } else {
                c = b;
            }
            const bool same = (a.Start == 0 && b.Start == 0 && c.Start == 0) ||
                              (a.Start == b.Start && a.End == b.End && a.Start == c.Start && a.End == c.End);
            errors += !same;
        }

        double lkInd = 1e99, lkDep = 1e99, walk = 1e99, t6walk = 1e99;
        for (int rep = 0; rep < 3; rep++) {
            t = Clock::now();
            for (unsigned long i = 0; i < N; i++) sink += table[keys[i]].Start;
            lkInd = std::min(lkInd, Seconds(t));
            t = Clock::now();
            unsigned prev = 0;
            for (unsigned long i = 0; i < N; i++) {
                unsigned zero = prev;
                asm("" : "+r"(zero));  // opaque: Start < 2^31, so this is 0, but the load must wait for prev
                prev = table[keys[i] ^ (zero >> 31)].Start;
            }
            sink += prev;
            lkDep = std::min(lkDep, Seconds(t));
            t = Clock::now();
            for (unsigned long i = 0; i < N; i++) sink += Walk_From_Root(R + i * W, k).Start;
            walk = std::min(walk, Seconds(t));
            if (k > 6) {
                t = Clock::now();
                for (unsigned long i = 0; i < N; i++) {
                    const char* r = R + i * W;
                    const Entry s = table6[Key(r, 6)];
                    sink += s.Start ? Walk(r, 6, s.Start, s.End, k).Start : 0;
                }
                t6walk = std::min(t6walk, Seconds(t));
            }
        }

        // Lazy: a fresh zeroed table per run (calloc'd pages are mapped on first touch).
        double lazy[3] = {0, 0, 0};
        unsigned long distinct = 0;
        const unsigned long lazyN[3] = {std::min(N, 100000ul), std::min(N, 1000000ul), N};
        for (int j = 0; j < 3; j++) {
            t = Clock::now();
            Entry* lt = static_cast<Entry*>(std::calloc(entries, sizeof(Entry)));
            unsigned long filled = 0;
            for (unsigned long i = 0; i < lazyN[j]; i++) {
                Entry& e = lt[keys[i]];
                if (!e.Start && !e.End) {
                    const char* r = R + i * W;
                    if (k > 6) {
                        const Entry s = table6[Key(r, 6)];
                        e = s.Start ? Walk(r, 6, s.Start, s.End, k) : Entry{0, 1};
                    } else {
                        e = Walk_From_Root(r, k);
                    }
                    filled++;
                }
                sink += e.Start;
            }
            lazy[j] = Seconds(t);
            std::free(lt);
            if (j == 2) distinct = filled;
        }

        const double ns = 1e9 / N;
        std::printf("%3d %9lu %9.1f %9.3f %9.3f | %8.1f %8.1f %8.1f %8.1f | %8.1f %8.1f %8.1f | %8.1f%% %6lu\n", k,
                    entries, entries * sizeof(Entry) / 1048576.0, buildLo, buildHi, lkInd * ns, lkDep * ns, walk * ns,
                    k > 6 ? t6walk * ns : 0.0, lazy[0] * 1e9 / lazyN[0], lazy[1] * 1e9 / lazyN[1], lazy[2] * ns,
                    100.0 * distinct / N, errors);
        std::fflush(stdout);
    }
    std::printf("(ns per read for lookup/walk/lazy columns; build in s; peak RSS %ld MB; sink %lu)\n", Peak_RSS_MB(),
                sink & 1);
    return 0;
}
