// Lossy dictionary of the genome's unique X-mers: does a confirmed hit make the start
// of a search cheaper than the lookup table + walk? Results: benchmark/unique_kmer_dict.md.
// Build from src/ (pulls in batman_zero_min.cpp for Get_SARange and the index loader):
//   g++ -std=c++17 -fpermissive -w -I. -I../index -O3 -funroll-loops -mpopcnt \
//       ../benchmark/unique_kmer_dict.cpp ../index/libfmindex_core.a -lm -o unique_kmer_dict
// Run:  taskset -c 2 ./unique_kmer_dict <genome prefix> <table k> <reads.fq>...
//       ./unique_kmer_dict <genome prefix> dump <X>   (unique X-mers as text, for checking)
//       ./unique_kmer_dict <genome prefix> build -    (only build both ways and compare)
//
// Unique X-mers: a depth-first walk over revfmi (the forward search's orientation)
// to depth 16. A range of one row at depth X (12, 14, 16) is a unique X-mer; its key
// is the X-mer with the first base in the low bits, as the lookup tables use.
// Single-row branches have one child, read from the BWT directly.
//
// Dictionary: 64-byte buckets of 8 {key, row} slots, bucket from a hash of the key,
// key stored whole, so a hit is exact (no false positives). A key whose bucket is
// full is dropped (the loss). Row 0 ($) marks an empty slot.
//
// Per read, from the read's first X bases (and the depth-k table for the first k):
//   base  : table, then Get_SARange to depth X                        (today)
//   A     : dictionary; on a hit, table then Start-only rank steps     (presence only)
//   B     : dictionary; on a hit, the stored row                       (presence + row)
//   misses in A and B fall back to base.
#define main orig_main
#include "batman_zero_min.cpp"
#undef main
#include <chrono>

namespace {

using Clock = std::chrono::steady_clock;
double Seconds(Clock::time_point t) { return std::chrono::duration<double>(Clock::now() - t).count(); }

constexpr int Xs[3] = {12, 14, 16};
constexpr int MAXX = 16;

struct Pair {
    unsigned Key;
    unsigned Row;
};
std::vector<Pair> uniques[3];

// Row's preceding character (its child's base), or -1 at the $ row.
inline int Bwt_Char(unsigned long row) {
    if (row == revfmi->inverseSa0) return -1;
    const unsigned long i = row - (row > revfmi->inverseSa0);
    return revfmi->bwtCode[i / 16] << ((i % 16) * 2) >> (BITS_IN_WORD - 2);
}

inline void Record(int level, unsigned label, unsigned long row) {
    for (int x = 0; x < 3; x++)
        if (level == Xs[x]) uniques[x].push_back({label, static_cast<unsigned>(row)});
}

// level = bases matched; label holds them, first base in the low bits.
void Dfs(unsigned long start, unsigned long end, int level, unsigned label) {
    if (start == end) {  // one row: one path down
        for (;;) {
            Record(level, label, start);
            if (level == MAXX) return;
            const int c = Bwt_Char(start);
            if (c < 0) return;
            start = revfmi->cumulativeFreq[c] + BWTOccValue(revfmi, start, c) + 1;
            label |= static_cast<unsigned>(c) << (2 * level);
            level++;
        }
    }
    if (level == MAXX) return;
    for (int b = 0; b < 4; b++) {
        unsigned long s = start, e = end;
        Get_SARange(b, s, e, revfmi);
        if (s) Dfs(s, e, level + 1, label | static_cast<unsigned>(b) << (2 * level));
    }
}

// Unique X-mers by sorting instead: every X-mer of the forward genome (.pac) as
// key << 32 | position, LSD radix sort on the key (11+11+10 bits), keep keys seen
// once. Gives the genome position, not the SA row. Positions must fit in 32 bits.
std::vector<Pair> Unique_By_Sort(int x, unsigned long n, double& sortSeconds) {
    std::vector<uint64_t> a(n - x + 1), b(a.size());
    auto t = Clock::now();
    unsigned key = 0;
    for (unsigned long p = 0; p < n; p++) {
        const unsigned g = (Pac[p / 4] >> (6 - 2 * (p % 4))) & 3;
        key = (key >> 2) | (g << (2 * (x - 1)));  // first base in the low bits
        if (x < 16) key &= (1u << (2 * x)) - 1;
        if (p + 1 >= static_cast<unsigned long>(x)) a[p + 1 - x] = static_cast<uint64_t>(key) << 32 | (p + 1 - x);
    }
    const int widths[3] = {11, 11, 10};
    int shift = 32;
    for (int w : widths) {
        std::vector<unsigned long> count((1u << w) + 1, 0);
        const uint64_t mask = (1u << w) - 1;
        for (uint64_t v : a) count[((v >> shift) & mask) + 1]++;
        for (std::size_t i = 1; i < count.size(); i++) count[i] += count[i - 1];
        for (uint64_t v : a) b[count[(v >> shift) & mask]++] = v;
        a.swap(b);
        shift += w;
    }
    std::vector<Pair> out;
    for (std::size_t i = 0; i < a.size();) {
        std::size_t j = i + 1;
        while (j < a.size() && (a[j] >> 32) == (a[i] >> 32)) j++;
        if (j == i + 1) out.push_back({static_cast<unsigned>(a[i] >> 32), static_cast<unsigned>(a[i])});
        i = j;
    }
    sortSeconds = Seconds(t);
    return out;
}

// Lookup table of depth k (as in lookup_table_scaling.cpp).
struct Entry {
    unsigned Start, End;
};
std::vector<Entry> table;
int tableK;
void Build(unsigned long s, unsigned long e, int level, unsigned label) {
    if (level == tableK) {
        table[label] = s ? Entry{static_cast<unsigned>(s), static_cast<unsigned>(e)} : Entry{0, 0};
        return;
    }
    for (int b = 0; b < 4; b++) {
        unsigned long s2 = s, e2 = e;
        if (s2) Get_SARange(b, s2, e2, revfmi);
        Build(s2, e2, level + 1, label | b << (2 * level));
    }
}

inline unsigned Hash(unsigned key) {  // murmur3 finalizer
    key ^= key >> 16; key *= 0x85ebca6bu; key ^= key >> 13; key *= 0xc2b2ae35u; key ^= key >> 16;
    return key;
}

struct Dict {
    struct alignas(64) Bucket { Pair Slot[8]; };
    std::vector<Bucket> buckets;
    unsigned long dropped = 0;
    void Init(unsigned long entries, double load) {
        buckets.assign(std::max(1ul, static_cast<unsigned long>(entries / load / 8)), Bucket{});
        dropped = 0;
    }
    unsigned long Bucket_Of(unsigned key) const { return (static_cast<unsigned long>(Hash(key)) * buckets.size()) >> 32; }
    void Insert(Pair p) {
        Bucket& b = buckets[Bucket_Of(p.Key)];
        for (Pair& s : b.Slot)
            if (!s.Row) { s = p; return; }
        dropped++;
    }
    unsigned Find(unsigned key) const {  // row, or 0
        const Bucket& b = buckets[Bucket_Of(key)];
        for (const Pair& s : b.Slot)
            if (s.Row && s.Key == key) return s.Row;
        return 0;
    }
    double MB() const { return buckets.size() * sizeof(Bucket) / 1048576.0; }
};

unsigned Key(const char* r, int n) {
    unsigned c = 0;
    for (int i = 0; i < n; i++) c |= static_cast<unsigned>(r[i]) << (2 * i);
    return c;
}

// Today: table, then Get_SARange to depth x. Returns Start (0 = absent); End in end.
inline unsigned long Walk(const char* r, int x, unsigned long& end) {
    const Entry t = table[Key(r, tableK)];
    unsigned long s = t.Start;
    end = t.End;
    for (int i = tableK; i < x && s; i++) Get_SARange(r[i], s, end, revfmi);
    return s;
}

// Start-only: valid only when r[0..x) is known to occur.
inline unsigned long Start_Only(const char* r, int x) {
    unsigned long s = table[Key(r, tableK)].Start;
    for (int i = tableK; i < x; i++) s = revfmi->cumulativeFreq[r[i]] + BWTOccValue(revfmi, s, r[i]) + 1;
    return s;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 4) {
        std::fprintf(stderr, "Usage: %s <genome prefix> <table k> <reads.fq>...\n       %s <genome prefix> dump <X>\n", argv[0], argv[0]);
        return 2;
    }
    std::string pre = argv[1], rb, ro, fb, fo;
    BuildIndexPaths(pre, rb, ro, fb, fo);
    revfmi = initFMI(rb.c_str(), ro.c_str());

    for (int i = 0; i < 3; i++)  // upper bounds, so the vectors never reallocate
        uniques[i].reserve(std::min(1ul << (2 * Xs[i]), static_cast<unsigned long>(revfmi->textLength)));
    auto t = Clock::now();
    for (int b = 0; b < 4; b++) {
        const Range r = InitialBaseRange(revfmi, b);
        if (r.Start) Dfs(r.Start, r.End, 1, b);
    }
    const double dfsSeconds = Seconds(t);

    // The same sets by sorting the .pac X-mers; keys must match the index walk exactly.
    double sortSeconds[3] = {0, 0, 0};
    unsigned long sortMismatch[3] = {0, 0, 0};
    if (LoadPac(pre + ".pac")) {
        for (int i = 0; i < 3; i++) {
            const std::vector<Pair> sorted = Unique_By_Sort(Xs[i], revfmi->textLength, sortSeconds[i]);
            std::vector<unsigned> walked(uniques[i].size());
            for (std::size_t j = 0; j < walked.size(); j++) walked[j] = uniques[i][j].Key;
            std::sort(walked.begin(), walked.end());
            sortMismatch[i] = sorted.size() != walked.size();
            for (std::size_t j = 0; !sortMismatch[i] && j < walked.size(); j++) sortMismatch[i] += sorted[j].Key != walked[j];
            // Each position really starts its X-mer.
            for (std::size_t j = 0; j < sorted.size(); j += 997) {
                unsigned k = 0;
                for (int q = 0; q < Xs[i]; q++) k |= ((Pac[(sorted[j].Row + q) / 4] >> (6 - 2 * ((sorted[j].Row + q) % 4))) & 3) << (2 * q);
                sortMismatch[i] += k != sorted[j].Key;
            }
        }
        std::vector<uint8_t>().swap(Pac);
    }

    if (std::strcmp(argv[2], "dump") == 0) {
        const int x = std::atoi(argv[3]);
        for (int i = 0; i < 3; i++) {
            if (Xs[i] != x) continue;
            for (const Pair& p : uniques[i]) {
                char s[MAXX + 1];
                for (int j = 0; j < x; j++) s[j] = "ACGT"[(p.Key >> (2 * j)) & 3];
                s[x] = 0;
                std::printf("%s\n", s);
            }
        }
        return 0;
    }

    std::printf("text=%lu  unique-X-mer walk to depth %d: %.1f s\n", revfmi->textLength, MAXX, dfsSeconds);
    for (int i = 0; i < 3; i++)
        std::printf("  X=%d  unique=%lu (%.1f%% of 4^X, %.1f%% of positions); by sorting the .pac: %.1f s, %s\n", Xs[i],
                    uniques[i].size(), 100.0 * uniques[i].size() / (1ul << (2 * Xs[i])),
                    100.0 * uniques[i].size() / revfmi->textLength, sortSeconds[i],
                    sortSeconds[i] == 0 ? "no .pac" : sortMismatch[i] ? "MISMATCH" : "same keys");
    if (std::strcmp(argv[2], "build") == 0) return 0;

    tableK = std::atoi(argv[2]);
    table.assign(1ul << (2 * tableK), Entry{});
    t = Clock::now();
    for (unsigned b = 0; b < 4; b++) {
        const Range r = InitialBaseRange(revfmi, b);
        Build(r.Start, r.End, 1, b);
    }
    std::printf("depth-%d table: %.2f s\n", tableK, Seconds(t));

    // Reads: first MAXX bases of each.
    std::vector<std::vector<char>> sets;
    std::vector<std::string> names;
    for (int a = 3; a < argc; a++) {
        std::ifstream in(argv[a]);
        std::string id, seq;
        std::vector<char> r;
        while (ReadNextRecord(in, id, seq)) {
            if (static_cast<int>(seq.size()) < MAXX) continue;
            char c[MAXX];
            bool ok = true;
            for (int i = 0; i < MAXX && ok; i++) ok = (c[i] = static_cast<char>(EncodeBase(seq[i]))) >= 0;
            if (ok) r.insert(r.end(), c, c + MAXX);
        }
        sets.push_back(std::move(r));
        const char* slash = std::strrchr(argv[a], '/');
        names.push_back(slash ? slash + 1 : argv[a]);
    }

    unsigned long sink = 0;
    const double loads[3] = {0.5, 0.8, 1.0};
    for (int xi = 0; xi < 3; xi++) {
        const int x = Xs[xi];
        if (x <= tableK || uniques[xi].empty()) continue;
        for (double load : loads) {
            Dict d;
            d.Init(uniques[xi].size(), load);
            t = Clock::now();
            for (const Pair& p : uniques[xi]) d.Insert(p);
            const double insertSeconds = Seconds(t);
            // Every unique X-mer not dropped is found, with its row.
            unsigned long found = 0, wrong = 0;
            for (const Pair& p : uniques[xi]) {
                const unsigned row = d.Find(p.Key);
                found += row != 0;
                wrong += row && row != p.Row;
            }
            std::printf("\nX=%d load=%.1f: %.0f MB (%.1f B/unique), insert %.2f s, dropped %lu (%.2f%%), found %lu, wrong %lu, lost %lu\n",
                        x, load, d.MB(), d.MB() * 1048576.0 / uniques[xi].size(), insertSeconds, d.dropped,
                        100.0 * d.dropped / uniques[xi].size(), found, wrong, uniques[xi].size() - d.dropped - found);
            std::printf("  %-14s %8s %7s %7s %7s %7s | %7s %7s %7s %6s | %7s %7s %7s %6s | %6s %6s\n", "reads", "n", "hit%",
                        "repeat%", "absent%", "drop%", "base", "A", "B", "query", "base_d", "A_d", "B_d", "qry_d",
                        "past_X", "errors");
            for (std::size_t si = 0; si < sets.size(); si++) {
                const char* R = sets[si].data();
                const unsigned long n = sets[si].size() / MAXX;
                // Classify and check, untimed.
                unsigned long hit = 0, repeat = 0, absent = 0, drop = 0, errors = 0, pastX = 0, pastN = 0;
                for (unsigned long i = 0; i < n; i++) {
                    const char* r = R + i * MAXX;
                    unsigned long end;
                    const unsigned long s = Walk(r, x, end);
                    const unsigned row = d.Find(Key(r, x));
                    if (row) {
                        hit++;
                        errors += !(s && end == s && s == row && Start_Only(r, x) == row);  // exact, unique, same row
                    } else if (!s) {
                        absent++;
                    } else if (end > s) {
                        repeat++;
                    } else {
                        drop++;
                    }
                    if (s && end > s) {  // steps past X until one row (or the 16 kept bases run out)
                        unsigned long ss = s, ee = end;
                        int lv = x;
                        while (ss && ee > ss && lv < MAXX) Get_SARange(r[lv++], ss, ee, revfmi);
                        pastX += lv - x;
                        pastN++;
                    }
                }
                // ns per read for base, A, B and the query alone; independent reads may
                // overlap in the CPU, dependent ones wait for the previous result
                // (as reads do in the search, where each has much more work).
                double tm[2][4];
                for (int dep = 0; dep < 2; dep++) {
                    for (int m = 0; m < 4; m++) {
                        tm[dep][m] = 1e99;
                        for (int rep = 0; rep < 3; rep++) {
                            t = Clock::now();
                            unsigned long prev = 0;
                            for (unsigned long i = 0; i < n; i++) {
                                unsigned long zero = dep ? prev : 0;
                                asm("" : "+r"(zero));
                                const char* r = R + i * MAXX + (zero >> 40);  // rows < 2^32: adds 0
                                unsigned long e, v;
                                if (m == 0) v = Walk(r, x, e);
                                else if (m == 1) v = d.Find(Key(r, x)) ? Start_Only(r, x) : Walk(r, x, e);
                                else if (m == 2) { const unsigned row = d.Find(Key(r, x)); v = row ? row : Walk(r, x, e); }
                                else v = d.Find(Key(r, x));
                                prev = v;
                                sink += v;
                            }
                            tm[dep][m] = std::min(tm[dep][m], Seconds(t));
                        }
                    }
                }
                const double ns = 1e9 / n;
                std::printf("  %-14s %8lu %6.1f%% %6.1f%% %6.1f%% %6.2f%% | %7.1f %7.1f %7.1f %6.1f | %7.1f %7.1f %7.1f %6.1f | %6.2f %6lu\n",
                            names[si].c_str(), n, 100.0 * hit / n, 100.0 * repeat / n, 100.0 * absent / n,
                            100.0 * drop / n, tm[0][0] * ns, tm[0][1] * ns, tm[0][2] * ns, tm[0][3] * ns,
                            tm[1][0] * ns, tm[1][1] * ns, tm[1][2] * ns, tm[1][3] * ns,
                            pastN ? static_cast<double>(pastX) / pastN : 0.0, errors);
                std::fflush(stdout);
            }
        }
    }
    std::printf("\n(ns per read; past_X = mean steps beyond X to one row for repeated X-mers, capped at %d bases; sink %lu)\n",
                MAXX, sink & 1);
    return 0;
}
