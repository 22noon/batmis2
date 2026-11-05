// Can prefetching overlap the cache misses of independent branches? Stand-alone
// experiment, no change to src/. Results: benchmark/branch_prefetch.md.
// Build from src/ (pulls in batman_zero_min.cpp):
//   g++ -std=c++17 -fpermissive -w -I. -I../index -O3 -funroll-loops -mpopcnt \
//       ../benchmark/branch_prefetch.cpp ../index/libfmindex_core.a -lm -o branch_prefetch
// Run:  taskset -c 2 ./branch_prefetch <genome prefix> <chains.fq> <reads.fq>...
//
// Part 1 (ceiling): B independent chains of Get_SARange steps (a read's depth-11 range,
//   then 8 of its bases), finished one after another, advanced round-robin, round-robin
//   with every chain's rows prefetched first, or "split" (all Starts, then prefetch all
//   End+1 rows, then all Ends).
// Part 2 (real branching): one-mismatch extension of each read's right half from its
//   exact left-half range, Count = 1 (every node branches until its mismatch is used):
//   V0 the real Extend_Right; V0c the same depth-first order with this file's step;
//   V1 V0c + prefetch of each child's rows when it is pushed; V2 level by level per
//   read (prefetch all live nodes, then step them); V3 V2 over R reads at once.
// Use_Ref is off: only FM-index steps are timed.
#define main orig_main
#include "batman_zero_min.cpp"
#undef main
#include <chrono>

namespace {

using Clock = std::chrono::steady_clock;
double Seconds(Clock::time_point t) { return std::chrono::duration<double>(Clock::now() - t).count(); }

std::vector<std::string> Read_Set(const char* path, unsigned long maxReads) {
    std::ifstream in(path);
    std::string id, seq;
    std::vector<std::string> out;
    while (out.size() < maxReads && ReadNextRecord(in, id, seq)) {
        std::string r(seq.size(), 0);
        bool ok = true;
        for (std::size_t i = 0; i < seq.size() && ok; i++) ok = (r[i] = static_cast<char>(EncodeBase(seq[i]))) >= 0;
        if (ok) out.push_back(std::move(r));
    }
    return out;
}

void Prefetch_Range(const SARange& r) {
    Prefetch_Row(revfmi, r.Start);
    if (!r.Skip && r.End != r.Start) Prefetch_Row(revfmi, r.End + 1);
}

// ---------------------------------------------------------------- part 1
struct Chain {
    unsigned long Start, End;
    const char* Bases;
};
constexpr int CHAIN_STEPS = 8;

unsigned long Sequential(std::vector<Chain>& c) {
    unsigned long steps = 0;
    for (Chain& x : c)
        for (int i = 0; i < CHAIN_STEPS && x.Start; i++, steps++) Get_SARange(x.Bases[i], x.Start, x.End, revfmi);
    return steps;
}

template <int MODE>  // 0 round-robin, 1 + prefetch, 2 split
unsigned long Batched(std::vector<Chain>& c, std::size_t B) {
    unsigned long steps = 0;
    for (std::size_t g = 0; g < c.size(); g += B) {
        const std::size_t e = std::min(c.size(), g + B);
        for (int i = 0; i < CHAIN_STEPS; i++) {
            if (MODE == 1) {
                for (std::size_t j = g; j < e; j++)
                    if (c[j].Start) { Prefetch_Row(revfmi, c[j].Start); if (c[j].End != c[j].Start) Prefetch_Row(revfmi, c[j].End + 1); }
            }
            if (MODE == 2) {
                // Starts in bulk (prefetching each chain's End+1 row as it goes), then Ends.
                unsigned long oldEnd[64];
                for (std::size_t j = g; j < e; j++) {
                    Chain& x = c[j];
                    if (!x.Start) continue;
                    oldEnd[j - g] = x.End;
                    Prefetch_Row(revfmi, x.End + 1);
                    const unsigned b = x.Bases[i];
                    x.Start = revfmi->cumulativeFreq[b] + BWTOccValue(revfmi, x.Start, b) + 1;
                }
                for (std::size_t j = g; j < e; j++) {
                    Chain& x = c[j];
                    if (!x.Start) continue;
                    const unsigned b = x.Bases[i];
                    x.End = revfmi->cumulativeFreq[b] + BWTOccValue(revfmi, oldEnd[j - g] + 1, b);
                    if (x.End < x.Start) x.Start = 0;
                    steps++;
                }
                continue;
            }
            for (std::size_t j = g; j < e; j++)
                if (c[j].Start) { Get_SARange(c[j].Bases[i], c[j].Start, c[j].End, revfmi); steps++; }
        }
    }
    return steps;
}

// ---------------------------------------------------------------- part 2
struct Hit {
    unsigned long Start, End;
    unsigned Pos, Char;
    bool operator<(const Hit& o) const { return std::tie(Start, End, Pos, Char) < std::tie(o.Start, o.End, o.Pos, o.Char); }
    bool operator==(const Hit& o) const { return Start == o.Start && End == o.End && Pos == o.Pos && Char == o.Char; }
};
Hit To_Hit(const SARange& r) { return {r.Start, r.Skip ? r.Start : r.End, r.Mismatch_Pos[0], r.Mismatch_Char}; }

struct Node {
    SARange R;
    unsigned Read;  // index into the batch (V3)
};

struct Stats {
    unsigned long steps = 0, branchings = 0;
};

// One step of a node, as Extend_Right / Extend_Right_Singleton do it (Use_Ref off,
// deferrals not recorded). Children go to Out; hits to Hits.
template <class Out>
inline void Step(const SARange& r, std::string_view tag, int offset, int count, Out&& push, std::vector<Hit>& hits, Stats& st) {
    const int len = static_cast<int>(tag.size());
    st.steps++;
    if (r.End == r.Start || r.Skip) {  // single row
        SARange t = r;
        if (!t.Skip && t.Start % SAINTERVAL == 0) { t.Skip = 1; t.End = t.Start; }
        if (t.Start == revfmi->inverseSa0) return;
        const unsigned long idx = t.Start - (t.Start > revfmi->inverseSa0);
        const unsigned now = revfmi->bwtCode[idx / 16] << ((idx % 16) * 2) >> (BITS_IN_WORD - 2);
        if (now == t.Forbidden_Base) return;
        const unsigned char expected = tag[t.Level - 1];
        if (now != expected && t.Mismatches == count) return;  // deferred in the real search
        t.Forbidden_Base = NO_BASE;
        const unsigned base = Read_Base(expected);  // on a match, rank with the read's base (see Read_Base)
        if (now == expected) {
            t.Start = revfmi->cumulativeFreq[base] + BWTOccValue(revfmi, t.Start, base) + 1;
        } else {
            asm volatile("");  // keeps the two calls apart: merged, the rank would wait for the BWT load
            t.Start = revfmi->cumulativeFreq[now] + BWTOccValue(revfmi, t.Start, now) + 1;
        }
        t.End = t.Start;  // still one row (the real loop never looks at End here)
        if (!t.Skip && t.Start % SAINTERVAL == 0) t.Skip = 1;
        if (expected != now) {
            t.Mismatch_Char |= now << t.Mismatches * 2;
            t.Mismatch_Pos[t.Mismatches] = offset + t.Level - 1;
            t.Mismatches++;
        }
        if (t.Level == len) {
            if (t.Mismatches == count) { if (!t.Skip) t.End = t.Start; hits.push_back(To_Hit(t)); }
            return;
        }
        t.Level++;
        push(t);
    } else if (r.Mismatches == count) {  // exact step only
        SARange t = r;
        Get_SARange(tag[t.Level - 1], t, revfmi);
        if (!t.Start) return;
        t.Forbidden_Base = NO_BASE;
        if (t.Level == len) { hits.push_back(To_Hit(t)); return; }
        t.Level++;
        push(t);
    } else {
        st.branchings++;
        Branch_Detect(r, revfmi, offset);
        SARange kids[4];
        unsigned nk = 0;
        for (int b = 0; b < 4; b++) {
            if (!Branch_Characters[b] || b == r.Forbidden_Base) continue;
            SARange t = r;
            t.Forbidden_Base = NO_BASE;
            t.Start = Branch_Ranges[b].Start;
            t.End = Branch_Ranges[b].End;
            if (tag[t.Level - 1] != b) {
                t.Mismatch_Char |= b << t.Mismatches * 2;
                t.Mismatch_Pos[t.Mismatches] = offset + t.Level - 1;
                t.Mismatches++;
            }
            if (t.Mismatches > count) continue;
            if (t.Level == len) {
                if (t.Mismatches == count) hits.push_back(To_Hit(t));
                continue;
            }
            t.Level++;
            kids[nk++] = t;
        }
        for (unsigned k = 0; k < nk; k++) push(kids[k]);
    }
}

}  // namespace

namespace {

struct Work {
    std::string_view Rh;
    SARange Start;  // exact range of the left half, Level 1 in the right half
};

// V0: the real Extend_Right.
void Real(const std::vector<Work>& w, int lh, std::vector<std::vector<Hit>>* hits) {
    ForwardResults out;
    for (std::size_t i = 0; i < w.size(); i++) {
        out.hits.clear(); out.mismatches.clear(); out.hitCount = 0;
        Extend_Right(w[i].Rh, lh, w[i].Start, 1, revfmi, out, 0);
        if (hits) for (const SARange& h : out.hits) (*hits)[i].push_back(To_Hit(h));
    }
}

// V0c / V1: Extend_Right's loop (single-row nodes walked by the real
// Extend_Right_Singleton); PF prefetches each child's rows when it is pushed.
template <bool PF>
void Depth_First(const std::vector<Work>& w, int lh, std::vector<std::vector<Hit>>& hits, Stats& st) {
    std::vector<SARange> stack;
    ForwardResults out;
    for (std::size_t i = 0; i < w.size(); i++) {
        const std::string_view tag = w[i].Rh;
        const int len = static_cast<int>(tag.size());
        out.hits.clear(); out.mismatches.clear(); out.hitCount = 0;
        stack.assign(1, w[i].Start);
        while (!stack.empty()) {
            SARange r = stack.back();
            stack.pop_back();
            st.steps++;
            if (r.End == r.Start || r.Skip) {
                Extend_Right_Singleton(tag, lh, r, 1, revfmi, out);
            } else if (r.Mismatches == 1) {
                SARange t = r;
                Get_SARange(tag[t.Level - 1], t, revfmi);
                if (!t.Start) continue;
                t.Forbidden_Base = NO_BASE;
                if (t.Level == len) { out.hits.push_back(t); continue; }
                t.Level++;
                if (PF) Prefetch_Range(t);
                stack.push_back(t);
            } else {
                st.branchings++;
                Branch_Detect(r, revfmi, lh);
                for (int b = 0; b < 4; b++) {
                    if (!Branch_Characters[b] || b == r.Forbidden_Base) continue;
                    SARange t = r;
                    t.Forbidden_Base = NO_BASE;
                    t.Start = Branch_Ranges[b].Start;
                    t.End = Branch_Ranges[b].End;
                    if (tag[t.Level - 1] != b) {
                        t.Mismatch_Char |= b << t.Mismatches * 2;
                        t.Mismatch_Pos[t.Mismatches] = lh + t.Level - 1;
                        t.Mismatches++;
                    }
                    if (t.Mismatches > 1) continue;
                    if (t.Level == len) { if (t.Mismatches == 1) out.hits.push_back(t); continue; }
                    t.Level++;
                    if (PF) Prefetch_Range(t);
                    stack.push_back(t);
                }
            }
        }
        for (const SARange& h : out.hits) hits[i].push_back(To_Hit(h));
    }
}

// V2 (R = 1) / V3: level by level over R reads; PF prefetches every live node's rows
// before stepping any of them.
template <bool PF>
void Frontier(const std::vector<Work>& w, int lh, std::size_t R, std::vector<std::vector<Hit>>& hits, Stats& st,
              double& widthSum, unsigned long& rounds, unsigned long& widthMax) {
    std::vector<Node> cur, next;
    for (std::size_t g = 0; g < w.size(); g += R) {
        cur.clear();
        for (std::size_t i = g; i < std::min(w.size(), g + R); i++) cur.push_back({w[i].Start, static_cast<unsigned>(i)});
        while (!cur.empty()) {
            widthSum += cur.size(); rounds++; widthMax = std::max<unsigned long>(widthMax, cur.size());
            if (PF) for (const Node& n : cur) Prefetch_Range(n.R);
            next.clear();
            for (const Node& n : cur)
                Step(n.R, w[n.Read].Rh, lh, 1, [&](const SARange& c) { next.push_back({c, n.Read}); }, hits[n.Read], st);
            cur.swap(next);
        }
    }
}

unsigned long Compare(std::vector<std::vector<Hit>>& a, std::vector<std::vector<Hit>>& b) {
    unsigned long bad = 0;
    for (std::size_t i = 0; i < a.size(); i++) {
        std::sort(a[i].begin(), a[i].end());
        std::sort(b[i].begin(), b[i].end());
        bad += !(a[i] == b[i]);
    }
    return bad;
}

void Tables_For(int length) {
    STRINGLENGTH = length;
    LOOKUP_DEPTH = 11;
    Build_Tables(revfmi);
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 4) {
        std::fprintf(stderr, "Usage: %s <genome prefix> <chains.fq> <reads.fq>...\n", argv[0]);
        return 2;
    }
    std::string pre = argv[1], rb, ro, fb, fo;
    BuildIndexPaths(pre, rb, ro, fb, fo);
    revfmi = initFMI(rb.c_str(), ro.c_str());
    Use_Ref = false;
    unsigned long sink = 0;

    // ---- part 1
    {
        const std::vector<std::string> reads = Read_Set(argv[2], 1000000);
        Tables_For(static_cast<int>(reads[0].size()));
        std::vector<Chain> base;
        for (const std::string& r : reads) {
            unsigned key = 0;
            for (int i = 0; i < LOOKUPSIZE; i++) key |= static_cast<unsigned>(r[i]) << (2 * i);
            base.push_back({Forward_Lookup[key].Start, Forward_Lookup[key].End, r.data() + LOOKUPSIZE});
        }
        std::vector<Chain> ref = base;
        const unsigned long steps = Sequential(ref);
        std::printf("Part 1: %zu chains of %d Get_SARange steps from depth %d (%lu steps); ns per step, best of 3\n",
                    base.size(), CHAIN_STEPS, LOOKUPSIZE, steps);
        std::printf("%4s %10s %10s %10s %10s   %s\n", "B", "sequential", "round-rob", "rr+prefch", "split", "check");
        for (std::size_t B : {1, 2, 4, 8, 16, 32}) {
            double t[4] = {1e99, 1e99, 1e99, 1e99};
            unsigned long bad = 0;
            for (int rep = 0; rep < 3; rep++) {
                for (int m = 0; m < 4; m++) {
                    std::vector<Chain> c = base;
                    auto t0 = Clock::now();
                    unsigned long s = m == 0 ? Sequential(c) : m == 1 ? Batched<0>(c, B) : m == 2 ? Batched<1>(c, B) : Batched<2>(c, B);
                    t[m] = std::min(t[m], Seconds(t0));
                    sink += s;
                    bad += s != steps;
                    for (std::size_t j = 0; j < c.size(); j++)
                        bad += c[j].Start != ref[j].Start || (c[j].Start && c[j].End != ref[j].End);
                }
            }
            std::printf("%4zu %10.1f %10.1f %10.1f %10.1f   %s\n", B, t[0] * 1e9 / steps, t[1] * 1e9 / steps,
                        t[2] * 1e9 / steps, t[3] * 1e9 / steps, bad ? "MISMATCH" : "same ranges");
            std::fflush(stdout);
        }
    }

    // ---- part 2
    std::printf("\nPart 2: one-mismatch extension of the right half (Count 1), ns per read, best of 3\n");
    std::printf("%-14s %8s %7s %7s %6s %6s | %8s %8s %8s %8s %8s %8s %8s | %s\n", "reads", "n", "steps", "branch",
                "width", "max", "V0 real", "V0c", "V1 pf", "V2 lvl", "V2 pf", "V3 R=4", "V3 R=16", "check");
    for (int a = 3; a < argc; a++) {
        const std::vector<std::string> reads = Read_Set(argv[a], 1000000);
        if (reads.empty()) continue;
        Tables_For(static_cast<int>(reads[0].size()));
        const int LH = STRINGLENGTH / 2;
        std::vector<Work> w;
        for (const std::string& r : reads) {
            if (static_cast<int>(r.size()) != STRINGLENGTH) continue;
            const std::string_view all(r);
            SARange range = Start_Range(all);
            if (!Search_Right(all.substr(0, LH), 0, range, revfmi) || !range.Start) continue;
            range.Level = 1;
            range.Mismatches = 0;
            range.Mismatch_Char = 0;
            w.push_back({all.substr(LH), range});
        }
        const std::size_t n = w.size();
        auto fresh = [&] { return std::vector<std::vector<Hit>>(n); };

        // Correctness and statistics, untimed.
        std::vector<std::vector<Hit>> h0 = fresh();
        Real(w, LH, &h0);
        unsigned long bad = 0;
        Stats st;
        double widthSum = 0;
        unsigned long rounds = 0, widthMax = 0;
        {
            auto h = fresh(); Stats s; Depth_First<false>(w, LH, h, s); bad += Compare(h0, h);
            h = fresh(); s = {}; Depth_First<true>(w, LH, h, s); bad += Compare(h0, h);
            h = fresh(); s = {}; Frontier<true>(w, LH, 1, h, s, widthSum, rounds, widthMax); bad += Compare(h0, h); st = s;
            double ws = 0; unsigned long rr = 0, wm = 0;
            h = fresh(); s = {}; Frontier<true>(w, LH, 16, h, s, ws, rr, wm); bad += Compare(h0, h);
        }

        double t[7];
        for (double& x : t) x = 1e99;
        for (int rep = 0; rep < 3; rep++) {
            for (int m = 0; m < 7; m++) {
                auto h = fresh();
                Stats s;
                double ws = 0; unsigned long rr = 0, wm = 0;
                auto t0 = Clock::now();
                switch (m) {
                    case 0: Real(w, LH, nullptr); break;
                    case 1: Depth_First<false>(w, LH, h, s); break;
                    case 2: Depth_First<true>(w, LH, h, s); break;
                    case 3: Frontier<false>(w, LH, 1, h, s, ws, rr, wm); break;
                    case 4: Frontier<true>(w, LH, 1, h, s, ws, rr, wm); break;
                    case 5: Frontier<true>(w, LH, 4, h, s, ws, rr, wm); break;
                    case 6: Frontier<true>(w, LH, 16, h, s, ws, rr, wm); break;
                }
                t[m] = std::min(t[m], Seconds(t0));
                sink += s.steps;
            }
        }
        const char* slash = std::strrchr(argv[a], '/');
        const double ns = 1e9 / n;
        std::printf("%-14s %8zu %7.1f %7.2f %6.2f %6lu | %8.0f %8.0f %8.0f %8.0f %8.0f %8.0f %8.0f | %s\n",
                    slash ? slash + 1 : argv[a], n, static_cast<double>(st.steps) / n,
                    static_cast<double>(st.branchings) / n, widthSum / rounds, widthMax, t[0] * ns, t[1] * ns,
                    t[2] * ns, t[3] * ns, t[4] * ns, t[5] * ns, t[6] * ns, bad ? "MISMATCH" : "same hits");
        std::fflush(stdout);
    }
    std::printf("(steps, branch = nodes stepped / branching nodes per read; width = mean live nodes per level of one read; sink %lu)\n",
                sink & 1);
    return 0;
}
