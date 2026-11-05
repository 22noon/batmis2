// Resolving a small SA range on the reference: locate every row of the range (its SA
// value), compare the rest of the read with the genome (.pac) at each occurrence, and
// take the final rows from a sampled inverse SA instead of walking the remaining bases.
// Stand-alone, no change to src/. Results: benchmark/sampled_end.md.
// Build from src/ (pulls in batman_zero_min.cpp):
//   g++ -std=c++17 -fpermissive -w -I. -I../index -O3 -funroll-loops -mpopcnt \
//       ../benchmark/sampled_end.cpp ../index/libfmindex_core.a -lm -o sampled_end
// Run:  taskset -c 2 ./sampled_end <genome prefix> <isa cache file> <lookup k> <reads.fq>...
//
// Range [S, E] after matching read[0..m), k = L - m bases left (revfmi):
//   - LF keeps the order of rows with the same preceding character, so the rows whose
//     occurrences match the rest map, in order, onto the final range: it holds exactly
//     the c matching rows and ends on E_f = ISA(P_last - k), P_last the text position
//     of the last matching row. No matching row: no hit. One: its genome position is
//     the hit (the located form today's search reports).
//   - Not true: "if S and E both match, so do the rows between them" (rows are ordered
//     by the text after the matched part, while the search adds bases in front of it).
// Variants (reference on; ranges that become one row finish on today's located path):
//   today : Search_Right.
//   W2E   : the range has 2 rows and a sampled End whose rest matches: locate S only.
//   Wn    : the range has 2..n rows (n = 2, 4, 8, 16): locate all rows together.
#define main orig_main
#include "batman_zero_min.cpp"
#undef main
#include <chrono>

namespace {

using Clock = std::chrono::steady_clock;
double Seconds(Clock::time_point t) { return std::chrono::duration<double>(Clock::now() - t).count(); }

unsigned long N;  // text length

inline int Bwt_Char(unsigned long row) {  // -1 at the $ row
    if (row == revfmi->inverseSa0) return -1;
    const unsigned long i = row - (row > revfmi->inverseSa0);
    return revfmi->bwtCode[i / 16] << ((i % 16) * 2) >> (BITS_IN_WORD - 2);
}
inline unsigned long Lf(unsigned long row, unsigned c) { return revfmi->cumulativeFreq[c] + BWTOccValue(revfmi, row, c) + 1; }

// Sampled inverse SA: Isa[k] = the row whose suffix starts at text position 8k.
std::vector<unsigned> Isa;

void Build_Isa(const std::string& cache) {
    Isa.assign(N / SAINTERVAL + 2, 0);
    {
        std::ifstream in(cache, std::ios::binary);
        if (in && in.read(reinterpret_cast<char*>(Isa.data()), Isa.size() * sizeof(unsigned))) return;
    }
    // Row 0 is the suffix "$" (position N); LF moves to position - 1.
    unsigned long row = 0;
    for (unsigned long p = N;; p--) {
        if (p % SAINTERVAL == 0) Isa[p / SAINTERVAL] = row;
        if (p == 0) break;
        row = Lf(row, Bwt_Char(row));
    }
    if (row != revfmi->inverseSa0) std::fprintf(stderr, "inverse SA walk ended on row %lu, not $ row %u\n", row, revfmi->inverseSa0);
    std::ofstream(cache, std::ios::binary).write(reinterpret_cast<const char*>(Isa.data()), Isa.size() * sizeof(unsigned));
}

// The row of the suffix starting at text position q, as BWTInverseSa: from the next
// sampled position, LF back (at most 7 steps).
inline unsigned long Row_Of(unsigned long q) {
    unsigned long k = (q + SAINTERVAL - 1) / SAINTERVAL, row, skip;
    if (k * SAINTERVAL > N) { row = 0; skip = N - q; }
    else { row = Isa[k]; skip = k * SAINTERVAL - q; }
    while (skip--) row = Lf(row, Bwt_Char(row));
    return row;
}

// The text position of row r: LF with its own characters to a sampled row, as
// BWTSaValue. Rows of a range are stepped together so their misses overlap.
void Locate_All(const unsigned long* rows, unsigned w, unsigned long* P, unsigned long& steps) {
    unsigned long r[64], t[64];
    unsigned pending = 0, idx[64];
    for (unsigned j = 0; j < w; j++) {
        r[j] = rows[j];
        t[j] = 0;
        idx[pending++] = j;
    }
    while (pending) {
        unsigned keep = 0;
        for (unsigned q = 0; q < pending; q++) {
            const unsigned j = idx[q];
            if (r[j] % SAINTERVAL == 0) { P[j] = revfmi->saValue[r[j] / SAINTERVAL] + t[j]; continue; }
            if (r[j] == revfmi->inverseSa0) { P[j] = t[j]; continue; }
            r[j] = Lf(r[j], Bwt_Char(r[j]));
            t[j]++;
            steps++;
            idx[keep++] = j;
        }
        pending = keep;
    }
}

// Genome start of the occurrence at text position P after m bases, if read[m..L) matches there.
inline bool Rest_Matches_At(unsigned long P, int m, int L, long& s) {
    s = static_cast<long>(N) - m - static_cast<long>(P);
    return s >= 0 && s + L <= static_cast<long>(N) && Genome_Mismatches(s, m, L) == 0;
}

struct Stats {
    unsigned long triggers = 0, widthSum = 0, matched = 0, noHit = 0, oneHit = 0, multiHit = 0, locateSteps = 0,
                  isaSteps = 0, basesSkipped = 0, assertFail = 0;
};

// Resolves R ([S, E] after m bases); P/known: E's position if already known (W2E).
bool Resolve(std::string_view read, SARange& R, const unsigned long* knownE, Stats* st) {
    const int L = static_cast<int>(read.size()), m = R.Level - 1, k = L - m;
    const unsigned w = R.End - R.Start + 1;
    unsigned long rows[64], P[64], steps = 0;
    for (unsigned j = 0; j < w; j++) rows[j] = R.Start + j;
    if (knownE) { Locate_All(rows, w - 1, P, steps); P[w - 1] = *knownE; }
    else Locate_All(rows, w, P, steps);
    unsigned c = 0, last = 0;
    long sLast = 0;
    for (unsigned j = 0; j < w; j++) {
        long s;
        if (Rest_Matches_At(P[j], m, L, s)) { c++; last = j; sLast = s; }
    }
    SARange walked = R;  // stats only: the range walked to the end
    if (st) {
        st->triggers++; st->widthSum += w; st->locateSteps += steps; st->basesSkipped += k;
        for (int i = m; i < L && walked.Start; i++) Get_SARange(read[i], walked, revfmi);
    }
    if (!c) {
        if (st) { st->noHit++; st->assertFail += walked.Start != 0; }
        R.Start = 0;
        return false;
    }
    if (c == 1) {  // one occurrence: report its genome position
        if (st) { st->oneHit++; st->assertFail += !walked.Start || walked.End != walked.Start; }
        R.Skip = 1; R.Pos = sLast; R.Start = R.End = rows[last]; R.Level = L;
        return true;
    }
    if (st) st->multiHit++;
    const unsigned long q = P[last] - k;  // last matching occurrence after k more bases
    const unsigned long Ef = Row_Of(q);
    if (st) {
        st->isaSteps += (q + SAINTERVAL - 1) / SAINTERVAL * SAINTERVAL - q;
        st->assertFail += walked.End != Ef || walked.Start != Ef - c + 1;
    }
    R.Start = Ef - c + 1; R.End = Ef; R.Level = L;
    return true;
}

// Width 2, End's occurrence confirmed on the reference (genome start sE), no inverse SA:
// Start follows the single-row extension (its character checked at each base; with
// LOCATED, at a sampled row its own rest is compared on the .pac). Exactly the rows that
// match survive, in order, so:
//   Start fails      -> only End's occurrence matches: a located hit at sE.
//   Start survives   -> both match: the final range is [S_f, S_f + 1].
template <bool LOCATED>
bool Resolve_W2S(std::string_view read, SARange& R, long sE, Stats* st) {
    const int L = static_cast<int>(read.size()), m = R.Level - 1;
    unsigned long S = R.Start;
    bool sConfirmed = false;  // S's rest already matched on the .pac
    SARange walked = R;       // stats only
    if (st) {
        st->triggers++; st->widthSum += 2;
        for (int i = m; i < L && walked.Start; i++) Get_SARange(read[i], walked, revfmi);
    }
    auto endOnly = [&](int i) {  // Start dropped out after i bases: End's occurrence is the hit
        if (st) { st->oneHit++; st->basesSkipped += L - i; st->assertFail += !walked.Start || walked.End != walked.Start; }
        R.Skip = 1; R.Pos = sE; R.Start = R.End = R.End; R.Level = L;
        return true;
    };
    for (int i = m; i < L; i++) {
        if (LOCATED && !sConfirmed && S % SAINTERVAL == 0) {
            const long s = static_cast<long>(N) - i - static_cast<long>(revfmi->saValue[S / SAINTERVAL]);
            if (s < 0 || s + L > static_cast<long>(N) || Genome_Mismatches(s, i, L)) return endOnly(i);
            sConfirmed = true;  // both occurrences match: walk on only for S's final row
        }
        const int c = Bwt_Char(S);
        const unsigned expected = static_cast<unsigned char>(read[i]);
        const unsigned base = Read_Base(expected);
        if (!sConfirmed && c != static_cast<int>(expected)) return endOnly(i);
        S = Lf(S, base);
        if (st) st->locateSteps++;
    }
    if (st) { st->multiHit++; st->assertFail += walked.Start != S || walked.End != S + 1; }
    R.Start = S; R.End = S + 1; R.Level = L;
    return true;
}

// Today's multi-row phase; N_MAX = 0 for today, -1 for W2E, -2 / -3 for W2S without /
// with the located check on Start, else the width limit n.
template <int N_MAX>
bool Search(std::string_view read, SARange& R, Stats* st) {
    const int L = static_cast<int>(read.size());
    if (R.Level > L) return R.Start != 0;
    if (!R.Start) return false;
    while (R.End != R.Start && !R.Skip) {
        const unsigned long w = R.End - R.Start + 1;
        if (N_MAX > 0 && w <= static_cast<unsigned long>(N_MAX)) return Resolve(read, R, nullptr, st);
        if (N_MAX < 0 && w == 2 && R.End % SAINTERVAL == 0) {
            long s;
            const unsigned long PE = revfmi->saValue[R.End / SAINTERVAL];
            if (Rest_Matches_At(PE, R.Level - 1, L, s)) {
                if (N_MAX == -1) return Resolve(read, R, &PE, st);
                if (N_MAX == -2) return Resolve_W2S<false>(read, R, s, st);
                return Resolve_W2S<true>(read, R, s, st);
            }
        }
        Get_SARange(read[R.Level - 1], R, revfmi);
        if (!R.Start) return false;
        if (R.Level == L) return true;
        R.Level++;
    }
    return Search_Right(read, 0, R, revfmi);
}

using Search_Fn = bool (*)(std::string_view, SARange&, Stats*);
const Search_Fn Variants[] = {Search<0>, Search<-2>, Search<-3>, Search<-1>, Search<2>, Search<4>};
const char* Names[] = {"today", "W2S", "W2S+loc", "W2E", "W2", "W4"};
constexpr int NV = 6;

}  // namespace

int main(int argc, char** argv) {
    if (argc < 5) {
        std::fprintf(stderr, "Usage: %s <genome prefix> <isa cache> <lookup k> <reads.fq>...\n", argv[0]);
        return 2;
    }
    std::string pre = argv[1], rb, ro, fb, fo, rs, fs;
    BuildIndexPaths(pre, rb, ro, fb, fo, &rs, &fs);
    revfmi = initFMI(rb.c_str(), ro.c_str(), rs.c_str());
    if (!LoadPac(pre + ".pac")) return 2;
    Use_Ref = true;
    N = Text_Length = revfmi->textLength;
    LOOKUP_DEPTH = std::atoi(argv[3]);

    auto t = Clock::now();
    Build_Isa(argv[2]);
    unsigned long isaBad = 0;
    for (unsigned long k = 0; k < 100000; k++) {  // ISA(SA(r)) = r on rows spread over the index
        const unsigned long r = 1 + (k * 2654435761ul) % N;
        isaBad += Row_Of(BWTSaValue(revfmi, r)) != r;
    }
    std::printf("text=%lu  sampled inverse SA: %.0f MB, ready in %.1f s, check on 100k rows: %s\n", N,
                Isa.size() * 4 / 1048576.0, Seconds(t), isaBad ? "MISMATCH" : "ok");

    for (int a = 4; a < argc; a++) {
        std::ifstream in(argv[a]);
        std::string id, seq;
        std::vector<std::string> reads;
        while (ReadNextRecord(in, id, seq) && reads.size() < 1000000) {
            std::string r(seq.size(), 0);
            bool ok = true;
            for (std::size_t i = 0; i < seq.size() && ok; i++) ok = (r[i] = static_cast<char>(EncodeBase(seq[i]))) >= 0;
            if (ok && (reads.empty() || r.size() == reads[0].size())) reads.push_back(std::move(r));
        }
        STRINGLENGTH = static_cast<int>(reads[0].size());
        Build_Tables(revfmi);

        // Correctness, statistics, and the split of reads by today's final range size.
        auto positions = [&](const SARange& x) {
            std::vector<long> p;
            if (x.Skip) { p.push_back(x.Pos); return p; }
            for (unsigned long row = x.Start; row <= x.End && p.size() < 5000; row++)
                p.push_back(static_cast<long>(N) - STRINGLENGTH - static_cast<long>(BWTSaValue(revfmi, row)));
            std::sort(p.begin(), p.end());
            return p;
        };
        Stats st[NV];
        unsigned long differ[NV] = {}, multiReads = 0, hitReads = 0;
        std::vector<const std::string*> sets[3];  // all, final one row (or no hit), final several rows
        for (const std::string& r : reads) {
            const std::string_view all(r);
            Pack_Read(all);
            SARange base = Start_Range(all);
            const bool f0 = Search<0>(all, base, nullptr);
            const std::vector<long> p0 = f0 ? positions(base) : std::vector<long>();
            const bool multi = f0 && !base.Skip && base.End > base.Start;
            multiReads += multi;
            hitReads += f0;
            sets[0].push_back(&r);
            sets[multi ? 2 : 1].push_back(&r);
            for (int v = 1; v < NV; v++) {
                SARange x = Start_Range(all);
                const bool f = Variants[v](all, x, &st[v]);
                differ[v] += f != f0 || (f && positions(x) != p0);
            }
        }

        double tm[NV][3];
        unsigned long sink = 0;
        for (int v = 0; v < NV; v++)
            for (int g = 0; g < 3; g++) {
                tm[v][g] = 1e99;
                if (sets[g].empty()) continue;
                for (int rep = 0; rep < 3; rep++) {
                    t = Clock::now();
                    for (const std::string* r : sets[g]) {
                        const std::string_view all(*r);
                        Pack_Read(all);
                        SARange x = Start_Range(all);
                        sink += Variants[v](all, x, nullptr) + x.Start;
                    }
                    tm[v][g] = std::min(tm[v][g], Seconds(t));
                }
            }

        const char* slash = std::strrchr(argv[a], '/');
        const unsigned long n = reads.size();
        std::printf("\n%s  (%lu reads of %d bp, lookup %d): hits %.1f%%, final range of several rows %.2f%%\n",
                    slash ? slash + 1 : argv[a], n, STRINGLENGTH, LOOKUPSIZE, 100.0 * hitReads / n, 100.0 * multiReads / n);
        std::printf("  %-6s %9s %9s %9s | %7s %6s %6s %6s %6s | %7s %7s %7s | %5s %5s\n", "", "ns all", "ns 1-row",
                    "ns multi", "trig%", "width", "none", "one", "many", "locate", "isa", "skipped", "diff", "assert");
        for (int v = 0; v < NV; v++) {
            const Stats& s = st[v];
            const double tr = s.triggers ? s.triggers : 1;
            std::printf("  %-6s %9.0f %9.0f %9.0f | %6.2f%% %6.2f %6.2f %6.2f %6.2f | %7.2f %7.2f %7.1f | %5lu %5lu\n",
                        Names[v], tm[v][0] * 1e9 / sets[0].size(), sets[1].empty() ? 0 : tm[v][1] * 1e9 / sets[1].size(),
                        sets[2].empty() ? 0 : tm[v][2] * 1e9 / sets[2].size(), 100.0 * s.triggers / n, s.widthSum / tr,
                        s.noHit / tr, s.oneHit / tr, s.multiHit / tr, s.locateSteps / tr,
                        s.multiHit ? (double)s.isaSteps / s.multiHit : 0.0, s.basesSkipped / tr, differ[v], s.assertFail);
        }
        std::fflush(stdout);
        if (sink == 42) std::printf(" ");
    }
    return 0;
}
