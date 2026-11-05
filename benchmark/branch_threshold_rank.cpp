// Branch_Detect cost by SA-range width: stored-count lookup versus direct counting.
// Results and the BRANCHTHRESHOLD decision: benchmark/branch_threshold.md.
// Build from src/ (pulls in batman_zero_min.cpp for Count_Range and the index loader):
//   g++ -std=c++17 -fpermissive -w -I. -I../index -O3 -funroll-loops -mpopcnt \
//       ../benchmark/branch_threshold_rank.cpp ../index/libfmindex_core.a -lm -o branch_threshold_rank
// Run:  taskset -c 2 ./branch_threshold_rank <genome prefix used with -g>
#define main orig_main
#include "batman_zero_min.cpp"
#undef main
#include <immintrin.h>
#include <chrono>
#include <random>
#include <vector>
// bench <prefix>: per-call cost of finding all four child ranges of an SA range of a given width.
//   lookup : BWTAllOccValueTwoIndex(Start, End+1)               (two stored-count lookups)
//   count  : BWTAllOccValue(Start) + Count_Range (scalar v2)    (one lookup + bit-plane count)
//   avx2   : BWTAllOccValue(Start) + Count_Range3 (AVX2 middle)
// Regimes: cold = random ranges over the whole index; hot = 256 ranges reused (L1/L2 resident).
//          indep = calls may overlap (memory-level parallelism); dep = each call waits for the last.

__attribute__((target("avx2,popcnt"))) static inline __m256i popcnt_bytes(__m256i v) {
    const __m256i lut = _mm256_setr_epi8(0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4, 0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4);
    const __m256i lo4 = _mm256_set1_epi8(0x0F);
    return _mm256_add_epi8(_mm256_shuffle_epi8(lut, _mm256_and_si256(v, lo4)),
                           _mm256_shuffle_epi8(lut, _mm256_and_si256(_mm256_srli_epi16(v, 4), lo4)));
}
__attribute__((target("avx2,popcnt"))) static inline unsigned long hsum(__m256i v) {
    __m128i s = _mm_add_epi64(_mm256_castsi256_si128(v), _mm256_extracti128_si256(v, 1));
    return (unsigned long)(_mm_cvtsi128_si64(s) + _mm_extract_epi64(s, 1));
}
__attribute__((target("avx2,popcnt"))) void Count_Range3(const BWT* fmi, unsigned long Start, unsigned long End, unsigned Count[4]) {
    const unsigned long inv = fmi->inverseSa0;
    const unsigned long First = Start - (Start > inv), Last = End + 1 - (End + 1 > inv);
    const unsigned* code = fmi->bwtCode;
    if (First >= Last) { Count[0] = Count[1] = Count[2] = Count[3] = 0; return; }
    const unsigned long W0 = First / 16, W1 = (Last - 1) / 16;
    if (W1 - W0 < 16) { Count_Range(fmi, Start, End, Count); return; }
    const unsigned head = ~0u >> (2 * (First % 16)), tail = ~0u << (30 - 2 * ((Last - 1) % 16));
    unsigned long ct = 0, gt = 0, t = 0;
    auto word32 = [&](unsigned w, unsigned m) {
        m &= 0x55555555u; const unsigned l = w & m, h = (w >> 1) & m;
        ct += __builtin_popcount(l); gt += __builtin_popcount(h); t += __builtin_popcount(l & h);
    };
    word32(code[W0], head);
    unsigned long w = W0 + 1;
    const __m256i m55 = _mm256_set1_epi32(0x55555555), zero = _mm256_setzero_si256();
    __m256i vct = zero, vgt = zero, vt = zero;
    for (; w + 8 <= W1; w += 8) {
        const __m256i x = _mm256_loadu_si256((const __m256i*)(code + w));
        const __m256i l = _mm256_and_si256(x, m55), h = _mm256_and_si256(_mm256_srli_epi32(x, 1), m55);
        vct = _mm256_add_epi64(vct, _mm256_sad_epu8(popcnt_bytes(l), zero));
        vgt = _mm256_add_epi64(vgt, _mm256_sad_epu8(popcnt_bytes(h), zero));
        vt  = _mm256_add_epi64(vt,  _mm256_sad_epu8(popcnt_bytes(_mm256_and_si256(l, h)), zero));
    }
    ct += hsum(vct); gt += hsum(vgt); t += hsum(vt);
    for (; w < W1; ++w) word32(code[w], ~0u);
    word32(code[W1], tail);
    const unsigned long n = Last - First;
    Count[3] = t; Count[1] = ct - t; Count[2] = gt - t; Count[0] = n - ct - gt + t;
}

template <int M> static inline void Children(BWT* f, unsigned long s, unsigned long e, unsigned* lo, unsigned* hi) {
    if (M == 0) { BWTAllOccValueTwoIndex(f, s, e + 1, lo, hi); return; }
    unsigned cnt[4];
    BWTAllOccValue(f, s, lo);
    if (M == 1) Count_Range(f, s, e, cnt); else Count_Range3(f, s, e, cnt);
    for (int c = 0; c < 4; c++) hi[c] = lo[c] + cnt[c];
}

template <int M> static double Time(BWT* f, const std::vector<unsigned long>& st, unsigned long width, int n, bool dep, unsigned long& sink) {
    unsigned long d = 0;
    auto t0 = std::chrono::steady_clock::now();
    for (int i = 0; i < n; i++) {
        unsigned int ALIGN_16 lo[4], hi[4];
        const unsigned long s = st[i % st.size()] + d;
        Children<M>(f, s, s + width - 1, lo, hi);
        sink += hi[0] - lo[0] + hi[3];
        if (dep) d = (lo[0] | hi[3]) >> 31;        // always 0 (counts < 2^31), but the CPU must wait for it
    }
    return std::chrono::duration<double, std::nano>(std::chrono::steady_clock::now() - t0).count() / n;
}

int main(int, char** argv) {
    std::string pre = argv[1], rb, ro, fb, fo; BuildIndexPaths(pre, rb, ro, fb, fo);
    BWT* f = initFMI(rb.c_str(), ro.c_str());
    const unsigned long T = f->textLength;
    std::mt19937_64 rng(7); unsigned long sink = 0, bad = 0;
    // correctness of all three on random ranges, including ones around $
    for (int i = 0; i < 1000000; i++) {
        unsigned long w = 1 + rng() % (i & 1 ? 200 : 20000), s = (i % 5 == 0 && f->inverseSa0 > w) ? f->inverseSa0 - rng() % w : rng() % (T - w);
        unsigned int ALIGN_16 a[4], b[4], c[4], x[4], y[4], z[4];
        Children<0>(f, s, s + w - 1, a, x); Children<1>(f, s, s + w - 1, b, y); Children<2>(f, s, s + w - 1, c, z);
        bad += memcmp(a, b, 16) || memcmp(a, c, 16) || memcmp(x, y, 16) || memcmp(x, z, 16);
    }
    printf("text=%lu  check: 1M ranges, mismatches=%lu\n\n", T, bad);
    const unsigned long widths[] = {1, 2, 4, 8, 16, 24, 32, 48, 64, 96, 128, 192, 256, 384, 512, 1024, 2048, 4096, 16384, 65536};
    for (int regime = 0; regime < 3; regime++) {
        const bool hot = regime == 2, dep = regime != 0;
        printf("%s\n%7s %9s %9s %9s %9s %9s\n", regime == 0 ? "COLD, independent calls (ns/call)" : regime == 1 ? "COLD, dependent calls (ns/call)" : "HOT (256 ranges reused), dependent calls (ns/call)",
               "rows", "lookup", "count", "avx2", "count/lk", "avx2/lk");
        for (unsigned long w : widths) {
            const int n = hot ? 2000000 : w >= 16384 ? 100000 : 1000000;
            std::vector<unsigned long> st(hot ? 256 : n);
            for (auto& s : st) s = rng() % (T - w - 2);
            double r[3] = {1e99, 1e99, 1e99};
            for (int rep = 0; rep < 3; rep++) {
                if (hot) { r[0] = std::min(r[0], Time<0>(f, st, w, n, dep, sink)); r[1] = std::min(r[1], Time<1>(f, st, w, n, dep, sink)); r[2] = std::min(r[2], Time<2>(f, st, w, n, dep, sink)); }
                else { // fresh positions per method so nothing is cached from the previous method
                    for (auto& s : st) s = rng() % (T - w - 2); r[0] = std::min(r[0], Time<0>(f, st, w, n, dep, sink));
                    for (auto& s : st) s = rng() % (T - w - 2); r[1] = std::min(r[1], Time<1>(f, st, w, n, dep, sink));
                    for (auto& s : st) s = rng() % (T - w - 2); r[2] = std::min(r[2], Time<2>(f, st, w, n, dep, sink));
                }
            }
            printf("%7lu %9.1f %9.1f %9.1f %8.2fx %8.2fx\n", w, r[0], r[1], r[2], r[1] / r[0], r[2] / r[0]);
        }
        printf("\n");
    }
    if (sink == 42) puts("");
}
