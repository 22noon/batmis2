// brute_force: every hit of every read within K mismatches, by scanning the whole genome.
// The ground truth for the all-hits comparison on genomes small enough to scan
// (SARS-CoV-2, E. coli, chr21).
//
//   brute_force genome.fa reads.fq K [threads] > hits.tsv
//
// A hit is a window of the read's length inside one sequence, holding only A/C/G/T, at
// Hamming distance <= K from the read (forward strand) or from its reverse complement
// (strand -). Windows that cover a non-ACGT base are skipped: the tools treat reference N
// differently (BatMis stores it as C, BWA picks a random base), so hits over N are left
// out of every comparison.
//
// Output: read, contig, pos0, strand, mismatches; tab-separated, in read order, as
// normalize_hits.py writes them except for its last two columns.
//
// Method: the genome and the reads are 2-bit packed; each window is compared 32 bases at a
// time with XOR and popcount, stopping once more than K mismatches are seen. Reads are
// split among threads.
#include <algorithm>
#include <atomic>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include <thread>
#include <vector>

typedef uint64_t u64;

struct Contig {
    std::string Name;
    std::vector<uint8_t> Code;  // 0-3, 4 for anything else
};

static uint8_t Encode(char c) {
    switch (c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default: return 4;
    }
}

// The contig as three bit-planes (bit 0 of each base, bit 1, and non-ACGT), 64 bases a word,
// base i in bit i % 64 of word i / 64; Get reads 64 bases from any start.
struct Packed {
    std::vector<u64> Lo, Hi, Bad;  // bit j of word w: base 64w + j
    size_t N = 0;
    explicit Packed(const std::vector<uint8_t>& C) : N(C.size()) {
        const size_t W = (N + 63) / 64 + 2;
        Lo.assign(W, 0); Hi.assign(W, 0); Bad.assign(W, 0);
        for (size_t i = 0; i < N; i++) {
            const u64 b = 1ULL << (i % 64);
            if (C[i] > 3) { Bad[i / 64] |= b; continue; }
            if (C[i] & 1) Lo[i / 64] |= b;
            if (C[i] & 2) Hi[i / 64] |= b;
        }
    }
    // 64 bits of plane P starting at base s.
    static inline u64 Get(const std::vector<u64>& P, size_t s) {
        const size_t w = s / 64, o = s % 64;
        return o ? (P[w] >> o) | (P[w + 1] << (64 - o)) : P[w];
    }
};

struct Read {
    std::string Name;
    std::vector<u64> Lo[2], Hi[2];  // strand 0 = forward, 1 = reverse complement
};

int main(int argc, char** argv) {
    if (argc < 4) {
        std::fprintf(stderr, "usage: brute_force genome.fa reads.fq K [threads]\n");
        return 1;
    }
    const int K = std::atoi(argv[3]);
    const unsigned Threads = argc > 4 ? std::atoi(argv[4]) : std::thread::hardware_concurrency();

    std::vector<Contig> Genome;
    {
        std::ifstream f(argv[1]);
        std::string line;
        while (std::getline(f, line)) {
            if (line.empty()) continue;
            if (line[0] == '>') {
                Genome.push_back({line.substr(1, line.find_first_of(" \t") - 1), {}});
            } else {
                for (char c : line) Genome.back().Code.push_back(Encode(c));
            }
        }
    }
    std::vector<Packed> Packs;
    for (auto& C : Genome) Packs.emplace_back(C.Code);

    std::vector<Read> Reads;
    size_t L = 0;
    {
        std::ifstream f(argv[2]);
        std::string h, s, p, q;
        while (std::getline(f, h) && std::getline(f, s) && std::getline(f, p) && std::getline(f, q)) {
            if (!L) L = s.size();
            if (s.size() != L) { std::fprintf(stderr, "reads must have one length\n"); return 1; }
            Read R;
            R.Name = h.substr(1, h.find_first_of(" \t") - 1);
            const size_t W = (L + 63) / 64;
            for (int st = 0; st < 2; st++) { R.Lo[st].assign(W, 0); R.Hi[st].assign(W, 0); }
            bool ok = true;
            for (size_t i = 0; i < L; i++) {
                const uint8_t c = Encode(s[i]);
                if (c > 3) { ok = false; break; }
                const uint8_t rc = 3 - c;
                const size_t j = L - 1 - i;
                if (c & 1) R.Lo[0][i / 64] |= 1ULL << (i % 64);
                if (c & 2) R.Hi[0][i / 64] |= 1ULL << (i % 64);
                if (rc & 1) R.Lo[1][j / 64] |= 1ULL << (j % 64);
                if (rc & 2) R.Hi[1][j / 64] |= 1ULL << (j % 64);
            }
            if (ok) Reads.push_back(std::move(R));
        }
    }

    const size_t W = (L + 63) / 64;
    const u64 LastMask = L % 64 ? (1ULL << (L % 64)) - 1 : ~0ULL;
    std::vector<std::string> Out(Reads.size());
    std::atomic<size_t> Next{0};
    auto Work = [&] {
        for (size_t r; (r = Next.fetch_add(1)) < Reads.size();) {
            const Read& R = Reads[r];
            std::string& O = Out[r];
            for (size_t c = 0; c < Genome.size(); c++) {
                const Packed& P = Packs[c];
                if (P.N < L) continue;
                for (size_t s = 0; s + L <= P.N; s++) {
                    for (int st = 0; st < 2; st++) {
                        int d = 0;
                        bool bad = false;
                        for (size_t w = 0; w < W && d <= K; w++) {
                            const u64 m = w + 1 == W ? LastMask : ~0ULL;
                            const size_t at = s + 64 * w;
                            if (Packed::Get(P.Bad, at) & m) { bad = true; break; }
                            const u64 x = ((Packed::Get(P.Lo, at) ^ R.Lo[st][w]) | (Packed::Get(P.Hi, at) ^ R.Hi[st][w])) & m;
                            d += __builtin_popcountll(x);
                        }
                        if (!bad && d <= K) {
                            O += R.Name + '\t' + Genome[c].Name + '\t' + std::to_string(s) + '\t' + (st ? '-' : '+') + '\t' +
                                 std::to_string(d) + '\n';
                        }
                    }
                }
            }
        }
    };
    std::vector<std::thread> T;
    for (unsigned t = 0; t < std::max(1u, Threads); t++) T.emplace_back(Work);
    for (auto& t : T) t.join();
    for (const auto& o : Out) std::fwrite(o.data(), 1, o.size(), stdout);
    return 0;
}
