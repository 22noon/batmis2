// Brute-force test of Genome_Map::Place: for every text position G and each read length
// L, Place must accept [G, G + L) exactly when it lies inside one sequence (not the
// padding sequence) and holds no non-ACGT base of the FASTA, and then give that
// sequence's SAM header place and the in-sequence position.
// Usage: test_genome_map <genome_prefix> <genome.fa>
#include <cctype>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

#include "genome_map.hpp"

int main(int argc, char** argv) {
    if (argc < 3) {
        std::fprintf(stderr, "Usage: %s <genome_prefix> <genome.fa>\n", argv[0]);
        return 2;
    }
    Genome_Map Map;
    std::string Error;
    if (!Map.Load(argv[1], Error)) {
        std::fprintf(stderr, "%s\n", Error.c_str());
        return 2;
    }
    // The text as the index has it: per position, its sequence and whether it is an N.
    std::ifstream In(argv[2]);
    std::string Line;
    std::vector<int> Seq_Of;
    std::vector<char> Is_N;
    int Seq = -1;
    while (std::getline(In, Line)) {
        if (!Line.empty() && Line[0] == '>') {
            Seq++;
            continue;
        }
        for (char C : Line) {
            if (C == '\r') continue;
            const char U = static_cast<char>(std::toupper(static_cast<unsigned char>(C)));
            Seq_Of.push_back(Seq);
            Is_N.push_back(!(U == 'A' || U == 'C' || U == 'G' || U == 'T'));
        }
    }
    const long Total = static_cast<long>(Seq_Of.size());
    std::vector<long> Seq_Start(Total);  // the text position where each position's sequence starts
    for (long G = 0; G < Total; G++) Seq_Start[G] = G > 0 && Seq_Of[G - 1] == Seq_Of[G] ? Seq_Start[G - 1] : G;
    unsigned long Checked = 0, Accepted = 0, Failures = 0;
    for (int L : {32, 75, 100}) {
        long Ns = 0;  // N count in [G, G + L)
        for (long G = 0; G < L && G < Total; G++) Ns += Is_N[G];
        for (long G = 0; G < Total; G++) {
            if (G > 0 && G + L - 1 < Total) Ns += Is_N[G + L - 1] - Is_N[G - 1];
            else if (G > 0) Ns -= Is_N[G - 1];
            const bool Inside = G + L <= Total && Seq_Of[G] == Seq_Of[G + L - 1] && Map.Tid_Of(Seq_Of[G]) >= 0;
            const bool Want = Inside && Ns == 0;
            int Tid = -1;
            int64_t Pos = -1;
            const bool Got = Map.Place(static_cast<uint64_t>(G), L, Tid, Pos);
            Checked++;
            Accepted += Got;
            bool Ok = Got == Want;
            if (Ok && Got) Ok = Tid == Map.Tid_Of(Seq_Of[G]) && Pos == G - Seq_Start[G];
            if (!Ok && Failures++ < 10)
                std::fprintf(stderr, "FAIL G=%ld L=%d: want %d (seq %d), got %d (tid %d pos %lld)\n", G, L, Want,
                             G < Total ? Seq_Of[G] : -1, Got, Tid, static_cast<long long>(Pos));
        }
    }
    std::printf("positions=%ld checked=%lu accepted=%lu failures=%lu\n", Total, Checked, Accepted, Failures);
    return Failures ? 1 : 0;
}
