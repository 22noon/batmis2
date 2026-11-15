// bi_check <genome_prefix> [patterns] [seed]: checks the bidirectional step.
//
// Each pattern (1-60 bases: from the text, many at its start or end, or random) is grown
// from a random point by a random mix of Extend_Left and Extend_Right. After every step:
//   - the Fw range equals the classic backward search of the pattern in Fw, and the Rev
//     range the classic search in Rev (Backward_Search, Forward_Search);
//   - the two have the same size, and it equals the occurrences counted in the text
//     (<genome_prefix>.bin) by brute force, for patterns of up to 12 bases;
//   - for ranges of up to 16 rows, the positions decoded from Fw and from Rev are the same
//     set, and the text holds the pattern at each;
//   - a step from a single row (the singleton step) or a narrow range (the narrow step,
//     Count_Range) gives the same Bi_Range as the general step (Extend_General);
//   - Children_Left / Children_Right give, for each base, the range Extend_Left /
//     Extend_Right gives.
// Then the inline Occ (Occ_All, Occ_One) equals BWTAllOccValue / BWTOccValue in both
// indexes, at every row of texts up to 4 Mb, else at 10M random rows and every row near
// the $ row and near 1000 random checkpoints. And every entry of a depth-8 lookup table (Build_Lookup; depth 6 for texts under
// 64 kb) equals the k-mer's range searched base by base.
// Prints one line per failure (the first 20) and a summary; exit code 1 on any failure.
#include <algorithm>
#include <cstdio>
#include <fstream>
#include <iterator>
#include <random>
#include <string>
#include <vector>

#include "bi_index.hpp"

using namespace bifm;

int main(int argc, char** argv) {
    if (argc < 2) {
        std::fprintf(stderr, "usage: %s <genome_prefix> [patterns] [seed]\n", argv[0]);
        return 2;
    }
    const std::string Prefix = argv[1];
    const int Patterns = argc > 2 ? std::atoi(argv[2]) : 2000;
    std::mt19937_64 Rng(argc > 3 ? std::strtoull(argv[3], nullptr, 10) : 1);

    std::ifstream Bin(Prefix + ".bin", std::ios::binary);
    std::string Text((std::istreambuf_iterator<char>(Bin)), std::istreambuf_iterator<char>());
    for (char& c : Text) c = static_cast<char>(Base_Code(c));
    Bi_Index Index(Prefix);
    const unsigned long N = Index.Text_Length();
    if (Text.size() != N) {
        std::fprintf(stderr, "%s.bin has %zu bases, the index %lu\n", Prefix.c_str(), Text.size(), N);
        return 2;
    }

    auto Pick = [&](unsigned long Hi) { return std::uniform_int_distribution<unsigned long>(0, Hi)(Rng); };
    unsigned long Steps = 0, Children_Checked = 0, Located = 0, Counted = 0, Singleton_Steps = 0, Narrow_Steps = 0, Failures = 0;
    auto Fail = [&](const std::string& What, const std::string& P) {
        if (++Failures <= 20) {
            std::string S;
            for (char c : P) S += "ACGT"[static_cast<int>(c)];
            std::fprintf(stderr, "FAIL %s: %s\n", What.c_str(), S.c_str());
        }
    };
    for (int t = 0; t < Patterns; t++) {
        const unsigned long Len = 1 + Pick(59);
        std::string P(Len, 0);
        const int Kind = static_cast<int>(Pick(9));
        if (Kind == 0) {  // random bases
            for (char& c : P) c = static_cast<char>(Pick(3));
        } else {  // from the text: at its start (1 in 9), at its end (1 in 9), anywhere
            if (Len > N) continue;
            const unsigned long At = Kind == 1 ? 0 : Kind == 2 ? N - Len : Pick(N - Len);
            P = Text.substr(At, Len);
        }
        // Grow P[i, j) from a random point to the whole pattern.
        unsigned long i = Pick(Len), j = i;
        Bi_Range R = Index.Empty_Pattern();
        while (i > 0 || j < Len) {
            const bool Left = j == Len || (i > 0 && Pick(1));
            const Bi_Range Prev = R;
            const int c = Left ? P[--i] : P[j++];
            {
                Bi_Range Out[4];
                Left ? Index.Children_Left(Prev, Out) : Index.Children_Right(Prev, Out);
                for (int b = 0; b < 4; b++) {
                    const Bi_Range E = Left ? Index.Extend_Left(Prev, b) : Index.Extend_Right(Prev, b);
                    if (E.Empty() != Out[b].Empty() ||
                        (!E.Empty() && (E.Fw_Start != Out[b].Fw_Start || E.Fw_End != Out[b].Fw_End || E.Rev_Start != Out[b].Rev_Start)))
                        Fail(std::string(Left ? "left" : "right") + " children differ from the extensions", P.substr(i, j - i));
                }
                Children_Checked++;
            }
            R = Left ? Index.Extend_Left(R, c) : Index.Extend_Right(R, c);
            Steps++;
            const std::string Q = P.substr(i, j - i);
            if (Prev.Size() < Index.Threshold + 1) {  // singleton or narrow step
                const Bi_Range G = Left ? Bi_Index::Extend_General(Index.Fw(), Prev.Fw(), Prev.Rev_Start, c, true)
                                        : Bi_Index::Extend_General(Index.Rev(), Prev.Rev(), Prev.Fw_Start, c, false);
                const char* Kind = Prev.Size() == 1 ? " singleton" : " narrow";
                if (G.Empty() != R.Empty() ||
                    (!R.Empty() && (G.Fw_Start != R.Fw_Start || G.Fw_End != R.Fw_End || G.Rev_Start != R.Rev_Start)))
                    Fail(std::string(Left ? "left" : "right") + Kind + " step differs from the general step", Q);
                (Prev.Size() == 1 ? Singleton_Steps : Narrow_Steps)++;
            }
            const Range Fw = Index.Backward_Search(Q), Rev = Index.Forward_Search(Q);
            if (Fw.Size() != Rev.Size()) Fail("classic Fw and Rev sizes differ", Q);
            if (R.Size() != Fw.Size() || (!R.Empty() && (R.Fw_Start != Fw.Start || R.Fw_End != Fw.End)))
                Fail(std::string(Left ? "left" : "right") + " step: Fw range", Q);
            if (!R.Empty() && R.Rev().Start != Rev.Start) Fail(std::string(Left ? "left" : "right") + " step: Rev range", Q);
            if (Q.size() <= 12) {  // brute-force count
                unsigned long Count = 0;
                for (std::size_t p = Text.find(Q); p != std::string::npos; p = Text.find(Q, p + 1)) Count++;
                if (Count != R.Size()) Fail("count against the text", Q);
                Counted++;
            }
            if (!R.Empty() && R.Size() <= 16) {
                std::vector<unsigned long> A, B;
                for (unsigned long Row = R.Fw_Start; Row <= R.Fw_End; Row++) A.push_back(Index.Locate_Fw(Row));
                for (unsigned long Row = R.Rev().Start; Row <= R.Rev().End; Row++) B.push_back(Index.Locate_Rev(Row, Q.size()));
                std::sort(A.begin(), A.end());
                std::sort(B.begin(), B.end());
                if (A != B) Fail("Fw and Rev positions differ", Q);
                for (unsigned long p : A)
                    if (Text.compare(p, Q.size(), Q) != 0) Fail("a position does not hold the pattern", Q);
                Located++;
            }
            if (R.Empty()) break;
        }
    }
    // The inline Occ against the library.
    unsigned long Occ_Rows = 0;
    for (const BWT* B : {Index.Fw(), Index.Rev()}) {
        auto Check_Row = [&](unsigned long Row) {
            if (Row > N + 1) return;
            unsigned Want[4], Got[4];
            BWTAllOccValue(B, static_cast<unsigned>(Row), Want);
            Bi_Index::Occ_All(B, Row, Got);
            for (int b = 0; b < 4; b++)
                if (Got[b] != Want[b] || Bi_Index::Occ_One(B, Row, b) != BWTOccValue(B, static_cast<unsigned>(Row), b) ||
                    Want[b] != BWTOccValue(B, static_cast<unsigned>(Row), b))
                    Fail("Occ at row " + std::to_string(Row), "");
            Occ_Rows++;
        };
        if (N <= 4000000) {
            for (unsigned long Row = 0; Row <= N + 1; Row++) Check_Row(Row);
        } else {
            for (int k = 0; k < 10000000; k++) Check_Row(Pick(N + 1));
            for (unsigned long Row = B->inverseSa0 > 300 ? B->inverseSa0 - 300 : 0; Row <= B->inverseSa0 + 300; Row++) Check_Row(Row);
            for (int k = 0; k < 1000; k++) {
                const unsigned long At = Pick((N + 1) / OCC_INTERVAL) * OCC_INTERVAL;
                for (unsigned long Row = At > 140 ? At - 140 : 0; Row <= At + 140; Row++) Check_Row(Row);
            }
        }
    }

    // The lookup table against base-by-base searches.
    const int K = N < 65536 ? 6 : 8;
    Index.Build_Lookup(K);
    unsigned long Entries = 0, Nonempty = 0;
    std::string Q(K, 0);
    for (unsigned long Key = 0; Key < 1ul << (2 * K); Key++) {
        for (int b = 0; b < K; b++) Q[b] = static_cast<char>(Key >> (2 * (K - 1 - b)) & 3);
        Bi_Range W = Index.Empty_Pattern();
        for (int b = K; b-- > 0 && !W.Empty();) W = Index.Extend_Left(W, Q[b]);
        const Bi_Range T = Index.Look_Up(Q.data());
        if (W.Empty() != T.Empty() ||
            (!W.Empty() && (W.Fw_Start != T.Fw_Start || W.Fw_End != T.Fw_End || W.Rev_Start != T.Rev_Start)))
            Fail("lookup entry", Q);
        Entries++;
        Nonempty += !T.Empty();
    }
    std::printf("patterns=%d steps=%lu singleton_steps=%lu narrow_steps=%lu children=%lu counted=%lu located=%lu "
                "lookup_entries=%lu nonempty=%lu occ_rows=%lu failures=%lu\n",
                Patterns, Steps, Singleton_Steps, Narrow_Steps, Children_Checked, Counted, Located, Entries, Nonempty,
                Occ_Rows, Failures);
    return Failures ? 1 : 0;
}
