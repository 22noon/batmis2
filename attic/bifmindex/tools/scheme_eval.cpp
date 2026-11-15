// scheme_eval: prototype of bidirectional search schemes for one mismatch class (step 1 of
// the scheme plan). Counts search nodes per read for a scheme set, to compare with
// batman_zero_min's steps for the same class, and prints the hits (pos:mismatches) for a
// brute-force check.
//
//   scheme_eval -g prefix -q reads [-p parts] [-k total] [-lh n] [--lookup k] [--ref] [--hits]
//
// The class: total -k mismatches with exactly -lh in LH (= the first half of the parts);
// e.g. -k 4 -lh 2 is -n 4's 2|2 class. The read is cut into -p parts (even; LH is the
// first p/2). Every distribution of counts over the parts in the class is a case; each case
// gets a seed (its longest run of exact parts, the most central on ties) and an order (from
// the seed, extend toward the neighbouring part with the fewer mismatches; ties toward the
// side with fewer mismatches left). Cases with the same seed and order form a group,
// searched once: a depth-first walk carrying the bitmask of its cases still possible
// (the counts of finished parts equal, the current part's not above); a mismatch is tried
// where some case still needs one. A sampled single row (either index) locates the node's
// one occurrence: a hit if the genome's mismatch distribution is a case of the group.
#include <algorithm>
#include <array>
#include <bitset>
#include <chrono>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "bi_index.hpp"
#include "bi_ref.hpp"

using namespace bifm;

namespace {

constexpr int MAXP = 16;
constexpr int MAXCASES = 1024;
using Mask = std::bitset<MAXCASES>;

struct Case {
    unsigned char C[MAXP];
};

struct Group {
    int Seed_Lo, Seed_Hi;      // parts [Seed_Lo, Seed_Hi) exact, searched first (leftward from Seed_Hi)
    std::vector<int> Order;    // the parts after the seed, in order
    std::vector<int> Cases;    // indices into the case list
};

int P = 8, K = 4, LH_Need = 2;
std::vector<Case> Cases;
std::vector<Group> Groups;

void Enumerate(std::vector<int>& Cur, int i, int Left_LH, int Left_RH) {
    if (i == P) {
        if (Left_LH == 0 && Left_RH == 0) {
            Case C{};
            for (int t = 0; t < P; t++) C.C[t] = static_cast<unsigned char>(Cur[t]);
            Cases.push_back(C);
        }
        return;
    }
    const bool In_LH = i < P / 2;
    const int Left = In_LH ? Left_LH : Left_RH;
    for (int c = 0; c <= Left; c++) {
        Cur[i] = c;
        Enumerate(Cur, i + 1, In_LH ? Left_LH - c : Left_LH, In_LH ? Left_RH : Left_RH - c);
    }
}

// Seeds by greedy set cover: runs [lo, hi) of exact parts, longest first, each time the
// run covering the most cases left; one group per seed, with one order for all its cases
// (at each turn toward the side whose next part has the fewer mismatches summed over them).
void Build_Groups() {
    std::vector<char> Covered(Cases.size(), 0);
    std::size_t Left = Cases.size();
    for (int Len = P; Len >= 1 && Left; Len--) {
        for (;;) {
            int Best_Lo = -1;
            std::size_t Best_N = 0;
            for (int lo = 0; lo + Len <= P; lo++) {
                std::size_t N = 0;
                for (std::size_t ci = 0; ci < Cases.size(); ci++) {
                    if (Covered[ci]) continue;
                    bool Exact = true;
                    for (int t = lo; t < lo + Len && Exact; t++) Exact = !Cases[ci].C[t];
                    N += Exact;
                }
                if (N > Best_N) Best_N = N, Best_Lo = lo;
            }
            if (!Best_N || (Len == 1 && false)) break;
            // Below Min_Seed parts only for cases with no longer exact run (all such are left now).
            Group G;
            G.Seed_Lo = Best_Lo, G.Seed_Hi = Best_Lo + Len;
            for (std::size_t ci = 0; ci < Cases.size(); ci++) {
                if (Covered[ci]) continue;
                bool Exact = true;
                for (int t = Best_Lo; t < Best_Lo + Len && Exact; t++) Exact = !Cases[ci].C[t];
                if (Exact) G.Cases.push_back(static_cast<int>(ci)), Covered[ci] = 1, Left--;
            }
            int L = G.Seed_Lo, R = G.Seed_Hi;
            while (L > 0 || R < P) {
                long Lc = 1L << 40, Rc = 1L << 40;
                if (L > 0) { Lc = 0; for (int ci : G.Cases) Lc += Cases[ci].C[L - 1]; }
                if (R < P) { Rc = 0; for (int ci : G.Cases) Rc += Cases[ci].C[R]; }
                if (Lc <= Rc) G.Order.push_back(--L);
                else G.Order.push_back(R++);
            }
            Groups.push_back(G);
        }
    }
}

struct Searcher {
    const Bi_Index& I;
    const Packed_Genome* G;
    std::string_view Rd;
    const uint64_t* Packed;
    std::size_t L;
    std::size_t Cut[MAXP + 1];
    const Group* Gr;
    unsigned long Nodes = 0, Located = 0;
    std::vector<std::pair<unsigned long, int>>* Out;  // hits (pos or row marker) with mismatches
    std::vector<Bi_Range>* Range_Hits;

    int Part_Of(std::size_t Pos) const {
        int t = 0;
        while (Pos >= Cut[t + 1]) t++;
        return t;
    }

    // Node: range of P'[i, j); Step: index into the group's order (-1: still in the seed);
    // Cnt: mismatches per part; M: cases still possible.
    void Visit(Bi_Range R, std::size_t i, std::size_t j, int Step, std::array<unsigned char, MAXP> Cnt, Mask M) {
        if (R.Empty() || M.none()) return;
        Nodes++;
        if (G && R.Size() == 1) {
            const bool Fw = R.Fw_Start % 8 == 0;
            if (Fw || R.Rev_Start % 8 == 0) {
                Located++;
                const long S = Fw ? Located_Search::Sa(I.Fw(), R.Fw_Start) - long(i)
                                  : long(I.Text_Length()) - Located_Search::Sa(I.Rev(), R.Rev_Start) - long(j);
                if (S < 0 || static_cast<unsigned long>(S) + L > I.Text_Length()) return;
                std::size_t At[MAXP];
                const unsigned D = G->Diff_At(static_cast<unsigned long>(S), Packed, L, At, MAXP);
                if (D != static_cast<unsigned>(K)) return;
                unsigned char Dist[MAXP] = {};
                for (unsigned k = 0; k < D; k++) Dist[Part_Of(At[k])]++;
                for (int ci : Gr->Cases)
                    if (std::memcmp(Dist, Cases[ci].C, P) == 0) {
                        Out->push_back({static_cast<unsigned long>(S), K});
                        return;
                    }
                return;
            }
        }
        // The part being extended: the seed's (exact) until it is done, then the order's.
        int Part;
        bool Left;
        if (Step < 0) {
            if (i > Cut[Gr->Seed_Lo]) Part = Part_Of(i - 1), Left = true;
            else Step = 0;
        }
        if (Step >= 0) {
            // Past finished parts: keep the cases whose count there matches.
            for (;;) {
                if (Step == static_cast<int>(Gr->Order.size())) {
                    Range_Hits->push_back(R);
                    return;
                }
                const int T = Gr->Order[Step];
                const bool Done = T < Gr->Seed_Lo ? i <= Cut[T] : j >= Cut[T + 1];
                if (!Done) break;
                Mask N;
                for (int b = 0; b < static_cast<int>(Gr->Cases.size()); b++)
                    if (M[b] && Cases[Gr->Cases[b]].C[T] == Cnt[T]) N[b] = true;
                M = N;
                if (M.none()) return;
                Step++;
            }
            Part = Gr->Order[Step];
            Left = Part < Gr->Seed_Lo;
        }
        const std::size_t Pos = Left ? i - 1 : j;
        const int Read_Base = static_cast<unsigned char>(Rd[Pos]);
        // A mismatch here if some case still needs one in Part; a match if some case can
        // still finish Part with what is left of it.
        bool Can = false, Match = false;
        const std::size_t Left_In_Part = Left ? i - Cut[Part] - 1 : Cut[Part + 1] - j - 1;
        Mask M_Mis, M_Match;
        for (int b = 0; b < static_cast<int>(Gr->Cases.size()); b++) {
            if (!M[b]) continue;
            const int Need = Cases[Gr->Cases[b]].C[Part] - Cnt[Part];
            if (Need >= 1 && static_cast<std::size_t>(Need - 1) <= Left_In_Part) Can = true, M_Mis[b] = true;
            if (Need >= 0 && static_cast<std::size_t>(Need) <= Left_In_Part) Match = true, M_Match[b] = true;
        }
        auto Child = [&](const Bi_Range& C, int b) {
            std::array<unsigned char, MAXP> N = Cnt;
            if (b != Read_Base) N[Part]++;
            Visit(C, Left ? i - 1 : i, Left ? j : j + 1, Step, N, b != Read_Base ? M_Mis : M_Match);
        };
        if (!Can) {
            if (Match) Child(Left ? I.Extend_Left(R, Read_Base) : I.Extend_Right(R, Read_Base), Read_Base);
            return;
        }
        Bi_Range C[4];
        if (Left) I.Children_Left(R, C);
        else I.Children_Right(R, C);
        for (int b = 0; b < 4; b++) {
            if (C[b].Empty()) continue;
            if (b == Read_Base && !Match) continue;
            Child(C[b], b);
        }
    }

    void Run(const Group& Grp) {
        Gr = &Grp;
        Mask M;
        for (std::size_t b = 0; b < Grp.Cases.size(); b++) M[b] = true;
        std::array<unsigned char, MAXP> Cnt{};
        // The seed from its right end leftward, from the table if long enough.
        std::size_t i = Cut[Grp.Seed_Hi], j = Cut[Grp.Seed_Hi];
        Bi_Range R = I.Empty_Pattern();
        const std::size_t K_T = static_cast<std::size_t>(I.Lookup_Depth());
        if (K_T && Cut[Grp.Seed_Hi] - Cut[Grp.Seed_Lo] >= K_T) i -= K_T, R = I.Look_Up(Rd.data() + i);
        Visit(R, i, j, -1, Cnt, M);
    }
};

bool Next_Record(std::istream& In, std::string& Id, std::string& Seq, std::string& Line, bool& Have) {
    while (Have || std::getline(In, Line)) {
        Have = false;
        if (Line.empty()) continue;
        if (Line[0] == '@') {
            Id = Line.substr(1);
            std::string A, B;
            return std::getline(In, Seq) && std::getline(In, A) && std::getline(In, B);
        }
        if (Line[0] == '>') {
            Id = Line.substr(1);
            Seq.clear();
            while (std::getline(In, Line)) {
                if (!Line.empty() && (Line[0] == '>' || Line[0] == '@')) {
                    Have = true;
                    break;
                }
                Seq += Line;
            }
            return true;
        }
    }
    return false;
}

}  // namespace

int main(int argc, char** argv) {
    std::string Prefix, Reads;
    int Lookup = 0;
    bool Ref = false, Print = false;
    for (int a = 1; a < argc; a++) {
        if (!std::strcmp(argv[a], "-g")) Prefix = argv[++a];
        else if (!std::strcmp(argv[a], "-q")) Reads = argv[++a];
        else if (!std::strcmp(argv[a], "-p")) P = std::atoi(argv[++a]);
        else if (!std::strcmp(argv[a], "-k")) K = std::atoi(argv[++a]);
        else if (!std::strcmp(argv[a], "-lh")) LH_Need = std::atoi(argv[++a]);
        else if (!std::strcmp(argv[a], "--lookup")) Lookup = std::atoi(argv[++a]);
        else if (!std::strcmp(argv[a], "--ref")) Ref = true;
        else if (!std::strcmp(argv[a], "--hits")) Print = true;
    }
    std::vector<int> Cur(P);
    Enumerate(Cur, 0, LH_Need, K - LH_Need);
    Build_Groups();
    // Each case in exactly one group, by construction; print the groups.
    std::fprintf(stderr, "cases=%zu groups=%zu\n", Cases.size(), Groups.size());
    for (const Group& G : Groups) {
        std::fprintf(stderr, "  seed [%d,%d) order", G.Seed_Lo, G.Seed_Hi);
        for (int t : G.Order) std::fprintf(stderr, " %d", t);
        std::fprintf(stderr, "  cases=%zu\n", G.Cases.size());
    }
    Bi_Index Index(Prefix);
    std::unique_ptr<Packed_Genome> Genome;
    if (Ref) Genome = std::make_unique<Packed_Genome>(Prefix);
    Index.Build_Lookup(Lookup);
    std::ifstream In(Reads);
    std::string Id, Seq, Line;
    bool Have = false;
    unsigned long Reads_N = 0, Nodes = 0, Located = 0, Hits = 0;
    std::vector<uint64_t> Packed;
    const auto Start = std::chrono::steady_clock::now();
    while (Next_Record(In, Id, Seq, Line, Have)) {
        std::string C(Seq.size(), 0);
        bool Ok = !Seq.empty();
        for (std::size_t t = 0; t < Seq.size() && Ok; t++) {
            const int b = Base_Code(Seq[t]);
            Ok = b >= 0;
            C[t] = static_cast<char>(b);
        }
        Reads_N++;
        std::vector<std::pair<unsigned long, int>> Out;
        std::vector<Bi_Range> Range_Hits;
        if (Ok && static_cast<int>(C.size()) >= P) {
            Packed.resize((C.size() + 31) / 32 + 1);
            Pack_Read(C, Packed.data());
            Searcher Sr{Index, Genome.get(), C, Packed.data(), C.size(), {}, nullptr, 0, 0, &Out, &Range_Hits};
            for (int t = 0; t <= P; t++) Sr.Cut[t] = C.size() * t / P;
            for (const Group& G : Groups) Sr.Run(G);
            Nodes += Sr.Nodes;
            Located += Sr.Located;
        }
        for (const Bi_Range& R : Range_Hits)
            for (unsigned long r = R.Fw_Start; r <= R.Fw_End; r++) Out.push_back({Index.Locate_Fw(r), K});
        Hits += Out.size();
        if (Print) {
            std::sort(Out.begin(), Out.end());
            std::printf("%s\t%zu", Id.substr(0, Id.find_first_of(" \t")).c_str(), Out.size());
            for (auto& [p, k] : Out) std::printf("\t%lu:%d", p, k);
            std::printf("\n");
        }
    }
    const double Sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - Start).count();
    std::fprintf(stderr, "reads=%lu nodes_per_read=%.1f located_per_read=%.2f hits=%lu seconds=%.3f\n", Reads_N,
                 double(Nodes) / Reads_N, double(Located) / Reads_N, Hits, Sec);
    return 0;
}
