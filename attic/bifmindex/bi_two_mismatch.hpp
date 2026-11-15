// bi_two_mismatch.hpp: two-mismatch search on the bidirectional index (batman_zero_min
// -n 2): one read at a time (Two_Mismatch) or B reads in lanes (Batch_Two_Mismatch).
//
// Pass 1 is the one-mismatch search (the exact and one-mismatch hits). Pass 2 finds the
// hits with exactly two mismatches, split as batman_zero_min's -n 2
// (Two_Mismatch_Thirds): with cut points a = L/3 < LH = L/2 < b = 2L/3, parts
// P0 = [0, a), P1 = [a, LH), P2 = [LH, b), P3 = [b, L), each hit falls in exactly one
// search (p the mismatch in LH, q in RH for 1|1):
//   A  (0|2): LH exact, then right over RH (two);
//   B  (2|0): RH exact, then left over LH (two);
//   C1 (1|1, p in P1): P0 exact, then right: P1 one, RH one;
//   C2 (1|1, p in P0, q in P2): P3 exact, then left: P2 one, P1 none, P0 one;
//   C3 (1|1, p in P0, q in P3): P1 P2 exact (left from b), then left over P0 (one), then
//      right over P3 (one).
// The halves seed both mismatches of one half; the thirds only the split ones. C3 turns
// around: batman_zero_min finds its seed in one index and re-walks it into the other
// (Mm_Start_Switch); here the same Bi_Range extends either way.
//
// A search is a depth-first walk: at a position that may take a mismatch, Children_*
// gives all four bases from one set of counts; elsewhere one Extend_* step. A part must
// end with its minimum count, and the walk stops once the mismatches left cannot be
// placed. With a packed genome, a single row whose SA value is sampled locates the node's
// only occurrence: a hit if the genome there differs from the read in exactly two places
// that fit the search's per-part counts (located singletons, as One_Mismatch).
// Pass 1's main paths end at a prefix (suffix) that does not occur or occurs at one read
// start: one inside a search's seed settles the search (its hits, if any, are there), and
// one that does not occur bounds where a mismatch must be (Bounds_Possible).
// The ranges of a node anchored at the read's start (i = 0) or end (j = L) are recorded
// by level, so a wide hit gets its range chain for locating (bi_locate.hpp).
#pragma once

#include <string_view>
#include <vector>

#include "bi_index.hpp"
#include "bi_locate.hpp"
#include "bi_mismatch.hpp"
#include "bi_ref.hpp"

namespace bifm {

// Pass 2's searches, nodes and tests, shared by the recursive walk (Two_Mismatch) and the
// stepped one (Mm2_Search).
class Pass2_Core {
public:
    // The searches run (bit k: A, B, C1, C2, C3), to compare; all by default.
    static inline unsigned Scheme_Mask = 31;
    // Settle searches from pass 1's prefix and suffix (off: walk them all, to compare).
    static inline bool Settle = true;
    // Prune by pass 1's failed prefix and suffix (Bounds_Possible; off: to compare).
    static inline bool Bound = true;
    unsigned long Settled = 0;                // searches settled by pass 1 (--timing)
    unsigned long Nodes = 0, Located_Nodes = 0;  // nodes entered, located (--timing)

protected:
    Pass2_Core(const Bi_Index& Index, const Packed_Genome* Genome, bool Record_Chains)
        : I(Index), G(Genome), Record(Record_Chains) {}

    // A search: the parts in order, the direction each extends, and each part's (P0 ..
    // P3) mismatch bounds; every search places exactly two.
    static constexpr int NP = 4, SEARCHES = 5;
    struct Scheme {
        int Order[NP];
        bool Left[NP];
        unsigned char Min[NP], Max[NP];
    };
    static constexpr bool R_ = false, L_ = true;
    static constexpr Scheme Schemes[SEARCHES] = {
        {{0, 1, 2, 3}, {R_, R_, R_, R_}, {0, 0, 0, 0}, {0, 0, 2, 2}},  // A
        {{3, 2, 1, 0}, {L_, L_, L_, L_}, {0, 0, 0, 0}, {2, 2, 0, 0}},  // B
        {{0, 1, 2, 3}, {R_, R_, R_, R_}, {0, 1, 0, 0}, {0, 1, 1, 1}},  // C1 (RH: one by the total)
        {{3, 2, 1, 0}, {L_, L_, L_, L_}, {1, 0, 1, 0}, {1, 0, 1, 0}},  // C2
        {{2, 1, 0, 3}, {L_, L_, L_, R_}, {1, 0, 0, 1}, {1, 0, 0, 1}},  // C3
    };

    // A node: R is the range of P'[i, j), the read with the mismatches placed so far.
    struct Node {
        Bi_Range R;
        std::size_t i, j;
        int Step;               // the part being extended (index into Order)
        unsigned char Cnt[NP];  // mismatches in each part
        unsigned Total;
        unsigned Sub_At[2];
        unsigned char Sub_Base[2];
        unsigned Rev_Low, Fw_Low;  // the first recorded chain levels (0: none yet)
    };

    // How Enter leaves a node: nothing more (pruned, a hit, or the end), its sampled single
    // row to locate (Locate_At), or its children to make (Children_Of, with Ex).
    enum Entered { NONE, LOCATE, EXPAND };
    struct Expansion {
        int T;          // the part
        bool Left;      // the direction
        std::size_t Pos;  // the read position extended
        bool Can, Must;   // a mismatch may / must go there
        int Read_Base;
    };

    // Sets up pass 2 for a read, after pass 1 (whose Prefix and Suffix are Pre, Suf).
    void Begin(std::string_view Read, const uint64_t* Read_Packed, unsigned long Max, Mm_Result& Result,
               Locate_Stats& St, const Mm_Search::Anchor& Pre, const Mm_Search::Anchor& Suf) {
        P = Read, Packed = Read_Packed, MaxHits = Max, Out = &Result, Stats = &St;
        L = P.size();
        Lo[0] = 0, Hi[0] = Lo[1] = L / 3, Hi[1] = Lo[2] = L / 2, Hi[2] = Lo[3] = 2 * L / 3, Hi[3] = L;
        if (Record && Rev_S.size() < L + 1) Rev_S.resize(L + 1), Rev_E.resize(L + 1), Fw_S.resize(L + 1), Fw_E.resize(L + 1);
        Prefix = Pre, Suffix = Suf;
        // A prefix or suffix that does not occur (or fits no read) bounds every hit: one
        // mismatch before Pre_Fail, one at or after Suf_Start (batman_zero_min's
        // LH_Fail / RH_Fail windows).
        Pre_Fail = Bound && Pre.Len && Pre.Pos < 0 ? Pre.Len : L;
        Suf_Start = Bound && Suf.Len && Suf.Pos < 0 ? L - Suf.Len : 0;
    }

    bool Full() const { return MaxHits && Out->Exact.Count() + Out->Hit_Count >= MaxHits; }

    // Whether pass 1 settles search k: its seed starts at the read's start (A, C1) or end
    // (B, C2) and holds the prefix (suffix) pass 1 ended at, so its hits can only be at
    // that one read start. Pos: that start, or -1 (none).
    bool Settles(int k, long& Pos) {
        const std::size_t Seed_Len[SEARCHES] = {Hi[1], L - Lo[2], Hi[0], L - Lo[3], 0};
        const Mm_Search::Anchor& An = k == 0 || k == 2 ? Prefix : Suffix;
        if (k == 4 || !Settle || !An.Len || An.Len > Seed_Len[k]) return false;
        Settled++;
        Pos = G ? An.Pos : -1;
        return true;
    }

    // Search k's root, at its first part's inner end; Lookup: P's index of the table entry
    // to start from (the seed is long enough), else npos.
    Node Root(int k, std::size_t& Lookup) const {
        const Scheme& Sc = Schemes[k];
        Node N{I.Empty_Pattern(), 0, 0, 0, {0, 0, 0, 0}, 0, {0, 0}, {0, 0}, 0, 0};
        const int First = Sc.Order[0];
        // The exact seed: the first parts in the first direction with no mismatch allowed.
        std::size_t Seed = 0;
        for (int s = 0; s < NP && Sc.Left[s] == Sc.Left[0] && Sc.Max[Sc.Order[s]] == 0; s++)
            Seed += Hi[Sc.Order[s]] - Lo[Sc.Order[s]];
        const std::size_t K = static_cast<std::size_t>(I.Lookup_Depth());
        N.i = N.j = Sc.Left[0] ? Hi[First] : Lo[First];
        Lookup = std::string_view::npos;
        if (K && Seed >= K) {
            if (Sc.Left[0]) N.i -= K;
            else N.j += K;
            Lookup = N.i;
        }
        return N;
    }

    // Takes node N of search Sc as far as it goes without memory: prunes it, records its
    // range, finds a sampled single row (LOCATE; Fw_Sampled), passes finished parts, makes
    // a whole-read node a hit, or works out its extension (EXPAND; Ex).
    Entered Enter(const Scheme& Sc, Node& N, bool& Fw_Sampled, Expansion& Ex) {
        if (N.R.Empty() || Full() || !Bounds_Possible(N)) return NONE;
        Nodes++;
        if (Record) {
            if (N.i == 0 && N.j > 0) {
                Rev_S[N.j] = N.R.Rev_Start, Rev_E[N.j] = N.R.Rev_Start + N.R.Size() - 1;
                if (!N.Rev_Low) N.Rev_Low = static_cast<unsigned>(N.j);
            }
            if (N.j == L && N.i < L) {
                Fw_S[L - N.i] = N.R.Fw_Start, Fw_E[L - N.i] = N.R.Fw_End;
                if (!N.Fw_Low) N.Fw_Low = static_cast<unsigned>(L - N.i);
            }
        }
        // A sampled single row: the node's only occurrence decides it.
        if (G && N.R.Size() == 1) {
            Fw_Sampled = N.R.Fw_Start % Located_Search::SA_INTERVAL == 0;
            if (Fw_Sampled || N.R.Rev_Start % Located_Search::SA_INTERVAL == 0) {
                Located_Nodes++;
                return LOCATE;
            }
        }
        // Past finished parts (each must have its minimum).
        for (;;) {
            if (N.Step == NP) {  // the whole read
                if (N.Total == 2) Add_Range_Hit(N);
                return NONE;
            }
            const int T = Sc.Order[N.Step];
            const bool Done = Sc.Left[N.Step] ? N.i <= Lo[T] : N.j >= Hi[T];  // the table may jump past it
            if (!Done) break;
            if (N.Cnt[T] < Sc.Min[T]) return NONE;
            N.Step++;
        }
        Ex.T = Sc.Order[N.Step];
        Ex.Left = Sc.Left[N.Step];
        Ex.Pos = Ex.Left ? N.i - 1 : N.j;
        // The mismatches still to place must fit in the positions left, and this part's
        // minimum in what is left of it.
        const std::size_t Remaining = L - (N.j - N.i), In_Part = Ex.Left ? N.i - Lo[Ex.T] : Hi[Ex.T] - N.j;
        if (2 - N.Total > Remaining || Sc.Min[Ex.T] > N.Cnt[Ex.T] + In_Part) return NONE;
        Ex.Can = N.Cnt[Ex.T] < Sc.Max[Ex.T] && N.Total < 2;
        Ex.Must = 2 - N.Total == Remaining || Sc.Min[Ex.T] == N.Cnt[Ex.T] + In_Part;  // no room for a match here
        if (Ex.Must && !Ex.Can) return NONE;
        Ex.Read_Base = static_cast<unsigned char>(P[Ex.Pos]);
        return EXPAND;
    }

    // N's children (one memory access: the step's rows), in base order, into C;
    // returns how many.
    int Children_Of(const Node& N, const Expansion& Ex, Node C[4]) const {
        auto Make = [&](Node& M, const Bi_Range& R, int b) {
            M = N;
            M.R = R;
            if (Ex.Left) M.i--;
            else M.j++;
            if (b != Ex.Read_Base) {
                M.Cnt[Ex.T]++;
                M.Sub_At[M.Total] = static_cast<unsigned>(Ex.Pos), M.Sub_Base[M.Total] = static_cast<unsigned char>(b);
                M.Total++;
            }
        };
        if (!Ex.Can) {
            Make(C[0], Ex.Left ? I.Extend_Left(N.R, Ex.Read_Base) : I.Extend_Right(N.R, Ex.Read_Base), Ex.Read_Base);
            return 1;
        }
        Bi_Range R[4];
        if (Ex.Left) I.Children_Left(N.R, R);
        else I.Children_Right(N.R, R);
        int n = 0;
        for (int b = 0; b < 4; b++)
            if (!R[b].Empty() && !(Ex.Must && b == Ex.Read_Base)) Make(C[n++], R[b], b);
        return n;
    }

    // The read start a located node gives (from its sampled row), and whether the read
    // fits there.
    long Start_Of(const Node& N, bool Fw_Sampled) const {
        return Fw_Sampled ? Located_Search::Sa(I.Fw(), N.R.Fw_Start) - long(N.i)
                          : long(I.Text_Length()) - Located_Search::Sa(I.Rev(), N.R.Rev_Start) - long(N.j);
    }
    bool Fits_Genome(long S) const { return S >= 0 && static_cast<unsigned long>(S) + L <= I.Text_Length(); }

    // A hit of Sc at S if the read there has exactly two mismatches that fit Sc.
    bool Check_At(const Scheme& Sc, long S) {
        std::size_t At[3];
        if (G->Diff_At(static_cast<unsigned long>(S), Packed, L, At, 3) != 2 || !Fits(Sc, At)) return false;
        Add({Bi_Range{0, 0, 0}, static_cast<unsigned long>(S), unsigned(At[0]), 0, -1, 2});
        return true;
    }
    // A located node at S: its hit, if any, with the located stats.
    void Located_At(const Scheme& Sc, long S) {
        if (!Fits_Genome(S)) return void(Stats->Off_Genome++);
        Check_At(Sc, S) ? Stats->Moving++ : Stats->Mismatched++;
    }

    // Whether N's hits can still have a mismatch before Pre_Fail and one at or after
    // Suf_Start: placed already, or a position left for it (positions [0, i) and [j, L)
    // are left, with Total < 2).
    bool Bounds_Possible(const Node& N) const {
        bool Low = false, High = false;
        for (unsigned s = 0; s < N.Total; s++) Low |= N.Sub_At[s] < Pre_Fail, High |= N.Sub_At[s] >= Suf_Start;
        const bool More = N.Total < 2;
        if (!Low && !(More && (N.i > 0 || N.j < Pre_Fail))) return false;
        if (!High && !(More && (N.j < L || N.i > Suf_Start))) return false;
        return true;
    }

    // Whether mismatches at At[0], At[1] fit the search's per-part counts.
    bool Fits(const Scheme& Sc, const std::size_t At[2]) const {
        unsigned Cnt[NP] = {0, 0, 0, 0};
        for (int k = 0; k < 2; k++) Cnt[At[k] < Hi[0] ? 0 : At[k] < Hi[1] ? 1 : At[k] < Hi[2] ? 2 : 3]++;
        for (int t = 0; t < NP; t++)
            if (Cnt[t] < Sc.Min[t] || Cnt[t] > Sc.Max[t]) return false;
        return true;
    }

    void Add(const Mm_Hit& H) {
        Out->Hits.push_back(H);
        Out->Hit_Count += H.Count();
    }

    // A whole-read node: a hit of R's rows, with its range chain if wide (the index
    // whose chain reaches lower).
    void Add_Range_Hit(const Node& N) {
        int Chain = -1;
        if (Record && N.R.Size() >= Locator::Chain_Min && (N.Rev_Low || N.Fw_Low)) {
            const bool Use_Rev = N.Rev_Low && (!N.Fw_Low || N.Rev_Low <= N.Fw_Low);
            Range_Chain& Ch = Out->Chains.emplace_back();
            Ch.Reset(L, Use_Rev, Use_Rev ? N.Rev_Low : N.Fw_Low);
            for (unsigned v = Ch.Low; v <= L; v++) {
                Ch.S[v] = Use_Rev ? Rev_S[v] : Fw_S[v];
                Ch.E[v] = Use_Rev ? Rev_E[v] : Fw_E[v];
            }
            for (unsigned s = 0; s < N.Total; s++) Ch.Substitute(N.Sub_At[s], N.Sub_Base[s]);
            Chain = static_cast<int>(Out->Chains.size() - 1);
        }
        Add({N.R, Exact_Result::NO_POS, N.Sub_At[0], N.Sub_Base[0], Chain, 2});
    }

    const Bi_Index& I;
    const Packed_Genome* G;
    const bool Record;
    std::string_view P;
    const uint64_t* Packed = nullptr;
    unsigned long MaxHits = 0;
    Mm_Result* Out = nullptr;
    Locate_Stats* Stats = nullptr;
    std::size_t L = 0, Lo[NP] = {0, 0, 0, 0}, Hi[NP] = {0, 0, 0, 0};
    std::size_t Pre_Fail = 0, Suf_Start = 0;  // Bounds_Possible
    Mm_Search::Anchor Prefix, Suffix;
    // Record: the current path's ranges by level, anchored at the start (Rev) or end (Fw).
    // A depth-first walk only ever overwrites levels deeper than the node it returns to.
    std::vector<unsigned long> Rev_S, Rev_E, Fw_S, Fw_E;
};

// The two-mismatch search one read at a time: pass 1 (One_Mismatch's plain loop), then
// pass 2 as a recursive walk.
class Two_Mismatch : public Pass2_Core {
public:
    // G null: no located singletons. Record: give wide unlocated hits their range chain.
    Two_Mismatch(const Bi_Index& Index, const Packed_Genome* Genome, bool Record_Chains = false)
        : Pass2_Core(Index, Genome, Record_Chains), Pass1(Index, Genome, Record_Chains) {}

    // P's hits with at most two mismatches into Out: pass 1's, then the two-mismatch ones
    // (Mm_Hit::Mismatches 2). MaxHits > 0 stops between the passes, the searches and the
    // nodes once the read has that many.
    void Search(std::string_view Read, const uint64_t* Read_Packed, unsigned long Max, Mm_Result& Result, Locate_Stats& St) {
        Pass1.Search(Read, Read_Packed, Max, Result, St);
        Begin(Read, Read_Packed, Max, Result, St, Pass1.Prefix, Pass1.Suffix);
        for (int k = 0; k < SEARCHES; k++) {
            if (Full()) return;
            if (!(Scheme_Mask >> k & 1)) continue;
            if (long Pos; Settles(k, Pos)) {
                if (Pos >= 0) Check_At(Schemes[k], Pos);
                continue;
            }
            std::size_t Lookup;
            Node N = Root(k, Lookup);
            if (Lookup != std::string_view::npos) N.R = I.Look_Up(P.data() + Lookup);
            Visit(Schemes[k], N);
        }
    }

private:
    void Visit(const Scheme& Sc, Node N) {
        bool Fw_Sampled;
        Expansion Ex;
        switch (Enter(Sc, N, Fw_Sampled, Ex)) {
            case NONE: return;
            case LOCATE: return Located_At(Sc, Start_Of(N, Fw_Sampled));
            case EXPAND: {
                Node C[4];
                const int n = Children_Of(N, Ex, C);
                for (int c = 0; c < n; c++) Visit(Sc, C[c]);
            }
        }
    }

    One_Mismatch Pass1;
};

// Pass 2 of one read, step by step (as Mm_Search): Start, then Prefetch and Step until
// Step returns true. The walk keeps its nodes on a stack, children pushed in reverse, so
// it visits them in Two_Mismatch's order, with the same hits.
class Mm2_Search : public Pass2_Core {
public:
    Mm2_Search(const Bi_Index& Index, const Packed_Genome* Genome, bool Record_Chains = false)
        : Pass2_Core(Index, Genome, Record_Chains) {}

    // Starts pass 2 of Read after its pass 1 (Pre, Suf: pass 1's Prefix and Suffix).
    // Returns true if it is already done.
    bool Start(std::string_view Read, const uint64_t* Read_Packed, unsigned long Max, Mm_Result& Result, Locate_Stats& St,
               const Mm_Search::Anchor& Pre, const Mm_Search::Anchor& Suf) {
        Begin(Read, Read_Packed, Max, Result, St, Pre, Suf);
        k = -1;
        Stack.clear();
        return Next_Search();
    }

    // Prefetches what the next Step reads.
    void Prefetch() const {
        switch (Ph) {
            case LOOKUP: I.Prefetch_Look_Up(P.data() + Lookup_At); break;
            case STEP: I.Prefetch_Step(Cur.R, Ex.Left); break;
            case SA:
                if (Fw_Sampled) __builtin_prefetch(I.Fw()->saValue + Cur.R.Fw_Start / Located_Search::SA_INTERVAL);
                else __builtin_prefetch(I.Rev()->saValue + Cur.R.Rev_Start / Located_Search::SA_INTERVAL);
                break;
            case VERIFY:
            case CHECK: G->Prefetch(static_cast<unsigned long>(S), L); break;
        }
    }

    // One step; true when the read's pass 2 has finished.
    bool Step() {
        switch (Ph) {
            case LOOKUP:
                Cur.R = I.Look_Up(P.data() + Lookup_At);
                return Enter_Cur();
            case STEP: {
                // The first child is entered at once (it would be the next popped); the
                // others wait on the stack, last first.
                Node C[4];
                const int n = Children_Of(Cur, Ex, C);
                if (!n) return Walk();
                for (int c = n; c-- > 1;) Stack.push_back(C[c]);
                Cur = C[0];
                return Enter_Cur();
            }
            case SA:
                S = Start_Of(Cur, Fw_Sampled);
                if (!Fits_Genome(S)) {
                    Stats->Off_Genome++;
                    return Walk();
                }
                Ph = VERIFY;
                return false;
            case VERIFY:
                Check_At(Schemes[k], S) ? Stats->Moving++ : Stats->Mismatched++;
                return Walk();
            case CHECK:
                Check_At(Schemes[k], S);
                return Next_Search();
        }
        return true;
    }

private:
    enum Phase { LOOKUP, STEP, SA, VERIFY, CHECK };

    // Enters Cur: its next phase if it needs memory, else on to the next node.
    bool Enter_Cur() {
        switch (Enter(Schemes[k], Cur, Fw_Sampled, Ex)) {
            case NONE: return Walk();
            case LOCATE: Ph = SA; return false;
            case EXPAND: Ph = STEP; return false;
        }
        return Walk();
    }

    // Pops nodes until one needs memory (its children, or its sampled row), or the
    // search's stack is empty (then the next search).
    bool Walk() {
        while (!Stack.empty()) {
            Cur = Stack.back();
            Stack.pop_back();
            switch (Enter(Schemes[k], Cur, Fw_Sampled, Ex)) {
                case NONE: break;
                case LOCATE: Ph = SA; return false;
                case EXPAND: Ph = STEP; return false;
            }
        }
        return Next_Search();
    }

    // The next search: settled (a genome check at most), from the table, or its root.
    bool Next_Search() {
        while (++k < SEARCHES) {
            if (Full()) return true;
            if (!(Scheme_Mask >> k & 1)) continue;
            if (long Pos; Settles(k, Pos)) {
                if (Pos < 0) continue;
                S = Pos;
                Ph = CHECK;
                return false;
            }
            Cur = Root(k, Lookup_At);
            if (Lookup_At != std::string_view::npos) {
                Ph = LOOKUP;
                return false;
            }
            return Enter_Cur();
        }
        return true;
    }

    int k = 0;  // the search
    std::vector<Node> Stack;
    Node Cur{};
    Expansion Ex{};
    Phase Ph = STEP;
    bool Fw_Sampled = false;
    long S = 0;  // SA, VERIFY, CHECK: the read start
    std::size_t Lookup_At = 0;
};

// Reads Codes (base codes; those with !Valid are skipped, with no hits) searched for two
// mismatches B at a time: each lane runs pass 1 (Mm_Search), then pass 2 (Mm2_Search);
// each round prefetches what every lane's next step reads, then steps every lane. Results[r]
// gets read r's hits, as Two_Mismatch::Search. Pass2 (if given) gets the lanes' pass 2
// counts.
inline void Batch_Two_Mismatch(const Bi_Index& I, const Packed_Genome* G, const std::vector<std::string>& Codes,
                               const std::vector<unsigned char>& Valid, unsigned B, unsigned long MaxHits,
                               std::vector<Mm_Result>& Results, Locate_Stats& St, bool Record = false,
                               unsigned long* Pass2_Counts = nullptr) {
    struct Lane {
        Mm_Search M1;
        Mm2_Search M2;
        bool Second;  // in pass 2
        std::size_t Read;
        std::vector<uint64_t> Packed;
    };
    Results.resize(Codes.size());
    std::vector<Lane> Lanes;
    Lanes.reserve(B ? B : 1);
    for (unsigned l = 0; l < (B ? B : 1); l++) {
        Lanes.push_back(Lane{Mm_Search(I, G), Mm2_Search(I, G, Record), false, 0, {}});
        Lanes.back().M1.Record = Record;
    }
    std::size_t Next = 0;
    // Pass 1 has finished in L: starts pass 2; true if that is done too.
    auto To_Pass2 = [&](Lane& L) {
        L.Second = true;
        return L.M2.Start(Codes[L.Read], L.Packed.data(), MaxHits, Results[L.Read], St, L.M1.Prefix, L.M1.Suffix);
    };
    // Starts the next valid read in L that needs a step; false if none is left.
    auto Fill = [&](Lane& L) {
        while (Next < Codes.size()) {
            const std::size_t r = Next++;
            Results[r].Clear();
            if (!Valid[r]) continue;
            const std::string& P = Codes[r];
            if (G) {
                L.Packed.resize((P.size() + 31) / 32);
                Pack_Read(P, L.Packed.data());
            }
            L.Read = r;
            L.Second = false;
            if (!L.M1.Start(P, L.Packed.data(), MaxHits, Results[r], St) || !To_Pass2(L)) return true;
        }
        return false;
    };
    // Lanes[Order[0, Active)] are searching (the searches hold references, so lanes stay put).
    std::vector<unsigned> Order(Lanes.size());
    unsigned Active = 0;
    while (Active < Lanes.size() && Fill(Lanes[Active])) Order[Active] = Active, Active++;
    while (Active) {
        for (unsigned l = 0; l < Active; l++) {
            const Lane& L = Lanes[Order[l]];
            L.Second ? L.M2.Prefetch() : L.M1.Prefetch();
        }
        for (unsigned l = 0; l < Active;) {
            Lane& L = Lanes[Order[l]];
            const bool Done = L.Second ? L.M2.Step() : L.M1.Step() && To_Pass2(L);
            if (!Done || Fill(L)) {
                l++;
                continue;
            }
            Order[l] = Order[--Active];  // no reads left: close the gap
        }
    }
    if (Pass2_Counts)
        for (const Lane& L : Lanes) Pass2_Counts[0] += L.M2.Nodes, Pass2_Counts[1] += L.M2.Located_Nodes, Pass2_Counts[2] += L.M2.Settled;
}

}  // namespace bifm
