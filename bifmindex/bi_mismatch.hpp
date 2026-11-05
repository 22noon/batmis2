// bi_mismatch.hpp: one-mismatch search on the bidirectional index (batman_zero_min -n 1).
//
// The read P (L bases) splits into LH = P[0, L/2) and RH = P[L/2, L). A hit with at most
// one mismatch has LH exact (0|0, 0|1) or RH exact (1|0):
//   - stage A: LH by Extend_Right, then RH base by base along the exact path. Each node's
//     Children_Right gives, from one set of counts, the exact child and the three others;
//     those (one mismatch placed) are kept, and after the exact path each is extended
//     exactly to the read's end (0|1).
//   - stage B: RH by Extend_Left, then LH leftward the same way with Children_Left (1|0).
//     If stage A's LH search failed at depth d (P[0, d) does not occur), a 1|0 hit's
//     mismatch is before d, so only children before d are kept.
// Both stages give each hit's range in both indexes, so every hit is located in Fw
// (batman_zero_min reports 0|x in revfmi and 1|0 in fwfmi).
//
// With a packed genome, a single row whose SA value is sampled, in either index, gives
// the read's start: Fw row of P[i, j) at SA_fw - i, Rev row at n - SA_rev - j. The
// matched part occurs only there, so the node's search ends with one comparison of the
// whole read against the genome: in stage A's half search and exact path, a hit if it
// differs nowhere (exact) or once in RH (0|1); in stage B's, once in LH (1|0); in a kept
// child's extension, exactly at the child's mismatch.
//
// Mm_Search is the search of one read cut into steps that each read one thing from
// memory (a lookup entry, a step's occurrence rows, an SA sample, the genome), so that
// Batch_One_Mismatch can interleave reads; One_Mismatch runs it one read at a time.
#pragma once

#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "bi_index.hpp"
#include "bi_locate.hpp"
#include "bi_ref.hpp"

namespace bifm {

// A one-mismatch hit: the read with base Base at At (not the read's) has range R, or,
// once located, its one position Pos.
struct Mm_Hit {
    Bi_Range R;
    unsigned long Pos = Exact_Result::NO_POS;
    unsigned At = 0;
    unsigned char Base = 0;
    int Chain = -1;  // with Record, a wide unlocated hit's range chain (Mm_Result::Chains)
    unsigned char Mismatches = 1;  // 2: a Two_Mismatch hit (At, Base: its first mismatch)
    bool Located() const { return Pos != Exact_Result::NO_POS; }
    unsigned long Count() const { return Located() ? 1 : R.Size(); }
};

// A read's hits: the exact ones, then the 0|1 hits, then the 1|0 ones.
struct Mm_Result {
    Exact_Result Exact;
    std::vector<Mm_Hit> Hits;
    unsigned long Hit_Count = 0;  // locations in Hits
    std::vector<Range_Chain> Chains;
    int Exact_Chain = -1;
    void Clear() {
        Exact = {Bi_Range{1, 0, 0}};
        Hits.clear();
        Hit_Count = 0;
        Chains.clear();
        Exact_Chain = -1;
    }
};

// One read's one-mismatch search, step by step: Start, then Prefetch and Step until Step
// returns true.
class Mm_Search {
public:
    // G null: no located singletons.
    Mm_Search(const Bi_Index& Index, const Packed_Genome* Genome) : I(Index), G(Genome) {}

    // Record each step's range (stage A: Rev ranges of P'[0, j); stage B: Fw ranges of
    // P'[i, L)), and give unlocated hits of Locator::Chain_Min rows or more their chain.
    bool Record = false;

    // Starts the search of Read (Read_Packed: Pack_Read(Read), with a genome) into Result.
    // Max > 0 stops between the exact hits, the 0|1 and the 1|0 searches, and between the
    // kept children, once the read has that many. Returns true if it is already done.
    bool Start(std::string_view Read, const uint64_t* Read_Packed, unsigned long Max, Mm_Result& Result,
               Locate_Stats& Stats) {
        P = Read, Packed = Read_Packed, MaxHits = Max, Out = &Result, St = &Stats;
        Out->Clear();
        L = P.size(), LH = L / 2, LH_Fail = LH;
        Prefix = Suffix = {};
        Pending = {};
        return Begin_Stage(A_MAIN);
    }

    // Prefetches what the next Step reads.
    void Prefetch() const {
        switch (Ph) {
            case LOOKUP: I.Prefetch_Look_Up(P.data() + (Left() ? L - K() : 0)); break;
            case STEP: I.Prefetch_Step(R, Left()); break;
            case SA:
                if (Fw_Sampled) __builtin_prefetch(I.Fw()->saValue + R.Fw_Start / Located_Search::SA_INTERVAL);
                else __builtin_prefetch(I.Rev()->saValue + R.Rev_Start / Located_Search::SA_INTERVAL);
                break;
            case VERIFY: G->Prefetch(static_cast<unsigned long>(S), L); break;
        }
    }

    // One step; true when the read's search has finished.
    bool Step() {
        switch (Ph) {
            case LOOKUP:
                if (Left()) i = L - K(), R = I.Look_Up(P.data() + i);
                else j = K(), R = I.Look_Up(P.data());
                return After();
            case STEP:
                switch (Part) {
                    case A_MAIN:
                        if (j < LH) {
                            R = I.Extend_Right(R, P[j++]);
                        } else {
                            Bi_Range C[4];
                            I.Children_Right(R, C);
                            Keep(C, j);
                            R = C[static_cast<int>(P[j++])];
                        }
                        break;
                    case B_MAIN:
                        if (i > LH) {
                            R = I.Extend_Left(R, P[--i]);
                        } else {
                            Bi_Range C[4];
                            I.Children_Left(R, C);
                            if (--i < LH_Fail) Keep(C, i);
                            R = C[static_cast<int>(P[i])];
                        }
                        break;
                    case A_KEPT: R = I.Extend_Right(R, P[j++]); break;
                    case B_KEPT: R = I.Extend_Left(R, P[--i]); break;
                }
                return After();
            case SA:
                S = Fw_Sampled ? Located_Search::Sa(I.Fw(), R.Fw_Start) - long(i)
                               : long(I.Text_Length()) - Located_Search::Sa(I.Rev(), R.Rev_Start) - long(j);
                if (S < 0 || static_cast<unsigned long>(S) + L > I.Text_Length()) {
                    St->Off_Genome++;
                    Learn(-1);
                    return End_Node();
                }
                Ph = VERIFY;
                return false;
            case VERIFY: {
                std::size_t First;
                const unsigned D = G->Diff(static_cast<unsigned long>(S), Packed, L, First);
                const unsigned long Pos = static_cast<unsigned long>(S);
                Learn(S);
                bool Hit = false;
                switch (Part) {
                    case A_MAIN:  // exact, or one mismatch in RH (one in LH is stage B's)
                        if (D == 0) Out->Exact.Pos = Pos, Out->Exact.R = {0, 0, 0}, Hit = true;
                        else if (D == 1 && First >= LH) Pending = {Bi_Range{0, 0, 0}, Pos, unsigned(First), 0}, Hit = true;
                        break;
                    case B_MAIN:  // one mismatch in LH
                        if (D == 1 && First < LH) Pending = {Bi_Range{0, 0, 0}, Pos, unsigned(First), 0}, Hit = true;
                        break;
                    default:  // the kept child's mismatch alone
                        if (D == 1) Add({Bi_Range{0, 0, 0}, Pos, Kept[k].At, Kept[k].Base}), Hit = true;
                }
                Hit ? St->Moving++ : St->Mismatched++;
                return End_Node();
            }
        }
        return true;
    }

    // What a stage's main path learned of the read's prefix P[0, Len) (stage A) or suffix
    // P[L - Len, L) (stage B): it does not occur, or fits no read there (Pos -1), or occurs
    // only at the read start Pos. Len 0: nothing. Two_Mismatch skips the searches it
    // settles.
    struct Anchor {
        std::size_t Len = 0;
        long Pos = -1;
    };
    Anchor Prefix, Suffix;

private:
    // Records what the main path's node ending now (R empty, or located at S) shows.
    void Learn(long At) {
        if (Part == A_MAIN) Prefix = {j, At};
        else if (Part == B_MAIN) Suffix = {L - i, At};
    }

    // A_MAIN: stage A's LH search and exact path (R is P[0, j)'s range); A_KEPT: a kept
    // child extended right. B_MAIN, B_KEPT: stage B's (R is P[i, L)'s).
    enum Part_Kind { A_MAIN, A_KEPT, B_MAIN, B_KEPT };
    enum Phase { LOOKUP, STEP, SA, VERIFY };

    // A child off the exact path: P with Base at At has range R (P[0, At] in stage A,
    // P[At, L) in stage B).
    struct Kept_Child {
        Bi_Range R;
        unsigned At;
        unsigned char Base;
    };

    bool Left() const { return Part == B_MAIN || Part == B_KEPT; }
    std::size_t K() const { return static_cast<std::size_t>(I.Lookup_Depth()); }
    bool Full() const { return MaxHits && Out->Exact.Count() + Out->Hit_Count >= MaxHits; }
    void Add(const Mm_Hit& H) {
        Out->Hits.push_back(H);
        Out->Hit_Count += H.Count();
    }
    void Keep(const Bi_Range C[4], std::size_t At) {
        for (int b = 0; b < 4; b++)
            if (b != P[At] && !C[b].Empty()) Kept.push_back({C[b], unsigned(At), static_cast<unsigned char>(b)});
    }

    // Starts stage A or B's half search, from the table if the half is long enough.
    bool Begin_Stage(Part_Kind Main) {
        Part = Main;
        Kept.clear();
        R = I.Empty_Pattern();
        i = j = Main == A_MAIN ? 0 : L;
        Path_Low = 1;
        if (K() && (Main == A_MAIN ? LH : L - LH) >= K()) {
            Path_Low = static_cast<unsigned>(K());
            Ph = LOOKUP;
            return false;
        }
        return After();
    }

    // After R changes: the node ends (empty, a sampled single row, the read's end), or the
    // next step. Returns true when the read's search has finished.
    bool After() {
        if (R.Empty()) {
            if (Part == A_MAIN && j <= LH) LH_Fail = j;  // with the table, P[0, K) already
            Learn(-1);
            return End_Node();
        }
        if (Part == B_MAIN && i == 0) return End_Node();  // exact: stage A's
        if (Record) Note();
        if (G && R.Size() == 1) {
            Fw_Sampled = R.Fw_Start % Located_Search::SA_INTERVAL == 0;
            if (Fw_Sampled || R.Rev_Start % Located_Search::SA_INTERVAL == 0) {
                Ph = SA;
                return false;
            }
        }
        if (Left() ? i == 0 : j == L) {
            const int Ch = Record && R.Size() >= Locator::Chain_Min ? Chain_Of() : -1;
            if (Part == A_MAIN) Out->Exact.R = R, Out->Exact_Chain = Ch;
            else Add({R, Exact_Result::NO_POS, Kept[k].At, Kept[k].Base, Ch});
            return End_Node();
        }
        Ph = STEP;
        return false;
    }

    // Records R at its level (bases matched) in the main path's or the kept extension's
    // chain, in the index the stage extends.
    void Note() {
        const unsigned Level = static_cast<unsigned>(j - i);
        const bool Main = Part == A_MAIN || Part == B_MAIN;
        std::vector<unsigned long>& S = Main ? Path_S : Ext_S;
        std::vector<unsigned long>& E = Main ? Path_E : Ext_E;
        if (S.size() < L + 1) S.resize(L + 1), E.resize(L + 1);
        if (Left()) S[Level] = R.Fw_Start, E[Level] = R.Fw_End;
        else S[Level] = R.Rev_Start, E[Level] = R.Rev_Start + R.Size() - 1;
    }

    // The chain of the hit R ends (the whole read): the main path's levels, then for a
    // kept child its own from the child's level. Returns its index in Out->Chains.
    int Chain_Of() {
        Range_Chain& Ch = Out->Chains.emplace_back();
        Ch.Rev = !Left();
        Ch.Low = Path_Low;
        Ch.S.assign(L + 1, 0);
        Ch.E.assign(L + 1, 0);
        unsigned Split = static_cast<unsigned>(L) + 1;  // first level from the extension
        if (Part == A_KEPT || Part == B_KEPT) {
            Ch.Substitute(Kept[k].At, Kept[k].Base);
            Split = Part == A_KEPT ? Kept[k].At + 1 : static_cast<unsigned>(L) - Kept[k].At;
        }
        for (unsigned v = Path_Low; v <= L; v++) {
            Ch.S[v] = v < Split ? Path_S[v] : Ext_S[v];
            Ch.E[v] = v < Split ? Path_E[v] : Ext_E[v];
        }
        return static_cast<int>(Out->Chains.size() - 1);
    }

    // The node has ended: its stage's next kept child, the next stage, or the end.
    bool End_Node() {
        if (Part == A_MAIN || Part == B_MAIN) {
            Part = Part == A_MAIN ? A_KEPT : B_KEPT;
            k = 0;
        } else {
            k++;
        }
        if (k < Kept.size() && !Full()) {
            R = Kept[k].R;
            if (Part == A_KEPT) i = 0, j = Kept[k].At + 1;
            else i = Kept[k].At, j = L;
            return After();
        }
        // The stage's kept children are done: then its located exact-path hit, if any.
        if (Pending.Located()) {
            if (!Full()) Add(Pending);
            Pending = {};
        }
        if (Part == A_KEPT && !Full()) return Begin_Stage(B_MAIN);
        return true;
    }

    const Bi_Index& I;
    const Packed_Genome* G;
    std::string_view P;
    const uint64_t* Packed = nullptr;
    unsigned long MaxHits = 0;
    Mm_Result* Out = nullptr;
    Locate_Stats* St = nullptr;
    std::size_t L = 0, LH = 0, LH_Fail = 0;
    Part_Kind Part = A_MAIN;
    Phase Ph = STEP;
    Bi_Range R;
    std::size_t i = 0, j = 0;  // R is the range of P[i, j) (changed at the kept child's At)
    std::vector<Kept_Child> Kept;
    std::size_t k = 0;  // the kept child being extended
    // Record: the main path's and the kept extension's ranges by level; Path_Low the
    // first recorded (the lookup depth, or 1).
    std::vector<unsigned long> Path_S, Path_E, Ext_S, Ext_E;
    unsigned Path_Low = 1;
    bool Fw_Sampled = false;
    long S = 0;      // SA, VERIFY: the read's start
    Mm_Hit Pending;  // a located hit of the stage's exact path, added after its kept children
};

// The search one read at a time: Mm_Search's work as a plain loop, the same steps in the
// same order (the same hits, chains and Prefix / Suffix), without the cost of stopping
// after each step (measured 9-14% of the search; README.md, "Two-mismatch search").
class One_Mismatch {
public:
    One_Mismatch(const Bi_Index& Index, const Packed_Genome* Genome, bool Record_Chains = false)
        : I(Index), G(Genome), Record(Record_Chains) {
        if (Genome && (Index.Fw()->saInterval != Located_Search::SA_INTERVAL ||
                       Index.Rev()->saInterval != Located_Search::SA_INTERVAL))
            throw std::runtime_error("the SA sample interval is not 8");
    }

    // As Mm_Search::Start and its steps.
    void Search(std::string_view Read, const uint64_t* Read_Packed, unsigned long MaxHits, Mm_Result& Result, Locate_Stats& Stats) {
        P = Read, Packed = Read_Packed, Out = &Result, St = &Stats;
        Out->Clear();
        L = P.size();
        const std::size_t LH = L / 2, K = static_cast<std::size_t>(I.Lookup_Depth());
        std::size_t LH_Fail = LH;
        Prefix = Suffix = {};
        Mm_Hit Pending;
        auto Full = [&] { return MaxHits && Out->Exact.Count() + Out->Hit_Count >= MaxHits; };
        if (Record && Path_S.size() < L + 1) Path_S.resize(L + 1), Path_E.resize(L + 1), Ext_S.resize(L + 1), Ext_E.resize(L + 1);

        // Stage A: LH exact, then the exact path over RH, keeping the other children.
        Kept.clear();
        Bi_Range R = I.Empty_Pattern();
        std::size_t j = 0;  // R is P[0, j)'s range
        Path_Low = 1;
        if (K && LH >= K) R = I.Look_Up(P.data()), j = K, Path_Low = static_cast<unsigned>(K);
        for (;;) {
            if (R.Empty()) {
                if (j <= LH) LH_Fail = j;
                Prefix = {j, -1};
                break;
            }
            if (Record) Note(Path_S, Path_E, j, R, false);
            if (long S; G && R.Size() == 1 && Sampled(R, 0, j, S)) {
                if (!Fits(S)) {
                    St->Off_Genome++;
                    Prefix = {j, -1};
                    break;
                }
                Prefix = {j, S};
                std::size_t First;
                const unsigned D = G->Diff(static_cast<unsigned long>(S), Packed, L, First);
                bool Hit = false;
                if (D == 0) Out->Exact.Pos = static_cast<unsigned long>(S), Out->Exact.R = {0, 0, 0}, Hit = true;
                else if (D == 1 && First >= LH) Pending = {Bi_Range{0, 0, 0}, static_cast<unsigned long>(S), unsigned(First), 0}, Hit = true;
                Hit ? St->Moving++ : St->Mismatched++;
                break;
            }
            if (j == L) {
                Out->Exact.R = R;
                Out->Exact_Chain = Record && R.Size() >= Locator::Chain_Min ? Chain_Of(false, nullptr) : -1;
                break;
            }
            if (j < LH) {
                R = I.Extend_Right(R, P[j++]);
                continue;
            }
            Bi_Range C[4];
            I.Children_Right(R, C);
            Keep(C, j);
            R = C[static_cast<int>(P[j++])];
        }
        for (const Kept_Child& Ch : Kept) {
            if (Full()) break;
            Extend_Kept(Ch, false);
        }
        if (Pending.Located()) {
            if (!Full()) Add(Pending);
            Pending = {};
        }
        if (Full()) return;

        // Stage B: RH exact, then the exact path over LH leftward, keeping children before
        // LH_Fail.
        Kept.clear();
        R = I.Empty_Pattern();
        std::size_t i = L;  // R is P[i, L)'s range
        Path_Low = 1;
        if (K && L - LH >= K) i = L - K, R = I.Look_Up(P.data() + i), Path_Low = static_cast<unsigned>(K);
        for (;;) {
            if (R.Empty()) {
                Suffix = {L - i, -1};
                break;
            }
            if (i == 0) break;  // exact: stage A's
            if (Record) Note(Path_S, Path_E, L - i, R, true);
            if (long S; G && R.Size() == 1 && Sampled(R, i, L, S)) {
                if (!Fits(S)) {
                    St->Off_Genome++;
                    Suffix = {L - i, -1};
                    break;
                }
                Suffix = {L - i, S};
                std::size_t First;
                const unsigned D = G->Diff(static_cast<unsigned long>(S), Packed, L, First);
                const bool Hit = D == 1 && First < LH;
                if (Hit) Pending = {Bi_Range{0, 0, 0}, static_cast<unsigned long>(S), unsigned(First), 0};
                Hit ? St->Moving++ : St->Mismatched++;
                break;
            }
            if (i > LH) {
                R = I.Extend_Left(R, P[--i]);
                continue;
            }
            Bi_Range C[4];
            I.Children_Left(R, C);
            if (--i < LH_Fail) Keep(C, i);
            R = C[static_cast<int>(P[i])];
        }
        for (const Kept_Child& Ch : Kept) {
            if (Full()) break;
            Extend_Kept(Ch, true);
        }
        if (Pending.Located() && !Full()) Add(Pending);
    }

    // As Mm_Search's: what the stages' main paths learned (Two_Mismatch).
    Mm_Search::Anchor Prefix, Suffix;

private:
    // A child off the exact path: P with Base at At has range R (P[0, At] in stage A,
    // P[At, L) in stage B).
    struct Kept_Child {
        Bi_Range R;
        unsigned At;
        unsigned char Base;
    };

    void Keep(const Bi_Range C[4], std::size_t At) {
        for (int b = 0; b < 4; b++)
            if (b != P[At] && !C[b].Empty()) Kept.push_back({C[b], unsigned(At), static_cast<unsigned char>(b)});
    }
    void Add(const Mm_Hit& H) {
        Out->Hits.push_back(H);
        Out->Hit_Count += H.Count();
    }
    bool Fits(long S) const { return S >= 0 && static_cast<unsigned long>(S) + L <= I.Text_Length(); }
    // With R one row (the range of P[i, j)): whether a row of it is sampled, and S, the
    // read's start from it.
    bool Sampled(const Bi_Range& R, std::size_t i, std::size_t j, long& S) const {
        if (R.Fw_Start % Located_Search::SA_INTERVAL == 0) return S = Located_Search::Sa(I.Fw(), R.Fw_Start) - long(i), true;
        if (R.Rev_Start % Located_Search::SA_INTERVAL == 0)
            return S = long(I.Text_Length()) - Located_Search::Sa(I.Rev(), R.Rev_Start) - long(j), true;
        return false;
    }
    // Records R as level Level in the index the stage extends (Rev going right).
    static void Note(std::vector<unsigned long>& S, std::vector<unsigned long>& E, std::size_t Level, const Bi_Range& R, bool Left) {
        if (Left) S[Level] = R.Fw_Start, E[Level] = R.Fw_End;
        else S[Level] = R.Rev_Start, E[Level] = R.Rev_Start + R.Size() - 1;
    }

    // Kept child K extended exactly to the read's end (Left: its start); the hit, located
    // at a sampled single row (only K's mismatch) or as a range.
    void Extend_Kept(const Kept_Child& K, bool Left) {
        Bi_Range R = K.R;
        std::size_t i = Left ? K.At : 0, j = Left ? L : K.At + 1;  // R is the range of P[i, j), changed
        for (;;) {
            if (Record) Note(Ext_S, Ext_E, j - i, R, Left);
            if (long S; G && R.Size() == 1 && Sampled(R, i, j, S)) {
                if (!Fits(S)) return void(St->Off_Genome++);
                std::size_t First;
                if (G->Diff(static_cast<unsigned long>(S), Packed, L, First) != 1) return void(St->Mismatched++);
                St->Moving++;
                Add({Bi_Range{0, 0, 0}, static_cast<unsigned long>(S), K.At, K.Base});
                return;
            }
            if (Left ? i == 0 : j == L) {
                const int Ch = Record && R.Size() >= Locator::Chain_Min ? Chain_Of(Left, &K) : -1;
                Add({R, Exact_Result::NO_POS, K.At, K.Base, Ch});
                return;
            }
            R = Left ? I.Extend_Left(R, P[--i]) : I.Extend_Right(R, P[j++]);
            if (R.Empty()) return;
        }
    }

    // The chain of a whole-read hit: the main path's levels, then for a kept child its
    // own from the child's level (as Mm_Search::Chain_Of). Returns its index in Out->Chains.
    int Chain_Of(bool Left, const Kept_Child* K) {
        Range_Chain& Ch = Out->Chains.emplace_back();
        Ch.Rev = !Left;
        Ch.Low = Path_Low;
        Ch.S.assign(L + 1, 0);
        Ch.E.assign(L + 1, 0);
        unsigned Split = static_cast<unsigned>(L) + 1;
        if (K) {
            Ch.Substitute(K->At, K->Base);
            Split = Left ? static_cast<unsigned>(L) - K->At : K->At + 1;
        }
        for (unsigned v = Path_Low; v <= L; v++) {
            Ch.S[v] = v < Split ? Path_S[v] : Ext_S[v];
            Ch.E[v] = v < Split ? Path_E[v] : Ext_E[v];
        }
        return static_cast<int>(Out->Chains.size() - 1);
    }

    const Bi_Index& I;
    const Packed_Genome* G;
    const bool Record;
    std::string_view P;
    const uint64_t* Packed = nullptr;
    Mm_Result* Out = nullptr;
    Locate_Stats* St = nullptr;
    std::size_t L = 0;
    std::vector<Kept_Child> Kept;
    std::vector<unsigned long> Path_S, Path_E, Ext_S, Ext_E;
    unsigned Path_Low = 1;
};

// Reads Codes (base codes; those with !Valid are skipped, with no hits) searched B at a
// time: each round prefetches what every lane's next step reads, then steps every lane,
// so their cache misses overlap; a finished lane takes the next read. Results[r] gets read
// r's hits, as One_Mismatch::Search.
inline void Batch_One_Mismatch(const Bi_Index& I, const Packed_Genome* G, const std::vector<std::string>& Codes,
                               const std::vector<unsigned char>& Valid, unsigned B, unsigned long MaxHits,
                               std::vector<Mm_Result>& Results, Locate_Stats& St, bool Record = false) {
    struct Lane {
        Mm_Search M;
        std::vector<uint64_t> Packed;
    };
    Results.resize(Codes.size());
    std::vector<Lane> Lanes(B ? B : 1, Lane{Mm_Search(I, G), {}});
    for (Lane& L : Lanes) L.M.Record = Record;
    std::size_t Next = 0;
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
            if (!L.M.Start(P, L.Packed.data(), MaxHits, Results[r], St)) return true;
        }
        return false;
    };
    // Lanes[Order[0, Active)] are searching (Mm_Search holds a reference, so lanes stay put).
    std::vector<unsigned> Order(Lanes.size());
    unsigned Active = 0;
    while (Active < Lanes.size() && Fill(Lanes[Active])) Order[Active] = Active, Active++;
    while (Active) {
        for (unsigned l = 0; l < Active; l++) Lanes[Order[l]].M.Prefetch();
        for (unsigned l = 0; l < Active;) {
            Lane& L = Lanes[Order[l]];
            if (!L.M.Step() || Fill(L)) {
                l++;
                continue;
            }
            Order[l] = Order[--Active];  // no reads left: close the gap
        }
    }
}

}  // namespace bifm
