// bi_locate.hpp: hit rows to text positions with the bidirectional index.
//
// A row's position comes from an LF walk to a sampled row (row % 8 == 0), about 7
// dependent cache misses (BWTSaValue). Two things help here:
//   - Pair walks. An occurrence of P[i, j) with one row in each index (a single-row
//     range) has two walks: LF in Fw moves its start left, LF in Rev moves its end
//     right. Stepped together, the first to reach a sample wins: about 3.3 rounds of two
//     independent misses instead of 7 serial ones. A narrow range is first split into
//     such pairs by Children_Left: a child's Rev rows are a sub-range of its parent's
//     (Rev rows do not move under Extend_Left), so a single-row child is a matched pair.
//   - The range chain. The search's Extend_Left steps give the Fw ranges of P[i, L) and
//     its Extend_Right steps the Rev ranges of P[0, j): the chain A_1 .. A_L that
//     batman_zero_min's Chain_Decode rebuilds (2 Occ a step) before walking it back.
//     Recorded during the search it costs nothing; walking it back is Chain_Decode's:
//     each hit row has an ancestor in every A_k (LF keeps order), and a sampled
//     ancestor gives the row's position with no walk (Chain_Backtrack).
// Locator queues walks (Add_*), resolving chains at once, then steps every walk a round
// at a time with prefetching (Run), so the walks of all reads overlap.
#pragma once

#include <algorithm>
#include <string_view>
#include <vector>

#include "bi_index.hpp"

namespace bifm {

// Which route resolved each position (bi_* --timing).
struct Locate_Counts {
    unsigned long Pairs = 0, Single_Walks = 0, Rounds = 0, Split_Nodes = 0, Chain_Hits = 0, Chain_Rows = 0,
                  Chain_Resolved = 0, Chain_Rebuilt = 0, Dollar = 0,
                  Walk_Steps = 0,  // LF steps (rounds) of every walk, summed
                  Serial_Rows = 0;  // rows of failed-trial chains, by the library's walk
};

class Locator {
public:
    static constexpr unsigned long NO_ROW = ~0UL;
    static constexpr unsigned long SA_INTERVAL = 8;
    static inline unsigned long Chain_Min = 24;     // hits this wide use their chain
    static inline unsigned Split_Depth = 8;         // narrow ranges split this deep at most
    static inline bool Split = false;               // narrow ranges split into pairs (measured no faster)
    static inline std::size_t Window = 1024;        // walks stepped at a time (Run)
    static inline bool Both_Ends = true;            // a pair walks both rows (off: Fw only, to compare)
    static inline unsigned long Chain_Scan_Per_Row = 4, Chain_Idle_Scans = 12, Chain_Trial = 256;

    explicit Locator(const Bi_Index& Index) : I(Index), N(Index.Text_Length()) {}

    // A walk to *Out: from Fw row Fw_Row (start = SA_fw - Fw_Off) and/or Rev row Rev_Row
    // (start = n - SA_rev - Rev_End); NO_ROW for a missing one.
    void Add_Pair(unsigned long Fw_Row, unsigned long Rev_Row, long Fw_Off, long Rev_End, unsigned long* Out) {
        if (!Both_Ends && Fw_Row != NO_ROW) Rev_Row = NO_ROW;
        Walk W{Fw_Row, Rev_Row, Fw_Off, Rev_End, 0, LF, Out};
        if (Next_Phase(W)) return;
        Walks.push_back(W);
        Fw_Row != NO_ROW && Rev_Row != NO_ROW ? C.Pairs++ : C.Single_Walks++;
    }

    // Every occurrence of R, the range of P[i, j) (Out[0, R.Size())): one pair, split
    // pairs (narrow, Split), or Fw walks.
    void Add_Range(const Bi_Range& R, long i, long j, unsigned long* Out) {
        if (R.Size() == 1) return Add_Pair(R.Fw_Start, R.Rev_Start, i, j, Out);
        if (Split && R.Size() < Chain_Min) return Split_Range(R, i, j, Out);
        for (unsigned long r = 0; r < R.Size(); r++) Add_Pair(R.Fw_Start + r, NO_ROW, i, j, Out + r);
    }

    // The first Count rows of a whole-read hit P' (L bases) at level L of chain Ch: those
    // whose chain ancestors are sampled at once, the rest as walks (Out[0, Count)).
    // Codes: the read. Fw_All: the hit's Fw rows, if Count is all of them (a failed trial
    // then walks those, as the library would).
    void Add_Chain(Range_Chain& Ch, std::string_view Codes, unsigned long Count, unsigned long* Out,
                   const Range* Fw_All = nullptr) {
        Chain_Backtrack(Ch, Codes, Count, Out, Fw_All);
    }

    // Steps every queued walk until all are done.
    void Run() {
        const BWT *Fw = I.Fw(), *Rv = I.Rev();
        // Window walks at a time, refilled from the queue as they finish: enough to overlap
        // their misses, few enough that walks stalled in a homopolymer run (whose rows
        // stay in cache) keep their locality.
        std::vector<Walk> Active;
        Active.reserve(Window);
        std::size_t Next = 0;
        while (Next < Walks.size() || !Active.empty()) {
            while (Active.size() < Window && Next < Walks.size()) Active.push_back(Walks[Next++]);
            C.Rounds++;
            for (const Walk& W : Active) {
                switch (W.Ph) {
                    case LF:
                        if (W.Fw != NO_ROW) Bi_Index::Prefetch_Row(Fw, W.Fw);
                        if (W.Rv != NO_ROW) Bi_Index::Prefetch_Row(Rv, W.Rv);
                        break;
                    case FW_SA: __builtin_prefetch(Fw->saValue + W.Fw / SA_INTERVAL); break;
                    case REV_SA: __builtin_prefetch(Rv->saValue + W.Rv / SA_INTERVAL); break;
                }
            }
            std::size_t Kept = 0;
            for (Walk& W : Active) {
                bool Done = true;
                switch (W.Ph) {
                    case LF:
                        if (W.Fw != NO_ROW) W.Fw = Lf(Fw, W.Fw);
                        if (W.Rv != NO_ROW) W.Rv = Lf(Rv, W.Rv);
                        W.Steps++;
                        Done = Next_Phase(W);
                        break;
                    case FW_SA: *W.Out = Fw->saValue[W.Fw / SA_INTERVAL] + W.Steps - W.Fw_Off; break;
                    case REV_SA: *W.Out = N - (Rv->saValue[W.Rv / SA_INTERVAL] + W.Steps) - W.Rev_End; break;
                }
                if (!Done) Active[Kept++] = W;
                else C.Walk_Steps += W.Steps;
            }
            Active.resize(Kept);
        }
        Walks.clear();
    }

    const Locate_Counts& Counts() const { return C; }

private:
    enum Phase { LF, FW_SA, REV_SA };
    struct Node {  // Split_Range: a range and how many bases it has been extended left
        Bi_Range R;
        unsigned Depth;
    };
    struct Row_Of {  // Chain_Backtrack
        unsigned long Row;  // the hit row's ancestor in the current level
        unsigned long J;    // the hit row's offset
    };
    struct Walk {
        unsigned long Fw, Rv;  // the rows, or NO_ROW
        long Fw_Off, Rev_End;
        unsigned long Steps;   // LF steps taken (both rows step together)
        Phase Ph;
        unsigned long* Out;
    };

    // LF of a row that is not the $ row.
    static unsigned long Lf(const BWT* B, unsigned long Row) {
        const unsigned c = Bi_Index::Bwt_Char(B, Row);
        return B->cumulativeFreq[c] + Bi_Index::Occ_Base(B, Row, static_cast<int>(c)) + 1;
    }

    // W's next phase once its rows have moved: done at a $ row (SA 0), a sampled row's
    // SA read, or another LF step. Returns true when done.
    bool Next_Phase(Walk& W) {
        const BWT *Fw = I.Fw(), *Rv = I.Rev();
        if (W.Fw != NO_ROW && W.Fw == Fw->inverseSa0) {
            *W.Out = W.Steps - W.Fw_Off;
            C.Dollar++;
            return true;
        }
        if (W.Rv != NO_ROW && W.Rv == Rv->inverseSa0) {
            *W.Out = N - W.Steps - W.Rev_End;
            C.Dollar++;
            return true;
        }
        if (W.Fw != NO_ROW && W.Fw % SA_INTERVAL == 0) W.Ph = FW_SA;
        else if (W.Rv != NO_ROW && W.Rv % SA_INTERVAL == 0) W.Ph = REV_SA;
        else W.Ph = LF;
        return false;
    }

    // Splits R (P[i, j)) by Children_Left into single-row children, matched pairs, down
    // to Split_Depth; wider children left there walk their Fw rows. Out gets the
    // occurrences in an order of its own.
    void Split_Range(const Bi_Range& R, long i, long j, unsigned long* Out) {
        std::vector<Node>& Stack = Split_Stack;
        Stack.clear();
        Stack.push_back({R, 0});
        unsigned long k = 0;
        const BWT* Fw = I.Fw();
        while (!Stack.empty()) {
            const Node Nd = Stack.back();
            Stack.pop_back();
            C.Split_Nodes++;
            const long Off = i - long(Nd.Depth);
            // The occurrence at T's start has no preceding base: Fw's $ row, SA 0.
            if (Fw->inverseSa0 >= Nd.R.Fw_Start && Fw->inverseSa0 <= Nd.R.Fw_End) Out[k++] = static_cast<unsigned long>(-Off), C.Dollar++;
            Bi_Range Ch[4];
            I.Children_Left(Nd.R, Ch);
            for (int c = 0; c < 4; c++) {
                if (Ch[c].Empty()) continue;
                if (Ch[c].Size() == 1) Add_Pair(Ch[c].Fw_Start, Ch[c].Rev_Start, Off - 1, j, Out + k++);
                else if (Nd.Depth + 1 < Split_Depth) Stack.push_back({Ch[c], Nd.Depth + 1});
                else for (unsigned long r = 0; r < Ch[c].Size(); r++) Add_Pair(Ch[c].Fw_Start + r, NO_ROW, Off - 1, j, Out + k++);
            }
        }
    }

    // The character the chain's index consumes at level k (1-based).
    static unsigned Consumed(const Range_Chain& Ch, std::string_view Codes, unsigned k) {
        const std::size_t L = Codes.size();
        const std::size_t p = Ch.Rev ? k - 1 : L - k;  // the read position
        for (unsigned s = 0; s < Ch.Subs; s++)
            if (p == Ch.Sub_At[s]) return Ch.Sub_Base[s];
        return static_cast<unsigned char>(Codes[p]);
    }

    // Fills the chain's levels 1 .. Low - 1 by backward steps from the empty pattern.
    void Rebuild(Range_Chain& Ch, std::string_view Codes) {
        const BWT* B = Ch.Rev ? I.Rev() : I.Fw();
        Range R{0, I.Text_Length()};
        for (unsigned k = 1; k < Ch.Low; k++) {
            const int c = static_cast<int>(Consumed(Ch, Codes, k));
            unsigned At_Start[4], At_End[4];
            Bi_Index::Occ_All(B, R.Start, At_Start);
            Bi_Index::Occ_All(B, R.End + 1, At_End);
            R = {B->cumulativeFreq[c] + At_Start[c] + 1, B->cumulativeFreq[c] + At_End[c]};
            Ch.S[k] = R.Start;
            Ch.E[k] = R.End;
        }
        C.Chain_Rebuilt += Ch.Low - 1;
        Ch.Low = 1;
    }

    // Chain_Decode (src/batman_zero_min.cpp) on a recorded chain: hit rows S[L] + 0 ..
    // Count - 1 resolved by sampled ancestors, the rest queued as walks of their own row.
    void Chain_Backtrack(Range_Chain& Ch, std::string_view Codes, unsigned long Count, unsigned long* Out,
                         const Range* Fw_All) {
        const BWT* F = Ch.Rev ? I.Rev() : I.Fw();
        const unsigned L = static_cast<unsigned>(Codes.size());
        const long Lc = static_cast<long>(L);
        // A hit row's SA value v (in F) as a start in T.
        auto Start = [&](unsigned long v) { return Ch.Rev ? N - v - Lc : v; };
        C.Chain_Hits++;
        C.Chain_Rows += Count;
        Done.assign(Count, 0);
        std::vector<Row_Of>& Rows = Chain_Rows;
        int Class_K[SA_INTERVAL];
        std::size_t Resolved_All = 0;
        auto Backtrack = [&](unsigned long J0, unsigned long J1, unsigned long Scale) {
            Rows.clear();
            for (unsigned long r = J0; r < J1; r++) Rows.push_back({Ch.S[L] + r, r});
            std::size_t Resolved = 0, Mark = Rows.size();
            unsigned long Idle_Scans = 0;
            for (unsigned k = L;;) {
                // The stretch A_a .. A_k of equal width: one level per class of S[k].
                const unsigned long Width = Ch.E[k] - Ch.S[k];
                unsigned a = k;
                for (;;) {
                    if (a == 1) break;
                    if (a - 1 < Ch.Low) Rebuild(Ch, Codes);
                    if (Ch.E[a - 1] - Ch.S[a - 1] != Width) break;
                    a--;
                }
                std::fill(Class_K, Class_K + SA_INTERVAL, -1);
                for (unsigned k2 = k + 1; k2-- > a;)
                    if (Class_K[Ch.S[k2] % SA_INTERVAL] < 0) Class_K[Ch.S[k2] % SA_INTERVAL] = static_cast<int>(k2);
                std::size_t Left = 0;
                for (const Row_Of& R : Rows) {
                    const unsigned long o = R.Row - Ch.S[k];
                    const int k2 = Class_K[(SA_INTERVAL - o % SA_INTERVAL) % SA_INTERVAL];
                    if (k2 >= 0) {
                        Out[R.J] = Start(F->saValue[(Ch.S[k2] + o) / SA_INTERVAL] - (L - static_cast<unsigned>(k2)));
                        Done[R.J] = 1;
                        Resolved++;
                    } else {
                        Rows[Left++] = {Ch.S[a] + o, R.J};
                    }
                }
                Rows.resize(Left);
                if (8 * Left <= 7 * Mark) {
                    Mark = Left;
                    Idle_Scans = 0;
                }
                if (Rows.empty() || a == 1) break;
                if (a - 1 < Ch.Low) Rebuild(Ch, Codes);
                if (Ch.E[a - 1] - Ch.S[a - 1] + 1 > Chain_Scan_Per_Row * Scale * Rows.size() + 4096) break;
                if (Idle_Scans++ >= Chain_Idle_Scans) break;  // ancestors stall (homopolymer)
                // A_{a-1} -> A_a: the t-th row of A_{a-1} with the consumed character is
                // row S[a] + t. A BWT word (16 characters) at a time: a mask of the
                // matching characters, popcounted; only words holding a wanted ancestor
                // are walked. Rows map to BWT positions skipping the $ row.
                const unsigned long Inv = F->inverseSa0;
                const unsigned long Lo = Ch.S[a - 1] - (Ch.S[a - 1] > Inv), Hi = Ch.E[a - 1] - (Ch.E[a - 1] >= Inv);
                const unsigned Pattern = Consumed(Ch, Codes, a) * 0x55555555u;
                unsigned long t = 0, w = Lo / 16;
                std::size_t x = 0;
                for (; w <= Hi / 16 && x < Rows.size(); w++) {
                    const unsigned y = F->bwtCode[w] ^ Pattern;
                    unsigned m = ~(y | y >> 1) & 0x55555555u;  // low bit of each matching pair; pair p at bit 30 - 2p
                    if (w == Lo / 16) m &= (1u << (31 - 2 * (Lo % 16))) - 1;
                    if (w == Hi / 16) m &= ~((1u << (30 - 2 * (Hi % 16))) - 1);
                    const unsigned Count_W = __builtin_popcount(m);
                    if (Ch.S[a] + t + Count_W <= Rows[x].Row) {
                        t += Count_W;
                        continue;
                    }
                    for (; m && x < Rows.size(); t++) {
                        const unsigned Bit = 31 - __builtin_clz(m);
                        m &= ~(1u << Bit);
                        if (Ch.S[a] + t == Rows[x].Row) {
                            const unsigned long Index = w * 16 + (30 - Bit) / 2;
                            Rows[x++].Row = Index + (Index >= Inv);
                        }
                    }
                }
                k = a - 1;
            }
            return Resolved;
        };
        // A sample first (the first Chain_Trial rows, whose ancestors come first in each
        // scan): the rest only if it resolves 7/8 of them.
        const unsigned long Trial = std::min(Count, Chain_Trial);
        const std::size_t Got = Backtrack(0, Trial, Count / std::max(Trial, 1UL));
        const bool Passed = Got >= Trial - Trial / 8;
        Resolved_All = Got;
        if (Passed && Trial < Count) Resolved_All += Backtrack(Trial, Count, 1);
        C.Chain_Resolved += Resolved_All;
        // The rest: after a passed trial, walks of their own row, stepped with the others.
        // A failed trial means ancestors that stall (a homopolymer run, whose rows stay in
        // cache), where the library's serial walk is faster per step: the rest take it, in
        // Fw if Fw_All gives all the rows there (the trial's few are then dropped, as Fw
        // and Rev rows do not pair up).
        if (!Passed && Fw_All && Ch.Rev) {
            C.Serial_Rows += Count;
            for (unsigned long r = 0; r < Count; r++) Out[r] = I.Locate_Fw(Fw_All->Start + r);
            return;
        }
        for (unsigned long r = 0; r < Count; r++) {
            if (Done[r]) continue;
            if (!Passed) {
                C.Serial_Rows++;
                Out[r] = Ch.Rev ? I.Locate_Rev(Ch.S[L] + r, L) : I.Locate_Fw(Ch.S[L] + r);
            } else if (Ch.Rev) {
                Add_Pair(NO_ROW, Ch.S[L] + r, 0, Lc, Out + r);
            } else {
                Add_Pair(Ch.S[L] + r, NO_ROW, 0, Lc, Out + r);
            }
        }
    }

    const Bi_Index& I;
    const unsigned long N;
    std::vector<Walk> Walks;
    Locate_Counts C;
    // Scratch.
    std::vector<unsigned char> Done;
    std::vector<Node> Split_Stack;
    std::vector<Row_Of> Chain_Rows;
};

}  // namespace bifm
