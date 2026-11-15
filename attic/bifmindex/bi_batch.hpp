// bi_batch.hpp: exact search of B reads at a time (batman_zero_min's Batch_Search for -n 0).
//
// Each lane holds one read's search, cut into steps that each read one thing from
// memory: the lookup entry, a step's occurrence rows, a sampled SA value, the genome.
// A round first prefetches what every lane's next step reads, then takes each lane's
// step, so the lanes' cache misses overlap. With Early, each lane instead prefetches its
// next step right after taking one, so every lane's prefetch has a whole round to arrive
// (with the round's prefetches first, the round's first lane has almost none). A finished
// lane takes the next read. The
// results are those of Bi_Index::Search_Left / Search_Right, or with a Located_Search,
// those of its searches.
#pragma once

#include <string>
#include <string_view>
#include <vector>

#include "bi_index.hpp"
#include "bi_locate.hpp"
#include "bi_ref.hpp"

namespace bifm {

// Reads Codes (base codes; those with !Valid are skipped, with an empty range) searched
// by Extend_Left from their last base (Left) or Extend_Right from their first, B at a
// time. Ranges[r] gets read r's range, or with Loc, Pos[r] its position once located
// (Exact_Result::NO_POS otherwise). With Chains, a read left with Locator::Chain_Min rows
// or more gets its recorded range chain (Fw going left, Rev going right) as
// (*Chains)[(*Chain_Of)[r]]; -1 otherwise.
inline void Batch_Exact(const Bi_Index& I, const Located_Search* Loc, const Packed_Genome* G,
                        const std::vector<std::string>& Codes, const std::vector<unsigned char>& Valid, bool Left,
                        unsigned B, std::vector<Bi_Range>& Ranges, std::vector<unsigned long>& Pos, Locate_Stats& St,
                        bool Early = false, std::vector<Range_Chain>* Chains = nullptr, std::vector<int>* Chain_Of = nullptr) {
    enum Phase { LOOKUP, STEP, SA_FIXED, SA_MOVING, VERIFY };
    struct Lane {
        std::size_t Read;
        std::string_view P;
        Bi_Range R;
        std::size_t At;  // R is the range of P[At, |P|) (Left) or P[0, At)
        bool Single;     // the fixed row has been tested
        Phase Ph;
        long S;          // VERIFY: the read's start
        std::vector<uint64_t> Packed;
        Range_Chain Ch;  // with Chains
    };
    const std::size_t K = static_cast<std::size_t>(I.Lookup_Depth());
    const unsigned long N = I.Text_Length();
    Ranges.assign(Codes.size(), Bi_Range{1, 0, 0});
    Pos.assign(Codes.size(), Exact_Result::NO_POS);
    if (Chains) Chains->clear(), Chain_Of->assign(Codes.size(), -1);

    std::size_t Next = 0;
    auto Finish = [&](Lane& L, const Exact_Result& X) {
        Ranges[L.Read] = X.Located() ? Bi_Range{0, 0, 0} : X.R;
        Pos[L.Read] = X.Pos;
        if (Chains && !X.Located() && X.R.Size() >= Locator::Chain_Min) {
            (*Chain_Of)[L.Read] = static_cast<int>(Chains->size());
            Chains->push_back(L.Ch);
        }
    };
    // After R changes: done, a sampled single row to locate, or the next step.
    // Returns true when the read is done.
    auto Check = [&](Lane& L) {
        if (L.R.Empty()) return Finish(L, {L.R}), true;
        if (Loc && L.R.Size() == 1) {
            if (!L.Single) {
                L.Single = true;
                if ((Left ? L.R.Rev_Start : L.R.Fw_Start) % Located_Search::SA_INTERVAL == 0) return L.Ph = SA_FIXED, false;
            }
            if ((Left ? L.R.Fw_Start : L.R.Rev_Start) % Located_Search::SA_INTERVAL == 0) return L.Ph = SA_MOVING, false;
        }
        if (Left ? L.At == 0 : L.At == L.P.size()) return Finish(L, {L.R}), true;
        L.Ph = STEP;
        return false;
    };
    // Starts the next valid read in L; false if none is left.
    auto Fill = [&](Lane& L) {
        while (Next < Codes.size()) {
            if (!Valid[Next]) {
                Next++;
                continue;
            }
            L.Read = Next++;
            L.P = Codes[L.Read];
            L.Single = false;
            if (Loc) {
                L.Packed.resize((L.P.size() + 31) / 32);
                Pack_Read(L.P, L.Packed.data());
            }
            if (K && L.P.size() >= K) {
                L.Ph = LOOKUP;
                return true;
            }
            L.R = I.Empty_Pattern();
            L.At = Left ? L.P.size() : 0;
            if (Chains) L.Ch.Reset(L.P.size(), !Left, 1);
            if (!Check(L)) return true;
        }
        return false;
    };
    auto Prefetch = [&](const Lane& L) {
        switch (L.Ph) {
            case LOOKUP: I.Prefetch_Look_Up(L.P.data() + (Left ? L.P.size() - K : 0)); break;
            case STEP: I.Prefetch_Step(L.R, Left); break;
            case SA_FIXED:
                if (Left) __builtin_prefetch(I.Rev()->saValue + L.R.Rev_Start / Located_Search::SA_INTERVAL);
                else __builtin_prefetch(I.Fw()->saValue + L.R.Fw_Start / Located_Search::SA_INTERVAL);
                break;
            case SA_MOVING:
                if (Left) __builtin_prefetch(I.Fw()->saValue + L.R.Fw_Start / Located_Search::SA_INTERVAL);
                else __builtin_prefetch(I.Rev()->saValue + L.R.Rev_Start / Located_Search::SA_INTERVAL);
                break;
            case VERIFY:
                if (L.S >= 0 && static_cast<unsigned long>(L.S) + L.P.size() <= N) G->Prefetch(L.S, L.P.size());
                break;
        }
    };
    // One step of L; true when the read is done.
    auto Step = [&](Lane& L) {
        const long Len = static_cast<long>(L.P.size());
        switch (L.Ph) {
            case LOOKUP:
                L.At = Left ? L.P.size() - K : K;
                L.R = I.Look_Up(L.P.data() + (Left ? L.At : 0));
                if (Chains) {
                    L.Ch.Reset(L.P.size(), !Left, static_cast<unsigned>(K));
                    if (!L.R.Empty()) L.Ch.Note(K, L.R);
                }
                return Check(L);
            case STEP:
                if (Left) L.R = I.Extend_Left(L.R, L.P[--L.At]);
                else L.R = I.Extend_Right(L.R, L.P[L.At++]);
                if (Chains && !L.R.Empty()) L.Ch.Note(Left ? L.P.size() - L.At : L.At, L.R);
                return Check(L);
            case SA_FIXED:  // Left: Rev's row of the read's suffix P[At, |P|); else Fw's of P[0, At)
                St.Fixed++;
                L.S = Left ? long(N) - Located_Search::Sa(I.Rev(), L.R.Rev_Start) - Len
                           : Located_Search::Sa(I.Fw(), L.R.Fw_Start);
                L.Ph = VERIFY;
                return false;
            case SA_MOVING:
                St.Moving++;
                L.S = Left ? Located_Search::Sa(I.Fw(), L.R.Fw_Start) - long(L.At)
                           : long(N) - Located_Search::Sa(I.Rev(), L.R.Rev_Start) - long(L.At);
                L.Ph = VERIFY;
                return false;
            case VERIFY:
                Finish(L, Loc->Verify(L.S, L.P, L.Packed.data(), St));
                return true;
        }
        return true;
    };

    std::vector<Lane> Lanes(B ? B : 1);
    unsigned Active = 0;
    while (Active < Lanes.size() && Fill(Lanes[Active])) Active++;
    if (Early)
        for (unsigned i = 0; i < Active; i++) Prefetch(Lanes[i]);
    while (Active) {
        if (!Early)
            for (unsigned i = 0; i < Active; i++) Prefetch(Lanes[i]);
        for (unsigned i = 0; i < Active;) {
            Lane& L = Lanes[i];
            if (!Step(L) || Fill(L)) {
                if (Early) Prefetch(L);
                i++;
                continue;
            }
            std::swap(L, Lanes[--Active]);  // no reads left: close the gap
        }
    }
}

}  // namespace bifm
