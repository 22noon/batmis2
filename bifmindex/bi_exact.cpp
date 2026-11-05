// bi_exact: exact search of reads, four ways, to compare the bidirectional step with the
// classic one-index searches.
//
//   bi_exact -g <genome_prefix> -q <reads.fa|fq> -mode M [-m N] [-o out] [--timing]
//
//   -mode backward   classic backward search in Fw, the read's last base first
//         forward    classic search in Rev, the read's first base first
//         bi-left    Bi_Index::Search_Left (Extend_Left)
//         bi-right   Bi_Index::Search_Right (Extend_Right)
//   -m N             positions printed per read (default 1, 0 = all)
//   --no-singleton   bi-*: no singleton step; single rows take the narrow step (or the
//                    general step with --threshold 0) (Bi_Index::Singleton)
//   --lookup k       bi-*: start from a table of every k-mer's range (Build_Lookup; its
//                    build is timed as table)
//   --ref            bi-*: finish a single row at its first sampled SA row against the
//                    packed genome (<genome_prefix>.pac; Located_Search)
//   --batch B        bi-*: search B reads at a time, interleaved (Batch_Exact; 0, the
//                    default: one at a time)
//   --early          with --batch: each lane prefetches its next step right after a step,
//                    instead of all lanes at the start of each round
//   --occ lib|fused  bi-*: Occ by the library (BWTAllOccValue) or inline (Occ_All,
//                    Occ_One; the default) (Bi_Index::Fused)
//   --threshold N    bi-*: ranges with End - Start < N take the narrow step (default 128,
//                    0 = never) (Bi_Index::Threshold)
//   --locate lib|bi  bi-*: locate each row by the library's walk (lib), or by the Locator
//                    (bi, the default: range chains and pair walks; bi_locate.hpp)
//   --chain-min N    hits of N rows or more locate by their range chain (default 24)
//   --one-end        pair walks walk the Fw row only (to compare)
//   --split          narrow ranges split into pairs (Children_Left) instead of walking Fw rows
//
// Output, one line per read: id, the number of occurrences, then the positions of N of
// them in T (0-based starts, ascending): the first N rows in Fw (forward: Rev), or with
// --locate bi, of a range chain's index (Rev for bi-right). A read with a base other than ACGT has 0. All reads are
// searched first (timed as search), then located (timed as locate).
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "bi_index.hpp"
#include "bi_locate.hpp"
#include "bi_batch.hpp"
#include "bi_ref.hpp"

using namespace bifm;
using Clock = std::chrono::steady_clock;

namespace {

// The next FASTA or FASTQ record (FASTA sequences may span lines).
bool Next_Record(std::istream& In, std::string& Id, std::string& Seq, std::string& Line, bool& Have_Line) {
    while (Have_Line || std::getline(In, Line)) {
        Have_Line = false;
        if (Line.empty()) continue;
        if (Line[0] == '@') {
            Id = Line.substr(1);
            std::string Plus, Qual;
            return std::getline(In, Seq) && std::getline(In, Plus) && std::getline(In, Qual);
        }
        if (Line[0] == '>') {
            Id = Line.substr(1);
            Seq.clear();
            while (std::getline(In, Line)) {
                if (!Line.empty() && (Line[0] == '>' || Line[0] == '@')) {
                    Have_Line = true;
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
    std::string Prefix, Reads_File, Mode, Out_File;
    unsigned long Max_Positions = 1, Threshold = 128;
    int Lookup = 0;
    unsigned Batch = 0;
    bool Timing = false, Singleton = true, Use_Ref = false, Early = false, Locate_Bi = true;
    for (int i = 1; i < argc; i++) {
        if (!std::strcmp(argv[i], "-g") && i + 1 < argc) Prefix = argv[++i];
        else if (!std::strcmp(argv[i], "-q") && i + 1 < argc) Reads_File = argv[++i];
        else if (!std::strcmp(argv[i], "-mode") && i + 1 < argc) Mode = argv[++i];
        else if (!std::strcmp(argv[i], "-m") && i + 1 < argc) Max_Positions = std::strtoul(argv[++i], nullptr, 10);
        else if (!std::strcmp(argv[i], "-o") && i + 1 < argc) Out_File = argv[++i];
        else if (!std::strcmp(argv[i], "--timing")) Timing = true;
        else if (!std::strcmp(argv[i], "--no-singleton")) Singleton = false;
        else if (!std::strcmp(argv[i], "--ref")) Use_Ref = true;
        else if (!std::strcmp(argv[i], "--early")) Early = true;
        else if (!std::strcmp(argv[i], "--locate") && i + 1 < argc) Locate_Bi = std::strcmp(argv[++i], "lib") != 0;
        else if (!std::strcmp(argv[i], "--chain-min") && i + 1 < argc) Locator::Chain_Min = std::strtoul(argv[++i], nullptr, 10);
        else if (!std::strcmp(argv[i], "--split")) Locator::Split = true;
        else if (!std::strcmp(argv[i], "--one-end")) Locator::Both_Ends = false;
        else if (!std::strcmp(argv[i], "--occ") && i + 1 < argc) Bi_Index::Fused = std::strcmp(argv[++i], "lib") != 0;
        else if (!std::strcmp(argv[i], "--batch") && i + 1 < argc) Batch = static_cast<unsigned>(std::atoi(argv[++i]));
        else if (!std::strcmp(argv[i], "--lookup") && i + 1 < argc) Lookup = std::atoi(argv[++i]);
        else if (!std::strcmp(argv[i], "--threshold") && i + 1 < argc) Threshold = std::strtoul(argv[++i], nullptr, 10);
        else Prefix.clear(), i = argc;
    }
    const int M = Mode == "backward" ? 0 : Mode == "forward" ? 1 : Mode == "bi-left" ? 2 : Mode == "bi-right" ? 3 : -1;
    if (Prefix.empty() || Reads_File.empty() || M < 0 || Lookup < 0 || Lookup > 14) {
        std::fprintf(stderr,
                     "usage: %s -g <genome_prefix> -q <reads.fa|fq> -mode backward|forward|bi-left|bi-right "
                     "[-m N] [-o out] [--timing] [--no-singleton] [--lookup k] [--ref] [--batch B] [--early] [--occ lib|fused] [--threshold N] [--locate lib|bi] [--chain-min N] [--split] [--one-end]\n",
                     argv[0]);
        return 2;
    }

    const auto Load_Start = Clock::now();
    Bi_Index Index(Prefix);
    Index.Singleton = Singleton;
    Index.Threshold = Threshold;
    std::unique_ptr<Packed_Genome> Genome;
    std::unique_ptr<Located_Search> Located;
    if (Use_Ref && M >= 2) {
        Genome = std::make_unique<Packed_Genome>(Prefix);
        Located = std::make_unique<Located_Search>(Index, *Genome);
    }
    const double Load_Seconds = std::chrono::duration<double>(Clock::now() - Load_Start).count();
    const auto Table_Start = Clock::now();
    if (M >= 2) Index.Build_Lookup(Lookup);
    const double Table_Seconds = std::chrono::duration<double>(Clock::now() - Table_Start).count();

    std::ifstream In(Reads_File);
    if (!In) {
        std::fprintf(stderr, "cannot open %s\n", Reads_File.c_str());
        return 2;
    }
    std::vector<std::string> Ids, Codes;
    std::vector<unsigned char> Valid;
    {
        std::string Id, Seq, Line;
        bool Have_Line = false;
        while (Next_Record(In, Id, Seq, Line, Have_Line)) {
            std::string C(Seq.size(), 0);
            bool Ok = !Seq.empty();
            for (std::size_t i = 0; i < Seq.size() && Ok; i++) {
                const int b = Base_Code(Seq[i]);
                Ok = b >= 0;
                C[i] = static_cast<char>(b);
            }
            Ids.push_back(Id.substr(0, Id.find_first_of(" \t")));
            Codes.push_back(std::move(C));
            Valid.push_back(Ok);
        }
    }

    // Search: each read's range (backward: Classic in Fw, forward: Classic in Rev, bi-*:
    // Ranges); with --ref, a read found at a sampled single row has its position instead.
    // With --locate bi, a bi-* read left with Locator::Chain_Min rows or more keeps its
    // range chain (Chains[Chain_Of[r]]).
    const bool Bi = M >= 2, Record = Bi && Locate_Bi;
    std::vector<Bi_Range> Ranges(Codes.size(), Bi_Range{1, 0, 0});
    std::vector<Range> Classic(Bi ? 0 : Codes.size(), Range{1, 0});
    std::vector<unsigned long> Located_Pos(Codes.size(), Exact_Result::NO_POS);
    std::vector<Range_Chain> Chains;
    std::vector<int> Chain_Of(Codes.size(), -1);
    Range_Chain Ch;
    std::vector<uint64_t> Packed;
    Locate_Stats Stats;
    unsigned long Bases = 0;
    const auto Search_Start = Clock::now();
    if (Batch && Bi) {
        for (std::size_t r = 0; r < Codes.size(); r++) Bases += Valid[r] ? Codes[r].size() : 0;
        Batch_Exact(Index, Located.get(), Genome.get(), Codes, Valid, M == 2, Batch, Ranges, Located_Pos, Stats, Early,
                    Record ? &Chains : nullptr, &Chain_Of);
    }
    for (std::size_t r = 0; r < Codes.size() && !(Batch && Bi); r++) {
        if (!Valid[r]) continue;
        const std::string& P = Codes[r];
        Bases += P.size();
        Range_Chain* Rec = Record ? &Ch : nullptr;
        if (Located) {
            Packed.resize((P.size() + 31) / 32);
            Pack_Read(P, Packed.data());
            const Exact_Result X = M == 2 ? Located->Search_Left(P, Packed.data(), Stats, Rec)
                                          : Located->Search_Right(P, Packed.data(), Stats, Rec);
            Ranges[r] = X.Located() ? Bi_Range{0, 0, 0} : X.R;
            Located_Pos[r] = X.Pos;
        } else {
            switch (M) {
                case 0: Classic[r] = Index.Backward_Search(P); break;
                case 1: Classic[r] = Index.Forward_Search(P); break;
                case 2: Ranges[r] = Index.Search_Left(P, Rec); break;
                default: Ranges[r] = Index.Search_Right(P, Rec); break;
            }
        }
        if (Rec && Located_Pos[r] == Exact_Result::NO_POS && Ranges[r].Size() >= Locator::Chain_Min) {
            Chain_Of[r] = static_cast<int>(Chains.size());
            Chains.push_back(Ch);
        }
    }
    const double Search_Seconds = std::chrono::duration<double>(Clock::now() - Search_Start).count();

    // Locate the first N rows of each read (all with -m 0): located reads directly; the
    // classic modes, and bi-* with --locate lib, by the library's walk per row (Fw rows;
    // forward: Rev rows); else by the Locator: chains (their own index's first rows), pair
    // walks and split ranges (all rows), or Fw walks (the first N), all stepped together.
    std::FILE* Out = Out_File.empty() ? stdout : std::fopen(Out_File.c_str(), "w");
    if (!Out) {
        std::fprintf(stderr, "cannot open %s\n", Out_File.c_str());
        return 2;
    }
    auto Size_Of = [&](std::size_t r) { return Bi ? Ranges[r].Size() : Classic[r].Size(); };
    std::vector<unsigned long> Begin(Codes.size() + 1, 0);
    for (std::size_t r = 0; r < Codes.size(); r++)
        Begin[r + 1] = Begin[r] + (Max_Positions ? std::min(Max_Positions, Size_Of(r)) : Size_Of(r));
    std::vector<unsigned long> Pos(Begin.back());
    Locator Loc(Index);
    const auto Locate_Start = Clock::now();
    for (std::size_t r = 0; r < Codes.size(); r++) {
        const unsigned long Want = Begin[r + 1] - Begin[r], L = Codes[r].size();
        unsigned long* O = Pos.data() + Begin[r];
        if (!Want) continue;
        if (Located_Pos[r] != Exact_Result::NO_POS) {
            *O = Located_Pos[r];
        } else if (!Bi) {
            for (unsigned long k = 0; k < Want; k++)
                O[k] = M == 1 ? Index.Locate_Rev(Classic[r].Start + k, L) : Index.Locate_Fw(Classic[r].Start + k);
        } else if (!Locate_Bi) {
            for (unsigned long k = 0; k < Want; k++) O[k] = Index.Locate_Fw(Ranges[r].Fw_Start + k);
        } else if (Chain_Of[r] >= 0) {
            const Range Fw_All = Ranges[r].Fw();
            Loc.Add_Chain(Chains[Chain_Of[r]], Codes[r], Want, O, Want == Ranges[r].Size() ? &Fw_All : nullptr);
        } else if (Want == Ranges[r].Size()) {
            Loc.Add_Range(Ranges[r], 0, static_cast<long>(L), O);
        } else {
            for (unsigned long k = 0; k < Want; k++) Loc.Add_Pair(Ranges[r].Fw_Start + k, Locator::NO_ROW, 0, long(L), O + k);
        }
    }
    Loc.Run();
    const double Locate_Seconds = std::chrono::duration<double>(Clock::now() - Locate_Start).count();
    unsigned long Hits = 0, Mapped = 0;
    for (std::size_t r = 0; r < Codes.size(); r++) {
        std::sort(Pos.begin() + Begin[r], Pos.begin() + Begin[r + 1]);
        Hits += Size_Of(r);
        Mapped += Size_Of(r) > 0;
        std::fprintf(Out, "%s\t%lu", Ids[r].c_str(), Size_Of(r));
        for (unsigned long k = Begin[r]; k < Begin[r + 1]; k++) std::fprintf(Out, "\t%lu", Pos[k]);
        std::fputc('\n', Out);
    }
    if (Out != stdout) std::fclose(Out);

    // The share of bidirectional steps taken from a single row, and from a narrow range
    // (untimed, the first 100k reads, extending right).
    unsigned long Steps = 0, Single = 0, Narrow = 0;
    for (std::size_t r = 0; Timing && r < Codes.size() && r < 100000; r++) {
        if (!Valid[r]) continue;
        Bi_Range R = Index.Empty_Pattern();
        for (std::size_t i = 0; i < Codes[r].size() && !R.Empty(); i++) {
            Steps++;
            Single += R.Size() == 1;
            Narrow += R.Size() > 1 && R.Size() <= Threshold;
            R = Index.Extend_Right(R, Codes[r][i]);
        }
    }
    if (Timing && Bi && Locate_Bi) {
        const Locate_Counts& C = Loc.Counts();
        std::fprintf(stderr,
                     "bi_exact locate: pairs=%lu single_walks=%lu rounds=%lu split_nodes=%lu chain_hits=%lu chain_rows=%lu "
                     "chain_resolved=%lu chain_rebuilt_levels=%lu dollar=%lu steps_per_walk=%.2f serial_rows=%lu\n",
                     C.Pairs, C.Single_Walks, C.Rounds, C.Split_Nodes, C.Chain_Hits, C.Chain_Rows, C.Chain_Resolved,
                     C.Chain_Rebuilt, C.Dollar, C.Pairs + C.Single_Walks ? double(C.Walk_Steps) / (C.Pairs + C.Single_Walks) : 0.0, C.Serial_Rows);
    }
    if (Timing && Located)
        std::fprintf(stderr, "bi_exact: located_by_fixed_row=%lu located_by_moving_row=%lu off_genome=%lu mismatched=%lu\n",
                     Stats.Fixed, Stats.Moving, Stats.Off_Genome, Stats.Mismatched);
    if (Timing)
        std::fprintf(stderr, "bi_exact: singleton_share=%.3f narrow_share=%.3f\n", Steps ? double(Single) / Steps : 0.0,
                     Steps ? double(Narrow) / Steps : 0.0);
    if (Timing)
        std::fprintf(stderr,
                     "bi_exact %s: load=%.3fs table=%.3fs search=%.3fs locate=%.3fs reads=%zu mapped=%lu occurrences=%lu "
                     "ns_per_read=%.1f ns_per_base=%.2f\n",
                     Mode.c_str(), Load_Seconds, Table_Seconds, Search_Seconds, Locate_Seconds, Codes.size(), Mapped, Hits,
                     Codes.empty() ? 0.0 : 1e9 * Search_Seconds / Codes.size(), Bases ? 1e9 * Search_Seconds / Bases : 0.0);
    return 0;
}
