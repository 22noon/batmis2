// bi_mm: one- and two-mismatch search on the bidirectional index (One_Mismatch,
// Two_Mismatch), to compare with batman_zero_min -n 1 / -n 2 and brute force.
//
//   bi_mm -g <genome_prefix> -q <reads.fa|fq> [-n 1|2] [-m N] [-o out] [--lookup k] [--ref]
//         [--batch B] [--counts] [--locate lib|bi] [--chain-min N] [--split] [--one-end] [--timing]
//
//   -n 1|2      the most mismatches (default 1); -n 2 (Two_Mismatch, Batch_Two_Mismatch)
//               prints the two-mismatch count after the one-mismatch one
//   -m N        stop a read's search once it has N hits (between stages and kept
//               children, as One_Mismatch::Search; default 0, all)
//   --lookup k  start the half searches from a table of every k-mer's range (0-14)
//   --ref       finish single rows at their first sampled SA row against the packed
//               genome (<genome_prefix>.pac)
//   --batch B   search B reads at a time, interleaved (Batch_One_Mismatch; 0, the
//               default: one at a time)
//   --occ lib|fused  Occ by the library (BWTAllOccValue) or inline (the default)
//   --counts    print the counts only, not the positions
//   --locate    lib: each unlocated row by the library's walk (BWTSaValue); bi (the
//               default): Locator's pair walks and range chains
//   --chain-min N  hits of N rows or more use their range chain (default 24)
//   --one-end        pair walks walk the Fw row only (to compare)
//   --split     narrow ranges split into pairs (Children_Left) instead of walking Fw rows
//
// Output, one line per read: id, the number of exact occurrences, the number with one
// mismatch (and with -n 2, two), then every position (0-based start in T) as pos:mismatches, ascending. A
// read with a base other than ACGT has 0 and 0. All reads are searched first (timed as
// search), then located (timed as locate).
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "bi_index.hpp"
#include "bi_locate.hpp"
#include "bi_mismatch.hpp"
#include "bi_two_mismatch.hpp"
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
    std::string Prefix, Reads_File, Out_File;
    unsigned long Max_Hits = 0;
    int Lookup = 0;
    unsigned Batch = 0;
    int Max_Mismatches = 1;
    unsigned long Threshold = 128;
    bool Timing = false, Use_Ref = false, Counts_Only = false, Locate_Bi = true;
    for (int i = 1; i < argc; i++) {
        if (!std::strcmp(argv[i], "-g") && i + 1 < argc) Prefix = argv[++i];
        else if (!std::strcmp(argv[i], "-q") && i + 1 < argc) Reads_File = argv[++i];
        else if (!std::strcmp(argv[i], "--schemes") && i + 1 < argc) Pass2_Core::Scheme_Mask = std::strtoul(argv[++i], nullptr, 0);
        else if (!std::strcmp(argv[i], "--no-settle")) Pass2_Core::Settle = false;
        else if (!std::strcmp(argv[i], "--no-bound")) Pass2_Core::Bound = false;
        else if (!std::strcmp(argv[i], "--occ") && i + 1 < argc) Bi_Index::Fused = std::strcmp(argv[++i], "lib") != 0;
        else if (!std::strcmp(argv[i], "--threshold") && i + 1 < argc) Threshold = std::strtoul(argv[++i], nullptr, 10);
        else if (!std::strcmp(argv[i], "-n") && i + 1 < argc) Max_Mismatches = std::atoi(argv[++i]);
        else if (!std::strcmp(argv[i], "-m") && i + 1 < argc) Max_Hits = std::strtoul(argv[++i], nullptr, 10);
        else if (!std::strcmp(argv[i], "-o") && i + 1 < argc) Out_File = argv[++i];
        else if (!std::strcmp(argv[i], "--timing")) Timing = true;
        else if (!std::strcmp(argv[i], "--ref")) Use_Ref = true;
        else if (!std::strcmp(argv[i], "--counts")) Counts_Only = true;
        else if (!std::strcmp(argv[i], "--locate") && i + 1 < argc) Locate_Bi = std::strcmp(argv[++i], "lib") != 0;
        else if (!std::strcmp(argv[i], "--chain-min") && i + 1 < argc) Locator::Chain_Min = std::strtoul(argv[++i], nullptr, 10);
        else if (!std::strcmp(argv[i], "--split")) Locator::Split = true;
        else if (!std::strcmp(argv[i], "--one-end")) Locator::Both_Ends = false;
        else if (!std::strcmp(argv[i], "--chain-idle") && i + 1 < argc) Locator::Chain_Idle_Scans = std::strtoul(argv[++i], nullptr, 10);
        else if (!std::strcmp(argv[i], "--window") && i + 1 < argc) Locator::Window = std::strtoul(argv[++i], nullptr, 10);
        else if (!std::strcmp(argv[i], "--batch") && i + 1 < argc) Batch = static_cast<unsigned>(std::atoi(argv[++i]));
        else if (!std::strcmp(argv[i], "--lookup") && i + 1 < argc) Lookup = std::atoi(argv[++i]);
        else Prefix.clear(), i = argc;
    }
    if (Prefix.empty() || Reads_File.empty() || Lookup < 0 || Lookup > 14 || Max_Mismatches < 1 || Max_Mismatches > 2) {
        std::fprintf(stderr,
                     "usage: %s -g <genome_prefix> -q <reads.fa|fq> [-n 1|2] [-m N] [-o out] [--lookup k] [--ref] [--batch B] [--counts] [--locate lib|bi] [--chain-min N] [--split] [--one-end] [--timing]\n",
                     argv[0]);
        return 2;
    }

    const auto Load_Start = Clock::now();
    Bi_Index Index(Prefix);
    Index.Threshold = Threshold;
    std::unique_ptr<Packed_Genome> Genome;
    if (Use_Ref) Genome = std::make_unique<Packed_Genome>(Prefix);
    One_Mismatch Search(Index, Genome.get(), Locate_Bi && !Counts_Only);
    Two_Mismatch Search2(Index, Genome.get(), Locate_Bi && !Counts_Only);

    const double Load_Seconds = std::chrono::duration<double>(Clock::now() - Load_Start).count();
    const auto Table_Start = Clock::now();
    Index.Build_Lookup(Lookup);
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

    std::vector<Mm_Result> Results(Codes.size());
    std::vector<uint64_t> Packed;
    Locate_Stats Stats;
    const auto Search_Start = Clock::now();
    unsigned long Pass2_Counts[3] = {0, 0, 0};  // -n 2 batched: nodes, located nodes, settled
    if (Batch && Max_Mismatches == 2)
        Batch_Two_Mismatch(Index, Genome.get(), Codes, Valid, Batch, Max_Hits, Results, Stats, Locate_Bi && !Counts_Only,
                           Pass2_Counts);
    else if (Batch)
        Batch_One_Mismatch(Index, Genome.get(), Codes, Valid, Batch, Max_Hits, Results, Stats, Locate_Bi && !Counts_Only);
    for (std::size_t r = 0; r < Codes.size() && !Batch; r++) {
        if (!Valid[r]) {
            Results[r].Clear();
            continue;
        }
        const std::string& P = Codes[r];
        if (Genome) {
            Packed.resize((P.size() + 31) / 32);
            Pack_Read(P, Packed.data());
        }
        if (Max_Mismatches == 2) Search2.Search(P, Packed.data(), Max_Hits, Results[r], Stats);
        else Search.Search(P, Packed.data(), Max_Hits, Results[r], Stats);
    }
    const double Search_Seconds = std::chrono::duration<double>(Clock::now() - Search_Start).count();

    std::FILE* Out = Out_File.empty() ? stdout : std::fopen(Out_File.c_str(), "w");
    if (!Out) {
        std::fprintf(stderr, "cannot open %s\n", Out_File.c_str());
        return 2;
    }
    // Locate: every hit's positions into one array (read r's from Begin[r]); located hits
    // directly, the rest by the library's walk per Fw row (lib), or by the Locator: pair
    // walks, split ranges, and chains, all stepped together.
    unsigned long Exact_Hits = 0, Mm_Hits = 0, Two_Hits = 0, Mapped = 0;
    std::vector<unsigned long> Begin(Codes.size() + 1, 0);
    for (std::size_t r = 0; r < Codes.size(); r++)
        Begin[r + 1] = Begin[r] + (Counts_Only ? 0 : Results[r].Exact.Count() + Results[r].Hit_Count);
    std::vector<unsigned long> Pos(Begin.back());
    std::vector<unsigned char> Mism(Begin.back());
    Locator Loc(Index);
    const auto Locate_Start = Clock::now();
    for (std::size_t r = 0; r < Codes.size() && !Counts_Only; r++) {
        Mm_Result& X = Results[r];
        unsigned long At = Begin[r];
        auto Add = [&](const Bi_Range& R, unsigned long Located, int Chain, int Mismatches) {
            const unsigned long Count = Located != Exact_Result::NO_POS ? 1 : R.Size();
            std::fill(Mism.begin() + At, Mism.begin() + At + Count, Mismatches);
            if (Located != Exact_Result::NO_POS) Pos[At] = Located;
            else if (!Locate_Bi) for (unsigned long k = 0; k < Count; k++) Pos[At + k] = Index.Locate_Fw(R.Fw_Start + k);
            else if (Chain >= 0) {
                const Range Fw_All = R.Fw();
                Loc.Add_Chain(X.Chains[Chain], Codes[r], Count, Pos.data() + At, &Fw_All);
            }
            else Loc.Add_Range(R, 0, static_cast<long>(Codes[r].size()), Pos.data() + At);
            At += Count;
        };
        if (X.Exact.Count()) Add(X.Exact.R, X.Exact.Pos, X.Exact_Chain, 0);
        for (const Mm_Hit& H : X.Hits) Add(H.R, H.Pos, H.Chain, H.Mismatches);
    }
    Loc.Run();
    const double Locate_Seconds = std::chrono::duration<double>(Clock::now() - Locate_Start).count();
    std::vector<std::pair<unsigned long, int>> Positions;
    for (std::size_t r = 0; r < Codes.size(); r++) {
        const Mm_Result& X = Results[r];
        Positions.clear();
        for (unsigned long k = Begin[r]; k < Begin[r + 1]; k++) Positions.push_back({Pos[k], Mism[k]});
        std::sort(Positions.begin(), Positions.end());
        Exact_Hits += X.Exact.Count();
        Mm_Hits += X.Hit_Count;
        Mapped += X.Exact.Count() + X.Hit_Count > 0;
        unsigned long Two = 0;
        for (const Mm_Hit& H : X.Hits) Two += H.Mismatches == 2 ? H.Count() : 0;
        Two_Hits += Two;
        std::fprintf(Out, "%s\t%lu\t%lu", Ids[r].c_str(), X.Exact.Count(), X.Hit_Count - Two);
        if (Max_Mismatches == 2) std::fprintf(Out, "\t%lu", Two);
        for (const auto& [p, k] : Positions) std::fprintf(Out, "\t%lu:%d", p, k);
        std::fputc('\n', Out);
    }
    if (Out != stdout) std::fclose(Out);

    if (Timing && Locate_Bi) {
        const Locate_Counts& C = Loc.Counts();
        std::fprintf(stderr,
                     "bi_mm locate: pairs=%lu single_walks=%lu rounds=%lu split_nodes=%lu chain_hits=%lu chain_rows=%lu "
                     "chain_resolved=%lu chain_rebuilt_levels=%lu dollar=%lu steps_per_walk=%.2f serial_rows=%lu\n",
                     C.Pairs, C.Single_Walks, C.Rounds, C.Split_Nodes, C.Chain_Hits, C.Chain_Rows, C.Chain_Resolved,
                     C.Chain_Rebuilt, C.Dollar, C.Pairs + C.Single_Walks ? double(C.Walk_Steps) / (C.Pairs + C.Single_Walks) : 0.0, C.Serial_Rows);
    }
    if (Timing && Max_Mismatches == 2)
        std::fprintf(stderr, "bi_mm: two_mismatch=%lu pass2_nodes=%lu pass2_located_nodes=%lu settled=%lu\n", Two_Hits,
                     Search2.Nodes + Pass2_Counts[0], Search2.Located_Nodes + Pass2_Counts[1], Search2.Settled + Pass2_Counts[2]);
    if (Timing && Genome)
        std::fprintf(stderr, "bi_mm: located=%lu off_genome=%lu mismatched=%lu\n", Stats.Moving, Stats.Off_Genome,
                     Stats.Mismatched);
    if (Timing)
        std::fprintf(stderr,
                     "bi_mm: load=%.3fs table=%.3fs search=%.3fs locate=%.3fs reads=%zu mapped=%lu exact=%lu one_mismatch=%lu "
                     "ns_per_read=%.1f\n",
                     Load_Seconds, Table_Seconds, Search_Seconds, Locate_Seconds, Codes.size(), Mapped, Exact_Hits, Mm_Hits,
                     Codes.empty() ? 0.0 : 1e9 * Search_Seconds / Codes.size());
    return 0;
}
