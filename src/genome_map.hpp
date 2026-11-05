// Genome_Map: the sequences of a BatMis index and their N runs, to place a hit of the
// concatenated text (forward-genome coordinate) on a sequence.
//
// The index (index/bwt_build_index.cpp) concatenates the FASTA sequences with no
// separator and stores each non-ACGT base as 'c'. Its side files:
//   <genome>.ann.location  "0", then a name line and a length line per sequence;
//   <genome>.N             a FASTA header line per sequence, then "start:end" (end
//                          exclusive, in-sequence) per run of non-ACGT bases.
// A reference shorter than the index minimum gets a padding sequence (Pad_Name) last.
//
// Place rejects a hit that leaves its sequence, lies on the padding sequence or covers
// an N (the index matched it against the 'c' stored there). The SAM header lists the
// sequences that can hold a hit: not the padding sequence, nor empty ones (SAM needs
// LN >= 1), so a hit's Tid is its sequence's place in that list.
//   Sequence lookup: the text is cut into blocks of 2^Shift bases, Shift chosen so a
//     block holds about one sequence start; First[b] is the sequence holding the
//     block's first base. A hit's sequence lies in [First[b], First[b + 1]]: usually
//     one candidate, else a binary search over that block's sequences only.
//   N runs: a bitmap with one bit per 4096-base block (any N in it) answers most hits
//     in a word or two; only a hit touching a marked block searches the run list
//     (sorted, global, non-overlapping).
#pragma once

#include <algorithm>
#include <cstdint>
#include <fstream>
#include <string>
#include <vector>

class Genome_Map {
public:
    static constexpr const char* Pad_Name = "__batmis_pad__";

    // Reads <Prefix>.ann.location and <Prefix>.N. Returns false with Error set if they
    // are missing or malformed.
    bool Load(const std::string& Prefix, std::string& Error) {
        Names.clear();
        Starts.clear();
        Lengths.clear();
        std::ifstream Loc(Prefix + ".ann.location");
        if (!Loc) return Fail(Error, "cannot open " + Prefix + ".ann.location");
        std::string Line, Name;
        if (!std::getline(Loc, Line) || Line != "0") return Fail(Error, Prefix + ".ann.location: no leading 0");
        uint64_t At = 0;
        while (std::getline(Loc, Name)) {
            if (Name.empty()) continue;
            if (!std::getline(Loc, Line)) return Fail(Error, Prefix + ".ann.location: no length for " + Name);
            const uint64_t Length = std::stoull(Line);
            Names.push_back(Name);
            Starts.push_back(At);
            Lengths.push_back(Length);
            At += Length;
        }
        if (Names.empty()) return Fail(Error, Prefix + ".ann.location: no sequences");
        Total = At;
        Pad = Names.back() == Pad_Name ? static_cast<int>(Names.size()) - 1 : -1;
        Tids.assign(Names.size(), -1);
        Header.clear();
        for (int I = 0; I < static_cast<int>(Names.size()); I++)
            if (I != Pad && Lengths[I] > 0) {
                Tids[I] = static_cast<int>(Header.size());
                Header.push_back(I);
            }
        Build_Blocks();
        return Load_N(Prefix + ".N", Error);
    }

    // Places [G, G + L) (G: forward-text coordinate): Tid is its sequence's place in the
    // SAM header and Pos its 0-based start there. False if it leaves the sequence, lies on
    // padding or covers an N.
    bool Place(uint64_t G, int L, int& Tid, int64_t& Pos) const {
        const uint64_t E = G + static_cast<uint64_t>(L);
        if (L <= 0 || E > Total) return false;
        const int I = Find(G);
        if (I == Pad || E > Starts[I] + Lengths[I] || Has_N(G, E)) return false;
        Tid = Tids[I];
        Pos = static_cast<int64_t>(G - Starts[I]);
        return true;
    }

    // The SAM header's sequences (Tid 0 .. Count() - 1).
    int Count() const { return static_cast<int>(Header.size()); }
    const std::string& Name(int Tid) const { return Names[Header[Tid]]; }
    uint64_t Length(int Tid) const { return Lengths[Header[Tid]]; }
    // The header place of the index's I-th sequence, or -1 (padding, empty).
    int Tid_Of(int I) const { return Tids[I]; }
    uint64_t Text_Length() const { return Total; }

private:
    static constexpr unsigned N_SHIFT = 12;  // N bitmap block: 4096 bases

    std::vector<std::string> Names;
    std::vector<uint64_t> Starts, Lengths;
    uint64_t Total = 0;
    int Pad = -1;
    std::vector<int> Tids, Header;  // sequence -> header place, and back
    unsigned Shift = 0;
    std::vector<uint32_t> First;              // per block, plus a sentinel
    std::vector<uint64_t> N_Bits;             // per 4096-base block: any N
    std::vector<uint64_t> Run_Starts, Run_Ends;

    static bool Fail(std::string& Error, const std::string& What) {
        Error = What;
        return false;
    }

    void Build_Blocks() {
        const uint64_t Per = std::max<uint64_t>(1, Total / Names.size());
        Shift = 0;
        while (Shift < 63 && (uint64_t{2} << Shift) <= Per) Shift++;
        Shift = std::min(24u, std::max(8u, Shift));
        const uint64_t Blocks = (Total >> Shift) + 1;
        First.assign(Blocks + 1, 0);
        uint32_t I = 0;
        for (uint64_t B = 0; B < Blocks; B++) {
            const uint64_t G = B << Shift;
            while (I + 1 < Names.size() && Starts[I + 1] <= G) I++;
            First[B] = I;
        }
        First[Blocks] = static_cast<uint32_t>(Names.size() - 1);
    }

    // The sequence holding G (G < Total): the last one starting at or before G, which
    // skips empty sequences.
    int Find(uint64_t G) const {
        const uint64_t B = G >> Shift;
        const uint32_t Lo = First[B], Hi = First[B + 1];
        if (Lo == Hi || Starts[Lo + 1] > G) return static_cast<int>(Lo);
        return static_cast<int>(std::upper_bound(Starts.begin() + Lo, Starts.begin() + Hi + 1, G) - Starts.begin()) - 1;
    }

    bool Load_N(const std::string& Name, std::string& Error) {
        Run_Starts.clear();
        Run_Ends.clear();
        N_Bits.assign(((Total >> N_SHIFT) + 64) / 64, 0);
        std::ifstream In(Name);
        if (!In) return Fail(Error, "cannot open " + Name);
        std::string Line;
        int Seq = -1;
        while (std::getline(In, Line)) {
            if (Line.empty()) continue;
            if (Line[0] == '>') {
                if (++Seq >= static_cast<int>(Names.size())) return Fail(Error, Name + ": more sequences than .ann.location");
                continue;
            }
            const std::size_t Colon = Line.find(':');
            if (Seq < 0 || Colon == std::string::npos) return Fail(Error, Name + ": bad line " + Line);
            const uint64_t S = Starts[Seq] + std::stoull(Line.substr(0, Colon));
            const uint64_t E = Starts[Seq] + std::stoull(Line.substr(Colon + 1));
            if (E <= S) continue;
            Run_Starts.push_back(S);
            Run_Ends.push_back(E);
            for (uint64_t B = S >> N_SHIFT; B <= (E - 1) >> N_SHIFT; B++) N_Bits[B / 64] |= uint64_t{1} << (B % 64);
        }
        return true;
    }

    // Whether [G, E) holds an N.
    bool Has_N(uint64_t G, uint64_t E) const {
        bool Marked = false;
        for (uint64_t B = G >> N_SHIFT; B <= (E - 1) >> N_SHIFT && !Marked; B++) Marked = N_Bits[B / 64] >> (B % 64) & 1;
        if (!Marked) return false;
        const std::size_t I = std::upper_bound(Run_Starts.begin(), Run_Starts.end(), E - 1) - Run_Starts.begin();
        return I > 0 && Run_Ends[I - 1] > G;
    }
};
