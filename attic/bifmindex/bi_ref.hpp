// bi_ref.hpp: exact search finished against the packed genome (batman_zero_min's located
// single-row paths, benchmark/pac_singleton.md).
//
// Once a pattern has one occurrence, its Bi_Range is a row in each index. A row whose SA
// value is sampled (row % 8 == 0) gives the occurrence's position in T at once:
//   - Fw row r of P = read[i, j): SA_fw[r] is where P starts, so the read starts at
//     SA_fw[r] - i;
//   - Rev row r' of P: P reversed starts at SA_rev[r'] in reversed T, so P starts at
//     n - SA_rev[r'] - (j - i), and the read at n - SA_rev[r'] - j.
// The read then matches if the genome there equals it, compared 32 bases at a time,
// with no more steps. A step moves only its own index's row (the singleton step), so the
// other row is tested once, when the pattern becomes a single row, and the moving row
// after each step.
#pragma once

#include <cstdint>
#include <cstring>
#include <fstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "bi_index.hpp"

namespace bifm {

// The forward genome (<prefix>.pac: 2 bits a base, first base in bits 7-6), with 16 bytes
// of padding so a word can be read past the end.
class Packed_Genome {
public:
    explicit Packed_Genome(const std::string& Prefix) {
        std::ifstream In(Prefix + ".pac", std::ios::binary | std::ios::ate);
        if (!In) throw std::runtime_error("cannot open " + Prefix + ".pac");
        const std::streamsize Size = In.tellg();
        In.seekg(0);
        Pac.assign(static_cast<std::size_t>(Size) + 16, 0);
        if (!In.read(reinterpret_cast<char*>(Pac.data()), Size)) throw std::runtime_error("cannot read " + Prefix + ".pac");
    }

    // The 32 bases from G0, first base in bits 63-62.
    uint64_t Word(unsigned long G0) const {
        const uint8_t* P = Pac.data() + G0 / 4;
        uint64_t G;
        std::memcpy(&G, P, 8);
        G = __builtin_bswap64(G);
        const unsigned Shift = 2 * (G0 % 4);
        if (Shift) G = (G << Shift) | (P[8] >> (8 - Shift));
        return G;
    }

    // Whether the genome from S holds the read of Length bases (Pack_Read).
    bool Matches(unsigned long S, const uint64_t* Packed, std::size_t Length) const {
        for (std::size_t W = 0; W * 32 < Length; W++) {
            uint64_t X = Word(S + W * 32) ^ Packed[W];
            const std::size_t Left = Length - W * 32;
            if (Left < 32) X &= ~(~0ULL >> (2 * Left));  // bases [0, Left) only
            if (X) return false;
        }
        return true;
    }

    // The number of bases where the genome from S differs from the read of Length bases,
    // and in First the first such base (Length if none).
    unsigned Diff(unsigned long S, const uint64_t* Packed, std::size_t Length, std::size_t& First) const {
        unsigned N = 0;
        First = Length;
        for (std::size_t W = 0; W * 32 < Length; W++) {
            uint64_t X = Word(S + W * 32) ^ Packed[W];
            const std::size_t Left = Length - W * 32;
            if (Left < 32) X &= ~(~0ULL >> (2 * Left));
            X = (X | X >> 1) & 0x5555555555555555ULL;  // a base's low bit: set where it differs
            if (X && First == Length) First = W * 32 + __builtin_clzll(X) / 2;
            N += __builtin_popcountll(X);
        }
        return N;
    }

    // The number of bases where the genome from S differs from the read of Length bases;
    // the first Max of them (in order) in At.
    unsigned Diff_At(unsigned long S, const uint64_t* Packed, std::size_t Length, std::size_t* At, unsigned Max) const {
        unsigned N = 0, Kept = 0;
        for (std::size_t W = 0; W * 32 < Length; W++) {
            uint64_t X = Word(S + W * 32) ^ Packed[W];
            const std::size_t Left = Length - W * 32;
            if (Left < 32) X &= ~(~0ULL >> (2 * Left));
            X = (X | X >> 1) & 0x5555555555555555ULL;  // a base's low bit: set where it differs
            N += __builtin_popcountll(X);
            for (; X && Kept < Max; X &= ~(1ULL << (63 - __builtin_clzll(X))))  // first base = highest bit
                At[Kept++] = W * 32 + __builtin_clzll(X) / 2;
        }
        return N;
    }

    // Prefetches the words Matches reads for a read of Length bases at S.
    void Prefetch(unsigned long S, std::size_t Length) const {
        __builtin_prefetch(Pac.data() + S / 4);
        __builtin_prefetch(Pac.data() + (S + Length) / 4);
    }

private:
    std::vector<uint8_t> Pac;
};

// A read's base codes, 32 bases a word, first base in bits 63-62 (as Packed_Genome::Word).
inline void Pack_Read(std::string_view Read, uint64_t* Packed) {
    std::memset(Packed, 0, (Read.size() + 31) / 32 * sizeof(uint64_t));
    for (std::size_t i = 0; i < Read.size(); i++)
        Packed[i / 32] |= static_cast<uint64_t>(static_cast<unsigned char>(Read[i])) << (62 - 2 * (i % 32));
}

// An exact search's result: the read's Bi_Range, or, once located, its one position.
struct Exact_Result {
    static constexpr unsigned long NO_POS = ~0UL;
    Bi_Range R;
    unsigned long Pos = NO_POS;  // located: the read's start in T (R is then one row, unused)
    bool Located() const { return Pos != NO_POS; }
    unsigned long Count() const { return Located() ? 1 : R.Size(); }
};

// Which row located the read (bi_exact --timing).
struct Locate_Stats {
    unsigned long Fixed = 0, Moving = 0, Off_Genome = 0, Mismatched = 0;
};

class Located_Search {
public:
    static constexpr unsigned long SA_INTERVAL = 8;

    Located_Search(const Bi_Index& Index, const Packed_Genome& Genome) : I(Index), G(Genome) {
        if (I.Fw()->saInterval != SA_INTERVAL || I.Rev()->saInterval != SA_INTERVAL)
            throw std::runtime_error("the SA sample interval is not 8");
    }

    // Search_Right (from the read's first base), finished at the first sampled single row.
    // With Ch, each step's Rev range is recorded (Bi_Index::Search_Right).
    Exact_Result Search_Right(std::string_view P, const uint64_t* Packed, Locate_Stats& St, Range_Chain* Ch = nullptr) const {
        const std::size_t L = P.size();
        Bi_Range R = I.Empty_Pattern();
        std::size_t j = 0;  // R is the range of P[0, j)
        if (I.Lookup_Depth() && L >= static_cast<std::size_t>(I.Lookup_Depth())) {
            R = I.Look_Up(P.data());
            j = I.Lookup_Depth();
        }
        if (Ch) Ch->Reset(L, true, j ? static_cast<unsigned>(j) : 1);
        if (Ch && j && !R.Empty()) Ch->Note(j, R);
        bool Single = false;
        for (;;) {
            if (R.Empty()) return {R};
            if (R.Size() == 1) {
                if (!Single) {  // the fixed Fw row, once
                    Single = true;
                    if (R.Fw_Start % SA_INTERVAL == 0) return St.Fixed++, Verify(Sa(I.Fw(), R.Fw_Start), P, Packed, St);
                }
                if (R.Rev_Start % SA_INTERVAL == 0)
                    return St.Moving++, Verify(long(I.Text_Length()) - Sa(I.Rev(), R.Rev_Start) - long(j), P, Packed, St);
            }
            if (j == L) return {R};
            R = I.Extend_Right(R, P[j++]);
            if (Ch && !R.Empty()) Ch->Note(j, R);
        }
    }

    // Search_Left (from the read's last base): the mirror, Rev's row fixed and Fw's moving;
    // with Ch, Fw ranges recorded.
    Exact_Result Search_Left(std::string_view P, const uint64_t* Packed, Locate_Stats& St, Range_Chain* Ch = nullptr) const {
        const std::size_t L = P.size();
        Bi_Range R = I.Empty_Pattern();
        std::size_t i = L;  // R is the range of P[i, L)
        if (I.Lookup_Depth() && L >= static_cast<std::size_t>(I.Lookup_Depth())) {
            i = L - I.Lookup_Depth();
            R = I.Look_Up(P.data() + i);
        }
        if (Ch) Ch->Reset(L, false, i == L ? 1 : static_cast<unsigned>(L - i));
        if (Ch && i < L && !R.Empty()) Ch->Note(L - i, R);
        bool Single = false;
        for (;;) {
            if (R.Empty()) return {R};
            if (R.Size() == 1) {
                if (!Single) {  // the fixed Rev row, once
                    Single = true;
                    if (R.Rev_Start % SA_INTERVAL == 0)
                        return St.Fixed++, Verify(long(I.Text_Length()) - Sa(I.Rev(), R.Rev_Start) - long(L), P, Packed, St);
                }
                if (R.Fw_Start % SA_INTERVAL == 0) return St.Moving++, Verify(Sa(I.Fw(), R.Fw_Start) - long(i), P, Packed, St);
            }
            if (i == 0) return {R};
            R = I.Extend_Left(R, P[--i]);
            if (Ch && !R.Empty()) Ch->Note(L - i, R);
        }
    }

    // The SA value of sampled Row in B.
    static long Sa(const BWT* B, unsigned long Row) { return static_cast<long>(B->saValue[Row / SA_INTERVAL]); }

    // The read at S, the only place its matched part occurs: found if it fits in T and the
    // genome there holds it.
    Exact_Result Verify(long S, std::string_view P, const uint64_t* Packed, Locate_Stats& St) const {
        if (S < 0 || static_cast<unsigned long>(S) + P.size() > I.Text_Length()) {
            St.Off_Genome++;
            return {Bi_Range{1, 0, 0}};
        }
        if (!G.Matches(static_cast<unsigned long>(S), Packed, P.size())) {
            St.Mismatched++;
            return {Bi_Range{1, 0, 0}};
        }
        Exact_Result X;
        X.R = {0, 0, 0};
        X.Pos = static_cast<unsigned long>(S);
        return X;
    }

private:
    const Bi_Index& I;
    const Packed_Genome& G;
};

}  // namespace bifm
