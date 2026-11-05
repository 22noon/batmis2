// Sam_Writer: SAM or BAM output through htslib, for ungapped hits of single reads.
//
// Header: @HD (unsorted), an @SQ per sequence of the Genome_Map (not the padding
// sequence) and @PG with the command line. Each mapped record has CIGAR <L>M, MAPQ 255
// (unavailable) and the tags NM:i (mismatches), MD:Z (reference bases at the mismatches)
// and NH:i (records written for the read); the first record of a read is primary and the
// rest are secondary (0x100). A read with no record is written unmapped (FLAG 4).
//
// Records are built by a Sam_Encoder, on any thread, into bytes: SAM lines (sam_format1),
// or BAM records laid out as bam_write1 writes them. Sam_Writer::Write then writes the
// bytes in order, as sam_write1 would have, so the file does not depend on how the
// records were split between encoders.
#pragma once

#include <htslib/bgzf.h>
#include <htslib/hfile.h>
#include <htslib/kstring.h>
#include <htslib/sam.h>

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <string>
#include <utility>
#include <vector>

#include "genome_map.hpp"

static_assert(__BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__, "Sam_Encoder lays out BAM records for a little-endian host");

class Sam_Writer;

class Sam_Encoder {
public:
    Sam_Encoder() = default;
    Sam_Encoder(const Sam_Encoder&) = delete;
    Sam_Encoder& operator=(const Sam_Encoder&) = delete;
    ~Sam_Encoder() {
        if (Record) bam_destroy1(Record);
        free(Line.s);
    }

    // A mapped record: Seq (ACGT) and Qual (FASTQ, '!'-based; empty for FASTA) already on
    // the reference strand; Mismatches holds (read position, reference base 0-3) pairs.
    bool Mapped(const std::string& Name, uint16_t Flag, int Tid, int64_t Pos, const std::string& Seq,
                const std::string& Qual, std::vector<std::pair<int, int>>& Mismatches, int Records) {
        const uint32_t Cigar = bam_cigar_gen(Seq.size(), BAM_CMATCH);
        if (!Set(Name, Flag, Tid, Pos, 255, 1, &Cigar, Seq, Qual)) return false;
        const int32_t NM = static_cast<int32_t>(Mismatches.size());
        const int32_t NH = Records;
        const std::string MD = Md_String(static_cast<int>(Seq.size()), Mismatches);
        return bam_aux_append(Record, "NM", 'i', sizeof NM, reinterpret_cast<const uint8_t*>(&NM)) == 0 &&
               bam_aux_append(Record, "MD", 'Z', static_cast<int>(MD.size() + 1), reinterpret_cast<const uint8_t*>(MD.c_str())) == 0 &&
               bam_aux_append(Record, "NH", 'i', sizeof NH, reinterpret_cast<const uint8_t*>(&NH)) == 0 && Add();
    }

    // An unmapped record (FLAG 4): the read as given.
    bool Unmapped(const std::string& Name, const std::string& Seq, const std::string& Qual) {
        return Set(Name, BAM_FUNMAP, -1, -1, 0, 0, nullptr, Seq, Qual) && Add();
    }

    // MD:Z for a read of Length bases with mismatches at (position, reference base 0-3),
    // e.g. "10A5C15"; sorts Mismatches by position.
    static std::string Md_String(int Length, std::vector<std::pair<int, int>>& Mismatches) {
        std::sort(Mismatches.begin(), Mismatches.end());
        std::string MD;
        int At = 0;
        for (const auto& M : Mismatches) {
            MD += std::to_string(M.first - At);
            MD += "ACGT"[M.second & 3];
            At = M.first + 1;
        }
        MD += std::to_string(Length - At);
        return MD;
    }

private:
    friend class Sam_Writer;
    const sam_hdr_t* Header = nullptr;  // set by Sam_Writer::Start
    bool Bam = false;
    std::string Bytes;  // the records built since Start
    bam1_t* Record = nullptr;
    std::string Phred;
    kstring_t Line = KS_INITIALIZE;

    // The record's core, sequence and quality (no tags yet).
    bool Set(const std::string& Name, uint16_t Flag, int Tid, int64_t Pos, uint8_t Mapq, size_t N_Cigar,
             const uint32_t* Cigar, const std::string& Seq, const std::string& Qual) {
        if (!Record && !(Record = bam_init1())) return false;
        const char* Q = nullptr;
        if (!Qual.empty()) {
            Phred.resize(Qual.size());
            for (std::size_t I = 0; I < Qual.size(); I++) Phred[I] = static_cast<char>(Qual[I] - 33);
            Q = Phred.data();
        }
        return bam_set1(Record, Name.size(), Name.c_str(), Flag, Tid, Pos, Mapq, N_Cigar, Cigar, -1, -1, 0, Seq.size(),
                        Seq.c_str(), Q, 32) >= 0;
    }

    // Appends Record to Bytes: a SAM line, or as bam_write1 writes it (block_len, the
    // core's 32 bytes, the data; one CIGAR operation, so never the long-CIGAR form).
    bool Add() {
        if (!Bam) {
            if (sam_format1(Header, Record, &Line) < 0) return false;
            Bytes.append(Line.s, Line.l);
            Bytes += '\n';
            return true;
        }
        const bam1_core_t& C = Record->core;
        if (C.l_qname - C.l_extranul > 255 || C.n_cigar > 0xffff) return false;
        const uint32_t Block_Len = Record->l_data - C.l_extranul + 32;
        const uint32_t X[8] = {static_cast<uint32_t>(C.tid),
                               static_cast<uint32_t>(C.pos),
                               static_cast<uint32_t>(C.bin) << 16 | static_cast<uint32_t>(C.qual) << 8 | (C.l_qname - C.l_extranul),
                               static_cast<uint32_t>(C.flag) << 16 | (C.n_cigar & 0xffff),
                               static_cast<uint32_t>(C.l_qseq),
                               static_cast<uint32_t>(C.mtid),
                               static_cast<uint32_t>(C.mpos),
                               static_cast<uint32_t>(C.isize)};
        Bytes.append(reinterpret_cast<const char*>(&Block_Len), 4);
        Bytes.append(reinterpret_cast<const char*>(X), 32);
        Bytes.append(reinterpret_cast<const char*>(Record->data), C.l_qname - C.l_extranul);
        Bytes.append(reinterpret_cast<const char*>(Record->data) + C.l_qname, Record->l_data - C.l_qname);
        return true;
    }
};

class Sam_Writer {
public:
    Sam_Writer() = default;
    Sam_Writer(const Sam_Writer&) = delete;
    Sam_Writer& operator=(const Sam_Writer&) = delete;
    ~Sam_Writer() { Close(); }

    // Opens Path ("-": stdout) as BAM or SAM and writes the header. False with Error set
    // on failure. BAM is compressed by Threads threads (htslib's pool) if more than one.
    bool Open(const std::string& Path, bool Bam, const Genome_Map& Map, const std::string& Command_Line,
              std::string& Error, int Threads = 1) {
        File = hts_open(Path.c_str(), Bam ? "wb" : "w");
        if (!File) return Fail(Error, "cannot open " + Path + " for writing");
        if (Bam && Threads > 1 && hts_set_threads(File, Threads) < 0) return Fail(Error, "cannot start the BAM threads");
        Header = sam_hdr_init();
        if (!Header || sam_hdr_add_line(Header, "HD", "VN", "1.6", "SO", "unsorted", NULL) < 0)
            return Fail(Error, "cannot build the SAM header");
        for (int I = 0; I < Map.Count(); I++) {
            const std::string Length = std::to_string(Map.Length(I));
            if (sam_hdr_add_line(Header, "SQ", "SN", Map.Name(I).c_str(), "LN", Length.c_str(), NULL) < 0)
                return Fail(Error, "cannot add @SQ " + Map.Name(I));
        }
        if (sam_hdr_add_line(Header, "PG", "ID", "batman_zero_min", "PN", "batman_zero_min", "CL", Command_Line.c_str(),
                             NULL) < 0)
            return Fail(Error, "cannot add @PG");
        if (sam_hdr_write(File, Header) < 0) return Fail(Error, "cannot write the SAM header");
        this->Bam = Bam;
        return true;
    }

    // Starts E on a new set of records for this file.
    void Start(Sam_Encoder& E) const {
        E.Header = Header;
        E.Bam = Bam;
        E.Bytes.clear();
    }

    // Writes E's records (Start, then Mapped / Unmapped), as sam_write1 would have.
    bool Write(const Sam_Encoder& E) {
        const std::string& B = E.Bytes;
        if (!Bam) return B.empty() || hwrite(File->fp.hfile, B.data(), B.size()) == static_cast<ssize_t>(B.size());
        for (std::size_t At = 0; At < B.size();) {  // bam_write1: no record across a block if it fits in one
            uint32_t Block_Len;
            std::memcpy(&Block_Len, B.data() + At, 4);
            const std::size_t N = 4 + static_cast<std::size_t>(Block_Len);
            if (bgzf_flush_try(File->fp.bgzf, N) < 0 || bgzf_write(File->fp.bgzf, B.data() + At, N) != static_cast<ssize_t>(N))
                return false;
            At += N;
        }
        return true;
    }

    // Flushes and closes; false if the output could not be completed.
    bool Close() {
        bool Ok = true;
        if (File) Ok = hts_close(File) == 0;
        if (Header) sam_hdr_destroy(Header);
        File = nullptr;
        Header = nullptr;
        return Ok;
    }

private:
    htsFile* File = nullptr;
    sam_hdr_t* Header = nullptr;
    bool Bam = false;

    static bool Fail(std::string& Error, const std::string& What) {
        Error = What;
        return false;
    }
};
