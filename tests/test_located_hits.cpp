// Differential test for located hits in the mapper: a single-row path that
// reaches a sampled SA row is finished against the packed genome (.pac) instead
// of walking the index. Runs the -n 0 search (Search_Right) and the -n 1 program
// (Batch_One_Mismatch) twice, walking (Use_Ref = false) and located (Use_Ref = true),
// and checks per read and pass that
//   - both give the same hits in the same order, at the same genome positions
//     (walked SA rows decoded with BWTSaValue) and with the same hit counts;
//   - each hit has the same Mismatches, Mismatch_Pos and Mismatch_Char;
//   - every position is at Hamming distance 0 (exact passes) or 1 (mismatch
//     passes) from the read, and the recorded mismatches are where the read and
//     the genome differ, with the genome's bases.
// It also runs the -n 1 to 7 programs (Batch_K: Batch_Search on Program_Of(k)) with
// one lane and with 16, reference off and on, and requires the same result lists
// (Result_List) for every read: hits and deferred ranges, field by field, the hits'
// stages, hit counts and the exact hit count. And with max_hits 0, -n k (k = 1 to 7)
// must report exactly
// the genome positions within Hamming distance k of each read (brute force over the
// genome), each once and with its distance.
// Usage: test_located_hits <genome_prefix> <reads.fa|fq> [max_hits (default 0 = all)]
// The mapper under test: batmis1.cpp, or batmis2.cpp with -DMAPPER_SOURCE='"batmis2.cpp"'.
#ifndef MAPPER_SOURCE
#define MAPPER_SOURCE "batmis1.cpp"
#endif
#define main mapper_main
#include MAPPER_SOURCE
#undef main

namespace {

struct Pass {
    const char* name;
    bool forwardIndex;  // rows in fwfmi (B) rather than revfmi (F)
    unsigned mismatches;
};

unsigned Genome_Base(unsigned long p) { return (Pac[p / 4] >> (6 - 2 * (p % 4))) & 3; }

// Genome start of each location of a hit, sorted.
std::vector<long> Positions(const SARange& hit, bool forwardIndex) {
    std::vector<long> pos;
    if (Use_Ref && hit.Skip) {
        pos.push_back(static_cast<long>(hit.Pos));
        return pos;
    }
    const unsigned long last = hit.Skip ? hit.Start : hit.End;
    for (unsigned long row = hit.Start; row <= last; row++) {
        if (forwardIndex) pos.push_back(BWTSaValue(fwfmi, row));
        else pos.push_back(static_cast<long>(Text_Length) - STRINGLENGTH - static_cast<long>(BWTSaValue(revfmi, row)));
    }
    std::sort(pos.begin(), pos.end());
    return pos;
}

unsigned long failures = 0;

void Fail(const std::string& id, const char* pass, const std::string& what) {
    if (failures++ < 20) std::fprintf(stderr, "FAIL %s %s: %s\n", id.c_str(), pass, what.c_str());
}

// Checks one location against the genome: distance, and the recorded mismatches
// (every row of a range shares them, as they come from one path through the index).
void Check_Location(const std::string& id, const Pass& p, const std::string& read, const SARange& hit, long s) {
    if (s < 0 || s + STRINGLENGTH > static_cast<long>(Text_Length)) {
        Fail(id, p.name, "position " + std::to_string(s) + " off the genome");
        return;
    }
    unsigned d = 0;
    for (int i = 0; i < STRINGLENGTH; i++) d += read[i] != static_cast<char>(Genome_Base(s + i));
    if (d != p.mismatches) Fail(id, p.name, "position " + std::to_string(s) + " at distance " + std::to_string(d));
    if (static_cast<unsigned>(hit.Mismatches) != d) {
        Fail(id, p.name, "Mismatches=" + std::to_string(hit.Mismatches) + " at distance " + std::to_string(d));
        return;
    }
    for (int k = 0; k < hit.Mismatches; k++) {
        const int i = hit.Mismatch_Pos[k];
        const unsigned g = Genome_Base(s + i), c = (hit.Mismatch_Char >> (2 * k)) & 3;
        if (i >= STRINGLENGTH || read[i] == static_cast<char>(g) || c != g)
            Fail(id, p.name, "mismatch " + std::to_string(k) + " recorded at read position " + std::to_string(i));
    }
}

void Compare(const std::string& id, const Pass& p, const std::string& read,
             const std::vector<SARange>& walked, const std::vector<SARange>& located,
             unsigned long& hits, unsigned long& locatedHits) {
    if (walked.size() != located.size()) {
        Fail(id, p.name, std::to_string(walked.size()) + " walked hits, " + std::to_string(located.size()) + " located");
        return;
    }
    for (std::size_t h = 0; h < walked.size(); h++) {
        const SARange& a = walked[h];
        const SARange& b = located[h];
        Use_Ref = false;
        const std::vector<long> pa = Positions(a, p.forwardIndex);
        Use_Ref = true;
        const std::vector<long> pb = Positions(b, p.forwardIndex);
        if (pa != pb) Fail(id, p.name, "hit " + std::to_string(h) + " at different genome positions");
        if (a.Mismatches != b.Mismatches || a.Mismatch_Char != b.Mismatch_Char ||
            std::memcmp(a.Mismatch_Pos, b.Mismatch_Pos, std::min<int>(a.Mismatches, MAX_MISMATCHES_BOUND)) != 0)
            Fail(id, p.name, "hit " + std::to_string(h) + " has different mismatch fields");
        for (long s : pb) Check_Location(id, p, read, b, s);
        hits += pb.size();
        locatedHits += b.Skip ? 1 : 0;
    }
}

// All fields a search sets (Stage only for reported hits: Same_Lists).
bool Same_Range(const SARange& a, const SARange& b) {
    return a.Start == b.Start && a.End == b.End && a.Level == b.Level && a.Mismatches == b.Mismatches &&
           a.Skip == b.Skip && a.Forbidden_Base == b.Forbidden_Base && a.Mismatch_Char == b.Mismatch_Char &&
           a.Pos == b.Pos &&
           std::memcmp(a.Mismatch_Pos, b.Mismatch_Pos, std::min<int>(std::max<int>(a.Mismatches, 0), MAX_MISMATCHES_BOUND)) == 0;
}

bool Same_List(const std::vector<SARange>& a, const std::vector<SARange>& b) {
    if (a.size() != b.size()) return false;
    for (std::size_t i = 0; i < a.size(); i++)
        if (!Same_Range(a[i], b[i])) return false;
    return true;
}

// Two reads' result lists (Result_List), field by field, with the hits' stages.
bool Same_Lists(const ForwardResults* a, const ForwardResults* b) {
    for (unsigned k = 0; k < MM_LISTS; k++)
        if (a[k].hitCount != b[k].hitCount || !Same_List(a[k].hits, b[k].hits) || !Same_List(a[k].mismatches, b[k].mismatches))
            return false;
    for (std::size_t i = 0; i < a[RL_HITS].hits.size(); i++)
        if (a[RL_HITS].hits[i].Stage != b[RL_HITS].hits[i].Stage) return false;
    return true;
}

// The mismatch hits of one stage.
std::vector<SARange> Stage_Hits(const ForwardResults* Res, int Stage) {
    std::vector<SARange> out;
    for (const SARange& h : Res[RL_HITS].hits)
        if (h.Stage == Stage) out.push_back(h);
    return out;
}

// The -n k program (1 to MAX_MISMATCH_SEARCH) on reads [0, n), B lanes.
void Batch_K(int k, const char* codes, unsigned n, unsigned B, unsigned maxHits, std::vector<ForwardResults>& res,
             std::vector<unsigned long>& readHits) {
    Batch_Search(codes, n, STRINGLENGTH, B, maxHits, Program_Of(k), &res, &readHits, nullptr, nullptr);
}

// -n k: one lane against 16 on every read (codes back to back), all lists.
void Check_Batched(const std::vector<std::string>& ids, const std::string& codes, unsigned maxHits, int k) {
    const unsigned n = static_cast<unsigned>(ids.size());
    if (!n || STRINGLENGTH / 2 <= LOOKUPSIZE) return;
    std::vector<ForwardResults> one(MM_LISTS * n), bat(MM_LISTS * n);
    std::vector<unsigned long> oneHits(n), batHits(n);
    for (int mode = 0; mode < 2; mode++) {
        Use_Ref = mode == 1;
        Batch_K(k, codes.data(), n, 1, maxHits, one, oneHits);
        Batch_K(k, codes.data(), n, 16, maxHits, bat, batHits);
        for (unsigned i = 0; i < n; i++)
            if (oneHits[i] != batHits[i] || !Same_Lists(&one[MM_LISTS * i], &bat[MM_LISTS * i]))
                Fail(ids[i], ("-n " + std::to_string(k) + (Use_Ref ? " batch ref" : " batch no-ref")).c_str(),
                     "B=1 and B=16 results differ");
    }
}

// -n k with no hit limit: the reported positions, by mismatch count, are exactly the
// genome positions within Hamming distance k (brute force). Reference off and on.
void Check_Complete(const std::vector<std::string>& ids, const std::string& codes, int k) {
    const unsigned n = static_cast<unsigned>(ids.size());
    if (!n || STRINGLENGTH / 2 <= LOOKUPSIZE) return;
    const Hit_Stage* Stages = Stages_Of(k);
    std::vector<ForwardResults> res(MM_LISTS * n);
    std::vector<unsigned long> readHits(n);
    const long last = static_cast<long>(Text_Length) - STRINGLENGTH;
    for (int mode = 0; mode < 2; mode++) {
        Use_Ref = mode == 1;
        Batch_K(k, codes.data(), n, 16, 0, res, readHits);
        for (unsigned i = 0; i < n; i++) {
            const std::string_view all(codes.data() + static_cast<std::size_t>(i) * STRINGLENGTH, STRINGLENGTH);
            const ForwardResults* r = &res[MM_LISTS * i];
            std::vector<long> got[MAX_MISMATCH_SEARCH + 1], want[MAX_MISMATCH_SEARCH + 1];
            for (const SARange& h : r[RL_EXACT].hits)
                for (long p : Positions(h, false)) got[0].push_back(p);
            for (const SARange& h : r[RL_HITS].hits)
                for (long p : Positions(h, Stages[h.Stage].Dir[0] == 'B')) got[Stages[h.Stage].Mismatches].push_back(p);
            for (long s = 0; s <= last; s++) {
                int d = 0;
                for (int j = 0; j < STRINGLENGTH && d <= k; j++) d += all[j] != static_cast<char>(Genome_Base(s + j));
                if (d <= k) want[d].push_back(s);
            }
            for (int d = 0; d <= k; d++) {
                std::sort(got[d].begin(), got[d].end());
                if (got[d] != want[d])
                    Fail(ids[i], ("-n " + std::to_string(k) + (Use_Ref ? " ref" : " no-ref")).c_str(),
                         std::to_string(got[d].size()) + " positions at distance " + std::to_string(d) + ", brute force " +
                             std::to_string(want[d].size()));
            }
        }
    }
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 3) {
        std::fprintf(stderr, "Usage: %s <genome_prefix> <reads.fa|fq> [max_hits]\n", argv[0]);
        return 2;
    }
    const std::string prefix = argv[1];
    const unsigned maxHits = argc > 3 ? static_cast<unsigned>(std::strtoul(argv[3], nullptr, 10)) : 0;
    std::string revBwt, revOcc, fwBwt, fwOcc, revSa, fwSa;
    BuildIndexPaths(prefix, revBwt, revOcc, fwBwt, fwOcc, &revSa, &fwSa);
    for (const std::string* name : {&revSa, &fwSa}) {
        if (!FileExists(*name)) {
            std::fprintf(stderr, "missing %s\n", name->c_str());
            return 2;
        }
    }
    revfmi = initFMI(revBwt.c_str(), revOcc.c_str(), revSa.c_str());
    fwfmi = BWTLoad(mmPool, fwBwt.c_str(), fwOcc.c_str(), fwSa.c_str(), nullptr, nullptr, nullptr);
    if (!Set_Sampling(revfmi, fwfmi) || !LoadPac(prefix + ".pac")) {
        std::fprintf(stderr, "needs %s.pac and SA samples at one power-of-two interval in both indexes\n", prefix.c_str());
        return 2;
    }
    Text_Length = revfmi->textLength;

    Read_File in(argv[2]);
    std::string id, seq, read;
    bool initialized = false;
    unsigned long reads = 0, hits = 0, locatedHits = 0;
    std::vector<std::string> batchIds;
    std::string batchCodes;
    const Pass exactPass{"-n0", false, 0}, pass00{"0|0", false, 0}, pass01{"0|1", false, 1}, pass10{"1|0", true, 1};
    while (ReadNextRecord(in, id, seq)) {
        if (!initialized) {
            STRINGLENGTH = static_cast<int>(seq.size());
            Build_Tables(revfmi);
            read.resize(STRINGLENGTH);
            initialized = true;
        }
        if (static_cast<int>(seq.size()) != STRINGLENGTH) continue;
        bool valid = true;
        for (int i = 0; i < STRINGLENGTH && valid; i++) {
            const int code = EncodeBase(seq[i]);
            valid = code >= 0;
            read[i] = static_cast<char>(code);
        }
        if (!valid) continue;
        reads++;
        batchIds.push_back(id);
        batchCodes += read;
        Pack_Read(read);

        const std::string_view all(read);
        std::vector<SARange> exact[2];
        for (int mode = 0; mode < 2; mode++) {
            Use_Ref = mode == 1;
            SARange r = Start_Range(all);
            if ((r.Level > STRINGLENGTH || Search_Right(all, 0, r, revfmi)) && r.Start != 0) exact[mode].push_back(r);
        }
        Compare(id, exactPass, read, exact[0], exact[1], hits, locatedHits);
    }
    // The -n 1 passes walked and located (the program's stages 0 and 1: 0|1 and 1|0).
    const unsigned n = static_cast<unsigned>(batchIds.size());
    if (n && STRINGLENGTH / 2 > LOOKUPSIZE) {
        std::vector<ForwardResults> walked(MM_LISTS * n), located(MM_LISTS * n);
        std::vector<unsigned long> walkedHits(n), locatedCount(n);
        Use_Ref = false;
        Batch_One_Mismatch(batchCodes.data(), n, STRINGLENGTH, 16, maxHits, walked, walkedHits);
        Use_Ref = true;
        Batch_One_Mismatch(batchCodes.data(), n, STRINGLENGTH, 16, maxHits, located, locatedCount);
        for (unsigned i = 0; i < n; i++) {
            const std::string read(batchCodes, static_cast<std::size_t>(i) * STRINGLENGTH, STRINGLENGTH);
            const ForwardResults *w = &walked[MM_LISTS * i], *l = &located[MM_LISTS * i];
            Compare(batchIds[i], pass00, read, w[RL_EXACT].hits, l[RL_EXACT].hits, hits, locatedHits);
            Compare(batchIds[i], pass01, read, Stage_Hits(w, 0), Stage_Hits(l, 0), hits, locatedHits);
            Compare(batchIds[i], pass10, read, Stage_Hits(w, 1), Stage_Hits(l, 1), hits, locatedHits);
            for (unsigned k : {RL_EXACT, RL_HITS})
                if (w[k].hitCount != l[k].hitCount) Fail(batchIds[i], "count", "different hit counts");
        }
    }
    for (int k = 1; k <= MAX_MISMATCH_SEARCH; k++) Check_Batched(batchIds, batchCodes, maxHits, k);
    if (maxHits == 0)
        for (int k = 1; k <= MAX_MISMATCH_SEARCH; k++) Check_Complete(batchIds, batchCodes, k);
    std::printf("reads=%lu length=%d max_hits=%u locations=%lu located_hits=%lu failures=%lu\n",
                reads, STRINGLENGTH, maxHits, hits, locatedHits, failures);
    return failures ? 1 : 0;
}
