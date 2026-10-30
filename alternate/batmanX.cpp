#include "config.h"

#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

extern "C" {
#include "BWT.h"
#include "MemManager.h"
}

#define SAINTERVAL 8

namespace {

static const int kSmallRangeThreshold = 8;
static const int kLookupSize = 6;

static const char *BuildSimdMode() {
    return "MMX/SSE2";
}

struct LookupRange {
    unsigned int start;
    unsigned int end;
};

struct PrefixLookupTable {
    std::vector<LookupRange> ranges;
};

struct SearchState {
    unsigned int start;
    unsigned int end;
    int level;
    int mismatches;
    int mismatchPos[2];
    unsigned int mismatchBits;
    bool skip;
};

struct ExactForwardState {
    unsigned int start;
    unsigned int end;
    int level;
    unsigned char skip;
};

struct OneMismatchStateLite {
    unsigned int start;
    unsigned int end;
    int level;
    unsigned char mismatches;
    unsigned char skip;
};

static unsigned int BaseToCode(char base) {
    switch (base) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default: return 0;
    }
}

static unsigned int BwtCharAt(const BWT *fmi, unsigned int pos) {
    unsigned int idx = pos;
    if (idx >= fmi->inverseSa0) {
        idx--;
    }
    unsigned int word = fmi->bwtCode[idx / 16];
    return (word << (((idx % 16) * 2))) >> (BITS_IN_WORD - 2);
}

static bool AdvanceRange(const BWT *fmi, SearchState &state, unsigned int nextBase) {
    if (state.start == 0) {
        return false;
    }
    unsigned int nextStart = fmi->cumulativeFreq[nextBase] + BWTOccValue(fmi, state.start, nextBase) + 1;
    unsigned int nextEnd = fmi->cumulativeFreq[nextBase] + BWTOccValue(fmi, state.end + 1, nextBase);
    if (nextEnd < nextStart) {
        state.start = 0;
        state.end = 0;
        return false;
    }
    state.start = nextStart;
    state.end = nextEnd;
    return true;
}

static bool AdvanceSmallRange(const BWT *fmi, SearchState &state, unsigned int expectedBase) {
    if (state.start == 0 || state.end < state.start) {
        return false;
    }

    unsigned int branchCounts[4] = {0, 0, 0, 0};
    const unsigned int first = (state.start + 1 >= fmi->inverseSa0) ? state.start : state.start + 1;
    const unsigned int last = (state.start + 1 >= fmi->inverseSa0) ? state.end : state.end + 1;

    for (unsigned long pos = first; pos <= last; ++pos) {
        const unsigned int idx = pos - 1;
        const unsigned int now = (fmi->bwtCode[idx / 16] << ((idx % 16) * 2)) >> (BITS_IN_WORD - 2);
        branchCounts[now]++;
    }

    if (branchCounts[expectedBase] == 0) {
        state.start = 0;
        state.end = 0;
        return false;
    }

    state.start = fmi->cumulativeFreq[expectedBase] + BWTOccValue(fmi, state.start, expectedBase) + 1;
    state.end = state.start + branchCounts[expectedBase] - 1;
    return true;
}

static bool AdvanceSingleBranch(const BWT *fmi, SearchState &state, unsigned int expectedBase) {
    if (state.start == 0) {
        return false;
    }

    const unsigned int now = BwtCharAt(fmi, state.start);
    if (now != expectedBase) {
        state.start = 0;
        state.end = 0;
        return false;
    }

    state.start = fmi->cumulativeFreq[now] + BWTOccValue(fmi, state.start, now) + 1;
    state.end = state.start;
    /*if (state.skip) {
        state.skip = true;
    } else if (state.start % SAINTERVAL == 0) {
        state.skip = true;
        state.end = state.start;
    }*/

    return true;
}

static void EmitHits(const std::vector<SearchState> &hits, unsigned int maxHits) {
    unsigned int count = 0;
    for (std::vector<SearchState>::const_iterator it = hits.begin(); it != hits.end() && count < maxHits; ++it) {
        const SearchState &s = *it;
        if (s.start == 0 || s.end < s.start) {
            continue;
        }
        std::cout << "SA[" << s.start << ", " << s.end << "] mismatch=" << s.mismatches << "\n";
        ++count;
    }
}

static PrefixLookupTable BuildPrefixLookupTable(const BWT *fmi) {
    PrefixLookupTable table;
    const std::size_t tableSize = static_cast<std::size_t>(1u) << (2 * kLookupSize);
    table.ranges.resize(tableSize);

    for (std::size_t key = 0; key < tableSize; ++key) {
        SearchState state = {1, fmi->textLength, 0, 0, {0, 0}, 0, false};
        bool ok = true;
        for (int level = 0; level < kLookupSize; ++level) {
            const unsigned int base = static_cast<unsigned int>((key >> (2 * level)) & 0x3u);
            if (!AdvanceRange(fmi, state, base)) {
                ok = false;
                break;
            }
        }

        if (ok) {
            table.ranges[key].start = state.start;
            table.ranges[key].end = state.end;
        } else {
            table.ranges[key].start = 0;
            table.ranges[key].end = 0;
        }
    }

    return table;
}

static bool SeedStateFromLookup(const PrefixLookupTable &table,
                                const std::string &query,
                                SearchState &state) {
    if (query.size() < static_cast<std::size_t>(kLookupSize)) {
        return false;
    }

    unsigned int key = 0;
    for (int i = 0; i < kLookupSize; ++i) {
        key |= (BaseToCode(query[static_cast<std::size_t>(i)]) << (2 * i));
    }

    if (key >= table.ranges.size()) {
        return false;
    }

    const LookupRange &seed = table.ranges[key];
    if (seed.start == 0 || seed.end < seed.start) {
        state.start = 0;
        state.end = 0;
        state.level = kLookupSize;
        return true;
    }

    state.start = seed.start;
    state.end = seed.end;
    state.level = kLookupSize;
    return true;
}

static bool SeedExactStateFromLookup(const PrefixLookupTable &table,
                                     const std::string &query,
                                     ExactForwardState &state) {
    if (query.size() < static_cast<std::size_t>(kLookupSize)) {
        return false;
    }

    unsigned int key = 0;
    for (int i = 0; i < kLookupSize; ++i) {
        key |= (BaseToCode(query[static_cast<std::size_t>(i)]) << (2 * i));
    }

    if (key >= table.ranges.size()) {
        return false;
    }

    const LookupRange &seed = table.ranges[key];
    if (seed.start == 0 || seed.end < seed.start) {
        state.start = 0;
        state.end = 0;
        state.level = kLookupSize + 1;
        state.skip = 0;
        return true;
    }

    state.start = seed.start;
    state.end = seed.end;
    state.level = kLookupSize + 1;
    state.skip = 0;
    return true;
}

static bool SearchExactForwardFast(const BWT *fmi,
                                   const PrefixLookupTable &lookupTable,
                                   const std::string &query,
                                   ExactForwardState &state) {
    if (query.empty()) {
        return false;
    }

    // For tiny queries, keep the generic path.
    if (!SeedExactStateFromLookup(lookupTable, query, state)) {
        return false;
    }
    if (state.start == 0) {
        return false;
    }

    for (;;) {
        if (state.end == state.start || state.skip) {
            if (state.start % SAINTERVAL == 0 && !state.skip) {
                state.skip++;
                state.end = state.start;
            }

            for (;;) {
                unsigned int index = state.start;
                if (index >= fmi->inverseSa0) {
                    index--;
                }

                const unsigned int now = (fmi->bwtCode[index / 16] << ((index % 16) * 2)) >> (BITS_IN_WORD - 2);
                const unsigned int expected = BaseToCode(query[static_cast<std::size_t>(state.level - 1)]);

                if (expected == now) {
                    state.start = fmi->cumulativeFreq[now] + BWTOccValue(fmi, state.start, now) + 1;
                    if (state.skip) {
                        state.skip++;
                    } else if (state.start % SAINTERVAL == 0) {
                        state.skip++;
                        state.end = state.start;
                    }

                    if (state.level == static_cast<int>(query.size())) {
                        state.end = state.start;
                        return true;
                    }
                    state.level++;
                    continue;
                }
                return false;
            }
        }

        const unsigned int expected = BaseToCode(query[static_cast<std::size_t>(state.level - 1)]);
        state.start = fmi->cumulativeFreq[expected] + BWTOccValue(fmi, state.start, expected) + 1;
        state.end = fmi->cumulativeFreq[expected] + BWTOccValue(fmi, state.end + 1, expected);
        if (state.end < state.start) {
            state.start = 0;
        }

        if (state.start != 0) {
            if (state.level == static_cast<int>(query.size())) {
                return true;
            }
            state.level++;
            continue;
        }
        return false;
    }
}

static bool SearchOneMismatchSingleHitFast(const BWT *fmi,
                                           const PrefixLookupTable &lookupTable,
                                           const std::string &query,
                                           const SearchState &root,
                                           OneMismatchStateLite &result) {
    SearchState seededRoot = root;
    SeedStateFromLookup(lookupTable, query, seededRoot);

    if (seededRoot.start == 0 || seededRoot.end < seededRoot.start) {
        return false;
    }

    std::vector<OneMismatchStateLite> stack;
    stack.reserve(query.size() * 4);
    stack.push_back({seededRoot.start, seededRoot.end, seededRoot.level, 0, 0});

    while (!stack.empty()) {
        OneMismatchStateLite state = stack.back();
        stack.pop_back();

        if (state.start == 0 || state.end < state.start) {
            continue;
        }
        if (state.level == static_cast<int>(query.size())) {
            result = state;
            return true;
        }

        const unsigned int expected = BaseToCode(query[static_cast<std::size_t>(state.level)]);

        if (state.end == state.start || state.skip) {
            OneMismatchStateLite exact = state;
            exact.level = state.level + 1;
            {
                const unsigned int now = BwtCharAt(fmi, exact.start);
                if (now == expected) {
                    exact.start = fmi->cumulativeFreq[now] + BWTOccValue(fmi, exact.start, now) + 1;
                    exact.end = exact.start;
                    stack.push_back(exact);
                }
            }

            if (state.mismatches >= 1) {
                continue;
            }

            for (unsigned int branch = 0; branch < 4; ++branch) {
                if (branch == expected) {
                    continue;
                }

                OneMismatchStateLite mismatch = state;
                mismatch.level = state.level + 1;
                mismatch.mismatches = static_cast<unsigned char>(state.mismatches + 1);

                const unsigned int now = BwtCharAt(fmi, mismatch.start);
                if (now != branch) {
                    continue;
                }
                mismatch.start = fmi->cumulativeFreq[now] + BWTOccValue(fmi, mismatch.start, now) + 1;
                mismatch.end = mismatch.start;
                stack.push_back(mismatch);
            }
            continue;
        }

        if ((state.end - state.start + 1) <= kSmallRangeThreshold) {
            unsigned int branchCounts[4] = {0, 0, 0, 0};
            const unsigned int first = (state.start + 1 >= fmi->inverseSa0) ? state.start : state.start + 1;
            const unsigned int last = (state.start + 1 >= fmi->inverseSa0) ? state.end : state.end + 1;
            for (unsigned long pos = first; pos <= last; ++pos) {
                const unsigned int idx = static_cast<unsigned int>(pos - 1);
                const unsigned int now = (fmi->bwtCode[idx / 16] << ((idx % 16) * 2)) >> (BITS_IN_WORD - 2);
                branchCounts[now]++;
            }

            if (branchCounts[expected] > 0) {
                OneMismatchStateLite exact = state;
                exact.level = state.level + 1;
                exact.start = fmi->cumulativeFreq[expected] + BWTOccValue(fmi, exact.start, expected) + 1;
                exact.end = exact.start + branchCounts[expected] - 1;
                stack.push_back(exact);
            }

            if (state.mismatches >= 1) {
                continue;
            }

            for (unsigned int branch = 0; branch < 4; ++branch) {
                if (branch == expected || branchCounts[branch] == 0) {
                    continue;
                }
                OneMismatchStateLite mismatch = state;
                mismatch.level = state.level + 1;
                mismatch.mismatches = static_cast<unsigned char>(state.mismatches + 1);
                mismatch.start = fmi->cumulativeFreq[branch] + BWTOccValue(fmi, mismatch.start, branch) + 1;
                mismatch.end = mismatch.start + branchCounts[branch] - 1;
                stack.push_back(mismatch);
            }
            continue;
        }

        {
            OneMismatchStateLite exact = state;
            exact.level = state.level + 1;
            exact.start = fmi->cumulativeFreq[expected] + BWTOccValue(fmi, exact.start, expected) + 1;
            exact.end = fmi->cumulativeFreq[expected] + BWTOccValue(fmi, exact.end + 1, expected);
            if (exact.end >= exact.start) {
                stack.push_back(exact);
            }
        }

        if (state.mismatches >= 1) {
            continue;
        }

        for (unsigned int branch = 0; branch < 4; ++branch) {
            if (branch == expected) {
                continue;
            }

            OneMismatchStateLite mismatch = state;
            mismatch.level = state.level + 1;
            mismatch.mismatches = static_cast<unsigned char>(state.mismatches + 1);
            mismatch.start = fmi->cumulativeFreq[branch] + BWTOccValue(fmi, mismatch.start, branch) + 1;
            mismatch.end = fmi->cumulativeFreq[branch] + BWTOccValue(fmi, mismatch.end + 1, branch);
            if (mismatch.end >= mismatch.start) {
                stack.push_back(mismatch);
            }
        }
    }

    return false;
}

static void SearchZeroMismatch(const BWT *fmi,
                              const std::string &query,
                              SearchState root,
                              std::vector<SearchState> &hits,
                              unsigned int maxHits) {
    std::vector<SearchState> stack;
    stack.push_back(root);

    while (!stack.empty() && hits.size() < maxHits) {
        SearchState state = stack.back();
        stack.pop_back();

        if (state.start == 0 || state.end < state.start) {
            continue;
        }
        if (state.level == static_cast<int>(query.size())) {
            hits.push_back(state);
            continue;
        }

        const unsigned int expected = BaseToCode(query[state.level]);

        if (state.end == state.start || state.skip) {
            SearchState next = state;
            next.level = state.level + 1;
            if (AdvanceSingleBranch(fmi, next, expected)) {
                stack.push_back(next);
            }
            continue;
        }

        if ((state.end - state.start + 1) <= kSmallRangeThreshold) {
            SearchState next = state;
            next.level = state.level + 1;
            if (AdvanceSmallRange(fmi, next, expected)) {
                stack.push_back(next);
            }
            continue;
        }

        SearchState exact = state;
        exact.level = state.level + 1;
        if (AdvanceRange(fmi, exact, expected)) {
            stack.push_back(exact);
        }
    }
}

static void SearchOneMismatch(const BWT *fmi,
                             const std::string &query,
                             SearchState root,
                             std::vector<SearchState> &hits,
                             unsigned int maxHits) {
    std::vector<SearchState> stack;
    stack.push_back(root);

    while (!stack.empty() && hits.size() < maxHits) {
        SearchState state = stack.back();
        stack.pop_back();

        if (state.start == 0 || state.end < state.start) {
            continue;
        }
        if (state.level == static_cast<int>(query.size())) {
            if (state.mismatches <= 1) {
                hits.push_back(state);
            }
            continue;
        }

        const unsigned int expected = BaseToCode(query[state.level]);

        if (state.end == state.start || state.skip) {
            SearchState exact = state;
            exact.level = state.level + 1;
            if (AdvanceSingleBranch(fmi, exact, expected)) {
                stack.push_back(exact);
            }
            if (state.mismatches >= 1) {
                continue;
            }
            for (unsigned int branch = 0; branch < 4; ++branch) {
                if (branch == expected) {
                    continue;
                }
                SearchState mismatch = state;
                mismatch.level = state.level + 1;
                mismatch.mismatchPos[mismatch.mismatches] = state.level;
                mismatch.mismatchBits |= (branch << (mismatch.mismatches * 2));
                mismatch.mismatches++;
                if (AdvanceSingleBranch(fmi, mismatch, branch)) {
                    stack.push_back(mismatch);
                }
            }
            continue;
        }

        if ((state.end - state.start + 1) <= kSmallRangeThreshold) {
            SearchState exact = state;
            exact.level = state.level + 1;
            if (AdvanceSmallRange(fmi, exact, expected)) {
                stack.push_back(exact);
            }
            if (state.mismatches >= 1) {
                continue;
            }
            for (unsigned int branch = 0; branch < 4; ++branch) {
                if (branch == expected) {
                    continue;
                }
                SearchState mismatch = state;
                mismatch.level = state.level + 1;
                mismatch.mismatchPos[mismatch.mismatches] = state.level;
                mismatch.mismatchBits |= (branch << (mismatch.mismatches * 2));
                mismatch.mismatches++;
                if (AdvanceSmallRange(fmi, mismatch, branch)) {
                    stack.push_back(mismatch);
                }
            }
            continue;
        }

        SearchState exact = state;
        exact.level = state.level + 1;
        if (AdvanceRange(fmi, exact, expected)) {
            stack.push_back(exact);
        }

        if (state.mismatches >= 1) {
            continue;
        }

        for (unsigned int branch = 0; branch < 4; ++branch) {
            if (branch == expected) {
                continue;
            }
            SearchState mismatch = state;
            mismatch.level = state.level + 1;
            mismatch.mismatchPos[mismatch.mismatches] = state.level;
            mismatch.mismatchBits |= (branch << (mismatch.mismatches * 2));
            mismatch.mismatches++;
            if (AdvanceRange(fmi, mismatch, branch)) {
                stack.push_back(mismatch);
            }
        }
    }
}

static std::string CanonicalizeQuery(const std::string &input) {
    std::string out;
    for (std::size_t i = 0; i < input.size(); ++i) {
        const char c = input[i];
        if (c == '\n' || c == '\r' || c == ' ' || c == '\t') {
            continue;
        }
        out.push_back(c);
    }
    return out;
}

static void SearchQuery(BWT *fmi,
                        const PrefixLookupTable &lookupTable,
                        const std::string &query,
                        const SearchState &root,
                        unsigned int maxHits,
                        int mismatchLimit,
                        bool forceExpensiveExtension,
                        bool zigzag,
                        bool quiet,
                        std::size_t queryIndex) {
    std::string currentQuery = CanonicalizeQuery(query);
    if (currentQuery.empty()) {
        return;
    }

    const bool useExactFastPath = !forceExpensiveExtension &&
                                  mismatchLimit == 0 &&
                                  currentQuery.size() >= static_cast<std::size_t>(kLookupSize);
    const bool useOneMismatchSingleHitFastPath = !forceExpensiveExtension &&
                                                 mismatchLimit == 1 &&
                                                 maxHits == 1 &&
                                                 currentQuery.size() >= static_cast<std::size_t>(kLookupSize);

    if (useExactFastPath) {
        ExactForwardState exactState = {0, 0, 0, 0};
        const bool found = SearchExactForwardFast(fmi, lookupTable, currentQuery, exactState);

        if (!quiet) {
            std::cout << "query[" << queryIndex << "]=" << currentQuery << "\n";
            if (zigzag) {
                std::cout << "zigzag\n";
            }
            if (found && exactState.start != 0 && maxHits > 0) {
                std::cout << "SA[" << exactState.start << ", " << exactState.end << "] mismatch=0\n";
            }
        }

        return;
    }

    if (useOneMismatchSingleHitFastPath) {
        OneMismatchStateLite hit = {0, 0, 0, 0, 0};
        const bool found = SearchOneMismatchSingleHitFast(fmi, lookupTable, currentQuery, root, hit);

        if (!quiet) {
            std::cout << "query[" << queryIndex << "]=" << currentQuery << "\n";
            if (zigzag) {
                std::cout << "zigzag\n";
            }
            if (found && hit.start != 0) {
                std::cout << "SA[" << hit.start << ", " << hit.end << "] mismatch=" << static_cast<int>(hit.mismatches) << "\n";
            }
        }

        return;
    }

    SearchState seededRoot = root;
    if (!forceExpensiveExtension) {
        SeedStateFromLookup(lookupTable, currentQuery, seededRoot);
    }

    std::vector<SearchState> hits;
    if (mismatchLimit == 0) {
        SearchZeroMismatch(fmi, currentQuery, seededRoot, hits, maxHits);
    } else {
        SearchOneMismatch(fmi, currentQuery, seededRoot, hits, maxHits);
    }

    if (!quiet) {
        std::cout << "query[" << queryIndex << "]=" << currentQuery << "\n";
        if (zigzag) {
            std::cout << "zigzag\n";
        }
        EmitHits(hits, maxHits);
    }

}

static void SearchFastqFile(BWT *fmi,
                           const PrefixLookupTable &lookupTable,
                           const std::string &fastqPath,
                           const SearchState &root,
                           unsigned int maxHits,
                           int mismatchLimit,
                           bool forceExpensiveExtension,
                           bool zigzag,
                           bool quiet) {
    std::ifstream in(fastqPath.c_str());
    if (!in) {
        std::cerr << "failed to open FASTQ file: " << fastqPath << "\n";
        return;
    }

    std::string header;
    std::string seq;
    std::string plus;
    std::string quality;
    std::size_t queryIndex = 0;

    while (std::getline(in, header)) {
        if (header.empty() || header[0] != '@') {
            continue;
        }
        if (!std::getline(in, seq)) {
            break;
        }
        if (!std::getline(in, plus)) {
            break;
        }
        if (!std::getline(in, quality)) {
            break;
        }

        SearchQuery(fmi, lookupTable, seq, root, maxHits, mismatchLimit, forceExpensiveExtension, zigzag, quiet, queryIndex++);
    }
}

static void Usage(const char *name) {
    std::cout << "Usage: " << name << " --bwt <bwt-code> --occ <occ-file> [--query <seq>] [--maxhits N] [--mismatch 0|1] [--zigzag]\n";
    std::cout << "       " << name << " [--quiet] ...\n";
    std::cout << "       " << name << " [--force-expensive-extension] ...\n";
    std::cout << "  The compact zero/one mismatch core keeps the direct small-range BWT walk optimization.\n";
    std::cout << "  --force-expensive-extension disables fast shortcuts and lookup seeding for benchmarking.\n";
    std::cout << "  Build SIMD mode: " << BuildSimdMode() << " (use --build-info for script-friendly output).\n";
}

}  // namespace

int main(int argc, char **argv) {
    std::string bwtPath;
    std::string occPath;
    std::string query;
    unsigned int maxHits = 10;
    int mismatchLimit = 1;
    bool zigzag = true;
    bool quiet = false;
    bool forceExpensiveExtension = false;

    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--build-info") {
            std::cout << "batmanX SIMD mode: " << BuildSimdMode() << "\n";
            return 0;
        }
        if (arg == "--bwt" || arg == "-b") {
            if (i + 1 < argc) bwtPath = argv[++i];
        } else if (arg == "--occ" || arg == "-o") {
            if (i + 1 < argc) occPath = argv[++i];
        } else if (arg == "--query" || arg == "-q") {
            if (i + 1 < argc) query = argv[++i];
        } else if (arg == "--maxhits" || arg == "-m") {
            if (i + 1 < argc) maxHits = static_cast<unsigned int>(std::atoi(argv[++i]));
        } else if (arg == "--mismatch" || arg == "-n") {
            if (i + 1 < argc) mismatchLimit = std::atoi(argv[++i]);
        } else if (arg == "--zigzag" || arg == "-z") {
            zigzag = true;
        } else if (arg == "--no-zigzag") {
            zigzag = false;
        } else if (arg == "--quiet") {
            quiet = true;
        } else if (arg == "--force-expensive-extension") {
            forceExpensiveExtension = true;
        } else if (arg == "--help" || arg == "-h") {
            Usage(argv[0]);
            return 0;
        }
    }

    if (bwtPath.empty() || occPath.empty()) {
        Usage(argv[0]);
        return 1;
    }
    if (query.empty()) {
        query = "ACGT";
    }

    if (mismatchLimit < 0 || mismatchLimit > 1) {
        mismatchLimit = 1;
    }

    MMMasterInitialize(3, 0, FALSE, NULL);
    MMPool *mmPool = MMPoolCreate(1 << 20);

    BWT *fmi = BWTLoad(mmPool, bwtPath.c_str(), occPath.c_str(), NULL, NULL, NULL, NULL);
    if (fmi == NULL) {
        std::cerr << "failed to load FM-index: " << bwtPath << " / " << occPath << "\n";
        return 2;
    }

    const unsigned int totalLength = fmi->textLength;
    SearchState root = {1, totalLength, 0, 0, {0, 0}, 0, false};
    PrefixLookupTable lookupTable = BuildPrefixLookupTable(fmi);

    const bool isFastqInput = (query.size() >= 3 && query.find(".fq") != std::string::npos) ||
                              (query.size() >= 5 && query.find(".fastq") != std::string::npos);

    if (isFastqInput) {
        SearchFastqFile(fmi, lookupTable, query, root, maxHits, mismatchLimit, forceExpensiveExtension, zigzag, quiet);
    } else {
        const std::string currentQuery = CanonicalizeQuery(query);
        if (currentQuery.empty()) {
            std::cerr << "query is empty after filtering\n";
            BWTFree(mmPool, fmi);
            return 1;
        }
        SearchQuery(fmi, lookupTable, currentQuery, root, maxHits, mismatchLimit, forceExpensiveExtension, zigzag, quiet, 0);
    }

    BWTFree(mmPool, fmi);
    return 0;
}
