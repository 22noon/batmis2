#include "exact_search.hpp"

#include <cassert>

namespace batmis {
namespace {
constexpr int kSaInterval = 8;
constexpr int kBranchThreshold = 0;
}

ExactSearch::ExactSearch(const BWT* index) : index_(index) {}

void ExactSearch::extend(unsigned long& start, unsigned long& end,
                         unsigned char base) const {
    start = index_->cumulativeFreq[base] + BWTOccValue(index_, start, base) + 1;
    end = index_->cumulativeFreq[base] + BWTOccValue(index_, end + 1, base);
    if (end < start) {
        start = 0;
    }
}

ExactSearch::SeedRange ExactSearch::initialBaseRange(int base) const {
    SeedRange range{};
    range.Start = index_->cumulativeFreq[base] + 1;
    range.End = index_->cumulativeFreq[base + 1];
    if (range.Start > range.End) range.Start = 0;
    return range;
}

void ExactSearch::buildSeed(SeedRange range, int level, int base) {
    range.Label |= static_cast<unsigned>(base) << (2 * (level - 1));
    if (range.Start) {
        unsigned long start = range.Start;
        unsigned long end = range.End;
        extend(start, end, static_cast<unsigned char>(base));
        range.Start = static_cast<unsigned>(start);
        range.End = static_cast<unsigned>(end);
    }

    if (level == lookupSize_) {
        assert(range.Start <= range.End);
        seedStart_[range.Label] = range.Start;
        seedEnd_[range.Label] = range.End;
        return;
    }

    for (int nextBase = 0; nextBase < 4; ++nextBase) {
        buildSeed(range, level + 1, nextBase);
    }
}

bool ExactSearch::prepare(int readLength) {
    if (!index_ || readLength <= 0 || readLength >= 4096) return false;
    readLength_ = readLength;
    lookupSize_ = readLength < 28 ? 3 : 6;

    const int seedCount = 1 << (2 * lookupSize_);
    for (int key = 0; key < seedCount; ++key) {
        seedStart_[key] = 0;
        seedEnd_[key] = 0;
    }

    for (int base = 0; base < 4; ++base) {
        SeedRange range = initialBaseRange(base);
        range.Label = static_cast<unsigned>(base);
        for (int nextBase = 0; nextBase < 4; ++nextBase) {
            buildSeed(range, 2, nextBase);
        }
    }
    return true;
}

bool ExactSearch::searchExact(const unsigned char* read, SearchRange& range) const {
    unsigned long first, last;
    unsigned int now;

    if (!range.Start) return false;

    for (;;) {
        if (range.End == range.Start || range.Skip) {
            if (range.Start % kSaInterval == 0 && !range.Skip) {
                ++range.Skip;
                range.End = range.Start;
            }

            for (;;) {
                unsigned long index = range.Start;
                if (index >= index_->inverseSa0) --index;
                now = index_->bwtCode[index / 16] << ((index % 16) * 2) >> (BITS_IN_WORD - 2);
                if (read[range.Level - 1] == now) {
                    range.Start = index_->cumulativeFreq[now] +
                                  BWTOccValue(index_, range.Start, now) + 1;
                    if (range.Skip) {
                        ++range.Skip;
                    } else if (range.Start % kSaInterval == 0) {
                        ++range.Skip;
                        range.End = range.Start;
                    }

                    if (range.Level == readLength_) {
                        range.End = range.Start;
                        return true;
                    }
                    ++range.Level;
                    continue;
                }
                return false;
            }
        } else {
            if (range.End - range.Start < kBranchThreshold) {
                unsigned int branchCharacters[4] = {0, 0, 0, 0};
                if (range.Start + 1 >= index_->inverseSa0) {
                    first = range.Start;
                    last = range.End;
                } else {
                    first = range.Start + 1;
                    last = range.End + 1;
                }

                for (unsigned long pos = first; pos <= last; ++pos) {
                    now = index_->bwtCode[(pos - 1) / 16] << (((pos - 1) % 16) * 2) >> (BITS_IN_WORD - 2);
                    ++branchCharacters[now];
                }

                now = read[range.Level - 1];
                if (branchCharacters[now]) {
                    range.Start = index_->cumulativeFreq[now] +
                                  BWTOccValue(index_, range.Start, now) + 1;
                    range.End = range.Start + branchCharacters[now] - 1;
                } else {
                    range.Start = 0;
                }
            } else {
                const unsigned char base = read[range.Level - 1];
                range.Start = index_->cumulativeFreq[base] +
                              BWTOccValue(index_, range.Start, base) + 1;
                range.End = index_->cumulativeFreq[base] +
                            BWTOccValue(index_, range.End + 1, base);
                if (range.End < range.Start) range.Start = 0;
            }

            if (range.Start != 0) {
                if (range.Level == readLength_) return true;
                ++range.Level;
                continue;
            }
            return false;
        }
    }
}

bool ExactSearch::search(const unsigned char* encodedRead, ExactResult& result) const {
    if (!encodedRead || lookupSize_ == 0) return false;
    if (readLength_ <= lookupSize_) {
        // The seed table is keyed on lookupSize_ bases: match a shorter read whole.
        const SeedRange first = initialBaseRange(encodedRead[0]);
        unsigned long start = first.Start, end = first.End;
        for (int i = 1; i < readLength_ && start; ++i) extend(start, end, encodedRead[i]);
        result.Start = start;
        result.End = end;
        result.Skip = 0;
        return start != 0;
    }
    unsigned key = 0;
    for (int i = 0; i < lookupSize_; ++i) {
        key |= static_cast<unsigned>(encodedRead[i]) << (2 * i);
    }

    SearchRange range{};
    range.Start = seedStart_[key];
    range.End = seedEnd_[key];
    range.Level = lookupSize_ + 1;
    range.Skip = 0;
    const bool found = searchExact(encodedRead, range);
    result.Start = range.Start;
    result.End = range.End;
    result.Skip = range.Skip;
    return found;
}

}  // namespace batmis
