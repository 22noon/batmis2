#include <array>
#include <cstdio>
#include <vector>

extern "C" {
#include "BWT.h"
#include "MemManager.h"
}

int main(int argc, char** argv) {
    if (argc != 3) return 2;
    MMMasterInitialize(3, 0, FALSE, nullptr);
    MMPool* pool = MMPoolCreate(524288);
    BWT* bwt = BWTLoad(pool, argv[1], argv[2], nullptr, nullptr, nullptr, nullptr);

    std::vector<std::array<unsigned, 4>> prefix(bwt->textLength + 1);
    prefix[0] = {0, 0, 0, 0};
    for (unsigned row = 0; row < bwt->textLength; ++row) {
        prefix[row + 1] = prefix[row];
        const unsigned base = (bwt->bwtCode[row / 16] >> (30 - 2 * (row % 16))) & 3;
        ++prefix[row + 1][base];
    }

    unsigned long long checks = 0;
    for (unsigned row = 0; row <= bwt->textLength; ++row) {
        const unsigned adjusted = row - (row > bwt->inverseSa0);
        for (unsigned base = 0; base < 4; ++base) {
            const unsigned expected = prefix[adjusted][base];
            const unsigned observed = BWTOccValue(bwt, row, base);
            if (observed != expected) {
                std::fprintf(stderr,
                    "rank mismatch row=%u adjusted=%u primary=%u base=%u expected=%u observed=%u\n",
                    row, adjusted, bwt->inverseSa0, base, expected, observed);
                return 1;
            }
            ++checks;
        }
    }
    std::printf("rank oracle passed rows=%u symbol_queries=%llu primary=%u\n",
                bwt->textLength + 1, checks, bwt->inverseSa0);
    return 0;
}
