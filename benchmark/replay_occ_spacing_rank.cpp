#include <chrono>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <vector>

extern "C" {
#include "BWT.h"
#include "MemManager.h"
}

struct TraceEntry {
    std::uint32_t row;
    std::uint8_t base;
    std::uint8_t origin;
    std::uint16_t reserved;
};
static_assert(sizeof(TraceEntry) == 8, "rank trace entry layout changed");

int main(int argc, char** argv) {
    if (argc != 4) return 2;
    MMMasterInitialize(3, 0, FALSE, nullptr);
    MMPool* pool = MMPoolCreate(524288);
    BWT* bwt = BWTLoad(pool, argv[1], argv[2], nullptr, nullptr, nullptr, nullptr);

    std::ifstream input(argv[3], std::ios::binary | std::ios::ate);
    const std::streamoff bytes = input.tellg();
    if (bytes < 0 || bytes % sizeof(TraceEntry)) return 3;
    input.seekg(0);
    std::vector<TraceEntry> trace(static_cast<std::size_t>(bytes) / sizeof(TraceEntry));
    input.read(reinterpret_cast<char*>(trace.data()), bytes);
    if (!input) return 4;

    volatile std::uint64_t checksum = 0;
    const auto begin = std::chrono::steady_clock::now();
    for (const TraceEntry& query : trace) checksum += BWTOccValue(bwt, query.row, query.base);
    const double seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - begin).count();
    std::printf("replay: queries=%zu checksum=%llu seconds=%.9f occ_minor_bytes=%u occ_major_bytes=%u\n",
        trace.size(), static_cast<unsigned long long>(checksum), seconds,
        bwt->occSizeInWord * 4, bwt->occMajorSizeInWord * 4);
    return 0;
}
