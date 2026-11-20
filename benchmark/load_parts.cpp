// load_parts: where batmis's index load time goes, and what mmap saves (mmap_loading.md).
//   load_parts <prefix>    (prefix.bwt/.fmv/.sa, rev<prefix>..., prefix.pac)
// Times, per index: BWTLoad without the SA; with it (fread, as the mappers do); mapping the
// .sa instead (MAP_PRIVATE, then touching every page, or only mapping); and reading the .pac.
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

extern "C" {
#include "MemManager.h"
#include "BWT.h"
}

using Clock = std::chrono::steady_clock;
static double Since(Clock::time_point t) { return std::chrono::duration<double>(Clock::now() - t).count(); }

// Maps a .sa file; saValue points at its 7th word, which the loader would overwrite with -1.
static unsigned* MapSa(const std::string& name, size_t& bytes, int flags) {
    int fd = open(name.c_str(), O_RDONLY);
    struct stat st;
    fstat(fd, &st);
    bytes = st.st_size;
    void* p = mmap(nullptr, bytes, PROT_READ | PROT_WRITE, MAP_PRIVATE | flags, fd, 0);
    close(fd);
    if (p == MAP_FAILED) { std::perror("mmap"); std::exit(1); }
    unsigned* sa = reinterpret_cast<unsigned*>(static_cast<char*>(p) + 24);
    sa[0] = static_cast<unsigned>(-1);
    return sa;
}

int main(int argc, char** argv) {
    const std::string prefix = argv[1];
    const std::string dir = prefix.substr(0, prefix.find_last_of('/') + 1), base = prefix.substr(prefix.find_last_of('/') + 1);
    MMMasterInitialize(3, 0, FALSE, nullptr);
    MMPool* pool = MMPoolCreate(524288);
    for (const std::string& p : {dir + "rev" + base, prefix}) {
        auto t = Clock::now();
        BWT* a = BWTLoad(pool, (p + ".bwt").c_str(), (p + ".fmv").c_str(), nullptr, nullptr, nullptr, nullptr);
        const double noSa = Since(t);
        t = Clock::now();
        BWT* b = BWTLoad(pool, (p + ".bwt").c_str(), (p + ".fmv").c_str(), (p + ".sa").c_str(), nullptr, nullptr, nullptr);
        const double withSa = Since(t);
        size_t bytes;
        t = Clock::now();
        unsigned* m = MapSa(p + ".sa", bytes, 0);
        const double mapOnly = Since(t);
        t = Clock::now();
        unsigned long sum = 0;
        for (size_t i = 0; i < bytes / 4 - 8; i += 1024) sum += m[i];  // one word per page
        const double touch = Since(t);
        t = Clock::now();
        unsigned* q = MapSa(p + ".sa", bytes, MAP_POPULATE);
        const double populate = Since(t);
        std::printf("%s: BWTLoad without SA %.3f s; with SA (fread) %.3f s, so SA %.3f s | mmap %.4f s + touching every page "
                    "%.3f s; MAP_POPULATE %.3f s (%zu MB, check %lu %u %u)\n",
                    base.c_str(), noSa, withSa, withSa - noSa, mapOnly, touch, populate, bytes >> 20, sum % 7, b->saValue[5], q[5]);
        (void)a;
    }
    auto t = Clock::now();
    std::ifstream in((prefix + ".pac").c_str(), std::ios::binary | std::ios::ate);
    const std::streamsize size = in.tellg();
    in.seekg(0);
    std::vector<unsigned char> pac(size + 16, 0);
    in.read(reinterpret_cast<char*>(pac.data()), size);
    std::printf(".pac read into a vector: %.3f s (%ld MB)\n", Since(t), size >> 20);
    return 0;
}
