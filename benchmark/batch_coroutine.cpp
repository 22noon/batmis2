// Batched exact search written as C++20 coroutines, against the hand-written lanes
// of Batch_Exact (src/batman_zero_min.cpp). Each read's search is ordinary sequential
// code that prefetches what its next step reads and suspends; a scheduler resumes B
// of them round-robin. If this costs about the same as the lanes, the one-mismatch
// search (much more control flow) can be batched the same way without rewriting it
// as a state machine. Results: benchmark/batched_reads.md.
// Build from src/ (pulls in batman_zero_min.cpp):
//   g++ -std=c++20 -fpermissive -w -I. -I../index -O3 -funroll-loops -mpopcnt \
//       ../benchmark/batch_coroutine.cpp ../index/libfmindex_core.a -lm -o batch_coroutine
// Run:  taskset -c 2 ./batch_coroutine <genome prefix> <reads.fq>... [--no-ref]
#define main orig_main
#include "batman_zero_min.cpp"
#undef main
#include <chrono>
#include <coroutine>

namespace {

using Clock = std::chrono::steady_clock;
double Seconds(Clock::time_point t) { return std::chrono::duration<double>(Clock::now() - t).count(); }

// Coroutine frames come from a free list, so starting a read allocates nothing.
struct Frame_Pool {
    std::vector<void*> free;
    std::size_t size = 0;
    void* Get(std::size_t n) {
        if (n <= size && !free.empty()) {
            void* p = free.back();
            free.pop_back();
            return p;
        }
        size = std::max(size, n);
        return ::operator new(size);
    }
    void Put(void* p) { free.push_back(p); }
} pool;

struct Task {
    struct promise_type {
        Task get_return_object() { return Task{std::coroutine_handle<promise_type>::from_promise(*this)}; }
        std::suspend_always initial_suspend() noexcept { return {}; }
        std::suspend_always final_suspend() noexcept { return {}; }
        void return_void() {}
        void unhandled_exception() { std::terminate(); }
        static void* operator new(std::size_t n) { return pool.Get(n); }
        static void operator delete(void* p) { pool.Put(p); }
    };
    std::coroutine_handle<promise_type> h;
};

// co_await Yield{}: the caller has issued its prefetches; let the other reads run.
struct Yield {
    bool await_ready() const noexcept { return false; }
    void await_suspend(std::coroutine_handle<>) const noexcept {}
    void await_resume() const noexcept {}
};

struct Result {
    SARange R;
    bool Found;
};

// Search_Right + Search_Right_Singleton from R (Start_Range), as Batch_Exact's lanes.
Task Search(std::string_view read, const uint64_t* packed, SARange R, Result& out) {
    const int len = static_cast<int>(read.size());
    out.Found = false;
    while (!(R.End == R.Start || R.Skip)) {
        Prefetch_Row(revfmi, R.Start);
        Prefetch_Row(revfmi, R.End + 1);
        co_await Yield{};
        Get_SARange(read[R.Level - 1], R, revfmi);
        if (!R.Start) { out.R = R; co_return; }
        if (R.Level == len) { out.R = R; out.Found = true; co_return; }
        R.Level++;
    }
    int done = R.Level - 1;  // bases matched when a sampled row is reached
    for (;;) {
        if (!R.Skip && R.Start % SAINTERVAL == 0) {
            if (Use_Ref) {
                __builtin_prefetch(revfmi->saValue + R.Start / SAINTERVAL);
                co_await Yield{};
                if (!Locate(revfmi, R, true, 0, len, done)) { out.R = R; co_return; }
                __builtin_prefetch(Pac.data() + (R.Pos + done) / 4);
                __builtin_prefetch(Pac.data() + (R.Pos + len) / 4);
                co_await Yield{};
                if (Rest_Mismatches(R, true, 0, len, done, packed) == 0) { R.Level = len; out.Found = true; }
                out.R = R;
                co_return;
            }
            R.Skip++;
            R.End = R.Start;
        }
        if (done == len) { R.End = R.Start; out.R = R; out.Found = true; co_return; }  // matched the last base
        Prefetch_Row(revfmi, R.Start);
        co_await Yield{};
        if (R.Start == revfmi->inverseSa0) { out.R = R; co_return; }
        const unsigned long index = R.Start - (R.Start > revfmi->inverseSa0);
        const unsigned now = revfmi->bwtCode[index / 16] << ((index % 16) * 2) >> (BITS_IN_WORD - 2);
        const unsigned expected = static_cast<unsigned char>(read[R.Level - 1]);
        const unsigned base = Read_Base(expected);
        if (expected != now) { out.R = R; co_return; }
        R.Start = revfmi->cumulativeFreq[base] + BWTOccValue(revfmi, R.Start, base) + 1;
        done = R.Level;
        if (R.Level < len) R.Level++;
    }
}

struct Slot {
    Task t;
    unsigned index;
    uint64_t* packed;  // the read's packed copy; stays put when slots are swapped, as the coroutine points at it
};

void Batch_Coroutine(const std::vector<std::string>& reads, unsigned B, std::vector<Result>& results) {
    std::vector<std::array<uint64_t, MAXTAG / 32 + 1>> buffers(B);
    std::vector<Slot> slots(B);
    for (unsigned i = 0; i < B; i++) slots[i].packed = buffers[i].data();
    unsigned active = 0, next = 0;
    auto start = [&](Slot& s) {
        while (next < reads.size()) {
            const unsigned i = next++;
            const std::string_view r(reads[i]);
            SARange R = Start_Range(r);
            if (R.Level > static_cast<int>(r.size()) || !R.Start) {
                results[i] = {R, R.Start != 0};
                continue;
            }
            if (Use_Ref) Pack_Read(r, s.packed);
            s.index = i;
            s.t = Search(r, s.packed, R, results[i]);
            return true;
        }
        return false;
    };
    for (unsigned i = 0; i < B && start(slots[active]); i++) active++;
    while (active) {
        for (unsigned i = 0; i < active;) {
            Slot& s = slots[i];
            s.t.h.resume();
            if (!s.t.h.done()) { i++; continue; }
            s.t.h.destroy();
            if (start(s)) { i++; continue; }
            std::swap(s, slots[--active]);
        }
    }
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 3) {
        std::fprintf(stderr, "Usage: %s <genome prefix> <reads.fq>... [--no-ref]\n", argv[0]);
        return 2;
    }
    bool useRef = true;
    std::vector<const char*> files;
    for (int a = 2; a < argc; a++) {
        if (std::strcmp(argv[a], "--no-ref") == 0) useRef = false;
        else files.push_back(argv[a]);
    }
    std::string pre = argv[1], rb, ro, fb, fo, rs, fs;
    BuildIndexPaths(pre, rb, ro, fb, fo, &rs, &fs);
    revfmi = initFMI(rb.c_str(), ro.c_str(), useRef ? rs.c_str() : nullptr);
    Use_Ref = useRef && LoadPac(pre + ".pac");
    Text_Length = revfmi->textLength;
    std::printf("%s; ns per read, best of 3 (search only)\n", Use_Ref ? "reference on" : "--no-ref");
    std::printf("%-14s %8s | %8s %8s %8s %8s | %8s %8s %8s %8s | %s\n", "reads", "n", "lane B=1", "lane 4", "lane 8",
                "lane 16", "coro 1", "coro 4", "coro 8", "coro 16", "check");
    for (const char* f : files) {
        std::ifstream in(f);
        std::string id, seq;
        std::vector<std::string> reads;
        std::string codes;
        while (ReadNextRecord(in, id, seq)) {
            std::string r(seq.size(), 0);
            bool ok = true;
            for (std::size_t i = 0; i < seq.size() && ok; i++) ok = (r[i] = static_cast<char>(EncodeBase(seq[i]))) >= 0;
            if (!ok || (!reads.empty() && r.size() != reads[0].size())) continue;
            codes += r;
            reads.push_back(std::move(r));
        }
        STRINGLENGTH = static_cast<int>(reads[0].size());
        Build_Tables(revfmi);
        const unsigned n = reads.size();
        const unsigned Bs[4] = {1, 4, 8, 16};
        double lane[4], coro[4];
        unsigned long bad = 0;
        std::vector<SARange> ranges(n);
        std::vector<unsigned char> found(n);
        std::vector<Result> res(n);
        for (int k = 0; k < 4; k++) {
            lane[k] = coro[k] = 1e99;
            for (int rep = 0; rep < 3; rep++) {
                auto t = Clock::now();
                Batch_Exact(codes.data(), n, STRINGLENGTH, Bs[k], ranges, found);
                lane[k] = std::min(lane[k], Seconds(t));
                t = Clock::now();
                Batch_Coroutine(reads, Bs[k], res);
                coro[k] = std::min(coro[k], Seconds(t));
            }
            for (unsigned i = 0; i < n; i++) {
                if (found[i] != res[i].Found) { bad++; continue; }
                if (!found[i]) continue;
                const SARange &a = ranges[i], &b = res[i].R;
                bad += Use_Ref && a.Skip ? !(b.Skip && a.Pos == b.Pos)
                                         : !(a.Start == b.Start && a.Skip == b.Skip && (a.Skip || a.End == b.End));
            }
        }
        const char* slash = std::strrchr(f, '/');
        const double ns = 1e9 / n;
        std::printf("%-14s %8u | %8.0f %8.0f %8.0f %8.0f | %8.0f %8.0f %8.0f %8.0f | %s\n", slash ? slash + 1 : f, n,
                    lane[0] * ns, lane[1] * ns, lane[2] * ns, lane[3] * ns, coro[0] * ns, coro[1] * ns, coro[2] * ns,
                    coro[3] * ns, bad ? "MISMATCH" : "same results");
        std::fflush(stdout);
    }
    return 0;
}
