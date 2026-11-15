// Sa_Decode against BWTSaValue: every row of both indexes (or Rows random rows), decoded
// in ranges of mixed widths (single rows, narrow and wide ranges, one range of all rows),
// must give BWTSaValue's value.
// Usage: test_sa_decode <genome_prefix> [Rows (default 0 = all)]
// The mapper under test: batmis1.cpp, or batmis2.cpp with -DMAPPER_SOURCE='"batmis2.cpp"'.
#ifndef MAPPER_SOURCE
#define MAPPER_SOURCE "batmis1.cpp"
#endif
#define main mapper_main
#include MAPPER_SOURCE
#undef main

#include <random>

int main(int argc, char** argv) {
    if (argc < 2) {
        std::fprintf(stderr, "Usage: %s <genome_prefix> [rows]\n", argv[0]);
        return 2;
    }
    std::string revBwt, revOcc, fwBwt, fwOcc, revSa, fwSa;
    BuildIndexPaths(argv[1], revBwt, revOcc, fwBwt, fwOcc, &revSa, &fwSa);
    revfmi = initFMI(revBwt.c_str(), revOcc.c_str(), revSa.c_str());
    fwfmi = BWTLoad(mmPool, fwBwt.c_str(), fwOcc.c_str(), fwSa.c_str(), nullptr, nullptr, nullptr);
    const unsigned long Rows = argc > 2 ? std::strtoul(argv[2], nullptr, 10) : 0;
    std::mt19937_64 Rng(7);
    unsigned long Checked = 0, Failures = 0;
    for (const ActiveBWT* F : {static_cast<const ActiveBWT*>(revfmi), static_cast<const ActiveBWT*>(fwfmi)}) {
        const unsigned long N = F->textLength + 1;  // rows 0 .. textLength
        std::vector<Sa_Request> Requests;
        std::size_t Out = 0;
        auto Add = [&](unsigned long S, unsigned long E) {
            Requests.push_back({F, S, E, Out});
            Out += E - S + 1;
        };
        if (Rows == 0) {
            Add(0, N - 1);  // one range of everything
            for (unsigned long S = 0; S < N;) {  // and again in mixed widths
                const unsigned long W = std::min<unsigned long>(N - S, std::vector<unsigned long>{1, 2, 5, 7, 8, 9, 16, 40, 300}[Rng() % 9]);
                Add(S, S + W - 1);
                S += W;
            }
        } else {
            for (unsigned long k = 0; k < Rows; k++) {
                const unsigned long W = std::vector<unsigned long>{1, 1, 3, 8, 20, 200}[Rng() % 6];
                const unsigned long S = Rng() % (N - W + 1);
                Add(S, S + W - 1);
            }
        }
        std::vector<unsigned> Values(Out);
        Sa_Decode(Requests, Values);
        for (const Sa_Request& R : Requests)
            for (unsigned long Row = R.Start; Row <= R.End; Row++) {
                const unsigned Want = BWTSaValue(F, static_cast<unsigned>(Row)), Got = Values[R.Out + (Row - R.Start)];
                Checked++;
                if (Got != Want && Failures++ < 10)
                    std::fprintf(stderr, "FAIL %s row %lu: %u, BWTSaValue %u\n", F == revfmi ? "revfmi" : "fwfmi", Row, Got, Want);
            }
    }
    std::printf("checked=%lu failures=%lu\n", Checked, Failures);
    return Failures ? 1 : 0;
}
