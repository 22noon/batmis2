# Branch_Detect: stored-count lookup versus direct counting (BRANCHTHRESHOLD)

## Summary

`Get_SARange` (`src/batman_zero_min.cpp`) narrows an SA range `[Start, End]`
to the child for one base. `Branch_Detect` finds all four child ranges at once.
Each can work in one of two ways:

- **lookup**: stored-count (rank) lookups at both `Start` and `End+1`.
  `Branch_Detect` uses `BWTAllOccValueTwoIndex`; `Get_SARange` uses two
  `BWTOccValue` calls, or `BWTOccValueTwoIndex` in the paired-rank build. The
  cost does not depend on the range width.
- **count**: one lookup at `Start` plus `Count_Range(Start, End)`, a direct
  bit-plane count of the BWT characters in the range. The cost grows with the
  width.

`BRANCHTHRESHOLD` in `src/batman_zero_min.cpp` picks between them, in both
functions. Ranges with
`End-Start < BRANCHTHRESHOLD` are counted; wider ranges use the lookup. The
default is **128**. Counting is faster than the lookup at every width up to
about 256 rows when the index data is not cached, and up to about 400 rows when
it is. Beyond that it falls behind quickly, and is about 16× slower at 65,536
rows. 128 keeps about a 2× margin below the cold break-even point.

On the random 500 Mb test genome the setting makes no measurable difference to
whole-program time (within ±3%; see below). Most `Branch_Detect` calls there
miss the cache, and cold counting is only 2–5% faster per call. The threshold
matters for cached data and for repetitive genomes, where wide ranges persist
deeper into the search. Those wide ranges must keep using the lookup, which is
why the threshold should not be raised above about 256.

Override at build time with `make DEFS=-DBATMAN_BRANCHTHRESHOLD=N`. `N=0`
disables counting entirely.

## Method

The benchmark (`branch_threshold_rank.cpp`) loads the 500 Mb random-genome
index and times one call of each method per SA range. It reports the best of 3
runs, in ns per call, pinned to one core (`taskset -c 2`) on an Intel Core
Ultra 7 165H, built with `-O3 -funroll-loops -mpopcnt`.

- **count** uses the production `Count_Range`. This is the scalar bit-plane
  version: low bit of each base pair = C or T, high bit = G or T, both = T,
  and A comes from the length. It handles whole words 64 bits at a time.
- **avx2** uses the same count with an AVX2 middle loop for ranges of 16 or
  more words. It is included for comparison and is not in production.
- **Cold** means random ranges over the whole index, so nearly every call
  misses the cache. **Hot** means 256 ranges reused, so the data is in L1/L2.
- **Dependent** means each call waits for the previous result, as in a search
  that walks down the tree. **Independent** lets the CPU overlap calls.

Before timing, 1M random ranges are checked, including ranges that cross the
`$` row (`inverseSa0`). All three methods agree on all of them.

## Results (ns per call)

| rows | cold dep. lookup | count | avx2 | hot dep. lookup | count | avx2 |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 163.6 | 158.4 | 158.7 | 31.4 | 21.5 | 21.2 |
| 4 | 161.7 | 160.8 | 158.8 | 30.8 | 21.4 | 21.3 |
| 16 | 165.3 | 161.2 | 159.9 | 29.5 | 21.0 | 21.7 |
| 32 | 162.6 | 159.9 | 159.2 | 29.5 | 21.2 | 21.8 |
| 64 | 169.3 | 162.4 | 158.9 | 29.9 | 21.8 | 23.2 |
| 128 | 168.1 | 164.7 | 161.2 | 29.8 | 23.4 | 24.9 |
| 256 | 172.9 | 172.6 | 172.3 | 30.1 | 27.0 | 32.8 |
| 384 | 174.9 | 179.5 | 179.8 | 31.4 | 30.4 | 35.1 |
| 512 | 173.9 | 185.1 | 186.1 | 32.8 | 34.2 | 37.8 |
| 1024 | 171.6 | 204.6 | 208.8 | 31.2 | 44.7 | 46.3 |
| 4096 | 173.4 | 351.9 | 353.7 | 29.6 | 113.0 | 101.7 |
| 65536 | 169.1 | 2804.1 | 2497.8 | 46.0 | 1476.2 | 1209.3 |

The break-even point is about 256 rows cold and about 400–500 rows hot.

## Why counting gains so little on cold data

`BWTAllOccValueTwoIndex` prefetches both positions before decoding either, so
its two cache misses overlap. In practice it costs about one miss latency, the
same as the single miss that counting still pays for `Start`. Counting therefore
saves only the decode work, about 8 ns per call. That saving is visible on hot
data (about 30% faster below 128 rows) but is lost in the ~160 ns miss on cold
data.

On cold data the AVX2 loop adds about 10–15% beyond 2,048 rows, where loading
the words from memory, not popcounting, sets the cost. Up to 256 rows it never
runs, because those ranges are under 16 words and fall back to the scalar code.

## Whole-program check

500 Mb genome, 2 runs each, mean search time:

| reads, options | T=0 (s) | T=128 (s) | change |
|---|---:|---:|---:|
| r32 `-n 1 -m 0` | 14.49 | 14.57 | +0.6% |
| r32 `-n 1 -m 1` | 9.42 | 9.17 | −2.7% |
| p32_mm1L `-n 1 -m 1` | 6.04 | 6.06 | +0.3% |
| p32_mm2 `-n 1 -m 0` | 4.03 | 4.12 | +2.3% |
| r24 `-n 1 -m 0` | 21.91 | 21.35 | −2.5% |

All differences are within run-to-run noise. On this genome every
`Branch_Detect` range for 32 bp reads is under 64 rows, and 99.2% are for
24 bp reads, so T=64 and T=128 take the same path almost everywhere.

Two build details matter:

- **`-mpopcnt` is required.** Without it, `__builtin_popcount` compiles to a
  libgcc call. With that build, T=64 was 1–9% slower than T=0. It is in the
  default `CXXFLAGS` of both `src/Makefile` and `index/Makefile`.
- **The older per-base `Count_Range`** (four compare-and-popcounts per word)
  was 1–8% slower at T=128 and should not be reinstated.

Output at T=128 is byte-identical to T=0. This was checked over 48 runs: three
read sets plus 12, 20 and 36 bp genome-end overhang sets, `-n 0` and
`-n 1 -m 0/1/5`, with standard and paired rank. Brute-force checks of the mutated and genome-end overhang reads
found no errors.
