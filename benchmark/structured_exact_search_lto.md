# LTO experiment: structured BatMIS exact search

## Result

LTO removed the `ExactSearch::search()` and `ExactSearch::searchExact()` call
boundaries, but it did **not** improve this workload. Across ten alternating
three-configuration blocks, LTO was slower than the otherwise identical
structured build in 8/10 paired search-time comparisons. The median paired
search difference was **+8.596 ms (+5.41%)**, where positive means LTO took
longer. Keep the existing non-LTO structured build and default build
unchanged.

This is evidence for the tested 100,000-read workload, not a claim about other
read lengths or machines. The index/rank implementation was not built with
LTO; the experiment isolates LTO to the two structured C++ translation units.

## Build configuration

Environment: GCC/G++ 13.3.0 (`Ubuntu 13.3.0-6ubuntu2~24.04`), Linux
6.17.0-19-generic, x86-64. The matched build flags were:

```text
-w -O3 -funroll-loops -maccumulate-outgoing-args -msse2 -g -O2
```

The final `-O2` is the effective optimization level. Both structured variants
compiled `src/batman_zero_structured.cpp` and `src/exact_search.cpp` with the
same flags; the LTO build added only `-flto` to both compilation and final
link. Both linked the same prebuilt `index/libfmindex_core.a` (SHA256
`2fea72bad1118ac6178f607e5e9cebcd21cd8e61e0daf68477126b7cd6a8c574`). The
minimized executable was the previously built matched-flags reference.

The source tree's existing `src/Makefile` and default targets were not changed
for this experiment. All binaries and intermediate files were placed under
`/tmp/batmis_structured_lto_experiment/`.

Build commands, from the repository root:

```sh
BASE=/tmp/batmis_structured_build_matched
BUILD=/tmp/batmis_structured_lto_experiment
mkdir -p "$BUILD/no_lto" "$BUILD/lto"
CXXFLAGS='-std=c++17 -fpermissive -w -O3 -funroll-loops -maccumulate-outgoing-args -msse2 -g -O2'
CORE="$BASE/libfmindex_core.a"
g++ $CXXFLAGS -I index src/batman_zero_structured.cpp src/exact_search.cpp \
  "$CORE" -lm -o "$BUILD/no_lto/batman_zero_structured"
g++ $CXXFLAGS -flto -I index src/batman_zero_structured.cpp src/exact_search.cpp \
  "$CORE" -lm -o "$BUILD/lto/batman_zero_structured"
g++ $CXXFLAGS -I src -I index tests/test_exact_search_intervals.cpp \
  src/exact_search.cpp "$CORE" -lm -o "$BUILD/no_lto/test_exact_search_intervals"
g++ $CXXFLAGS -flto -I src -I index tests/test_exact_search_intervals.cpp \
  src/exact_search.cpp "$CORE" -lm -o "$BUILD/lto/test_exact_search_intervals"
```

`$BUILD/{no_lto,lto}/batman_zero_min` and `bwt_build_index` were symlinks to
the matched reference binaries. Hashes:

| Binary | SHA256 |
|---|---|
| Minimized A | `85f3120b2c6b993216c8be8502e0b585135ff3601d051eedb377811a6a8dcabb` |
| Structured B, no LTO | `080f67de2b0c619fe8a702a079d05cb0c467bcef41be46664b2f02dd48113a83` |
| Structured C, LTO | `602724df664f86fd082d6d0b572293b327c4b071dbe20b1b0f3fa965e95a4122` |

The no-LTO hash exactly matches the prior structured benchmark binary. The
index library and data hashes are recorded in the raw JSONL file.

## Generated-code inspection

Commands:

```sh
nm -C "$BUILD/no_lto/batman_zero_structured" | grep 'ExactSearch::'
nm -C "$BUILD/lto/batman_zero_structured" | grep 'ExactSearch::'
objdump -d -C "$BUILD/no_lto/batman_zero_structured"
objdump -d -C "$BUILD/lto/batman_zero_structured"
```

In B, `nm` shows separate `ExactSearch::search()` and
`ExactSearch::searchExact()` symbols. Disassembly shows `main` calling
`ExactSearch::search()`, which in turn calls `searchExact()`. In C, neither
method remains as a symbol and the corresponding calls are absent; the seed
lookup and extension loop are emitted in `main`. The inlined path still calls
the C `BWTOccValue()` rank function twice for the general multi-row extension;
the index archive was not LTO-enabled, so rank implementation boundaries
remain. The LTO executable's text size was 158,567 bytes versus 160,251 bytes
for B. The smaller text and removed C++ calls did not translate to faster
search in these measurements. No source-level inlining changes were made.

## Correctness

Commands:

```sh
tests/test_structured_exact_search.sh \
  /tmp/batmis_structured_lto_experiment/no_lto \
  /tmp/batmis_structured_lto_experiment/no_lto/bwt_build_index \
  /tmp/batmis_structured_lto_experiment/no_lto/test_exact_search_intervals
tests/test_structured_exact_search.sh \
  /tmp/batmis_structured_lto_experiment/lto \
  /tmp/batmis_structured_lto_experiment/lto/bwt_build_index \
  /tmp/batmis_structured_lto_experiment/lto/test_exact_search_intervals
```

Both passed. Each ran the existing generated-index differential suite over
5,592 queries (lengths 1–12), comparing exact interval results against the
minimized baseline; the helper checks interval bounds and `Skip`. Oracle
checks for lengths 4–12 also passed. The established limitation for read
lengths 1–3 remains preserved; these lengths are compared with the baseline,
not asserted as general exact-search correctness.

The 100,000-read TSV outputs from A, B, and C were byte-identical, each with
SHA256 `4412d1193a1aebc82e7266b2b101580f67dcb475ad3016ddd438868a6fdf6e48`.
All timed variants reported 100,000 records and 75,000 hits on every run.

## Benchmark method and data

Workload: `/tmp/zero_min_exp1/reads_100k.fq`, 100,000 reads of 100 bases;
genome prefix `/tmp/zero_min_exp1/ref.fa`, using the same reverse BatMIS index
files as the prior matched benchmark. Command per run:

```sh
BINARY -q /tmp/zero_min_exp1/reads_100k.fq \
  -g /tmp/zero_min_exp1/ref.fa -m 0 --quiet --timing
```

No SA load, output, or location. Ten blocks ran all three variants, rotating
through all six execution-order permutations. The runner is
`benchmark/run_structured_lto_benchmark.py`; it records external wall time,
child user CPU, application stage timings, hit counts, binary/data hashes and
load averages. The raw 30-run records are in
`benchmark/structured_exact_search_lto_trials.jsonl` (SHA256
`da84d801e4bdbee2e82e75e98076e771ac437ee4d23443460b61b6bafb02c845`).
Load average's 1-minute value ranged from 2.91 to 2.99 over the experiment.

### Per-configuration summary

Times are seconds; ranges are min–max across ten trials.

| Metric | A: minimized median (range) | B: structured median (range) | C: structured LTO median (range) |
|---|---:|---:|---:|
| Search | 0.162129 (0.146209–0.181100) | 0.155687 (0.136258–0.171418) | 0.165097 (0.146786–0.182148) |
| Internal total | 0.181868 (0.163943–0.203295) | 0.175681 (0.153706–0.193122) | 0.184929 (0.164954–0.204098) |
| External wall | 0.184378 (0.166540–0.206273) | 0.178104 (0.155928–0.195596) | 0.187426 (0.167578–0.206556) |
| User CPU | 0.179903 (0.160083–0.200856) | 0.173794 (0.148458–0.191612) | 0.183413 (0.165108–0.202181) |

### Paired B-to-C differences

Differences are `C - B`; positive values mean LTO was slower. The percentage
is the median of the ten per-block percentage differences using B as the
denominator.

| Metric | Median paired difference | Median paired change | Blocks favoring LTO |
|---|---:|---:|---:|
| Search | +0.008596 s | +5.41% | 2/10 |
| Internal total | +0.008458 s | +4.69% | 3/10 |
| External wall | +0.008642 s | +4.68% | 3/10 |
| User CPU | +0.009952 s | +5.69% | 3/10 |

Search paired differences by block, in seconds: `+0.017961, +0.004406,
-0.001429, +0.010730, +0.000751, +0.006824, -0.002260, +0.012661,
+0.010367, +0.010528`. Thus 8/10 blocks favored B. A versus B is included as
context: B was faster in 8/10 blocks, with a median paired search change of
−5.18%; this is consistent with the earlier structured-search result and is
not an LTO comparison.

## Conclusion

LTO successfully inlined the structured C++ search entry points across the
translation-unit boundary, but that did not help end-to-end search here. C
was slower than B in most paired blocks, with strongly overlapping ranges.
This does not establish why LTO regressed; the experiment shows that removing
these calls alone was insufficient to improve the measured workload. Keep LTO
opt-in and do not change the default build. A useful next experiment would be
to profile the B and C binaries at the same workload to see whether the inlined
loop shifts time into rank calls, branch behavior, or another component before
considering any build-default change.
