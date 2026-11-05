# Structured single-index exact mapper

## Design and scope

The original `src/batman_zero_min.cpp` remains unchanged and is the comparison
baseline. The new `batman_zero_structured` executable keeps its command line,
input parsing, index selection, timing stages, and output format. Search logic
and fixed-length seed state live in `src/exact_search.cpp`, declared in
`src/exact_search.hpp`.

Responsibilities are intentionally narrow:

- `batman_zero_structured.cpp` owns the existing BatMIS memory-pool/index load
  and handles FASTA/FASTQ input, timing, and TSV output.
- `ExactSearch` is a non-owning view of the loaded `BWT`; it builds the existing
  3- or 6-base preindex once for the first read length and runs the specialized
  exact-search loop. Its two 4096-entry `std::array`s contain inclusive SA-row
  bounds. Search has no per-extension allocation.
- `ExactResult` carries the SA interval and the sampled-singleton `Skip` state
  used by the minimized mapper. Output continues to report SA rows, not
  reference coordinates. The current mapper does not load an SA or recover
  coordinates; coordinate recovery remains a separate follow-up component if
  output requirements expand.

The extension operation remains BatMIS's native backward extension:

```text
start = cumulativeFreq[base] + BWTOccValue(start, base) + 1
end   = cumulativeFreq[base] + BWTOccValue(end + 1, base)
```

The native BWT uses inclusive SA-row bounds. `BWTOccValue()` adjusts ranks for
the omitted sentinel row at `inverseSa0`; the search code retains its existing
packed-symbol read and SA-sampling behavior (`saInterval` 8). No BWA index,
rank implementation, or index builder was introduced or changed.

Seed behavior is preserved: the initial seed is 3 bases for reads shorter than
28 and 6 bases otherwise; subsequent bases are extended from the existing
preindex interval. The parser still fixes the read length from the first
record, reports later length mismatches as misses, and rejects non-ACGT
sequences. The known length 1–3 behavior is preserved: those cases are compared
with the baseline, but are not claimed as a correct general exact-search mode.

## Files added or changed

- `src/exact_search.hpp`, `src/exact_search.cpp`: fixed-size seed workspace,
  lightweight interval/result types, and exact-search algorithm.
- `src/batman_zero_structured.cpp`: structured executable/driver; baseline
  `src/batman_zero_min.cpp` is unchanged.
- `src/Makefile`: explicit `batman_zero_structured` and `test-structured`
  targets; default `all` still builds the original minimized executable.
- `tests/test_structured_exact_search.sh` and
  `tests/test_exact_search_intervals.cpp`: generated-index CLI and exact
  interval differential/oracle regressions.
- `benchmark/run_structured_exact_benchmark.py` and
  `benchmark/structured_exact_search_trials.jsonl`: repeatable harness and raw
  trial data.

## Correctness checks

The regression script creates a 2,000-base single-contig reference, builds its
BatMIS forward/reverse indexes, and compares the baseline and structured TSV
byte-for-byte for read lengths 1–12. The query batches include every possible
ACGT string for lengths 1–6, reference substrings (including boundary
substrings), repeated patterns, and absent patterns.

For lengths 4–12, the structured hit counts also match a naïve substring-count
oracle for every generated query. For lengths 1–3, only byte-identical behavior
to the baseline is asserted, preserving the documented limitation.

Command and result:

```sh
tests/test_structured_exact_search.sh \
  /tmp/batmis_structured_build_matched \
  /tmp/batmis_structured_build_matched/bwt_build_index \
  /tmp/batmis_structured_build_matched/test_exact_search_intervals
```

Result: passed. The batches contained 5,592 queries total. The C++ differential
test compared the found flag, inclusive interval start/end (including the
baseline's internal empty insertion bounds), and sampled `Skip` state for
every query. The CLI outputs were byte-identical as well.
Additionally, 5,508 queries at lengths 4–12 were checked against the naïve
oracle. The 1–3 batches matched the baseline only; they are not claimed as
correct general exact-search behavior. The test uses the baseline's fixed
zero-padded query buffer to preserve that short-read behavior.

The same matched-build binaries were also run on the complete 100,000-read
benchmark input with output enabled (outside timed trials). Their 100,000-row
TSV files were byte-identical, SHA256
`4412d1193a1aebc82e7266b2b101580f67dcb475ad3016ddd438868a6fdf6e48`.

The BatMIS source tree had no pre-existing registered CTest suite. The existing
CTest suite in the sibling `minimal_bwt` workspace was also run:
`ctest --test-dir /home/chandana.tennakoon@pirbright.ac.uk/git/batmis/minimal_bwt/build
--output-on-failure` — 11/11 passed. Those CTests exercise the sibling BWA-backed
implementation and are not a substitute for the BatMIS differential test.

## Benchmark method

The old and new executable were built from the current sources using the same
GCC/G++ 13.3.0 toolchain and the flags used by the preceding controlled index
comparison:

```text
-w -O3 -funroll-loops -maccumulate-outgoing-args -msse2 -g -O2
```

The final `-O2` is the effective optimization level. The seven native index C
objects and both executables were built in `/tmp/batmis_structured_build_matched`
so no generated objects or archives were written into the source tree. Both
variants linked the same `libfmindex_core.a`.

Workload and files:

- 100,000 100-base reads from `/tmp/zero_min_exp1/reads_100k.fq`.
- BatMIS reverse BWT/occurrence files at
  `/tmp/zero_min_exp1/revref.fa.bwt` and `.fmv`, queried via genome prefix
  `/tmp/zero_min_exp1/ref.fa`.
- `-m 0 --quiet --timing`: exact only, count all hits, no SA load, location, or
  output formatting.
- Ten paired trials; execution order alternated baseline-first and
  structured-first. Both variants reported 100,000 records and 75,000 hits in
  every run.
- Host: Linux 6.17, Intel Core Ultra 7 165H, 22 logical CPUs. Load average was
  2.10–2.20 during the final trials.

Reproduction commands (run from the repository root):

```sh
BUILD=/tmp/batmis_structured_build_matched
mkdir -p "$BUILD"
CFLAGS='-w -O3 -funroll-loops -maccumulate-outgoing-args -msse2 -g -O2'
for source in BWT.c BWTConstruct.c MiscUtilities.c MemManager.c TextConverter.c QSufSort.c DNACount.c; do
  gcc $CFLAGS -I index -c "index/$source" -o "$BUILD/${source%.c}.o"
done
ar rcs "$BUILD/libfmindex_core.a" "$BUILD"/*.o
g++ -std=c++17 -fpermissive $CFLAGS -I index src/batman_zero_min.cpp "$BUILD/libfmindex_core.a" -lm -o "$BUILD/batman_zero_min"
g++ -std=c++17 -fpermissive $CFLAGS -I index src/batman_zero_structured.cpp src/exact_search.cpp "$BUILD/libfmindex_core.a" -lm -o "$BUILD/batman_zero_structured"
g++ -std=c++17 -fpermissive $CFLAGS -I src -I index tests/test_exact_search_intervals.cpp src/exact_search.cpp "$BUILD/libfmindex_core.a" -lm -o "$BUILD/test_exact_search_intervals"
g++ -std=c++17 $CFLAGS -I index index/bwt_build_index.cpp "$BUILD/libfmindex_core.a" -lm -o "$BUILD/bwt_build_index"
tests/test_structured_exact_search.sh "$BUILD" "$BUILD/bwt_build_index" "$BUILD/test_exact_search_intervals"
python3 benchmark/run_structured_exact_benchmark.py \
  --baseline "$BUILD/batman_zero_min" \
  --structured "$BUILD/batman_zero_structured" \
  --reads /tmp/zero_min_exp1/reads_100k.fq \
  --genome-prefix /tmp/zero_min_exp1/ref.fa \
  --trials 10 --output benchmark/structured_exact_search_trials.jsonl
```

The raw trial records contain exact command lines, binary/input/index hashes,
per-stage times, external wall times, user CPU times, and load samples. Binary
SHA256 values for this run:

| Executable | SHA256 |
|---|---|
| baseline `batman_zero_min` | `85f3120b2c6b993216c8be8502e0b585135ff3601d051eedb377811a6a8dcabb` |
| `batman_zero_structured` | `080f67de2b0c619fe8a702a079d05cb0c467bcef41be46664b2f02dd48113a83` |

Workload identities: reads `a6adffaf0ca4718a0a00b3d1b8378503b39d31130608d1a684a030b9a347f4fd`,
reference `820adc6c0960567688d20e8cc4c775c4c6cfedfab893755d76e50938db940b6a`,
reverse BWT `140087efc5f2260ded757e61c64a1222edf368aae19b73bbd01d79c40eb397ad`,
reverse FMV `4cca722cb7dea4d6b8a929e22e552a40d25233819a0e2ca556b3c62614e67533`.
The raw JSONL SHA256 is
`1dcc677d644adf837e29b65531a3f194bc04b85166cb1f84d37de8ca49c868dc`.

## Results

Times are seconds. Paired differences are structured minus baseline; negative
values favor the structured version. The percentage is the median of the ten
per-pair percentage changes, not the ratio of independent medians.

| Stage | Baseline median (range) | Structured median (range) | Paired median difference | Median paired change | Pairs favoring structured |
|---|---:|---:|---:|---:|---:|
| Search | 0.142497 (0.129804–0.148624) | 0.135419 (0.120017–0.152121) | -0.008295 | -5.87% | 7/10 |
| Internal total | 0.160237 (0.145818–0.167261) | 0.153116 (0.135833–0.171215) | -0.008064 | -5.07% | 7/10 |
| External wall | 0.163069 (0.147926–0.170436) | 0.155703 (0.138088–0.174041) | -0.008184 | -4.99% | 7/10 |
| User CPU | 0.159197 (0.142662–0.168031) | 0.152912 (0.131673–0.170582) | -0.008262 | -5.18% | 7/10 |

Input parsing medians were 0.012657 s baseline and 0.012634 s structured.
Preindex construction medians were 0.000257 s and 0.000229 s, respectively;
index loading medians were 0.000684 s and 0.000461 s. Output formatting was
zero in both modes. Search time starts before base encoding and ends after the
interval and hit-count calculation; parsing, index loading, seed-table setup,
and output are separately timed. Search-time paired differences (structured
minus baseline, seconds) were `-0.015490, -0.002241, +0.004186, +0.006180,
-0.023402, -0.012400, +0.013511, -0.014180, -0.017368, -0.004189`.

The structured version was faster in seven of ten pairs, with median paired
search reduction of 5.87%; three pairs favored baseline and the timing ranges
overlap substantially. This is suggestive evidence for the tested workload,
not a general performance guarantee or statistical-significance claim. The
cause of the observed timing difference is unresolved: the search state and
compilation-unit boundaries changed together, and no profiling/disassembly
experiment isolated them. The rank primitive and interval arithmetic are
unchanged.

## Recommendation and limitations

Keep `batman_zero_min` as the default and retain `batman_zero_structured` as a
separately selectable implementation while gathering broader workload
evidence. The exact interval results agree with the minimized baseline, and
the measured 100k workload shows a modest, mostly consistent speed reduction.
Do not attribute that difference to the BatMIS rank/index implementation:
both versions use the same index, rank calls, seed strategy, read order, and
matching semantics.

The structured executable intentionally preserves the baseline's interval
reporting rather than loading sampled SA values and emitting coordinates. The
known read-length 1–3 limitation is unchanged. The next useful step is to
exercise the same binary pair on additional read lengths and datasets before
considering a default change or adding a separate optional location component.
