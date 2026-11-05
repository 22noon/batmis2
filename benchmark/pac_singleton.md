# Finishing single-row paths from a sampled SA row (.sa + .pac)

> **Since the program refactor** (`five_mismatch.md`, "Program builder"): the hand-written
> per-read searches named below (`One_Mismatch_Search` … `Five_Mismatch_Search`,
> `Extend_Windows`, `Switch_Index`, the `Program_*` arrays, the `SEG_*` / `WIN_*` enums) are
> gone. `Mm_Builder` builds the same programs step for step, and `--batch 1` runs them
> on one lane of the engine.

## Summary

Once an SA range narrows to a single row, `batman_zero_min` extends it one base
at a time. Each step is an LF mapping, which costs one dependent cache miss.
When such a path lands on a row whose SA value is stored (`Start % 8 == 0`),
the program now stops walking:

1. It reads the stored SA value, which gives the read's position in the genome.
2. It compares the rest of the read with the 2-bit packed genome (`.pac`). For
   each 32 bases: XOR, fold each base's 2 bits into 1, popcount.

The remaining walk becomes two memory accesses. Paths reach a sampled row after
about 8 steps on average, so the saving grows with read length. The search is
21–25% faster on 32 bp reads, 42% on 40 bp, 58–68% on 75 bp and 66–76% on
100 bp. Results are unchanged: every hit is found at the same genome position,
in the same order.

The cost is memory: the SA samples (250 MB per index for a 500 Mb genome) plus
the `.pac` (125 MB). `--no-ref` keeps the old walk-only behaviour and memory
use, with byte-identical output and the same speed.

## What changes in `src/batman_zero_min.cpp`

- **`Locate`** turns a sampled row into the read's forward-genome start:
  - in `revfmi`, which indexes the reversed genome, the row stands for
    `read[0..M)`, so the start is `n − M − SA`;
  - in `fwfmi`, the row stands for the read's last `M` bases, so the start is
    `SA − (L − M)`.
  - A read that would run off the genome is dropped, as the walk drops it on
    reaching `$`.
- **`Genome_Mismatches`** does the XOR / popcount comparison against the `.pac`.
  The read is packed once per read (`Pack_Read`).
- **The four singleton walks.** `Search_Right/Left_Singleton` and
  `Extend_Right/Left_Singleton` switch to the comparison at a sampled row,
  either on entry or after a step.
  - `Extend_Located` keeps the walk's rules: a hit at exactly `Count`
    mismatches, and above `Count` the tag is deferred to the next pass, where
    the comparison is redone.
  - On a hit, `Record_Mismatches` appends the remaining mismatches to
    `Mismatch_Pos` / `Mismatch_Char` in the order the walk would have met them
    (ascending read position on the right, descending on the left). A located
    hit therefore carries the same mismatch fields as a walked one.
  - A located tag keeps `Skip` set and its genome start in `SARange::Pos`. The
    full set of fields and how they pass from search to extension is documented
    at `SARange`. The search halves hand a located tag to the extension halves,
    and the 0|0 pass hands it to the 0|1 resume, by this position alone.
- **`One_Mismatch_Search`** holds the `-n 1` scan (0|0, the 0|1 resume, 1|0).
  It was moved out of `main` so the test runs the same code; `main` only formats
  the output.
- **Output.** A located hit prints `id P pos pos 1` (plus `1` for a
  one-mismatch hit), where `pos` is the 0-based genome start. Ranges that are
  not located still print `F` / `B` SA rows.
- **Loading.** The program loads `rev<g>.sa`, `<g>.sa` (with `-n 1`) and
  `<g>.pac` by default. If one is missing or the SA interval is not 8, it warns
  and falls back to the walk-only path. `--no-ref` skips them.
- **`Skip` is now only a flag.** It used to count steps since the sampled row,
  meant for a later decode stage (sampled row's SA value + Skip). Nothing reads
  the count today, and as an `unsigned char` it would wrap to 0 after 255 steps
  on long reads. A row-based decode can recover it the way `BWTSaValue` does,
  stepping from the final row until it reaches a sampled row. Alternatively,
  store the `Level` at the moment of sampling: the count is the final `Level`
  minus that. The located path makes it unnecessary, because it records the
  genome position itself.

## Verification

`make test-located` (in `src/`) runs `tests/test_located_hits.sh`:
- **Toy reference:** two contigs with a repeated segment. Reads at 12, 20, 36
  and 75 bp: exact, one mismatch in either half, two mismatches, and reads
  crossing the contig junction or overhanging the genome ends.
- **`tests/test_located_hits.cpp`:** runs the `-n 0` search and
  `One_Mismatch_Search` with the reference off and on, for `-m 0` and `-m 1`.
  Per pass, it requires the same hits in the same order, at the same genome
  positions (walked rows decoded with `BWTSaValue`), with the same hit counts and
  the same `Mismatches` / `Mismatch_Pos` / `Mismatch_Char`. Every position must
  be at the right Hamming distance from the read, and the recorded mismatches
  must be where the read and the genome differ.
- **End to end:** compares per-read hit counts of `batman_zero_min` with and
  without `--no-ref`.
- **Sensitivity:** disabling `Record_Mismatches`, or shifting `Locate`'s
  position by one, makes the test fail (354 and 709 failures).

Checks run while developing it:

- **`--no-ref`:** output byte-identical to the previous binary over 48 runs.
  These were the small-genome sets (reads.fa, mut12, mut36, 12/20/36 bp
  genome-end overhangs), with `-n 0` and `-n 1 -m 0/1/5`, standard and
  paired-rank builds.
- **Default mode, same 48 runs:** every output line was decoded to genome
  positions (`F` rows as `n − L − SA_rev`, `B` rows as `SA_fw`, using full
  `BWTSaValue`). The decoded lists are identical to the old binary's.
- **Every decoded position** was checked against the genome: Hamming distance 0
  on exact lines, 1 on one-mismatch lines, no position repeated for a read. The
  brute-force counts (`check.py`) also match.
- **500 Mb genome:** the same decode-and-compare check on 20k-read samples of
  r32, p32_mm1L, p40_mm2, p75_mm0, p75_mm1R, p100_mm1L and p100_mm2, with
  `-n 0 -m 0`, `-n 1 -m 0` and `-n 1 -m 1`. All 21 runs were identical and every
  position checked out.

## Timing

500 Mb genome (20 contigs, about 5% repeat families). Search time in seconds,
mean of 2 runs, pinned to one core (`taskset -c 2`), Intel Core Ultra 7 165H.
"Old" is the build before this change.

| reads | options | old | new | change | `--no-ref` |
|---|---|---:|---:|---:|---:|
| r24 (2M) | `-n 1 -m 0` | 22.33 | 18.76 | −16.0% | 19.14 |
| r32 (2M) | `-n 1 -m 0` | 15.10 | 11.35 | −24.8% | 14.80 |
| r32 (2M) | `-n 1 -m 1` | 9.47 | 7.32 | −22.7% | 9.50 |
| p32_mm1L (1M) | `-n 1 -m 1` | 6.20 | 4.88 | −21.2% | 6.14 |
| p32_mm2 (1M) | `-n 1 -m 0` | 4.20 | 3.96 | −5.7% | 4.20 |
| p40_mm0 (1M) | `-n 1 -m 0` | 12.08 | 7.05 | −41.7% | 12.08 |
| p75_mm0 (500k) | `-n 0 -m 1` | 5.57 | 1.77 | −68.2% | 5.61 |
| p75_mm0 | `-n 1 -m 0` | 11.90 | 3.88 | −67.4% | 11.94 |
| p75_mm1L | `-n 1 -m 0` | 7.31 | 3.07 | −58.0% | 7.19 |
| p75_mm1R | `-n 1 -m 0` | 7.23 | 3.07 | −57.6% | 7.17 |
| p75_mm2 | `-n 1 -m 0` | 4.24 | 2.59 | −39.0% | 4.29 |
| p100_mm0 (500k) | `-n 0 -m 1` | 7.46 | 1.80 | −75.9% | 7.46 |
| p100_mm0 | `-n 1 -m 0` | 15.86 | 3.93 | −75.2% | 15.80 |
| p100_mm1L | `-n 1 -m 0` | 9.62 | 3.21 | −66.6% | 9.58 |
| p100_mm1R | `-n 1 -m 0` | 9.62 | 3.25 | −66.2% | 9.59 |
| p100_mm2 | `-n 1 -m 0` | 5.66 | 2.78 | −51.0% | 5.66 |

- **`--no-ref`** is within noise of the old build, except r24, where it was
  14% faster in both runs. That gain is not from the reference path. In the old
  build GCC inlined `Extend_Left` into the read loop in `main`; the larger
  singleton code made it stop doing so.
  - Forcing `Extend_Left` out of line in the old source alone takes r24 from
    21.3–21.6 s to 18.3–18.4 s.
  - Neither the `Skip` change (21.3–21.7 s) nor the larger `SARange` (21.1–21.3 s)
    made a difference.
  - Both `Extend_Left` and `Extend_Right` are now marked `noinline`, so the
    gain no longer depends on GCC's inlining limits.
  - The other read sets barely use the one-mismatch left extension, which is
    why only r24 (12 bp halves) showed it.
- **Gain grows with read length:** at 75–100 bp, walks that would continue for
  tens of steps end after about 8.
- **Smallest gains** are on reads with no hit (p32_mm2, p75/p100_mm2 under
  `-n 1`). Most of their paths die within a few steps, before reaching a
  sampled row.

## Memory and load time

Peak resident memory on the 500 Mb genome (`/usr/bin/time -v`). Index load times
are with the files in the page cache.

| mode | `-n 0` | `-n 1` |
|---|---:|---:|
| default (`.sa` + `.pac`) | 495 MB, 0.19 s | 868 MB, 0.31 s |
| `--no-ref` | 138 MB, 0.05 s | 273 MB, 0.09 s |

The SA samples take `4n/8` bytes per index and the `.pac` takes `n/4` bytes.
For a 3 Gb human genome that is about 1.5 GB per index plus 0.75 GB, so
`--no-ref` stays useful on small machines.

## Tried and reverted: settling located paths at once

**The idea.** A located path's full Hamming distance D is known at its first
comparison (`Extend_Located`). When D > Count, it is deferred unchanged, and
each later round of its chain resumes and compares it again. Instead:
- **D > n** (beyond `-n`): drop it, since it can never be a hit. This applies at
  any `-m` and keeps the output byte-identical.
- **D ≤ n with `-m 0`:** record it as a hit at once, carried to the round with
  Count = D on a per-chain list, so it keeps its stage.

**Correct.** The `-n 4` brute-force oracle passed on the 1 Mb and 2 Mb
repetitive genomes, with 0 duplicates, and batch 1 equalled batch 16.

**Not worth it.** 500 Mb, 20k mixed reads, `-m 0`, `--batch 16`, settling off →
on:

| reads | `-n` | lane steps per read | located re-entries per read | carried per read: D ≤ n / D > n |
|---|---|---|---|---|
| 100 bp | 2 | 119.0 → 118.7 | 2.39 → 2.29 | 0.07 / 0.10 |
| 100 bp | 3 | 160.8 → 159.9 | 3.29 → 3.01 | 0.20 / 0.21 |
| 100 bp | 4 | 283.7 → 282.4 | 3.66 → 3.26 | 0.29 / 0.22 |
| 75 bp | 4 | 411.4 → 410.6 | 3.52 → 3.24 | 0.20 / 0.15 |
| 32 bp | 3 | 924 → 911 | 16.8 → 12.5 | 0.03 / 14.5 |
| 32 bp | 4 | 8312 → 8269 | 143 → 129 | 0.03 / 141 |

- **Little to save:** 0.1–1.4% of lane steps, and search time within noise.
- **Few tags are carried.** A located tag is carried only if it passes Count in
  a round before the last. On long reads that is 0.1–0.5 tags per read.
- **Most tags beyond `-n` come from the final round,** whose rejects were already
  dropped at the stage's end.
- **The 32 bp `-n 3` / `-n 4` cost is elsewhere:** in the wide seeds (thousands
  of steps per read), not in carried located tags.
- **The cost of keeping it:** search parameters set at every entry point, an
  extra list threaded through the chains, and a `-m 0` record order that differs
  within a stage. So it was reverted.
