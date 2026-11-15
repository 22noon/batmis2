# Bidirectional FM index (step 1: exact search)

A bidirectional FM index built on the BatMis index pair, developed here on its own,
without touching `src/`.

- **The goal:** let a match grow left or right at any point, without re-walking it in the
  other index, as `batman_zero_min`'s `Mm_Start_Switch` does.
- **Step 1:** the basic operations, unoptimised, checked and timed against the classic
  one-index searches.

## The index

`bwt_build_index` writes, from the same text T (`<g>.bin`):
- `<g>.bwt/.fmv/.sa`, the FM index of T (Fw);
- `rev<g>.bwt/.fmv/.sa`, the FM index of T reversed (Rev).

That is the pair a bidirectional index needs, so there is no new format. `Bi_Index`
(`bi_index.hpp`) loads both through `../index/libfmindex_core.a`.

**Conventions** (as `index/BWT.h`):
- rows run from 0 to n (n = `textLength`), and ranges are inclusive;
- the `$` row is `inverseSa0`, and its character is not stored;
- `Occ(i, c)` counts c in rows [0, i), skipping `$`.

**A pattern P has two ranges of the same size:**
- **Fw:** T's suffixes that start with P;
- **Rev:** reversed T's suffixes that start with P reversed.

`Bi_Range` keeps {Fw_Start, Fw_End, Rev_Start}. The empty pattern is [0, n] in both.

**`Extend_Left(R, c)`** (P → cP):
- **Fw:** the classic backward step,
  Start = C[c] + Occ(Start, c) + 1, End = C[c] + Occ(End + 1, c).
- **Rev:** P reversed's rows are ordered by the base that follows it in reversed T, which
  is the base before P in T: the Fw BWT character of the matching row. So the new Rev
  range starts after the rows whose character is smaller than c:
  `Rev_Start += [Fw's $ row in Fw's range] + Σ_{b<c} (Occ(End + 1, b) − Occ(Start, b))`.
- **The `$` row** is the occurrence of P at T's start. Its reversed suffix is followed by
  the end of reversed T, so it sorts first.
- **Cost:** `BWTAllOccValue` at both ends gives all four counts.

**`Extend_Right(R, c)`** (P → Pc) is the mirror: a backward step in Rev, with Fw shifted
by Rev's counts.

`Locate_Fw(row)` is `BWTSaValue` in Fw. `Locate_Rev(row, |P|)` = n − SA_rev(row) − |P|.

## Programs

- **`bi_check <prefix> [patterns] [seed]`:** grows random patterns of 1–60 bases from a
  random point, by a random mix of `Extend_Left` and `Extend_Right`. The patterns come
  from the text (many at its start or end, the `$` cases) or are random strings. After
  every step it checks:
  - Fw equals the classic backward search in Fw, and Rev the classic search in Rev;
  - the sizes are equal, and equal the brute-force count in the text (patterns of up to 12
    bases);
  - for up to 16 rows, the positions decoded from Fw and from Rev are the same set, and
    the text holds the pattern at each.
- **`bi_exact -g <prefix> -q reads -mode M [-m N] [-o out] [--timing]`:** exact search,
  with M one of the following (the `bi-*` modes also take `--lookup k`, `--ref`,
  `--batch B`, `--threshold N` and `--no-singleton`; see below):
  - `backward`: classic, Fw only;
  - `forward`: classic, Rev only;
  - `bi-left`: `Extend_Left` from the last base;
  - `bi-right`: `Extend_Right` from the first base.

  Output: one line per read, with the read id, the occurrence count and the first N
  positions (0-based in T, ascending). With `--timing` it prints search and locate
  seconds, and ns per read and per base.

- **`bi_mm -g <prefix> -q reads [-m N] [-o out] [--lookup k] [--ref] [--batch B] [--counts] [--timing]`:**
  the one-mismatch search (below). `-m N` stops a read's search at N hits (0, the
  default: all). Output: one line per read, with the read id, the exact count, the
  one-mismatch count, then every position as `pos:mismatches` (`--counts`: the counts
  only).
- **Locate options** (`bi_exact` bi-* modes and `bi_mm`; see "Locate"): `--locate lib|bi`,
  `--chain-min N`, `--split`, `--one-end`.

```sh
make          # bi_check, bi_exact, bi_mm (builds ../index/libfmindex_core.a)
make test     # toy index: bi_check, every bi_exact mode and bi_mm against brute force (test.sh)
```

## Correctness

- **`bi_check`:** no failure.
  - 1 Mb random genome: 55,860 steps.
  - 2 Mb repeat-rich genome: 28,242 steps.
  - `test.sh` toy reference: 86,213 steps. It has three contigs, the shortest 300 bases,
    with repeats.
- **`test.sh`:** every `bi_exact` mode, at `-m 0`, equals brute force over `.bin`. The
  reads are taken from the text, have one mismatch, are random, or sit at the text's two
  ends, at lengths 8–75, plus a read with an N.
- **Against `batman_zero_min`:** 500 Mb genome (20 sequences), 100k 100 bp reads.
  - `bi_exact -mode bi-right -m 0` gives the same positions as
    `batman_zero_min -n 0 -m 0 --sam --forward-only` for every read (28,996 occurrences),
    once SAM's positions are moved to the concatenated text.

## Timing (step 1, unoptimised)

**Setup:**
- 2M reads that occur exactly in the 500 Mb random genome, at 100 bp and at 32 bp.
- One P-core (CPU 1), search seconds only (no loading, no locating), `-m 1`.
- **`batman_zero_min`:** `-n 0 --batch 1 --lookup 3 -t 1`, unbatched with only a 3-base
  table, with and without `--no-ref`.

**Noise:** other users kept the load average at about 6, and the run's wall time stretched
to six hours. In round 1 at 100 bp, the first three modes took about 63 s, twice round
2's time. Those are left out; the table gives round 2 at 100 bp, and both rounds at 32 bp.

| | 100 bp: per read | per base | 32 bp: per read | per base |
|---|---:|---:|---:|---:|
| classic backward (Fw) | 16.9 µs | 169 ns | 5.5–8.1 µs | 172–252 ns |
| classic forward (Rev) | 17.9 µs | 179 ns | 6.4–7.9 µs | 199–246 ns |
| bi-left | 18.3 µs | 183 ns | 8.9–9.1 µs | 277–283 ns |
| bi-right | 16.3–17.4 µs | 163–174 ns | 10.4–10.6 µs | 325–331 ns |
| `batman_zero_min --no-ref` | 15.8–17.9 µs | | 7.5–8.2 µs | |
| `batman_zero_min` (reference) | 3.7–4.7 µs | | 5.9–6.5 µs | |

- **At 100 bp, keeping the second range costs almost nothing.** After the first 15 or so
  bases the range is a single row, so the two `BWTAllOccValue` calls at Start and
  End + 1 read the same occurrence block as the classic step's two `BWTOccValue` calls.
- **At 32 bp, it costs 1.2–1.9× the classic step.** More of each read is matched while
  the range is still wide, where Start and End + 1 are in different blocks. There each
  step adds the all-base count (`BWTAllOccValue`, SSE) to the same two cache misses.
- **The classic searches match `batman_zero_min --no-ref`** (unbatched: per-read lookups,
  a singleton walk). So these timings are a fair baseline.
- **With the reference, `batman_zero_min` is 4× faster at 100 bp.** A path that reaches a
  sampled SA row is finished by comparing with the packed genome
  (`../benchmark/pac_singleton.md`). The bidirectional index can use the same shortcut.

## Singleton step

**Once a match is a single row, the other index's row does not move.** Take Fw = [r, r]
(one occurrence, so Rev = [r′, r′]) and extend left with c:
- **Fw:** the occurrence extends only if r's BWT character is c. The `$` row is P at T's
  start, which extends with nothing. Then Fw = LF(r) = C[c] + Occ(r, c) + 1.
- **Rev:** the general step adds `[$ in range] + Σ_{b<c} count_b` over [r, r]. Both terms
  are 0, because the row's character is c. So r′ is unchanged.
- **Extending right** is the mirror: Rev moves and Fw stays.

`Extend_Left` and `Extend_Right` take this path for single rows (`Singleton_Step`). The
other index's row is copied, not computed. `bi_exact --no-singleton` turns it off, and
`bi_check` compares every single-row step with `Extend_General`: 93,224 steps on three
references, all equal.

**Share of steps taken from a single row** (exact reads): 82% at 100 bp and 50% at 32 bp,
on the 500 Mb genome.

**Which Occ to use matters more than the saved update.** Dependent chains of 100 LF steps
from random rows, one P-core, ns per step (`scratchpad/bi/micro/lfbench.cpp`):

| step | ns |
|---|---:|
| BWT character + `BWTOccValue` | 400–418 |
| `BWTOccValueOnSpot` | 384–419 |
| BWT character + `BWTAllOccValue` ×2 (the general step's lookups) | 293–326 |
| BWT character + `BWTAllOccValue` ×1 | 239–269 |

`BWTOccValue` (one base) is about 40% slower than `BWTAllOccValue` (all four), so
`Singleton_Step` uses `BWTAllOccValue`. The first version used `BWTOccValue`, and was
20% slower than the general step.

**Exact search, 2M reads, one P-core, two rounds, µs per read.** The machine ran slower
than in the step 1 timing; compare within this table.

| | 100 bp | 32 bp |
|---|---:|---:|
| classic backward | 25.5–25.6 | 7.3–7.4 |
| bi-left, general step only | 28.3–28.9 | 8.0–8.4 |
| bi-right, general step only | 30.4–31.1 | 8.7–9.0 |
| bi-left, singleton step | 26.2–26.8 | 7.7–7.9 |
| bi-right, singleton step | 26.9–27.3 | 8.0–8.2 |
| `batman_zero_min --no-ref`, unbatched | 24.5–24.7 | 7.2–7.4 |

- **The singleton step saves 7–12% at 100 bp** and 4–8% at 32 bp, against the general
  step.
- **The bidirectional search is now within 3–7% of classic backward at 100 bp**, and
  5–11% at 32 bp. The rest is the general step's cost on the first, wide steps of each
  read (18% of the steps at 100 bp, 50% at 32 bp).
- **The classic step needs no BWT character read,** since the read gives the base. Its
  two `BWTOccValue` calls fall in the same block, so a real search does not show the gap
  the chain benchmark shows.

## batman_zero_min's exact-search optimisations, ported

`batman_zero_min` (`../src`) stays fixed as the baseline. Each of its `-n 0`
optimisations maps onto a bidirectional form:

| `batman_zero_min` | Bidirectional | Where |
|---|---|---|
| Singleton walk | singleton step (above) | `Bi_Index::Extend` |
| `Count_Range` below `BRANCHTHRESHOLD` | narrow step | `Extend_Narrow` |
| `--lookup k`: two tables, one per index | one table of `Bi_Range`s, both directions | `Build_Lookup`, `Look_Up` |
| located single rows (`.pac`) | either row of a singleton | `bi_ref.hpp` |
| `--batch B` | lanes with per-step prefetch | `bi_batch.hpp` |

### Narrow step

**The counts give both ranges.** For a range below `Threshold` rows (default 128),
`Count_Range` counts each base's BWT characters directly: Count[c] is the new range's
size, and `[$ in range] + Σ_{b<c} Count[b]` is the other index's offset. One
`BWTAllOccValue` at Start gives the new start. So the second range costs nothing beyond
the classic narrow step. `Children_Left` / `Children_Right` return all four extensions
from the same counts; `Build_Lookup` uses them, and so will the mismatch search.

**Timing:** no measurable effect at any threshold (0, 32, 64, 128, 256: within ±5%
noise, load average about 6). Only 6.5% of the steps at 100 bp, and 15% at 32 bp, start
from a narrow range on this genome; `../benchmark/branch_threshold.md` found the same
for `batman_zero_min`. It matters on repetitive genomes, where ranges stay wide longer.

### Lookup table

**One table serves both directions.** Entry k-mer P (key: P's bases, first in the high
bits) holds P's `Bi_Range`, {Fw_Start, Fw_End, Rev_Start}, 16 bytes, so one cache line
per lookup. `Search_Right` starts from the read's first k bases, `Search_Left` from its
last k. The same table can start a seed anywhere in the read.
- **Build:** a depth-first walk by `Children_Left`, 0.73–0.84 s at k = 11.
- **Size:** 4^k × 16 bytes, 64 MB at k = 11 (`batman_zero_min`'s two tables:
  2 × 4^k × 8 bytes).

### Located singletons

**A singleton has two rows, and either can be sampled.** Fw row r of read[i, j) gives
the read's start SA_fw[r] − i; Rev row r′ gives n − SA_rev[r′] − j. The read then matches
if the packed genome there holds it (32 bases at a time).
- **Extending right,** the Fw row never moves: it is tested once, when the pattern
  becomes a single row, and the Rev row after each step. Extending left is the mirror.
- **The fixed row locates 14% of the reads:** on the 100k mixed reads, 11.4k of 80k
  located reads (68.5k by the moving row). `batman_zero_min` has only the moving row.
- `bi_exact --ref`; it needs SA sampling every 8 rows in both indexes.

### Batched lanes

`Batch_Exact` (`bi_batch.hpp`) runs B reads at once. Each lane's search is cut into steps
that read one thing from memory:
- the lookup entry;
- a step's occurrence rows (`Prefetch_Step`: the own row for a singleton; Start and the
  last `bwtCode` word for a narrow range; Start and End + 1 for a wide one);
- a sampled SA value;
- the genome words.

A round prefetches every lane's next step, then takes each lane's step, as
`batman_zero_min`'s `Batch_Search`. `bi_exact --batch B`. With `--early`, each lane
prefetches its next step right after taking one, so that every lane's prefetch has a whole
round to arrive.

### Inline Occ and the batched profile

**First timing, library Occ** (2M exact reads, one P-core, `--lookup 11`, with the
reference; search seconds, two rounds, load average 5.6–9.6):
- **Unbatched:** the bidirectional search equalled `batman_zero_min` at 100 bp
  (13.3–16.5 s against 13.9–17.2 s) and at 32 bp (12.1–20.4 against 12.7–16.5).
- **Batched:** it gained little. At B = 8, 100 bp: 13.0–14.2 s against 9.0–10.6 s.

**Profile** (`perf record`, B = 8; a `-DBI_PROFILE` build keeps the three step kinds out
of line):

| % of samples | 100 bp | 32 bp |
|---|---:|---:|
| `Singleton_Step` (own code: the BWT character load) | 17 | 16 |
| `BWTAllOccValue` + `BWTDecodeAll` (nearly all from single rows) | 38 | 33 |
| `Extend_Narrow` | 3.1 | 2.8 |
| `Extend_General` (wide ranges) | 0.06 | 0.05 |
| `Batch_Exact` (lanes; stalls on the SA, genome and lookup loads) | 25 | 21 |

- **Wide ranges barely occur** once the k = 11 table starts the search.
- **Narrow ranges are cheap.**
- **The cost was the single-row step:** a four-base library decode (`BWTDecodeAll`, SSE
  through partition tables) where one base is needed. The step itself needs no more
  memory than the classic walk: LLC misses per read were within 10% of
  `batman_zero_min`'s.

**Inline Occ** (`Occ_All`, `Occ_One` in `bi_index.hpp`):
- **Method:** the nearest checkpoint (major + 16-bit minor count, as `BWT.c`'s
  `BWTAllOccValueExplicit`), plus or minus a direct count of the at most 128 characters
  to the row. All four bases use `Count_Range`'s bit-planes. One base is two
  `popcnt` a word, on `~(word ^ c's pattern)`.
- **No tables, and inlined:** the single-row and narrow steps call `Occ_One`, and wide
  steps call `Occ_All`.
- **Check:** `bi_check` compares both with the library at every row of the toy, 1 Mb and
  2 Mb indexes, and at 20.6M rows of the 500 Mb pair. All equal.
- **Effect:** 10% fewer instructions, and search 6–15% faster than with the library Occ.
- `bi_exact --occ lib` keeps the library calls.

**SIMD for the other index's rank:**
- **All four children** (`Children_*`: table build, and later branching): the four counts
  are one SSE vector. Their exclusive prefix sum `[0, n0, n0+n1, n0+n1+n2]` (two shifts
  and adds), plus Other_Start and `[$ in range]`, gives all four other-index starts.
- **One base:** the sum of at most three counts stays scalar; it does not show in the
  profile.

**Timing now** (inline Occ, `--lookup 11`, with the reference; search seconds, two
rounds, load average 5.7–6.5):

| | B = 1 | 4 | 8 | 16 | 32 |
|---|---:|---:|---:|---:|---:|
| 100 bp bi-right | 16.9–17.8 | 12.0–12.2 | 11.8–12.0 | 11.9 | 11.5–11.6 |
| 100 bp bi-right `--early` | 16.6–17.3 | 11.8–12.2 | 11.5–12.1 | 11.5–12.2 | 11.5–12.1 |
| 100 bp `batman_zero_min` | 18.1–18.9 | 12.0–12.2 | 12.2 | 11.5–12.0 | |
| 32 bp bi-right | 12.1–12.7 | 8.7–9.8 | 8.8–9.3 | 8.5–8.6 | 8.3–9.0 |
| 32 bp bi-right `--early` | 12.9–13.9 | 8.6–9.9 | 8.8–9.8 | 8.4–8.6 | 8.8–9.4 |
| 32 bp `batman_zero_min` | 13.3–13.4 | 8.5–8.7 | 8.7–11.5 | 9.4–9.5 | |

- **The bidirectional search now matches `batman_zero_min` at every B:**
  - unbatched, it is 5–10% faster, since the fixed row locates 14% of the reads a few
    steps sooner;
  - batched, they are equal within noise.
- **`--early` prefetch makes no measurable difference,** so it stays off by default.
  Neither does B beyond 8.

### Correctness

- **`bi_check`** (toy, 1 Mb, 2 Mb repeat genomes): every narrow step equals
  `Extend_General`; every step's `Children_*` equal the four `Extend_*`; every lookup
  entry (4^6 or 4^8) equals the k-mer searched base by base; the inline Occ equals the
  library's (above).
- **`test.sh`:** 25 `bi_exact` runs against brute force: each mode; the bidirectional ones
  with `--threshold 0`, `--lookup 6` and `10` (longer than the 8-base reads), `--ref`,
  `--batch 1–16` with and without `--ref` and `--lookup`, `--early`, and `--occ lib`.
- **500 Mb genome, 100k mixed reads, `-m 0`:** 12 configurations (both directions ×
  lookup, ref, batch, and combinations) give identical output, and
  `bi-right --ref --lookup 11 --batch 8` equals `batman_zero_min` for every read.

The same inline Occ was tried in `batman_zero_min` (`src/fm_occ.hpp`): identical output,
6–9% fewer instructions but no fewer cycles, and `-n 1` 10–15% slower. There it stays off
by default (`-DBATMAN_INLINE_OCC`; `benchmark/inline_occ.md`).

## One-mismatch search (`bi_mismatch.hpp`, `bi_mm`)

The hits of `batman_zero_min -n 1`, on one `Bi_Range` per node. LH = P[0, L/2),
RH = P[L/2, L).

- **Stage A (0|0, 0|1):** LH by `Extend_Right`, then RH along the exact path. Each node's
  `Children_Right` gives, from one set of counts, the exact child and the other three
  (one mismatch placed); those are kept, and after the path each is extended exactly to
  the read's end. This is `batman_zero_min`'s deferral (`Branch_Detect`, then the 0|1
  resume), without a second pass.
- **Stage B (1|0):** RH by `Extend_Left`, then LH leftward the same way with
  `Children_Left`. If stage A's LH search failed at depth d, a 1|0 hit's mismatch is
  before d, so only children before d are kept.
- **One index for every hit:** both stages give each hit's range in Fw (and Rev).
  `batman_zero_min` reports 0|x hits in revfmi and 1|0 hits in fwfmi. One lookup table
  serves both half searches.
- **Located singletons** (`--ref`), everywhere: a single row whose Fw or Rev SA value is
  sampled gives the read's start (`SA_fw − i` or `n − SA_rev − j` for P[i, j)), and one
  comparison of the whole read with the genome (`Packed_Genome::Diff`) ends the node:
  - stage A's half search and path: a hit if it differs nowhere (exact) or once in RH;
  - stage B's: a hit if it differs once in LH;
  - a kept child's extension: a hit if it differs only at the child's mismatch.

  Locating during the half searches matters: on 500 Mb a 50-base half is unique after
  about 15 bases, and without it the remaining 35 single-row steps of each half took
  100 bp searches from 10 s to 3.6 s (500k reads).
- **Steps** (`Mm_Search`): the search is cut into steps that each read one thing from
  memory (lookup entry, step rows, SA sample, genome), so `Batch_One_Mismatch` runs B
  reads in lanes as `Batch_Exact` does.
- **One read at a time** (`One_Mismatch`, `--batch 0`): the same steps in the same order
  as a plain loop, with the same hits, chains and prefix/suffix anchors. Byte-identical to
  `--batch 1` and `8` in 48 runs (output and located counts). Running the steps one at a
  time cost 9–14% of the search: 3.43 s against 2.61–2.67 s at 100 bp, 2.76 against 2.52–2.60
  at 32 bp, with `batman_zero_min --batch 1` at 3.02 and 2.48–2.57.
- **`-m N`:** checked between the exact hits, the 0|1 and 1|0 searches, and the kept
  children, so the reported hits are those of the earlier stages first, as
  `batman_zero_min`'s.

**Per read** (500k reads, 500 Mb, `--ref --lookup 11`): 100 bp: 2 lookups, 22 steps
(10 single-row, 11 narrow, 0.5 wide), 1.7 SA samples and genome comparisons, 0.14 kept
children; 32 bp: 2, 18 (8, 9, 0.5), 1.5, 0.5.

### Correctness

- **Brute force** (`test.sh`, and a 71 kb toy with an `AC` repeat): reads of 1–100 bp,
  exact, one or two mismatches, random, at the text's ends; `-m 0` equals every position
  with at most one mismatch, one read at a time and batched (B = 1–16), with and without
  the table and `--ref`.
- **`batman_zero_min -n 1 -m 0 --sam --forward-only`**, 100k mixed reads, 500 Mb: the same
  (position, NM) set for every read at 32 bp, and at 100 bp but one read, whose 1|0 hit
  crosses the chr18/chr19 boundary (a hit in T that SAM cannot hold).
- **The steps against the first version** (a plain loop per read): byte-identical in 40
  runs (100/32 bp × `--ref`/lookup combinations × `-m 0/1/3` × B = 0, 1, 8, 16).

### Timing

Search seconds, 500k reads, 500 Mb genome, `--lookup 11`, `--ref`, one CPU (`taskset -c
1`, a shared machine: about ±10%). `batman_zero_min` at `-n 1 -t 1`.

| run | bz B = 1 | bz B = 8 | bi, plain loop | bi steps, one lane | bi B = 8 | bi B = 16 |
|---|---:|---:|---:|---:|---:|---:|
| 100 bp `-m 0` | 4.5–4.8 | 2.9–3.0 | 2.9–3.2 | 3.9–4.4 | 3.1–3.6 | 3.3–3.6 |
| 100 bp `-m 1` | 3.5 | 2.2–2.7 | 2.8 | 3.2–3.6 | 2.4–2.7 | 2.5–2.6 |
| 32 bp `-m 0` | 3.7–3.9 | 2.6–2.9 | 3.4 | 3.5–3.6 | 2.9–3.6 | 2.9 |
| 32 bp `-m 1` | 2.4–2.6 | 1.6–1.8 | 2.2 | 2.3–2.4 | 1.8 | 1.7–1.9 |

- **One read at a time, the bidirectional search is faster** (the plain loop, the first
  version, before the steps: 25–35% faster at 100 bp, 10% at 32 bp). It has the same hits
  with fewer steps: no second index, no deferral pass, and the located shortcut in both
  rows.
- **Batched, it is level with `batman_zero_min` at 32 bp and 10–20% behind at 100 bp.**
- **The steps cost about 1 s unbatched** (one lane against the plain loop; `--batch 0`
  is the plain loop again since the two-mismatch work: 15% more
  instructions, 25% more cycles, the same cache misses), and batching wins back less than
  for `batman_zero_min` (1.25× against 1.5×). B beyond 4–8 gains nothing, so the cache
  misses are hidden and the cost is the work per step: about 400 instructions a step
  (240 in `Extend`), not yet looked into.

## Locate (`bi_locate.hpp`)

A hit row's position comes from an LF walk to a sampled row (`row % 8 == 0`): about 7
dependent misses per row with the library's `BWTSaValue` (`--locate lib`). `Locator`
(`--locate bi`, the default) uses what the bidirectional index adds:

- **Pair walks (single rows).** A single-row range is one occurrence with a row in each
  index. LF in Fw moves its start left; LF in Rev moves its end right
  (`start = n − SA_rev − j`). Stepped together, the first to reach a sample wins:
  - **expected rounds:** about 3.3 instead of 7 (p = 15/64 a round);
  - **measured:** 3.26 rounds per pair, against 6.98 steps per walk with `--one-end`.
- **The range chain (wide hits, from `--chain-min` rows, default 24).** A search records
  each step's range in the index it steps: Fw ranges of P′[i, L) going left, Rev ranges of
  P′[0, j) going right. Both `Mm_Search` stages, the kept extensions, `Search_Left/Right`,
  `Located_Search` and the `Batch_Exact` lanes record it. That is the chain `batman_zero_min`'s
  `Chain_Decode` rebuilds at about 2L Occ per hit before walking it back (which is why it
  stayed off in `src/`); here it costs one 16-byte store per step.
  - **The backtrack is `Chain_Decode`'s:** stretches by `s_k mod 8`, word-wise scans at
    width changes, and the 256-row trial.
  - **Levels below the lookup table's depth** are rebuilt only if the backtrack gets there.
  - **Rows left after a passed trial** join the walks.
  - **A failed trial** (ancestors that stall: poly-C hits from N runs) walks the hit's Fw
    rows with the library's serial walk. Those rows stay in cache, where its per-step cost
    beats the batched walks. Walking the Rev rows instead was 9% slower there.
- **All walks of all reads** step together, at most 1024 at a time, with each walk's next
  rows prefetched.
- **Split** (`--split`, off): narrow ranges split into matched pairs by `Children_Left`.
  A child's Rev rows are a sub-range of its parent's, so a single-row child is a pair. It
  measured no faster, since the splitting runs unbatched.
- **`bi_exact -m N`** (N above 0) locates the first N rows of the chain's index: Rev for
  `bi-right` hits, as `batman_zero_min` reports from revfmi. With `-m 0` the output is
  unchanged.

**Correctness:**
- **`test.sh`:** every route against brute force, in both programs: `--locate lib`, chains on
  every hit of 2+ rows (`--chain-min 2`, with the table's rebuilt levels), `--split` and
  `--one-end`.
- **Against `--locate lib`, `-m 0`:** byte-identical in 30 runs with the final code, and in
  50 `bi_exact` runs before the failed-trial fallback changed:
  - 500 Mb (20k mixed reads, 100 and 32 bp), `rx` (50k reads, 32 and 100 bp) and `rb`;
  - `bi_exact` in both directions and `bi_mm`, ±`--ref`, ±table, batched;
  - and the rep200 poly-C reads.

**Locate seconds** (`taskset -c 5`, two rounds; `rx`: 20 Mb, 40% near-exact repeats;
`rb`: copies that differ just upstream of the reads):

| run | `lib` | `bi` | resolved by chain |
|---|---:|---:|---:|
| 500 Mb, 500k × 100 bp, `bi_mm`, no `--ref` | 0.84–0.97 | 0.45–0.47 | 944 of 944 |
| 500 Mb, 500k × 100 bp, `bi_mm --ref` | 0.10–0.12 | 0.07–0.10 | |
| 500 Mb, 500k × 32 bp, `bi_mm`, no `--ref` | 1.22–1.26 | 0.64–0.65 | 270k of 281k |
| 500 Mb, 500k × 32 bp, `bi_mm --ref` | 0.44 | 0.23–0.27 | |
| `rx` 100 bp, `bi_mm`, no `--ref` / `--ref` | 0.21–0.22 / 0.17–0.19 | 0.14–0.15 / 0.09–0.10 | 615.0k of 615.1k |
| `rx` 32 bp, `bi_mm`, no `--ref` / `--ref` | 0.28–0.31 / 0.31–0.36 | 0.23–0.27 / 0.20–0.21 | 1.04M of 1.08M |
| `rx` 100 bp, `bi_exact -m 0` left / right | 0.07 / 0.09 | 0.02–0.03 / 0.02–0.04 | 100% |
| `rx` 32 bp, `bi_exact -m 0` left / right | 0.14 / 0.14–0.16 | 0.07–0.08 / 0.06–0.08 | 97% |
| `rb` 100 bp, `bi_mm --ref` | 0.07 | 0.006–0.023 | 100% |
| rep200 poly-C, 500 × 100 bp (9.7M rows) | 29.3–30.5 | 28.9–29.6 | 0% (serial) |

- **Single rows:** 1.9× faster. Batching the walks gives most of it (`--one-end`: 0.70
  against 1.3); the second end adds 5–20%. In throughput a pair's two misses a round cost
  about as much as one walk's seven; the second end mainly halves latency.
- **Wide hits:** 1.2–3× faster on `rx`, and 3–12× on `rb`, the case where
  `batman_zero_min`'s `Sa_Decode` spans help least.
  - **Reference:** `batman_zero_min`'s SAM decode phase (`Sa_Decode`, stats build,
    `-n 1 -m 0 --forward-only`) on the same `rx` reads takes 0.45–0.47 s for 658k rows
    (100 bp) and 0.65–0.66 s for 1.16M (32 bp). That is about 0.7 and 0.56 µs a row,
    against 0.13 and 0.17 µs here with `--ref`. It is a different program and phase, so
    only indicative.
- **Poly-C:** level with the library walk, but `Sa_Decode`'s spans (one LF for a whole
  range of one character) do this case several times faster. They are the missing piece
  if N-run hits matter here.

## Two-mismatch search (`bi_two_mismatch.hpp`, `bi_mm -n 2`)

**Pass 1** is the one-mismatch search, unchanged. **Pass 2** finds the hits with exactly
two mismatches.

**The split.** Cut points a = L/3 < LH = L/2 < b = 2L/3 give parts P0 = [0, a),
P1 = [a, LH), P2 = [LH, b), P3 = [b, L). Each hit falls in exactly one search, as in
`batman_zero_min`'s `-n 2` (p is the mismatch in LH, q the one in RH):

| search | seed (exact) | then | hits |
|---|---|---|---|
| A | LH | right over RH, two | 0\|2 |
| B | RH | left over LH, two | 2\|0 |
| C1 | P0 | right: P1 one, RH one | 1\|1, p in P1 |
| C2 | P3 | left: P2 one, P1 none, P0 one | 1\|1, p in P0, q in P2 |
| C3 | P1 P2 | left over P0 one, then right over P3 one | 1\|1, p in P0, q in P3 |

**C3 turns around.** `batman_zero_min` finds that seed in one index and re-walks it into
the other (`Mm_Start_Switch`). Here the same `Bi_Range` extends either way.

**How a search walks** (depth first):
- At a position that may take a mismatch, `Children_*` gives all four bases from one set
  of counts; elsewhere it takes one `Extend_*` step.
- Each part must end with its minimum. The walk stops once the mismatches left no longer
  fit.
- **Located singletons:** a sampled single row is the node's only occurrence. It is a hit
  if the genome there differs in exactly two places that fit the search's per-part counts
  (`Packed_Genome::Diff_At`).
- **Chains:** nodes anchored at the read's start or end record their ranges, so wide hits
  get a range chain (C3's, after its turn, is in Rev from level b).

**Pass 1 settles most of pass 2.** Pass 1's LH search ends at a prefix P[0, j) that
either does not occur (or fits no read) or occurs only at read start S. Stage B gives the
same for a suffix (`Mm_Search::Prefix`, `Suffix`).
- **When the prefix settles a search:** if it is no longer than the seed of A (LH) or C1
  (P0), every hit of that search is at S. One genome comparison settles it, with no walk.
  The suffix does the same for B and C2 (`--no-settle` walks them all).
- **Coverage:** on 500 Mb, 95% of A, B, C1 and C2 at 100 bp (1.91M of 2M) are settled; 23%
  at 32 bp, whose thirds (10–11 bases) are shorter than the prefix that becomes unique.
  C3, whose seed is the middle, always walks.

**Correctness:**
- **`test.sh`:** 8 runs against brute force (± table, `--ref`, `--no-settle`, chains,
  `--locate lib`), reads of 1–75 bp.
- **`batman_zero_min -n 2 -m 0 --sam --forward-only`**, 20k mixed reads, 500 Mb: the same
  (position, NM) set for every read at 100 and 32 bp, with and without settling and the
  table.
- **A bug the 32 bp comparison caught:** where the lookup table jumped past the end of the
  first part, that part never ended, and mismatches were refused.

**Bounds from pass 1's failures.** A prefix P[0, d) that does not occur (or fits no read)
means every hit has a mismatch before d; a failed suffix P[L − e, L) means one at or after
L − e (`batman_zero_min`'s LH_Fail / RH_Fail windows). `Bounds_Possible` prunes a node
that can no longer place them (`--no-bound` to compare).
- The output is unchanged.
- On genome reads it prunes almost nothing (6.38M to 6.37M nodes at 32 bp), since their
  prefixes and suffixes occur. It is aimed at reads with no near match, which were not
  timed.

**Timing** (search seconds, 500k reads, 500 Mb, `--lookup 11`, one CPU, three interleaved
rounds, no decode in either: `bi_mm --counts`, `batman_zero_min` text output to
`/dev/null`):

| run | bz `-n 1` | bi `-n 1` | bz `-n 2` | bi `-n 2` |
|---|---:|---:|---:|---:|
| 100 bp `-m 0` | 2.93–3.20 | 2.63–2.83 | 4.94–5.18 | 4.93–5.20 |
| 32 bp `-m 0` | 2.58–2.65 | 2.64–2.73 | 15.0–16.3 | 17.3–17.8 |

- **100 bp: level.** Pass 1 is 10% faster, pass 2 about 15% slower.
- **32 bp: pass 2 is 13–15% slower.** Both take about the same steps per read
  (`batman_zero_min` 142, of which about 94 are in its 1|1 window searches; bi about 150,
  of which 113 are in C1–C3), so the difference is the cost per step: about 290 ns
  against 210.
  - **Profile:** `Count_Codes` 41% (narrow ranges and Occ, mostly cache-miss stalls),
    singleton steps 15%, the walk 15%; 4× `batman_zero_min`'s branch misses (the
    recursive walk).
  - **The narrow threshold** (128, 64, 32, 16, 0) changes nothing; 0 is 10% worse.
**Batched** (`Batch_Two_Mismatch`, `--batch B`): each lane runs pass 1 (`Mm_Search`), then
pass 2 as `Mm2_Search`. That is the same walk with an explicit stack, one memory access a
step (lookup entry, a node's children, SA sample, genome check, including the settling
checks).
- **Visiting order:** children are pushed last first and the first is entered at once,
  so nodes are visited in the recursive walk's order.
- **Byte-identical to `--batch 0`** in 84 runs (B = 1, 8, 16; output, located and pass-2
  counts), and `test.sh` runs it against brute force.

Search seconds, 500k reads, quiet machine, three interleaved rounds:

| run | bz B = 1 | bz B = 8 | bz `auto` | bi B = 0 | bi B = 4 | bi B = 8 | bi B = 16 |
|---|---:|---:|---:|---:|---:|---:|---:|
| 100 bp `-m 0` | 4.18–4.27 | 2.76–2.80 | 2.78–2.85 | 4.11–4.35 | 3.10–3.26 | 3.15–3.21 | 3.14–3.24 |
| 32 bp `-m 0` | 12.09–12.17 | 8.82–9.51 | 8.83–9.42 | 12.89–13.48 | 11.28–11.75 | 11.09–11.68 | 11.17–11.88 |

- **Batched, bi is 12–15% behind at 100 bp and about 25% behind at 32 bp.** Its batching
  gain is 1.3× and 1.15×, against `batman_zero_min`'s 1.5× and 1.35×.
- **The counters show a compute gap, not a memory one.** At 32 bp, B = 8, search only, an
  earlier run, against `batman_zero_min`:
  - +31% cycles and +29% instructions;
  - 2.25× the branch misses;
  - 12% fewer LLC misses.

  With the misses overlapped, the work per step decides it.
- **Profile, 32 bp, B = 8:** `Count_Codes` 39% of cycles and 22% of branch misses (mostly
  stalls on the BWT words; the library Occ is slower: 1.19 s against 1.15 s on 50k
  reads), singleton steps 14%. The walk itself (`Enter`, `Children_Of`, the lane loop) is
  about a third of the instructions.
- **Entering the first child directly** (no push and pop of 64-byte nodes) gained 4–5%.
  The rest is spread: branch misses cost about 12% of cycles in all, and a branch-free
  `Count_Codes` would save about 2.5%.

- **Before settling,** the five seed searches made pass 2 2.6× slower. An earlier
  three-search split on thirds alone (T1 exact over T2 T3, …) was slower again at 32 bp
  (47 s against 21 s), as its 10-base seeds branch over about 500 rows.
- **A first timing on a loaded machine** (load average 6–7) had bi 15% slower at 100 bp;
  that was pass 1 run as steps, now the plain loop.

## batman_zero_min ported to bidirectional ranges (`bi_batman.cpp`)

A copy of `src/batman_zero_min.cpp` (same options, I/O, SAM output, threads and decode)
whose ranges also carry `Other_Start`: the start of the same string's range in the other
index.
- **Where it changes** (single-row steps leave it alone):
  - a step moves it past the range's `$` row and its rows with a smaller base
    (`Get_SARange`, `Branch_Detect`);
  - a 1-base range starts at the same row in both indexes;
  - an index switch ends at the seed's own `Start`.
- **The one extra cost:** wide steps (128+ rows) take all four bases' Occ at both ends
  instead of the added base's alone.
- **The O(1) switch** (default; `--switch-walk` for the old way): an index switch
  (`Mm_Start_Switch`) re-walks the seed's string in the other index. With `Other_Start`
  that range is known at once (`Switched_Range`), so the cross step extends it directly.

**Correctness:**
- **`bi_batman_verify`** (`-DBI_VERIFY_OTHER`) checks each range of up to 32 rows: the
  other index's rows from `Other_Start` hold the same occurrences, a string length apart.
  - About 100M ranges checked (`-n 1–7`, 100/32 bp, and `rx`), with no failures.
  - Every switch walk ended at the seed's `Other_Start` (4.9M checks).
- **Byte-identical to `batman_zero_min`:** 128 configurations with the walk (`-n 0–7` ×
  text/SAM × `--batch 1 -t 1` / `-t 22` × ±`--no-ref` × 100/32 bp, 20k reads), and 96 with
  the O(1) switch (`-n 2–7`).

**Timing** (search seconds, 500k reads, 500 Mb, `--lookup 11`, CPU 5, three rounds):

| run | bz | port, walk | port, O(1) switch |
|---|---:|---:|---:|
| 100 bp `-n 2` B = 1 | 4.22–4.43 | 4.27–4.37 | 4.21–4.27 |
| 100 bp `-n 2` B = 8 | 2.83–2.92 | 2.89–2.95 | 2.91–2.97 |
| 100 bp `-n 3` B = 1 | 5.86–6.37 | 5.85–6.23 | 5.74–6.04 |
| 100 bp `-n 3` B = 8 | 3.98–4.09 | 4.12–4.28 | 4.10–4.27 |
| 100 bp `-n 4` B = 1 | 7.75–7.96 | 7.98–8.14 | 7.80–8.02 |
| 100 bp `-n 4` B = 8 | 6.27–6.34 | 6.43–6.50 | 6.30–6.33 |
| 32 bp `-n 2` B = 1 | 11.74–11.87 | 12.11–12.24 | 12.00–12.14 |
| 32 bp `-n 2` B = 8 | 8.86–9.05 | 9.21–9.42 | 9.21–9.25 |
| 32 bp `-n 3` B = 1 | 35.0–35.2 | 35.8–36.7 | 34.4–35.3 |
| 32 bp `-n 3` B = 8 | 26.3–26.5 | 28.5–28.8 | 27.6–27.7 |
| 32 bp `-n 4` B = 1 | 258–318 | 277–379 | 251–285 |
| 32 bp `-n 4` B = 8 | 214–225 | 224–243 | 206–216 |

- **Carrying `Other_Start` (the walk column) costs 0–3% one read at a time and 2–8%
  batched.** An earlier run at `-n 0–3` agrees.
- **The O(1) switch** gains 1–4% at `-n 2–3` and 7–11% at 32 bp `-n 4`, where seeds
  switch most. There the port is 3–7% faster than `batman_zero_min`; at 100 bp it is level.
  The remaining gap is `-n 2–3` batched (up to 5%).
- **Compared with the hand-written searches** (`bi_mm`), the port's lanes are level with
  `batman_zero_min` batched, where `bi_mm -n 2` was 12–25% behind. The cost per step was
  the code, not the index.

### Locating from either row (`Both_Rows`; `--own-row` to compare)

A single-row range is one occurrence with a row in each index. `batman_zero_min` locates it
only when its own row is sampled. The port also takes the other row:
- **Where the path becomes single** (`Mm_Search_Enter`, `Mm_Single_Enter`), both rows are
  tested. The other row never moves on a single-row path, so after that only the moving own
  row is.
- **The other row needs to know where the string sits in the read.** Each range carries
  `Anchor`: the read position where its string starts (`revfmi` ranges) or ends (`fwfmi`).
  It is set where a search starts its segment or a switch makes a range. `fwfmi`'s SA value
  gives the read start `SA − Anchor`; `revfmi`'s gives `n − SA − Anchor`.
- **`-n 0` now loads `fwfmi` too**, with its SA samples, when the reference is on (its
  backward table is not built).

**Correctness:**
- **The verifying build caught a bug.** Segments no longer than the lookup key took their
  range straight from the table (`L.Cur = R`, also kept by `OP_SEARCH_KEEP`) without an
  `Anchor`. At 32 bp `-n 4`, 7,472 of 3.9M other-row locates were 16 bases off. After the
  fix, about 104M other-row locates (`-n 0–7`, 100/32 bp, `rx`) all equal the own row's walk.
- **SAM `-m 0`, records sorted, header byte-identical to `batman_zero_min`:** 27
  configurations (`-n 0–5` × 100/32 bp × `--batch 1 -t 1` / `-t 22`, 100 bp `-n 7`,
  `--no-ref`, `rx`). A path located earlier is found in a different order, so with `-m N`
  the first N hits may differ (any N hits are correct).

**Timing** (search seconds, O(1) switch in both port columns, CPU 5, three rounds; load
about 2–2.5, so rounds vary up to 10%):

| run | bz | port `--own-row` | port, both rows |
|---|---:|---:|---:|
| 100 bp `-n 0` B = 1 | 1.19–1.20 | 1.22–1.26 | 1.17–1.18 |
| 100 bp `-n 0` B = 8 | 0.63–0.65 | 0.66–0.70 | 0.65 |
| 100 bp `-n 1` B = 1 | 2.61–2.80 | 2.63–2.71 | 2.58–2.70 |
| 100 bp `-n 1` B = 8 | 1.50–1.58 | 1.55–1.61 | 1.50–1.64 |
| 100 bp `-n 2` B = 1 | 4.20–4.32 | 4.15–4.38 | 3.99–4.08 |
| 100 bp `-n 2` B = 8 | 2.71–2.76 | 2.85–3.14 | 2.68–2.82 |
| 100 bp `-n 3` B = 1 | 5.76–6.49 | 5.77–6.62 | 5.57–6.42 |
| 100 bp `-n 3` B = 8 | 4.47–4.84 | 4.62–4.97 | 4.54–4.75 |
| 100 bp `-n 4` B = 1 | 7.71–8.64 | 7.98–8.92 | 7.46–8.45 |
| 100 bp `-n 4` B = 8 | 6.62–7.23 | 6.90–7.38 | 6.61–7.60 |
| 32 bp `-n 0` B = 1 | 1.04–1.06 | 1.05–1.07 | 1.02–1.07 |
| 32 bp `-n 0` B = 8 | 0.53–0.56 | 0.58–0.62 | 0.58–0.59 |
| 32 bp `-n 1` B = 1 | 2.41–2.45 | 2.40 | 2.30–2.35 |
| 32 bp `-n 1` B = 8 | 1.51–1.53 | 1.60–1.72 | 1.56–1.57 |
| 32 bp `-n 2` B = 1 | 13.8–14.1 | 14.2–14.5 | 14.4–14.8 |
| 32 bp `-n 2` B = 8 | 9.0–10.6 | 10.0–11.3 | 9.6–11.5 |
| 32 bp `-n 3` B = 1 | 35.0–35.4 | 34.5–35.1 | 34.6–34.9 |
| 32 bp `-n 3` B = 8 | 26.6–29.0 | 27.8–31.3 | 28.4–30.7 |
| 32 bp `-n 4` B = 8 (one round) | 239 | 232 | 226 |

- **One read at a time, the other row gains 2–5%** at 100 bp `-n 0`, `-n 2` and `-n 4`
  and 32 bp `-n 1`. There the port is 2–5% faster than `batman_zero_min`, where it was
  level.
- **Batched, the gain is within noise.** With the walks' misses overlapped, ending a walk a
  few steps sooner saves little time.
- **32 bp `-n 2` is 1–3% slower with both rows** one read at a time. That is the one
  setting where the extra test does not pay.
- **Overall:** the port with both advantages is level with `batman_zero_min` or up to 5%
  faster one read at a time, 3–7% faster at 32 bp `-n 4`, and level to about 5% behind
  batched at 32 bp `-n 1–3`.

## Search schemes for `-n 4`'s 2|2 stage (tried, not adopted)

Could plans built for free direction changes (search schemes: Lam et al.'s 2BWT,
Kucherov et al. 2016) beat `batman_zero_min`'s half-based plans? The check was a
prototype first, with a go/no-go gate: at least 25% fewer steps before porting
(`tools/scheme_eval.cpp`, built by hand against `bi_index.hpp` / `bi_ref.hpp`).

**Where the time is.** At 32 bp, `-n 3` and `-n 4` are dominated by one kind of stage: seeds
of one half with mismatches, grown from a short exact piece (5–8 bases, thousands of
rows on 500 Mb), then crossed over the other half.
- **`-n 3`:** 300 of 440 steps a read are four explorations (an exact quarter, then the
  next quarter with exactly one mismatch). Any scheme set needs the same four, for the
  cases with no two adjacent exact quarters, so the ceiling was about 10%. Not
  prototyped.
- **`-n 4`:** 80% of 3,720 steps a read are the 2|2 stage (2,180 seed steps, 460
  crossing, 340 switch).

**The prototype**, for the 2|2 class:
- **Cases:** the read in eighths, every distribution with 2 mismatches in each half (100
  cases).
- **Search:** each case seeded from a run of exact eighths, extended toward the fewer
  mismatches first; cases sharing a seed searched once, a bitmask of the cases still
  possible at each node, located singletons from either row.
- **Correctness:** the same 2.2M hits as brute force on the toy index.

| 2\|2 stage, steps per read (10k reads, 500 Mb) | 32 bp | 100 bp |
|---|---:|---:|
| `batman_zero_min` (steps 48–59) | about 2,980 (2,640 without switches) | about 55 |
| prototype, each case's longest exact run (30 seeds) | 4,594 | 360 |
| prototype, seeds by greedy set cover (21 seeds) | 3,729 | 250 |

**Not adopted.** The prototype does more work than `batman_zero_min`, at both lengths.
- **Every seed costs a search per read.** At 100 bp that is 21–30 seeds × about 12 steps.
  `batman_zero_min`'s handful of seeds reuse the ranges of earlier stages
  (`OP_SEARCH_KEEP`, the chains) and are bounded by pass 1's failure depths.
- **Cases with no adjacent exact pair** still need 4–5-base seeds, as
  `batman_zero_min`'s thirds do.

The bidirectional index's gains stay the ones in `bi_batman`: the O(1) switch (7–11% at
32 bp `-n 4`) and locating from either row.

## Next steps

1. **Search schemes:** tried for `-n 4`'s 2|2 stage and not adopted (above). A
   scheme design would have to keep `batman_zero_min`'s reuse of earlier stages to compete.
2. **`bi_mm` and `bi_two_mismatch.hpp`** stay as references until the port has those.
