# Six and seven mismatches (`-n 6`, `-n 7`): a proof of concept

## Summary

`batman_zero_min -n 6` and `-n 7` report hits with up to six or seven mismatches. They
are built by the same rule as `-n 1..5` (`five_mismatch.md`, "The seed pattern"). The
rule needed one new piece: seeds of a half with exactly three mismatches.

- **Correct on every check:**
  - **Brute force:** all hits at `-m 0`, the best hits at `-m 1` and `-m 3`, no
    duplicate record, NM equal to the mismatches in MD.
  - **Read sets:** 100 bp and 32 bp, every split i|j up to 7. The read sets place the
    three seed-side mismatches over every combination of the half's parts.
  - **Genomes:** a random genome and a repeat-rich one.
  - **Engine and threads:** `--batch 1`, `--batch 16` and `-t 4` give identical output.
- **`-n 0..5` output is unchanged**, byte for byte.
- **Cost at 100 bp** (reads with 0–7 mismatches): per read, `-n 6` takes 1.7× and `-n 7`
  2.5× the time of `-n 5`.
- **Smallest seed:** an eighth of the read, 12 bases at 100 bp and 4 at 32 bp. With three
  mismatches in a half, no seed choice inside the half can do better.

This is a proof of concept. It is tested for accuracy at 100 bp and 32 bp, and timed at
100 bp only.

## What `-n 6` and `-n 7` add

A hit splits as i | j, i mismatches in LH and j in RH:
- **i ≤ j:** grows from LH seeds with exactly i mismatches.
- **i > j:** grows from RH seeds with exactly j mismatches.
- **Later rounds** resume the previous round's rejects (chains).

| `-n` | new cases | from |
|---|---|---|
| 6 | 0\|6, 6\|0, 1\|5, 5\|1, 2\|4, 4\|2 | chains |
| 6 | **3\|3** | LH seeds, s = 3 (new) |
| 7 | 0\|7, 7\|0, 1\|6, 6\|1, 2\|5, 5\|2, 3\|4 | chains (3\|4 resumes 3\|3's rejects) |
| 7 | **4\|3** | RH seeds, s = 3 (the mirror) |

`Mm_Builder::Build`'s round loop already called `Family(LH, 3, 6)` and `Family(RH, 3, 7)`,
so `-n 6/7` needed:
- **`Mm_Builder::Three_Seeds`**, the s = 3 seeds (below);
- **two more chain lists**: `RL_CHAINS_END = RL_FWD + 7`. `-n 7` keeps 0|·, ·|0, 1|·,
  ·|1, 2|·, ·|2 and 3|·;
- **`RL_SEED_4`**, a second list for seeds extended in two steps. `OP_CROSS` does not
  consume its input list, so the two switch cases cannot share `RL_SEED_3`;
- **limits:**
  - `MAX_MISMATCH_SEARCH = 7`;
  - the programs and zig-zag programs built up to 7;
  - `MM_PROGRAM_MAX = 320` (the `-n 7` zig-zag program has 307 steps);
  - `-n 0..7` on the command line.

The engine needed no new ops or window kinds. Each read now has 13 result lists instead
of 10.

| `-n` | steps | zig-zag steps | stages |
|---|---:|---:|---:|
| 5 | 88 | 179 | 20 |
| 6 | 119 | 241 | 27 |
| 7 | 152 | 307 | 35 |

## The s = 3 seeds: the `-n 3` search on the half

The pattern so far is that a half's s-mismatch seeds are the `-n s` search run on that
half. The half is split, counting from its outer end:
- **quarters:** Q_out = [0, ½) and Q_in = [½, 1), the read's quarters;
- **eighths:** E1 = [0, ¼), E2 = [¼, ½), E3 = [½, ¾) and E4 = [¾, 1).

p is the outer mismatch (the first in LH, the last in RH), and r the inner one.

| (Q_out, Q_in) | case | seed (exact) | extension | windows |
|---|---|---|---|---|
| (0, 3) | | Q_out (the half scan's kept range) | inward over Q_in | {Q_in: 3}, skipped once the half's exact scan fails inside Q_out |
| (3, 0) | | Q_in | outward over Q_out | {Q_out: 3} |
| (1, 2) | p ∈ E2 | E1 | inward over [¼, 1) | {E2: 1 narrowed, Q_in: 2} |
| (1, 2) | p ∈ E1 | E2 | outward over E1 {E1: 1 narrowed}, switch, inward over Q_in {Q_in: 2} | |
| (2, 1) | r ∈ E3 | E4 | outward over [0, ¾) | {Q_out: 2, E3: 1} |
| (2, 1) | r ∈ E4 | E3 | inward over E4 {E4: 1}, switch, outward over Q_out {Q_out: 2} | |

- **Complete and disjoint:** the six cases split every placement of exactly three
  mismatches in the half.
- **Cheaper part first:** the two-step cases extend over the part with one mismatch
  first, then switch indexes (`OP_CROSS` with a switch) for the part with two.
- **Failure bound:** read[0, LH_Fail) does not occur, so p < LH_Fail. Only a window that
  holds p alone is narrowed: E2 in the first (1, 2) case, E1 in the second.
- **RH:** the same cases mirrored (`Seeds(PART_RH, 3)`), for 4|3.

## Smallest seed in the worst case

**No seed inside a half can be longer.** Three mismatches can sit at a half's quarter
points, leaving no exact stretch longer than a quarter of the half, which is L/8. Five of
the six cases reach that bound; (0, 3) and (3, 0) have a quarter of the read.

Smallest seed per stage, in bases (fractions floored as `Part_Position` does):

| stage | L = 100 | L = 75 | L = 40 | L = 32 |
|---|---:|---:|---:|---:|
| 0\|k, k\|0: the half | 50 | 37 | 20 | 16 |
| 1\|k, k\|1: a quarter | 25 | 18 | 10 | 8 |
| 2\|k, k\|2, both in one quarter | 25 | 18 | 10 | 8 |
| 2\|k, k\|2, one in each: a sixth (`-n 4/5`'s minimum) | 16 | 12 | 6 | 5 |
| 3\|k, k\|3, all three in one quarter | 25 | 18 | 10 | 8 |
| **3\|k, k\|3, otherwise: an eighth** | **12** | **9** | **5** | **4** |

- **`-n 7` has the same minimum as `-n 6`.** Its new seeded stage, 4|3, uses RH's s = 3
  seeds; the rest are chains.
- **Crossing the halves could do a little better.** Six mismatches anywhere in the read
  leave an exact stretch of at least L/7 (14 at 100 bp). That would need a read-level
  case split, like `-n 2`'s thirds for 1|1, and it gains nothing at 32 bp.
- **BatMis 3.00 stops at `-n 5`**, so there is nothing to compare with.
- **What it means at 32 bp:** a 4-base seed matches a large share of the index's rows.
  At 32 bp, `-n 6/7` are for accuracy checks rather than throughput.

## Correctness

**Targeted read sets** (`scratchpad/n7/make_parts7.py`):
- every split i|j with i + j ≤ 7, 6 reads each;
- for 3|3, 3|4 (LH) and 4|3 (RH), every placement of the three seed-side mismatches over
  the half's six parts (cut at ¼, ⅓, ½, ⅔ and ¾), 2 reads each;
- 552 reads at 100 bp and 504 at 32 bp;
- on `g1m` (1 Mb random genome) and `rep2m` (2 Mb with repeats).

**Matrix** (`scratchpad/n7/matrix7.sh`): `-n 5, 6, 7` × `-m 0, 1, 3` × with and without
`--no-ref` × 100 bp and 32 bp × both genomes. That is 72 combinations, each run with
`--batch 1`, `--batch 16` and `-t 4`.
- **All 72 pass:** identical output, no wrong or missing hit against brute force
  (`oracle4.py`), no duplicate record and no bad NM.

**Unchanged:** `-n 0..5` output on 20k 100 bp and 5k 32 bp reads is byte-identical to
the build before this change (40 runs: text and SAM, `-t 1` and `-t 4`).

**Test suites:** `make test-located`, `make test-sam` and `make test-structured` pass.
- `test-located` compares `--batch 1` with 16, and checks every hit against brute force,
  for k = 1..7. Its reads now include reads with 3–7 mismatches.
- `test-sam`'s reads have up to 7 planted mismatches, and its loops run `-n 0..7`.
- `test-located` now takes about 9 minutes and `test-sam` about 2.5 minutes, mostly in
  `-n 6/7` on the 12- and 20-base reads.

## Timing (100 bp)

`mm07_100.fq` has 100k reads from the 500 Mb random genome: 0–7 planted mismatches
(about 11k reads each) plus 10% random reads, half reverse-complemented. Text output,
search seconds, 2 runs each. The load average from other users was 6–9, so ±25% at
`-t 22`.

| | `-n 5` | `-n 6` | `-n 7` |
|---|---:|---:|---:|
| `-t 1`, 20k reads, `-m 0` | 1.52 s (76 µs/read) | 2.56 s (128 µs) | 3.85 s (192 µs) |
| `-t 22`, 100k reads, `-m 0` | 0.37–0.61 s | 0.97–1.09 s | 1.44–1.71 s |
| `-t 22`, 100k reads, `-m 1` | 0.43–0.56 s | 0.57–0.72 s | 0.60–0.92 s |
| hits reported (`-m 0`) | 66,687 | 80,597 | 94,410 |

- **`-n 6` costs 1.7× `-n 5` per read, and `-n 7` 2.5×** (`-t 1`).
- **Correct mismatch counts:** in SAM on the first 1,000 reads, each planted read's best
  hit has its planted count, or one fewer where a substitution happens to match another
  location. The random reads stay unmapped.
- **On a set of mostly exact reads** (`shuf100.fq`: 0–2 mismatches), `-n 7` adds little:
  at `-m 1` all three counts take 0.11–0.17 s, since the run stops at the first hits.

32 bp is checked for accuracy only (above); it was not timed.
