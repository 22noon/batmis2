# Five-mismatch search (`-n 5`), and the seed pattern behind `-n 1..5`

## Summary

`batman_zero_min -n 5` reports hits with up to five mismatches. Like every
mismatch program, it is built by one rule (`Mm_Builder`, "Program builder"
below) and runs on the engine for any `--batch B`, with byte-identical output.

- **One rule picks every seed.** From `-n 1` to `-n 5` (below), so `-n 5` needs
  just one new seed family: RH with exactly two mismatches, the mirror of
  `-n 4`'s LH family.
- **Seeds are never shorter than BatMis 3.00's,** stage by stage and for every
  placement of the mismatches.
- **`-n 0..4` output is unchanged** by adding `-n 5`.

## The seed pattern

A hit splits as i | j: i mismatches in LH, j in RH.
- **i ≤ j:** it grows from LH seeds with exactly i mismatches, extended right
  over RH. Round j = i crosses RH; each later round j + 1 resumes the previous
  round's rejects (a *chain*).
- **i > j:** the mirror. RH seeds with exactly j mismatches, extended left over
  LH.

**Seeds of a half with exactly s mismatches:**

| s | how | used for |
|---|---|---|
| 0 | the half itself, exact | 0\|k, k\|0 |
| 1 | the half's quarters: one is exact, the other holds the mismatch | 1\|k (k ≥ 1), k\|1 (k ≥ 2) |
| 2 | both in one quarter: the other quarter exact. One in each: the half's thirds, one exact (cases A, B, C1–C3) | 2\|k (k ≥ 2), k\|2 (k ≥ 3) |

**Cases by `-n`:**

| `-n` | new cases | seeds |
|---|---|---|
| 1 | 0\|1, 1\|0 | exact halves |
| 2 | 0\|2, 2\|0 (chains); 1\|1 | read thirds (`two_mismatch.md`) |
| 3 | 0\|3, 3\|0 (chains); 1\|1, 1\|2 (chain); 2\|1 | LH and RH s = 1 |
| 4 | 0\|4, 4\|0, 1\|3, 3\|1 (chains); 2\|2 | LH s = 2 (`four_mismatch.md`) |
| 5 | 0\|5, 5\|0, 1\|4, 4\|1, 2\|3 (chains); 3\|2 | RH s = 2 (new) |

**`-n 5` covers every split.** With i + j ≤ 5:
- from LH seeds: (0, 0..5), (1, 1..4), (2, 2..3);
- from RH seeds: (1..5, 0), (2..4, 1), (3, 2).

**`-n 6` and `-n 7`** add LH and RH s = 3 seeds, for 3|3 and 4|3
(`seven_mismatch.md`).

## Stages

The `-n 4` stages run first. Their five count-4 rounds now keep their rejects:
0|4 in `RL_FWD`, 4|0 in `RL_BWD`, 1|3 in `RL_CROSS`, 3|1 in `RL_CROSS2`, and 2|2
in the new `RL_CROSS22`. Then:

| # | stage | from | extension | list |
|---|---|---|---|---|
| 18 | 0\|5 | 0\|4 rejects | resume right over RH, Count 5 | `F` |
| 19 | 5\|0 | 4\|0 rejects | resume left over LH, Count 5 | `B` |
| 20 | 1\|4 | 1\|3 rejects | resume right over RH, Count 5 | `F` |
| 21 | 4\|1 | 3\|1 rejects | resume left over LH, Count 5 | `B` |
| 22 | 2\|3 | 2\|2 rejects | resume right over RH, Count 5 | `F` |
| 23 | 3\|2 | RH two-mismatch seeds | cross left over LH, Count 5 | `B` |

With `-m`, the run stops at `MaxHits` after the 18–22 group, as between
`-n 4`'s groups.

## The 3|2 seeds: RH's quarters first, then RH's thirds

This is the mirror of `four_mismatch.md`'s 2|2 seeds. Q3 = [LH, c3) and
Q4 = [c3, L) are RH's quarters. R1 = [LH, r1), R2 = [r1, r2) and R3 = [r2, L)
are its thirds, with r1 = LH + R/3, r2 = LH + 2R/3 and R = L − LH, so
r1 ≤ c3 ≤ r2. p < q are the two RH mismatches:

| case | M | seed (exact) | extension | windows |
|---|---|---|---|---|
| A′ | both in Q4 | Q3, revfmi | right over Q4 | {Q4: 2} |
| B′ | both in Q3 | Q4, fwfmi (the RH scan's kept range) | left over Q3 | {Q3: 2} |
| C1′ | q ∈ [c3, r2), p ∈ Q3 | R3, fwfmi | left over [LH, r2) | {[c3, r2): 1, Q3: 1} |
| C2′ | q ∈ R3, p ∈ [r1, c3) | R1, revfmi | right over [r1, L) | {R3: 1, [r1, c3): 1} |
| C3′ | q ∈ R3, p ∈ R1 | R2, fwfmi | left over R1 {R1: 1}, switch to revfmi, right over R3 {R3: 1} | |

- **Seed lists:** B′ and C1′ seeds (fwfmi) cross LH directly. The A′, C2′ and C3′
  seeds (revfmi) are switched to fwfmi first.
- **Failure bound:** read[RH_Fail, L) does not occur, so q ≥ RH_Fail.
  - Only a window that holds q alone is narrowed: C1′'s [c3, r2), and C2′'s and
    C3′'s R3.
  - B′ is skipped when RH_Fail ≥ c3, and C1′ when RH_Fail ≥ r2.

## Seed lengths against BatMis 3.00

BatMis's `-n 5` seeds k|2 the same way as its 2|2, mirrored onto RH:
- both mismatches in one RH quarter: the other quarter (`Possible_20`,
  `Possible_02`);
- one in each quarter: RHQR's halves, eighths of the read (`Search_XL01`,
  `Search_Backwards_XL10`).

Smallest seed per stage, in bases:

| stage | L = 100: BatMis | L = 100: port | L = 40: BatMis | L = 40: port |
|---|---|---|---|---|
| 0\|k, k\|0 | 50 | 50 | 20 | 20 |
| `-n 2` 1\|1 | 25 | 33 | 10 | 13 |
| 1\|k, k\|1 | 25 | 25 | 10 | 10 |
| 2\|2, 3\|2: both in one quarter | 25 | 25 | 10 | 10 |
| 2\|2, 3\|2: one in each quarter | 12 | 16 | 5 | 6 |

## Correctness

- **Brute force** (`oracle4.py`, SAM with `--forward-only`, K = 5):
  - **Method:** every genome position within Hamming distance 5, by distance.
    Also checked: no record repeated (same read, strand, sequence, position,
    MD), and NM equal to the mismatches in MD.
  - **Genomes:** the 1 Mb random genome and the 2 Mb repetitive genome.
  - **Reads:** 100 and 75 bp first, then 32 and 36 bp; 2 reads for each
    placement (as a multiset) of up to five mismatches over LH's parts H1,
    [h1, c1), [c1, h2), H3 and RH's R1, [r1, c3), [c3, r2), R3. That is every
    case A–C3 and A′–C3′ (1,287 placements).
  - **Settings:** both reference modes, `-m` 0, 1 and 3, `--batch` 1 and 16.
  - **Result:** 48 runs, every read matches, 0 duplicates, and batch 1 equals
    batch 16 byte for byte.
- **Fault injection:** each fault is caught:
  - case A′ dropped: 123 / 160 reads wrong;
  - C2′'s R3 window widened to [c3, L), so that it overlaps C1′: 44 / 147
    duplicates, which only the duplicate check sees.
- **Test suites:** `make test-located` (batched against per read, all lists;
  brute force) and `make test-sam` (SAM oracle, batch / `--no-ref` / BAM
  equivalence) now run up to `-n 5`. `make test-structured` passes.
- **Regression:** `-n 0..3` byte-identical (96 runs, 500 Mb); `-n 4` the same
  hit set for every read.

## Timing

Wall seconds including index load, against BatMis 3.00 (`batman` plus
`batdecode`): 500 Mb genome, 100k reads (25k each with 0, 1 in LH, 1 in RH and 2
mismatches), `--batch 16`, pinned to one core, load average 3.8–7.2. With
`-mall`, the hit sets are identical to the original's for every read: 230,195
hits at 100 bp, 247,188 at 75 bp and 251,909 at 32 bp.

| reads | `-m` | original | port |
|---|---|---:|---:|
| 100 bp | 1 | 7.5 | 3.2 |
| 100 bp | all | 24.8 | 12.5 |
| 75 bp | 1 | 6.3 | 2.9 |
| 75 bp | all | 38.0 | 16.9 |
| 32 bp | 1 | 4.9 | 3.7 |
| 32 bp | all | 618.1 | 469.6 |

- **Long reads:** 2.0–2.3× faster than the original.
- **32 bp, all hits:** 1.3× faster. As at `-n 4` (`four_mismatch.md`, Timing),
  the cost is in enumerating two-mismatch variants of a 16-base half, and in
  switching seeds between indexes.

## Program builder

The rule above is now the code. `Mm_Builder` (`src/batman_zero_min.cpp`) builds
each `-n k` program from it, and replaces the hand-written programs and per-read
searches.
- **`One_Mismatch`:** the half scans, 0|1 and 1|0.
- **`Two_Mismatch_Thirds`:** `-n 2`'s 1|1 by the read's thirds.
- **The rounds:** each chain resumed with Count c, then `Family(side, s, c)`,
  the round's new seeds crossing the other half.
- **`Seeds(side, s)`:** the `-n s` search on a half, kept to its exactly-s
  cases. One function serves both halves.

**Segments are fractions** of a part of the read (`Segment_Table`). In a half, a
fraction u counts from its outer end and sits at LH·u in LH, or LH + R·(1 − u)
in RH. So RH's cases are LH's with the same fractions. These give exactly the
earlier boundaries (c1, h1, h2, c3, r1, r2).

**Windows** (`Window_Table`) are segments with a count and one of three kinds:
- `WIN_FREE`: not bounded;
- `WIN_NARROW`: it holds the half's outer mismatch alone, so it is cut by
  `LH_Fail` / `RH_Fail`;
- `WIN_AFTER_SEED`: its seed is the exact outer part before it, so the case
  cannot hold past the failure bound.

**The same programs, step for step.** A dump of every program was identical
before and after:
- `-n 0` to `5`, plain and zig-zag;
- each step's segments resolved for 13 read lengths (12–150);
- each window set evaluated for every `LH_Fail` and `RH_Fail`.

`--batch 1` (text output) now runs the programs on one lane of the engine, and
the per-read searches are gone. `batman_zero_min.cpp` went from 4,240 to 3,275
lines.

**Checks:**
- Output byte-identical to the build before: 216 runs (`-n 0` to `5`, 32 / 75 /
  100 bp, `-m 0/1`, both reference modes, SAM and text at `--batch` 1 and 16).
- `make test-located`, `make test-sam` and `make test-structured` pass.
  `test_located_hits` now compares one lane against 16 (all lists) and checks
  `-n 1` to `5` against brute force.

**Speed:** search seconds, 500 Mb, 20k reads, `-m 0`, two runs each, load 5–8.
- **`--batch 16`:** unchanged within noise. At 32 bp: `-n 4` 21.4, 21.7 → 21.9,
  21.8; `-n 5` 42.4, 42.9 → 43.6, 43.2.
- **`--batch 1`**, where the hand-written per-read code ran before:

| reads | `-n` | per-read code | engine, one lane |
|---|---|---|---|
| 100 bp | 2–4 | equal within noise | |
| 100 bp | 5 | 1.07, 1.18, 1.31, 1.26 | 1.08, 1.30, 1.31, 1.38 |
| 32 bp | 4 | 19.2, 20.0, 21.0, 20.8 | 20.8, 21.4, 21.6, 21.9 |
| 32 bp | 5 | 39.1, 39.2, 46.8, 48.7 | 42.0, 43.9, 49.5, 49.7 |

- **The cost:** one lane of the engine is 3–9% slower than the per-read code at
  32 bp `-n 4/5`, and the same at long reads. It is not the prefetch: skipping
  it at B = 1 does not help.
- **A side finding:** at 32 bp `-n 4/5`, `--batch 16` is no faster than one lane.
  Those runs are dominated by wide windowed extensions, whose many branches per
  step already overlap their misses within one read.
