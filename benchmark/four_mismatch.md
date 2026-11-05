# Four-mismatch search (`-n 4`)

## Summary

`batman_zero_min -n 4` reports hits with up to four mismatches. Its program is
built by `Mm_Builder` (`five_mismatch.md`, "Program builder") and runs on the
engine for any `--batch B`; one lane and many give byte-identical output.

**What comes from `-n 3`:** four mismatches split over the halves LH | RH give
the new cases 0|4, 4|0, 1|3, 3|1 and 2|2.
- **0|4, 4|0, 1|3 and 3|1 are resumed rejects:** each continues a chain one
  round further. Under `-n 4` the 0|3, 3|0, 1|2 and 2|1 stages keep their rejects
  (`RL_FWD`, `RL_BWD`, `RL_CROSS`, and the new `RL_CROSS2`, the k|1 chain); under
  `-n 3` they drop them as before.
- **Only 2|2 needs new seeds:** LH with exactly two mismatches. Both in one
  quarter: the other quarter is the seed, as in BatMis. One in each quarter:
  LH's thirds, the `-n 2` scheme applied to LH, so every seed is at least as
  long as BatMis's for the same placement (section below).
- **`-n 0..3` output is unchanged:** 96 runs on the 500 Mb samples, byte for byte.

## Stages

The `-n 3` stages (`three_mismatch.md`) run first, then:

| # | stage | from | extension | list |
|---|---|---|---|---|
| 13 | 0\|4 | 0\|3 rejects | resume right over RH, Count 4 | `F` |
| 14 | 4\|0 | 3\|0 rejects | resume left over LH, Count 4 | `B` |
| 15 | 1\|3 | 1\|2 rejects | resume right over RH, Count 4 | `F` |
| 16 | 3\|1 | 2\|1 rejects | resume left over LH, Count 4 | `B` |
| 17 | 2\|2 | LH two-mismatch seeds | cross right over RH, Count 4 | `F` |

**Order:** with `-m`, the fewest mismatches come first. The run stops between
stage groups at `MaxHits`: the group 0|4, 4|0, 1|3, 3|1, then 2|2.

## The 2|2 seeds: quarters first, then LH's thirds

Q1 = [0, c1) and Q2 = [c1, LH) are LH's quarters; H1 = [0, h1), H2 = [h1, h2),
H3 = [h2, LH) are its thirds, with h1 = LH/3 and h2 = 2·LH/3 (`Sixths`), so
h1 ≤ c1 ≤ h2. p < q are the two LH mismatches. Each placement is in exactly one
case:

| case | M | seed (exact) | extension | windows |
|---|---|---|---|---|
| A | both in Q1 | Q2, fwfmi | left over Q1 | {Q1: 2} |
| B | both in Q2 | Q1, revfmi (the LH scan's kept range) | right over Q2 | {Q2: 2} |
| C1 | p ∈ [h1, c1), q ∈ Q2 | H1, revfmi | right over [h1, LH) | {[h1, c1): 1, Q2: 1} |
| C2 | p ∈ H1, q ∈ [c1, h2) | H3, fwfmi | left over [0, h2) | {H1: 1, [c1, h2): 1} |
| C3 | p ∈ H1, q ∈ H3 | H2, revfmi | right over H3 {H3: 1}, switch to fwfmi, left over H1 {H1: 1} | |

- **Complete and disjoint:** A and B take the placements with both in one
  quarter; C1–C3 the rest. C1 has p ≥ h1; C2 has p < h1 and q < h2; C3 has p < h1
  and q ≥ h2. Every position outside a case's windows must match.
- **Seed lists:** B and C1 seeds (revfmi) cross RH directly. The A, C2 and C3
  seeds (fwfmi) are first switched to revfmi. `Mm_Start_Switch` applies all of a
  seed's mismatches.
- **Failure bound:** read[0, LH_Fail) does not occur, so p < LH_Fail; q is not
  bounded.
  - Only a window that holds p alone is narrowed: C1's [h1, c1), and C2's and
    C3's H1.
  - B is skipped when LH_Fail ≤ c1, and C1 when LH_Fail ≤ h1.
- **Code:** `Mm_Builder::Seeds(PART_LH, 2)`, with window kinds `WIN_FREE`,
  `WIN_NARROW` and `WIN_AFTER_SEED` (`Window_Set`). The engine needed no new ops.

### Seed lengths against BatMis

BatMis 3.00 (`alternate/batman.cpp`) seeds 2|2 as follows:
- both mismatches in one LH quarter: from the other quarter (`Left_Mishits`,
  `Right_Mishits`);
- one in each quarter: from LHQL's halves, eighths of the read (`Search_01LX`,
  `Search_10LX`).

The port's seed is at least as long in every case:

| placement | L = 100: BatMis | L = 100: port | L = 40: BatMis | L = 40: port | L = 32: BatMis | L = 32: port |
|---|---|---|---|---|---|---|
| both in one quarter (A, B) | 25 | 25 | 10 | 10 | 8 | 8 |
| one in each (C1, C2, C3) | 12–13 | 16–17 | 5 | 6–7 | 4 | 5–6 |

**The first version was not always longer.** It used LH's thirds for every
placement (cases a, b1, b2, c). That gave a sixth (16 bases at 100 bp, 5 at
32 bp) where BatMis uses a quarter (25, or 8). Both mismatches in [25, 50) at
100 bp is one example.

## Each hit is found once

A hit's mismatch set S splits as i in LH and j in RH:
- **Stages are disjoint:**
  - 0|k and k|0 have the other half exact;
  - 1|k grows from exactly-one-in-LH seeds, k|1 (k ≥ 2) from exactly-one-in-RH
    seeds;
  - 2|2 grows from exactly-two-in-LH seeds, and cases A, B, C1, C2 and C3 are
    disjoint on M.
- **Within a stage, each string is reached once:**
  - an extension walks distinct strings depth first;
  - a resumed reject restarts at its level and skips the branch already taken
    (`Forbidden_Base`);
  - a chain list holds only its latest round.

**Checked empirically:** every run of the matrix below also checks that no
record repeats (same read, strand, sequence, position and MD), and that NM
equals the mismatches in MD.
- **0 duplicates** in every run.
- **The check is needed:** a fault that widens C2's H1 window to Q1, so that it
  overlaps C1, gives 31 / 120 duplicate records. The position-set comparison
  alone passes it.

## Correctness

- **Brute force** (`oracle4.py`, SAM with `--forward-only`):
  - **Method:** every genome position within Hamming distance 4, by distance;
    with `-m`, every record at its distance and the best distance present.
  - **Genomes:** the 1 Mb random genome, and a 2 Mb repetitive one (40 families
    at 0–3% divergence).
  - **Reads:** 75 and 100 bp first, then 32 and 36 bp; 3 reads for each
    placement (as a multiset) of up to four mismatches over LH's parts H1,
    [h1, c1), [c1, h2), H3 and RH's mirror parts (`make_parts.py`, 495
    placements).
  - **Settings:** both reference modes, `-m` 0, 1 and 3, `--batch` 1 and 16.
  - **Result:** 48 runs, every read matches, 0 duplicates, and batch 1 equals
    batch 16 byte for byte.
- **Fault injection:** each fault is caught:
  - C1's window shifted by one base: 6 / 6 reads wrong;
  - case A dropped: 90 / 120 reads wrong;
  - C2 overlapping C1: 31 / 120 duplicates.
- **Against the thirds-only version:** the same hit set for every read (20k
  reads per length at 32, 75 and 100 bp, 500 Mb, both reference modes). Only
  the order within the 2|2 stage can differ. `-n 0..3` stayed byte-identical
  (96 runs).
- **`make test-located`:** now also checks `-n 4`: batched (B = 1, 16) against
  per read, all lists field by field, and per read against brute force.
- **`make test-sam`:** the SAM oracle and the batch / `--no-ref` / BAM
  equivalence now include `-n 4`, with reads of up to four mismatches.
- **Against the original BatMis 3.00** (500 Mb, 100k mixed reads, `-mall`):
  identical hit sets for every read (216,526 hits at 100 bp, 239,258 at 75 bp,
  250,551 at 32 bp).
  With `-m 1` only the choice among equally good hits differs, as at `-n 3`.

## Timing

**Against the original** (wall seconds including index load): 500 Mb genome,
100k reads (25k each with 0, 1 in LH, 1 in RH and 2 mismatches), `--batch 16`,
pinned to one core, load average 4.8–7.2. The original is `batman` plus
`batdecode` (BatMis 3.00, built with `-fcommon` only).

| reads | `-m` | original | port |
|---|---|---:|---:|
| 100 bp | 1 | 7.7 | 3.1 |
| 100 bp | all | 19.5 | 9.0 |
| 75 bp | 1 | 6.3 | 3.1 |
| 75 bp | all | 25.7 | 11.9 |
| 32 bp | 1 | 5.0 | 3.6 |
| 32 bp | all | 345.9 | 235.7 |

**Quarters first against thirds only:** 20k reads, `-m 0`, search seconds, the
two builds back to back (`-DBATMAN_BATCH_STATS`).

| reads | lane steps per read | search |
|---|---|---|
| 100 bp | 283.7 → 283.5 | 1.44 → 1.48 |
| 75 bp | 411.4 → 352.2 (−14%) | 2.48 → 2.01 (−19%) |
| 32 bp | 8312 → 8198 (−1%) | 46.9 → 46.5 |

- **Long reads:** 2.1–2.5× faster than the original. At 75 bp the quarter seeds
  save a fifth of the search.
- **32 bp, all hits:** 1.5× faster than the original (1.44 GB peak), but the new
  partition barely helps there. The 2|2 stage is about 3,500 of 4,100 lane steps
  per strand:
  - **Windowed seed extensions:** about 2,300 steps. On a 16-base LH, every
    two-mismatch variant that occurs in the genome is listed, whichever part is
    the seed.
  - **Index switches of the fwfmi seeds:** about 700 steps. Each one walks LH's
    16 bases, one rank per base.
- **Next step for short reads:** a single-row fwfmi seed could be located (LF to
  a sampled row, about 4 steps on average) and its RH compared with the
  reference, instead of switching it to revfmi and extending.
