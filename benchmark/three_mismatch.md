# Three-mismatch search (`-n 3`)

> **Since the program refactor** (`five_mismatch.md`, "Program builder"): the hand-written
> per-read searches named below (`One_Mismatch_Search` … `Five_Mismatch_Search`,
> `Extend_Windows`, `Switch_Index`, the `Program_*` arrays, the `SEG_*` / `WIN_*` enums) are
> gone. `Mm_Builder` builds the same programs step for step, and `--batch 1` runs them
> on one lane of the engine.

## Summary

`batman_zero_min -n 3` reports hits with up to three mismatches. It exists both
per read (`Three_Mismatch_Search`) and batched (`Program_Three_Mismatch`,
`--batch B`), and the two give byte-identical output.

**Why quarters, not thirds:** three thirds can't guarantee an exact piece with
three mismatches (`two_mismatch.md`, section 6).
- The read's halves LH | RH are cut into quarters: Q1 = [0, c1), Q2 = [c1, LH),
  Q3 = [LH, c3), Q4 = [c3, L), with c1 = LH/2 and c3 = LH + (L − LH)/2.
- p is a mismatch in LH, q one in RH.

**What comes from earlier rounds:**
- Every stage but two either resumes an earlier round's rejects or uses a range
  the half scans kept.
- Only Q2 and Q3 are searched anew.
- The one-mismatch seeds use the failure bounds d and e
  (`two_mismatch.md`, section 4).

## Stages

"Count c" is the existing extension: a hit at exactly c mismatches, rejects
(deferred ranges) above c.

| # | stage | seed | extension | output (list) |
|---|---|---|---|---|
| 1 | 0\|0 | LH exact, revfmi; **keeps Q1's range** | right over RH, Count 0 | hits, rejects |
| 2 | 0\|1 | 0\|0 rejects | resume, Count 1 | hits, rejects (`TM_RIGHT1`) |
| 3 | 1\|0 | RH exact, fwfmi; **keeps Q4's range** | left over LH, Count 1 | hits, rejects (`TM_LEFT1`) |
| 4 | 0\|2 | 0\|1 rejects | resume, Count 2 | hits, rejects (`TM_RIGHT2`) |
| 5 | 2\|0 | 1\|0 rejects | resume, Count 2 | hits, rejects (`TM_LEFT2`) |
| 6a | LH seeds, p ∈ Q2 | **kept Q1**; skipped if d ≤ c1 | right over Q2, window {Q2: 1} | seeds (`TM_SEED_LR`, revfmi) |
| 6b | LH seeds, p ∈ Q1 | Q2 exact, fwfmi (**new search**) | left over Q1, window {Q1: 1} | seeds (`TM_SEED_LF`, fwfmi) |
| 7 | 1\|1 | the seeds; 6b's switched to revfmi (Start-only walk) unless located | right over RH, Count 2 | hits, rejects (`TM_CROSS_11`) |
| 8 | 0\|3 | 0\|2 rejects | resume, Count 3 | hits (`TM_RIGHT3`) |
| 9 | 3\|0 | 2\|0 rejects | resume, Count 3 | hits (`TM_LEFT3`) |
| 10 | 1\|2 | 1\|1 rejects | resume, Count 3 | hits (`TM_CROSS_12`) |
| 11a | RH seeds, q ∈ Q3 | **kept Q4**; skipped if e ≥ c3 | left over Q3, window {Q3: 1} | seeds (`TM_SEED_RF`, fwfmi) |
| 11b | RH seeds, q ∈ Q4 | Q3 exact, revfmi (**new search**) | right over Q4, window {Q4: 1} | seeds (`TM_SEED_RR`, revfmi) |
| 12 | 2\|1 | the seeds; 11b's switched to fwfmi unless located | left over LH, Count 3 (exactly two in LH) | hits (`TM_CROSS_21`) |

**Each hit is found once:**
- Stages 1–5, 8 and 9 have one half exact.
- 7 and 10 have exactly one mismatch in LH: the seeds have exactly one, and
  6a / 6b split on which quarter holds it.
- 12 has exactly one in RH and two in LH.

**Order:** the stages run, and report, in the order 0|1, 1|0, 0|2, 2|0, 1|1, 0|3,
3|0, 1|2, 2|1. So with `-m` the fewest mismatches come first, and the run stops
between stages at `MaxHits` (`Format_Three_Mismatch`, output count 3 for the
last four).

**Seed windows:** the one-mismatch seeds use the windowed extension with the
failure bounds.
- An LH window ends by d, an RH window starts at e or later.
- 6a is skipped when d ≤ c1 (Q1 then cannot be exact), and 11a when e ≥ c3.

**Lists:** the rolling lists of `two_mismatch.md`, section 2.
- The 0|k chain keeps its rejects in `RL_FWD`, the k|0 chain in `RL_BWD`, and
  the 1|k chain in `RL_CROSS`.
- 6a / 11a seeds go in `RL_SEED_1`, 6b / 11b in `RL_SEED_2`, cleared at each
  stage's end.
- That is seven lists per read, against 18 with one list per stage.
- **Under `-n 4`** (`four_mismatch.md`), 0|3, 3|0, 1|2 and 2|1 keep their
  rejects (2|1's in `RL_CROSS2`) for the four-mismatch round. `-n 3` drops them,
  as here.
- **The seed rule behind `-n 1..5`,** and a stage-by-stage comparison with
  BatMis's seed lengths: `five_mismatch.md`.

## Shared code

**Per read:**
- `Three_Mismatch_Search` reuses `One_Mismatch_Search` (with `Half_Scan` keeping
  Q1 and Q4), `Extend_Right` / `Extend_Left` (resumes and crosses),
  `Extend_Windows` (seeds) and `Switch_Index`.
- Its `Resume`, `Cross` and `Exact` lambdas are the stage shapes above.

**Batched:**
- `Program_Three_Mismatch` uses only existing ops (`OP_SEARCH_KEEP`,
  `OP_RESUME`, `OP_SEARCH` with window and kept slot, windowed `OP_EXTEND`,
  `OP_CROSS` with and without switch).
- The engine needed no new step kinds for `-n 3`.

## Correctness

- **Brute force** (`oracle3.py`):
  - **Method:** reported positions, by mismatch count, equal every genome
    position within Hamming distance 3, with no repeats.
  - **Genomes:** the 1 Mb random genome and the 49 kb repetitive toy.
  - **Reads:** 32, 36, 75 and 100 bp; 40 reads for each of the 35 placements
    of up to three mismatches over the quarters (`make_quarters.py`).
  - **Settings:** both reference modes, `--lookup` 3, 6, 11 and 14.
  - **Result:** 64 runs, every read matches. Batched (`--batch 16`) is
    byte-identical to per read in all of them.
- **500 Mb samples** (20k reads each of p32_mm2, p75_mm2, p100_mm2):
  - every reported position is at its reported distance, with no repeats;
  - batched (`--batch` 16 and 32) equals per read at `-m` 0, 1 and 3, both
    reference modes (36 runs).
- **`-n 0/1/2`:** byte-identical to the previous build (96 runs on 500 Mb).
- **`make test-located`:** now checks `-n 3` too: batched against per read, all
  lists field by field, and per read against a brute-force scan.
  - **Fault injection:** each fault is caught (20 failures each):
    - the 2|1 Count set to 2 for one seed list;
    - WIN_Q2 shifted by one base;
    - Q1 kept one base short.

## Memory accesses per read

`-n 3 -m 0 --lookup 11 --batch 16`, 20k reads, 500 Mb (`-DBATMAN_BATCH_STATS`):

| reads | total | of which the four seed extensions | `-n 2` total |
|---|---:|---:|---:|
| p32_mm2 | 434.1 | 316 | 130.2 |
| p75_mm2 | 65.8 | 18 | 41.4 |
| p100_mm2 | 68.0 | 14 | 44.2 |

- **Long reads:** `-n 3` costs about 24 more memory accesses per read than
  `-n 2`. The kept Q1 / Q4 ranges and the skipped 6a / 11a (44% of reads each)
  keep the seed cost low.
- **32 bp:** the four one-mismatch seed extensions are 73% of the work, about 79
  steps each. An 8-base quarter has thousands of rows on 500 Mb, the width
  problem quarters always had.
- **Not counted:** the exact searches of segments no longer than the lookup key
  (Q1–Q4 on 32 bp with `--lookup 11`) run inline in the program, not as counted
  steps. Their cost shows in the time, not in these counts.

## Timing

Input + search seconds, mean of 3 runs, `--lookup 11`, reference on, 500 Mb
genome, load average 1.1–1.8. The `-n 2` column is the same build, with kept
T1/T3.

| reads | `-m` | `-n 3`, batch 1 | `-n 3`, `--batch 16` | `-n 2`, `--batch 16` |
|---|---:|---:|---:|---:|
| p32_mm2 (1M) | 0 | 67.5 | 49.3 | 16.1 |
| p32_mm2 | 1 | 18.2 | 13.2 | 7.0 |
| p75_mm2 (500k) | 0 | 4.75 | 3.26 | 2.13 |
| p75_mm2 | 1 | 2.78 | 1.78 | 1.57 |
| p100_mm2 (500k) | 0 | 5.03 | 3.30 | 2.30 |
| p100_mm2 | 1 | 3.07 | 1.91 | 1.74 |

- **Long reads:** `-n 3` costs 1.4–1.5× `-n 2` at `-m 0`, and about 1.1× at
  `-m 1`.
- **32 bp:** `-n 3` costs 3× `-n 2` at `-m 0`, because of the wide 8-base
  quarter seeds (Memory accesses above). That is where `-n 3` should improve
  next.
