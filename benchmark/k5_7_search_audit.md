# Redundant work at K = 5–7: an audit of the search program

## Summary

- **No repeated hits, and no repeated exact walks remain.** The cases are disjoint
  (`seven_mismatch.md`). Since the shared seed walks (`shared_walks.md`), each exact search
  that still runs has its own origin and direction, and no served search fell back to
  searching: none of them ran index steps of its own on 100,000 reads.
- **One repeat is left in the program itself: the same seed extended over the same
  segment in two or three rounds.** The 1|k, 2|k and 3|k families extend the half's two
  quarter seeds with exactly 1, 2 and 3 mismatches, in rounds 2, 4 and 6 (and the
  mirrors in 3, 5 and 7). Those 8 to 12 extensions take 8.7–9.9% of the all-hits search
  and 4.2–4.7% at `-m 1`. Doing each seed's extension once would save at most about 5% (all hits)
  and 1.5% (`-m 1`).
- **The same locus is compared with the genome again and again, but that costs little.**
  At `-n 7`, 9 of a read's 14 genome compares are of a read strand and locus already
  compared. The repeats cost about 2% of the search, because the first compare has
  already brought the genome words into cache.
- **The chains and crossings do no rework.** Resumed ranges continue where the previous
  round deferred them. All of them together, with the 0|0 extensions, take 9–10% of the search.
- **The cost is the number of distinct seeds.** The 6, 8 or 10 exact searches a strand
  that still run take 54–56% of the all-hits search and 67–68% at `-m 1`. Each costs
  about the same whatever its length (750–1,000 cycles a strand on E. coli). These
  searches are not repeats. They can only be removed with a scheme that needs fewer
  distinct seeds.

## Method

- **Program:** `-n 5`, `-n 6` and `-n 7` dumped from `Program_Of` (a scratch tool that
  includes `batmis1.cpp`). Segments are given at 150 and 100 bp. The `-n 5` and `-n 6`
  programs are steps 0–86 and 0–117 of the `-n 7` program (plus the end), the same
  except for which lists keep rejects. So one table covers all three.
- **Cost:** batmis1 with `-DBATMAN_BATCH_STATS` and rdtsc around each lane step. The run
  maps 100,000 simulated 150 bp reads on E. coli (0–7 substitutions, 10% random), SAM,
  both strands, on one P-core. Each step's cost is the sum over both strands, as a share of
  all the search's cycles, including the prefetch pass.
- **Repeated loci:** a second scratch build times every genome compare
  (`Genome_Mismatches`, `First_Mismatch`). It keys each compare by the read strand's
  packed bases and the locus.

## Partitions

The read splits into halves LH = [0, L/2) and RH = [L/2, L). Each half's parts are
fractions counted from its outer end (read start for LH, read end for RH). Positions are
floored as `Part_Position` does.

| part of a half | fraction from the outer end | LH at 150 bp | LH at 100 bp | RH at 150 bp | RH at 100 bp | used by |
|---|---|---|---|---|---|---|
| the half | 0 – 1 | [0,75) | [0,50) | [75,150) | [50,100) | 0\|j, j\|0 walks; the other half of every crossing |
| Q_out | 0 – ½ | [0,37) | [0,25) | [112,150) | [75,100) | seed for s = 1, B and (0,3); extension for s = 1, A and (3,0) |
| Q_in | ½ – 1 | [37,75) | [25,50) | [75,112) | [50,75) | seed for s = 1, A and (3,0); extension for s = 1, B and (0,3) |
| H_out | 0 – ⅓ | [0,25) | [0,16) | [125,150) | [83,100) | seed C1; narrowed window for C2 and C3 |
| H_mid | ⅓ – ⅔ | [25,50) | [16,33) | [100,125) | [66,83) | seed C3 |
| H_in | ⅔ – 1 | [50,75) | [33,50) | [75,100) | [50,66) | seed C2; extension C3 |
| [⅓, ½), [½, ⅔) | | [25,37), [37,50) | [16,25), [25,33) | [112,125), [100,112) | [75,83), [66,75) | windows of C1 and C2 |
| E1 | 0 – ¼ | [0,18) | [0,12) | [131,150) | [87,100) | seed (1,2) p ∈ E2; extension (1,2) p ∈ E1 |
| E2 | ¼ – ½ | [18,37) | [12,25) | [112,131) | [75,87) | seed (1,2) p ∈ E1 |
| E3 | ½ – ¾ | [37,56) | [25,37) | [93,112) | [62,75) | seed (2,1) r ∈ E4 |
| E4 | ¾ – 1 | [56,75) | [37,50) | [75,93) | [50,62) | seed (2,1) r ∈ E3; extension (2,1) r ∈ E4 |

The cases, by how many mismatches the seed half has (s):

| s | rounds | cases of the seed half | smallest seed (150 / 100 bp) |
|---|---|---|---|
| 0 | 1 (exact halves), then chains | the half exact | 75 / 50 |
| 1 | 2 (LH), 3 (RH) | Q_out exact + 1 in Q_in; Q_in exact + 1 in Q_out | 37 / 25 |
| 2 | 4 (LH), 5 (RH) | A (both in Q_out), B (both in Q_in), C1, C2, C3 (thirds) | 25 / 16 |
| 3 | 6 (LH), 7 (RH) | (0,3), (3,0), (1,2) ×2, (2,1) ×2 (eighths) | 18 / 12 |

The other half takes the rest of the budget in the family's crossing (round C = 2s or
2s + 1, at most C mismatches in all), then in later rounds as a chain of resumed rejects.
`-n 5` ends after round 5, `-n 6` after round 6.

## Every scan and extension (`-n 7`; `-n 5` is steps 0–86, `-n 6` steps 0–117)

- **Steps:** step numbers are those of the single-strand program. Stop checks and stage
  ends are left out.
- **Windows:** "exactly N in" a segment.
  - *narrow*: the window is cut at the exact half's failure depth.
  - *after-seed*: the case is skipped once that depth falls inside its seed.
- **Exact seed obtained:**
  - *walk*: searched, and saves the states of the later searches it serves.
  - *served*: taken from a walk, with no index steps.
  - *searched*: a search of its own.
- **Cost:** share of all `-n 7` search cycles at 150 bp, with all hits and with `-m 1`. A
  served search costs 0 and its step is charged nothing. Windowed searches the windows
  skip cost nothing either.

| step | round | case | op | segment (150 / 100 bp) | windows: exactly N in | exact seed obtained | all % | -m 1 % |
|---:|---|---|---|---|---|---|---:|---:|
| 0 | 1 | 0\|j: LH exact | exact → | [0,75) / [0,50) |  | walk; keeps [0,37); serves [0,37), [0,25), [0,18) | 5.7 | 14.1 |
| 1 | 1 | 0\|0: RH exact from the LH walk | extend → | [75,150) / [50,100) |  |  | 0.5 | 1.0 |
| 3 | 1 | 0\|1: the 0\|0 rejects, resumed | resume →, ≤ 1 in all | [75,150) / [50,100) |  |  | 0.2 | 0.5 |
| 6 | 1 | j\|0: RH exact | exact ← | [75,150) / [50,100) |  | walk; keeps [112,150); serves [112,150), [125,150), [131,150) | 6.1 | 12.9 |
| 7 | 1 | 1\|0: LH with one mismatch | extend ← | [0,75) / [0,50) |  |  | 0.6 | 0.5 |
| 10 | 2 | chain → RH, list 2 | resume →, ≤ 2 in all | [75,150) / [50,100) |  |  | 0.1 | 0.2 |
| 12 | 2 | chain ← LH, list 3 | resume ←, ≤ 2 in all | [0,75) / [0,50) |  |  | 0.1 | 0.2 |
| 15 | 2 | LH seeds, exactly 1: Q_out seed (mismatch in Q_in) | exact → | [0,37) / [0,25) | [37,75)/[25,50) = 1 (narrow) | served by step 0's walk | 0.0 | 0.0 |
| 16 | 2 |  | extend → | [37,75) / [25,50) | [37,75)/[25,50) = 1 (narrow) |  | 0.7 | 0.8 |
| 17 | 2 | LH seeds, exactly 1: Q_in seed (mismatch in Q_out) | exact ← | [37,75) / [25,50) |  | walk; serves [37,75), [50,75), [56,75) | 5.6 | 10.5 |
| 18 | 2 |  | extend ← | [0,37) / [0,25) | [0,37)/[0,25) = 1 (narrow) |  | 0.9 | 0.9 |
| 19 | 2 | cross → RH (seeds list 9) | cross →, ≤ 2 in all | [75,150) / [50,100) |  |  | 0.3 | 0.5 |
| 20 | 2 | cross → RH (seeds list 10) | cross → (switch), ≤ 2 in all | [75,150) / [50,100) |  |  | 0.3 | 0.5 |
| 23 | 3 | chain → RH, list 2 | resume →, ≤ 3 in all | [75,150) / [50,100) |  |  | 0.1 | 0.1 |
| 25 | 3 | chain ← LH, list 3 | resume ←, ≤ 3 in all | [0,75) / [0,50) |  |  | 0.1 | 0.1 |
| 27 | 3 | chain → RH, list 4 | resume →, ≤ 3 in all | [75,150) / [50,100) |  |  | 0.2 | 0.5 |
| 30 | 3 | RH seeds, exactly 1: Q_out seed (mismatch in Q_in) | exact ← | [112,150) / [75,100) | [75,112)/[50,75) = 1 (narrow) | served by step 6's walk | 0.0 | 0.0 |
| 31 | 3 |  | extend ← | [75,112) / [50,75) | [75,112)/[50,75) = 1 (narrow) |  | 0.7 | 0.5 |
| 32 | 3 | RH seeds, exactly 1: Q_in seed (mismatch in Q_out) | exact → | [75,112) / [50,75) |  | walk; serves [75,112), [75,100), [75,93) | 5.8 | 8.9 |
| 33 | 3 |  | extend → | [112,150) / [75,100) | [112,150)/[75,100) = 1 (narrow) |  | 0.9 | 0.6 |
| 34 | 3 | cross ← LH (seeds list 9) | cross ←, ≤ 3 in all | [0,75) / [0,50) |  |  | 0.3 | 0.3 |
| 35 | 3 | cross ← LH (seeds list 10) | cross ← (switch), ≤ 3 in all | [0,75) / [0,50) |  |  | 0.4 | 0.3 |
| 38 | 4 | chain → RH, list 2 | resume →, ≤ 4 in all | [75,150) / [50,100) |  |  | 0.0 | 0.1 |
| 40 | 4 | chain ← LH, list 3 | resume ←, ≤ 4 in all | [0,75) / [0,50) |  |  | 0.0 | 0.1 |
| 42 | 4 | chain → RH, list 4 | resume →, ≤ 4 in all | [75,150) / [50,100) |  |  | 0.1 | 0.3 |
| 44 | 4 | chain ← LH, list 5 | resume ←, ≤ 4 in all | [0,75) / [0,50) |  |  | 0.2 | 0.3 |
| 47 | 4 | LH seeds, exactly 2: A: both in Q_out | exact ← | [37,75) / [25,50) |  | served by step 17's walk | 0.0 | 0.0 |
| 48 | 4 |  | extend ← | [0,37) / [0,25) | [0,37)/[0,25) = 2 |  | 0.7 | 0.4 |
| 49 | 4 | LH seeds, exactly 2: B: both in Q_in | exact → | [0,37) / [0,25) | [37,75)/[25,50) = 2 (after-seed) | served by step 0's walk | 0.0 | 0.0 |
| 50 | 4 |  | extend → | [37,75) / [25,50) | [37,75)/[25,50) = 2 (after-seed) |  | 0.7 | 0.4 |
| 51 | 4 | LH seeds, exactly 2: C1: p in [1/3,1/2) | exact → | [0,25) / [0,16) | [25,37)/[16,25) = 1 (narrow), [37,75)/[25,50) = 1 | served by step 0's walk | 0.0 | 0.0 |
| 52 | 4 |  | extend → | [25,75) / [16,50) | [25,37)/[16,25) = 1 (narrow), [37,75)/[25,50) = 1 |  | 0.8 | 0.7 |
| 53 | 4 | LH seeds, exactly 2: C2: p in H_out, other in [1/2,2/3) | exact ← | [50,75) / [33,50) | [0,25)/[0,16) = 1 (narrow), [37,50)/[25,33) = 1 | served by step 17's walk | 0.0 | 0.0 |
| 54 | 4 |  | extend ← | [0,50) / [0,33) | [0,25)/[0,16) = 1 (narrow), [37,50)/[25,33) = 1 |  | 0.8 | 0.6 |
| 55 | 4 | LH seeds, exactly 2: C3: p in H_out, other in H_in | exact → | [25,50) / [16,33) | [0,25)/[0,16) = 1 (narrow) | searched | 5.4 | 6.9 |
| 56 | 4 |  | extend → | [50,75) / [33,50) | [50,75)/[33,50) = 1 |  | 1.1 | 0.7 |
| 57 | 4 |  | cross ← (switch) | [0,25) / [0,16) | [0,25)/[0,16) = 1 (narrow) |  | 0.4 | 0.3 |
| 58 | 4 | cross → RH (seeds list 9) | cross →, ≤ 4 in all | [75,150) / [50,100) |  |  | 0.2 | 0.3 |
| 59 | 4 | cross → RH (seeds list 10) | cross → (switch), ≤ 4 in all | [75,150) / [50,100) |  |  | 0.4 | 0.5 |
| 62 | 5 | chain → RH, list 2 | resume →, ≤ 5 in all | [75,150) / [50,100) |  |  | 0.0 | 0.0 |
| 64 | 5 | chain ← LH, list 3 | resume ←, ≤ 5 in all | [0,75) / [0,50) |  |  | 0.0 | 0.0 |
| 66 | 5 | chain → RH, list 4 | resume →, ≤ 5 in all | [75,150) / [50,100) |  |  | 0.1 | 0.2 |
| 68 | 5 | chain ← LH, list 5 | resume ←, ≤ 5 in all | [0,75) / [0,50) |  |  | 0.1 | 0.2 |
| 70 | 5 | chain → RH, list 6 | resume →, ≤ 5 in all | [75,150) / [50,100) |  |  | 0.2 | 0.4 |
| 73 | 5 | RH seeds, exactly 2: A: both in Q_out | exact → | [75,112) / [50,75) |  | served by step 32's walk | 0.0 | 0.0 |
| 74 | 5 |  | extend → | [112,150) / [75,100) | [112,150)/[75,100) = 2 |  | 0.7 | 0.2 |
| 75 | 5 | RH seeds, exactly 2: B: both in Q_in | exact ← | [112,150) / [75,100) | [75,112)/[50,75) = 2 (after-seed) | served by step 6's walk | 0.0 | 0.0 |
| 76 | 5 |  | extend ← | [75,112) / [50,75) | [75,112)/[50,75) = 2 (after-seed) |  | 0.6 | 0.2 |
| 77 | 5 | RH seeds, exactly 2: C1: p in [1/3,1/2) | exact ← | [125,150) / [83,100) | [112,125)/[75,83) = 1 (narrow), [75,112)/[50,75) = 1 | served by step 6's walk | 0.0 | 0.0 |
| 78 | 5 |  | extend ← | [75,125) / [50,83) | [112,125)/[75,83) = 1 (narrow), [75,112)/[50,75) = 1 |  | 0.8 | 0.4 |
| 79 | 5 | RH seeds, exactly 2: C2: p in H_out, other in [1/2,2/3) | exact → | [75,100) / [50,66) | [125,150)/[83,100) = 1 (narrow), [100,112)/[66,75) = 1 | served by step 32's walk | 0.0 | 0.0 |
| 80 | 5 |  | extend → | [100,150) / [66,100) | [125,150)/[83,100) = 1 (narrow), [100,112)/[66,75) = 1 |  | 0.8 | 0.4 |
| 81 | 5 | RH seeds, exactly 2: C3: p in H_out, other in H_in | exact ← | [100,125) / [66,83) | [125,150)/[83,100) = 1 (narrow) | searched | 5.7 | 5.4 |
| 82 | 5 |  | extend ← | [75,100) / [50,66) | [75,100)/[50,66) = 1 |  | 1.1 | 0.5 |
| 83 | 5 |  | cross → (switch) | [125,150) / [83,100) | [125,150)/[83,100) = 1 (narrow) |  | 0.4 | 0.2 |
| 84 | 5 | cross ← LH (seeds list 9) | cross ←, ≤ 5 in all | [0,75) / [0,50) |  |  | 0.2 | 0.2 |
| 85 | 5 | cross ← LH (seeds list 10) | cross ← (switch), ≤ 5 in all | [0,75) / [0,50) |  |  | 0.4 | 0.3 |
| 88 | 6 | chain → RH, list 2 | resume →, ≤ 6 in all | [75,150) / [50,100) |  |  | 0.0 | 0.0 |
| 90 | 6 | chain ← LH, list 3 | resume ←, ≤ 6 in all | [0,75) / [0,50) |  |  | 0.0 | 0.0 |
| 92 | 6 | chain → RH, list 4 | resume →, ≤ 6 in all | [75,150) / [50,100) |  |  | 0.0 | 0.1 |
| 94 | 6 | chain ← LH, list 5 | resume ←, ≤ 6 in all | [0,75) / [0,50) |  |  | 0.0 | 0.1 |
| 96 | 6 | chain → RH, list 6 | resume →, ≤ 6 in all | [75,150) / [50,100) |  |  | 0.1 | 0.2 |
| 98 | 6 | chain ← LH, list 7 | resume ←, ≤ 6 in all | [0,75) / [0,50) |  |  | 0.1 | 0.2 |
| 101 | 6 | LH seeds, exactly 3: (0,3) | exact → | [0,37) / [0,25) | [37,75)/[25,50) = 3 (after-seed) | served by step 0's walk | 0.0 | 0.0 |
| 102 | 6 |  | extend → | [37,75) / [25,50) | [37,75)/[25,50) = 3 (after-seed) |  | 0.7 | 0.1 |
| 103 | 6 | LH seeds, exactly 3: (3,0) | exact ← | [37,75) / [25,50) |  | served by step 17's walk | 0.0 | 0.0 |
| 104 | 6 |  | extend ← | [0,37) / [0,25) | [0,37)/[0,25) = 3 |  | 0.7 | 0.1 |
| 105 | 6 | LH seeds, exactly 3: (1,2), p in E2 | exact → | [0,18) / [0,12) | [18,37)/[12,25) = 1 (narrow), [37,75)/[25,50) = 2 | served by step 0's walk | 0.0 | 0.0 |
| 106 | 6 |  | extend → | [18,75) / [12,50) | [18,37)/[12,25) = 1 (narrow), [37,75)/[25,50) = 2 |  | 1.1 | 0.3 |
| 107 | 6 | LH seeds, exactly 3: (1,2), p in E1 | exact ← | [18,37) / [12,25) | [0,18)/[0,12) = 1 (narrow) | searched | 4.5 | 2.8 |
| 108 | 6 |  | extend ← | [0,18) / [0,12) | [0,18)/[0,12) = 1 (narrow) |  | 0.9 | 0.2 |
| 109 | 6 |  | cross → (switch) | [37,75) / [25,50) | [37,75)/[25,50) = 2 |  | 0.3 | 0.2 |
| 110 | 6 | LH seeds, exactly 3: (2,1), r in E3 | exact ← | [56,75) / [37,50) | [0,37)/[0,25) = 2, [37,56)/[25,37) = 1 | served by step 17's walk | 0.0 | 0.0 |
| 111 | 6 |  | extend ← | [0,56) / [0,37) | [0,37)/[0,25) = 2, [37,56)/[25,37) = 1 |  | 0.9 | 0.3 |
| 112 | 6 | LH seeds, exactly 3: (2,1), r in E4 | exact → | [37,56) / [25,37) |  | searched | 5.6 | 3.4 |
| 113 | 6 |  | extend → | [56,75) / [37,50) | [56,75)/[37,50) = 1 |  | 1.3 | 0.3 |
| 114 | 6 |  | cross ← (switch) | [0,37) / [0,25) | [0,37)/[0,25) = 2 |  | 0.4 | 0.2 |
| 115 | 6 | cross → RH (seeds list 9) | cross →, ≤ 6 in all | [75,150) / [50,100) |  |  | 0.3 | 0.2 |
| 116 | 6 | cross → RH (seeds list 10) | cross → (switch), ≤ 6 in all | [75,150) / [50,100) |  |  | 0.3 | 0.2 |
| 119 | 7 | chain → RH, list 2 | resume →, ≤ 7 in all | [75,150) / [50,100) |  |  | 0.0 | 0.0 |
| 121 | 7 | chain ← LH, list 3 | resume ←, ≤ 7 in all | [0,75) / [0,50) |  |  | 0.0 | 0.0 |
| 123 | 7 | chain → RH, list 4 | resume →, ≤ 7 in all | [75,150) / [50,100) |  |  | 0.0 | 0.0 |
| 125 | 7 | chain ← LH, list 5 | resume ←, ≤ 7 in all | [0,75) / [0,50) |  |  | 0.0 | 0.0 |
| 127 | 7 | chain → RH, list 6 | resume →, ≤ 7 in all | [75,150) / [50,100) |  |  | 0.1 | 0.1 |
| 129 | 7 | chain ← LH, list 7 | resume ←, ≤ 7 in all | [0,75) / [0,50) |  |  | 0.1 | 0.1 |
| 131 | 7 | chain → RH, list 8 | resume →, ≤ 7 in all | [75,150) / [50,100) |  |  | 0.1 | 0.1 |
| 134 | 7 | RH seeds, exactly 3: (0,3) | exact ← | [112,150) / [75,100) | [75,112)/[50,75) = 3 (after-seed) | served by step 6's walk | 0.0 | 0.0 |
| 135 | 7 |  | extend ← | [75,112) / [50,75) | [75,112)/[50,75) = 3 (after-seed) |  | 0.7 | 0.0 |
| 136 | 7 | RH seeds, exactly 3: (3,0) | exact → | [75,112) / [50,75) |  | served by step 32's walk | 0.0 | 0.0 |
| 137 | 7 |  | extend → | [112,150) / [75,100) | [112,150)/[75,100) = 3 |  | 0.7 | 0.0 |
| 138 | 7 | RH seeds, exactly 3: (1,2), p in E2 | exact ← | [131,150) / [87,100) | [112,131)/[75,87) = 1 (narrow), [75,112)/[50,75) = 2 | served by step 6's walk | 0.0 | 0.0 |
| 139 | 7 |  | extend ← | [75,131) / [50,87) | [112,131)/[75,87) = 1 (narrow), [75,112)/[50,75) = 2 |  | 1.0 | 0.1 |
| 140 | 7 | RH seeds, exactly 3: (1,2), p in E1 | exact → | [112,131) / [75,87) | [131,150)/[87,100) = 1 (narrow) | searched | 4.5 | 1.3 |
| 141 | 7 |  | extend → | [131,150) / [87,100) | [131,150)/[87,100) = 1 (narrow) |  | 0.9 | 0.1 |
| 142 | 7 |  | cross ← (switch) | [75,112) / [50,75) | [75,112)/[50,75) = 2 |  | 0.4 | 0.0 |
| 143 | 7 | RH seeds, exactly 3: (2,1), r in E3 | exact → | [75,93) / [50,62) | [112,150)/[75,100) = 2, [93,112)/[62,75) = 1 | served by step 32's walk | 0.0 | 0.0 |
| 144 | 7 |  | extend → | [93,150) / [62,100) | [112,150)/[75,100) = 2, [93,112)/[62,75) = 1 |  | 1.0 | 0.1 |
| 145 | 7 | RH seeds, exactly 3: (2,1), r in E4 | exact ← | [93,112) / [62,75) |  | searched | 4.7 | 1.4 |
| 146 | 7 |  | extend ← | [75,93) / [50,62) | [75,93)/[50,62) = 1 |  | 1.0 | 0.1 |
| 147 | 7 |  | cross → (switch) | [112,150) / [75,100) | [112,150)/[75,100) = 2 |  | 0.4 | 0.0 |
| 148 | 7 | cross ← LH (seeds list 9) | cross ←, ≤ 7 in all | [0,75) / [0,50) |  |  | 0.2 | 0.1 |
| 149 | 7 | cross ← LH (seeds list 10) | cross ← (switch), ≤ 7 in all | [0,75) / [0,50) |  |  | 0.2 | 0.1 |

prefetch pass: 15.4% (all), 13.3% (-m 1)

## Where the cycles go, by K

Share of all search cycles, 150 bp, both strands:

| | `-n 5` all | `-n 5 -m 1` | `-n 6` all | `-n 6 -m 1` | `-n 7` all | `-n 7 -m 1` |
|---|---:|---:|---:|---:|---:|---:|
| cycles per read | 19,237 | 10,232 | 26,135 | 11,686 | 33,220 | 12,407 |
| exact searches still run, per strand | 6 | 6 | 8 | 8 | 10 | 10 |
| exact walks (halves, quarter walks) | 38.6% | 54.6% | 29.0% | 48.4% | 23.2% | 46.4% |
| other exact searches (middle seeds) | 17.6% | 13.5% | 25.5% | 18.9% | 30.4% | 21.2% |
| windowed extensions (seed cases) | 18.7% | 8.5% | 20.7% | 9.0% | 22.0% | 8.9% |
| … of which the same seed over the same segment | 9.9% | 4.7% | 9.1% | 4.6% | 8.7% | 4.2% |
| crossings and switch crossings | 5.3% | 4.0% | 5.5% | 4.4% | 5.7% | 4.5% |
| chains (resumes) and the 0\|0 extensions | 4.5% | 5.6% | 3.8% | 5.6% | 3.4% | 5.7% |
| prefetch pass | 15.3% | 13.7% | 15.4% | 13.7% | 15.3% | 13.3% |

## Findings

### 1. Repeated exact walks: none left

- **Every served search was served.** None of the 18 served searches at `-n 7` ever
  ran index steps of its own.
- **The 10 searches that run all differ in origin or direction:**
  - the two halves;
  - the two inner quarter walks, [37,75) leftward and [75,112) rightward;
  - the six middle seeds H_mid, E2 and E3 of each half.
- **None of them lies on another's walk.** A one-direction index cannot reuse a walk
  from another origin. For example, [37,56) rightward starts at 37, and no other
  search walks rightward from 37.
- **Their length hardly matters.** A search costs 750–1,000 cycles a strand, from the
  18-base E2 seed up to the 75-base half. On E. coli a row is single after about 11
  bases, and then one locate and one genome compare finish the search.

### 2. The same seed extended in several rounds (removable, at most about 5%)

| seed | extended over | rounds (exactly N mismatches) | `-n 7` all | `-n 7 -m 1` |
|---|---|---|---:|---:|
| LH Q_out [0,37) → | Q_in [37,75) | 2 (N=1, narrow), 4 (B, N=2), 6 ((0,3), N=3) | 2.1% | 1.3% |
| LH Q_in [37,75) ← | Q_out [0,37) | 2 (N=1, narrow), 4 (A, N=2), 6 ((3,0), N=3) | 2.2% | 1.3% |
| RH Q_out [112,150) ← | Q_in [75,112) | 3, 5, 7 | 2.0% | 0.8% |
| RH Q_in [75,112) → | Q_out [112,150) | 3, 5, 7 | 2.3% | 0.9% |

- **What repeats:** each round extends the same seed range over the same quarter. Only
  the exact mismatch count its window asks for differs. The rows are nearly always
  located by then (0.4 located extensions per read for each step), so each repeat is
  mostly one more genome compare and its bookkeeping. The N = 3 tree includes the N = 1
  and N = 2 trees, except where the window's reach check prunes them near the end.
- **How it could be removed:** extend once, in the first round, and keep what goes over
  that round's count for the later rounds. That is how the plain extensions already
  defer (`Extend_Located`: "above Count, deferred unchanged"). This would take a
  windowed form of that deferral and a list for each of the four seeds.
- **What it would save:** at most the cost of all but one extension in each group:
  4.7%, 5.2% and 5.5% of the all-hits search at `-n 5/6/7`, and 1.3–1.5% at `-m 1`.
  With `-m 1` the later rounds seldom run, so there is little to save.

### 3. The same locus compared again (cheap)

| | `-n 5` all | `-n 5 -m 1` | `-n 6` all | `-n 6 -m 1` | `-n 7` all | `-n 7 -m 1` |
|---|---:|---:|---:|---:|---:|---:|
| genome compares per read | 9.0 | 6.0 | 11.4 | 6.5 | 13.7 | 6.8 |
| … of a read strand and locus already compared | 5.7 | 3.4 | 7.4 | 3.8 | 9.0 | 3.9 |
| cycles of those repeats, share of the search | 2.2% | 2.4% | 2.1% | 2.3% | 2.1% | 2.3% |

- **Why the same locus comes back:** on the strand where a read maps, each seed that
  matches exactly lands on its true locus, and each case compares the read there.
- **Why the repeats are cheap:** the first compare at a locus averages about 1,300
  cycles, because the genome words are not yet in cache. The repeats average about 80
  cycles.
- **Not worth removing:** a per-lane cache of compared loci would save under 2%. Even
  then each search would still have to walk and locate to learn the locus.
- **Caveat:** the timer adds its own overhead (search time +17% at `-n 7`). Counts are
  exact; cycle shares are approximate.

### 4. Chains and crossings: no rework

- **Resumes continue where the last round stopped.** A resume takes only the ranges the
  previous round deferred at its count, and picks each up at its saved depth. The pops
  per read fall round by round. The → RH chain from the exact LH half takes 0.32 pops per
  read (0|1), then 0.15, 0.07, 0.03, 0.02, 0.01 and 0.00 in rounds 2 to 7.
- **Each family's crossings take only that family's seeds.** The crossings of a family
  run over its own seed lists (9 and 10). The seeds' rejects then become a new chain.
- **Total:** 9–10% of the search for the chains, the crossings and the 0|0 extensions
  together, at every K.

### 5. Both strands

The reverse-complement strand runs the whole program again. With an index of one strand
this is inherent, so the audit does not count it as redundant. On the strand where the
read does not map, most searches end before they locate. Per read, a search locates on
  about 1.05 of its two strands.

## What would close the gap

- **The remaining lever is the number of distinct seeds.** Exact searches cost the same
  at any length on E. coli, so the 6 middle seeds (`-n 7`) cost as much as the halves.
- **Fewer, longer seeds:** Columba's search schemes cover the same mismatch placements
  with fewer searches, each extended both ways in a bidirectional index.
- **Human may differ:** on the 3 Gb human genome a row becomes single later (about 16
  bases). There the 12-base eighths at 100 bp likely stay wide to the end and cost
  more than a half. This is not measured.
- **The prefetch pass:** at 13–15% it is the largest single item that is not a search.
  It was already noted in `shared_walks.md`.

## Reproducing

Scratch tools, not in the tree:
- **`program_dump.cpp`** (includes `batmis1.cpp` via `MAPPER_SOURCE`): prints a
  program's steps at a read length.
- **`patch.py`**: the per-step rdtsc profiling build. **`patch_loci.py`**: the
  genome-compare timer.
- **`audit.py`**: joins a program dump with the two-strand profile.
- **`tables.py`**: writes the table above.

The profiles use the same reads, genome and core as `shared_walks.md`.
