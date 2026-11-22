# Same count: dropping a half's seed cases once a seed is as rare as the half

## Summary

- **The rule.** Let H be a half of the read (LH or RH) and X a substring of it.
  - **Counts:** every occurrence of H contains X at the same offset, so |SA(X)| ≥ |SA(H)|.
  - **Equal counts:** if they are equal, every occurrence of X lies inside an occurrence of
    H. So wherever X matches exactly, H matches exactly too.
  - **Consequence:** every seed case of H (the 1|k, 2|k and 3|k families, and -n 2's thirds
    cases (a) and (b)) needs at least one mismatch in H. Once its exact seed is as rare as
    the exact half, the case has no hit.
- **Stopping mid-search.** A search's range only shrinks as it extends, and every longer
  prefix of the seed is still inside H. So the seed's search stops as not found at the
  first depth where its range size equals |SA(H)|.
  - **Where it stops:** for a read whose half occurs once, at the depth the range becomes a
    single row, on average 12.4 bases on E. coli and 14.6 on chr21.
  - **What is skipped:** the rest of the walk, the extension, the crossing and the switch.
- **Where it applies:** only when the half matched exactly (|SA(H)| ≥ 1). With mismatches
  in both halves nothing changes, and on the strand a read does not map on the half
  almost never matches, so that strand costs nothing.
- **Speed** on default-built indexes (every SA row sampled), `-n 7`, all hits, both mappers:

  | reads | E. coli, 150 bp | chr21, 150 bp | chr21, 100 bp |
  |---|---:|---:|---:|
  | exact | 32–34% faster | 12–17% | 3–5% |
  | one mismatch | 15–17% | 4–9% | 2–4% |
  | all mismatches in one half (0\|7, 7\|0) | 17–21% | 8–12% | 2–7% |
  | clustered | 13–17% | 7–10% | 0–4% |
  | mismatches in both halves | no change | no change | no change |
  | random (no hit) | no change | no change | no change |

  At `-n 5` exact reads are 28–29% faster on E. coli. In cycles (`perf stat`), random
  reads cost 0.5–1%.
- **Same output:** SAM output is identical to HEAD (9fa3780) and to `--no-same-count`.
  - **Mappers and indexes:** batmis1 and batmis2; E. coli on indexes sampled at every row
    and at every 8th row, the latter also with `--known`; chr21 at every row, 150 and
    100 bp.
  - **Settings:** `-n 1..7`, all hits and `-m 1`, batch 16 and 1.
  - **Reads:** 20 read layouts and 100k simulated reads.
  - **`make test`** also passes.
- **On by default.** `--no-same-count` turns it off, for timing.

## How it works (`src/batmis1.cpp`, `src/batmis2.cpp`)

- **Marking the seeds:** `Mm_Step::Half` marks the `OP_SEARCH` of each seed case with its
  half (1 LH, 2 RH). `Mm_Builder::Mark_Half` sets it for every search that `Seeds` and
  `Three_Seeds` add, and the -n 2 thirds cases (a) and (b) are marked by hand. Each of
  these seeds lies inside its half, and its case needs a mismatch there.
- **The half's count:** when the exact search of LH or RH finishes, `Mm_Half_Done` keeps its
  range size: `LH_Count` and `RH_Count`, per strand (`Mm_Strand`, swapped by
  `Mm_Set_Strand`). It is 0 if the half did not match.
- **At a marked `OP_SEARCH`, with C the half's count (C > 0):**
  - **A served seed** (a kept range or a shared walk), or one the lookup table already
    covers: if its size is C, the case is dropped (`Found = false`, as the window skip
    does). Its extension then does not run, and its seed lists stay empty for the crossings.
  - **An own search:** `Same_Target = C`. `Mm_Search_Step` (MM_WIDE) ends the search as not
    found when the range size reaches C.
  - **Shared walks:** the search's pending walk points become failed. They are longer
    prefixes from the same origin inside the same half, so they are at least as rare, and
    every search they serve is a seed case of that half.
  - **Windows:** the failure depths for the windows come only from the half searches
    (`Mm_Search_Failed`), so they are unaffected.

Counts per read, batmis1, `-n 7`, all hits, default-built indexes (`-DBATMAN_BATCH_STATS`,
line `same_count:`):

| reads | genome | seed cases dropped, served | own searches stopped | mean depth |
|---|---|---:|---:|---:|
| exact (u0) | E. coli 150 bp | 10.1 | 8.0 | 12.4 |
| 2 mismatches (u2) | E. coli 150 bp | 2.5 | 2.0 | 12.4 |
| 7 mismatches (u7) | E. coli 150 bp | 0.07 | 0.06 | 12.3 |
| 0\|7 (s0_7) | E. coli 150 bp | 5.0 | 4.0 | 12.4 |
| exact (u0) | chr21 150 bp | 9.8 | 7.6 | 14.6 |
| exact (u0) | chr21 100 bp | 7.9 | 4.5 | 14.0 |

With both halves exact, nearly all of the read's 18 seed cases on the mapping strand are
dropped. At 100 bp the eighth seeds are only 12 bases long. Many are rarer than the half
only when they are complete, or never, so fewer of them are dropped.

## Timings

**Setup:**
- **Reads:** 20,000 reads per layout, simulated from the genome with planted
  substitutions.
- **Indexes:** built by the current `bwt_build_index` with its default sampling, which
  stores every SA row for both genomes.
- **Run:** all hits, Hamming, one P-core, SAM output.
- **BatMis time:** search + output_format.
- **Columba time:** its own mapping time (Columba 1.3, `-e K -a all`).
- **Repeats:** median of 5 interleaved runs (3 on chr21).

The machine's speed drifted between sessions, so compare times within a table only.

**Layouts:**
- **uM:** M mismatches at random positions.
- **sI_J:** I mismatches in the left half and J in the right.
- **q3_4:** three mismatches at the left half's quarter points and four in the right half.
- **spread7:** one mismatch near each eighth of the read.
- **cluster7:** seven mismatches within 20 bases.
- **random:** random sequence, with no hit.

### E. coli, 150 bp, `-n 7`

| reads | batmis1 HEAD (s) | batmis1 new (s) | new / HEAD | batmis2 new / HEAD | Columba (s) | Columba / batmis1 new |
|---|---:|---:|---:|---:|---:|---:|
| u0 | 0.387 | 0.257 | 0.66 | 0.68 | 0.504 | 1.96 |
| u1 | 0.361 | 0.299 | 0.83 | 0.85 | 0.490 | 1.64 |
| u2 | 0.349 | 0.325 | 0.93 | 0.93 | 0.483 | 1.49 |
| u4 | 0.311 | 0.313 | 1.01 | 1.00 | 0.450 | 1.44 |
| u7 | 0.271 | 0.276 | 1.02 | 1.03 | 0.402 | 1.46 |
| s0_7 | 0.317 | 0.258 | 0.81 | 0.83 | 0.421 | 1.63 |
| s7_0 | 0.330 | 0.262 | 0.79 | 0.82 | 0.424 | 1.62 |
| s3_4 | 0.255 | 0.261 | 1.02 | 1.03 | 0.411 | 1.57 |
| q3_4 | 0.248 | 0.248 | 1.00 | 1.02 | 0.371 | 1.50 |
| spread7 | 0.207 | 0.208 | 1.00 | 1.01 | 0.228 | 1.10 |
| cluster7 | 0.322 | 0.268 | 0.83 | 0.87 | 0.382 | 1.43 |
| random | 0.149 | 0.152 | 1.02 | 1.01 | 0.210 | 1.38 |

At `-n 5`, batmis1 new / HEAD: u0 0.71, u1 0.85, u2 0.94, u4 1.01, s0_7 0.81, random 1.01.
batmis2 gives 0.72, 0.87, 0.95, 0.99, 0.82 and 1.02.

### chr21, `-n 7`

| reads | batmis1 new / HEAD, 150 bp | batmis2, 150 bp | batmis1, 100 bp | batmis2, 100 bp | Columba / batmis1 new, 150 bp |
|---|---:|---:|---:|---:|---:|
| u0 | 0.83 | 0.88 | 0.95 | 0.97 | 2.04 |
| u1 | 0.91 | 0.96 | 0.96 | 0.98 | 1.92 |
| u2 | 0.95 | 0.97 | 0.98 | 0.99 | 1.82 |
| u4 | 0.98 | 1.01 | 1.01 | 0.98 | 1.92 |
| u7 | 0.99 | 1.02 | 1.00 | 1.03 | 2.01 |
| s0_7 | 0.88 | 0.91 | 0.93 | 0.98 | 1.93 |
| s7_0 | 0.88 | 0.92 | – | – | 1.89 |
| cluster7 | 0.90 | 0.93 | 0.97 | 1.00 | 1.93 |
| random | 1.01 | 1.04 | 1.04 | 1.04 | 1.96 |

Differences of 2–3% on chr21 are within run-to-run noise. On 100 bp chr21 reads, Columba
takes 4–16 times as long as batmis1 at `-n 7`.

## Quarters as super-seeds

- **The general form:** a case's exact seed X lies inside an earlier exact string Y whose
  count is known, and |X| = |Y| ≥ 1. If the case needs a mismatch inside Y (it can only be
  in Y \ X), the case has no hit.
- **Super-seeds:** besides the halves, the program knows four exact counts when the seed
  cases run:
  - **The outer quarters** [0,37) and [112,150) at 150 bp, kept by the half searches in
    round 1 (`Mm_Keep_Rest`).
  - **The inner quarters** [37,75) and [75,112), searched whole as the s=1 seeds in rounds 2
    and 3 (steps marked `Records`).
  - **Nothing else:** no other searched string (the thirds and eighths on the walks, or the
    middle seeds) contains another case's seed together with the region where that case
    needs its mismatch.

Every seed case of `-n 7` at 150 bp. Steps are the single-strand program's; RH mirrors are
in brackets.

| step (RH) | round | case | seed | needs a mismatch in | dropped when the seed is as rare as |
|---|---|---|---|---|---|
| 15 (30) | 2 (3) | s=1, in Q_in | Q_out [0,37) | Q_in | the half |
| 17 (32) | 2 (3) | s=1, in Q_out | Q_in [37,75) | Q_out | the half |
| 47 (73) | 4 (5) | A, 2 in Q_out | Q_in | Q_out | the half |
| 49 (75) | 4 (5) | B, 2 in Q_in | Q_out | Q_in | the half |
| 51 (77) | 4 (5) | C1, p in [25,37) | H_out [0,25) | [25,37) | the half, or **Q_out** |
| 53 (79) | 4 (5) | C2, one in [37,50) | H_in [50,75) | [37,50) | the half, or **Q_in** |
| 55 (81) | 4 (5) | C3 | H_mid [25,50) | H_out, H_in | the half |
| 101 (134) | 6 (7) | (0,3) | Q_out | Q_in | the half |
| 103 (136) | 6 (7) | (3,0) | Q_in | Q_out | the half |
| 105 (138) | 6 (7) | (1,2), p in E2 | E1 [0,18) | E2 | the half, or **Q_out** |
| 107 (140) | 6 (7) | (1,2), p in E1 | E2 [18,37) | E1 | the half, or **Q_out** |
| 110 (143) | 6 (7) | (2,1), r in E3 | E4 [56,75) | E3 | the half, or **Q_in** |
| 112 (145) | 6 (7) | (2,1), r in E4 | E3 [37,56) | E4 | the half, or **Q_in** |

- **How many cases the quarters add:** 4 at `-n 5` (C1 and C2 of each half), 8 at `-n 6`,
  12 at `-n 7`. The same containments hold at 100 bp (`Part_Position`'s fractions).
- **Where they help:** only when the half has a mismatch but the quarter is exact. Reads
  with both halves exact were already handled by the half rule.
- **Implementation:**
  - `Mm_Step::Super` names the quarter.
  - The quarter's count is used in place of the half's when it is known. X ⊂ Q keeps
    |X| ≥ |Q| ≥ |H|, so a search can meet |Q| but never pass it.
  - If the inner quarter's own search was stopped by the half rule, its count is |H|.
  - None of these searches is a shared walk, which the builder checks.
  - The sizes are kept per strand in `Mm_Lane::Same_Count`, indexed by the strand and not
    swapped. Swapping them at every strand switch cost random reads 2.4%.

Seed cases dropped by a quarter alone, per read, batmis1, `-n 7`, E. coli 150 bp: u1 3.1,
u2 3.8, u4 3.1, u7 1.5, s1_6 3.1, s3_4 1.1, q3_4 0.3. None on exact or random reads.

Search + output time, new / the half rule alone (e07b661), default-built indexes, `-n 7`:

| reads | batmis1 E. coli 150 bp | batmis2 E. coli 150 bp | batmis1 chr21 150 bp | batmis2 chr21 150 bp | batmis1 chr21 100 bp |
|---|---:|---:|---:|---:|---:|
| u1 | 0.95 | 0.97 | 0.99 | 0.99 | 1.00 |
| u2 | 0.95 | 0.97 | 0.97 | 0.97 | 1.01 |
| u4 | 0.95 | 0.96 | 0.99 | 0.99 | 1.00 |
| u7 | 0.99 | 1.00 | 1.00 | 0.99 | 1.00 |
| s1_6 | 0.95 | 0.97 | 0.98 | 0.98 | 1.01 |
| s2_5 | 0.98 | 0.98 | 0.99 | 1.00 | – |
| cluster7 | 0.96 | 0.98 | 0.97 | 1.00 | – |
| u0, s0_7, s3_4, q3_4 | 1.00–1.01 | 1.00–1.01 | 0.99–1.02 | 1.00–1.01 | – |
| spread7 | 1.02 | 1.02 | 1.00 | 1.01 | – |
| random | 1.01 | 1.04 | 1.00 | 1.02 | 1.01 |

- **At `-n 5`:** batmis1 u2 0.98, s1_6 0.97.
- **In cycles (`perf stat`):**
  - **u2:** batmis1 −5%, batmis2 −4%.
  - **Random and spread7, batmis1:** within 1%.
  - **Random and spread7, batmis2:** 1.3–3% higher at the same instruction count, and
    1.3–1.7% with aligned builds (`-falign-functions=64 -falign-loops=64`). That is within
    the 3% that repeated runs of one binary vary, so it is code layout, not work.
- **At 100 bp, no gain:** the 12-base eighths rarely become as rare as their quarter.

## The exact read as the super-seed (RH's k|0 cases)

- **The idea:** R = LH|RH matches exactly. If RH is as rare as something that already
  pins it to R, then extending RH leftward over LH can only lead back to the exact hit,
  so the cases that keep RH exact can be dropped.
- **The target must be |SA(R)|, not |SA(LH)|.**
  - Let S_k = LH[k..] + RH be RH's leftward extension, so RH ⊂ S_k ⊂ R.
  - |S_k| = |R| puts every exact S_k inside an exact R. Any alignment with S_k exact is
    then the exact hit, and every case needs a mismatch somewhere.
  - S_k does not contain LH, so |S_k| = |LH| says nothing on its own. For example, take LH
    3 times, R once and RH 3 times: S_k can keep count 3, and 2 of those loci can be real
    k|0 hits.
  - The two targets agree when |LH| = |R|, e.g. both unique.
- **When |R| is known:** at `OP_COUNT_EXACT`, after the 0|0 extension and before RH's
  search, per strand (`Mm_Lane::Exact_Count`).

What it prunes:

| seed exact | cases | new? |
|---|---|---|
| RH, when \|RH\| = \|R\| at the end of its search | the 1\|0 extension over LH and its whole chain, 2\|0 … K\|0 | yes |
| T2, the -n 2 middle third, which spans both halves | case (c) | yes, `-n 2` |
| any seed X inside a half H, quarters included | H's seed cases | no: \|X\| ≥ \|H\| ≥ \|R\|, so X meets \|H\| first and the half rule drops it |
| LH | 0\|k | can't: it runs before \|R\| is known, and it finds R |

**What it saved:**
- **Firing rate:** RH's k|0 cases are dropped about once per exact read (exact_rh 1.01 on
  E. coli, 1.03 on chr21), on the strand that has the hit.
- **Cycles (`perf stat`), new / 1421df9, E. coli 150 bp:**
  - **Exact reads:** at `-n 7` batmis1 0.97, batmis2 0.98; at `-n 2` batmis1 0.89,
    batmis2 0.91.
  - **Random reads:** 1.00–1.02 at +0.2% instructions, which is code layout as above.
- **The middle third's share at `-n 2`:**
  - **Firing rate:** T2 is dropped 1.00 times per exact read on E. coli and 1.02 on chr21
    (`exact_seed`). That is 0 at `-n 7`, which has no seed outside a half.
  - **Exact-read cycles, batmis1, relative to 1421df9:**
    - RH drop alone: 0.95 on both genomes.
    - With T2 too: 0.90 on both.
  - So T2 gives about half the `-n 2` gain.
- **Search and output time, exact reads:**
  - **`-n 7`:** E. coli 150 bp 0.97–0.98, chr21 150 bp 0.98–0.99, chr21 100 bp 1.00.
  - **Every other layout:** within ±2%.
- **The rest of the idea, not built:**
  - **Stopping RH's own search at |R|:** this would also save its locate and genome
    compare, but the kept quarter and the shared walks need the search to reach the
    quarter's depth.
  - **Testing the exact branch inside the 1|0 extension** (RH repeated, so S_k reaches |R|
    inside LH): this is the branch test below, which saved under 1%.

## Pruning extension branches by the same rule (not adopted)

- **The rule is valid in extensions too.** Take a branch of a seed case's windowed
  extension that is still exact. Its string is a substring of the half: the seed plus the
  bases matched. If its range size has fallen to |SA(H)|, every way it can continue inside
  H follows H's own bases, so the case's required mismatch in H cannot come. The branch
  can be dropped.
- **Why equality is the trigger:** an exact substring is never rarer than its superstring,
  so the test is equality, not "smaller".
- **Branches with mismatches:** these are not read substrings, and no known count bounds
  them.
- **Measured:** it dropped 0.16 branches per exact read on E. coli at 150 bp, 1.6 on chr21
  at 150 bp, and 7.7 on chr21 at 100 bp.
- **What it saved:** under 1% of instructions, and nothing in cycles. The branches it drops
  end soon anyway, once their windows can no longer reach their counts.
- **Left out** of the hot loop. A note at `Range_Count` in the code points here.

## Reproducing

Scratch scripts, not in the tree, as in `known_positions.md`:
- **Reads:** `make_classes.py genome.fa outdir 20000 [length]`.
- **Timings:** `final_d.sh mapper genome K rounds layouts…`, with `SUF=_100` for 100 bp.
- **Output checks:** `check_f.sh` with `FLAG=--no-same-count`.
- **Counts:** a `-DBATMAN_BATCH_STATS` build prints `same_count:` with `--timing`.
