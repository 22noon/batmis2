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
