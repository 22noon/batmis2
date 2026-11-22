# Known positions: locating single-row paths from the read's earlier hits

## Summary

- **Off by default; `--known` turns it on** in both mappers. It pays only where single-row
  paths walk far to a sampled row, which the current index builder avoids for genomes up
  to 134 Mb by sampling every row. See "Why the SA sampling decides it" below.
- **The rule.**
  - **Why it holds:** a single-row path's string occurs exactly once in the genome. Suppose
    the read, with the path's own mismatches and no others, matches the genome over that
    string's interval at a position already located for this read strand. Then that one
    occurrence is there.
  - **What it saves:** the path is located at once (`Known_Find`). It does not take LF steps
    to a sampled row, or look up the SA sample.
  - **Why the output cannot change:** the compare that follows (`MM_LOCATED`) is the one
    the walk would have reached. So the path's mismatches, failure depth, deferrals and hits
    are the same.
- **Speed on E. coli, index sampled at every 8th row, all hits, 150 bp, `-n 7`:**
  - **Faster, exact, near-exact or clustered reads:** batmis1 10–15%, batmis2 5–12%.
  - **Faster, spread mismatches:** 3–4%.
  - **Slower:** 2–4% for reads that do not map.
- **Speed on chr21, index sampled at every 8th row:** about neutral.
  - **Faster:** 1–2% for exact and clustered reads.
  - **Slower:** 1–2% for spread mismatches, and 2.5–5% for reads that do not map.
  - **Likely why (not measured separately):** a row becomes single later (about 14 bases
    rather than 11), so less of the walk is left to skip. And the index is not in cache, so
    the batch's prefetching already hides much of the walk.
- **On an index sampled at every row** (the builder's default up to 134 Mb), no gain at
  all: 1–4% slower on E. coli, in every layout.
- **Same output:** SAM output is identical to HEAD (9cefed5) with the table on and off.
  - **Mappers and genomes:** batmis1 and batmis2, E. coli and chr21.
  - **Settings:** `-n 1..7`, all hits and `-m 1`, batch 16 and 1.
  - **Reads:** the 20 read layouts below and 100k simulated reads.
  - **`make test`** also passes.
- **Against Columba:** with the table on, batmis1 is ahead of Columba on the E. coli index
  sampled at every 8th row in every layout at `-n 7`. Without the table, Columba was level
  or slightly ahead there on exact and clustered reads. On the index sampled at every row,
  batmis1 is ahead without it (u0: 0.30 s, Columba 0.37 s).

## How it works (`src/batmis1.cpp`, `src/batmis2.cpp`)

- **`SARange::Anchor`:** the read position where the path's string starts. That is its first
  base for revfmi paths, and one past its last for fwfmi paths. Together with the task's
  Offset and Level, it gives the string's interval.
  - batmis2 already kept this field for its other-index locate.
  - batmis1 now sets it where a search or index switch starts, in one spare byte of the
    struct, with `NO_ANCHOR` meaning unknown. It is copied with the range everywhere
    else.
- **The table:** `Known_Table`, per lane and kept outside `Mm_Lane`. It holds up to 8 read
  starts per strand that a Locate accepted.
  - **What is added:** a position is kept only if the whole read has at most `-n`
    mismatches there.
  - **Why the filter matters:** without it, the strand a read does not map on adds the
    positions where its seeds occur by chance. Each later lookup would then read the
    genome there for nothing.
  - **When it is added:** at `MM_LOCATED`, after `Mm_Prefetch` has fetched the genome
    words. Adding at `MM_LOCATE` took a full memory miss per check and made chr21 2–7%
    slower.
- **The lookup:** where a path becomes single and is not yet located (`Mm_Search_Enter`,
  `Mm_Single_Enter`). For each known position, it compares the genome over the interval
  with the read and the path's mismatches: the same count, at the same read positions,
  with the same genome bases.
  - **Mismatches included:** this covers mismatch branches and switched seeds too, not
    only exact seeds. A path that comes back to the read's position by any route is
    located there.
- **The guard:** the lookup and the add are out of line behind a one-byte test
  (`Known_Here`), so the walk's code is unchanged when the table is empty.
- **`--known`:** turns it on. Off, the lookup and the add stop at their first test. That
  build's search time is within 1% of HEAD's.

Counts per read, batmis1, `-n 7`, all hits (`-DBATMAN_BATCH_STATS`, line `known:`):

| reads | genome | lookups | located from the table | add checks | positions added |
|---|---|---:|---:|---:|---:|
| exact (u0) | E. coli | 9.9 | 9.6 | 2.0 | 1.0 |
| 7 random mismatches (u7) | E. coli | 6.4 | 4.5 | – | 1.0 |
| random sequence | E. coli | 0 | 0 | 1.6 | 0 |
| exact (u0) | chr21 | 27.9 | 10.9 | 14.6 | 1.2 |
| 7 random mismatches (u7) | chr21 | 9.4 | 4.1 | 5.5 | 1.0 |

On chr21, seeds also reach repeat copies of the read's sequence, where the read has more
than `-n` mismatches. Those positions are checked but not added.

## Timings

**Setup:**
- **Reads:** 20,000 reads of 150 bp per layout, simulated from the genome with planted
  substitutions.
- **Run:** all hits, Hamming, one P-core (`taskset -c 2`), SAM output.
- **BatMis time:** search + output_format from `--timing`.
- **Columba time:** its "finished mapping" time. This is Columba 1.3, default strategy,
  `-e K -a all`.
- **Repeats:** median of 5 interleaved runs (3 on chr21). "new / HEAD" < 1 means faster.

**Layouts:**
- **uM:** M mismatches at random positions.
- **sI_J:** I mismatches in the left half and J in the right.
- **q3_4:** three mismatches at the left half's quarter points and four in the right half.
- **spread7:** one mismatch near each eighth of the read.
- **cluster7:** seven mismatches within 20 bases.
- **random:** random sequence, with no hit.

### batmis1, E. coli, `-n 7`

| reads | HEAD (s) | new (s) | new / HEAD | Columba (s) | Columba / new |
|---|---:|---:|---:|---:|---:|
| u0 | 0.400 | 0.338 | 0.85 | 0.369 | 1.09 |
| u1 | 0.381 | 0.329 | 0.86 | 0.355 | 1.08 |
| u2 | 0.365 | 0.322 | 0.88 | 0.355 | 1.10 |
| u4 | 0.318 | 0.291 | 0.92 | 0.350 | 1.20 |
| u7 | 0.269 | 0.257 | 0.96 | 0.322 | 1.25 |
| s0_7 | 0.324 | 0.305 | 0.94 | 0.325 | 1.07 |
| s3_4 | 0.253 | 0.242 | 0.96 | 0.310 | 1.28 |
| s7_0 | 0.322 | 0.297 | 0.92 | 0.330 | 1.11 |
| q3_4 | 0.237 | 0.228 | 0.96 | 0.297 | 1.30 |
| spread7 | 0.231 | 0.223 | 0.97 | 0.278 | 1.25 |
| cluster7 | 0.364 | 0.327 | 0.90 | 0.357 | 1.09 |
| random | 0.142 | 0.147 | 1.04 | 0.216 | 1.47 |

At `-n 5`: u0 0.89, u2 0.93, u4 0.96, random 1.02. Columba / new is 1.47–1.85.

### batmis2, E. coli, `-n 7`

| reads | HEAD (s) | new (s) | new / HEAD | Columba (s) | Columba / new |
|---|---:|---:|---:|---:|---:|
| u0 | 0.425 | 0.374 | 0.88 | 0.372 | 0.99 |
| u2 | 0.386 | 0.349 | 0.90 | 0.365 | 1.05 |
| u4 | 0.346 | 0.319 | 0.92 | 0.346 | 1.08 |
| u7 | 0.299 | 0.289 | 0.97 | 0.314 | 1.09 |
| s0_7 | 0.352 | 0.327 | 0.93 | 0.329 | 1.01 |
| s7_0 | 0.353 | 0.330 | 0.93 | 0.330 | 1.00 |
| s3_4 | 0.277 | 0.266 | 0.96 | 0.314 | 1.18 |
| cluster7 | 0.399 | 0.380 | 0.95 | 0.354 | 0.93 |
| random | 0.166 | 0.169 | 1.02 | 0.220 | 1.30 |

At `-n 5`: u0 0.90, u2 0.93, u4 0.95, random 1.01.

### chr21, `-n 7`

| reads | batmis1 new / HEAD | batmis2 new / HEAD | Columba / batmis1 new |
|---|---:|---:|---:|
| u0 | 0.99 | 0.99 | 1.52 |
| u2 | 0.99 | 0.98 | 1.64 |
| u4 | 0.99 | 0.99 | 1.76 |
| u7 | 1.01 | 1.01 | 1.91 |
| s0_7 | 0.99 | 1.06 | 1.58 |
| s3_4 | 1.01 | 1.02 | 2.04 |
| cluster7 | 0.98 | 0.99 | 1.62 |
| spread7 | 1.02 | 1.00 | 2.07 |
| random | 1.05 | 1.03 | 2.08 |

On chr21, differences of about 2% are within run-to-run noise. In whole-process cycles
(`perf stat`), batmis1 is 3% faster than HEAD on u0 and 2% faster on u7, and 1.3% slower on
random reads.

## Why the SA sampling decides it

- **What the table saves:** the LF steps of a single-row path to the next sampled row,
  where `Locate` reads the stored SA value.
- **Every 8th row sampled:** that walk is about 8 steps (geometric), or up to the end of the
  segment. This is what the timings above use: the bench indexes were built on 2026-11-17,
  before the 1 GB budget (1c7cf9e).
- **Every row sampled:** a path is located the moment it becomes single, so there is
  nothing to skip. The table's add checks and lookups are then pure cost.

Search time (s), E. coli index sampled at every row, `-n 7`, all hits (min of 5 interleaved
runs):

| reads | batmis1 HEAD | batmis1 `--known` | batmis2 HEAD | batmis2 `--known` |
|---|---:|---:|---:|---:|
| u0 | 0.303 | 0.310 | 0.339 | 0.344 |
| u2 | 0.285 | 0.289 | 0.321 | 0.335 |
| u7 | 0.218 | 0.224 | 0.276 | 0.288 |
| cluster7 | 0.289 | 0.294 | 0.362 | 0.371 |
| random | 0.138 | 0.140 | 0.180 | 0.184 |

- **Dense sampling beats the table.** Dense sampling removes the walk for every path. The
  table removes it only at positions the read already located. So HEAD on this index
  (u0: 0.303 s) is faster than the table on the sparse one (0.338 s).
- **Where it could still help: large genomes.** The 1 GB budget keeps every 8th row there:
  for 3 Gb, every row would be about 12 GB an index. chr21 at every 8th row, the closest
  stand-in measured here, is within ±2%. Human, with larger misses per step but more
  repeat copies, is untested.

## Ideas checked but not adopted

- **Skipping an own exact search whose prefix contains a known unique interval.**
  - **The idea:** if a string I is unique and is at position p, any longer string that
    contains I can occur only at p. So the search could be replaced by a compare at p.
  - **How often it applies:** counted at every own `OP_SEARCH`, with the interval recorded
    where each path became single. Only 1.5 of the 16 own searches per read have such an
    interval in their prefix on E. coli u0, 0.8 on u7, and 0.4 on chr21 u0.
  - **Why so few:** the own searches start from distinct origins (audit finding 1), and a
    unique interval sits near the origin of the search that found it.
  - **What it would save:** about 1% of the search. Not worth a second path through the
    walks.
- **Locating every row of a small range (Columba's in-text verification).**
  - **Already done in effect:** BatMis's single-row path already follows the genome's base
    and records mismatches. A range of k rows has at most k live paths, which share steps
    until they split.
  - **No cheaper:** locating a row still needs LF steps to a sampled row, so locating all k
    rows costs about as much as the walk.
  - **The real lever** is the distance to a sample: denser SA samples, batmis2's second
    index, or this table.
- **Deciding a range from its top and bottom rows.**
  - **Same count is not enough:** two rows with n mismatches each say nothing about the
    rows between them. Only a common continuation of the two (their LCP) covers the rows
    between.
  - **Wrong side:** in a one-direction index, rows are sorted by the text on the side
    opposite to the extension. So the rule decides the side already passed, not the one
    ahead.
  - **Little to gain:** unique located seeds already skip the index switch.
  - **Left out:** a possible later aid for repeat-rich genomes.

## Reproducing

Scratch scripts, not in the tree:
- **Reads:** `make_classes.py genome.fa outdir 20000` writes the layouts (seed 7).
- **Timings:** `final.sh mapper genome K rounds layouts…` gives the tables. It runs HEAD,
  this build and Columba, interleaved.
- **Output checks:** `check.sh old new genome readsdir "K…" "m…" [args]` compares sorted
  SAM records of HEAD, this build and this build with the table's option flipped.
- **Counts:** a `-DBATMAN_BATCH_STATS` build prints the `known:` line with `--timing`.
