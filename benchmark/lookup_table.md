# Start-of-search lookup table: depth, and pre-built vs on the fly

## Summary

Every search in `batman_zero_min` starts from a table that holds the SA range of
each k-base prefix (`LOOKUPSIZE`). The table is built eagerly at startup, one
table per direction. k was fixed at 3 for reads under 28 bp and 6 otherwise. It
is now the `--lookup <k>` option (3–14, default 6; see "Setting the depth").

- **Deeper tables pay off.** Each extra level saves one `Get_SARange` step per
  lookup, about 150–180 ns once the ranges stop fitting in cache (from k ≈ 9 on
  the 500 Mb genome). There are about 2 lookups per read under `-n 1` and 1
  under `-n 0`.
- **Build cost grows 4× per level.** About 80 ns per entry per table: 0.08 s at
  k = 10, 0.3 s at k = 11, 1.2 s at k = 12, 5 s at k = 13, 20 s at k = 14. Memory
  is 8·4^k bytes per table: 8 MB at k = 10, 32 MB at 11, 128 MB at 12.
- **The best depth is set by the read count.** Going from k−1 to k pays off
  once reads ≳ 4^k / 3. On the 500 Mb genome:
  - **k = 11 for 2M reads:** r32 search + build 11.1 s → 9.1 s (−18%), r24
    18.4 s → 15.9 s (−14%).
  - **k = 10 for 500k reads:** p100 2.8 s instead of 3.2 s (−12%), `-n 0`
    1.73 s instead of 1.82 s.
  - **k = 12 is past the break-even** at these read counts: search keeps
    improving, but the 3.1 s build eats the gain.
- **Building on the fly loses.** Filling the table on first use (lazy) never
  beats eager building at a well-chosen depth.
  - Wherever lazy filling would save build time (k ≥ 11 with 2M reads), almost
    every prefix is used only once: 79% at k = 11, 94% at k = 12. So nearly
    every lookup pays the full walk, plus a page fault on the first touch of
    the table.
  - Below that depth the eager build is already cheap (0.08 s at k = 10).
  - To reach a given depth, lazy costs at least as much as the cheaper of
    eager building and plain walking (see the table below).
- **Table layout.** Start and End are now stored together in one array of
  `{Start, End}` pairs, sized 4^k. The old layout was two arrays of 2·4^k each.
  Output is byte-identical and timing is unchanged at k = 6. The layout matters
  at large k: a lookup is one cache miss instead of two, and memory is a quarter
  of what it was.

## Setting the depth

`--lookup <k>` sets the depth, 3–14, default 6.
- **Clamp.** The depth in use is `max(3, min(k, LH − 1))`. The `-n 1` halves
  search on from base k+1, so the key must be shorter than the left half. Reads
  under 8 bp keep the minimum 3, as before.
- **Short reads.** The old "3 for reads under 28 bp" rule is gone: r24 is faster
  at every k from 6 to 11. Output does not depend on k.
- **Choosing k.** From the measurements below, a good choice is
  `k ≈ round(log4(3 · reads))`, kept below about log4(genome size):
  - 10 for ~500k reads;
  - 11 for ~2M reads;
  - 12 from ~10M reads.

  Memory is 8·4^k bytes per index: 32 MB at 11, 128 MB at 12, 2 GB at 14.

## What was measured

`benchmark/lookup_table_scaling.cpp` covers k = 3…14 on the 500 Mb random genome
(revfmi, using the read's first bases as the forward table does). It uses 2M
r32 reads and 500k p100 reads, pinned with `taskset -c 2`, best of 3. For every
read and depth it checks that the eager table, the lazy table, the walk from the
root and the walk from the k=6 table give the same range: 0 disagreements.

What each column measures:

| column | what it measures |
|---|---|
| build | eager 4-ary recursion, one `Get_SARange` per node, both key layouts (first base in the low or high bits; no difference) |
| lk dep | one table read, each waiting for the previous one (as in the search) |
| t6walk | today: k=6 table, then k−6 steps |
| walk | no table: k steps from the first base |
| lazy | fill on first use (k=6 table + walk on a miss), zero pages from `calloc`; per read, all 2M reads |
| distinct | prefixes seen / reads |

Results for r32, 2M reads (ns per read unless noted):

| k | MB | build (s) | lk dep | t6walk | walk | lazy | distinct |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 6 | 0.03 | 0.000 | 1.5 | – | 77 | 0.9 | 0.2% |
| 8 | 0.5 | 0.004 | 3.9 | 52 | 147 | 3.1 | 3.3% |
| 9 | 2 | 0.016 | 8.0 | 170 | 275 | 27 | 13% |
| 10 | 8 | 0.076 | 15 | 332 | 462 | 172 | 45% |
| 11 | 32 | 0.311 | 104 | 500 | 620 | 478 | 79% |
| 12 | 128 | 1.22 | 139 | 683 | 783 | 732 | 94% |
| 13 | 512 | 5.08 | 149 | 882 | 961 | 1042 | 98% |
| 14 | 2048 | 19.8 | 183 | 1026 | 1114 | 1579 | 99% |

Cost per read of reaching depth k, with 2M reads:

| k | eager (build/N + lookup) | lazy | t6walk |
|---:|---:|---:|---:|
| 10 | 53 | 172 | 332 |
| 11 | 260 | 478 | 500 |
| 12 | 749 | 732 | 683 |
| 13 | 2689 | 1042 | 882 |

Eager is cheapest up to k = 11. From k = 12 walking is cheapest, and lazy is
never the best. p100 (500k reads) shows the same picture, with every crossover
one level lower.

Notes:
- **Lookup time jumps between k = 10 and k = 11.** The table outgrows the L2/L3
  cache there: 8 MB → 32 MB.
- **Walking costs grow by about 150–180 ns per level from k ≈ 9.** At that depth
  the ranges' stored-count blocks are spread across the whole index. Levels up
  to ~7 cost ~15–25 ns each, because they hit the same few cached blocks. This is
  why tables of depth ≤ 8 gain little.

## Table fill, and build time on larger genomes

The build skips the subtree under an empty prefix, so a sparse table would
understate the build cost for a larger genome. On the 500 Mb genome the table is
full through k = 13. The table below shows the forward table only; the scratch
tool was `fill.cpp`.

| k | leaves filled | nodes visited | build (s) | ns per node | mean rows per leaf | single-row leaves |
|---:|---:|---:|---:|---:|---:|---:|
| 8 | 100% | 100% | 0.005 | 56 | 7629 | 0% |
| 10 | 100% | 100% | 0.074 | 53 | 477 | 0% |
| 11 | 100% | 100% | 0.296 | 53 | 119 | 0% |
| 12 | 100% | 100% | 1.18 | 53 | 30 | 0% |
| 13 | 99.9% | 100% | 4.52 | 51 | 7.5 | 0.6% |
| 14 | 83% | 99.9% | 17.7 | 50 | 2.2 | 36% |

So the measured build times are for a complete table, with no pruning.

The cost per node is flat, at about 50–56 ns, whatever the range width:
- It covers both `Get_SARange` paths: direct counting below 128 rows, and the
  stored-count lookups from 128 rows up (k ≤ 10 here).
- Each node's four children start from the same range as their parent, so they
  read the same stored-count blocks. Only the first child misses the cache.

Build time is therefore about 53 ns × (4/3)·4^k per table, independent of genome
size, as long as the genome is large enough to fill the table (n ≫ 4^k). For a
3 Gb genome the table is full through k ≈ 14, so the times above hold. The one
thing not measured here is the TLB and page-walk cost of a 6× larger index,
which would add to both the build and the walk it saves.

A genome too small to fill the table builds faster, because empty subtrees are
skipped. But most of its deep entries are then empty or single rows, which is
what the `log4(genome size)` cap is for.

## Whole program

Depths were set at build time for this sweep (a `BATMAN_LOOKUPSIZE` define, since
replaced by `--lookup`), clamped to `k ≤ LH − 1`, so r24 runs k = 12 as 11.

The table below is search + table build in seconds, mean of 2 runs, pinned to
one core, on the 500 Mb genome:
- "before" is the build before this change (two arrays, 2·4^k).
- "default" is the new layout with the old 3/6 rule.

| reads, options | before | default | k=3 | k=6 | k=8 | k=10 | k=11 | k=12 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| r24 (2M) `-n 1 -m 0` | 18.43 | 18.57 | 18.30 | 17.91 | 17.04 | 16.01 | **15.85** | (=11) |
| r32 (2M) `-n 1 -m 0` | 11.31 | 11.11 | 11.25 | 11.01 | 10.31 | 9.25 | **9.10** | 10.98 |
| p100_mm1L (500k) `-n 1 -m 0` | 3.25 | 3.21 | 3.27 | 3.21 | 3.00 | **2.83** | 3.19 | 5.52 |
| p100_mm0 (500k) `-n 0 -m 0` | 1.85 | 1.82 | 1.86 | 1.83 | 1.80 | **1.73** | 1.83 | 3.58 |

Table build time on its own (both directions under `-n 1`):

| k | 8 | 10 | 11 | 12 |
|---|---:|---:|---:|---:|
| build (s) | 0.011 | 0.17 | 0.64 | 3.1 |

Search time on its own keeps falling with k. r32 runs 11.0 s at k = 6, 9.1 s at
10, 8.5 s at 11 and 7.9 s at 12: about 0.3 µs per read per level from k = 9,
which is two lookups at ~160 ns. Peak RSS (`-n 1`) is 868 MB at k ≤ 8, 884 MB at
k = 10, 932 MB at 11 and 1124 MB at 12.

## Break-even rule

Going from k−1 to k costs about 0.75 · 4^k · 80 ns per table. It saves about
170 ns per lookup, and there is one lookup per table per read. That evens out at
4^k ≈ 3 · reads, so the best k is round(log4(3 · reads)):
- 500k reads → 10.2 → **10** (measured best: 10)
- 2M reads → 11.2 → **11** (measured best: 11)
- 10M reads → 12.4 → 12
- 50M reads → 13.5 → 13; build ≈ 10 s, 1 GB

The constants are for this machine and the 500 Mb random genome. A repetitive
genome keeps wide ranges deeper, which makes walking cheaper per level (cached
blocks) but leaves more levels to walk. The rule should be rechecked on a real
genome before the cap is raised above 12.

## Correctness

- Output is byte-identical to the build before the change:
  - **150 runs at fixed depths:** builds at the default depth and at k = 3, 8,
    10 and 12, on 20k-read samples of r24, r32, p32_mm1L, p75_mm1L and
    p100_mm2, with `-n 0` and `-n 1`, `-m 0` and `-m 1`, with and without
    `--no-ref`.
  - **120 runs with `--lookup`:** default, 3, 10 and 12 on the same sets.
  - **60 runs on 1, 2, 3, 4 and 7 bp reads:** default and `--lookup 14`.
- `make test-located` and `make test-structured` pass, at the default depth and
  at depth 10.
- With the default `--lookup 6`, r24 takes 19.2 s against 19.9 s for the old
  build (3 runs each, machine not idle), so the change from k = 3 to 6 for short
  reads is within noise.
