# Prefetching for branching: overlapping independent cache misses

## Summary

Every search step (`Get_SARange`, or one `BWTOccValue` on a single-row path)
waits for a cache miss of ~155 ns on the 500 Mb index. The idea tested here:
find the ranges of several independent branches in bulk, prefetch their rows
for the next step together, and then compute them, so the misses overlap
instead of queueing. This is a stand-alone experiment (`branch_prefetch.cpp`);
nothing in `src/` was changed.

- **Overlapping works.** With 4 or more independent chains of steps in flight,
  a step costs 86 ns instead of 155 ns.
  - Plain round-robin over the chains gets there through out-of-order execution
    alone. An explicit `Prefetch_Row` adds nothing on top.
  - The literal form of the idea does best: compute all Starts, prefetch all
    End+1 rows, then compute the Ends. It reaches 74 ns at 8 chains and 66 ns
    at 32, 2.4× faster than one chain at a time.
- **But one read's branches give almost nothing to overlap.**
  - In the one-mismatch extension of the right half, a read has on average
    1.03–1.08 live branches per level. Its 16-base exact left half has usually
    narrowed the range to a single row already, so only 0.5–0.75 nodes per read
    branch at all.
  - Prefetching each child's rows when it is pushed (V1) changes nothing, within
    ±5%.
  - Stepping one read's branches level by level (V2) is slower, from step
    overhead (below).
- **Batching across reads gains 17–21%.** Stepping the live nodes of 16 reads
  together, after prefetching all their rows (V3), beats the real
  `Extend_Right`:

  | reads | real | 16 reads at a time |
  |---|---:|---:|
  | p32_mm1R | 2716 ns | 2151 ns (−21%) |
  | p32_mm0 | 2734 ns | 2236 ns (−18%) |
  | r32 | 2561 ns | 2133 ns (−17%) |
  | p100_mm1R | 8105 ns | 6636 ns (−18%) |

  That is with a step function that costs ~60% more per step than the real
  single-row loop (V2 vs V0). Part 1 suggests up to ~2× is available with a
  lean batched step.

**Verdict:**
- Prefetching sibling branches within a read is not worth doing on this data.
- The parallelism to exploit is across reads: interleave the searches of a
  batch of reads (~8–16) and keep each step lean. This applies to every pass,
  exact search included, since all of them are dependent chains of misses.
- It means restructuring the per-read loop into a batched state machine. It is
  not a local change to `Extend_Right`.

## Part 1: how far overlapping goes

The chains are 1M real r32 reads. Each starts from its depth-11 table range and
takes 8 `Get_SARange` steps along its own bases: 6.87M steps in total. Times
are ns per step, best of 3, pinned to one core. Every variant ends each chain
on the same range as the sequential run.

| chains in flight (B) | sequential | round-robin | + prefetch | split (Starts, prefetch End+1, Ends) |
|---:|---:|---:|---:|---:|
| 1 | 157 | 156 | 155 | 157 |
| 2 | 155 | 98 | 98 | 102 |
| 4 | 155 | 86 | 87 | 82 |
| 8 | 156 | 86 | 86 | 74 |
| 16 | 156 | 87 | 87 | 69 |
| 32 | 156 | 86 | 86 | 66 |

- **Round-robin flattens at ~86 ns from 4 chains.** `Get_SARange` is branchy and
  calls into the index library. Only about two to four steps' worth of
  instructions fit in the reorder window, and a prefetch issued in the same
  window cannot reach further.
- **The split loop keeps more misses in flight.** Its passes are short loops of
  a single kind of work, so it keeps improving up to 32 chains.

## Part 2: the real branching search

Each read's right half is extended with one mismatch allowed (`Count = 1`),
starting from its exact left-half range (left half searched with `Start_Range`
and `Search_Right`, depth-11 table). This is the 0|0 and 0|1 search in one pass.
`Use_Ref` is off, so only FM-index steps are timed.

| variant | what it does |
|---|---|
| V0 | the real `Extend_Right` |
| V0c | the same loop copied, single-row nodes walked by the real `Extend_Right_Singleton` |
| V1 | V0c + `Prefetch_Row` of each child's rows when it is pushed |
| V2 | one read level by level: every live node prefetched, then stepped once |
| V3 | V2 with the live nodes of 4 or 16 reads together |

Results, in ns per read, best of 3:

| reads | n | steps / read | branching nodes / read | live nodes per level (mean / max) | V0 | V0c | V1 | V2 | V2 + pf | V3 R=4 | V3 R=16 |
|---|---:|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|
| p32_mm1R | 1M | 16.6 | 0.47 | 1.04 / 15 | 2716 | 2702 | 2578 | 4251 | 4238 | 2212 | **2151** |
| p32_mm0 | 1M | 17.3 | 0.75 | 1.08 / 14 | 2734 | 2621 | 2629 | 4292 | 4293 | 2334 | **2236** |
| r32 | 731k | 16.1 | 0.64 | 1.07 / 14 | 2561 | 2476 | 2472 | 4046 | 4070 | 2195 | **2133** |
| p100_mm1R | 500k | 51.6 | 0.61 | 1.03 / 17 | 8105 | 8114 | 8113 | 13834 | 14023 | 6651 | **6636** |

Every variant finds the same hits as the real `Extend_Right` for every read:
the same (Start, End, mismatch position, mismatch base) set.

Notes:
- **Why V2 is slower.** V2 and V3 step every node through a general per-node
  function (`Step`), which pushes nodes through vectors. It costs about 255 ns
  per step against 160 ns for the real single-row loop, while doing the same
  number of rank calls (16.0 per read, counted with a linker wrap). With about
  one live node per read there is nothing to overlap, so V2 pays that overhead
  in full.
- **V3 hides the overhead.** It overlaps the misses of 16 reads, which more than
  pays for the slower step. A batched version with a lean step, like part 1's
  split loop, should get closer to the ~70 ns per step ceiling.
- **The data shapes this result.** On a more repetitive genome, or with shorter
  reads, left halves stay multi-row longer and more nodes branch. Within-read
  prefetching would then have more to work with. On this genome, branching
  happens at fewer than one node per read.
