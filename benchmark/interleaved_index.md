# Interleaving the Occ checkpoints with the BWT (the BWA layout)

BWA stores each block's Occ counts next to that block's BWT characters, so one rank reads
one cache line. This asks whether the BatMis index would gain from the same layout.
Stand-alone experiment (`interleaved_occ.cpp`); nothing in `src/` or `index/` was changed.

## Summary

- **Interleaving gains little here.** On a 500 Mb index, a 64-byte interleaved block makes
  a rank step 0–10% faster at up to 8 interleaved walks and 15–25% faster at 16–32, against
  the present layout with working prefetches. With threads the gain is within the noise.
  It needs a new index format and changes to every direct BWT read. Not worth doing now.
- **Why so little:** BatMis's two reads per rank (the BWT line and the count line) are
  independent, so they are already in flight together. Interleaving removes a miss that
  was not being waited for.
- **The experiment found a larger problem: the mappers' prefetches are not in the
  binaries.** GCC deletes every call to `Prefetch_Row` and `Mm_Prefetch` in `batmis1` and
  `batmis2`. Forcing the two functions inline (two attributes, no other change) cuts
  batmis2's batched search time by 33–45% on one core and 23–30% on 22 threads, with
  byte-identical output in all 76 runs compared. **Applied to both mappers on
  2026-11-15: see `prefetch_fix.md`.**

## The layouts

From `index/BWT.h` and `index/BWT.c`, for N bases:

| array | contents | size | 500 Mb | 3.1 Gb |
|---|---|---|---:|---:|
| `bwtCode` | 2 bits per base | N / 4 B | 125 MB | 775 MB |
| `occValue` | 4 × 16-bit counts per checkpoint, one every 256 rows | N / 32 B | 15.6 MB | 97 MB |
| `occValueMajor` | 4 × 32-bit counts every 65,536 rows | N / 16,384 B | 122 KB | 757 KB |

- An Occ takes the nearest checkpoint and counts at most 128 characters towards the row.
- A checkpoint interval is 64 bytes of `bwtCode`, and `bwtCode` is 64-byte aligned
  (checked at load), so the characters counted are in the row's own cache line.
- One Occ therefore reads one BWT line and one count line. The major table stays in cache.
  `BWTDecode` also loads 16 bytes of the next BWT line for a quarter of the rows, without
  using them.

Tested against it:

| name | layout | bits per base | size vs now | most characters counted |
|---|---|---:|---:|---:|
| LIB | the library's arrays and `BWTOccValue` | 2.25 | — | 128 |
| SEP | the same three arrays, counted with popcounts | 2.25 | — | 128 |
| B | 64-byte block: 4 × 32-bit counts + 192 bases | 2.67 | +19% | 96 |
| C | 64-byte block: 4 × 16-bit counts + 224 bases; major table every 256 blocks | 2.29 | +1.6% | 128 |

- In B and C the counts hold Occ at the block's character 96, and a row counts forward or
  backward from there.
- SEP uses the popcount code of B and C, so SEP against B/C compares layouts only.
- BWA's own block (4 × 64-bit counts + 128 bases) would be 4 bits per base, +78%. It was
  not tested: B is the same idea with counts half the size.

## Method

- **Index:** a 500 Mb random genome, forward index, built with `bwt_build_index`. The
  earlier write-ups' 500 Mb index was no longer on the machine.
- **Step:** an LF walk, the mappers' single-row step: read the row's BWT character, one
  Occ for it, go to that row. On this index every step lands on an unrelated row.
- **Chains:** 1–32 walks stepped round-robin, each next row prefetched a round ahead, as
  the batched lanes do.
- **Check:** all layouts equal `BWTOccValue` on 20M random rows plus the rows around `$`,
  both ends and the block edges, for all four bases (0 differences), and every timed run
  ends on the same rows.
- **Machine:** Core Ultra 7 165H, GCC 13.3 `-O3 -mpopcnt`. One-core runs are pinned to a
  P-core, load about 1, best of 3 × 20M steps. Thread runs are not pinned.
- `perf` is not usable on this machine (`perf_event_paranoid=4`), so there are no miss
  counts.

## Results: rank steps

One core, ns per step. LIB0 is the library with no prefetch, which is what the mappers
run today (next section).

| chains | LIB0 | LIB | SEP | B | C |
|---:|---:|---:|---:|---:|---:|
| 1 | 264 | 264 | 152 | 144 | 143 |
| 2 | 134 | 77 | 76 | 76 | 75 |
| 4 | 73 | 43 | 40 | 39 | 40 |
| 8 | 74 | 22.5 | 22.3 | 20.4 | 20.6 |
| 16 | 75 | 16.3 | 17.8 | 13.7 | 15.7 |
| 32 | 75 | 15.9 | 17.0 | 12.1 | 13.2 |

- **Up to 8 chains** the layouts are within 10% of each other. The time is the misses,
  and both layouts wait for one per step.
- **At 16–32 chains** B is 16–24% faster than LIB and C 4–17%. Here the limit is how many
  lines can be in flight, and the interleaved block asks for half as many.
- **One chain, 264 against 152 ns:** the library indexes the count array by the base
  (`occValue[... + character]`). In an LF walk the base comes from the BWT line, so the
  count line is not requested until the BWT line arrives, and the two misses run one after
  the other. SEP loads all four counts and then selects, so they overlap. This does not
  affect the mappers' searches, where the base comes from the read, nor `Sa_Decode`, which
  uses `BWTAllOccValue`.
- **Huge pages** (`--thp`, 590 of 730 MB backed by 2 MB pages) changed nothing: 16
  chains took 17.6 / 14.8 / 15.6 ns for LIB / B / C. TLB misses are not what
  interleaving would save.

Threads, 16 chains each, total million steps per second, best of 3:

| threads | LIB0 | LIB | SEP | B | C |
|---:|---:|---:|---:|---:|---:|
| 6 | 64 | 228 | 155 | 190 | 194 |
| 22 | 114 | 286 | 206 | 322 | 300 |

- B and C are 5–12% ahead of LIB on 22 threads and behind it on 6. The threads were not
  pinned on a CPU with three kinds of core, so I would not read more into this than "no
  clear gain".
- When everything is in cache (a 2 Mb index), the library's SSE2 count is faster than the
  popcount code: 10 ns against 13–16 ns a step.

## The prefetches that are not there

In the first runs LIB did not improve past 4 chains (74 ns) while SEP reached 17 ns on
the same arrays. The cause was in the compiler output, not the memory:

- `Lib::prefetch` was a copy of batmis2's `Prefetch_Row`. Its `Walk` contained no
  prefetch instruction at all.
- A function that only prefetches has no side effect that GCC models. When GCC does not
  inline it, it treats the call as dead and removes it. With two prefetches the function
  was inlined and kept; with three or four it was not, and all of them went.
- **The binaries of 2026-11-05 had the same fault.** `src/batmis2` held 7 prefetch
  instructions of its own, all the `saValue` prefetch of `Sa_Decode`'s walks (unrolled).
  There is no `Prefetch_Row` symbol and no call to one. `Mm_Prefetch`, which is called
  once per lane per round, is gone with everything in it, including its `.pac` and
  `saValue` prefetches. `src/batmis1` is the same (26 prefetch instructions in each
  binary, 19 of them inside the index library and htslib).
- This agrees with `branch_prefetch.md`, which found that an explicit `Prefetch_Row`
  "adds nothing on top" of round-robin and that steps flatten at about 86 ns. The batched
  search has been getting its overlap from out-of-order execution alone.

**Test on the mapper.** A copy of `batmis2.cpp` with `__attribute__((always_inline))` on
`Prefetch_Row` and `Mm_Prefetch`, built with the Makefile's flags next to an unchanged
copy (71 against 26 prefetch instructions). 1M reads sampled from the genome with 1%
substitutions, `--lookup 11 -m 1`, text output, both indexes with the reference.

Search seconds on one P-core, best of 2, unchanged → with prefetches:

| reads | `-n` | B = 1 | B = 4 | B = 16 | change at B = 16 |
|---|---:|---|---|---|---:|
| 100 bp | 0 | 1.86 → 1.87 | 1.04 → 0.85 | 1.19 → 0.71 | −40% |
| 100 bp | 1 | 3.02 → 2.98 | 1.93 → 1.50 | 1.90 → 1.28 | −33% |
| 100 bp | 2 | 3.77 → 3.73 | 2.33 → 1.82 | 2.32 → 1.53 | −34% |
| 100 bp | 3 | 3.77 → 3.72 | 2.43 → 1.86 | 2.37 → 1.59 | −33% |
| 32 bp | 0 | 1.48 → 1.45 | 0.84 → 0.62 | 0.83 → 0.50 | −39% |
| 32 bp | 1 | 2.29 → 2.19 | 1.42 → 1.06 | 1.41 → 0.89 | −37% |
| 32 bp | 2 | 11.8 → 11.0 | 8.9 → 6.3 | 8.8 → 5.3 | −40% |
| 32 bp | 3 | 34.4 → 35.3 | 30.5 → 19.8 | 31.0 → 17.0 | −45% |

`-t 22`, `-n 2`, median of 3:

| reads | B = 4 | B = 16 |
|---|---|---|
| 100 bp | 0.315 → 0.229 (−27%) | 0.321 → 0.224 (−30%) |
| 32 bp | 1.367 → 1.048 (−23%) | 1.429 → 1.058 (−26%) |

- The output files were byte-identical in all 64 one-core runs and all 12 threaded ones.
- B = 1 is unchanged, as it should be: one read has nothing to prefetch ahead of.
- `--batch auto` gave the same gains (for example 100 bp `-n 2`: 2.42 → 1.61 s).
- With working prefetches B = 16 beats B = 4 everywhere. Without them the two are level.
  The best-B findings in `batched_reads.md` and `threads.md` would need redoing.
- The `test` suites were not run on the changed copy.

## What I would do

1. **Fix the prefetches.** Done in both mappers, with a build check
   (`prefetch_fix.md`).
2. **Leave the index layout alone for now.** Repeat this experiment only
   if rank steps at high B still dominate. Two cases this run did not cover could change
   the answer:
   - a human-size genome, where the count table is 97 MB and no part of it stays in L3
     (at 500 Mb its 15.6 MB partly fits the 24 MB L3, which favours the present layout);
   - a real, repetitive genome, where ranges stay wide and steps are `Get_SARange` pairs
     instead of single-row walks.
3. If it is ever done, layout C is the one: +1.6% memory, one line per rank, and the
   existing files can be converted at load in one pass.

## Reproducing

```sh
g++ -std=c++17 -fpermissive -w -O3 -mpopcnt -pthread -Iindex \
    -o interleaved_occ benchmark/interleaved_occ.cpp index/libfmindex_core.a
taskset -c 2 ./interleaved_occ <index prefix>                       # one core, chains 1-32
taskset -c 2 ./interleaved_occ <index prefix> --thp                 # with huge pages
./interleaved_occ <index prefix> --threads 6,22 --chains 1,4,16 --steps 10000000
```
