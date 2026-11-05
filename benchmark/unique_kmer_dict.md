# Lossy dictionary of unique X-mers

## Summary

The idea is a dictionary of the genome's X-mers that occur exactly once.
- **Lossy:** some entries may be dropped, so a lookup can miss.
- **Exact keys:** a hit always confirms the X-mer occurs, once.

A search starts by looking up the read's first X bases. On a hit the read's range
at depth X is known to be a single row; on a miss the search runs as today. Two
uses of a hit were measured on the 500 Mb random genome, against today's start
(the depth-11 lookup table, then `Get_SARange` to depth X):

- **A: presence only, then a Start-only walk.** One rank per step instead of two,
  and End = Start since the X-mer is unique. **No gain:** within ±5% of today,
  and often slower.
  - The second rank was already almost free: below 128 rows `Get_SARange`
    counts instead of doing a second lookup, and above it the two lookups
    overlap.
  - Each step still waits for one cache miss.
  - The dictionary query adds another miss.
- **B: the entry also stores where the X-mer is** (its SA row, or its genome
  position). A hit skips the walk entirely.
  - **At X = 16 the cost of reaching depth 16 roughly halves:** r32
    968 → 482 ns, p100_mm0 1008 → 441 ns, p75_mm1L 948 → 545 ns (load 0.5).
    Reads whose 16-mer is repeated or absent pay the extra ~100 ns query.
  - **At X = 14 the gain is ~10%:** only 16–25% of reads hit, because most
    14-mers are repeated in 500 Mb.
- **Building by sorting is ~15× faster than walking the index.** Sorting the
  `.pac`'s X-mers took 10.6 s for X = 16; walking the FM index to depth 16 took
  159 s. Both give exactly the same unique sets. Sorting needs no index, so it
  fits naturally at index-build time.
- **Memory is the obstacle.** At X = 16 the 500 Mb genome has 428M unique
  16-mers. The dictionary takes 6.5 GB at load 0.5 and 3.3 GB at load 1.0
  (13% of keys dropped). For comparison, `batman_zero_min -n 1` peaks at
  868 MB with both indexes, the SA samples and the `.pac` loaded. Storing only part of each key could bring an entry to ~5 bytes
  (2 GB here).

**Verdict:**
- Version A is not worth pursuing.
- Version B works and saves ~500 ns per lookup at X = 16, and more in the real
  search: a hit also skips the single-row walk to a sampled SA row, which these
  numbers don't include.
- It costs memory comparable to or larger than the index itself.
- The next step, if wanted, is a whole-program prototype of B. It would store
  the genome position and pass hits to the existing `.pac` located path. It
  would run with the compact layout, measured on the full read sets.

## How it was measured

The benchmark is `benchmark/unique_kmer_dict.cpp`. It runs on the 500 Mb random
genome, pinned to one core, best of 3, starting from a depth-11 lookup table.

**Finding the unique X-mers:**
- **Index walk:** a depth-first walk over revfmi to depth 16. It records each
  range of exactly one row at depth X, together with its row.
- **Sorting:** every X-mer of the `.pac` packed as `key << 32 | position`,
  sorted with a 3-pass LSD radix sort on the key, keeping keys seen once.
- **Check on a 49 kb repetitive toy genome:** both methods give exactly the
  same unique X-mers as a Python brute-force count, at X = 12, 14 and 16. On
  500 Mb the two methods give identical key sets.

**Dictionary:**
- **Layout:** 64-byte buckets of 8 `{key, row}` slots. The bucket comes from a
  hash of the key, and the slot stores the whole key, so a hit is always exact.
- **Losses:** a key whose bucket is full is dropped. That loses 0.6–0.8% of keys
  at load 0.5, 6–7% at 0.8 and 13–14% at 1.0.
- **Correctness:** checked for every read. There were no wrong hits: each hit's
  X-mer walks to a range of exactly one row, the stored row equals the walked
  row, and A's Start-only walk ends on that row. Every stored key was found
  again with its row.

**Timing:** each read is timed two ways, as independent reads and with each read
waiting for the previous one. The two agree to within a few percent, so overlap
between reads does not distort the comparison.

## Results (ns per read to reach depth X, dependent timing)

Unique X-mers in the 500 Mb genome:

| X | unique | % of positions | sort build | index-walk build |
|---:|---:|---:|---:|---:|
| 12 | 0 | 0% | 5.3 s | (one walk to depth 16 for all X: 159 s) |
| 14 | 80M | 16% | 7.7 s | |
| 16 | 428M | 86% | 10.6 s | |

On 500 Mb every 12-mer occurs about 30 times, so a 12-mer dictionary is empty.

**X = 16, load 0.5 (6.5 GB, 0.6% dropped):**

| reads | hit | repeated | absent | today | A | B | query alone |
|---|---:|---:|---:|---:|---:|---:|---:|
| r32 (2M) | 63% | 10% | 27% | 911 | 939 | **507** | 121 |
| p100_mm0 (500k) | 85% | 14% | 0% | 1018 | 1030 | **439** | 185 |
| p32_mm2 (1M) | 28% | 4% | 68% | 849 | 853 | **631** | 111 |
| p75_mm1L (500k) | 53% | 8% | 39% | 944 | 925 | **535** | 115 |

**X = 16, load 1.0 (3.3 GB, 13% dropped):** B takes 590 / 494 / 694 / 590 ns for
the same four sets.

**X = 14, load 0.5 (1.2 GB):** hits on 16–25% of reads. Today takes ~665 ns, A
~690 ns and B 600–630 ns.

The "absent" share is the reads whose first 16 bases are not in the genome. It
is high for the mismatch sets, whose first 16 bases often contain a mismatch.

## Before taking B further

- **Genome.** A random genome overstates uniqueness. A real genome has more
  repeats, so it has fewer unique X-mers and a lower hit rate.
- **Key size.** On 3 Gb, X = 16 is short of unique: with 4.3G possible 16-mers
  for 3G positions, about half the positions would start a unique 16-mer.
  X ≈ 18 needs 36-bit keys.
- **What to store.**
  - The genome position plugs straight into the `.pac` located path (`Locate`,
    `Extend_Located`), which already works from a position.
  - The SA row is needed only under `--no-ref`.
- **Memory layout.** The bucket index already fixes part of each key, so a slot
  could keep only the remaining key bits plus a 32-bit position, about 5 bytes
  per entry. Even so, 3 Gb at X = 18 means about 10–15 GB.
- **Read halves.** X ≤ LH lets the `-n 1` halves use the dictionary too. Half of
  a 32 bp read is exactly 16 bases.
