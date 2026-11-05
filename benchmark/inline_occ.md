# The inline Occ in batman_zero_min

`bifmindex/` replaced the library's table-based Occ (`BWTOccValue` / `BWTAllOccValue`,
which decode through `BWTDecode` / `BWTDecodeAll`) with an inline count. That saved 10% of
instructions and 6–15% of search time there (`bifmindex/README.md`, "Inline Occ and the
batched profile"). The same code, in `src/fm_occ.hpp`, was tried in `batman_zero_min`.
**Result: identical output, no speed-up; the library Occ stays the default.**
`-DBATMAN_INLINE_OCC` (`make batman_zero_min_inlineocc`) selects the inline one.

## Method

- **Occ:** the nearest checkpoint (`occValueMajor` plus the 16-bit minor count), plus or
  minus a direct count of the at most 128 characters between it and the row:
  - one base: `popcnt(Y & Y >> 1 & 0x55…)` with `Y = ~(word ^ c·0x55…)`;
  - all bases: bit-planes CT, GT and T, three popcounts per word.
- **Call sites:** every `BWTOccValue` (exact search, the mismatch lane walks, singleton
  walks) and `BWTAllOccValue` (`Branch_Detect`, `Sa_Decode`). `Count_Range` calls
  `fm_occ::Count_Codes` in both builds.

## Correctness

- **`make test-occ`:** every row of a 40-contig toy index (172,644 rows) and a 1.2 Mb
  genome (2,400,004 rows), both indexes, equals the library's.
- **500 Mb genome:** 20,563,202 rows (random, around the `$` row and around checkpoints),
  0 differences.
- **Output:** byte-identical to the pre-change binary in 100 configurations: `-n 0..5` ×
  text/SAM × (`--batch 1 -t 1`, `-t 22`) × ref/`--no-ref` × 100/32 bp mixed reads, plus
  `-n 6/7` on 20k reads.
- **Suites:** `make test-structured`, `make test-sam`, `make test-located` pass.

## Timing

500 Mb genome, `-t 1 --lookup 11 -m 1`, one P-core (`taskset -c 1`); 2M reads for
`-n 0/1`, 500k for `-n 3/5`. Search seconds, round 1 / round 2 (alternating order). The
machine was shared (load 5–18), so differences under ~10% are noise.

| run | library Occ | inline Occ |
|---|---|---|
| 100 bp `-n 0` B = 1 | 14.3 / 13.8 | 14.0 / 14.4 |
| 100 bp `-n 0` B = 8 | 9.1 / 8.0 | 9.0 / 8.3 |
| 100 bp `-n 1` B = 8 | 14.8 / 14.8 | 17.2 / 15.8 |
| 100 bp `-n 3` B = 8 | 4.3 / 4.9 | 4.4 / 4.8 |
| 100 bp `-n 5` B = 8 | 5.3 / 5.1 | 4.7 / 4.8 |
| 32 bp `-n 0` B = 1 | 9.4 / 9.4 | 9.3 / 9.2 |
| 32 bp `-n 0` B = 8 | 5.9 / 6.2 | 6.5 / 6.3 |
| 32 bp `-n 1` B = 8 | 11.4 / 9.9 | 12.8 / 11.1 |
| 32 bp `-n 3` B = 8 | 7.2 / 6.4 | 7.6 / 6.9 |
| 32 bp `-n 5` B = 8 | 7.5 / 6.9 | 7.0 / 6.7 |
| 100 bp `-n 3`, `-t 22`, total | 4.17 | 4.42 |

`perf stat` (user space):

| run | instructions, library → inline | cycles, library → inline |
|---|---|---|
| 100 bp `-n 0` B = 8 | 1.85e10 → 1.74e10 (−6%) | 1.45e10 → 1.44e10 |
| 100 bp `-n 3` B = 8 | 1.17e10 → 1.06e10 (−9%) | 1.0e10 → 1.0e10 |
| 32 bp `-n 1` B = 8 (×2) | 2.26e10 → 2.11e10 (−7%) | 1.83–1.90e10 → 1.94–1.97e10 |

At 32 bp `-n 1`, LLC load misses also fell (7.1e7 → 6.7e7), yet search took 10.1–10.3 s
against 11.3–11.6 s.

## Conclusion

- `batman_zero_min` is bound by memory latency: the inline count removes 6–9% of the
  instructions but no cycles.
- `-n 1` is consistently 10–15% slower with it, in every pair. The cause was not
  profiled; fewer instructions and misses with more cycles suggests a longer dependent
  chain in the one-base Occ of the mismatch lane walks.
- The single-index version is the fixed baseline, so it keeps the library Occ.
  `bifmindex/` keeps the inline Occ, where Occ is a larger share of the work.
