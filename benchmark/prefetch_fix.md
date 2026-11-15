# The prefetches the compiler dropped (fixed 2026-11-15)

## Summary

- Until 2026-11-15 the `batmis1` and `batmis2` binaries did not contain the prefetches
  that the batched search issues. GCC had removed every call to `Prefetch_Row` and
  `Mm_Prefetch`. The output was right; only the batched search was slower than it
  should have been.
- Both functions are now declared `PREFETCH_ONLY` (`src/fm_occ.hpp`), which forces them
  inline. Batched search is 35–42% faster on one core. Output is unchanged.
- `tests/check_prefetch.sh`, run by `make` and `make test`, fails the build if the
  prefetches go missing again.
- **Every batched timing in the other write-ups was measured without prefetches.** See
  "What this changes" below.

## Cause

- `__builtin_prefetch` has no side effect that GCC models. A function made only of
  prefetches therefore looks free of effects, and a call to it whose result is unused is
  dead code. GCC removes such a call when it does not inline the function. Nothing is
  reported.
- `Prefetch_Row` (four prefetches, many call sites) and `Mm_Prefetch` (a switch of them,
  called once per lane per round) were both `inline`, and GCC 13.3 `-O3` chose not to
  inline either. The binaries had no `Prefetch_Row` symbol and no call to one.
- What remained: 7 prefetch instructions in each mapper's own code, all the `saValue`
  prefetch written directly in `Sa_Decode`'s walk loop.
- It was found with `benchmark/interleaved_occ.cpp`, where a copy of `Prefetch_Row` did
  not speed up the library's walks and a two-prefetch version, small enough to be
  inlined, did (`interleaved_index.md`).
- `branch_prefetch.md` had seen the symptom: an explicit `Prefetch_Row` "adds nothing on
  top" of round-robin, and steps flatten near 86 ns. That flattening was out-of-order
  execution alone.

## The fix

```cpp
// src/fm_occ.hpp
#define PREFETCH_ONLY __attribute__((always_inline)) inline

// src/batmis1.cpp, src/batmis2.cpp
PREFETCH_ONLY void Prefetch_Row(const ActiveBWT* fmi, unsigned long Index) { ... }
PREFETCH_ONLY void Mm_Prefetch(const Mm_Lane& L) { ... }
```

Inlined, the prefetches sit in callers that have other effects (the lane loop,
`Sa_Decode`), and stay. The mappers' own prefetch instructions went from 7 to 51
(`batmis1`) and 52 (`batmis2`).

## Keeping them

A dropped prefetch changes no output and fails no test, so it has to be checked in the
binary. `tests/check_prefetch.sh <binary> <source>` checks:

1. **Source:** every function named `*Prefetch*` is declared `PREFETCH_ONLY`.
2. **Binary:** the mapper's own code (its anonymous namespace, so not the index library
   or htslib) holds at least 20 prefetch instructions: 7 without the fix, about 50 with.

`make` (target `check-prefetch`, part of `all`) and `make test` run it for both mappers.
It fails on the pre-fix binaries.

Rules for new code:
- A function that only prefetches must be `PREFETCH_ONLY`, and so must a function that
  only calls such functions. Name it `*Prefetch*` so that the check sees it.
- A prefetch written directly in a function that does other work is safe.
- After adding prefetch sites, look at the count the check prints: it should rise.
- The floor of 20 is for the Makefile's `-O3`. For a build where the count is
  legitimately lower (`-O0`, for instance), set `PREFETCH_MIN`, or `PREFETCH_MIN=0` to
  skip the binary check.
- With another compiler or GCC version, run the check before trusting any timing.

## Measurements

500 Mb random genome, 1M reads sampled from it with 1% substitutions, `--lookup 11 -m 1
--batch 16 -t 1`, pinned to one P-core (Core Ultra 7 165H). Search seconds, the binaries
of 2026-11-05 → the fixed ones, two runs each (the faster shown; the pairs agree within
2%).

| reads | `-n` | batmis1 | batmis2 |
|---|---:|---|---|
| 100 bp | 0 | 0.886 → 0.557 (−37%) | 0.945 → 0.578 (−39%) |
| 100 bp | 2 | 1.931 → 1.246 (−35%) | 2.011 → 1.313 (−35%) |
| 32 bp | 0 | 0.798 → 0.468 (−41%) | 0.820 → 0.494 (−40%) |
| 32 bp | 2 | 8.139 → 4.767 (−41%) | 8.749 → 5.299 (−39%) |

The fuller grid, on a scratch build of batmis2 before the fix went into `src/`
(`interleaved_index.md`):
- **B = 1:** unchanged (one read has nothing to prefetch ahead of).
- **B = 4:** 19–35% less search time.
- **B = 16:** 33–45% less, for `-n 0..3` at 100 and 32 bp.
- **`-t 22`, `-n 2`:** 23–30% less.

**Same output:**
- text output byte-identical, old against new, in all 16 runs above;
- SAM at `-t 22 -n 2`, 1M reads: identical for both mappers apart from the command line
  in `@PG`;
- `make test` passes (both suites and `test_batmis2.sh`).

## What this changes in the other write-ups

Their numbers stand as measured, but they were measured without prefetches. Not redone:

- **Best B** (`batched_reads.md`, `threads.md`): with prefetches B = 16 beats B = 4 in
  every case tried here; before, the two were level and threads preferred B = 2–4.
  `--batch auto` picks B by timing, so it follows by itself.
- **Batched timings everywhere** (`recap.md`, the `-n` write-ups): expect about a third
  less search time at B ≥ 4.
- **Experiments judged under the old behaviour** could come out differently now that
  memory waits are shorter and instructions matter more: the inline Occ
  (`inline_occ.md`), checkpoint spacing (`occurrence_checkpoint_spacing.md`), range
  chains, and the interleaved index (`interleaved_index.md` already compares against
  working prefetches).
- **batmis1 against batmis2:** both got the fix, so the comparison is still like for
  like. On these runs batmis1 is 4–11% ahead in search time at `-m 1`, before and after.
