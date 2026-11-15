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
  every case tried here; before, the two were level and threads preferred B = 2–4. The
  default changed with it (next section).
- **Batched timings everywhere** (`recap.md`, the `-n` write-ups): expect about a third
  less search time at B ≥ 4.
- **Experiments judged under the old behaviour** could come out differently now that
  memory waits are shorter and instructions matter more: the inline Occ
  (`inline_occ.md`), checkpoint spacing (`occurrence_checkpoint_spacing.md`), range
  chains, and the interleaved index (`interleaved_index.md` already compares against
  working prefetches).
- **batmis1 against batmis2:** both got the fix, so the comparison is still like for
  like. On these runs batmis1 is 4–11% ahead in search time at `-m 1`, before and after.

## The batch size after the fix

`--batch auto` was the default. It still chose well, but its trial became the expensive
part, so the default is now a fixed size and auto is a check on it.

**Fixed sizes.** batmis2, 1M reads, one P-core, search seconds (a build was running
alongside, so these are noisier than the tables above):

| reads | `-n` | B = 1 | 2 | 4 | 8 | 16 | 32 | 64 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 100 bp | 0 | 1.72 | 1.02 | 0.75 | 0.68 | 0.67 | 0.68 | 0.69 |
| 100 bp | 2 | 3.54 | 2.35 | 1.79 | 1.52 | 1.84 | 1.55 | 1.55 |
| 32 bp | 0 | 1.88 | 0.95 | 0.71 | 0.58 | 0.55 | 0.58 | 0.58 |
| 32 bp | 2 | 12.6 | 9.1 | 7.4 | 6.3 | 6.3 | 6.4 | 6.7 |

The time is flat from 8 to 64, so a size anywhere in that range is within a few percent
of the best. (The 1.84 at 100 bp `-n 2`, B = 16 is the build's noise: 1.24–1.25 s in the
runs below.)

**The old trial.** Auto tried B = 1, 2, 4, 8, 16 three times each on 4,096-read chunks:
61k reads per thread, 1.35M at `-t 22`. Without prefetches the sizes differed little and
that cost nothing. With them the B = 1 and B = 2 chunks run 1.5–2.5× slower than the
best:
- **`-t 22`, 1M reads, 100 bp `-n 2`:** no thread finished its trial; 0.225–0.230 s
  against 0.188–0.202 s at fixed B = 4–32.
- **`-t 22`, 4M reads:** 0.99–1.04 s against 0.88–0.98 s fixed.
- **One thread:** within 2–5% of the best fixed size.

**Now:**
- **Default:** `--batch 16` (`DEFAULT_BATCH`), with no trial.
- **`BATMIS_BATCH`:** an environment variable, `auto` or 1–64, that sets the size when
  `--batch` is not given. For a machine where 16 turns out not to suit.
- **`--batch auto`:** tries B = 4, 8, 16 only (nine chunks a thread), and prints its
  measurements to stderr, with or without `--timing`:

  ```
  batch auto: ns/read 4:1529 8:1277 16:1236; chosen 16; default 16 is +0.0% against the fastest (16)
  batch auto: ns/read 4:10054 8:9372 16:9894; chosen 4:5 8:12 16:5 threads; default 16 is +5.6% against the fastest (8)
  ```

  The first is one thread; the second is `-t 22`, where each time is the median over the
  threads and `chosen` counts the threads that settled on each size. An input too short
  for the trial still prints the times it has.

**Checking a new machine.** Run the workload a few times with `--batch auto` and read the
last figure. If the default is repeatedly more than a few percent from the fastest, set
`BATMIS_BATCH` there. If the fastest is at an end of the range (4 or 16), time
`--batch 2` or `--batch 32` by hand as well. The threaded figures move by 5% or so from
run to run, so one run is not enough. Only this machine (Core Ultra 7 165H) has been
measured.

The output does not depend on the batch size: the default, `--batch 16`, `--batch auto`,
`BATMIS_BATCH=4` and `-t 22 --batch auto` gave byte-identical files.
