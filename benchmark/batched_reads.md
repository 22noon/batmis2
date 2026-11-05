# Batched reads: interleaving the searches of several reads

> **Since the program refactor** (`five_mismatch.md`, "Program builder"): the hand-written
> per-read searches named below (`One_Mismatch_Search` … `Five_Mismatch_Search`,
> `Extend_Windows`, `Switch_Index`, the `Program_*` arrays, the `SEG_*` / `WIN_*` enums) are
> gone. `Mm_Builder` builds the same programs step for step, and `--batch 1` runs them
> on one lane of the engine.

## Summary

Each search step waits for one cache miss of ~155 ns, and a read's steps depend
on each other. `benchmark/branch_prefetch.md` showed that the independent work
is across reads, not within one. `batman_zero_min -n 0 --batch B` now searches B
reads at a time ("lanes"). Each round:
1. prefetches what every lane's next step reads;
2. steps every lane once;
3. hands a lane the next read when its read finishes.

Each lane follows `Search_Right` and `Search_Right_Singleton` exactly, so the
output is byte-identical to `--batch 1`.

- **Exact search is 1.7–2.6× faster.** The table below is input + search time on
  the 500 Mb genome, mean of 2 runs, pinned to one core. Encoding the reads moves
  from search to input in batched mode, so the two are added together.

  | reads | `--lookup` | `--batch 1` | `--batch 16` | speed-up |
  |---|---:|---:|---:|---:|
  | r32 (2M) | 6 | 5.10 s | 3.00 s | 1.70× |
  | r32 | 11 | 4.00 s | 2.19 s | 1.83× |
  | r32 `--no-ref` | 11 | 5.26 s | 2.54 s | 2.07× |
  | p100_mm0 (500k) | 6 | 1.88 s | 1.13 s | 1.66× |
  | p100_mm0 `--no-ref` | 6 | 7.69 s | 3.04 s | 2.53× |
  | p75_mm0 (500k) | 11 | 1.48 s | 0.81 s | 1.83× |
  | p75_mm0 `--no-ref` | 11 | 5.20 s | 2.08 s | 2.50× |

  - **Combined with a deeper table:** r32 goes from 5.10 s (today's default) to
    2.19 s, 2.3× overall.
  - **The gain is largest under `--no-ref`.** Long single-row walks there are
    exactly the dependent chains that batching overlaps.
- **B = 4 is enough.** Times are flat from 4 to 32 lanes, about the same point
  where part 1 of the prefetch experiment flattened.
- **Splitting wide steps** (all Start ranks first, then End ranks) gained
  nothing in the full program, so it was dropped.
- **Coroutines are a viable way to batch `-n 1`.**
  - The same batched exact search was written as C++20 coroutines
    (`benchmark/batch_coroutine.cpp`): each read is plain sequential code that
    prefetches, suspends, and is resumed round-robin.
  - It gives identical results and costs only 0–3% more than the hand-written
    lanes with the reference, and 4–8% more under `--no-ref`.
  - The one-mismatch search, with three passes, depth-first stacks, deferral
    and early exit, could therefore keep its current structure and gain
    suspension points, instead of becoming a hand-written state machine.

`--batch` now defaults to `auto`, which picks B on the fly ("Choosing B" below). It
covers every `-n` and SAM output.

## How it works

The exact-search lanes described in this section's first version (`Batch_Exact`
with its own `LANE_*` states) were later merged into the one engine described
under "One engine: the pass program" below. The timings above were measured with
that first version, and the merged engine times the same.

`src/batman_zero_min.cpp`:
- **`Run_Batched`:** used by `main` whenever `--batch` > 1.
  - It reads records in blocks of 4096 and encodes them. Invalid reads and reads
    of the wrong length get their `F 0 0 0` line in place, as before.
  - It searches the block, then writes the lines in input order.
- **Exact search:** each lane follows `Search_Right` and its single-row walk,
  one memory access per step:

  | step | prefetched before the step | does |
  |---|---|---|
  | wide | rows Start and End+1 | `Get_SARange` |
  | single row | row Start | one single-row step (`Read_Base`, one `BWTOccValue`) |
  | locate | `saValue[Start/8]` | `Locate` |
  | compare | the `.pac` bytes of the rest of the read | `Rest_Mismatches` |

- **`Packed` parameter:** `Pack_Read`, `Mismatch_Bits`, `Genome_Mismatches`,
  `Rest_Mismatches`, `Search_Located`, `Extend_Located` and `Record_Mismatches`
  take a `Packed` read, defaulting to the global `Read_Packed`, so existing
  callers are unchanged.

## Lanes vs coroutines (search only, ns per read, best of 3)

| reads | lanes B=1 | lanes 4 | lanes 16 | coroutines B=1 | coroutines 4 | coroutines 16 |
|---|---:|---:|---:|---:|---:|---:|
| r32 | 2349 | 1380 | 1303 | 2349 | 1320 | 1324 |
| p100_mm0 | 3393 | 1813 | 1853 | 3398 | 1876 | 1856 |
| r32 `--no-ref` | 2989 | 1550 | 1499 | 3007 | 1658 | 1623 |
| p100_mm0 `--no-ref` | 14703 | 5692 | 5677 | 14775 | 5959 | 5977 |

The coroutine frames come from a free list, so starting a read allocates
nothing. One pitfall turned up while building this: a suspended coroutine keeps
pointers into its slot. Slots must not be moved, or must point at storage that
stays put. That bug showed up as wrong results on the last few reads of a set.

## Verification

- **Output byte-identical to the build before batching:** 240 runs. That covers
  `--batch` 1, 4, 16 and 32, `--lookup` 6 and 11, `-m 0` and `-m 1`, with and
  without `--no-ref`, on 20k-read samples of r24, r32, p32_mm1L, p75_mm1L,
  p100_mm2 and p100_mm0. It includes `-n 1` with `--batch 16`, which ignores
  the flag.
- **1–7 bp reads:** output identical.
- **Tests:** `make test-located` now also compares `-n 0` with `--batch 16`
  against one read at a time, byte for byte. It passes, as does
  `make test-structured`.
- **Coroutines vs lanes:** identical ranges, positions and found flags on every
  read of r32 and p100_mm0, with and without the reference.

## One-mismatch search (`-n 1 --batch B`)

`Batch_One_Mismatch` runs `One_Mismatch_Search` as a per-read state machine, B
reads at a time. It uses hand-written lanes, chosen over coroutines.

**Tasks.** A lane runs one task at a time:
- **Searches:** a search of a half (`Search_Right` / `Search_Left` and their
  single-row walks);
- **Extensions:** a depth-first extension (`Extend_Right` / `Extend_Left`,
  their single-row walks and the located path), with the read's own stack.

**Steps.** Each step does one memory access, and is prefetched in the round
before:

| step | what it does |
|---|---|
| pop | take the next node of the depth-first stack |
| wide / threshold | `Get_SARange` |
| branch | `Branch_Detect` |
| single row | one step of the single-row walk |
| locate | the SA sample |
| located compare | the `.pac` |

The steps are the original loop bodies cut at each access, in the same order.

**Pass driver.** `Mm_Advance` runs:
1. 0|0 search, then 0|0 extension;
2. the 0|1 resume, one task per deferred parent;
3. 1|0 search, then 1|0 extension;

with the `MaxHits` checks where `One_Mismatch_Search` has them. So each read's
hits, deferred ranges and their order are unchanged.

**Per-read state:**
- `Branch_Detect`'s globals are written and read within one step;
- the packed read is per lane (`Search_Located`, `Extend_Located` and
  `Record_Mismatches` now take `Packed` too);
- `main` and the block runner share `Format_One_Mismatch` for the output lines.

**Timing.** Input + search seconds, mean of 2 runs, on the 500 Mb genome. The
machine had a load average of about 2–3, so expect ±5% noise.
- `old` is today's per-read code (`--batch 1`).
- `lanes at 1` is the state machine with a single lane (a scratch build), to
  show its own overhead.

| reads, options | old | lanes at 1 | `--batch 4` | `--batch 16` | speed-up |
|---|---:|---:|---:|---:|---:|
| r24 (2M), lookup 6, `-m 0` | 19.34 | 19.82 | 14.82 | 14.68 | 1.32× |
| r24, lookup 11, `-m 0` | 16.71 | 17.22 | 12.51 | 12.27 | 1.36× |
| r32 (2M), lookup 6, `-m 0` | 12.42 | 12.71 | 8.18 | 8.29 | 1.52× |
| r32, lookup 11, `-m 0` | 9.83 | 9.10 | 5.77 | 5.60 | 1.75× |
| r32, lookup 11, `-m 1` | 6.57 | 6.32 | 3.58 | 3.50 | 1.87× |
| r32, lookup 6, `--no-ref` | 15.84 | 15.79 | 9.55 | 9.59 | 1.66× |
| p32_mm1L (1M), lookup 11, `-m 0` | 3.89 | 3.89 | 2.31 | 2.28 | 1.71× |
| p75_mm1L (500k), lookup 11, `-m 0` | 2.61 | 2.53 | 1.47 | 1.48 | 1.78× |
| p75_mm1L, lookup 6, `--no-ref` | 7.65 | 7.63 | 3.80 | 3.80 | 2.01× |
| p100_mm1L (500k), lookup 11, `-m 0` | 2.93 | 2.90 | 1.72 | 1.70 | 1.73× |
| p100_mm1L, lookup 6, `--no-ref` | 10.70 | 10.29 | 4.95 | 4.91 | 2.18× |

- **The gain is 1.3–2.2×.** It is flat from 4 lanes, as for the exact search,
  and largest under `--no-ref`, where single-row walks are long.
- **r24 gains least (1.3–1.5×).** Its 12-base halves branch at many nodes
  (7 multi-row exact-remainder steps per read; see `sampled_end.md`). Its
  steps are therefore more often `Branch_Detect`, which already overlaps two
  rank lookups.
- **The lanes themselves cost nothing measurable.** At one lane they are
  within ±5% of today's code.
- **With the deeper lookup table:** r32 `-n 1` goes from 12.42 s (lookup 6, one
  read at a time) to 5.60 s (lookup 11, `--batch 16`), 2.2× overall.

**Verification.**
- **Output byte-identical to the build before batching:** 308 runs.
  - `-n 1`, `-m 0` and `-m 1`, `--batch` 1, 4, 16 and 32, `--lookup` 6 and 11,
    with and without `--no-ref`.
  - 20k-read samples of r24, r32, p32_mm1L, p75_mm1L, p75_mm1R, p100_mm1L and
    p100_mm2, and of the repetitive genome's x32_mm0 and x100_mm1.
  - The 1–7 bp reads.
- **`make test-located` field-by-field check.**
  - For every read, the lanes (B = 1 and 16, reference off and on) must give
    the same three result lists as `One_Mismatch_Search`: hits and deferred
    ranges compared field by field, hit counts, and the exact hit count.
  - It also compares the output of `-n 0` and `-n 1` at `--batch 16` against
    `--batch 1`, byte for byte.
  - Dropping one deferral from the lane engine made it fail on every
    affected read.
- `make test-structured` passes.

## Mixed reads

**The question.** Reads in one run differ widely in work:
- an exact read finishes after 0|0, or right after its exact search with `-m 1`;
- one-mismatch reads run all three passes;
- reads with no hit die early.

Lanes are refilled independently as each read finishes, with no barrier per
round, so a mixture should cost what its parts cost.

**The test.** Mixtures were built from the 500 Mb sets:
- **mix32:** p32_mm0, mm1L, mm1R and mm2, 250k each (1M reads);
- **mix100:** the p100 equivalents, 125k each;
- **xmix100:** the repetitive genome's x100_mm0 and x100_mm1, 200k each.

Each was run interleaved and shuffled, and compared with the same reads run as
separate pure files. Settings: `-n 1`, `--lookup 11`, reference on; input +
search seconds.

| reads | `-m` | mix, batch 1 | sum of parts, batch 1 | mix, batch 16 | sum of parts, batch 16 |
|---|---:|---:|---:|---:|---:|
| mix32 | 0 | 3.85 | 3.84 | 2.26 | 2.26 |
| mix32 | 1 | 3.05 | 2.96 | 1.70 | 1.74 |
| mix32 shuffled | 0 | 3.89 | 3.84 | 2.41 | 2.26 |
| mix100 | 0 | 2.75 | 2.96 | 1.55 | 1.76 |
| mix100 shuffled | 0 | 3.14 | 2.96 | 1.71 | 1.76 |
| xmix100 | 0 | 2.95 | 3.07 | 1.78 | 1.79 |
| xmix100 | 1 | 1.97 | 1.80 | 1.11 | 1.12 |

- **A mix costs the sum of its parts**, within ±10% run-to-run noise, batched
  or not. Interleaving or shuffling makes no difference.
- **Occupancy stays at 15.7–15.97 of 16 lanes** (`-DBATMAN_BATCH_STATS`,
  `--batch 16`).
  - Reads with long single-row walks or deep extensions do not hold the other
    lanes back.
  - The only idle time is the drain at the end of each 4096-read block: 0.7%
    of rounds for mix32, 1.5% for mix100, 2.3% for xmix100.
- **Larger blocks don't help.** 65536-read blocks bring occupancy to 15.98 but
  are no faster, and slightly slower in most runs, because the per-block
  result arrays grow 16-fold. 4096 stays.
- **Where the steps go:** per read, mix32 at `--lookup 11`, `--batch 16`:

  | stage | steps per read |
  |---|---:|
  | 0|0 search | 4.0 wide + 0.6 single-row + 0.25 locate/compare |
  | 0|0 extension | 4.3 |
  | 0|1 resume | 2.4 |
  | 1|0 search | 4.8 |
  | 1|0 extension | 5.5 |
  | total | 21.8 |

  The two exact half searches are 44% of the steps, and the extensions the
  rest.

## One engine: the pass program

Both searches now run on one engine. Each read runs a short pass program (`Mm_Program`) of five operations:

| operation | meaning |
|---|---|
| `OP_SEARCH(segment, direction)` | exact search of a read segment: forward table + revfmi to the right, backward table + fwfmi to the left, then wide steps, single-row walk and located shortcut |
| `OP_EXTEND(segment, direction, count, list)` | one depth-first extension of the last search's range over a segment, with `count` mismatches, hits and deferred ranges into result list `list` |
| `OP_RESUME(from, segment, direction, count, list)` | an `OP_EXTEND` for each range list `from` deferred, in order |
| `OP_COUNT_EXACT(list)` | the read's exact hit count |
| `OP_STOP_IF_MAXHITS` | the early exit |

The two programs:
- **`-n 0`:** `SEARCH(read, right)`.
- **`-n 1`:** `SEARCH(left half, right) → EXTEND(right half, right, 0, list 0) →
  COUNT_EXACT(list 0) → RESUME(list 0, right half, right, 1, list 1) →
  STOP_IF_MAXHITS → SEARCH(right half, left) → EXTEND(left half, left, 1,
  list 2)`.

**`MaxHits` budgets.** An extension that stops at `MaxHits` gets
MaxHits − (exact hits + the other lists' hits) as its budget, counted against
its own list. That reproduces `One_Mismatch_Search`'s two budgets exactly.

**Checked.**
- **Output:** byte-identical to the build before batching in 852 runs.
  - `-n 0` and `-n 1`, `-m 0` and `-m 1`, `--batch` 1, 4, 16 and 32,
    `--lookup` 6 and 11, both reference modes.
  - 13 read sets: the 500 Mb, repetitive and mixed samples.
  - The 1–7 bp reads.
- **Tests:** `make test-located` (field by field) and `make test-structured`
  pass.
- **Timing:** within noise of the two separate engines it replaced: r32 `-n 0`
  2.17–2.47 s against 2.19–2.22 s, `-n 1` 5.27–5.34 s against 5.32–5.50 s, and
  p100 within ±4%.

### Higher mismatch counts

A k-mismatch search in the same style splits the read into segments and uses the
pigeonhole rule. If the read has at most k mismatches, at least one of k + 1
segments matches exactly. Each segment in turn:
- is searched exactly (`OP_SEARCH`);
- has its neighbours extended with a growing budget (`OP_EXTEND` with count
  0, …);
- has its deferred ranges resumed with a larger count (`OP_RESUME`);
- with early exits between passes (`OP_STOP_IF_MAXHITS`).

The engine already supports everything per task: direction, segment offset,
`Count`, deferral with `Forbidden_Base`, located hits, and `MaxHits` budgets.
What a k-mismatch program would still need:
- **Segments beyond the two halves.** `Mm_Segment` knows 0 / 1 / 2. More
  segments need a table of (offset, length) per scheme.
- **More result lists than three.** `MM_LISTS`, plus `Run_Batched`'s output
  formatting by list and direction.
- **An extension that continues across several segments**, with a different
  budget per segment (as BatMis's multi-mismatch schemes do). That is either a
  chain of `OP_EXTEND`s, each starting from the previous one's hits (a new
  "from list" form of `OP_EXTEND`), or a per-level budget in the extension task.

None of these changes the lane step functions. They extend the program table and
the interpreter (`Mm_Run`).

## Next

- **Default batch size:** make batching the default, with B around 8.
- **One path for `--batch 1`:** retire the per-read path? The engine matches it
  at one lane, so `--batch 1` could run through it too, leaving one
  implementation. The per-read functions are still used by the tests as the
  reference.
- **k-mismatch programs:** when the k-mismatch scheme is designed, add it as a
  program, together with the segment table and the chained-extension form above.

## Choosing B (`--batch auto`, the default)

**What sets the best lane count.** It is the number of cache misses one core can
have in flight, not memory size.
- **Why lanes help:** each search step waits on a dependent cache miss, and B
  lanes overlap B of them.
- **Where it stops:** the gain ends once the core's line-fill buffers are
  busy. Their number differs between CPUs, and between this machine's P-cores
  and E-cores.
- **Lane state is small** (stack, read, result lists): even 64 lanes cost
  little memory.
- **The workload decides most.** Long single-row walks gain most (`-n 0`,
  `--no-ref`). Branching extensions already have several misses in flight per
  read, so more lanes only add cache pressure (`-n 4`).

**Sweep:** search time relative to B = 1, 20k reads (mostly exact) on the
500 Mb genome, `-m 0`, load 4–5, pinned to one core. P is CPU 2 (P-core), E is
CPU 14 (E-core), LP is CPU 20 (low-power E-core). Ranges cover B = 2 to 64.

| workload | P | E | LP |
|---|---|---|---|
| `-n 0` | 0.59–0.78 | 0.84–0.97 | 0.75–1.12 |
| `--no-ref -n 1` | 0.60–0.71 | 0.72–0.81 | 0.75–0.79 |
| `-n 1` | 0.79–0.96 | 0.84–1.00 | 0.75–1.06 |
| `-n 2`, `-n 3` | 0.87–1.04 | 0.89–1.11 | 0.74–1.19 |
| `-n 4`, 100 bp | 1.17–1.32 | 1.00–1.15 | 1.05–1.45 |
| SAM `-n 2` | 0.95–1.05 | 0.96–1.18 | 0.83–1.06 |

No fixed B suits every case. B = 16 costs up to 30% at `-n 4` on the P-core,
and B = 1 gives up 1.4–1.7× on `-n 0` and `--no-ref`.

**`--batch auto`** (`Batch_Chooser`):
- **The trial:** the first blocks (4096 reads each) try the candidates
  B = 1, 2, 4, 8 and 16, one block each. There are 3 rounds, in alternating
  order, which cancels drift such as caches warming up.
- **Timing:** a block's first 256 reads only warm up; the rest is timed.
- **Choice:** the lowest median time per read is kept for the rest of the run,
  or, within 3% of it, the one with the fewest lanes.
- **Cost:** the trial takes 15 blocks, about 61k reads. On a short input it ends
  before choosing, and `--timing` then reports `batch=auto:-`. An explicit
  `--batch N` skips it.

**Two first attempts that failed:**
- **Alternating every few hundred reads.** Chunks of 256–512 reads, round-robin
  over B values, never reached steady state. In the trial, B = 1 measured 39–42
  µs per read against 37–39 for B = 4–8. Run whole, B = 1 is 12–20% faster than
  either at `-n 4`. Every switch disturbs caches and branch prediction, and the
  single lane loses most.
- **A trial on the first reads only.** It assumes those reads are typical. An
  input sorted by read type, like this project's mixed sets (exact reads first),
  misleads it. Real read files are not sorted that way.

**Evaluation:** 400k reads (100k each with 0, 1 in LH, 1 in RH and 2 mismatches,
shuffled), search seconds, one run each (±5%).

| workload | P: auto / best fixed | P: auto / B = 1 | P: B chosen | E: auto / best fixed | E: auto / B = 1 | E: B chosen |
|---|---|---|---|---|---|---|
| `-n 0`, 100 bp | 0.94 | 0.80 | 2 | 1.07 | 0.88 | 2 |
| `-n 0`, 32 bp | 1.07 | 0.83 | 2 | 1.02 | 0.95 | 4 |
| `-n 1`, 100 bp | 1.01 | 0.90 | 16 | 0.97 | 0.89 | 8 |
| `-n 2`, 100 bp | 1.05 | 1.01 | 16 | 1.05 | 0.93 | 4 |
| `-n 2`, 32 bp | 1.06 | 1.02 | 8 | 1.03 | 0.97 | 2 |
| `-n 3`, 100 bp | 1.01 | 0.98 | 2 | 1.06 | 1.01 | 8 |
| `-n 3`, 32 bp | 0.99 | 0.94 | 2 | 1.09 | 0.97 | 2 |
| `-n 4`, 100 bp | 1.03 | 1.03 | **1** | 1.05 | 1.05 | 2 |
| `--no-ref -n 1` | 1.04 | 0.69 | 4 | 0.77 | 0.65 | 4 |
| SAM `-n 2` | 1.09 | 0.99 | 16 | 0.98 | 0.95 | 2 |

- **Close to the best fixed B:** `auto` stays within 0.94–1.09 of it on the
  P-core (single runs, ±5%).
- **It avoids the bad cases:** at `-n 4` on the P-core it chose B = 1. Every
  B ≥ 2 takes 17.0–17.1 s there, against 14.3 s for B = 1, so B = 16 would cost
  19%.
- **It keeps the large gains:** 0.65–0.69 of the B = 1 time with `--no-ref`,
  and 0.80–0.95 at `-n 0`.
- **Output is identical** for every B, and all runs now use the engine. The `-n 0`
  per-read loop in `main` is gone.

**With threads (`-t`, `threads.md`):**
- **Each thread has its own chooser,** timed on its own 4,096-read chunks while the
  other threads run.
- **The best B is lower under load.** Two hyperthreads on one core share its fill
  buffers, so B is about 2–4 instead of 4–16, and lanes gain little on a hyperthread
  pair.
- **One thread per core keeps the gain** (0.51–0.57 of the B = 1 time with `--no-ref`).
- **`auto` stays close:** within 1.00–1.12× of the best fixed B at 6, 12 and 22
  threads.
