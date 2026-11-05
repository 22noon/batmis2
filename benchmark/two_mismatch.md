# Two-mismatch search (`-n 2`)

> **Since the program refactor** (`five_mismatch.md`, "Program builder"): the hand-written
> per-read searches named below (`One_Mismatch_Search` … `Five_Mismatch_Search`,
> `Extend_Windows`, `Switch_Index`, the `Program_*` arrays, the `SEG_*` / `WIN_*` enums) are
> gone. `Mm_Builder` builds the same programs step for step, and `--batch 1` runs them
> on one lane of the engine.

## Summary

`batman_zero_min -n 2` reports hits with up to two mismatches. It exists both
per read (`Two_Mismatch_Search`) and batched (`Program_Two_Mismatch`,
`--batch B`). The two share the 1|1 extension and switch code.

The read's halves LH | RH are the pigeonholes. p and q denote a 1|1 hit's
mismatches in LH and RH, and T1 = [0, a), T2 = [a, b), T3 = [b, L) are the read's
thirds, with a = L/3 and b = 2L/3, so a < LH < b.

| case | how it is found | list |
|---|---|---|
| 0\|0, 0\|1, 1\|0 | `One_Mismatch_Search`, unchanged | exact, `F`·1, `B`·1 |
| 0\|2 | the 0\|1 pass's rejects (deferred ranges) extended over RH with Count 2 | `F`·2 |
| 2\|0 | the 1\|0 pass's rejects extended over LH with Count 2 | `B`·2 |
| 1\|1 (a), p ∈ [a, LH) | T1 exact, extended right over [a, L) with windows {[a, LH): 1, [LH, L): 1} | `F`·2 |
| 1\|1 (b), p ∈ T1, q ∈ [LH, b) | T3 exact, extended left over [0, b) with windows {[LH, b): 1, [0, a): 1} | `B`·2 |
| 1\|1 (c), p ∈ T1, q ∈ T3 | T2 exact, extended right over T3 with {T3: 1}; each such seed switched to fwfmi (unless located), then extended left over T1 with {T1: 1} | `B`·2 |

**The 1|1 partition** (your suggestion: reuse the 0|1 / 1|0-style extensions, and
keep the cases apart by which positions each may change):
- **Coverage:** every 1|1 hit has p in T1 or in [a, LH), and q in [LH, b) or in
  T3. Each of the three cases takes a disjoint part of those combinations, and
  the third it seeds on is exact in every hit of that part.
- **No duplicates:** a hit falls in exactly one case, because the cases'
  conditions on p and q are disjoint.
- **No overlap with 0|2 / 2|0:** every case requires one mismatch in each half.

**Mismatch windows** replace the single mismatch count in these extensions
(`Mm_Windows`):
- A node may take a mismatch only at a position whose window has room
  (`Window_Can_Change`).
- It is dropped as soon as a window can no longer reach its count with the
  positions left (`Window_Keep`).
- Positions outside every window must match.
- Located nodes check all windows on the `.pac` mismatches (`Windows_Hold`).
- Nothing is deferred: 1|1 has no later round at k = 2.

**Bounds from the exact half searches** (section 4): where the exact LH and RH
searches failed narrows every window, and skips (a) or (b) outright when their
windows cannot hold.

**The index switch** (case (c)) walks Start alone (`Mm_Start_Switch`):
- The seed's string, read[a, L) with q's mismatch applied, occurs (the seed
  proves it). So every range on the way is non-empty, and the Start-only step
  S′ = C[c] + Occ(c, S) + 1 is exact.
- A string has the same number of occurrences in both indexes, so
  End = Start + (seed width) − 1.
- There are no End ranks, `Count_Range` calls or character checks. It replaces
  the earlier full re-search.

**Against the earlier quarter seeding:**
- **32 bp:** 1.4–1.8× faster. Thirds are 10–11 bases instead of 8, so they start
  from about 16× fewer rows.
- **75 and 100 bp:** on their own, the thirds were up to 17% slower at `-m 0`,
  because of a third exact search.
- **With the bounds** (section 4), long reads match or beat quarters: 2–9%
  faster, except p75 `-m 0 --batch 16`, which is within 1% (Timing below).

## 1. Against the original BatMis (`alternate/batman.cpp`, two-mismatch block)

**0|2 and 2|0: the same algorithm.**
- **The original:** keeps the threshold rejects of the 0|1 forward scan
  (`Mismatches_Forward`) and of the 1|0 backward scan (`Mismatches_Backward`),
  and calls `Search_Forwards(…, 2, …)` / `Search_Backwards(…, 2, …)` on them.
  That is "resume the rejects with Count 2".
- **Two differences in order, not in the hits:**
  - **Reject order:** it resumes the rejects last-first.
  - **Hits at the last base:** it prints the 2-mismatch hits whose second
    mismatch is on the read's last base first
    (`Two_Mismatches_At_End[_Forward]`). Here that path is a threshold reject
    at the last level, and resuming it finds the same hit, in its stage.

**1|1: a different partition, the same hits.**
- **The original:** seeds RH's quarters (`Search_Backwards_X10`,
  `Search_X01` + `Reverse`), then extends over LH with Count 2.
- **Here:** seeds the thirds, with windows.
- **Both are complete for 1|1.** Neither reports a hit twice: the original by
  the quarter cases being disjoint, this one by the windows.
- **What the thirds change:** the seeds are longer than quarters, and two of
  the three cases need no index switch.

## 2. Bookkeeping: is "rejects in, hits + rejects out" enough?

**The original:** a separate array per category, some read and appended in place
(the 0|2 scan reads [0, P) of `Mismatches_Forward` and appends its 0|3
candidates after P).

**The stage model.** Each stage reads the previous round's rejects of its own
chain, and writes hits and new rejects.
- **Two chains:** forward (0|0 → 0|1 → 0|2) and backward (1|0 → 2|0).
- **They interleave by stage,** so the backward rejects must still be alive
  while 0|2 runs. One set of three structures is therefore not enough.
- **Enough, for any k:** an input and an output reject list per chain, one hit
  list per read tagged with (index F/B, mismatches), and a seed list for seeded
  cases such as 1|1 (c).

**Now (rolling lists, `Result_List`):** seven lists per read, for `-n 1` to 3
alike.

| list | holds |
|---|---|
| `RL_EXACT` | 0\|0 hits; its rejects are the forward chain's first round |
| `RL_HITS` | every mismatch hit, in stage order, each tagged with its stage (`SARange::Stage`); also the running stage's rejects |
| `RL_FWD`, `RL_BWD`, `RL_CROSS` | each chain's latest round of rejects (0\|k, k\|0, 1\|k) |
| `RL_SEED_1`, `RL_SEED_2` | the running stage's one-mismatch seeds |

- **Ending a stage** (`End_Stage`, `OP_END_STAGE`):
  - tags the stage's new hits;
  - swaps its rejects into its chain's list, which drops the round just read;
  - clears the seeds.
- **Each chain holds one round at a time,** plus the shared scratch list, so the
  list count no longer grows with k.
- **Output:** `Format_Hits` walks `RL_HITS` once and takes each hit's index and
  mismatch count from the search's stage table (`One_Stages`, `Two_Stages`,
  `Three_Stages`).
- **Order and `MaxHits` are unchanged:**
  - the stages already ran in output order;
  - an extension into `RL_HITS` gets the budget `MaxHits` minus the exact
    hits, and `RL_HITS.hitCount` already counts every earlier mismatch hit.
- **The per-read searches stop at the same points as their programs.** A
  program gained one `OP_STOP_IF_MAXHITS` after the shared passes, so both
  leave identical lists, which `make test-located` checks field by field,
  stages included.
- **Checks:**
  - byte-identical output to the per-stage layout (140 runs on 500 Mb: `-n 1/2/3`,
    `-m` 0/1/3, batch 1/16, both reference modes);
  - brute force at distances 2 and 3 (64 runs).
- **Time:** unchanged within noise on 500 Mb (mean of 3 runs, `-n 1/2/3`,
  batch 1 and 16). Changes ran from −3% to +2%, except one outlier: −8% for
  p100 `-n 2 -m 1 --batch 16`.

## 3. Shared routines

**Sequential.**
- **`Search_Segment(segment, offset, direction)`:** exact search of any
  segment. It walks from the end base when the segment is shorter than the
  lookup key (thirds can be).
- **`Extend_Windows` and `Switch_Index`:** run the batched engine's own window
  and switch steps (`Mm_Window_Step`, `Mm_Switch_Step`) on a private lane, so
  the two paths share one implementation of the new code.
- **`Read_Segment_Of`:** the read's segments (halves, thirds, [a, L), [0, b)).

**Batched.**
- **`MM_ONE_MISMATCH_PASSES`** (the `-n 1` passes) is shared by both programs.
- **1|1** is `OP_SEARCH` on a third, `OP_EXTEND` with a window set, and
  `OP_CROSS` (seed → switch → windowed extension) for case (c).
- **Extension tasks carry a task kind:** search, extend, windowed extend or
  switch.

**`-n 0` / `-n 1` are unaffected.** Their output is byte-identical to before.
Their timing was within noise when the shared pieces were introduced (see the
first `-n 2` table in `batched_reads.md`).

## 4. Pruning 1|1 with the exact half searches

**The 0|1 / 1|0 rejects can't help 1|1.** Every reject has one half exact, and a
1|1 hit has a mismatch in both halves, so no reject leads to one.

**Where the exact searches failed can.**
- **d:** if the 0|0 pass's exact LH search (revfmi) fails at depth d, read[0, d)
  does not occur, so every hit has a mismatch in [0, d). For 1|1 that is p, so
  p < d.
- **e:** if the 1|0 pass's exact RH search (fwfmi) fails at position e, read[e, L)
  does not occur, so q ≥ e.
- A search that succeeds gives no bound (d = L, e = 0).

**What the bounds change** (`Window_Set`): every LH window ends by d, and every
RH window starts at e or later.

| case | windows after the bounds | skipped when |
|---|---|---|
| (a) | [a, min(LH, d)): 1, [max(LH, e), L): 1 | d ≤ a |
| (b) | [0, min(a, d)): 1, [max(LH, e), b): 1 | e ≥ b |
| (c) | T3: [max(b, e), L): 1, then T1: [0, min(a, d)): 1 | never |

- **Skipping a case saves its exact third search.** On long reads that search is
  most of the case's cost.
- **Narrowing the windows** drops branches that the index would reject only
  later. For example, a mismatch at p in [d, a) during a leftward extension
  leaves read[0, d) to match exactly. That can't succeed, but the walk would
  find out only at position 0.
- **No hits are lost or duplicated:** the bounds hold for every 1|1 hit, and the
  windows only lose positions where no hit can have its mismatch.

**Where d and e come from.** The search stops in one of four ways, and each
leaves the failure depth in `Level`:
- **In the lookup table:** the table range is empty, so the first k bases don't
  occur. d = k, a valid but looser bound.
- **On a wide range or a single row:** the base that failed.
- **At the `$` row:** the base that could not be extended.
- **After `Locate`:** the first mismatch of the `.pac` compare
  (`First_Mismatch`, `Located_Matches`). The located row was the only
  occurrence of the prefix, so the prefix up to that base does not occur. If
  the read would run off the genome, the bound is the whole segment.

`One_Mismatch_Search` returns d and e to `Two_Mismatch_Search`. In the batched
engine, the `OP_SEARCH` steps of LH and RH set them on the lane, and `OP_SEARCH`
with a window set skips a case whose windows cannot hold.

**Effect on 500 Mb** (`--batch 16`, `-m 0`, 20k reads):
- **75 and 100 bp:** (a) and (b) are each skipped for 55–57% of reads, and the
  windowed extensions also shrink a little.
- **1|1 memory accesses per read:** p75 34.3 → 27.3, p100 37.9 → 29.3. That is
  almost back to the quarter scheme (25.7, 28.1).
- **32 bp:** almost unchanged (113.6 → 112.8). d is rarely ≤ a = 10 there,
  because a 10-base prefix of a read always occurs in 500 Mb.

**Time** (`-n 2`, input + search seconds, mean of 2 runs, `--lookup 11`,
reference on, 500 Mb genome, both builds in the same session, load average
about 2):

| reads | `-m` | thirds, batch 1 | + bounds, batch 1 | thirds, `--batch 16` | + bounds, `--batch 16` |
|---|---:|---:|---:|---:|---:|
| r32 | 0 | 51.8 | 49.6 | 36.4 | 35.8 |
| r32 | 1 | 8.89 | 8.85 | 5.75 | 5.75 |
| p32_mm2 | 0 | 23.1 | 22.7 | 17.0 | 17.1 |
| p32_mm2 | 1 | 10.1 | 10.1 | 7.35 | 7.26 |
| p75_mm2 | 0 | 3.38 | 3.29 | 2.46 | 2.28 (−7%) |
| p75_mm2 | 1 | 2.57 | 2.54 | 1.70 | 1.61 (−5%) |
| p100_mm2 | 0 | 3.71 | 3.61 | 2.69 | 2.46 (−8%) |
| p100_mm2 | 1 | 2.81 | 2.76 | 1.88 | 1.73 (−8%) |

- **Batched:** long reads gain 5–8%.
- **Per read:** the gain is 1–3%, because the skipped searches are a smaller
  share of the time when each step waits for its own cache miss.
- **32 bp:** within noise.
- **`-n 0` / `-n 1`:** within noise (r32 and p100, batch 1 and 16). They now
  use `First_Mismatch` in place of the mismatch count for a located search.

## Correctness

- **Brute force:** reported positions, split by mismatch count, are exactly the
  genome positions within Hamming distance 2, on every read. No position is
  repeated.
  - **Genomes:** a 1 Mb random genome with planted near-repeats, and the 49 kb
    repetitive toy.
  - **Reads:** 32, 36, 75 and 100 bp (the thirds' boundaries fall differently
    against LH), with mismatches in every placement: 0|1, 1|0, 0|2, 2|0 and
    all four (p, q) region pairs of 1|1.
  - **Settings:** both reference modes, `--lookup` 6 and 11.
  - **Size:** 48 runs of 540–1350 reads each.
- **500 Mb samples** (20k reads each of r32, p32_mm2, p75_mm2, p100_mm2,
  p32_mm1L): every reported position is at its reported distance. The counts per
  mismatch class are identical to the quarter scheme's.
- **Batched vs sequential:** byte-identical in every oracle run at `--batch 16`,
  and in 60 runs on 500 Mb (batch 1, 16 and 32, `-m` 0 and 1, both reference
  modes).
- **`-n 0` / `-n 1`:** byte-identical to the build before batching (852 runs).
- **Bounds (section 4):**
  - **Brute force:** matches every read in 64 runs: the 32 runs above, plus the
    same matrix at `--lookup` 3 and 14, where the search can fail inside the
    lookup table.
  - **Batched vs sequential:** byte-identical in all of them.
  - **500 Mb samples:** `-n` 0, 1 and 2 are byte-identical to the build
    without bounds (120 runs: 5 read sets, both reference modes, `-m` 0 / 1,
    batch 1 / 16).
  - **Fault injection:** each fault fails `make test-located`:
    - d one too small: many failures;
    - e one too large: 20 failures;
    - the located depth one too small: 5 failures.
- **`make test-located`:**
  - batched against per-read `-n 2`, all lists field by field;
  - per-read `-n 2` against a brute-force scan.

  A fault in each window set (shifting one window by one base) is caught, with
  88–126 failures each. `make test-structured` passes.
- **Two bugs were found by the oracle during the work:**
  - **Argument order:** `Search_Segment(Read_Segment_Of(All, Seg, Offset), Offset, …)`
    read `Offset` before the same call set it. C++ leaves argument order
    unspecified, so located compares used the previous segment's offset. It
    showed only with the reference on.
  - **A false shortcut in the switch:** I had added "locate at the first
    sampled row when the seed is a single row". During a Start-only walk the
    range is wider than the seed's until the last base, so Start is not yet the
    seed's row. The shortcut was removed.

## Timing

`-n 2`, input + search seconds, mean of 2 runs, `--lookup 11`, reference on,
500 Mb genome. The machine had a load average of about 2, so expect ±5% noise.

| reads | `-m` | quarters, batch 1 | thirds, batch 1 | quarters, `--batch 16` | thirds, `--batch 16` | thirds vs quarters |
|---|---:|---:|---:|---:|---:|---:|
| r32 (2M) | 0 | 74.1 | 52.5 | 54.4 | 39.4 | 1.38–1.41× |
| r32 | 1 | 12.9 | 9.7 | 9.27 | 6.59 | 1.33–1.41× |
| p32_mm2 (1M) | 0 | 35.5 | 24.3 | 25.3 | 18.1 | 1.39–1.46× |
| p32_mm2 | 1 | 18.7 | 10.7 | 13.1 | 7.40 | 1.76× |
| p75_mm2 (500k) | 0 | 3.39 | 3.40 | 2.33 | 2.49 | 0.93–1.00× |
| p75_mm2 | 1 | 2.75 | 2.69 | 1.81 | 1.76 | 1.02× |
| p100_mm2 (500k) | 0 | 3.71 | 3.79 | 2.57 | 3.11 | 0.83–0.98× |
| p100_mm2 | 1 | 3.03 | 2.88 | 1.94 | 1.92 | 1.01–1.05× |

Memory accesses per read (`-DBATMAN_BATCH_STATS`, 20k reads, `--batch 16`):

| reads | stage | quarters | thirds | thirds + bounds |
|---|---|---:|---:|---:|
| p32_mm2 | up to 2\|0 | 17.5 | 17.5 | 17.5 |
| | 1\|1 | 194.1 | 113.6 | 112.8 |
| | total | 211.6 | 131.1 | 130.2 |
| p75_mm2 | up to 2\|0 | 23.2 | 23.2 | 23.2 |
| | 1\|1 | 25.7 | 34.3 | 27.3 |
| | total | 49.0 | 57.5 | 50.5 |
| p100_mm2 | up to 2\|0 | 25.5 | 25.5 | 25.5 |
| | 1\|1 | 28.1 | 37.9 | 29.3 |
| | total | 53.6 | 63.3 | 54.7 |

The thirds column was measured again for this table. The first measurement came
from a stats build made before the located shortcut was removed from the switch,
and it differed slightly.

**Why 32 bp gains.** The seed is a 10–11-base third rather than an 8-base
quarter, so the windowed extension starts from far fewer rows.

**Why long reads lost a little before the bounds.** Long-read 1|1 is dominated
by its three exact third searches, about 9 steps each (wide, then single row,
then locate). The quarter scheme had two such searches. The bounds skip (a) and
(b) for about 55% of reads.

**Quarters against thirds + bounds, measured side by side** (mean of 3 runs,
load average under 0.6):

| reads | `-m` | quarters, batch 1 | thirds + bounds, batch 1 | quarters, `--batch 16` | thirds + bounds, `--batch 16` |
|---|---:|---:|---:|---:|---:|
| p75_mm2 | 0 | 3.31 | 3.24 (−2%) | 2.23 | 2.26 (+1%) |
| p75_mm2 | 1 | 2.71 | 2.52 (−7%) | 1.73 | 1.59 (−8%) |
| p100_mm2 | 0 | 3.65 | 3.59 (−2%) | 2.49 | 2.45 (−2%) |
| p100_mm2 | 1 | 2.96 | 2.76 (−7%) | 1.90 | 1.72 (−9%) |

Thirds + bounds still make slightly more memory accesses per read than quarters
on long reads (p75 50.5 against 49.0, p100 54.7 against 53.6), so at `-m 0` the
two are close.

## 5. T1 and T3 kept from the half scans

**The idea:** the exact searches of (a) and (b) repeat work the half scans already
did.
- **T1:** the 0|0 pass's LH scan (revfmi) passes depth a, so T1's range is a
  state of that scan.
- **T3:** the 1|0 pass's RH scan (fwfmi) passes depth L − b, which gives T3's
  range.
- Keeping both leaves T2 as the only new exact search.

**How:** `Search_Right_Stored` / `Search_Left_Stored`, and `OP_SEARCH_KEEP` in the
batched engine.
- The scan searches the third first and keeps its range (`Half_Scan`), then
  goes on from that state over the rest of the half.
- A located range compares the rest with the genome (`Located_Matches`).
- Every state on the way is the one the whole search would reach, so the 0|0 and
  1|0 results are unchanged.
- The kept range is what a fresh search of the third would give: the same SA
  range, and located at the same row.
- (a) and (b) use it (`OP_SEARCH` with `Stored`) in place of their search.

**When there is no kept range:**
- **The lookup table already covers the third** (k > a, e.g. `--lookup 11` on
  32 bp): (a) / (b) search as before.
- **The scan failed inside the third:** the bound (section 4) already skips that
  case.

**Memory accesses per read** (`--batch 16`, `-m 0`, 20k reads, 500 Mb):

| reads | `--lookup` | bounds only | + kept T1/T3 |
|---|---:|---:|---:|
| p100_mm2 | 11 | 54.7 | 44.2 (−19%) |
| p75_mm2 | 11 | 50.5 | 41.4 (−18%) |
| p100_mm2 | 6 | 74.2 | 59.2 (−20%) |
| p75_mm2 | 6 | 70.0 | 56.5 (−19%) |
| p32_mm2 | 6 | 154.7 | 145.7 (−6%) |
| p32_mm2 | 11 | 130.2 | 130.2 (table covers T1) |

**The scans themselves** take about 0.7 more steps per read: the located compare
is split at the third's end. Against quarters (49.0 and 53.6), thirds now need
16–18% fewer memory accesses on long reads.

**Checks:**
- **Brute force:** matches every read in 64 runs (`--lookup` 3, 6, 11, 14).
- **Byte-identical to the build without kept ranges:** `-n 0/1/2` on 500 Mb
  (120 runs, `--lookup` 6 and 11).
- **Fault injection:** each fault fails `make test-located` (20 failures each):
  - keeping T1 one base short (per read);
  - keeping RH instead of T3 (batched).

**Time** (`-n 2`, seconds, mean of 3 runs, `--lookup 11`, load average 1.1–1.8):

| reads | `-m` | bounds, batch 1 | + kept, batch 1 | bounds, `--batch 16` | + kept, `--batch 16` |
|---|---:|---:|---:|---:|---:|
| p75_mm2 | 0 | 3.32 | 3.28 | 2.33 | 2.13 (−8%) |
| p75_mm2 | 1 | 2.53 | 2.54 | 1.58 | 1.57 |
| p100_mm2 | 0 | 3.68 | 3.57 | 2.50 | 2.30 (−8%) |
| p100_mm2 | 1 | 2.79 | 2.80 | 1.72 | 1.74 |
| p32_mm2 | 0 | 22.9 | 22.0 | 17.2 | 16.1 (−7%) |
| p32_mm2 | 1 | 10.0 | 9.69 | 7.16 | 7.10 |

- **Long reads at `-m 0`:** 8% faster batched. At `-m 1` there is no change,
  since most reads stop before 1|1.
- **p32 at `-m 0`:** also 4–7% faster, although its memory accesses are
  unchanged (the table covers T1). The build also has the `-n 3` code, so this
  is probably code layout, not this change.

## 6. Why the thirds stop at k = 2

Three pieces with two mismatches always leave one piece exact. With three, they
don't. The thirds cut at LH give four pieces: T1, T2L = [a, LH), T2R = [LH, b)
and T3. Raising one window of the 1|1 extensions to 2 covers three of the four
1|2 cases:

| 1\|2 case (p; q1, q2) | from | exact |
|---|---|---|
| p ∈ T2L; q's anywhere in RH | (a), RH window → 2 | T1 |
| p ∈ T1; both q's in T2R | (b), T2R window → 2 | T3 |
| p ∈ T1; both q's in T3 | (c), T3 window → 2 | T2 |
| p ∈ T1; q1 ∈ T2R, q2 ∈ T3 | **none** | only T2L, a sixth of the read |

- **2|1 has the mirror gap:** p1 ∈ T1, p2 ∈ T2L, q ∈ T3, where only T2R is exact.
- **Filling the gaps** needs a one-mismatch seed, which is what the quarter
  scheme is built on.
- **So the rule is one more piece than the number of mismatches:** `-n 2` uses
  thirds and `-n 3` uses quarters (`three_mismatch.md`).

## Not done yet (from the earlier discussion)

- **Streaming seeds:** extend each (c) seed over T1 as soon as it is found,
  instead of collecting `TM_SEED_C`. That saves the list (memory and copies),
  and work after `MaxHits` is reached. It saves no steps at `-m 0`.
- **Rolling reject lists:** see Bookkeeping. `-n 3` went ahead without them
  (`three_mismatch.md`).
