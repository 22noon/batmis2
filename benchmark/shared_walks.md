# Shared seed walks

## Summary

- **The high-K search time is the exact seed searches.** Profiled at 150 bp, `-n 7`, all
  hits: exact seed searches take 64% of the search's cycles, mismatch extensions 16%,
  crossings into the other half 4%, the chains 1.5% and the prefetch pass 15%.
- **Many of those searches repeat work.** At `-n 7` the program searches the same segment
  up to three times per strand (in the 1|k, 2|k and 3|k families), and other seeds are
  prefixes, in search order, of a longer walk from the same origin: [0,18), [0,25) and
  [0,37) all lie on the exact half's walk over [0,75).
- **Now an earlier walk serves them.** The walk saves the state each later search would end
  in, and that search takes it instead of searching. Output is byte-identical.
- **Search time at K = 5-7 falls by 26-32% for all hits and 14-23% for first hit**, in both
  mappers, at 100 and 150 bp. K <= 3 is unchanged (nothing to share).

## Profile (before)

E. coli, 100,000 simulated 150 bp reads (equal shares with 0-7 substitutions, 10% random),
SAM, both strands, `-m 0`, one P-core. Cycles counted per program step (rdtsc around each
lane step, in a scratch build with `-DBATMAN_BATCH_STATS`):

| part of the search | share of cycles |
|---|---:|
| the seed families 2\|k, k\|2, 3\|k, k\|3 | 65% |
| the exact halves and 1\|k, k\|1 | 18% |
| the chains (resumed rejects), all rounds | 1.5% |
| the prefetch pass over the lanes | 15.5% |

By operation: exact searches 64%, extensions 16%, crossings 4%. Each exact seed search costs
about 6 index steps per strand plus a locate and one genome compare (on a 4.6 Mb genome a
row is single after about 11 bases); there are 28 per strand at `-n 7`. The seeds are not
the problem: the smallest, an eighth of the read (18 bases at 150 bp), is already specific.

Searches that were repeats or prefixes of an earlier walk, per strand at `-n 7`:

| served search | from the walk |
|---|---|
| [37,75) leftward, twice | [37,75) leftward (1\|k) |
| [50,75), [56,75) leftward | [37,75) leftward |
| [75,112) rightward, twice | [75,112) rightward (k\|1) |
| [75,100), [75,93) rightward | [75,112) rightward |
| [0,37) (already a kept range), [0,25), [0,18) rightward | the exact half [0,75) |
| [112,150) (kept range), [125,150), [131,150) leftward | the exact half [75,150) |

They cost 28% of the cycles (22% at first hit).

## How it works (`src/batmis1.cpp`, `src/batmis2.cpp`, identical)

- **`Mm_Builder::Share_Walks`**, after each program is built: an `OP_SEARCH` of a segment
  other than the read or a half is a *consumer* if an earlier walk (`OP_SEARCH` that is not
  a consumer, or `OP_SEARCH_KEEP`) covers it: same direction, same origin, and the
  consumer's string no longer, at every read length from 16 to 300 (`Covers`). A walk that
  always runs (no window, no kept range) is preferred. Each served segment gets a slot
  (`Walk_Use`); each walk the list of segments it serves (`Walk_Make`, `Walk_Table`).
  Searches of the halves stay as they are: their failure depths bound the windows.
- **The walk saves** (`Walk_Arm`, `Walk_Capture`, `Walk_Settle`, `Walk_Fail`): at start,
  its points that lie on it at this read length and that a step will reach (beyond the
  lookup table's bases) are pending. At each depth it then saves the state a search of that
  length would end in: the range (wide), the row with `End` set as a finished single-row
  search sets it, or, once located, the located tag. A failed step fails the pending points
  (they all contain the failing base). A located walk compares the rest with the genome in
  one go; its first mismatch settles which points match.
- **The consumer** takes a settled slot (found: the saved state; failed: no range) after its
  window check and kept range, as before. Unknown (the walk was skipped, served itself, or
  the geometry does not hold at this length): it searches, as before.
- **Why the output is unchanged:** a prefix search follows the same steps as the walk up to
  its depth: the same lookup start, the same range steps, the same sampled rows and the
  same locate (`Locate`'s position depends only on the origin and depth). Nothing after a
  search reads its `Level`, the one field that differs.
- **Per lane:** up to 16 slots a strand (`WALK_SLOTS`), up to 8 points a walk
  (`WALK_POINTS`).

Exact searches per strand, and how many a walk serves:

| `-n` | 2 | 3 | 4 | 5 | 6 | 7 |
|---|---:|---:|---:|---:|---:|---:|
| exact searches | 5 | 6 | 11 | 16 | 22 | 28 |
| served by a walk | 2 | 2 | 6 | 10 | 14 | 18 |
| walks serving them | 2 | 2 | 3 | 4 | 4 | 4 |
| slots | 2 | 2 | 5 | 8 | 10 | 12 |

At `-n 3` the two served searches make no measurable difference to the time.

## Timing

Same reads and machine, `--sam -o /dev/null`, `--timing` search seconds, median of 3.
Before → after:

| reads | mapper | `-m` | K = 3 | K = 5 | K = 6 | K = 7 |
|---|---|---|---|---|---|---|
| 150 bp | batmis1 | all | 0.39 → 0.39 | 0.98 → 0.71 | 1.35 → 1.11 | 1.97 → 1.37 |
| 150 bp | batmis1 | first | 0.32 → 0.34 | 0.57 → 0.47 | 0.67 → 0.51 | 0.68 → 0.54 |
| 150 bp | batmis2 | all | 0.49 → 0.50 | 1.21 → 0.88 | 1.67 → 1.19 | 2.16 → 1.47 |
| 150 bp | batmis2 | first | 0.36 → 0.36 | 0.62 → 0.50 | 0.71 → 0.58 | 0.77 → 0.60 |
| 100 bp | batmis1 | all | 0.42 → 0.42 | 1.04 → 0.76 | 1.52 → 1.10 | 2.00 → 1.45 |
| 100 bp | batmis1 | first | 0.30 → 0.30 | 0.53 → 0.43 | 0.63 → 0.51 | 0.69 → 0.55 |
| 100 bp | batmis2 | all | 0.48 → 0.48 | 1.16 → 0.85 | 1.70 → 1.21 | 2.20 → 1.58 |
| 100 bp | batmis2 | first | 0.34 → 0.34 | 0.60 → 0.48 | 0.71 → 0.56 | 0.76 → 0.66 |

The K = 3 differences (up to 5%) are run-to-run noise on this machine.

## Correctness

- **Byte-identical SAM** (without `@PG`) against the build before, for both mappers at
  `-n 2..7`, `-m 0` and `-m 1`, on 20,000 reads of 150 bp.
- **Other modes, byte-identical too:** 672 comparisons in all, both mappers, `-n 4..7`,
  `-m 0, 1, 5`, reads of 32, 40, 100 and 150 bp, each with `--sam`, `--sam --no-ref`,
  `--sam --forward-only`, text output, `--sam --batch 1`, `--sam -t 4` and
  `--sam --lookup 11`. No difference.
- **`make test` passes** with both SA sampling settings (the 1 GB budget, and every 8th row).

## What is left

- **The remaining cost is the number of cases.** 10 of `-n 7`'s 28 exact searches per
  strand still search, each with its own extension. Columba's search schemes cover the same
  mismatch placements with fewer searches; that is the next gap to look at.
- **The prefetch pass** takes 15% of the cycles, more than address arithmetic should.
