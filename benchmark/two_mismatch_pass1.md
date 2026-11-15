# What pass 1 tells pass 2 in `-n 2` (batmis1 / batmis2)

Two questions, investigated on 2026-11-05:
1. Why did the hand-written bidirectional `bi_mm -n 2` (`attic/bifmindex`) gain so much
   from settling pass-2 searches with pass 1's results, and does batmis2 need it?
2. Can a prefix or suffix that does not occur stretch pass 2's seeds?

**Answers:**
1. batmis1 and batmis2 already get settling's effect, through located kept ranges and
   resumed rejects.
2. Stretching seeds is valid but would save about 3% at 100 bp and nothing at 32 bp.
   Neither is implemented.

## Pass 1's searches and what they leave

Pass 1 searches the read's left half LH exactly, from the read's start, and its right half
RH from its end. The search is exact, so at depth j it holds one SA range: every
occurrence of P[0, j), one row each (two occurrences are two rows, not two paths). It ends
early in one of two ways:
- **Empty range:** P[0, j) does not occur. This is the failure depth (`LH_Fail`; for the
  suffix, `RH_Fail`).
- **A single row whose SA value is sampled:** P[0, j) occurs exactly once, at read start S
  (located). Adding bases can only remove occurrences, so P[0, j) stays unique for any
  longer j.

A range of two or more rows tells nothing yet; the search goes on. If the half ends with
several rows (a repeat), pass 1 learns nothing about that side.

## 1. Settling

**What it is** (in `bi_mm`, `bi_two_mismatch.hpp`): a pass-2 search seeded at the read's
start (A: LH exact; C1: the first third exact) whose seed contains pass 1's ended prefix
P[0, j) is decided without being walked.
- **Prefix absent:** the search has no hits.
- **Prefix unique at S:** every hit is at S, and one genome comparison decides it.

The same holds for the suffix with the RH-seeded searches.
- **100 bp:** 95% of these searches were settled. The prefix is unique after about 15
  bases, well inside a 33-base third.
- **32 bp:** 23%. A 10-base third is too short to be unique.

**Why batmis1 and batmis2 do not need it.** Their `-n 2` program reaches the same decision
on its own:
- **0|2 and 2|0** resume pass 1's deferred rejects, which are located when pass 1 located
  them, so each is one genome comparison (`Extend_Located`).
- **1|1 cases (a) and (b)** start from the first and last thirds' ranges, kept from the
  pass-1 half searches (`OP_SEARCH_KEEP`). When pass 1 located the read within that third,
  the kept range is a located single row: one genome comparison. That is the same
  condition as settling.
- **Failure depths** already narrow the 1|1 windows, or skip a case outright
  (`Window_Set`, "skipped" in the counts).

`bi_mm` needed settling because it searched each pass-2 seed afresh. Nothing in settling
is bidirectional. Locating from either row (batmis2) only finds S a few bases sooner, so
slightly more kept ranges are located.

Steps per read, batmis2 `-n 2 -m 0 --batch 1`, 50k reads, 500 Mb, `--lookup 11`
(`BATMAN_BATCH_STATS`):

| stage | starts from | 100 bp | 32 bp |
|---|---|---:|---:|
| pass 1 (exact, 0\|1, 1\|0) | half searches | about 31 | about 25 |
| 0\|2, 2\|0 | pass 1's rejects | 1.5 | 1.5 |
| 1\|1 (a) | kept first third | 3.6 | 49 |
| 1\|1 (b) | kept last third | 3.6 | 32 |
| 1\|1 (c) | new search of the middle third | 15.5 | 36 |
| total | | 55 | 143 |

What remains:
- **100 bp, case (c):** its seed, the middle third, is neither a prefix nor a suffix, so no
  pass-1 range or bound covers it. Its hits have a mismatch in the first third, so the
  prefix says nothing about where they start.
- **32 bp, all three cases:** the thirds are 10–11 bases, about 500 rows each on 500 Mb,
  never located. Each case branches with a one-mismatch window from a wide range: the
  short-seed limit also met with search schemes (`attic/bifmindex/README.md`).

## 2. Stretching seeds with failure depths

**The idea.** If P[0, j) does not occur (j ≤ LH), any hit differs from the read in
[0, j). A 1|1 hit has exactly one LH mismatch p, so p < j and the hit is exact on
[j, LH). Likewise, if the suffix from L − e does not occur, q ≥ L − e. When both fail,
every 1|1 hit is exact on [j, L − e), and cases (a), (b) and (c) become one search:
1. the seed [j, L − e) exact;
2. extend left over [0, j) with exactly one mismatch;
3. extend right over [L − e, L) with exactly one.

Starting in the middle suits batmis2's bidirectional ranges. batmis1 and batmis2 already
use j and e, but only as window bounds; the seeds stay the fixed thirds.

**Measured** (a probe build of batmis2 recording each read's failure depths and its steps
in the 1|1 cases, same reads):

| reads | 100 bp share | 100 bp 1\|1 steps/read | 32 bp share | 32 bp 1\|1 steps/read |
|---|---:|---:|---:|---:|
| no half failed | 26% | 29.9 | 54% | 118.8 |
| prefix failed only | 34% | 19.9 | 19% | 116.3 |
| suffix failed only | 26% | 22.7 | 19% | 117.2 |
| both failed | 14% | 16.3 | 8% | 116.9 |

For the reads where both failed, the seed length L − e − j:
- **100 bp:** 79% have 32 bases or more (1–7: 87 reads, 8–15: 261, 16–31: 1,158,
  32–63: 4,393, 64+: 1,149).
- **32 bp:** 7 bases or fewer (≤ 0: 446 reads, 1–7: 3,582), shorter than the thirds.

**Not implemented.** The idea applies only to reads where both halves fail.
- **100 bp:** these reads hold 10% of the 1|1 steps. Cutting their 16 steps to about 5 would
  save about 1.5 of 55 steps a read, about 3%.
- **32 bp:** nothing.
- **Where the 1|1 steps go:** mostly to reads whose halves both occur. They are typically
  reads with an exact or one-mismatch hit, still checked for 1|1 hits elsewhere at `-m 0`,
  and nothing bounds them.
- **Real data** with more errors or unmappable reads would have more both-failed reads.
  The gain stays limited to those reads (about 10 steps each at 100 bp).
