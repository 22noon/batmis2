# Resolving a small SA range on the reference

## Summary

**The question.** In the exact search the range [S, E] has matched read[0..m),
with k = L − m bases left. Can the rest of the search be replaced by checks on
the reference (`.pac`) plus a sampled inverse SA (ISA)? Two cases:
- **General:** a range of w rows.
- **Width 2 with a sampled End:** the End's occurrence matches the rest of the
  read, so only Start is left to settle.

**Correct.** The final range follows from the reference alone:
- LF keeps the order of rows that share a preceding character. So the rows of
  [S, E] whose occurrences match the rest of the read map, in order, onto the
  final range.
- That range holds exactly the c matching rows and ends on
  E_f = ISA(P_last − k), where P_last is the text position of the last
  matching row.
- So a range is resolved completely by three steps:
  1. locate each row (LF with its own characters to a sampled row, at most 7
     steps);
  2. compare each occurrence's rest with the `.pac`;
  3. one ISA lookup (at most 7 LF steps).
- No match means no hit, and one match is a located hit at its genome position.

**Checked.** Results are identical to today's search for every read:
- **Genomes:** the 500 Mb random genome (r32, p100_mm0, p75_mm1L) and a 100 Mb
  repetitive genome (47% repeats; exact 32, 75 and 100 bp reads and
  one-mismatch 100 bp reads).
- **Lookup depths:** 6 and 11.
- **Variants:** all of them (W2E, W2, W4, W8, W16).
- **Assertions:** every resolution was checked against the range walked to the
  end (final rows and match count). There were no failures.

**What it does not do.** Two matching end rows say nothing about the rows
between them, so the range size is not kept. Each row must be checked, or the
range walked.

**Faster?** Only the width-2 case with a sampled End (W2E) is close to break-even,
and it wins only on short reads in a repetitive genome. Resolving wider ranges is
always slower.

ns per read, all reads, best of 3:

| genome, reads, lookup | today | W2E | W2 | W4 | W8 |
|---|---:|---:|---:|---:|---:|
| repetitive, x32_mm0, 6 | 2177 | **1838** | 2073 | 2521 | 3532 |
| repetitive, x32_mm0, 11 | 1894 | **1628** | 1862 | 2805 | 3000 |
| repetitive, x75_mm0, 11 | 2114 | 2081 | 2345 | 3091 | 3796 |
| repetitive, x100_mm0, 11 | 2319 | 2173 | 2446 | 2793 | 3705 |
| repetitive, x100_mm1, 11 | 2059 | 2162 | 2611 | 2963 | 4622 |
| 500 Mb, r32, 11 | 1873 | 1880 | 2947 | 4057 | 5838 |
| 500 Mb, p100_mm0, 11 | 2897 | 2948 | 3813 | 4939 | 6756 |
| 500 Mb, p75_mm1L, 11 | 1714 | 1727 | 2816 | 4036 | 5777 |

The repetitive-genome sets are 200k reads, and their timings vary by ±10%
between groups. The W2E gains of 7–16% on x32 and x100 are larger than that, but
the x75 and x100_mm1 differences are within it.

## Why it does not pay more

- **Resolution costs about one locate per row.** That is ~7 dependent steps,
  averaging 7.9 per width-2 resolution with E known (only S to locate) and 14
  when both rows are located. Then come one `.pac` compare per row and, for two
  or more matches, ~3.5 ISA steps. Even with the rows' locates overlapped, Wn
  costs grow with width: 20, 36 and 56 locate steps for W4, W8 and W16.
- **Today is cheap in the same place.** A range of 2–16 rows usually narrows
  to one row within a few bases. Today's located path then walks to the next
  sampled row and compares once, which is about one locate. So in the common
  case resolution repeats today's work, plus locates for rows that would have
  dropped out anyway.
- **Repeats are where resolution saves the walk.** These are ranges that stay
  multi-row to the end, where today walks every remaining base. But on the
  random genome they are wider than 16 rows (p100 reads ending multi-row take
  14.5 µs today and no variant helps), so locating each row costs more than
  walking. Ranges of exactly 2 copies that persist, as for x32 reads in the
  repetitive genome (13% end multi-row), are the case W2E wins.
- **Reads without hits gain nothing.** For one-mismatch reads (x100_mm1,
  p75_mm1L), the check finds no match, and the rest was going to fail quickly
  anyway.

## The width-2 case in detail

With E sampled and E's rest matching:
1. E_f = ISA(P_E − k).
2. S is located (at most 7 steps; 7.9 steps per resolution on average,
   including the steps where S was already sampled) and its rest compared.
   - Match: the final range is [E_f − 1, E_f].
   - No match: the hit is E's genome position, as a located hit, with no ISA
     needed.
3. Resolution is complete; there is no walk.

How often W2E fires:
- on 2–6% of reads;
- it skips 14–84 bases;
- most triggers end with one match (0.82–0.99), and the rest with two (an ISA
  lookup of 3.4–3.6 steps).

## Width 2 without an inverse SA: End on the reference, Start as a single row

**The setup.** The range has 2 rows, E is sampled, and E's occurrence matches
the rest of the read on the `.pac`. Then Start follows the normal single-row
extension:
- its BWT character is checked against the read at each base;
- in `W2S+loc`, at a sampled row its own rest is also compared on the `.pac`,
  as the located path does.

Exactly the rows that match survive, in order, so:
- **Start drops out:** only E's occurrence matches. The hit is E's genome
  position, reported at once as a located hit.
- **Start survives to the end:** both rows match, and the final range is
  [S_f, S_f + 1], with no inverse SA needed. Start's final row comes first and
  the range has exactly two rows.

**Correct.** Every read gave the same genome positions as today, and every
trigger was checked against the range walked to the end. There were no
failures, on both genomes and at lookup 6 and 11.

**Speed: neutral**, within about ±5% everywhere. The table shows ns per read at
lookup 11:

| reads | today | W2S | W2S+loc | triggers | Start drops out | Start steps |
|---|---:|---:|---:|---:|---:|---:|
| repetitive x32_mm0 | 2123 | 2053 | 2082 | 4.9% | 82% | 1.6–1.8 |
| repetitive x75_mm0 | 2672 | 2773 | 2602 | 6.0% | 96% | 2.6–3.9 |
| repetitive x100_mm0 | 2844 | 3050 | 3113 | 5.9% | 99% | 2.2–3.9 |
| 500 Mb r32 | 1884 | 1857 | 1879 | 2.2% | 99% | 0.5 |
| 500 Mb p100_mm0 | 2976 | 3079 | 3051 | 4.3% | 98% | 1.1–1.4 |
| 500 Mb p75_mm1L | 1682 | 1660 | 1668 | 0.0% | 61% | 19–21 |

Why it is neutral:
- **Start usually drops out within 0.5–4 bases.** Today the range narrows to
  E's single row at that same base.
- **What today does next:** walks E's row to its next sampled row (on average
  about 4 more steps, at most 7) and compares once on the `.pac`.
- **What W2S saves:** only those few steps, on the 2–6% of reads that trigger.
- **What W2S costs:** the End check (an SA sample + a `.pac` compare) at every
  width-2 range with a sampled End, including the ones whose rest does not
  match.

The two roughly cancel. When Start survives (both copies match, 1–18% of
triggers), W2S walks Start to the end. That is the same walk as today, one step
per base.

## The extension (mismatch) search

- **Exact-remainder nodes** (Mismatches == Count: the 0|0 pass and the
  threshold branch of `Extend_Right`/`Extend_Left`): the resolution is the
  same, with the same economics.
- **Nodes that may still take a mismatch:** the argument fails. A matching row
  there is one hit, not a boundary of the final ranges, because the remainder
  can branch.

## Costs and recommendation

- **Inverse SA memory:** n/8 × 4 bytes per index, the same as the SA samples:
  238 MB for 500 Mb, ~1.4 GB per index for 3 Gb.
- **Building it:** one LF walk over the text took 144 s for 500 Mb and 12 s for
  100 Mb. At index build it could instead be written straight from the full SA
  (`BWTSaveInverseSa`, `index/BWTConstruct.c`).
- **Recommendation:** not worth adding to `batman_zero_min` on these numbers.
  - **The best case** is W2E on short reads of a repetitive genome: 14–16%
    faster on x32.
  - **The cost:** a second sampled array per index and a new index file.
  - **Elsewhere** it is neutral or slower.
  - **What to revisit:** the case worth re-measuring is a real genome (human),
    whose many 2–3-copy segmental duplicates would trigger W2E more often than
    this synthetic one.

## How it was measured

`benchmark/sampled_end.cpp` (stand-alone, nothing in `src/` changed):
- **today:** `Search_Right`, with the reference on.
- **W2E:** resolve when the range has 2 rows and a sampled End whose rest
  matches.
- **Wn:** resolve the first time the range has 2..n rows (n = 2, 4, 8, 16).
  Rows are located together, round-robin, so their misses overlap.
- **Results compared** per read as genome positions (a located hit's `Pos`, or
  each final row's SA value). Every resolution was asserted against the range
  walked with `Get_SARange`.
- **Timing:** reads are also split by whether today's final range has one row
  or several.
- **Sampled ISA:** built by one LF walk and cached. Checked as ISA(SA(r)) = r on
  100k rows per genome.
- **Repetitive genome:** `scratchpad/rep/scripts/gen_rep.py`. It has 100 Mb,
  with 45% interspersed copies of 20 families of 300 bp and 5 of 6 kb at 2–20%
  divergence, plus 3% tandem repeats. Reads come half from repeat copies and
  half from unique sequence.
- **An earlier version** of this experiment assumed the range size is kept when
  both end rows match. On a 49 kb repetitive toy genome, 143 of 20k reads came
  out wrong, and the assumption was dropped.
