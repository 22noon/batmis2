# Columba against BatMis at exact matching (K = 0)

The cluster pilot showed Columba ahead of batmis1/batmis2 at K = 0. The hypothesis to test:
Columba samples the suffix array more densely (every 4 positions, BatMis every 8 rows), and
its FM index is generally faster.

## Summary

- **The pilot gap is start-up, not search.** On the pilot's E. coli set (50k reads) BatMis's
  search is 4.5× faster than Columba's (0.033 s against 0.148 s). Its total was larger
  because the pilot config gives it the depth-11 lookup table above 40,000 reads, and
  building that for both indexes takes 0.14 s, four times its whole search. With the default
  depth (6) BatMis's total is 0.078 s against Columba's 0.25 s.
- **On a large genome BatMis is 5–9× faster than Columba at K = 0**, with the same hits
  (500 Mb, 1M reads: search 1.20 s against 6.86 s on one core, 0.13 s against 1.24 s on 16
  threads). Columba uses 2–3.4× the memory.
- **The sampling hypothesis is right about Columba itself**: at SA sparseness 8 it is 1.4×
  slower, at 16 2.3× slower. And BatMis does need more steps to locate a row: 7.0 LF steps
  on average (row sampling every 8, unbounded) against Columba's 1.5 (text sampling every 4,
  at most 3). But it is not what makes Columba fast: Columba is the slower one.
- **The rank structure makes no difference.** Columba's interleaved rank9 bit vectors and
  BatMis's checkpoint-and-count Occ take the same time per LF step (within 3%) once misses
  overlap; both read two cache lines per Occ.
- **What makes BatMis faster:** batching reads (2.7×: BatMis without it, 3.27 s, is still 2×
  faster than Columba) and finishing single-row paths against the packed genome (2.1× with
  batching, 3.3× without; Columba's equivalent, its in-text verification, matters to it in
  the same way).
- **What BatMis could take from Columba:** text-position SA sampling would cut its locate
  walk from 7 to 1.5 LF steps. Estimated gain at K = 0, 100 bp: about 10% of search time
  batched, up to a quarter unbatched; more for all-hits runs on repetitive genomes, where
  locating dominates. It needs a new SA format, a marker bit vector and changes to the
  locate code; not attempted here.

## Set-up

- **Genome:** the 500 Mb random genome; BatMis index from `bwt_build_index`, Columba 2.0
  (Vanilla, 32-bit) indexes at SA sparseness 4 (its default), 8 and 16: 70 s and 3.0–3.5 GB
  to build each, 1.5–1.8 GB on disk (BatMis: 1.4 GB, both directions).
- **Reads:** 1M reads of 100 bp and 1M of 32 bp sampled from the genome with 1%
  substitutions.
- **Runs:** K = 0, first hit (BatMis `-m 1`, Columba `-a best -I 100`) and all hits (`-m 0`;
  `-a all -e 0`); SAM output; BatMis `--lookup 11`, batch 16 (its defaults otherwise).
- **Machine:** Core Ultra 7 165H. One P-core (`taskset -c 2`) or 16 threads; median of 3 runs.
- **Times:** "search" is each tool's own figure: BatMis `--timing` `search`, Columba's
  "mapping threads finished mapping in". "Total" is wall time from `/usr/bin/time`.
- **Same output:** all five configurations report identical hit sets (470,901 hits at
  100 bp, 785,815 at 32 bp; compared with `paper/normalize_hits.py`).

## The pilot

E. coli, 50k reads of 100 bp, K = 0 first hit, one core:

| | BatMis2, lookup 11 | BatMis2, lookup 6 | Columba |
|---|---:|---:|---:|
| seed table build | 0.14 s | 0.0002 s | 0.09 s (10-mer table) |
| search | 0.033 s | 0.041 s | 0.148 s |
| total | 0.21 s | 0.078 s | 0.25 s |

`config_pilot.yaml` sets `large_reads: 40000`, which turns on `--lookup 11` for its 50k-read
sets. That table pays off from about 500k reads (`lookup_table.md`), which is what
`config.yaml` uses.

## 500 Mb genome, K = 0

| reads | mode | threads | batmis1 search / total | batmis2 search / total | Columba s=4 | s=8 | s=16 |
|---|---|---:|---|---|---|---|---|
| 100 bp | first | 1 | **1.20** / 2.62 | 1.23 / 2.63 | 7.02 / 7.94 | 9.76 / 10.66 | 15.89 / 16.71 |
| 100 bp | all | 1 | **1.20** / 2.55 | 1.25 / 2.66 | 6.86 / 7.79 | 9.79 / 10.73 | 15.94 / 16.85 |
| 32 bp | first | 1 | **1.02** / 2.27 | 1.09 / 2.37 | 6.55 / 7.43 | 9.43 / 10.33 | 15.89 / 16.77 |
| 32 bp | all | 1 | **1.06** / 2.27 | 1.09 / 2.36 | 6.79 / 7.72 | 9.65 / 10.55 | 15.65 / 16.59 |
| 100 bp | first | 16 | **0.13** / 1.20 | 0.13 / 1.29 | 1.17 / 2.24 | 1.19 / 2.09 | 1.46 / 2.36 |
| 100 bp | all | 16 | **0.13** / 1.20 | 0.13 / 1.29 | 1.24 / 2.24 | 1.32 / 2.24 | 1.50 / 2.41 |
| 32 bp | first | 16 | **0.11** / 1.09 | 0.12 / 1.20 | 0.82 / 1.81 | 2.54 / 3.51 | 1.45 / 2.40 |
| 32 bp | all | 16 | **0.11** / 1.08 | 0.11 / 1.16 | 0.80 / 1.71 | 0.94 / 1.85 | 1.43 / 2.37 |

Seconds, median of 3. Peak memory: batmis1 546 MB, batmis2 935 MB (one core; 830–1,290 MB
on 16 threads), Columba 1,853 / 1,615 / 1,494 MB at s = 4 / 8 / 16. The 32 bp first-hit s=8
figure on 16 threads (2.54 s) is out of line with its neighbours; the machine was otherwise
idle, so it is probably Columba's thread scheduling; it does not change any conclusion.

On 16 threads both tools' totals are mostly index loading and output: BatMis's search is
0.13 s of 1.2 s.

## The factors one by one

**SA sampling: locate steps.** `locate_steps.cpp` walks 1M random rows of the BatMis index to
a sampled row:

| | BatMis (rows divisible by 8) | Columba s = 4 (text positions) |
|---|---:|---:|
| mean LF steps | 7.0 | 1.5 |
| worst | 103 | 3 |
| share needing more than 8 steps | 30% | 0 |
| time per locate, one at a time | 2,190 ns | (about 1.5 × 290 ns) |
| time per locate, 16 interleaved | 396 ns | |

batmis2's either-row rule (`--own-row` off) takes the first sampled row of either index,
which roughly halves the walk; batmis1 has only its own rows. Columba's sparseness sweep
(table above) shows how much locating costs Columba: +43% from s = 4 to 8 and +130% to 16.

**Rank structure.** Columba's layout (`BitvecIntl`, rank9: 4 bit vectors interleaved per 64
rows, counts in a separate array) added to `interleaved_occ.cpp` as `R9`; LF steps on the
500 Mb BWT (ns per step; 0 differences against `BWTOccValue` on 20M rows):

| | BatMis library (LIB) | Columba rank9 (R9) |
|---|---:|---:|
| 1 core, 1 walk | 264 | 287 |
| 1 core, 4 walks interleaved | 41.2 | 41.8 |
| 1 core, 16 walks | 16.0 | 15.5 |
| 22 threads, 16 walks each | 44.2 | 44.3 |

Both read two cache lines per Occ, so they cost the same; rank9 spends 4 bits per base
against BatMis's 2.25.

**Seeds and verification.** Columba, 100 bp, first hit, one core, s = 4 (search seconds;
nodes visited per read):

| `-i` (in-text switch) | 1 | 4 (default) | 50 |
|---|---:|---:|---:|
| `-K 6` | 5.98 (29.4) | 6.76 (27.8) | 16.86 (24.0) |
| `-K 10` (default) | 5.82 | 6.69 | 16.97 |
| `-K 12` | 5.75 | 6.74 | 16.90 |

The seed table's depth does not matter (a 12-mer table only adds 4.5 s of building). Leaving
the index early for the text matters a lot: verifying in the text from a range of 50 is 2.5×
slower than from 4.

BatMis2, the same reads (search seconds):

| | batch 16 (default) | batch 1 |
|---|---:|---:|
| reference loaded, lookup 11 | 1.20 | 3.27 |
| reference loaded, lookup 6 | 1.62 | |
| `--no-ref`, lookup 11 | 2.55 | 10.76 |

## Why BatMis is faster

1. **Batching: 2.7×.** BatMis interleaves 16 reads so that their cache misses overlap;
   Columba searches one read at a time per thread. Even unbatched, BatMis (3.27 s) is 2×
   faster than Columba (6.7 s).
2. **Finishing against the reference: 2.1× batched, 3.3× unbatched.** Once a path is a single
   row, BatMis checks the rest of the read against `.pac` as soon as it reaches a sampled row.
   Columba does the same with its in-text verification, and it matters to Columba as much (`-i`).
3. **Less work per read at K = 0.** Columba runs its search-scheme machinery even at K = 0
   (about 28 search nodes per read); BatMis's exact search is one backward pass from its
   lookup range.

What Columba does better is locate: 1.5 LF steps against 7.

## What BatMis could adopt

- **Text-position SA sampling.** Store the SA value of every text position divisible by s,
  with a bit vector marking those rows, so a locate takes at most s − 1 LF steps.
  - **Estimate at K = 0, 100 bp:** about 470k located reads × (396 − about 90) ns ≈ 0.14 s of
    1.2 s batched (12%); unbatched 470k × (2.2 − 0.45) µs ≈ 0.8 s of 3.3 s.
  - **More in all-hits runs on repetitive genomes**, where `Sa_Decode` dominates.
  - **Cost:** a new `.sa` format and builder output, a rank on the marker vector per step, and
    changes to `Sa_Decode`, the located paths, and batmis2's either-row locate.
  - Worth a stand-alone trial with `locate_steps.cpp` extended before touching the mappers.
- **Pilot config:** raise `large_reads` in `config_pilot.yaml` to 500,000, so the pilot uses
  BatMis's default lookup depth like the real config would for sets that small.

## Not covered

- A real, repetitive genome (human): many hits per read make locating a larger share for
  both tools. The paper benchmark's human runs will show it.
- K > 0, where Columba's search schemes and BatMis's mismatch search differ in kind.

## Files

- `locate_steps.cpp`: LF steps to a sampled row on a BatMis index.
- `interleaved_occ.cpp`: now also times Columba's rank9 layout (`R9`).

## Denser SA sampling in BatMis (implemented)

Following the locate-step result, `bwt_build_index` now chooses the SA sampling interval from
a memory budget for the samples of both indexes (`--sa-budget`, default 1 GB; `--sa-interval`
to force one): every row up to 134 Mb, every 2nd up to 268 Mb, every 4th up to 537 Mb, every
8th (as before) above. batmis1 and batmis2 read the interval from the index
(`index/README.md`). Same hits on every interval: `make test` runs every suite on both
samplings, and E. coli at K = 0–4, both mappers, with and without the reference, gives
identical SAM hit sets; original BatMis 3.00 also reads the dense index correctly.

**Effect**, 1M reads, `--lookup 11`, median of 3; search = batmis's own figure:

| genome (interval) | reads | K | mode | search, 1 core: dense / 8 | wall, 1 core | wall, 16 threads | peak memory, 1 core |
|---|---|---:|---|---|---|---|---|
| E. coli (1 vs 8) | 100 bp | 0 | first | 0.65 / 0.72 (−10%) | 1.35 / 1.43 | 0.56 / 0.55 | 110 / 79 MB |
| E. coli (1 vs 8) | 100 bp | 2 | all | 2.43 / 2.76 (−12%) | 3.44 / 3.82 | 1.07 / 1.08 | 161 / 130 MB |
| chr21 (1 vs 8) | 100 bp | 0 | all | 0.86 / 0.93 (−7%) | 2.36 / 1.73 | 0.79 / 0.68 | 463 / 152 MB |
| chr21 (1 vs 8) | 100 bp | 2 | all | 3.71 / 4.15 (−11%) | 5.12 / 5.49 | 1.56 / 1.51 | 523 / 211 MB |
| 500 Mb (4 vs 8) | 100 bp | 0 | all | 1.08 / 1.18 (−8%) | 2.62 / 2.53 | 1.39 / 1.21 | 1,412 / 935 MB |
| 500 Mb (4 vs 8) | 100 bp | 2 | all | 4.36 / 4.71 (−7%) | 6.68 / 6.84 | 2.41 / 2.30 | 1,464 / 987 MB |
| 500 Mb (4 vs 8) | 32 bp | 2 | all | 22.8 / 21.9 (+4%) | 24.9 / 23.8 | 4.99 / 4.71 | 1,466 / 989 MB |

batmis2 shown; batmis1 behaves the same. Over all 64 pairs (both mappers, K = 0 and 2, both
modes, 1 and 16 threads): search −4% to −14%, except 32 bp K = 2 all hits, −5% to +4%.

- **Search is 5–13% faster** in nearly every case, on one core and on 16 threads, about the
  estimate above (12%): single-row paths stop at once and are checked against the reference.
- **32 bp, K = 2, all hits gains nothing** (up to 4% slower): short reads leave little walk to
  save, and the SA table, twice as large, is read for each of the run's million hits.
- **Wall time hardly improves.** Index loading grows with the samples (500 Mb: 0.33 → 0.50 s;
  chr21: 0.04 → 0.16 s), mostly allocating and filling the larger array rather than reading
  the file (0.09 against 0.04 s from page cache). On 16 threads, where search is a tenth of a
  second, the run is up to 0.2 s slower; on one core the gains and the load roughly cancel
  on short runs and the search gain wins on long ones (E. coli and chr21 K = 2).
- **Memory:** the SA samples grow by 8n/s − n bytes (500 Mb: +477 MB; chr21: +312 MB).

So the budget buys about a tenth of the search time for small and mid-sized genomes, where
runs are short anyway; human (interval 8) is unchanged.

**Update:** the mappers now map the SA samples instead of reading them (`mmap_loading.md`). The
500 Mb index then loads in 0.12–0.13 s at every 4th row and at every 8th, so the denser
sampling no longer costs load time, and every run is 0.1–0.5 s shorter.
