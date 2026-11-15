# Recap: the port against the original BatMis 3.00

> **Naming (2026-11-05):** `batman_zero_min` here is now **batmis1** (`src/batmis1.cpp`);
> its bidirectional version is batmis2 (`src/batmis2.cpp`, see the top-level `README.md`).
> Code these write-ups mention but outside batmis1 (the structured search, the inline Occ,
> `Chain_Decode`, the BWA variant) is in `attic/`.

> **Prefetches (2026-11-15):** the batched timings here and in the other write-ups were
> measured with binaries from which the compiler had dropped the prefetches. Fixed in both
> mappers; batched search is 35–42% faster (`prefetch_fix.md`). The numbers below are not
> redone.

## What changed

Each item has its own write-up in this directory.

**Search:**
- **Located single-row paths** (`pac_singleton.md`): a path that reaches a
  sampled SA row is finished against the packed genome. Search is 21–76%
  faster, depending on read length.
- **Start-of-search lookup table** (`lookup_table.md`): `--lookup k`, with
  depth 10–11 for 0.5–2M reads.
- **Batched reads** (`batched_reads.md`, `branch_prefetch.md`): `--batch B`
  interleaves B reads so their cache misses overlap. Exact search is 1.7–2.6×
  faster, with byte-identical output. The default, `--batch auto`, times
  B = 1–16 on the first blocks and keeps the fastest: the best B depends on the
  core and the workload (B = 1 at 100 bp `-n 4` on a P-core).
- **Inline Occ, tried and left off** (`inline_occ.md`): counting the BWT
  characters with popcounts instead of the library's decode tables gives identical
  output and 6–9% fewer instructions, but no fewer cycles, and `-n 1` is 10–15%
  slower. The library Occ stays the default; `-DBATMAN_INLINE_OCC` selects it.
- **Threads** (`threads.md`): `-t N`, byte-identical output for any N.
  - **Speed:** 2M reads at `-t 22` run 6.3–9.7× faster end to end than one thread
    (6.3× for SAM, where writing the records is still serial).
  - **Input:** a reader thread reads the next block while this one is searched, and
    the workers encode the reads.
  - **Lanes:** threads lower the best B to about 2–4, because hyperthreads share a
    core's misses in flight.
- **Gzip input** (`threads.md`): `-q` reads plain or gzip FASTA/FASTQ, or stdin.
  Inflating costs about 3 s per 2M reads, on a thread of its own.
- **Parallel SAM/BAM records** (`threads.md`): the workers build the records, and the
  main thread only writes the bytes.
  - SAM output at `-t 22` fell from 6.5 s to 1.1 s per 2M reads.
  - BAM stays limited by zlib compression.
- **Rank lookup vs direct count** (`branch_threshold.md`), and occurrence
  checkpoint spacing (`occurrence_checkpoint_spacing.md`): the defaults are
  confirmed.
- **Mismatch search, by halves and their parts:**
  - `-n 2` (`two_mismatch.md`): rolling reject lists, and thirds for 1|1;
  - `-n 3` (`three_mismatch.md`): quarters, with kept seed ranges;
  - `-n 4` (`four_mismatch.md`): for 2|2, LH's quarters when both mismatches
    share one, else LH's thirds; the other new cases resume existing reject
    chains;
  - `-n 5` (`five_mismatch.md`): the mirror on RH for 3|2. One rule picks every
    seed, and each is at least as long as BatMis's.
  - `-n 6`, `-n 7` (`seven_mismatch.md`, proof of concept): a half's three-mismatch
    seeds, for 3|3 and 4|3.
    - Smallest seed: an eighth of the read (12 bases at 100 bp, 4 at 32 bp).
    - At 100 bp, `-n 7` costs 2.5× `-n 5` per read.
  - Each hit is found once, and per read and batched give the same output.

**Output:**
- **SAM / BAM** (`sam_output.md`): htslib output, both strands, NM / MD / NH.
- **`Sa_Decode`:** SA rows are decoded in groups and spans instead of one LF
  walk per row. This removes the poly-C (N-run) stall.
- **`-m` counts dropped locations:** locations dropped because they cover an N
  run still count toward `-m`, so `-m` never decodes past its first rows.
- **`Chain_Decode`** (opt-in, `-DSA_CHAIN`): range-chain backtracking. Correct,
  but no faster, so it is off by default.

## Setup

- **Machine:** shared Linux server. Each run pinned to one core; load average
  5–8, so single runs are noisy, about ±10%.
- **Original:** BatMis 3.00 from `BatMis-3.00.tar.gz`, `./configure; make`
  (with `-fcommon` added for current GCC). Its time is `batman` (binary hits)
  plus `batdecode` (SAM). It reads the port's index files unchanged.
- **Port:** `batman_zero_min --batch 16 --sam`.
- **Times:** wall seconds, including index loading. The original starts in about
  1.45 s, the port in 0.7 s (`-n 0`) to 1.2 s (`-n ≥ 1`).

## Mismatch search (500 Mb random genome)

**Reads:** 100k per length: 25k each with no mismatch, one in LH, one in RH, and
two.

**Same answers:**
- With all hits (`-mall` / `-m 0`), the original and the port report identical
  hit sets for every read, in all 30 combinations of length, `-n` and reference
  mode.
- With `-m 1`, the same reads map with the same number of hits; only the choice
  among equally good hits differs.

Original (`batman` + `batdecode`) → port, with the port's `--no-ref` time in
brackets:

| length | `-n` | `-m 1` | all hits |
|---|---|---|---|
| 100 bp | 0 | 4.5 → 1.7 (2.4) | 4.5 → 1.7 (2.8) |
| 100 bp | 1 | 6.5 → 2.6 (5.4) | 8.6 → 3.2 (7.5) |
| 100 bp | 2 | 7.5 → 2.7 (6.1) | 12.7 → 4.8 (10.1) |
| 100 bp | 3 | 7.0 → 2.5 (5.3) | 14.7 → 6.0 (13.0) |
| 100 bp | 4 | 7.1 → 2.8 (3.5) | 17.9 → 8.3 (12.0) |
| 75 bp | 0 | 3.4 → 1.6 (2.2) | 3.5 → 1.6 (2.2) |
| 75 bp | 1 | 5.3 → 2.6 (4.2) | 7.1 → 3.2 (5.3) |
| 75 bp | 2 | 6.1 → 2.7 (4.4) | 10.1 → 4.6 (8.2) |
| 75 bp | 3 | 6.2 → 2.8 (4.4) | 12.5 → 6.0 (10.5) |
| 75 bp | 4 | 5.6 → 2.5 (3.1) | 21.4 → 12.7 (15.0) |
| 32 bp | 0 | 2.8 → 1.7 (1.6) | 2.5 → 1.6 (1.6) |
| 32 bp | 1 | 3.5 → 2.7 (2.9) | 4.6 → 3.5 (3.7) |
| 32 bp | 2 | 5.4 → 3.4 (3.8) | 28.1 → 11.5 (14.1) |
| 32 bp | 3 | 5.3 → 3.7 (4.6) | 49.5 → 27.2 (40.9) |
| 32 bp | 4 | 4.1 → 3.1 | 279.5 → 216.5 |

The `--no-ref` column was measured in a separate run. For `-n 0..3`, my `-n 4`
tests overlapped part of that run, so those `--no-ref` times run a little high.

**Findings:**
- **75–100 bp:** the port is 2–3× faster with the reference, and 1.3–2× faster
  without it.
- **32 bp:**
  - `-n 0/1` are dominated by index loading;
  - `-n 2` all hits is 2.4× faster;
  - `-n 3` / `-n 4` all hits are only 1.3–1.8× faster. Their 5–8-base seeds
    match too many rows; this is the next thing to improve.
- **Memory:** the port's peak is 0.5–0.9 GB with the reference (SA samples,
  `.pac`, both indexes), 0.4–0.8 GB with `--no-ref`, and 1.45 GB for 32 bp
  `-n 4` all hits. The original's is 0.28 GB in `batman`, then 0.77 GB in
  `batdecode`.

## Decoding wide hits (rep200: 200 Mb, 40% repeats, 200 N runs)

The reads match millions of rows. SAM goes to `/dev/null`; `>300` means the run
hit the 300 s limit.

| reads | `-n` | `-m` | port now | port before `Sa_Decode` | original |
|---|---|---|---:|---:|---:|
| 8 all-C reads, 100 bp | 0 | all | 22.1 | 213.6 | 217.8 |
| 8 all-C reads, 100 bp | 0 | 10 | 0.21 | 184.3 | 0.36 |
| 2k reads, 32 bp | 0 | all | 11.5 | 116.3 | 134.6 |
| 2k reads, 32 bp | 2 | all | 82.9 | >300 | >300 |
| 2k reads, 100 bp | 0 | all | 25.2 | 268.9 | >300 |
| 2k reads, 100 bp | 2 | all | 95.2 | >300 | >300 |
| 2k reads, 100 bp | 3 | all | 94.5 | >300 | >300 |
| any of the above | any | 10 | 0.18–0.48 | 120–280 or >300 | 0.4–0.7 |

- **All hits:** `Sa_Decode` decodes about 10× faster than both the original's
  `batdecode` and the port's old one-walk-per-row decoder.
- **`-m 10`:** counting dropped N-run locations toward `-m` removed the old
  port's long stalls. Both tools now take under a second.

## Correctness

- **Tests:** `make test-located`, `make test-sam` and `make test-structured`
  pass, with `-n` up to 4.
- **Brute force:** for `-n 4`, all genome positions within distance 4 on a 1 Mb
  random genome and a 2 Mb repetitive one, at 32–100 bp. No location is
  reported twice with the same mismatch pattern.
- **Regression:** adding `-n 4` left `-n 0..3` output byte-identical (96 runs).
