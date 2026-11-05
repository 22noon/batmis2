# SAM / BAM output

## Summary

`batman_zero_min --sam` / `--bam` writes alignments through htslib
(`htslib-1.24/`, built as `libhts.a` by `src/Makefile`). `-o` ending in `.sam` or
`.bam` also picks the format; without `-o` the output goes to stdout. The
existing text output stays the default and is byte-identical to before.

> **Since threads** (`threads.md`): records are built by a `Sam_Encoder` per worker thread
> (`Encode_Sam_Read`, formerly `Write_Sam_Read`). `Sam_Writer::Write` then writes their
> bytes in order, as `sam_write1` would have, so the files are byte-identical.

- **Both strands:** every read and its reverse complement are searched.
  - A hit of the reverse complement is written with FLAG 16, SEQ
    reverse-complemented and QUAL reversed.
  - `--forward-only` searches the reads as given only.
- **One record per location:** CIGAR `<L>M` and MAPQ 255 (not computed).
  - **Tags:** `NM:i` (mismatches), `MD:Z` (reference bases at the mismatches)
    and `NH:i` (records written for the read).
  - **Order:** fewest mismatches first; the first record is primary, the rest
    secondary (0x100).
  - **No placed location:** a read with none, or with a non-ACGT base, gets one
    unmapped record (FLAG 4).
- **Header:** `@HD VN:1.6 SO:unsorted`, an `@SQ` per sequence, and `@PG` with
  the command line.
- **`-m`** counts a read's locations over both strands.
- **Dropped locations:** a location that leaves its sequence, covers an N, or
  lies on the short-reference padding sequence is not written.
  - It still counts toward `-m`, as if written. So a read can get fewer than
    `-m` records, or none (unmapped).
  - This keeps `-m` from decoding past the first `-m` rows (see "Decoding SA
    rows").

## Placing a hit (`src/genome_map.hpp`)

**The index's text** (`index/bwt_build_index.cpp`) is the FASTA sequences joined
with no separator. Each non-ACGT base is stored as `c`.
- `<genome>.ann.location` gives the sequences' names and lengths.
- `<genome>.N` gives each sequence's runs of non-ACGT bases.

`Genome_Map::Place(G, L)` turns a text position into (sequence, position), or
rejects it.

**Finding the sequence: a block table.**
- **Blocks:** the text is cut into blocks of 2^k bases, with k chosen so a
  block holds about one sequence start (k = ⌊log2(text / sequences)⌋, in
  [8, 24]).
- **`First[b]`:** the sequence holding block b's first base.
- **Lookup:** a hit in block b belongs to a sequence in `[First[b], First[b+1]]`.
  That is usually a single candidate; otherwise a binary search runs over that
  block's sequences only.
- **Sizes:** 3 Gb with 25 sequences needs about 180 entries; a million
  scaffolds, about 1.5M (6 MB).
- **The boundary check** is then one comparison with the sequence's end.

**N runs: a bitmap filter with an exact check behind it.**
- **The filter:** one bit per 4096-base block that holds any N, 92 KB for
  3 Gb. A read touches one or two blocks, so most hits are cleared with a bit
  test or two.
- **The exact check:** only hits touching a marked block search the sorted run
  list. The last run starting before the hit's end overlaps it iff it ends
  after the hit's start.
- **Why not one bit per base:** that would cost 375 MB for 3 Gb.

**The SAM header** lists only sequences that can hold a hit. Empty sequences
(SAM needs LN ≥ 1) and the padding sequence are left out, so a hit's tid is its
sequence's place in that list.

**`tests/test_genome_map.cpp`** checks `Place` against a brute-force scan of
every text position for L = 32, 75 and 100:
- a 610 kb reference with 8 sequences: N runs at starts, ends and middles, an
  IUPAC base, sequences of 90, 40 and 0 bases;
- a 2.76 Mb reference with 6,000 sequences, mostly 0–120 bases (many per
  block);
- the 1 Mb and toy genomes.

All agree.

**Dropped hits still count toward `-m` in the search.** So with `-m 1`, a read
whose first hit leaves its sequence can be reported unmapped although another
hit exists. With `-m 0` nothing is lost.

## Positions without the reference (`Mark_Sampled`)

**With the reference,** a single-row path is "located" at its first sampled SA
row: its genome position `Pos` is read there, and the rest is compared on
`.pac`.

**Without it (`--no-ref`),** the walk goes on to the end of the read. Before
this change, a hit's position then needed an SA walk at output (`BWTSaValue`:
up to 7 LF steps from the hit's row).

**The trick:** with SAM/BAM output, `.sa` is loaded even with `--no-ref`, and
`Track_Pos` is on.
- **At the first sampled row** a single-row path meets, `Mark_Sampled` records
  `Pos` from the sample, with `Locate`'s formula (or `NO_POS` if the read would
  run off the text).
- **The walk then goes on as before.** Every later step, mismatches included,
  follows the same occurrence, so `Pos` stays the read's start.
- **At output** such a hit costs nothing. Only multi-row ranges, and paths that
  never met a sampled row (e.g. just after an index switch), are decoded row by
  row.
- **Coverage:** `Mark_Sampled` replaces the no-reference branch at all 13
  sampled-row sites: per-read and batched searches, extensions and windowed
  extensions.
- **With the text output,** `Track_Pos` is off and nothing changes.

**Check:** `--no-ref` SAM is byte-identical to SAM with the reference in every
oracle run. Shifting the matched-base count by one in `Mark_Sampled` produces 38
differences in the matrix below.

## Strands and `-m`: zig-zag search

**Both strands of a read run in one lane, in zig-zag order** (`Zig_Zag`,
`Batch_Both_Strands`):
- 0 mismatches on the read, then on its reverse complement; then 1 mismatch on
  the read, then on the reverse complement; and so on.
- `-m` is counted over both strands, and the read stops as soon as it is
  reached.
- So a read whose reverse complement matches exactly stops after the exact
  stage, without searching the forward strand any deeper. The other way round,
  the reverse complement is never searched past the level where the forward
  strand filled `-m`.

**How the two-strand program is built:** generated from the one-strand program.
- **Groups:** it is cut at each `OP_STOP_IF_MAXHITS`, and after
  `OP_COUNT_EXACT`.
- **Order:** each group runs on strand 0, then a stop check, then on strand 1,
  then a stop check.
- **Steps carry their strand** (`Mm_Step::Strand`). The lane keeps the active
  strand's state in its usual fields, and `Mm_Set_Strand` swaps it, between
  steps only, with the other strand's (`Mm_Strand`): read, packed read,
  lists, exact count, failure bounds, kept ranges.
- **`Mm_Used`** counts both strands.
- **Budget guard:** a budget-limited `OP_EXTEND` is skipped once the hits
  reach `-m`, because a budget of 0 would mean "no limit".

`-n 2` runs in these groups: LH search with 0|0 and the exact count; 0|1; RH
search with 1|0; 0|2 and 2|0; 1|1 (a); 1|1 (b); 1|1 (c).

**Why the hits are still the best:**
- Each group holds one mismatch count. The first group holds only exact hits,
  which are never budget-limited.
- Within a strand the groups keep their order.
- So when the read stops after a group, no hit with fewer mismatches is left
  on either strand.
- **Skipping the reverse strand's half:** after a stop following the forward
  strand's half of a group, the reverse strand could only have added ties.
  Ties go to the forward strand, as in the final merge.
- **The merge stays:** `Write_Sam_Read` still sorts by mismatch count and cuts
  at `-m`, so `-m 0` output is byte-identical to searching the strands
  separately (54 runs on the test references).

**Memory accesses per read** (`-n 2 -m 1 --batch 16`, 20k reads, 500 Mb):

| reads | planted on | both strands, separately (`--forward-only` + reverse) | zig-zag |
|---|---|---:|---:|
| p32_mm0 | forward | 19.1 + 129.7 | 19.1 |
| p32_mm0 | reverse | 129.7 + 19.1 | 28.6 |
| p100_mm0 | forward | 24.1 + 28.2 | 24.1 |
| p100_mm0 | reverse | 28.2 + 24.1 | 33.5 |
| p32_mm2 | forward | 67.3 + 129.8 | 107.9 |
| p32_mm2 | reverse | 129.8 + 67.3 | 127.5 |
| p100_mm2 | forward | 42.7 + 28.2 | 61.5 |
| p100_mm2 | reverse | 28.2 + 42.7 | 63.6 |

- **Exact reads** cost the same as a one-strand search when they lie on the
  forward strand. On the reverse strand they cost one failed exact search
  more.
- **Two-mismatch reads** need both strands searched through levels 0 and 1
  before the stop, whichever strand they lie on. That is the price of a
  best-first answer.

**One engine for SAM:** all SAM/BAM output now runs on the batched engine.
`--batch 1`, the default, is one lane. The text output keeps the per-read
search functions.

## Correctness

- **`sam_oracle.py`:**
  - **every record:** SEQ (the read, or its reverse complement for FLAG 16),
    QUAL (reversed for FLAG 16), CIGAR, a window inside one sequence with no N,
    `NM` equal to the Hamming distance, `MD` recomputed from the genome, `NH`,
    and the primary / secondary flags;
  - **with `-m 0`:** each read's (sequence, position, strand) set equals brute
    force over both strands, at distance ≤ `-n`, windows inside one sequence
    with no N;
  - **with `-m 1`:** the primary is at the best brute-force distance.
- **The matrix** (`sam_matrix.sh`):
  - **References:** the multi-sequence reference, the 1 Mb random genome and
    the toy.
  - **Reads:** 200 each at L = 32, 50 and 100, either strand, 0–3 mismatches,
    including reads across sequence boundaries and over N runs, random reads,
    and reads with an N.
  - **Settings:** `-n` 0–3, `-m` 0, 1 and 3.
  - **Consistency:** each run is also compared, byte for byte apart from
    `@PG`, with `--no-ref`, `--batch 16`, both together, and BAM decoded by
    htslib's `test_view`.
  - **Result** (zig-zag build): 108 oracle runs. 107 pass, and in all 108
    the `--no-ref`, batched and BAM output matches `--batch 1` (432
    comparisons).
  - **The `-m` check** now also requires min(`-m`, valid hits) records for
    every read with no invalid window at distance ≤ `-n`.
  - **The one failure** is the `-m` limit described above.
    - **What happens:** at `-n 3 -m 1`, the read's best hit (2 mismatches)
      covers an N, so it is dropped. With `-m 1` the search had already
      stopped there, so its 26 valid 3-mismatch hits were never found, and the
      read is written unmapped.
    - **With `-m 0`** all 26 are written correctly.
- **500 Mb samples** (2,000 reads each of p32_mm2 and p100_mm2, `-n` 2 and 3,
  `-m 5 --batch 16`): every record checks against the genome.
- **Fault injection:** `MD` built without sorting the mismatches by position
  fails the oracle in every run with `-n 3`.
- **Text output:** byte-identical to the previous build (48 runs on 500 Mb:
  `-n` 0–3, both reference modes, batch 1 / 16).
- **Tests:** `make test-located`, `make test-structured` and the new
  `make test-sam` pass.
  - `make test-sam` runs `test_genome_map` and the SAM oracle on a generated
    reference: -n 0–3, -m 0 / 1, `--no-ref`, `--batch 16`, and BAM.

## Timing

Search and output seconds, mean of 2 runs, `-n 2 -m 1 --batch 16`, 500 Mb
genome, load average 2–2.8. p100_mm2 has 500k reads of 100 bp, p32_mm2 1M of
32 bp.

| reads | ref | output | search | output |
|---|---|---|---:|---:|
| p100_mm2 | on | text | 2.06 | 0.06 |
| | | SAM, `--forward-only` | 2.07 | 0.41 |
| | | SAM | 3.68 | 0.43 |
| | | BAM | 3.80 | 2.48 |
| | `--no-ref` | text | 5.13 | 0.07 |
| | | SAM | 6.94 | 0.45 |
| | | SAM without the position trick | 6.65 | 1.62 |
| p32_mm2 | on | text | 7.38 | 0.12 |
| | | SAM, `--forward-only` | 7.38 | 1.13 |
| | | SAM | 22.2 | 1.22 |
| | | BAM | 22.6 | 3.63 |
| | `--no-ref` | SAM | 22.4 | 1.31 |
| | | SAM without the position trick | 21.7 | 3.35 |

- **Writing SAM** costs 0.8–1.1 µs per read.
- **BAM** adds single-threaded BGZF compression: 2.4 s per 500k reads.
  `hts_set_threads` would spread it over cores.
- **The position trick** (`--no-ref`) cuts SAM output from 1.62 to 0.45 s
  (p100) and from 3.35 to 1.31 s (p32), because fewer SA walks are needed.
  - Reading the SA sample during the search costs 0.3–0.7 s.
  - The net saving is about 0.9 s per run: 11% of p100's total, 4% of p32's.
- **The reverse-complement search costs more than writing.** The rows above
  search each strand on its own, so a strand with no hit runs every stage.

**Zig-zag** (search + output seconds, mean of 2 runs, `-n 2 -m 1`, SAM,
load average 2.6):

| reads | strands | batch | forward only | strands separately | zig-zag |
|---|---|---:|---:|---:|---:|
| p100_mm2 (500k) | forward-planted | 1 | 3.93 | 6.48 | 5.79 (−11%) |
| | | 16 | 2.71 | 4.39 | 3.92 (−11%) |
| | 50/50 mix | 1 | (3.30) | 6.47 | 5.84 (−10%) |
| | | 16 | (2.28) | 4.34 | 4.04 (−7%) |
| p32_mm2 (1M) | forward-planted | 1 | 12.3 | 32.4 | 19.5 (1.66×) |
| | | 16 | 9.26 | 24.0 | 14.5 (1.65×) |
| | 50/50 mix | 1 | (16.7) | 32.5 | 20.9 (1.55×) |
| | | 16 | (12.2) | 24.1 | 15.6 (1.54×) |

(Forward-only on the mixes misses half the reads; it is shown for scale only.)
- **32 bp:** zig-zag removes most of the second strand's cost. 2-mismatch reads
  still search both strands through 0 and 1 mismatches.
- **100 bp:** the gain is smaller, because there a failed search is short
  either way.

## Decoding SA rows

A hit range lists SA rows. A row's text position comes from walking it with LF
(`BWTSaValue`) until it reaches a sampled row (every 8th, `saInterval`).
- **The walk is serial:** about 7 dependent cache misses per row, so it matters
  when hits are wide: repeats at `-m 0` or a large `-m`.
- **Located hits need no walk:** with the reference, or with `--no-ref`'s
  `Mark_Sampled` trick, a hit's position is known already.

**Benchmark genome** (`rep200`, scratchpad): 10 × 20 Mb, 40% repeat copies
(300 bp / 6 kb families at 5–20% divergence), and 200 N runs of up to 50 kb.

**Finding: poly-C from N runs.**
- **Where the rows come from.** The index stores N as `c`, so an N run is a
  poly-C run. A few all-C reads (3 per 2000) match millions of rows there.
  Every one of those rows is decoded and then dropped by the N filter.
  `hg38` has about 150 Mb of N.
- **Why walks stall there.** In a homopolymer run, LF shifts a whole block of
  rows by a constant stride. When that stride is a multiple of 8, the same rows
  miss the samples step after step: 148 LF steps per row instead of about 8.

### `Sa_Decode` (`src/batman_zero_min.cpp`)

The SAM block now runs in two phases:
1. collect each read's rows to decode, capped at `-m` locations
   (`Sam_Requests`);
2. decode them all together (`Sa_Decode`), then write.

**Units of work, stepped one round at a time:**
- **Span** (a request of at least 16 rows): rows Lo..Hi plus a mask of the
  resolved classes `(row − Lo) mod 8`.
  - Each step reads its sampled class straight from `saValue`.
  - If every row has the same BWT character (`BWTAllOccValue` at both ends), LF
    shifts the whole span at O(1) cost, however long it stalls.
  - Otherwise it splits into runs of one character; runs of 16+ rows stay
    spans.
- **Group:** the unresolved rows of a narrower request or run, in row order.
  - Sampled rows are read off.
  - The rest map to up to 4 child groups: one `BWTAllOccValue`, then counting
    the BWT characters between neighbouring rows (a fresh lookup across a gap
    over 64 rows).
- **Walk:** a single row, stepped like `BWTSaValue`. The walks of a round are
  interleaved with prefetching, so their cache misses overlap.

**Replay** (`tests/bench_sa_decode.cpp`; decode requests recorded with
`BATMAN_DECODE_DUMP`, 2M rows each, seconds):

| workload | `BWTSaValue` | walks only | groups | groups + spans |
|---|---:|---:|---:|---:|
| random SA ranges (1–1000 rows) | 4.14 | 2.64 | 2.24 | 2.3–2.7 |
| poly-C hit (100 bp) | 16.3 | 18.1 | 7.6 | 1.2–1.6 |
| rep200 32 bp hits | 18.6 | 19.5 | 7.9 | 1.3–1.6 |

- The span threshold (8, 16 or 64 rows) and the group gap threshold (16, 64 or
  256) are within noise; the defaults are 16 and 64.
- Giving groups no single-row walks costs about 50% on random ranges.

**End to end** (rep200, `-n 0 -m 0`, SAM, seconds of output):

| reads | previous | `Sa_Decode` |
|---|---:|---:|
| 3 all-C 32 bp reads (14.5M rows) | 101 | 10 |
| same, `--no-ref --batch 16` | 94 | 11 |
| 2000 × 100 bp, `--batch 16` (53M rows) | 359 | 36 |

The SAM output is byte-identical to the previous build on all of these, and on
108 runs on gm/g1m/toy.

### Dropped locations count toward `-m`

**Before:** when the first `-m` rows were N-run locations and got dropped, the
writer decoded further rows one at a time until it had `-m` records. On rep200,
`-n 1 -m 10` spent **656 s** walking millions of poly-C rows, all of them
dropped.

**Now:** a dropped location counts as written, so only the first `-m` rows are
ever decoded.
- The same run takes **0.22 s** and gives identical output: those reads stay
  unmapped either way.
- On the gm/g1m/toy matrix, one read changed: an N-run hit now uses up `-m 3`
  (`gm L32 -n 3`).
- **Follow-up:** replace the genome's N's with random bases at index time to
  make such hits rare.

### Decoding from the hit's range chain (`Chain_Decode`, opt-in: `-DSA_CHAIN`)

**Idea (user's).** A hit range A_L is the SA range of the hit string P' (the
read with its recorded mismatches) in its index.
- **Rebuild:** exact backward search of P' rebuilds the chain A_1 … A_L (2 Occ
  lookups per step), whatever route the search took.
- **Ancestors:** LF keeps order, so the t-th row of A_{k−1} whose BWT character
  is the k-th one consumed maps to row `s_k + t` of A_k. Every hit row therefore
  has an exact ancestor in each A_k: a different row each time, with SA value
  (hit's) + L − k. A sampled ancestor resolves the hit row with no walk.
- **Stretches** (steps that keep the width): the ancestors are contiguous. Each
  distinct `s_k mod 8` covers another eighth of the hit, so a stretch is taken
  in one go.
- **Boundaries** (width changes): one word-wise scan of A_{k−1}'s BWT (a
  popcount per 16 characters) finds the ancestors by counting.
- **Tail shortcut:** widths never grow along the chain, so once one equals the
  hit's width, a step needs only one rank.
- **Not the nested offset form.** `start + Σ_{b<c} cnt_b` is for appending
  within the same SA (bidirectional right extension). There the hit's rows keep
  their row numbers, so it finds no new samples.

**Limits:**
- Backtracking stops when every row is resolved, when a scan would pass
  4 × rows left + 4096, or after 12 scans without resolving 1/8 of the rows
  left.
- A hit is first tried on its first 256 rows, which is cheap: their ancestors
  come first in each scan.
- A hit resolving under 7/8 goes to `Sa_Decode` whole. Poly-C hits do this:
  their ancestors stall as their walks do.

**Correctness:** `tests/test_chain_decode.cpp` (in `make test-sam`) checks every
row the chain resolves against `BWTSaValue` and against the genome text. It
uses random hits with up to 2 mismatches on both strands and both indexes:
5M+ rows, 0 failures. SAM output is byte-identical with and without the chain.

**Measured: no gain, so the chain is off by default.**
- **Coverage is high:** it resolves 95–100% of the rows of wide hits.
- **The time is not better.** Decode phase only (`decode: seconds`, stats
  build, best of 3, 50k reads):

| genome / reads | `Sa_Decode` | with chain |
|---|---:|---:|
| `rx` (exact-ish repeats), 32 bp `-n 0 -m 0` | 0.24 | 0.30 |
| `rx`, 100 bp `-n 2 -m 0` | 0.66 | 1.21 |
| `rb` (copies differ just upstream of the reads) | 0.016 | 0.017 |

  - **Spans already do the same thing forward.** Exact copies share their
    upstream sequence too, so the forward walk shifts all of a hit's rows as one
    span.
  - **Where the copies diverge just upstream** (`rb`), the rows fall into
    groups that still share lookups.
- **Decoding is no longer the cost:** it is about 0.1 µs per row, against about
  1 µs to write each SAM record.
- **The rebuild costs about 2L lookups per hit.** The tail shortcut saved only
  about 4 per hit on `rx`, because copies that differ inside the read keep
  narrowing the range until near its end.

### Index-format options (not done)

- **Sampling by text position** (every k-th position, plus a marked-row
  bitvector with rank) bounds a walk at k − 1 steps. It also removes the
  homopolymer stall, but needs a new `.sa` format.
- **`SAINTERVAL` 4** halves walks and doubles `.sa`.
