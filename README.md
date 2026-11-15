# BatMis: batmis1 and batmis2

Two short-read mappers on the same search, for comparison:

- **batmis1** (`src/batmis1.cpp`) searches a pair of one-direction FM indexes, the
  reverse index `rev<genome>` and the forward index `<genome>`. Each SA range lives in one
  index. A seed found in one direction is re-walked in the other index before it extends
  the other way, and a single row is located from its own index's sampled rows.
- **batmis2** (`src/batmis2.cpp`) is batmis1 with bidirectional SA ranges: each range also
  keeps the start of the same string's range in the other index. It uses this in two
  places:
  - **index switches** read the switched range off the seed, in O(1), instead of
    re-walking it (`--switch-walk` restores batmis1's walk);
  - **a single row is located** from either index's sampled row (`--own-row` restores
    batmis1's own-row locate).

  With both flags, batmis2's output is byte-identical to batmis1's. Without them, it finds
  the same hits; only which N hits `-m N` reports can differ.

Both report hits with up to `-n` mismatches (0–7), fewest first, as text or SAM/BAM.
They share the index format, I/O, threads (`-t`), batched search (`--batch`) and SA
decoding.

## Layout

| path | contents |
|---|---|
| `src/` | `batmis1.cpp`, `batmis2.cpp`; shared `fm_occ.hpp` (Occ), `genome_map.hpp` (sequences, N runs), `sam_writer.hpp` (SAM/BAM via htslib); `Makefile` |
| `index/` | the FM index library and builder (`bwt_build_index`), shared by both |
| `htslib-1.24/` | htslib, built from source for SAM/BAM output |
| `tests/` | the test suites (below) |
| `benchmark/` | the write-ups of how batmis1 was built and measured. They name batmis1 `batman_zero_min`, its name until 2026-11-05. `two_mismatch_pass1.md` covers what pass 1 tells pass 2 in `-n 2` (both mappers). |
| `attic/` | code kept for reference, out of the build (below) |

## Build and test

```sh
cd src
make            # batmis1 and batmis2 (builds index/ and htslib's library as needed)
                # then checks that their prefetches are in the binaries
make test       # both suites, about 5 minutes
make install    # batmis1, batmis2 and bwt_build_index to ../bin
```

An index is built once per genome with `index/bwt_build_index genome.fa`. Then:

```sh
batmis1 -g genome.fa -q reads.fq -n 2 -m 10 --sam -t 8 -o hits.sam
batmis2 -g genome.fa -q reads.fq -n 2 -m 10 --sam -t 8 -o hits.sam
```

**Batch size.** Both search 16 reads at a time (`--batch`), which does not affect the
output. To check that on another machine, run a workload a few times with `--batch auto`:
it times 4, 8 and 16, and prints the times and how far 16 is from the fastest. If another
size is clearly better there, set it with `--batch` or the environment variable
`BATMIS_BATCH` (`benchmark/prefetch_fix.md`).

**Tests** (`make test`; `make test-batmis1` and `make test-batmis2` run them separately).
Each mapper runs:
- **`test_located_hits.sh`:** the batched lanes against one read at a time, field by field;
  located paths against walked ones; and every hit set for `-n 1` to `-n 7` with all hits
  against brute force over the genome.
- **`test_sam_output.sh`:** SAM/BAM records checked against the genome; `-m 0` against
  brute force on both strands (`sam_oracle.py`); `--no-ref`, `--batch`, threads and BAM
  output; `Sa_Decode` on every row.

`make` and `make test` also run **`check_prefetch.sh`** on both binaries: the compiler
can drop the batched search's prefetches without any change in output
(`benchmark/prefetch_fix.md`).

batmis2 also runs **`test_batmis2.sh`**:
- with `--switch-walk --own-row`, byte-identical to batmis1, at `-n 0` to `-n 7` in text and
  SAM, batched, with threads and with `--no-ref`;
- as is, the same `-m 0` hits;
- `batmis2_verify` (built with `-DBATMIS2_VERIFY`) checks every range's other-index start and
  every other-row locate against the index. A deliberately broken `Other_Start` fails it.

## attic/

Kept for reference, not built.

- **`attic/bifmindex/`:** the bidirectional index work that led to batmis2: a hand-written
  bidirectional index (`bi_index.hpp`); exact, one- and two-mismatch searches (`bi_exact`,
  `bi_mm`); the range-chain `Locator`; the search-scheme prototype (`tools/scheme_eval.cpp`).
  Its `README.md` records the measurements, including the port to batmis2.
- **`attic/src/batmis2_range_chains.cpp`:** batmis2 with range chains recorded during the
  search (each anchored range by level; wide hits keep their chain), walked back before
  `Sa_Decode` for SAM/BAM output.
  - **Correct:** SAM byte-identical in 80 runs; recorded levels checked against rebuilt
    chains.
  - **Effect:** it resolved 93–100% of wide-hit rows, but decoding took the same time,
    within ±5% of `Sa_Decode` alone (`rx`, `rb`, `-n 1–2`, `-m 0`). The walk back (about
    10 BWT scans a hit, plus the levels under the lookup table) costs what `Sa_Decode`'s
    spans do.
  - **Removed 2026-11-05.**
- **`attic/src/`, `attic/tests/`:**
  - the structured exact search;
  - the inline-Occ variant of `fm_occ.hpp` and its test;
  - the opt-in `Chain_Decode` test;
  - SA-decode benchmarks.
- **`attic/alternate/`, `attic/awfm/`, `attic/doc/`, `attic/tmp/`:** the original BatMis
  sources, the BWA-index and AWFM variants, and their notes.
