# batmisAWFM: AwFm-based Zero-Mismatch Search

This document explains how to build and run the AwFm-backed program variant and its index builder.

## What This Adds

- `batmisAWFM`: zero-mismatch search binary using AwFm `.awfmi` index files.
- `batmisAWFM_index_builder`: FASTA -> `.awfmi` index builder utility.

Note: This flow uses native AwFm index files and does not directly consume BatMis `rev*.bwt` / `rev*.fmv` files.

## Prerequisites

- GCC/G++ with C++17 support
- CMake
- OpenMP-capable toolchain

## 1) Build AwFm dependency

From the AwFm repo:

```bash
cd ~/git/batmis/AvxWindowFmIndex
git submodule update --init --recursive
cmake -S . -B build
cmake --build build -j"$(nproc)"
```

This produces required libraries under:

- `~/git/batmis/AvxWindowFmIndex/build/build/libawfmindex_static.a`
- `~/git/batmis/AvxWindowFmIndex/build/build/libdivsufsort64.a`
- `~/git/batmis/AvxWindowFmIndex/build/build/libfastavector_static.a`

## 2) Build batmisAWFM and index builder

From this repository:

```bash
cd ~/git/batmis/batmis_port/BatMis-3.00/src
make batmisAWFM_index_builder
make batmisAWFM
```

Optional compatibility target (same program code):

```bash
make batman_zero_mismatch_indexX
```

## 3) Build an AwFm index from FASTA

```bash
cd ~/git/batmis/batmis_port/BatMis-3.00/src
./batmisAWFM_index_builder \
  -f /path/to/reference.fa \
  -o /path/to/reference.awfmi
```

Useful tuning options:

- `-k <int>`: k-mer seed table length (default `5` for DNA)
- `-s <int>`: suffix array compression ratio (default `8`)
- `--keep-sa-in-memory`: keep compressed SA in memory for faster queries
- `--store-seq`: store original sequence in the index

Show builder help:

```bash
./batmisAWFM_index_builder --help
```

## 4) Run batmisAWFM

```bash
cd ~/git/batmis/batmis_port/BatMis-3.00/src
./batmisAWFM \
  -q /path/to/reads.fq \
  -x /path/to/reference.awfmi \
  -m 1 \
  -o hits.txt
```

Arguments:

- `-q <reads.fa/fq>`: input reads (FASTA or FASTQ)
- `-x <index.awfmi>`: AwFm index file
- `-m <maxhits>`: cap reported hit count per read (default `1`)
- `-o <hits.txt>`: optional output path (default stdout)
- `--quiet` / `-Q`: suppress per-read output
- `--timing`: print timing summary to stderr

Legacy alias accepted for index path:

- `-b <index.awfmi>` (equivalent to `-x`)

## 5) One-command helper script

You can build/reuse index and run search in one command:

```bash
cd ~/git/batmis/batmis_port/BatMis-3.00
./scripts/run_batmisAWFM.sh -f /path/to/reference.fa -q /path/to/reads.fq
```

Common variants:

```bash
# Force index rebuild and write hits to file
./scripts/run_batmisAWFM.sh \
  -f /path/to/reference.fa \
  -q /path/to/reads.fq \
  -x /path/to/reference.awfmi \
  -o hits.txt \
  -m 5 \
  --rebuild-index

# Pass index-builder tuning options
./scripts/run_batmisAWFM.sh \
  -f /path/to/reference.fa \
  -q /path/to/reads.fq \
  -k 6 \
  -s 12 \
  --keep-sa-in-memory
```

Script help:

```bash
./scripts/run_batmisAWFM.sh --help
```

## Output format

One line per read:

```text
query_id\tF\tSA_start\tSA_end\thit_count
```

## Quick smoke test

```bash
./batmisAWFM --help
./batmisAWFM_index_builder --help
```

If both commands print usage text, build and linkage are working.

## Basic correctness test

Run the included regression test that builds a tiny reference/index and validates
one known hit plus one known miss:

```bash
cd ~/git/batmis/batmis_port/BatMis-3.00
./scripts/test_batmisAWFM_basic.sh
```

Expected result:

```text
[test][PASS] batmisAWFM basic exact-match behavior is correct.
```
