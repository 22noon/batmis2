# Threads (`-t N`)

## Summary

`batman_zero_min -t N` searches with N threads (default 1; `-t 0` uses all CPUs). The
output is byte-identical for any N, any `--batch` and the build before threads.

- **End to end** (2M 100 bp reads, 500 Mb genome, index load excluded, mean of 2 runs):

  | workload | before threads | `-t 6` (1 per P-core) | `-t 22` (all CPUs) | speed-up at 22 |
  |---|---:|---:|---:|---:|
  | `-n 0 -m 0` | 14.69 s | 2.98 s | 1.70 s | 8.6× |
  | `-n 2 -m 0` | 52.29 s | 10.84 s | 5.37 s | 9.7× |
  | `-n 2 -m 3 --sam` | 76.99 s | 18.36 s | 12.18 s | 6.3× |
  | `--no-ref -n 1 -m 0` | 67.94 s | 14.32 s | 7.81 s | 8.7× |

- **`-t 1` costs nothing.** It is within ±2% of the build before threads, the same as
  run-to-run noise.
- **Threads and lanes compete for the same thing**, the cache misses a core can keep in
  flight. So with threads the best `--batch` falls from 4–16 to about 2–4, and two
  hyperthreads on one core gain little from lanes.
- **Input is read ahead.** A reader thread reads block k+1 while block k is searched
  (double buffering), and the workers encode the reads. At `-t 22` this cut `-n 0`
  from 3.44 s to 1.70 s.
- **SAM records are built on the workers** (`Sam_Encoder`, below). The main thread only
  writes their bytes. At `-t 22`, SAM output fell from 6.5 s to 1.1 s.
- **What limits it now:**
  - **BAM:** compressing about 4–5 s per 2M reads with htslib's zlib.
  - **Gzip input:** inflating about 3 s per 2M reads.
  - **The reader thread itself:** about 1.4 s per 2M reads.

Machine: Intel Core Ultra 7 165H. It has 6 P-cores with 2 hyperthreads each (CPUs 0–11),
8 E-cores (12–19) and 2 low-power E-cores (20–21). Other users kept the load average at
5–12 during these runs, so differences under about 10% are noise.

## How it works

`Run_Batched` (`src/batman_zero_min.cpp`):

1. **Blocks.** A block holds 2 × N chunks of 4,096 reads, or a single chunk of 4,096 with
   one thread, as before.
2. **Reader** (`Block_Reader`). With N > 1, a thread of its own reads the next block's
   records into a second buffer while the current block is processed. It only splits
   records, then hands the block over at the boundary.
3. **Encode.** Each thread checks a contiguous part of the block: a read is searched if it
   has the first read's length and only ACGT. Then a prefix sum gives each searched read
   its slot, and each thread encodes its part's reads and their reverse complements
   through a 256-entry table (`Base_Codes`). The table replaced the per-base
   locale-aware `toupper` in `EncodeBase`. That alone took 1 s off the 2M-read input at
   `-t 1`.
4. **Search.** The threads take 4,096-read chunks from an atomic counter and run
   `Batch_Search` on `[First, Last)`. They write into per-read slots that nothing else
   touches, so the search needs no locks.
5. **Output.**
   - Each thread formats its part's text lines. For SAM/BAM it builds its part's
     candidates, decodes their SA rows (`Sa_Decode` per part) and encodes the records
     (`Sam_Encoder`).
   - The main thread then writes the parts' bytes in order.
   - BAM compression uses htslib's thread pool (`hts_set_threads`).

The worker threads are a persistent pool (`Thread_Pool`). It wakes them once per step of
each block.

**Per-thread state:**
- the lanes (`static thread_local` in `Batch_Search`);
- `Branch_Detect`'s branch arrays;
- the scratch vectors of `Sa_Decode`, `Chain_Requests`, `Sam_Requests` and
  `Encode_Sam_Read`, and each part's `Sam_Encoder`;
- the `BATMAN_BATCH_STATS` counters, which are added to the main thread's when a worker
  exits.

The program table is built once, thread-safely (a function-local static).

**`--timing`:**
- `input_parse` is the time spent waiting for the reader, plus encoding.
- `input_read` is the reader's own time. It overlaps the rest when N > 1, and is part of
  `input_parse` when N = 1.

## Does threading change batching?

**Lanes and threads both hide memory latency, and they share the hardware that does it.**
Each core has a fixed number of line-fill buffers (cache misses in flight), and its two
hyperthreads share them. The memory system's latency also rises once every core is
issuing misses.

The 100k-read sweep below gives search time at the best B as a fraction of B = 1. Lower
means lanes help more.

| | `--no-ref -n 1` | `-n 0` | SAM `-n 2` | `-n 2` | `-n 4` |
|---|---:|---:|---:|---:|---:|
| 1 thread | 0.52 | 0.64 | 0.68 | 0.85 | ~1.0 |
| 2 threads, separate cores | 0.51 | 0.70 | 0.71 | 0.89 | 1.0 (B = 1) |
| 2 threads, same core (HT) | 0.86 | 0.89 | 0.95 | 1.0 (B = 1) | 1.0 (B = 1) |
| 6 threads, 1 per P-core | 0.57 | 0.61 | 0.77 | 0.87 | 1.0 (B = 1) |
| 12 threads (P-cores, HT) | 0.78 | 0.87 | 0.96 | 0.89 | 1.0 (B = 1) |
| 8 threads (E-cores) | 0.82 | 0.84 | 0.89 | 1.0 (B = 1) | 1.0 (B = 1) |
| 22 threads | 0.75 | 0.84 | 0.91 | 0.85 | 0.97 |

- **One thread per core keeps the lane gain.** Each thread has its core's fill buffers
  to itself.
- **Hyperthreads take most of it away.** Two threads on a core already overlap each
  other's misses, which is what lanes do.
- **Under full load, B = 2–4 is best,** and the spread between B values is small: B ≥ 2
  is within about 10% of the best.
- **At `-n 4`, B = 1 stays best.** Its branching extensions already keep misses in
  flight.

**Search speed-up over 1 thread** (100k reads, best B): 5.1–6.7× with 6 threads (one per
P-core), 6.1–9.5× with 12, and 6.6–15.6× with 22.

### `--batch auto` with threads

Each thread has its own `Batch_Chooser`, which times the candidates on its own chunks
while the other threads run. So each choice reflects the load and the core type that
thread runs on.

On 2M reads (2 runs each), `auto` came within 1.00–1.12× of the best fixed B in all
12 combinations of `-n 0`, `-n 2`, `--no-ref -n 1` and `-n 4` with 6, 12 and 22 threads.

| | `auto` / B = 1 |
|---|---:|
| `--no-ref -n 1`, `-t 6` | 0.67 |
| `-n 0`, `-t 6` | 0.73 |
| `-n 4`, `-t 12` | 1.07 (B = 1 best) |

**Limits:**
- **The threads' choices scatter**, e.g. `batch=auto:[4,2,4,4,2,2]`. That is mostly
  noise between options that are within a few percent of each other.
- **The trial grows with N.** It is 15 chunks per thread, about 1.35M reads at
  `-t 22`, and smaller inputs end during it. That costs little, since under load the
  candidates are close. One trial shared by all threads would make it independent of N.

## Input: double buffering and encoding on the workers

Before this change, the main thread read and encoded each block while the workers
waited. At `-t 22`, that serial input step (2.3–2.5 s per 2M reads) was longer than the
search itself (0.9 s at `-n 0`).

| 2M reads, `-t 22` | input wait | reader (overlapped) | search | total − load |
|---|---:|---:|---:|---:|
| `-n 0`, threads only | 2.31 s | – | 0.91 s | 3.44 s |
| `-n 0`, + reader and encoding | 0.26 s | 1.42 s | 1.02 s | 1.70 s |
| `-n 2`, threads only | 2.36 s | – | 4.09 s | 6.92 s |
| `-n 2`, + reader and encoding | 0.14 s | 1.66 s | 4.66 s | 5.37 s |

- **Memory:** only the input side is doubled (record IDs, sequences, qualities), about
  27 MB at `-t 22`. The search block's result lists are not doubled, because block k's
  search ends before block k+1's begins.
- **Why search time rises a little:** the reader thread shares the cores with the
  workers.

## Gzip input

`-q` takes plain or gzip-compressed FASTA/FASTQ, or `-` for stdin. It also reads BGZF
files and multi-member files.
- **How:** `Read_File` reads through zlib (`gzread`), which passes plain files through
  unchanged. It splits lines from 1 MB buffers of its own. That replaced
  `std::ifstream` + `std::getline` and the FASTA `seekg` look-ahead.
- **Inflate thread:** with N > 1, a thread of its own inflates into a ring of four
  buffers while the reader thread splits the lines of the one before. With `-t 1`,
  `Get_Line` reads each buffer itself, as before.
- **Errors:** truncated or corrupt compressed data stops the run with an error and exit
  code 1.
- **Plain input got a little faster:** the reader's time fell from 0.95 s to 0.80 s per
  2M reads.

2M 100 bp reads, `-n 0 -m 0`, mean of 2 runs (the 443 MB FASTQ is 91 MB gzipped):

| | reader | input wait | search | total |
|---|---:|---:|---:|---:|
| `-t 1`, plain | 0.80 s | 1.41 s | 8.79 s | 11.23 s |
| `-t 1`, gzip | 3.80 s | 4.40 s | 8.84 s | 14.17 s |
| `-t 22`, plain | 0.84 s | 0.15 s | 0.77 s | 1.76 s |
| `-t 22`, gzip | 3.62 s | 2.79 s | 0.72 s | 4.28 s |

- **Inflating costs about 3 s per 2M reads**, on the reader thread.
- **It is hidden when the search takes longer than that,** e.g. `-n 2` at `-t 22`
  (5 s).
- **For fast workloads at high `-t`, the reader becomes the limit.**

**With the inflate thread** (2M reads, `-n 0`, `-t 22`, mean of 2 runs, load 6–24):

| | reader | total − index load |
|---|---:|---:|
| reader thread inflates | 4.99 s | 5.29 s |
| inflate thread | 4.19 s | 4.57 s |

- **Gain:** the line splitting (about 0.8 s) now overlaps inflating. 14% of the run, as
  expected.
- **No faster plain-gzip path without a new library.** zlib's inflate is serial, and
  libdeflate cannot stream a single-member file. ISA-L or zlib-ng would inflate 2–3×
  faster, but each is a new dependency.

## SAM/BAM output: records built on the workers

**The two writers** (`src/sam_writer.hpp`):
- **`Sam_Encoder`** (one per thread) builds records into bytes. For SAM it uses
  `sam_format1` plus a newline. For BAM it lays them out as `bam_write1` does (block_len,
  the 32-byte core, the data).
- **`Sam_Writer::Write`** writes an encoder's bytes:
  - SAM: one `hwrite`.
  - BAM: per record, `bgzf_flush_try` + `bgzf_write`, as `bam_write1`, so the BGZF
    blocks are the same.

`Encode_Sam_Read` (formerly `Write_Sam_Read`) runs in each part's job, right after
`Sa_Decode`. The main thread only writes the parts in order. The header is only read
(`sam_hdr_tid2name`).

**Output is byte-identical.** 228 runs of raw SAM and BAM bytes, including `@PG`, match
the previous build: `-t` 1, 4, 22; `-n 0..5`; `-m 0/1/5`; `--no-ref`; `--forward-only`;
`--quiet`; plain and gzip input; tiny and empty inputs.

2M reads, `-n 2 -m 3`, `-t 22`, mean of 2 runs (load 6–25, so ±15%):

| | output step, before | output step, after | total − load, before | total − load, after |
|---|---:|---:|---:|---:|
| SAM | 6.52 s | 1.12 s | 18.35 s | 13.31 s |
| BAM | 6.10 s | 4.46 s | 16.92 s | 17.32 s |

- **SAM:** building and formatting the records takes 1.2 s across the workers. Writing
  the 584 MB on the main thread takes 1.0 s. Both were measured in an instrumented build.
- **BAM is limited by compression,** not by the hand-over. The writes take about 4.9 s:
  the main thread waits on htslib's compression pool for 488 MB of BAM data.
  - System `samtools` 1.19, which uses libdeflate, needs 23 CPU-seconds for the same data
    (2.2 s on 22 threads).
  - The bundled htslib compresses with plain zlib, which is slower.
  - Faster BAM output needs libdeflate in htslib, which is a new dependency, or a lower
    compression level.
- **`-t 1`:** unchanged (103.0 s against 102.4 s for SAM, within noise).

## Verification

- **Output byte-identical** to the build before threads in 400 runs (`@PG` excluded):
  - `-t` 1, 2, 4, 8, 22 × `--batch` auto, 1, 16;
  - `-n 0..5`, text, SAM and BAM, `--forward-only`, with and without `--no-ref`;
  - 100k- and 20k-read sets, a 5k-read 32 bp set, a 100-read input and an empty one.
- **ThreadSanitizer:** no reports on 20k, 1k, 100 and 0 reads. That covers `-n 2`,
  `-n 5`, SAM, BAM and `--no-ref`, at `-t 4` with the default chunks. It also ran with
  16-read chunks (`-DBATMAN_BATCH_BLOCK=16`), so that each run spans many blocks and the
  reader's hand-over is exercised, including inputs that end exactly on a block
  boundary.
- **Test suites:** `make test-located`, `make test-sam` and `make test-structured` pass.
  They now also compare `-t 1`, `-t 4`, `-t 3 --batch 4`, `-t 3 --no-ref`,
  `-t 2 --forward-only` and threaded BAM output. These run on
  `batman_zero_min_chunk16`, a build with 16-read chunks, so that the tests' few hundred
  reads are split across blocks, chunks and threads. `test-located` also checks gzip input.
- **Gzip:**
  - Identical output from plain, gzip, BGZF and two-member gzip FASTQ, from FASTA and
    multi-line FASTA (plain and gzip), and from gzip on stdin. Checked for text and SAM,
    at `-t 1` and `-t 4`.
  - Truncated, corrupt and missing files fail with exit code 1.

## Next

- **BAM compression** (about 4–5 s per 2M reads at `-t 22`): build htslib with libdeflate,
  or offer a lower compression level. A writer thread would only overlap the writes with
  the next block's search; the cores are already busy compressing.
- **The reader** (1.4–1.7 s per 2M reads) is the next limit for `-n 0`. FASTQ's four-line
  records make it possible to split the reading across threads.
- **A shared `--batch auto` trial**, so that its length does not grow with N.
- **Gzip input:**
  - BGZF input could be inflated in parallel (htslib's `bgzf_mt`).
  - Plain gzip can only be inflated serially, though ISA-L or zlib-ng would do it 2–3×
    faster. Both are new dependencies.
