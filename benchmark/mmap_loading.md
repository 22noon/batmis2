# Mapping the index files instead of reading them

## Summary

- batmis1 and batmis2 now map the SA samples (`.sa`) and the reference (`.pac`) with `mmap`,
  read-only and prefaulted, instead of reading them into fresh allocations. `--no-mmap`
  restores reading.
- **Index loading falls by 70–90%** (500 Mb genome, every 4th SA row: 0.50 → 0.13 s for
  batmis2; chr21, every row: 0.16 → 0.03 s), saving 0.1–0.5 s of every run's wall time:
  up to 46% of a short run and 20–32% of the 16-thread runs on the 500 Mb genome. On
  E. coli, whose files are small, there is nothing to save.
- **Search time, peak memory and output are unchanged.** `make test` passes, and mapped and
  read loading give identical SAM hits in 24 comparisons.
- It also removes the load cost of the denser SA sampling (`columba_comparison.md`): mapped,
  the 500 Mb index loads in 0.12–0.13 s at every 4th row and at every 8th.

## Why

`fread` copies every byte of a file from the page cache into memory the program has just
allocated, touching each page twice. A read-only private mapping shares the page cache's
pages instead. `MAP_POPULATE` prefaults the whole mapping at once, so that the search does not
take page faults later, and reads the file in order if it is not cached.

Measured on the 500 Mb genome's 476 MB `.sa` (`load_parts.cpp`, page cache warm, one core):

| | seconds |
|---|---:|
| `BWTLoad` reading the samples (`fread` into a new allocation) | 0.163 |
| mapped read-only, every page then touched | 0.003–0.009 |
| mapped read-only with `MAP_POPULATE` | 0.006 |
| mapped writable (`MAP_PRIVATE`) with `MAP_POPULATE` | 0.140 (it copies every page) |

So the mapping has to stay read-only. The library expects `saValue[0]` to be −1, and in the
file that word holds the text length, so the one page holding it is mapped again writable:
a private copy of 4 KB.

## Implementation (`src/batmis1.cpp`, `src/batmis2.cpp`, identical)

- `Map_File`: maps a file read-only with `MAP_POPULATE` over an anonymous region a little
  larger, so the bytes past its end read as zeros (`.pac` is used with 16 bytes of padding).
- `Map_Sa`: checks the `.sa` header against the index (inverse SA0, cumulative counts), as
  `BWTLoad` does; maps the file; remaps the first page writable and sets `saValue[0] = -1`;
  sets `saInterval` and the sample count; builds the boundary values as `BWTLoad` does.
- `Load_Index`: `BWTLoad` without the SA, then `Map_Sa`; with `--no-mmap`, `BWTLoad` with it,
  as before.
- `LoadPac`: maps `.pac`, or with `--no-mmap` reads it into a buffer. `Pac` is now a small
  read-only view (`data()`, `[]`) of either.
- The BWT and occurrence counts (`.bwt`, `.fmv`) are still read by `BWTLoad`: the library
  allocates them a few words longer than the file and writes into that tail, which a
  read-only mapping does not allow. They are about 0.05 s per 500 Mb index.

## Timing

Page cache warm, `--lookup 11`, median of 3; 1M reads of 100 bp; one core pinned, or 16
threads. Wall time and batmis's own index-load time, read (`--no-mmap`) → mapped:

| index (SA interval) | run | batmis1 wall | batmis2 wall | batmis2 index load |
|---|---|---|---|---|
| 500 Mb (4) | 1 read | 0.82 → 0.56 | 0.98 → **0.53** | 0.48 → 0.11 |
| 500 Mb (4) | K = 0 first, 1 core | 2.41 → 2.17 | 2.60 → 2.17 | 0.49 → 0.13 |
| 500 Mb (4) | K = 2 all, 1 core | 6.33 → 5.87 | 6.31 → 5.81 | 0.50 → 0.12 |
| 500 Mb (4) | K = 0 first, 16 threads | 1.22 → 0.98 | 1.40 → **0.95** | 0.50 → 0.13 |
| 500 Mb (4) | K = 2 all, 16 threads | 2.41 → 1.97 | 2.31 → 1.82 | 0.50 → 0.15 |
| 500 Mb (8) | 1 read | 0.74 → 0.55 | 0.77 → 0.52 | 0.32 → 0.11 |
| 500 Mb (8) | K = 0 first, 16 threads | 1.15 → 0.99 | 1.21 → 0.97 | 0.34 → 0.13 |
| chr21 (1) | 1 read | 0.43 → 0.35 | 0.39 → 0.24 | 0.15 → 0.02 |
| chr21 (1) | K = 2 all, 1 core | 5.03 → 4.72 | 4.91 → 4.77 | 0.16 → 0.03 |
| chr21 (1) | K = 0 first, 16 threads | 0.79 → 0.76 | 0.82 → 0.66 | 0.16 → 0.03 |
| E. coli (1) | any | within ±0.05 s | within ±0.03 s | 0.03 → 0.02 |

Search time is the same both ways (within 3%, the run-to-run noise), and so is peak memory:
mapped pages still count as resident, but they are the page cache's, shared by every process
on the machine that maps the same index, and the kernel can drop and reread them under
memory pressure instead of swapping.

## Caveats

- **Cold page cache:** `MAP_POPULATE` reads the file as `fread` would, so nothing is lost,
  but the gain is the copy only. Not measured (dropping the page cache needs root).
- **Rebuilding an index in place** while a mapper is using it can make the mapper crash
  (SIGBUS) where reading would have kept the old copy. Build into another directory, or run
  with `--no-mmap`.
- **Filesystems:** `mmap` works on local disks and on Lustre/NFS; where it does not, use
  `--no-mmap`.

## Files

- `src/batmis1.cpp`, `src/batmis2.cpp`: the loaders above, and `--no-mmap`.
- `load_parts.cpp`: the per-file load times above
  (`g++ -O2 -std=c++17 -fpermissive -I../index load_parts.cpp ../index/libfmindex_core.a`).
