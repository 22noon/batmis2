# BatMis FM index builder

`bwt_build_index` builds the forward and reverse FM indexes that BatMis
searches from a FASTA reference in a single step. It replaces the old
`scripts/build_index` pipeline (`filter`, `reverse`, `bwtformatdb` run twice,
and `ann2loc.pl`), and produces byte-identical index files.

This directory is self-contained: it does not use any file from `../src`.

## Building

```sh
make            # libfmindex_core.a and bwt_build_index
make clean
```

`libfmindex_core.a` is the FM index library (BWT construction, occurrence
counting, SA sampling, memory manager). `../src/Makefile` links
`batman_zero_min` against it and calls `make -C ../index libfmindex_core.a`
itself, so building `src/` also builds the library.

Override the compilers or flags in the usual way, e.g.
`make CC=clang CXX=clang++ CFLAGS="-O2"`.

`make -C ../src install` copies `bwt_build_index` to `../bin` together with the
searchers, if it has been built.

## Usage

```
bwt_build_index [options] <genome.fa> [outdir]

  -f, --forward-only   build only the forward index
  -r, --reverse-only   build only the reverse index
  -o, --outdir DIR     write index files to DIR
  -p, --progress       show BWT construction progress
  -b, --sa-budget SIZE memory for the SA samples of both indexes (default 1G;
                       K, M, G suffixes): every row sampled while it fits, then every
                       2nd or 4th, every 8th above (8n/s bytes for n bases)
  -s, --sa-interval N  sample every N-th SA row (1, 2, 4 or 8), overriding the budget
  -h, --help           show this help
```

**SA sampling.** The `.sa` files hold the suffix-array value of every s-th row. Denser
sampling makes locating a row cheaper (about s - 1 LF steps on average, each a cache miss on
a large genome) and lets the mappers finish single-row search paths sooner, at 8n/s bytes for
both indexes. The builder picks the densest s of 1, 2, 4, 8 within the budget (1 GB by
default) and never goes sparser than 8:

| bases | interval | samples, both indexes |
|---|---:|---:|
| up to 134 M (viruses, bacteria, one human chromosome) | 1 | up to 1 GB |
| up to 268 M | 2 | up to 1 GB |
| up to 537 M | 4 | up to 1 GB |
| more (human: 3.1 G) | 8 | 8n/8 = 3.1 GB for human |

Both indexes always get the same interval, which is stored in each `.sa`; batmis1 and
batmis2 read it from there. `--sa-interval 8` reproduces the earlier files byte for byte.

By default both the forward and the reverse index are built. The output
directory can be given with `-o` or as the second argument, but not both; it
is created if it does not exist. Without one, files are written next to the
FASTA file.

Examples:

```sh
bwt_build_index hg38.fa                 # writes hg38.fa.* and revhg38.fa.* beside hg38.fa
bwt_build_index -o idx hg38.fa          # writes idx/hg38.fa.*, idx/revhg38.fa.*
bwt_build_index -r hg38.fa idx          # reverse index only (plus the shared files)
```

The searchers take the genome prefix `<outdir>/<genome.fa>`, for example
`batman_zero_min -g idx/hg38.fa -q reads.fa`, and find the `rev` files from it.

## Output files

For an input `<genome.fa>`, with `G = <outdir>/<genome.fa>` and
`R = <outdir>/rev<genome.fa>`:

| File | Built | Contents |
|---|---|---|
| `G.N` | always | Sequence headers, each followed by the `start:end` runs of non-ACGT characters (0-based offsets within that sequence) |
| `G.bin` | always | The concatenated sequences as plain text, lower-case `a/c/g/t`, no headers or newlines |
| `G.ann.location` | always | `0`, then one name line and one length line per sequence (name = header up to the first whitespace) |
| `G.pac` | always | The forward text packed 2 bits per base (see below) |
| `G.bwt`, `G.fmv`, `G.sa` | forward | Forward FM index: BWT, occurrence values, SA sampled every s rows (s by the memory budget, above) |
| `R.bwt`, `R.fmv`, `R.sa` | reverse | The same for the reversed text |

With `--reverse-only`, the forward index files are not built, but `G.pac` and
the other shared files are still written.

The reverse text's `.pac` is only needed while building and is deleted
afterwards.

### Text conversion

The reference is converted as the old `filter` did:
- Upper and lower case `A/C/G/T` become `a/c/g/t`.
- Every other character (N, IUPAC codes and so on) becomes `c`, and its
  position is recorded in `G.N`.
- Sequences are concatenated with no separator. Use `G.ann.location` to map
  text offsets back to sequence names.

### `.pac` format

Bases are coded `a=0, c=1, g=2, t=3` and packed four to a byte, high bits
first. The file ends with one byte holding the number of bases in the last
data byte. When the length is a multiple of 4, an extra zero byte comes before
it. This is the layout written by `SaveBytePacked()` in `TextConverter.c` and
read back by `BWTIncConstructFromPacked()`.

## Limits and small references

- The total length after filtering must be below 2^32 bases.
- Very small references are handled automatically:
  - **Under about 1.5 kb:** the library's workspace estimate is too small,
    and it fails with `BWTIncCreate() : targetNBit is too low!`. The builder
    raises `targetNBit` for these inputs and prints
    `Small reference: targetNBit raised to X`. The index is unchanged.
  - **Under 16 bases:** the BWT library crashes or hangs on texts this short.
    The builder appends a dummy sequence named `__batmis_pad__`, made only of
    `a`, to bring the text up to 16 bases. It is added last, so the
    coordinates of the real sequences don't change, and it appears in
    `G.N` and `G.ann.location` like any other sequence. **Ignore hits on
    `__batmis_pad__` when processing results.**

## Memory and I/O

- The forward text is packed in memory while the FASTA is read, so the input
  is read only once.
- The reverse text is packed by memory-mapping `G.bin` and reading it
  backwards.
  - If `mmap` fails, it falls back to reading 16 MB blocks from the end of
    the file.
  - Setting the environment variable `BWT_BUILD_INDEX_NO_MMAP` forces the
    fallback.
- The two indexes are built one after the other, so peak memory is that of a
  single index build.

## Build parameters

These are the `bwtformatdb` defaults; all but the SA sampling are fixed in the source:

| Parameter | Value |
|---|---|
| Memory pool size | 2097152 |
| `targetNBit` | 2.5 (raised for small references) |
| Initial / incremental build size | 10000000 / 10000000 |
| SA sampling interval | 1, 2, 4 or 8, by the memory budget (`--sa-budget`, `--sa-interval`) |

The `COLOR`, `CTOT`, `GTOA` and random-base options of the old tools are not
supported.
