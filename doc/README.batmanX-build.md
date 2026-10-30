# batmanX Build Instructions (SSE/MMX and scalar)

This document shows how to build `batmanX` in two modes:

- SSE/MMX enabled (default)
- SSE/MMX disabled (scalar)

For FM-index profile separation (legacy/full vs non-compressed FM-only), see `README.fmindex-profiles.md`.

The commands below were validated in this repository.

## Prerequisites

- `autotools` installed (`autoreconf` available)
- `gcc` and `g++`
- `make`
- `objdump` (optional, for binary-level verification)

## 1) SSE/MMX-enabled build (default)

From repository root:

```bash
cd /home/chandana.tennakoon@pirbright.ac.uk/git/batmis/batmis_port/BatMis-3.00
autoreconf -fi
./configure
```

Build `batmanX` directly:

```bash
OPT_CFLAGS="-w -O3 -funroll-loops -maccumulate-outgoing-args -msse2"
OPT_CXXFLAGS="-w -O3 -funroll-loops -maccumulate-outgoing-args -msse2"

rm -f src/*.o batmanX && \
gcc $OPT_CFLAGS -I. -Isrc -c src/BWT.c -o src/BWT.o && \
gcc $OPT_CFLAGS -I. -Isrc -c src/BWTConstruct.c -o src/BWTConstruct.o && \
gcc $OPT_CFLAGS -I. -Isrc -c src/MiscUtilities.c -o src/MiscUtilities.o && \
gcc $OPT_CFLAGS -I. -Isrc -c src/MemManager.c -o src/MemManager.o && \
gcc $OPT_CFLAGS -I. -Isrc -c src/TextConverter.c -o src/TextConverter.o && \
gcc $OPT_CFLAGS -I. -Isrc -c src/r250.c -o src/r250.o && \
gcc $OPT_CFLAGS -I. -Isrc -c src/QSufSort.c -o src/QSufSort.o && \
gcc $OPT_CFLAGS -I. -Isrc -c src/iniparser.c -o src/iniparser.o && \
gcc $OPT_CFLAGS -I. -Isrc -c src/inistrlib.c -o src/inistrlib.o && \
gcc $OPT_CFLAGS -I. -Isrc -c src/dictionary.c -o src/dictionary.o && \
gcc $OPT_CFLAGS -I. -Isrc -c src/DNACount.c -o src/DNACount.o && \
gcc $OPT_CFLAGS -I. -Isrc -c src/Timing.c -o src/Timing.o && \
gcc $OPT_CFLAGS -I. -Isrc -c src/Socket.c -o src/Socket.o && \
gcc $OPT_CFLAGS -I. -Isrc -c src/HSP.c -o src/HSP.o && \
gcc $OPT_CFLAGS -I. -Isrc -c src/HSPstatistic.c -o src/HSPstatistic.o && \
gcc $OPT_CFLAGS -I. -Isrc -c src/karlin.c -o src/karlin.o && \
g++ $OPT_CXXFLAGS -std=c++17 -fpermissive -I. -Isrc -c src/batmanX.cpp -o src/batmanX.o && \
g++ $OPT_CXXFLAGS -std=c++17 -fpermissive -o batmanX src/batmanX.o src/BWT.o src/BWTConstruct.o src/MiscUtilities.o src/MemManager.o src/TextConverter.o src/r250.o src/QSufSort.o src/iniparser.o src/inistrlib.o src/dictionary.o src/DNACount.o src/Timing.o src/Socket.o src/HSP.o src/HSPstatistic.o src/karlin.o -lm
```

Verify mode:

```bash
./batmanX --build-info
```

Expected output:

```text
batmanX SIMD mode: MMX/SSE2
```

Optional checks:

```bash
grep -E "^#define MMX" config.h
grep -E "^SIMD_CFLAGS" src/Makefile
objdump -d ./batmanX | grep -Ec "\<(movdqa|movdqu|pand|psrld|paddd|pcmpeqd)\>"
```

Expected indicators:

- `#define MMX 1` exists in `config.h`
- `SIMD_CFLAGS = -msse2` in `src/Makefile`
- instruction count is typically high (non-zero, often hundreds)

## 2) SSE/MMX-disabled build (scalar)

From repository root:

```bash
cd /home/chandana.tennakoon@pirbright.ac.uk/git/batmis/batmis_port/BatMis-3.00
./configure --disable-mmx
```

Rebuild `batmanX` with the same optimization profile, but without `-msse2`:

```bash
OPT_CFLAGS="-w -O3 -funroll-loops -maccumulate-outgoing-args"
OPT_CXXFLAGS="-w -O3 -funroll-loops -maccumulate-outgoing-args"

rm -f src/*.o batmanX && \
gcc $OPT_CFLAGS -I. -Isrc -c src/BWT.c -o src/BWT.o && \
gcc $OPT_CFLAGS -I. -Isrc -c src/BWTConstruct.c -o src/BWTConstruct.o && \
gcc $OPT_CFLAGS -I. -Isrc -c src/MiscUtilities.c -o src/MiscUtilities.o && \
gcc $OPT_CFLAGS -I. -Isrc -c src/MemManager.c -o src/MemManager.o && \
gcc $OPT_CFLAGS -I. -Isrc -c src/TextConverter.c -o src/TextConverter.o && \
gcc $OPT_CFLAGS -I. -Isrc -c src/r250.c -o src/r250.o && \
gcc $OPT_CFLAGS -I. -Isrc -c src/QSufSort.c -o src/QSufSort.o && \
gcc $OPT_CFLAGS -I. -Isrc -c src/iniparser.c -o src/iniparser.o && \
gcc $OPT_CFLAGS -I. -Isrc -c src/inistrlib.c -o src/inistrlib.o && \
gcc $OPT_CFLAGS -I. -Isrc -c src/dictionary.c -o src/dictionary.o && \
gcc $OPT_CFLAGS -I. -Isrc -c src/DNACount.c -o src/DNACount.o && \
gcc $OPT_CFLAGS -I. -Isrc -c src/Timing.c -o src/Timing.o && \
gcc $OPT_CFLAGS -I. -Isrc -c src/Socket.c -o src/Socket.o && \
gcc $OPT_CFLAGS -I. -Isrc -c src/HSP.c -o src/HSP.o && \
gcc $OPT_CFLAGS -I. -Isrc -c src/HSPstatistic.c -o src/HSPstatistic.o && \
gcc $OPT_CFLAGS -I. -Isrc -c src/karlin.c -o src/karlin.o && \
g++ $OPT_CXXFLAGS -std=c++17 -fpermissive -I. -Isrc -c src/batmanX.cpp -o src/batmanX.o && \
g++ $OPT_CXXFLAGS -std=c++17 -fpermissive -o batmanX src/batmanX.o src/BWT.o src/BWTConstruct.o src/MiscUtilities.o src/MemManager.o src/TextConverter.o src/r250.o src/QSufSort.o src/iniparser.o src/inistrlib.o src/dictionary.o src/DNACount.o src/Timing.o src/Socket.o src/HSP.o src/HSPstatistic.o src/karlin.o -lm
```

Verify mode:

```bash
./batmanX --build-info
```

Expected output:

```text
batmanX SIMD mode: scalar
```

Optional checks:

```bash
grep -E "^#define MMX" config.h || true
grep -E "^SIMD_CFLAGS" src/Makefile
objdump -d ./batmanX | grep -Ec "\<(movdqa|movdqu|pand|psrld|paddd|pcmpeqd)\>"
```

Expected indicators:

- no `#define MMX` in `config.h`
- `SIMD_CFLAGS` is empty in `src/Makefile`
- SSE instruction count drops significantly

## Notes

- `make all` may fail for other targets in this codebase due to unrelated link-order issues in legacy rules. The direct `batmanX` compile/link command above is the reliable path for `batmanX` itself.
- Warnings in legacy C files are expected and pre-existing.

## Minimal Legacy FM-index Build (no HSP/Karlin/Socket/Timing objects)

To build a minimal legacy FM-index archive that excludes non-core modules, run:

```bash
cd /home/chandana.tennakoon@pirbright.ac.uk/git/batmis/batmis_port/BatMis-3.00/src
make clean
make libfmindex_core.a
```

This target builds `BWT.c` in `FMINDEX_ONLY` mode and archives only FM-index core objects:

- `BWT_fmindex.o`
- `BWTConstruct.o`
- `MiscUtilities.o`
- `MemManager.o`
- `TextConverter.o`
- `r250.o`
- `QSufSort.o`
- `DNACount.o`

Excluded from this minimal archive:

- `HSP.o`
- `HSPstatistic.o`
- `karlin.o`
- `Socket.o`
- `Timing.o`

Behavior note:

- Approximate-match scoring paths are intentionally unavailable in `FMINDEX_ONLY` builds. Exact FM-index routines are retained.

### Minimal target smoke test

Build the dedicated minimal executable and run a quick sanity check:

```bash
cd /home/chandana.tennakoon@pirbright.ac.uk/git/batmis/batmis_port/BatMis-3.00/src
make clean
make batman_zero_min_fmindex
./batman_zero_min_fmindex --help | head -n 8
```

Verify that the minimal archive does not carry non-core module references:

```bash
ar t libfmindex_core.a
nm -u BWT_fmindex.o | grep -E "HSP|Socket|karlin|stat_" || true
```

Expected result:

- `batman_zero_min_fmindex` builds successfully.
- `libfmindex_core.a` contains only FM-index core objects.
- No unresolved references to `HSP`, `Socket`, `karlin`, or `stat_` symbols from `BWT_fmindex.o`.

## FM-index Function Purpose Map (legacy BWT path)

The legacy FM-index flow is split between index construction/save and index load/query.

Construction and save path:

- `BWTGenerateOccValueFromBwt()` in `src/BWTConstruct.c`
	- Derives `occValue` and `occValueMajor` tables from BWT code.
- `BWTSaveBwtCodeAndOcc()` in `src/BWTConstruct.c`
	- Writes BWT code and occurrence tables to disk.
- `BWTSaveSaValue()` in `src/BWTConstruct.c`
	- Writes sampled suffix-array values.
- `BWTSaveInverseSa()` in `src/BWTConstruct.c`
	- Writes inverse-SA samples.

Load and query path:

- `BWTLoad()` in `src/BWT.c`
	- Loads BWT/Occ/SA files into runtime structures.
- `BWTDecode()`, `BWTDecodeAll()` in `src/BWT.c`
	- Decodes packed BWT blocks and character counts.
- `BWTOccValue()`, `BWTAllOccValue()` in `src/BWT.c`
	- Computes rank/occurrence values for backward search.
- `BWTPsiMinusValue()`, `BWTSaValue()`, `BWTInverseSa()` in `src/BWT.c`
	- Core LF/SA/inverse-SA accessors used by FM-index traversal.
- `BWTBackwardSearch()` and related helpers in `src/BWT.c`
	- Performs exact FM-index backward search.

## FM-index Memory Allocation Map (old FM-index)

Primary allocation sites are in `BWTCreate()` and `BWTLoad()` in `src/BWT.c`.

- `bwtCode`
	- Size basis: `BWTResidentSizeInWord(textLength)` (+ alignment padding in load path).
	- Role: packed BWT symbols.
- `occValue`
	- Size basis: minor occurrence table size derived from text length and `OCC_INTERVAL`.
	- Role: rank acceleration table.
- `occValueMajor`
	- Size basis: major occurrence table with interval `OCC_INTERVAL_MAJOR`.
	- Role: coarse rank checkpoints.
- `saValue`
	- Size basis: sampled SA array (`textLength / saInterval` style layout).
	- Role: resolve text positions after FM steps.
- `inverseSa`
	- Size basis: sampled inverse-SA (`textLength / inverseSaInterval` layout).
	- Role: faster SA lookups for some query paths.
- `cachedSaIndex`
	- Size basis: either `ALPHABET_SIZE + 1` or loaded cache table size.
	- Role: cached SA range boundaries for short prefixes.
- `decodeTable`
	- Size basis: `DNA_OCC_CNT_TABLE_SIZE_IN_WORD`.
	- Role: fast packed-word occurrence decoding.

## Where adler32/crc32/r250/inf*/trees are used

Build placement in `src/Makefile.am`:

- `r250.c`
	- In FM core object list (`FMINDEX_CORE_OBJS`), required by legacy FM-index path.
- `adler32.c`, `crc32.c`, `infback.c`, `inffast.c`, `inflate.c`, `inftrees.c`, `trees.c`
	- Built into `libz.a` (zlib bundle), not part of `FMINDEX_CORE_OBJS`.

Runtime usage classification:

- FM-index core (`BWT.c`/`BWTConstruct.c`) reads and writes raw binary index artifacts using stdio (`fopen64`, `fread`, `fwrite`).
- zlib files support legacy gzip/BAT result I/O paths (for example in `batdecode.cpp` and legacy `batman.cpp`), not FM-index construction/loading.

## Can zlib-era files be removed for non-compressed FM-index?

Yes, for FM-index-only builds.

- Safe for non-compressed FM-index target:
	- Keep `r250.c`.
	- Exclude zlib objects by building `batman_zero_min_fmindex` / `libfmindex_core.a` only.
- Not safe if you still need legacy compressed BAT result tooling:
	- Keep `libz.a` and its object set.

Practical recommendation:

- Use dedicated FM-only target for non-compressed FM-index workflows.
- Keep zlib path isolated for legacy BAT compatibility, rather than mixing both dependency surfaces.

## Build profile commands

From `src/`:

```bash
make profiles-help
make legacy-full
make fmindex-only-lib
make fmindex-only
```

Profile intent:

- `legacy-full`: builds legacy utility set and zlib static library (`reverse`, `filter`, `libz.a`).
- `fmindex-only-lib`: builds only `libfmindex_core.a`.
- `fmindex-only`: builds `batman_zero_min_fmindex` against FM core only.

## FM-index symbol index (line anchors)

Core load/query symbols in `src/BWT.c`:

- `BWTCreate` at line 67
- `BWTLoad` at line 95
- `BWTDecode` at line 323
- `BWTDecodeAll` at line 442
- `BWTOccValue` at line 547
- `BWTAllOccValue` at line 635
- `BWTResidentSizeInWord` at line 845
- `BWTFileSizeInWord` at line 851
- `BWTPsiMinusValue` at line 885
- `BWTSaValue` at line 923
- `BWTInverseSa` at line 949
- `BWTBackwardSearch` at line 1645

Construction/save symbols in `src/BWTConstruct.c`:

- `BWTGenerateOccValueFromBwt` at line 802
- `BWTSaveBwtCodeAndOcc` at line 994
- `BWTSaveSaValue` at line 1126
- `BWTSaveInverseSa` at line 1195
