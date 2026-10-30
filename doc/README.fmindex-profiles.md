# FM-index Build Profiles

This guide explains the build profiles that separate non-compressed FM-index artifacts from legacy zlib-based tooling.

## Scope

- Non-compressed FM-index profile:
  - Builds FM-index core library and minimal executable.
  - Does not require zlib object linkage for FM-index artifacts.
- Legacy/full profile:
  - Builds legacy helper utilities and zlib static library.

## Profile Targets

Run from `src/`:

```bash
make profiles-help
make legacy-full
make fmindex-only-lib
make fmindex-only
```

### Target behavior

- `make legacy-full`
  - Builds: `reverse`, `filter`, `libz.a`
  - Intended for legacy utility path and compatibility
- `make fmindex-only-lib`
  - Builds: `libfmindex_core.a`
  - FM-index core static archive only
- `make fmindex-only`
  - Builds: `batman_zero_min_fmindex`
  - Minimal non-compressed FM-index executable

## FM-index Core Contents

`libfmindex_core.a` contains:

- `BWT_fmindex.o`
- `BWTConstruct.o`
- `MiscUtilities.o`
- `MemManager.o`
- `TextConverter.o`
- `r250.o`
- `QSufSort.o`
- `DNACount.o`

Excluded from FM-index core archive:

- `HSP.o`
- `HSPstatistic.o`
- `karlin.o`
- `Socket.o`
- `Timing.o`

## Dependency Boundary

- FM-index artifact load/save path uses raw binary I/O in BWT modules.
- zlib object set (for example adler32/crc32/inflate/inftrees/trees) is kept in `libz.a` for legacy compressed result tooling compatibility.
- `r250.c` remains part of FM-index core.

## Quick Validation

Run from `src/`:

```bash
make clean
make fmindex-only-lib
make fmindex-only
ar t libfmindex_core.a
nm -u BWT_fmindex.o | grep -E "HSP|Socket|karlin|stat_" || true
```

Expected:

- FM-index targets build successfully.
- `libfmindex_core.a` shows FM-index core objects listed above.
- No unresolved `HSP`, `Socket`, `karlin`, or `stat_` references from `BWT_fmindex.o`.

## Notes

- Environment in this workspace does not provide `autoreconf`; keep manual edits synchronized across:
  - `src/Makefile.am`
  - `src/Makefile.in`
  - `src/Makefile`
- If you remove zlib sources globally, legacy compressed result tooling may break. Keep profile separation unless you explicitly drop that compatibility.