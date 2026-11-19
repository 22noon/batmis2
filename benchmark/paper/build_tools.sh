#!/bin/bash
# build_tools.sh: builds every mapper of the paper benchmark with the same compiler and
# flags, on the machine that will run the benchmark, and writes $PREFIX/manifest.json.
#
#   PREFIX=/local/scratch/bench/tools ORIG_TARBALL=/path/BatMis-3.00.tar.gz \
#       benchmark/paper/build_tools.sh [tool ...]
#
# Tools: batmis original bwa bwa_stock bowtie bowtie2 columba mrsfast
# (default: all). Each goes to $PREFIX/bin/<tool>/.
#
#   OPT       optimisation flags for every tool (default -O3 -march=native)
#   JOBS      make -j (default: all CPUs)
#
# BWA is built from the unmodified release, two ways:
#   bwa         OPT -mpopcnt (the flag stated explicitly);
#   bwa_stock   its own Makefile flags (-O3, no -march), for reference.
# The manifest records how many popcnt instructions each binary holds.
set -euo pipefail

HERE=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$HERE/../.." && pwd)
PREFIX=${PREFIX:?set PREFIX, the directory to build the tools in}
OPT=${OPT:--O3 -march=native}
JOBS=${JOBS:-$(nproc)}
SRC=$PREFIX/src
BIN=$PREFIX/bin
mkdir -p "$SRC" "$BIN"

BWA_TAG=v0.7.19
BOWTIE_TAG=v1.3.1
BOWTIE2_TAG=v2.5.4
COLUMBA_TAG=v2.0
MRSFAST_TAG=v3.4.2

TOOLS=${*:-batmis original bwa bwa_stock bowtie bowtie2 columba mrsfast}

log() { echo "[build_tools] $*" >&2; }

clone() {  # clone <github repo> <tag> <dir>
    if [ ! -d "$SRC/$3" ]; then
        git -c advice.detachedHead=false clone -q --depth 1 -b "$2" "https://github.com/$1" "$SRC/$3"
    fi
}

build_batmis() {
    log "batmis (this repository, $(git -C "$REPO" rev-parse --short HEAD))"
    make -C "$REPO/index" clean >/dev/null
    make -C "$REPO/index" -j"$JOBS" CFLAGS="-w $OPT -funroll-loops" libfmindex_core.a bwt_build_index
    make -C "$REPO/src" clean >/dev/null
    make -C "$REPO/src" -j"$JOBS" CXXFLAGS="$OPT -funroll-loops -pthread" batmis1 batmis2
    mkdir -p "$BIN/batmis"
    cp "$REPO/src/batmis1" "$REPO/src/batmis2" "$REPO/index/bwt_build_index" "$BIN/batmis/"
    (cd "$REPO/src" && ../tests/check_prefetch.sh "$BIN/batmis/batmis1" batmis1.cpp &&
        ../tests/check_prefetch.sh "$BIN/batmis/batmis2" batmis2.cpp)
}

build_original() {
    local tarball=${ORIG_TARBALL:?set ORIG_TARBALL to BatMis-3.00.tar.gz}
    log "original BatMis 3.00 from $tarball"
    rm -rf "$SRC/original" && mkdir -p "$SRC/original"
    tar -xzf "$tarball" -C "$SRC/original" --strip-components=1
    # One guard against undefined behaviour, no change of results: Location_To_Genome looks a
    # hit up with map::lower_bound, and for a hit in the reference's last sequence (every hit,
    # in a one-sequence genome) that is end(), whose key it reads before stepping back.
    # Optimised GCC 13 builds crash on it; the guard skips the read.
    python3 - "$SRC/original/src/batdecode.cpp" <<'EOF'
import sys
p = sys.argv[1]; s = open(p).read()
old = "\tif (I->first != Location) I--;"
assert s.count(old) == 1, "batdecode.cpp: Location_To_Genome not as expected"
open(p, "w").write(s.replace(old, "\tif (I == Annotations.end() || I->first != Location) I--;"))
EOF
    # -fcommon: the 2012 sources define globals in headers (GCC >= 10 defaults to -fno-common).
    (cd "$SRC/original" && ./configure -q CFLAGS="$OPT -fcommon -w" CXXFLAGS="$OPT -fcommon -w -fpermissive" LIBS="-lm" &&
        make -j"$JOBS" >/dev/null)
    mkdir -p "$BIN/original"
    cp "$SRC/original/src/batman" "$SRC/original/src/batdecode" "$SRC/original/src/bwtformatdb" "$BIN/original/"
    cp "$SRC/original/scripts/"* "$BIN/original/" 2>/dev/null || true
}

build_bwa_variant() {  # build_bwa_variant <name> <CFLAGS or "" for the Makefile's own>
    log "$1"
    rm -rf "$SRC/$1" && cp -r "$SRC/bwa_release" "$SRC/$1"
    if [ -n "$2" ]; then make -C "$SRC/$1" -j"$JOBS" CFLAGS="-g -Wall -Wno-unused-function $2" >/dev/null
    else make -C "$SRC/$1" -j"$JOBS" >/dev/null; fi
    mkdir -p "$BIN/$1" && cp "$SRC/$1/bwa" "$BIN/$1/"
    echo "${2:-$(sed -n 's/^CFLAGS=[[:space:]]*//p' "$SRC/$1/Makefile") (Makefile)}" > "$BIN/$1/BUILD_FLAGS"
}

build_bowtie() {
    clone BenLangmead/bowtie "$BOWTIE_TAG" bowtie
    log "bowtie $BOWTIE_TAG"
    make -C "$SRC/bowtie" -j"$JOBS" POPCNT_CAPABILITY=1 RELEASE_FLAGS="$OPT" \
        bowtie-align-s bowtie-align-l bowtie-build-s bowtie-build-l >/dev/null
    mkdir -p "$BIN/bowtie"
    cp "$SRC/bowtie/"bowtie{,-build,-align-s,-align-l,-build-s,-build-l} "$BIN/bowtie/"
}

build_bowtie2() {
    clone BenLangmead/bowtie2 "$BOWTIE2_TAG" bowtie2
    log "bowtie2 $BOWTIE2_TAG"
    make -C "$SRC/bowtie2" -j"$JOBS" POPCNT_CAPABILITY=1 RELEASE_FLAGS="$OPT -funroll-loops" \
        bowtie2-align-s bowtie2-align-l bowtie2-build-s bowtie2-build-l >/dev/null
    mkdir -p "$BIN/bowtie2"
    cp "$SRC/bowtie2/"bowtie2{,-build,-align-s,-align-l,-build-s,-build-l} "$BIN/bowtie2/"
}

build_columba() {
    clone biointec/columba "$COLUMBA_TAG" columba
    log "columba $COLUMBA_TAG (Vanilla, 32-bit: texts below 2^32)"
    rm -rf "$SRC/columba/build_bench" && mkdir -p "$SRC/columba/build_bench"
    (cd "$SRC/columba/build_bench" &&
        cmake .. -DRUN_LENGTH_COMPRESSION=OFF -DTHIRTY_TWO=ON -DCMAKE_BUILD_TYPE=Release \
            -DCMAKE_CXX_FLAGS_RELEASE="$OPT -DNDEBUG" >/dev/null &&
        make -j"$JOBS" >/dev/null)
    mkdir -p "$BIN/columba"
    cp "$SRC/columba/build_bench/columba" "$SRC/columba/build_bench/columba_build" "$BIN/columba/"
}

build_mrsfast() {
    clone sfu-compbio/mrsfast "$MRSFAST_TAG" mrsfast
    log "mrsfast $MRSFAST_TAG"
    # Its Makefile adds -O2 in OPTIMIZE_FLAGS and -msse4.2 when the CPU has it; OPT is
    # appended after them, so its -O3 -march=native wins.
    grep -q BENCH_OPT "$SRC/mrsfast/Makefile" ||
        sed -i 's/^OPTIMIZE_FLAGS:$/OPTIMIZE_FLAGS:\n\t$(eval CFLAGS = $(CFLAGS) $(BENCH_OPT))/' "$SRC/mrsfast/Makefile"
    make -C "$SRC/mrsfast" clean >/dev/null
    # -fcommon and -no-pie: the 2017 sources define globals in headers and are compiled
    # with -fno-pic, which current GCC (-fno-common, PIE by default) rejects.
    make -C "$SRC/mrsfast" BENCH_OPT="$OPT -fcommon" LDFLAGS=-no-pie >/dev/null
    mkdir -p "$BIN/mrsfast" && cp "$SRC/mrsfast/mrsfast" "$BIN/mrsfast/"
}

for t in $TOOLS; do
    case $t in
        batmis) build_batmis ;;
        original) build_original ;;
        bwa|bwa_stock)
            if [ ! -d "$SRC/bwa_release" ]; then clone lh3/bwa "$BWA_TAG" bwa_release; fi
            case $t in
                bwa) build_bwa_variant bwa "$OPT -mpopcnt" ;;
                bwa_stock) build_bwa_variant bwa_stock "" ;;
            esac ;;
        bowtie) build_bowtie ;;
        bowtie2) build_bowtie2 ;;
        columba) build_columba ;;
        mrsfast) build_mrsfast ;;
        *) echo "unknown tool $t" >&2; exit 1 ;;
    esac
done

python3 "$HERE/manifest.py" "$PREFIX" "$OPT" > "$PREFIX/manifest.json"
log "manifest: $PREFIX/manifest.json"
python3 - "$PREFIX/manifest.json" <<'EOF'
import json, sys
m = json.load(open(sys.argv[1]))
for name, b in sorted(m["binaries"].items()):
    print(f"  {name:32s} popcnt={b['popcnt']:<6d} prefetch={b['prefetch']}")
EOF
