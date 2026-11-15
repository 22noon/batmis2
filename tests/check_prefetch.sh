#!/bin/sh
# check_prefetch.sh <binary> <source>: fails if the mapper's prefetches are not in the
# binary.
#
# GCC removes a call to a function that only prefetches when it does not inline it
# (src/fm_occ.hpp, PREFETCH_ONLY), and nothing else shows it: the output is the same and
# only the batched search gets slower. So this checks two things:
#   1. source: every function named *Prefetch* is declared PREFETCH_ONLY;
#   2. binary: the mapper's own code (its anonymous namespace; not the index library or
#      htslib) holds at least MIN prefetch instructions. With the prefetch functions
#      dropped it holds 7 (Sa_Decode's, unrolled); with them, about 45 (GCC 13.3 -O3).
# An unoptimised build (-O0) has one instruction per call site and fails this; set
# PREFETCH_MIN to change the floor, or PREFETCH_MIN=0 to skip.
set -e
BIN=$1
SRC=$2
MIN=${PREFETCH_MIN:-20}

BAD=$(grep -nE '^[A-Za-z_ ]*\b[A-Za-z_]*Prefetch[A-Za-z_]*\(.*\) *\{' "$SRC" | grep -v PREFETCH_ONLY || true)
if [ -n "$BAD" ]; then
    echo "check_prefetch: $SRC: prefetch function not declared PREFETCH_ONLY:" >&2
    echo "$BAD" >&2
    exit 1
fi

N=$(objdump -d -C --no-show-raw-insn "$BIN" | awk '
    /^[0-9a-f]+ <.*>:$/ { own = index($0, "(anonymous namespace)") > 0 }
    own && /\tprefetch/ { n++ }
    END { print n + 0 }')
if [ "$N" -lt "$MIN" ]; then
    echo "check_prefetch: $BIN has $N prefetch instructions of its own, expected at least $MIN:" >&2
    echo "  the compiler has dropped them (see src/fm_occ.hpp, PREFETCH_ONLY)" >&2
    exit 1
fi
echo "check_prefetch: $BIN: $N prefetch instructions"
