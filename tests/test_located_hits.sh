#!/usr/bin/env bash
# Located hits (single-row paths finished against .pac from a sampled SA row)
# must match walking the index: runs test_located_hits on a toy reference, then
# compares the mapper ($MAPPER: batmis1, the default, or batmis2) end to end with and without --no-ref, --batch 16, threads (-t)
# and gzip input.
set -euo pipefail

BUILD_DIR=${1:?usage: test_located_hits.sh BUILD_DIR [BWT_BUILD_INDEX] [LOCATED_TEST]}
INDEX_BUILDER=${2:-$BUILD_DIR/bwt_build_index}
LOCATED_TEST=${3:-$BUILD_DIR/test_located_hits}
MAPPER=${MAPPER:-batmis1}
TMP_DIR=$(mktemp -d)
trap 'rm -rf "$TMP_DIR"' EXIT

python3 - "$TMP_DIR" <<'PY'
import random
import sys

out = sys.argv[1]
rng = random.Random(11)
ACGT = "ACGT"

def rand_seq(n):
    return "".join(rng.choice(ACGT) for _ in range(n))

def mutate(s, positions):
    s = list(s)
    for p in positions:
        s[p] = rng.choice([c for c in ACGT if c != s[p]])
    return "".join(s)

# Two contigs with a repeated segment, so that some reads hit several rows.
contigs = [rand_seq(20000), rand_seq(15000)]
segment = rand_seq(150)
for _ in range(8):
    c = rng.randrange(2)
    p = rng.randrange(len(contigs[c]) - 150)
    copy = mutate(segment, rng.sample(range(150), rng.randrange(3)))
    contigs[c] = contigs[c][:p] + copy + contigs[c][p + 150:]
with open(f"{out}/ref.fa", "w") as f:
    for i, s in enumerate(contigs):
        f.write(f">chr{i + 1}\n")
        for j in range(0, len(s), 80):
            f.write(s[j:j + 80] + "\n")
genome = "".join(contigs)

for L in (12, 20, 36, 75):
    lh = L // 2
    with open(f"{out}/reads{L}.fa", "w") as f:
        n = 0
        def put(seq):
            global n
            f.write(f">q{n}\n{seq}\n")
            n += 1
        for _ in range(150):
            p = rng.randrange(len(genome) - L)
            r = genome[p:p + L]
            put(r)
            put(mutate(r, [rng.randrange(lh)]))
            put(mutate(r, [rng.randrange(lh, L)]))
            put(mutate(r, rng.sample(range(L), 2)))
            put(mutate(r, rng.sample(range(L), rng.randrange(3, 8))))  # for -n 3 to 7
        # Reads across the end of each contig, and overhanging the genome ends.
        for k in range(1, L):
            put(genome[len(contigs[0]) - k:len(contigs[0]) - k + L])
            put(rand_seq(L - k) + genome[:k])
            put(genome[len(genome) - k:] + rand_seq(L - k))
            put(mutate(genome[len(genome) - k:] + rand_seq(L - k), [rng.randrange(L)]))
PY

"$INDEX_BUILDER" -o "$TMP_DIR" "$TMP_DIR/ref.fa" >"$TMP_DIR/index.log" 2>&1

for length in 12 20 36 75; do
    reads="$TMP_DIR/reads${length}.fa"
    for max_hits in 0 1; do
        "$LOCATED_TEST" "$TMP_DIR/ref.fa" "$reads" "$max_hits"
    done
    # End to end: the same number of locations per read with and without the reference.
    for opts in "-n 0 -m 0" "-n 1 -m 0" "-n 1 -m 1"; do
        # shellcheck disable=SC2086
        "$BUILD_DIR/$MAPPER" -q "$reads" -g "$TMP_DIR/ref.fa" $opts -o "$TMP_DIR/ref.tsv"
        # shellcheck disable=SC2086
        "$BUILD_DIR/$MAPPER" -q "$reads" -g "$TMP_DIR/ref.fa" $opts --no-ref -o "$TMP_DIR/noref.tsv"
        for f in ref noref; do
            awk -F'\t' '{n[$1] += $5} END {for (q in n) print q, n[q]}' "$TMP_DIR/$f.tsv" | sort >"$TMP_DIR/$f.counts"
        done
        if ! cmp -s "$TMP_DIR/ref.counts" "$TMP_DIR/noref.counts"; then
            echo "length $length ($opts): hit counts differ with and without --no-ref" >&2
            exit 1
        fi
    done
    # Batched search: the same output, byte for byte, as one read at a time.
    for ref in "" "--no-ref"; do
        for n in 0 1 2 3 4 5 6 7; do
            for m in 0 1; do
                # shellcheck disable=SC2086
                "$BUILD_DIR/$MAPPER" -q "$reads" -g "$TMP_DIR/ref.fa" -n $n -m $m $ref -o "$TMP_DIR/one.tsv"
                # shellcheck disable=SC2086
                "$BUILD_DIR/$MAPPER" -q "$reads" -g "$TMP_DIR/ref.fa" -n $n -m $m $ref --batch 16 -o "$TMP_DIR/batch.tsv"
                if ! cmp -s "$TMP_DIR/one.tsv" "$TMP_DIR/batch.tsv"; then
                    echo "length $length (-n $n -m $m $ref): --batch 16 output differs" >&2
                    exit 1
                fi
                # Threads: 16-read chunks, so that the reads are shared out.
                for opts in "-t 1" "-t 4" "-t 3 --batch 4"; do
                    # shellcheck disable=SC2086
                    "$BUILD_DIR/${MAPPER}_chunk16" -q "$reads" -g "$TMP_DIR/ref.fa" -n $n -m $m $ref $opts \
                        -o "$TMP_DIR/threads.tsv"
                    if ! cmp -s "$TMP_DIR/one.tsv" "$TMP_DIR/threads.tsv"; then
                        echo "length $length (-n $n -m $m $ref $opts, 16-read chunks): output differs" >&2
                        exit 1
                    fi
                done
                # The reads gzip-compressed
                if [ $m = 0 ]; then
                    gzip -c "$reads" > "$TMP_DIR/reads.gz"
                    # shellcheck disable=SC2086
                    "$BUILD_DIR/$MAPPER" -q "$TMP_DIR/reads.gz" -g "$TMP_DIR/ref.fa" -n $n -m $m $ref -t 2 \
                        -o "$TMP_DIR/gz.tsv"
                    if ! cmp -s "$TMP_DIR/one.tsv" "$TMP_DIR/gz.tsv"; then
                        echo "length $length (-n $n -m $m $ref): output differs for gzip input" >&2
                        exit 1
                    fi
                fi
            done
        done
    done
done

echo "located-hit differential test passed"
