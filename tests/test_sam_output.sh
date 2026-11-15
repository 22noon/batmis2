#!/usr/bin/env bash
# SAM/BAM output of the mapper ($MAPPER: batmis1, the default, or batmis2) on a multi-sequence reference with N runs, short
# and empty sequences:
#   - Genome_Map::Place against a brute-force scan of every position (test_genome_map);
#   - every record against the genome (SEQ, QUAL, CIGAR, NM, MD, NH, flags), and with
#     -m 0 each read's hits against brute force on both strands (sam_oracle.py);
#   - the same output with --no-ref (positions recorded at sampled rows), --batch 16,
#     threads (-t, with 16-read chunks: <mapper>_chunk16), and as BAM (decoded
#     with htslib's test_view, if built);
#   - Sa_Decode against BWTSaValue on every row (test_sa_decode, if given), of the
#     reference and of a genome with long N runs, homopolymers and tandem repeats.
set -euo pipefail

BUILD_DIR=${1:?usage: test_sam_output.sh BUILD_DIR BWT_BUILD_INDEX GENOME_MAP_TEST [TEST_VIEW [SA_DECODE_TEST]]}
INDEX_BUILDER=${2:?}
MAP_TEST=${3:?}
TEST_VIEW=${4:-}
SA_TEST=${5:-}
TESTS=$(cd "$(dirname "$0")" && pwd)
MAPPER=${MAPPER:-batmis1}
TMP_DIR=$(mktemp -d)
trap 'rm -rf "$TMP_DIR"' EXIT

python3 - "$TMP_DIR" <<'PY'
import random, sys
out = sys.argv[1]
rng = random.Random(5)
ACGT = "ACGT"
def rand_seq(n): return "".join(rng.choice(ACGT) for _ in range(n))
shared = rand_seq(400)
seqs = []
for i, n in enumerate([30000, 20000, 90, 0, 15000, 40, 12000]):
    s = list(rand_seq(n))
    if n > 1000:
        for _ in range(4):  # near-repeats across sequences
            p = rng.randrange(n - 400); blk = list(shared)
            for _ in range(rng.randrange(3)): blk[rng.randrange(400)] = rng.choice(ACGT)
            s[p:p + 400] = blk
        for _ in range(4):
            p = rng.randrange(n); l = rng.randrange(1, 60); s[p:p + l] = "N" * len(s[p:p + l])
        s[:5] = "NNNNN"; s[-3:] = "nnn"
    seqs.append((f"chr{i + 1} test", "".join(s)))
with open(f"{out}/ref.fa", "w") as f:
    for name, s in seqs:
        f.write(f">{name}\n")
        for k in range(0, len(s), 60): f.write(s[k:k + 60] + "\n")
text = "".join(s for _, s in seqs).upper()
COMP = str.maketrans("ACGTN", "TGCAN")
for L in (32, 50, 75):
    with open(f"{out}/reads_{L}.fq", "w") as f:
        for i in range(150):
            if rng.random() < 0.08: r = rand_seq(L)
            else:
                p = rng.randrange(len(text) - L); r = ["C" if c not in ACGT else c for c in text[p:p + L]]
                for j in rng.sample(range(L), rng.randint(0, 7)): r[j] = rng.choice([b for b in ACGT if b != r[j]])
                r = "".join(r)
                if rng.random() < 0.5: r = r.translate(COMP)[::-1]
            if rng.random() < 0.03: j = rng.randrange(L); r = r[:j] + "N" + r[j + 1:]
            q = "".join(chr(33 + rng.randint(2, 40)) for _ in range(L))
            f.write(f"@q{i} planted\n{r}\n+\n{q}\n")
PY

"$INDEX_BUILDER" "$TMP_DIR/ref.fa" > /dev/null
"$MAP_TEST" "$TMP_DIR/ref.fa" "$TMP_DIR/ref.fa"

if [ -n "$SA_TEST" ]; then
    python3 - "$TMP_DIR/homo.fa" <<'PY'
import random, sys
rng = random.Random(9)
parts = []
for i in range(40):
    parts.append("".join(rng.choice("ACGT") for _ in range(rng.randrange(50, 3000))))
    kind = rng.randrange(3)
    if kind == 0: parts.append("N" * rng.randrange(20, 4000))
    elif kind == 1: parts.append(rng.choice("ACGT") * rng.randrange(20, 2000))
    else: parts.append("".join(rng.choice("ACGT") for _ in range(rng.randrange(2, 7))) * rng.randrange(10, 400))
with open(sys.argv[1], "w") as f:
    f.write(">homo\n" + "".join(parts) + "\n")
PY
    "$INDEX_BUILDER" "$TMP_DIR/homo.fa" > /dev/null
    for G in ref homo; do
        Result=$("$SA_TEST" "$TMP_DIR/$G.fa")
        echo "Sa_Decode, $G: $Result"
    done
fi

body() { grep -v '^@PG' "$1"; }
failed=0
for L in 32 50 75; do
    reads="$TMP_DIR/reads_$L.fq"
    for n in 0 1 2 3 4 5 6 7; do
        for m in 0 1; do
            "$BUILD_DIR/$MAPPER" -q "$reads" -g "$TMP_DIR/ref.fa" -n $n -m $m -o "$TMP_DIR/a.sam"
            result=$(python3 "$TESTS/sam_oracle.py" "$TMP_DIR/ref.fa" "$reads" "$TMP_DIR/a.sam" $n $m)
            case "$result" in
                *"bad=0"*) ;;
                *) echo "length $L (-n $n -m $m): $result" >&2; failed=1 ;;
            esac
            for opts in "--no-ref" "--batch 16" "--no-ref --batch 16"; do
                "$BUILD_DIR/$MAPPER" -q "$reads" -g "$TMP_DIR/ref.fa" -n $n -m $m $opts -o "$TMP_DIR/b.sam"
                if ! cmp -s <(body "$TMP_DIR/a.sam") <(body "$TMP_DIR/b.sam"); then
                    echo "length $L (-n $n -m $m $opts): SAM differs" >&2
                    failed=1
                fi
            done
            # Threads: 16-read chunks, so that the reads are shared out.
            for opts in "-t 4" "-t 3 --no-ref" "-t 2 --forward-only"; do
                # shellcheck disable=SC2086
                "$BUILD_DIR/${MAPPER}_chunk16" -q "$reads" -g "$TMP_DIR/ref.fa" -n $n -m $m $opts -o "$TMP_DIR/b.sam"
                if [ "$opts" = "-t 2 --forward-only" ]; then
                    "$BUILD_DIR/$MAPPER" -q "$reads" -g "$TMP_DIR/ref.fa" -n $n -m $m --forward-only -o "$TMP_DIR/f.sam"
                    expected="$TMP_DIR/f.sam"
                else
                    expected="$TMP_DIR/a.sam"
                fi
                if ! cmp -s <(body "$expected") <(body "$TMP_DIR/b.sam"); then
                    echo "length $L (-n $n -m $m $opts, 16-read chunks): SAM differs" >&2
                    failed=1
                fi
            done
            if [ -n "$TEST_VIEW" ] && [ -x "$TEST_VIEW" ]; then
                "$BUILD_DIR/$MAPPER" -q "$reads" -g "$TMP_DIR/ref.fa" -n $n -m $m -o "$TMP_DIR/c.bam"
                if ! cmp -s <(body "$TMP_DIR/a.sam") <("$TEST_VIEW" "$TMP_DIR/c.bam" | grep -v '^@PG'); then
                    echo "length $L (-n $n -m $m): BAM differs from SAM" >&2
                    failed=1
                fi
                # BAM compressed by htslib's threads (-t)
                "$BUILD_DIR/${MAPPER}_chunk16" -q "$reads" -g "$TMP_DIR/ref.fa" -n $n -m $m -t 4 -o "$TMP_DIR/c.bam"
                if ! cmp -s <(body "$TMP_DIR/a.sam") <("$TEST_VIEW" "$TMP_DIR/c.bam" | grep -v '^@PG'); then
                    echo "length $L (-n $n -m $m -t 4): BAM differs from SAM" >&2
                    failed=1
                fi
            fi
        done
    done
done
[ $failed = 0 ] && echo "SAM output test passed"
exit $failed
