#!/usr/bin/env bash
# test_batmis2.sh BUILD_DIR BWT_BUILD_INDEX: batmis2 against batmis1 on a toy reference
# (several sequences, near-repeats, an N run) with reads of 32, 50 and 75 bases and 0-4
# mismatches:
#   - with --switch-walk --own-row (batmis1's switch and locate), the output is
#     byte-identical to batmis1's (SAM: but for the @PG line), for -n 0 to 7, text and SAM,
#     -m 5 and -m 0, one read at a time, batched, threads, and --no-ref;
#   - as is, the same hits: -m 0 SAM records equal once sorted (a hit located from the
#     other row is found in another order, so with -m N the first N may differ);
#   - batmis2_verify (-DBATMIS2_VERIFY) finds no range whose other-index start, and no
#     other-row locate, disagrees with the other index.
set -euo pipefail
BUILD_DIR=${1:?usage: test_batmis2.sh BUILD_DIR BWT_BUILD_INDEX}
INDEX_BUILDER=${2:?}
# SA_INTERVAL=N builds the test indexes with SA samples every N rows (bwt_build_index
# --sa-interval); unset, the builder's memory budget decides (every row, for these sizes).
SA_OPT=${SA_INTERVAL:+--sa-interval $SA_INTERVAL}
B1=$BUILD_DIR/batmis1
B2=$BUILD_DIR/batmis2
V2=$BUILD_DIR/batmis2_verify
TMP_DIR=$(mktemp -d)
trap 'rm -rf "$TMP_DIR"' EXIT

python3 - "$TMP_DIR" <<'PY'
import random, sys
out = sys.argv[1]
rng = random.Random(23)
ACGT = "ACGT"
def rand_seq(n): return "".join(rng.choice(ACGT) for _ in range(n))
repeat = rand_seq(500)
seqs = []
for n in (60000, 35000, 400, 20000):
    s = list(rand_seq(n))
    for _ in range(n // 8000 + 1):  # near-repeats
        p = rng.randrange(max(1, n - 500)); r = list(repeat)
        for _ in range(rng.randrange(3)): r[rng.randrange(500)] = rng.choice(ACGT)
        s[p:p + 500] = r[:len(s) - p]
    seqs.append("".join(s)[:n])
seqs[3] = seqs[3][:9000] + "N" * 800 + seqs[3][9800:]
with open(f"{out}/ref.fa", "w") as f:
    for i, s in enumerate(seqs):
        f.write(f">s{i}\n")
        for k in range(0, len(s), 70): f.write(s[k:k + 70] + "\n")
text = "".join(s.replace("N", "C") for s in seqs)
def comp(s): return s[::-1].translate(str.maketrans("ACGT", "TGCA"))
for L in (32, 50, 75):
    with open(f"{out}/reads_{L}.fq", "w") as f:
        for k in range(300):
            p = rng.randrange(len(text) - L); s = list(text[p:p + L])
            for j in rng.sample(range(L), rng.choice((0, 0, 1, 1, 2, 3, 4))): s[j] = rng.choice([c for c in ACGT if c != s[j]])
            s = "".join(s)
            if k % 2: s = comp(s)
            if k % 50 == 7: s = rand_seq(L)
            f.write(f"@r{L}_{k}\n{s}\n+\n{'I' * L}\n")
PY
"$INDEX_BUILDER" $SA_OPT "$TMP_DIR/ref.fa" >"$TMP_DIR/index.log" 2>&1

failed=0
body() { grep -v '^@PG' "$1"; }
records() { grep -v '^@' "$1" | sort; }
for L in 32 50 75; do
    reads="$TMP_DIR/reads_$L.fq"
    for n in 0 1 2 3 4 5 6 7; do
        for opts in "-m 5" "-m 0 --batch 1" "-m 5 --batch 8" "-m 0 -t 3" "-m 5 --no-ref" "-m 5 --sam" "-m 0 --sam -t 3" \
                    "-m 0 --sam --no-ref"; do
            "$B1" -q "$reads" -g "$TMP_DIR/ref.fa" -n $n $opts -o "$TMP_DIR/one.out" 2>/dev/null
            "$B2" -q "$reads" -g "$TMP_DIR/ref.fa" -n $n $opts --switch-walk --own-row -o "$TMP_DIR/two.out" 2>/dev/null
            if ! cmp -s <(body "$TMP_DIR/one.out") <(body "$TMP_DIR/two.out"); then
                echo "batmis2 --switch-walk --own-row differs from batmis1: L=$L -n $n $opts" >&2
                failed=1
            fi
        done
        [ $n -le 5 ] || continue
        for opts in "" "-t 3" "--batch 1"; do
            "$B1" -q "$reads" -g "$TMP_DIR/ref.fa" -n $n -m 0 --sam $opts -o "$TMP_DIR/one.sam" 2>/dev/null
            "$B2" -q "$reads" -g "$TMP_DIR/ref.fa" -n $n -m 0 --sam $opts -o "$TMP_DIR/two.sam" 2>/dev/null
            if ! cmp -s <(records "$TMP_DIR/one.sam") <(records "$TMP_DIR/two.sam"); then
                echo "batmis2's hits differ from batmis1's: L=$L -n $n $opts" >&2
                failed=1
            fi
        done
        result=$("$V2" -q "$reads" -g "$TMP_DIR/ref.fa" -n $n -m 0 --batch 4 -o /dev/null 2>&1 | grep 'Verify_Other: checked' || true)
        if ! [[ "$result" =~ failed=0$ ]]; then
            echo "batmis2_verify: L=$L -n $n: ${result:-no report}" >&2
            failed=1
        fi
    done
done
[ $failed = 0 ] && echo "batmis2 test passed"
exit $failed
