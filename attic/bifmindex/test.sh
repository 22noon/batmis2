#!/usr/bin/env bash
# test.sh BUILD_DIR INDEX_BUILDER: builds a toy index (three contigs, repeats, a short
# contig) with bwt_build_index, then
#   - bi_check: the bidirectional step against the classic searches and the text;
#   - bi_exact: every mode (-m 0) against brute force over the indexed text (.bin), on
#     reads from the text, with one mismatch, random, and at the text's two ends;
#   - bi_mm: the one-mismatch search (-m 0) against brute force, on the same reads and
#     1-5 base ones;
#   - bi_mm -n 2: the two-mismatch search, likewise;
#   all with each locate route (bi_locate.hpp).
set -euo pipefail
BUILD_DIR=${1:-.}
INDEX_BUILDER=${2:-../index/bwt_build_index}
TMP_DIR=$(mktemp -d)
trap 'rm -rf "$TMP_DIR"' EXIT

python3 - "$TMP_DIR" <<'PY'
import random, sys
out = sys.argv[1]
rng = random.Random(11)
def rand_seq(n): return "".join(rng.choice("ACGT") for _ in range(n))
repeat = rand_seq(300)
contigs = []
for n in (40000, 25000, 300):
    s = list(rand_seq(n))
    for _ in range(n // 5000 + 1):  # copies of a repeat, some with a change
        p = rng.randrange(max(1, n - 300))
        r = list(repeat)
        if rng.random() < 0.5: r[rng.randrange(300)] = rng.choice("ACGT")
        s[p:p + 300] = r[:len(s) - p]
    contigs.append("".join(s)[:n])
with open(f"{out}/ref.fa", "w") as f:
    for i, c in enumerate(contigs):
        f.write(f">c{i}\n")
        for k in range(0, len(c), 70): f.write(c[k:k + 70] + "\n")
PY

"$INDEX_BUILDER" -o "$TMP_DIR" "$TMP_DIR/ref.fa" >"$TMP_DIR/index.log" 2>&1
"$BUILD_DIR/bi_check" "$TMP_DIR/ref.fa" 3000 5

python3 - "$TMP_DIR" <<'PY'
import random, sys
out = sys.argv[1]
rng = random.Random(12)
text = open(f"{out}/ref.fa.bin").read().upper()
n = len(text)
with open(f"{out}/reads.fa", "w") as f:
    k = 0
    def put(s):
        global k
        f.write(f">r{k}\n{s}\n"); k += 1
    for L in (8, 12, 20, 36, 75):
        for _ in range(200):
            p = rng.randrange(n - L); s = text[p:p + L]
            put(s)
            j = rng.randrange(L); put(s[:j] + rng.choice([b for b in "ACGT" if b != s[j]]) + s[j + 1:])
            put("".join(rng.choice("ACGT") for _ in range(L)))
        put(text[:L]); put(text[n - L:])  # the text's ends
    put("ACGTNACGT")
# Brute force: every occurrence, as bi_exact prints with -m 0.
with open(f"{out}/want.txt", "w") as w:
    for line in open(f"{out}/reads.fa"):
        if line.startswith(">"): rid = line[1:].strip(); continue
        s = line.strip()
        pos = [] if "N" in s else [p for p in range(n - len(s) + 1) if text.startswith(s, p)]
        w.write("\t".join([rid, str(len(pos))] + [str(p) for p in pos]) + "\n")
# bi_mm's reads: the same and 1-5 base ones. Brute force: every position where the read
# has at most one mismatch (one half of it occurs exactly), as bi_mm prints with -m 0.
with open(f"{out}/mm_reads.fa", "w") as f:
    f.write(open(f"{out}/reads.fa").read())
    for L in (1, 2, 3, 4, 5):
        for k in range(20): f.write(f">s{L}_{k}\n{''.join(rng.choice('ACGT') for _ in range(L))}\n")
with open(f"{out}/want_mm.txt", "w") as w:
    for line in open(f"{out}/mm_reads.fa"):
        if line.startswith(">"): rid = line[1:].strip(); continue
        s = line.strip(); L = len(s); H = L // 2
        if "N" in s:
            w.write(f"{rid}\t0\t0\n")
            continue
        cand = set(range(n - L + 1)) if H == 0 else set()
        for part, off in ((s[:H], 0), (s[H:], H)) if H else ():
            q = text.find(part)
            while q >= 0:
                if 0 <= q - off <= n - L: cand.add(q - off)
                q = text.find(part, q + 1)
        hits = [(p, d) for p in sorted(cand) for d in [sum(a != b for a, b in zip(s, text[p:p + L]))] if d <= 1]
        e = sum(1 for _, d in hits if d == 0)
        w.write(f"{rid}\t{e}\t{len(hits) - e}" + "".join(f"\t{p}:{d}" for p, d in hits) + "\n")
# bi_mm -n 2: at most two mismatches (one third of the read occurs exactly).
with open(f"{out}/want_mm2.txt", "w") as w:
    for line in open(f"{out}/mm_reads.fa"):
        if line.startswith(">"): rid = line[1:].strip(); continue
        s = line.strip(); L = len(s); a, b = L // 3, 2 * L // 3
        if "N" in s:
            w.write(f"{rid}\t0\t0\t0\n")
            continue
        cand = set(range(n - L + 1)) if a == 0 else set()
        for lo, hi in ((0, a), (a, b), (b, L)) if a else ():
            q = text.find(s[lo:hi])
            while q >= 0:
                if 0 <= q - lo <= n - L: cand.add(q - lo)
                q = text.find(s[lo:hi], q + 1)
        hits = [(p, d) for p in sorted(cand) for d in [sum(x != y for x, y in zip(s, text[p:p + L]))] if d <= 2]
        c = [sum(1 for _, d in hits if d == k) for k in range(3)]
        w.write(f"{rid}\t{c[0]}\t{c[1]}\t{c[2]}" + "".join(f"\t{p}:{d}" for p, d in hits) + "\n")
PY

failed=0
# The bidirectional modes also with the general step only, with lookup tables (depth 10
# is longer than the 8-base reads, which then search from the empty pattern), and
# finished against the packed genome (--ref), batched, and with the library Occ; and each
# locate route: the library walk, range chains on every hit of 2+ rows (with the table's
# rebuilt levels), split ranges, one-end walks.
for run in "backward" "forward" "bi-left" "bi-right" "bi-left --threshold 0" "bi-right --threshold 0" \
           "bi-left --lookup 6" "bi-right --lookup 6" "bi-left --lookup 10" "bi-right --lookup 10" \
           "bi-left --ref" "bi-right --ref" "bi-left --ref --lookup 6" "bi-right --ref --lookup 6" \
           "bi-left --batch 1" "bi-right --batch 5" "bi-left --batch 7 --lookup 6" "bi-right --batch 16 --lookup 10" \
           "bi-left --batch 4 --ref" "bi-right --batch 3 --ref --lookup 6" "bi-left --batch 16 --ref --lookup 6" \
           "bi-left --occ lib" "bi-right --occ lib --batch 5 --ref --lookup 6" \
           "bi-left --batch 6 --early --ref --lookup 6" "bi-right --batch 2 --early" \
           "bi-right --locate lib" "bi-left --chain-min 2" "bi-right --chain-min 2 --lookup 6" \
           "bi-left --batch 5 --chain-min 2 --lookup 6 --ref" "bi-right --batch 3 --chain-min 2 --lookup 10" \
           "bi-left --split --one-end" "bi-right --split --lookup 6"; do
    "$BUILD_DIR/bi_exact" -g "$TMP_DIR/ref.fa" -q "$TMP_DIR/reads.fa" -mode $run -m 0 -o "$TMP_DIR/out.txt"
    if ! cmp -s "$TMP_DIR/want.txt" "$TMP_DIR/out.txt"; then
        echo "bi_exact -mode $run differs from brute force" >&2
        diff "$TMP_DIR/want.txt" "$TMP_DIR/out.txt" | head -5 >&2
        failed=1
    fi
done
# The one-mismatch search: one read at a time and batched, with and without the table and
# the packed genome, and the locate routes as above.
for run in "" "--lookup 6" "--ref" "--ref --lookup 6" "--batch 1" "--batch 8 --ref --lookup 6" "--batch 3 --lookup 10" \
           "--batch 16 --ref" "--locate lib" "--chain-min 2 --lookup 6" "--chain-min 2 --batch 8 --ref --lookup 6" \
           "--split" "--one-end --chain-min 2"; do
    "$BUILD_DIR/bi_mm" -g "$TMP_DIR/ref.fa" -q "$TMP_DIR/mm_reads.fa" $run -o "$TMP_DIR/out.txt"
    if ! cmp -s "$TMP_DIR/want_mm.txt" "$TMP_DIR/out.txt"; then
        echo "bi_mm $run differs from brute force" >&2
        diff "$TMP_DIR/want_mm.txt" "$TMP_DIR/out.txt" | head -5 >&2
        failed=1
    fi
done
# The two-mismatch search: with and without the table (depth 10: from the empty pattern
# for short parts), the packed genome, pass 1's settling and bounds, chains, the library
# walk, and batched.
for run in "" "--lookup 6" "--ref" "--ref --lookup 6" "--ref --lookup 10 --no-settle" "--lookup 6 --no-settle" \
           "--ref --lookup 6 --chain-min 2" "--locate lib" "--batch 1" "--batch 8 --ref --lookup 6" \
           "--batch 5 --lookup 10 --no-settle --no-bound" "--batch 16 --ref --chain-min 2"; do
    "$BUILD_DIR/bi_mm" -g "$TMP_DIR/ref.fa" -q "$TMP_DIR/mm_reads.fa" -n 2 $run -o "$TMP_DIR/out.txt"
    if ! cmp -s "$TMP_DIR/want_mm2.txt" "$TMP_DIR/out.txt"; then
        echo "bi_mm -n 2 $run differs from brute force" >&2
        diff "$TMP_DIR/want_mm2.txt" "$TMP_DIR/out.txt" | head -5 >&2
        failed=1
    fi
done
[ $failed = 0 ] && echo "bidirectional index test passed"
exit $failed
