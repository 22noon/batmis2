#!/usr/bin/env python3
"""compare_hits.py: scores the hits of the correctness runs of one read set.

  compare_hits.py --genome g.fa --reads set.fq --out scores.json [--truth brute.tsv]
                  --run tool:k:mode:hits.tsv [--run ...]

Hits are normalize_hits.py's TSV (read, contig, pos0, strand, nm, exact). Every reported
hit is checked against the genome: it is valid for K if its CIGAR is exact, its window lies
inside the contig with only A/C/G/T, and the read (or its reverse complement) differs from
it in at most K places.

The truth for K is either the brute-force hits (--truth: read, contig, pos0, strand,
mismatches) with mismatches <= K, or, without it (human), the union of the valid hits that
any exact tool reported in any all-hits run, with mismatches <= K. Approximate tools (k -1)
are not part of the union.

Scores per run:
  all mode:    reported, valid, invalid (not a Hamming hit), over_k (a true alignment with
               more than K mismatches: Columba's best mode can report them), truth, missed (truth hits not reported), and the
               reads with a hit in the truth against those the tool reported;
  first mode:  reads with a truth hit, reported reads, reads whose hit is invalid, reads with
               no hit reported though one exists (missed), and reads whose hit has more
               mismatches than the read's best (not_best);
  both:        for simulated reads (named sim<i>_k<k>_<contig>_<pos0>_<strand>), the share of
               reads whose true origin is among the reported hits, by k; approximate tools
               may place it within 5 bases.
"""
import argparse, collections, json, re
import numpy as np

CODE = np.full(256, 4, dtype=np.uint8)
for i, c in enumerate(b"ACGT"):
    CODE[c] = i
    CODE[c + 32] = i
SIM = re.compile(r"^sim(\d+)_k(\d+)_(.+)_(\d+)_([+-])$")


def load_genome(path):
    names, parts, cur = [], [], []
    with open(path, "rb") as f:
        for line in f:
            if line.startswith(b">"):
                if names:
                    parts.append(b"".join(cur))
                names.append(line[1:].split()[0].decode())
                cur = []
            else:
                cur.append(line.strip())
    if names:
        parts.append(b"".join(cur))
    offset, lengths, pos = {}, {}, 0
    for n, p in zip(names, parts):
        offset[n], lengths[n] = pos, len(p)
        pos += len(p)
    text = CODE[np.frombuffer(b"".join(parts), dtype=np.uint8)]
    return text, offset, lengths


def load_reads(path):
    names, seqs = [], []
    with open(path, "rb") as f:
        while True:
            h = f.readline()
            if not h:
                break
            s = f.readline().strip()
            f.readline(); f.readline()
            names.append(h[1:].split()[0].decode())
            seqs.append(s)
    L = len(seqs[0])
    fwd = CODE[np.frombuffer(b"".join(seqs), dtype=np.uint8)].reshape(-1, L)
    rc = np.where(fwd < 4, 3 - fwd, 4)[:, ::-1]
    return names, {n: i for i, n in enumerate(names)}, fwd, rc


def read_hits(path):
    out = []
    with open(path) as f:
        for line in f:
            r, c, p, s, nm, ex = line.rstrip("\n").split("\t")
            out.append((r, c, int(p), s, int(ex)))
    return out


def mismatches(hits, idx, fwd, rc, text, offset, lengths):
    """Mismatches of each (read, contig, pos0, strand, exact) hit; -1 if it is not a Hamming
    hit (clipped or gapped, or outside its contig), -2 if its window covers a non-ACGT base
    (tools that index N as some base can report these; they are left out, not counted wrong)."""
    L = fwd.shape[1]
    n = len(hits)
    mm = np.full(n, -1, dtype=np.int64)
    ok = np.zeros(n, dtype=bool)
    ridx = np.zeros(n, dtype=np.int64)
    gpos = np.zeros(n, dtype=np.int64)
    minus = np.zeros(n, dtype=bool)
    for i, (r, c, p, s, ex) in enumerate(hits):
        if ex and c in offset and r in idx and 0 <= p and p + L <= lengths[c]:
            ok[i] = True
            ridx[i] = idx[r]
            gpos[i] = offset[c] + p
            minus[i] = s == "-"
    sel = np.nonzero(ok)[0]
    ar = np.arange(L)
    for a in range(0, len(sel), 200000):
        j = sel[a:a + 200000]
        win = text[gpos[j][:, None] + ar]
        rd = np.where(minus[j][:, None], rc[ridx[j]], fwd[ridx[j]])
        d = (win != rd).sum(1)
        d[(win > 3).any(1)] = -2
        mm[j] = d
    return mm


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--reads", required=True)
    ap.add_argument("--truth")
    ap.add_argument("--run", action="append", default=[])
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    names, idx, fwd, rc = load_reads(a.reads)
    text, offset, lengths = load_genome(a.genome)
    origin = {}
    for n in names:
        m = SIM.match(n)
        if m:
            origin[n] = (m.group(3), int(m.group(4)), m.group(5), int(m.group(2)))

    runs = []
    for spec in a.run:
        tool, k, mode, path = spec.split(":", 3)
        hits = read_hits(path)
        mm = mismatches(hits, idx, fwd, rc, text, offset, lengths)
        runs.append({"tool": tool, "k": int(k), "mode": mode, "hits": hits, "mm": mm})

    # Truth: hit key -> mismatches.
    truth = {}
    if a.truth:
        with open(a.truth) as f:
            for line in f:
                r, c, p, s, d = line.rstrip("\n").split("\t")
                truth[(r, c, int(p), s)] = int(d)
        source = "brute_force"
    else:
        for R in runs:
            if R["k"] < 0 or R["mode"] != "all":
                continue
            for h, d in zip(R["hits"], R["mm"]):
                if 0 <= d <= R["k"]:
                    truth[h[:4]] = int(d)
        source = "verified_union"

    def truth_for(k):
        return {h: d for h, d in truth.items() if d <= k}

    best_cache = {}

    def best_of(k):  # read -> fewest mismatches among its truth hits
        if k not in best_cache:
            b = {}
            for h, d in truth_for(k).items():
                b[h[0]] = min(d, b.get(h[0], 99))
            best_cache[k] = b
        return best_cache[k]

    scores = []
    for R in runs:
        k, hits, mm = R["k"], R["hits"], R["mm"]
        s = {"tool": R["tool"], "k": k, "mode": R["mode"], "truth_source": source, "reported": len(hits)}
        if k >= 0:
            valid = {h[:4]: int(d) for h, d in zip(hits, mm) if 0 <= d <= k}
            s["valid"] = len(valid)
            s["over_k"] = int(((mm > k)).sum())  # true alignments with more than K mismatches
            s["invalid"] = int((mm == -1).sum())  # not a Hamming hit at all
            s["over_n"] = int((mm == -2).sum())  # a window over reference N: left out
            T = truth_for(k)
            if R["mode"] == "all":
                s["truth"] = len(T)
                s["missed"] = sum(1 for h in T if h not in valid)
                s["extra_not_in_truth"] = sum(1 for h in valid if h not in T)
                s["truth_reads"] = len({h[0] for h in T})
                s["reported_reads"] = len({h[0] for h in valid})
            else:
                best = best_of(k)
                got = collections.defaultdict(lambda: 99)
                bad, over = set(), set()
                for h, d in zip(hits, mm):
                    if 0 <= d <= k:
                        got[h[0]] = min(got[h[0]], int(d))
                    elif d > k:
                        over.add(h[0])
                    elif d == -1:
                        bad.add(h[0])
                s["truth_reads"] = len(best)
                s["reported_reads"] = len(set(got) | bad | over)
                s["reads_invalid_only"] = len(bad - set(got))
                s["reads_over_k_only"] = len(over - set(got) - bad)
                s["missed_reads"] = sum(1 for r in best if r not in got)
                s["not_best"] = sum(1 for r, d in got.items() if r in best and d > best[r])
        # True origin, by k of the simulation.
        if origin:
            found = collections.Counter()
            total = collections.Counter()
            near = 5 if k < 0 else 0
            rep = collections.defaultdict(list)
            for h, d in zip(hits, mm):
                rep[h[0]].append(h)
            for n, (c, p, st, kk) in origin.items():
                total[kk] += 1
                if any(h[1] == c and h[3] == st and abs(h[2] - p) <= near for h in rep.get(n, ())):
                    found[kk] += 1
            s["origin_found_by_k"] = {str(kk): [found[kk], total[kk]] for kk in sorted(total)}
        scores.append(s)
    json.dump({"reads": a.reads, "genome": a.genome, "scores": scores}, open(a.out, "w"), indent=1)
    for s in scores:
        keys = ["missed", "invalid", "over_k", "not_best", "missed_reads"]
        print(f"  {s['tool']:10s} k={s['k']:2d} {s['mode']:5s} " + " ".join(f"{x}={s[x]}" for x in keys if x in s))


if __name__ == "__main__":
    main()
