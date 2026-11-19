#!/usr/bin/env python3
"""discrepancies.py: why the exact tools' hits differ from the truth, for one read set.

  discrepancies.py --genome g.fa --reads set.fq --out DIR [--truth brute.tsv]
                   --run tool:k:mode:hits.tsv [--run ...] [--mrsfast-ws 12]
  discrepancies.py ... --read NAME        one read: the truth and every tool's hits, aligned

Inputs are those of compare_hits.py (whose loaders it uses), and the truth is chosen the
same way: the brute-force hits if --truth is given, else the union of the valid hits of
every exact tool's all-hits runs. Approximate tools (k -1) are not analysed.

Every missed or extra hit gets one cause, the first that applies:

extra (reported, not a truth hit within K):
  over_n             the window covers a reference N (tools that index N as a base)
  clipped_or_gapped  the CIGAR is not one full-length match
  off_contig         unknown contig, or the window runs off its end
  over_k             a true alignment with more than K mismatches
  not_in_truth       a valid hit the truth lacks (with brute force: a truth bug; never with
                     the union, which contains every valid exact hit)
missed (a truth hit within K the tool did not report):
  window_limit       mrsFAST, whose pigeonhole needs K + 1 whole hash windows: L < (K+1)*ws
  overlapping        the tool reports another hit of the read, same contig and strand,
                     less than L bases away (overlapping occurrences in tandem repeats)
  many_hits_read     the read has more than 100 truth hits and the tool reports some
  contig_edge        the window lies within L bases of its contig's end
  low_complexity     the window's dinucleotide entropy is under 3 bits (of 4)
  unexplained        none of the above: the rows to look at by hand
First-hit runs: a read with a truth hit but no valid reported hit is a missed_read (with the
cause of its best truth hit); a read whose hit has more mismatches than its best is not_best.

Writes, in DIR:
  summary.tsv                  tool, k, mode, kind, cause, count
  pairwise.tsv                 k, tool A, tool B: valid hits of A's all-hits run that B's lacks
  <tool>_k<K>_<mode>.tsv       every discrepancy: kind, cause, read, contig, pos0, strand,
                               mismatches, the read's truth hits
  examples.txt                 up to 10 per tool, run and cause, aligned against the genome
"""
import argparse, collections, itertools, math, os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from compare_hits import SIM, load_genome, load_reads, mismatches, read_hits  # noqa: E402

BASES = "ACGTN"
MANY_HITS = 100
LOW_ENTROPY_BITS = 3.0
EXAMPLES = 10


def window(text, offset, c, p, L):
    return "".join(BASES[x] for x in text[offset[c] + p:offset[c] + p + L])


def read_as_aligned(seq, strand):
    """The read as it lies on the forward strand."""
    if strand == "+":
        return seq
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def entropy2(s):
    pairs = collections.Counter(s[i:i + 2] for i in range(len(s) - 1))
    n = sum(pairs.values())
    return -sum(v / n * math.log2(v / n) for v in pairs.values()) if n else 0.0


def aligned_view(seq, text, offset, lengths, hit, L):
    r, c, p, s = hit[:4]
    if c not in offset or p < 0 or p + L > lengths[c]:
        return [f"    {c}:{p}{s}  (outside the contig)"]
    g = window(text, offset, c, p, L)
    a = read_as_aligned(seq, s)
    marks = "".join(" " if x == y else ("N" if y == "N" else "*") for x, y in zip(a, g))
    return [f"    {c}:{p}{s}  {sum(m != ' ' for m in marks)} mismatches",
            f"      genome {g}", f"      read   {a}", f"             {marks}"]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--reads", required=True)
    ap.add_argument("--truth")
    ap.add_argument("--run", action="append", default=[])
    ap.add_argument("--out")
    ap.add_argument("--read", help="show one read instead of writing the analysis")
    ap.add_argument("--mrsfast-ws", type=int, default=12)
    a = ap.parse_args()

    names, idx, fwd, rc = load_reads(a.reads)
    L = fwd.shape[1]
    seqs = {n: "".join(BASES[x] for x in fwd[i]) for n, i in idx.items()} if a.read else None
    text, offset, lengths = load_genome(a.genome)

    runs = []
    for spec in a.run:
        tool, k, mode, path = spec.split(":", 3)
        if int(k) < 0:
            continue  # approximate tools
        hits = read_hits(path)
        runs.append({"tool": tool, "k": int(k), "mode": mode, "hits": hits,
                     "mm": mismatches(hits, idx, fwd, rc, text, offset, lengths)})

    truth = {}
    if a.truth:
        with open(a.truth) as f:
            for line in f:
                r, c, p, s, d = line.rstrip("\n").split("\t")
                if not a.read or r == a.read:
                    truth[(r, c, int(p), s)] = int(d)
    else:
        for R in runs:
            if R["mode"] == "all":
                for h, d in zip(R["hits"], R["mm"]):
                    if 0 <= d <= R["k"]:
                        truth[h[:4]] = int(d)
    truth_by_read = collections.defaultdict(list)
    for h, d in truth.items():
        truth_by_read[h[0]].append((h, d))

    if a.read:
        show_read(a.read, seqs, truth_by_read, runs, text, offset, lengths, L)
        return

    os.makedirs(a.out, exist_ok=True)
    summary = collections.Counter()
    examples = []

    def missed_cause(tool, k, hit, reported_valid):
        r, c, p, s = hit
        if tool == "mrsfast" and L < (k + 1) * a.mrsfast_ws:
            return "window_limit"
        for (rr, cc, pp, ss) in reported_valid.get(r, ()):
            if cc == c and ss == s and abs(pp - p) < L:
                return "overlapping"
        n_truth = sum(1 for _, d in truth_by_read[r] if d <= k)
        if n_truth > MANY_HITS and reported_valid.get(r):
            return "many_hits_read"
        if p < L or p + 2 * L > lengths.get(c, 0):
            return "contig_edge"
        if entropy2(window(text, offset, c, p, L)) < LOW_ENTROPY_BITS:
            return "low_complexity"
        return "unexplained"

    def extra_cause(hit, ex, d, k):
        if d == -2:
            return "over_n"
        if d == -1:
            return "clipped_or_gapped" if not ex else "off_contig"
        if d > k:
            return "over_k"
        return "not_in_truth"

    for R in runs:
        tool, k, mode = R["tool"], R["k"], R["mode"]
        rows = []
        valid = {}
        reported_valid = collections.defaultdict(list)
        for h, d in zip(R["hits"], R["mm"]):
            if 0 <= d <= k:
                valid[h[:4]] = int(d)
                reported_valid[h[0]].append(h[:4])
        n_truth = {r: sum(1 for _, d in hs if d <= k) for r, hs in truth_by_read.items()}
        # Extra hits, both modes.
        for h, d in zip(R["hits"], R["mm"]):
            key = h[:4]
            if key in truth and truth[key] <= k:
                continue
            cause = extra_cause(key, h[4], int(d), k)
            rows.append(("extra", cause, *key, int(d), n_truth.get(h[0], 0)))
        if mode == "all":
            for key, d in truth.items():
                if d <= k and key not in valid:
                    rows.append(("missed", missed_cause(tool, k, key, reported_valid), *key, d, n_truth[key[0]]))
        else:
            for r, hs in truth_by_read.items():
                hs = [(h, d) for h, d in hs if d <= k]
                if not hs:
                    continue
                best_h, best = min(hs, key=lambda x: x[1])
                got = [valid[h] for h in reported_valid.get(r, ())]
                if not got:
                    rows.append(("missed_read", missed_cause(tool, k, best_h, reported_valid), *best_h, best, len(hs)))
                elif min(got) > best:
                    rows.append(("not_best", "not_best", *best_h, best, len(hs)))
        rows.sort()
        with open(os.path.join(a.out, f"{tool}_k{k}_{mode}.tsv"), "w") as o:
            o.write("kind\tcause\tread\tcontig\tpos0\tstrand\tmismatches\tread_truth_hits\n")
            for row in rows:
                o.write("\t".join(map(str, row)) + "\n")
        per_cause = collections.defaultdict(list)
        for row in rows:
            summary[(tool, k, mode, row[0], row[1])] += 1
            per_cause[(row[0], row[1])].append(row)
        for (kind, cause), rs in sorted(per_cause.items()):
            examples.append((tool, k, mode, kind, cause, len(rs), rs[:EXAMPLES], dict(reported_valid)))

    with open(os.path.join(a.out, "summary.tsv"), "w") as o:
        o.write("tool\tk\tmode\tkind\tcause\tcount\n")
        for (tool, k, mode, kind, cause), n in sorted(summary.items()):
            o.write(f"{tool}\t{k}\t{mode}\t{kind}\t{cause}\t{n}\n")

    # Pairwise: valid hits of A's all-hits run missing from B's, per K.
    by_k = collections.defaultdict(dict)
    for R in runs:
        if R["mode"] == "all":
            by_k[R["k"]][R["tool"]] = {h[:4] for h, d in zip(R["hits"], R["mm"]) if 0 <= d <= R["k"]}
    with open(os.path.join(a.out, "pairwise.tsv"), "w") as o:
        o.write("k\ttool_a\ttool_b\ta_not_b\n")
        for k, sets in sorted(by_k.items()):
            for ta, tb in itertools.permutations(sorted(sets), 2):
                o.write(f"{k}\t{ta}\t{tb}\t{len(sets[ta] - sets[tb])}\n")

    # Examples, aligned.
    need = {row[2] for *_, rs, _ in examples for row in rs}
    ex_seqs = {n: "".join(BASES[x] for x in fwd[idx[n]]) for n in need if n in idx}
    with open(os.path.join(a.out, "examples.txt"), "w") as o:
        for tool, k, mode, kind, cause, n, rs, rep in examples:
            o.write(f"=== {tool} K={k} {mode}: {kind} / {cause} ({n})\n")
            for row in rs:
                r = row[2]
                o.write(f"  read {r}  truth hits within K: {row[7]}  reported by {tool}: {len(rep.get(r, []))}\n")
                m = SIM.match(r)
                if m:
                    o.write(f"    true origin {m.group(3)}:{m.group(4)}{m.group(5)} with {m.group(2)} substitutions\n")
                if r in ex_seqs:
                    o.write("\n".join(aligned_view(ex_seqs[r], text, offset, lengths, row[2:6], L)) + "\n")
                near = sorted(rep.get(r, []), key=lambda h: (h[1], abs(h[2] - row[4])))[:3]
                if near:
                    o.write("    the tool's nearest hits: " + ", ".join(f"{c}:{p}{s}" for _, c, p, s in near) + "\n")
            o.write("\n")
    total = sum(summary.values())
    unexplained = sum(n for key, n in summary.items() if key[4] == "unexplained")
    print(f"  {a.out}: {total} discrepancies, {unexplained} unexplained", file=sys.stderr)


def show_read(name, seqs, truth_by_read, runs, text, offset, lengths, L):
    if name not in seqs:
        sys.exit(f"read {name} is not in the read set")
    seq = seqs[name]
    print(f"read {name}\n  {seq}")
    m = SIM.match(name)
    if m:
        print(f"  true origin {m.group(3)}:{m.group(4)}{m.group(5)}, {m.group(2)} substitutions")
    truth = sorted(truth_by_read.get(name, []), key=lambda x: (x[1], x[0]))
    print(f"\ntruth: {len(truth)} hits")
    hits_of = {}
    for R in runs:
        hits_of[(R["tool"], R["k"], R["mode"])] = {h[:4]: int(d) for h, d in zip(R["hits"], R["mm"]) if h[0] == name}
    print("\nper run: reported / truth hits within K (missed, extra)")
    for (tool, k, mode), got in sorted(hits_of.items()):
        tk = {h for h, d in truth if d <= k}
        valid = {h for h, d in got.items() if 0 <= d <= k}
        print(f"  {tool:10s} K={k} {mode:5s} {len(got):6d} / {len(tk):6d}  (missed {len(tk - valid)}, extra {len(set(got) - tk)})")
    print("\ntruth hits (first 20), and the runs that miss each:")
    for h, d in truth[:20]:
        missing = [f"{t}/K{k}" for (t, k, mode), got in sorted(hits_of.items())
                   if mode == "all" and d <= k and not (0 <= got.get(h, -9) <= k)]
        for line in aligned_view(seq, text, offset, lengths, h, L):
            print(line)
        if missing:
            print(f"      missed by: {', '.join(missing)}")
    if len(truth) > 20:
        print(f"  ... {len(truth) - 20} more")


if __name__ == "__main__":
    main()
