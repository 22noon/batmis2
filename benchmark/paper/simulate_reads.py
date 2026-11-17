#!/usr/bin/env python3
"""simulate_reads.py: single-end reads with a known origin, for the paper benchmark.

  simulate_reads.py genome.fa out.fq --length L --reads N --mode stratified|realistic
                    [--max-k 7] [--random-fraction 0.1] [--indel-rate 0] [--seed S]

Origins are uniform over the windows of length L (L + 8 with indels) that lie inside one
sequence and hold only ACGT; each read is taken from the forward or the reverse strand with
equal chance.

stratified: read i gets exactly k = i mod (max_k + 1) substitutions, at distinct positions,
    each to one of the three other bases (uniformly). It tests sensitivity per mismatch count.
realistic: each base is substituted with a probability that rises linearly along the read
    from 0.1% to 1.5% (an Illumina-like profile); with --indel-rate r, each base is also
    deleted, or followed by an inserted base, with probability r/2 each (approximate group).

A --random-fraction of the reads are uniform random sequences, which match nowhere at
small k: they measure the cost of a read without hits.

Read names: sim<i>_k<k>_<contig>_<pos0>_<strand>, where pos0 is the 0-based leftmost
reference position of the read's window and strand is + or -; random reads are rand<i>.
For realistic reads k is the number of substitutions (indels are not counted).
Qualities are all 'I'.
"""
import argparse, sys
import numpy as np

CODE = np.full(256, 4, dtype=np.uint8)
for i, c in enumerate(b"ACGT"):
    CODE[c] = i
    CODE[c + 32] = i
BASES = np.frombuffer(b"ACGT", dtype=np.uint8)


def read_fasta(path):
    names, seqs, cur = [], [], []
    with open(path, "rb") as f:
        for line in f:
            if line.startswith(b">"):
                if names:
                    seqs.append(b"".join(cur))
                names.append(line[1:].split()[0].decode())
                cur = []
            else:
                cur.append(line.strip())
    if names:
        seqs.append(b"".join(cur))
    return names, seqs


def valid_starts(seq_codes, window):
    """Start positions whose window of the given length is all ACGT."""
    bad = (seq_codes > 3).astype(np.int64)
    if len(seq_codes) < window:
        return np.zeros(0, dtype=np.int64)
    c = np.concatenate([[0], np.cumsum(bad)])
    n = len(seq_codes) - window + 1
    return np.nonzero(c[window:window + n] - c[:n] == 0)[0]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("genome")
    ap.add_argument("out")
    ap.add_argument("--length", type=int, required=True)
    ap.add_argument("--reads", type=int, required=True)
    ap.add_argument("--mode", choices=["stratified", "realistic"], required=True)
    ap.add_argument("--max-k", type=int, default=7)
    ap.add_argument("--random-fraction", type=float, default=0.1)
    ap.add_argument("--indel-rate", type=float, default=0.0)
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    L = a.length
    window = L + (8 if a.indel_rate > 0 else 0)

    names, seqs = read_fasta(a.genome)
    codes = [CODE[np.frombuffer(s, dtype=np.uint8)] for s in seqs]
    starts = [valid_starts(c, window) for c in codes]
    weights = np.array([len(s) for s in starts], dtype=np.float64)
    if weights.sum() == 0:
        sys.exit("no valid window in the genome")
    weights /= weights.sum()

    n_random = int(round(a.reads * a.random_fraction))
    n_sim = a.reads - n_random
    contig = rng.choice(len(names), size=n_sim, p=weights)
    comp = np.array([3, 2, 1, 0], dtype=np.uint8)
    qual = b"I" * L
    sub_prob = np.linspace(0.001, 0.015, L)
    out = open(a.out, "wb")
    for i in range(n_sim):
        ci = contig[i]
        pos = int(starts[ci][rng.integers(len(starts[ci]))])
        win = codes[ci][pos:pos + window].copy()
        strand = "+" if rng.random() < 0.5 else "-"
        if strand == "-":
            win = comp[win[::-1]]
        if a.mode == "stratified":
            r = win[:L].copy()
            k = i % (a.max_k + 1)
            where = rng.choice(L, size=k, replace=False)
            r[where] = (r[where] + rng.integers(1, 4, size=k)) % 4
        else:
            out_bases, k, j = [], 0, 0
            while len(out_bases) < L and j < window:
                if a.indel_rate > 0:
                    u = rng.random()
                    if u < a.indel_rate / 2:  # deletion: skip a reference base
                        j += 1
                        continue
                    if u < a.indel_rate:  # insertion before this base
                        out_bases.append(int(rng.integers(4)))
                        continue
                b = int(win[j])
                if rng.random() < sub_prob[len(out_bases)]:
                    b = (b + int(rng.integers(1, 4))) % 4
                    k += 1
                out_bases.append(b)
                j += 1
            r = np.array(out_bases[:L], dtype=np.uint8)
            if len(r) < L:  # ran out of window (very high indel rate): pad from the reference
                r = np.concatenate([r, win[len(r):L]])
        name = f"sim{i}_k{k}_{names[ci]}_{pos}_{strand}".encode()
        out.write(b"@" + name + b"\n" + BASES[r].tobytes() + b"\n+\n" + qual + b"\n")
    for i in range(n_random):
        r = rng.integers(4, size=L).astype(np.uint8)
        out.write(b"@rand%d\n" % i + BASES[r].tobytes() + b"\n+\n" + qual + b"\n")
    out.close()


if __name__ == "__main__":
    main()
