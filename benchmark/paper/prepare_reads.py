#!/usr/bin/env python3
"""prepare_reads.py: the same fixed-length, ACGT-only read sets for every tool.

  prepare_reads.py <in.fastq[.gz] | URL> <out_prefix> --lengths 40,100,150
                   [--fraction 1.0] [--max-reads N]

BatMis maps reads of one length with only A, C, G, T, so every tool is given that input:
each read is kept with probability --fraction (decided by a hash of its name, so the
sample is the same on every run), trimmed to its first L bases for each length L, and
dropped from a length's set if it is shorter than L or holds another character there.
Writes <out_prefix>.<L>.fq per length and <out_prefix>.json with the counts. A URL is
streamed (curl) and not stored.
"""
import argparse, gzip, hashlib, json, subprocess, sys

ACGT = set(b"ACGTacgt")


def open_input(src):
    if src.startswith(("http://", "https://", "ftp://")):
        p = subprocess.Popen(["curl", "-sSfL", "--retry", "5", src], stdout=subprocess.PIPE)
        return gzip.open(p.stdout) if src.endswith(".gz") else p.stdout
    return gzip.open(src) if src.endswith(".gz") else open(src, "rb")


def keep(name, fraction):
    if fraction >= 1.0:
        return True
    h = int.from_bytes(hashlib.blake2b(name, digest_size=8).digest(), "little")
    return h < fraction * 2**64


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("src")
    ap.add_argument("out_prefix")
    ap.add_argument("--lengths", required=True)
    ap.add_argument("--fraction", type=float, default=1.0)
    ap.add_argument("--max-reads", type=int, default=0, help="stop after this many sampled reads (0: all)")
    a = ap.parse_args()
    lengths = [int(x) for x in a.lengths.split(",")]
    outs = {L: open(f"{a.out_prefix}.{L}.fq", "wb") for L in lengths}
    counts = {"input": 0, "sampled": 0, **{f"kept_{L}": 0 for L in lengths},
              **{f"short_{L}": 0 for L in lengths}, **{f"non_acgt_{L}": 0 for L in lengths}}
    f = open_input(a.src)
    while True:
        h = f.readline()
        if not h:
            break
        s = f.readline().rstrip(b"\r\n")
        f.readline()
        q = f.readline().rstrip(b"\r\n")
        counts["input"] += 1
        name = h[1:].split()[0]
        if not keep(name, a.fraction):
            continue
        counts["sampled"] += 1
        for L in lengths:
            if len(s) < L:
                counts[f"short_{L}"] += 1
            elif not ACGT.issuperset(s[:L]):
                counts[f"non_acgt_{L}"] += 1
            else:
                outs[L].write(b"@" + name + b"\n" + s[:L].upper() + b"\n+\n" + q[:L] + b"\n")
                counts[f"kept_{L}"] += 1
        if a.max_reads and counts["sampled"] >= a.max_reads:
            break
    for o in outs.values():
        o.close()
    counts["source"] = a.src
    counts["fraction"] = a.fraction
    json.dump(counts, open(f"{a.out_prefix}.json", "w"), indent=1)
    print(json.dumps(counts), file=sys.stderr)


if __name__ == "__main__":
    main()
