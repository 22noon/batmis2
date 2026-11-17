#!/usr/bin/env python3
"""normalize_hits.py: one tool's SAM output as a sorted, de-duplicated hit list.

  normalize_hits.py in.sam out.tsv

Every tool in the benchmark writes SAM (original BatMis through batdecode). Each mapped
record, primary or secondary, is a hit; supplementary records (chimeric parts) are not.
BWA's extra hits in the XA tag (chr,+-pos,CIGAR,NM;...) are hits too.

Output columns, tab-separated, sorted, one line per distinct hit:
  read  contig  pos0  strand  nm  exact
pos0 is the 0-based leftmost reference position of the whole read: the record's POS less
the bases soft- or hard-clipped at its left end, so a clipped hit still names the read's
window. nm is the NM tag (-1 if absent). exact is 1 if the CIGAR is a single <L>M or <L>=/X
mix (no clip, no indel), else 0: only exact hits count in the Hamming comparisons.
"""
import re, subprocess, sys

CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")


def parse_cigar(cigar):
    ops = CIGAR.findall(cigar)
    lead = 0
    for n, op in ops:
        if op in "SH":
            lead += int(n)
        else:
            break
    exact = all(op in "M=X" for _, op in ops)
    return lead, int(exact)


def main():
    src, dst = sys.argv[1], sys.argv[2]
    hits = set()
    with open(src, "rb") as f:
        for raw in f:
            if raw.startswith(b"@"):
                continue
            t = raw.rstrip(b"\n").split(b"\t")
            flag = int(t[1])
            if flag & 4 or flag & 2048:
                continue
            name = t[0].decode()
            nm = -1
            xa = None
            for tag in t[11:]:
                if tag.startswith(b"NM:i:"):
                    nm = int(tag[5:])
                elif tag.startswith(b"XA:Z:"):
                    xa = tag[5:].decode()
            lead, exact = parse_cigar(t[5].decode())
            hits.add((name, t[2].decode(), int(t[3]) - 1 - lead, "-" if flag & 16 else "+", nm, exact))
            if xa:
                for alt in xa.rstrip(";").split(";"):
                    chrom, spos, cig, anm = alt.split(",")
                    lead, exact = parse_cigar(cig)
                    hits.add((name, chrom, abs(int(spos)) - 1 - lead, spos[0], int(anm), exact))
    with open(dst + ".tmp", "w") as o:
        for h in hits:
            o.write("\t".join(map(str, h)) + "\n")
    # External sort: the hit lists of human all-hit runs do not fit comfortably in memory twice.
    subprocess.run(["sort", "-u", "-S", "25%", "-o", dst, dst + ".tmp"], check=True,
                   env={"LC_ALL": "C", "PATH": "/usr/bin:/bin"})
    subprocess.run(["rm", "-f", dst + ".tmp"])


if __name__ == "__main__":
    main()
