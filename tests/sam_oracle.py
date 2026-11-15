# sam_oracle.py genome.fa reads.fq out.sam K maxhits [forward_only]
# Checks every SAM record against the genome, and each read's hits against brute force on
# both strands over the index's text (the sequences joined, non-ACGT stored as C): the
# valid hits are windows inside one sequence, without N, at Hamming distance <= K.
#   maxhits 0: the records are exactly the valid hits;
#   maxhits > 0: the primary is at the best valid distance, and when no window at
#     distance <= K is invalid (those still count toward -m in the search), there are
#     min(maxhits, valid hits) records;
#   maxhits -1: the records only.
import sys, re, collections, numpy as np
fa, fq, samf, K, M = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]), int(sys.argv[5])
fwd_only = len(sys.argv) > 6 and sys.argv[6] == "1"
seqs, names, cur = {}, [], None
for l in open(fa):
    if l.startswith(">"):
        cur = l[1:].split()[0]; names.append(cur); seqs[cur] = []
    else: seqs[cur].append(l.strip().upper())
seqs = {k: "".join(v) for k, v in seqs.items()}
raw = "".join(seqs[n] for n in names)
text = np.frombuffer(raw.encode(), dtype=np.uint8).copy()
is_n = ~np.isin(text, np.frombuffer(b"ACGT", dtype=np.uint8))
text[is_n] = ord("C")  # as the index stores them
seq_of = np.concatenate([np.full(len(seqs[n]), i, dtype=np.int32) for i, n in enumerate(names)] or [np.zeros(0, np.int32)])
starts = np.cumsum([0] + [len(seqs[n]) for n in names])
ncs = np.concatenate([[0], np.cumsum(is_n)])
COMP = str.maketrans("ACGT", "TGCA")
def rc(s): return s.translate(COMP)[::-1]
reads = collections.OrderedDict()
with open(fq) as f:
    while True:
        h = f.readline()
        if not h: break
        s = f.readline().strip(); f.readline(); q = f.readline().strip()
        reads[h[1:].split()[0]] = (s, q)
recs = collections.defaultdict(list); bad = 0; msgs = []
def err(m):
    global bad
    bad += 1
    if len(msgs) < 15: msgs.append(m)
for l in open(samf):
    if l.startswith("@"): continue
    f = l.rstrip("\n").split("\t")
    recs[f[0]].append(f)
for rid, (s, q) in reads.items():
    R = recs.get(rid, [])
    if not R: err(f"{rid}: no record"); continue
    L = len(s)
    if R[0][1] == "4":
        if len(R) != 1 or R[0][9] != s: err(f"{rid}: bad unmapped record")
        got = set()
    else:
        got = set()
        for k, f in enumerate(R):
            flag, rname, pos = int(f[1]), f[2], int(f[3]) - 1
            if bool(flag & 0x100) != (k > 0): err(f"{rid}: primary/secondary flags")
            rev = bool(flag & 16)
            want_seq = rc(s.upper()) if rev else s.upper()
            if f[9] != want_seq: err(f"{rid}: SEQ"); continue
            if q and f[10] != (q[::-1] if rev else q): err(f"{rid}: QUAL")
            if f[5] != f"{L}M": err(f"{rid}: CIGAR")
            tags = dict((t.split(":")[0], t.split(":", 2)[2]) for t in f[11:])
            ref = seqs[rname][pos:pos + L]
            d = sum(a != b for a, b in zip(ref, f[9]))
            if len(ref) != L or "N" in ref: err(f"{rid}: window off sequence or on N"); continue
            if int(tags["NM"]) != d or d > K: err(f"{rid}: NM {tags['NM']} distance {d}")
            md, run = "", 0
            for a, b in zip(ref, f[9]):
                if a == b: run += 1
                else: md += f"{run}{a}"; run = 0
            md += str(run)
            if tags["MD"] != md: err(f"{rid}: MD {tags['MD']} want {md}")
            if int(tags["NH"]) != len(R): err(f"{rid}: NH")
            got.add((rname, pos, rev))
        if len(got) != len(R): err(f"{rid}: duplicate records")
    if not set(s.upper()) <= set("ACGT"):
        if got: err(f"{rid}: read with N mapped")
        continue
    if M < 0: continue  # records only
    want, invalid = {}, 0
    n = len(text) - L + 1
    if n > 0:
        valid = (seq_of[:n] == seq_of[L - 1:L - 1 + n]) & ((ncs[L:L + n] - ncs[:n]) == 0)
        for rev in ([False] if fwd_only else [False, True]):
            r = np.frombuffer((rc(s.upper()) if rev else s.upper()).encode(), dtype=np.uint8)
            d = np.zeros(n, dtype=np.int32)
            for i in range(L): d += text[i:i + n] != r[i]
            near = d <= K
            invalid += int((near & ~valid).sum())
            for p in np.flatnonzero(near & valid):
                t = int(seq_of[p]); want[(names[t], int(p - starts[t]), rev)] = int(d[p])
    if M == 0:
        if got != set(want): err(f"{rid}: {len(got)} hits, brute force {len(want)}: missing {sorted(set(want)-got)[:3]} extra {sorted(got-set(want))[:3]}")
    else:
        if len(R) > M: err(f"{rid}: {len(R)} records > -m {M}")
        if not invalid and len(got) != min(M, len(want)): err(f"{rid}: {len(got)} records, want min({M}, {len(want)})")
        if want and not got: err(f"{rid}: unmapped, brute force best {min(want.values())}")
        if got and R[0][1] != "4":
            best = min(want.values()) if want else None
            nm0 = int(dict((t.split(":")[0], t.split(":", 2)[2]) for t in R[0][11:])["NM"])
            if best is None or nm0 != best: err(f"{rid}: primary NM {nm0}, best {best}")
print(f"reads={len(reads)} records={sum(len(v) for v in recs.values())} bad={bad}")
for m in msgs: print("  " + m)
