"""tools.py: how each mapper of the paper benchmark builds an index and maps reads.

Every function returns a list of phases, each (phase name, argv, stdout file or None).
A run executes the phases in order, each under /usr/bin/time; the run's output is a SAM
file that normalize_hits.py reads.

Modes:
  first  one hit per read (BatMis -m 1: a hit with the fewest mismatches)
  all    every hit with at most K mismatches

Exact tools (every hit within K mismatches, Hamming distance):
  batmis1, batmis2  -n K -m 1|0 --sam (-t)
  original          BatMis 3.00: batman -n K -m 1|-mall, then batdecode to SAM (1 thread)
  bwa_aln           made exhaustive: no seed (-l 1000), no gaps (-o 0 -e 0), -k K;
                    all: -N (non-iterative, every hit) -m 20000000, samse -n 1000000 (XA)
  bowtie            -v K (K <= 3); first: -k 1 --best; all: -a
  columba           -m hamming; all: -a all -e K; first: -a best -I <identity of K>
                    (Columba has no one-hit mode: best mode reports every best-stratum hit)
  mrsfast           -e K; first: --best; all: default (all mappings)
bwa_aln_stock runs bwa_aln's commands with BWA built by its own Makefile flags.

Approximate tools (heuristic, gapped; no K):
  bwa_mem           first: default; all: -a
  bowtie2           --end-to-end --very-sensitive; first: default; all: -k 100
"""
import math, os

EXACT = ["batmis1", "batmis2", "original", "bwa_aln", "bowtie", "columba", "mrsfast"]
SUPPLEMENTARY = ["bwa_aln_stock"]
APPROXIMATE = ["bwa_mem", "bowtie2"]

MAX_K = {"bowtie": 3, "original": 5, "bwa_aln": 5, "bwa_aln_stock": 5}
THREADED = {"original": False}

# The index each tool reads, by directory under <index>/<genome>/.
INDEX_OF = {"batmis1": "batmis", "batmis2": "batmis", "original": "batmis", "bwa_aln": "bwa",
            "bwa_aln_stock": "bwa", "bwa_mem": "bwa", "bowtie": "bowtie",
            "bowtie2": "bowtie2", "columba": "columba", "mrsfast": "mrsfast"}
INDEXES = ["batmis", "bwa", "bowtie", "bowtie2", "columba", "mrsfast"]


def accepts(tool, k):
    return k <= MAX_K.get(tool, 7)


def build_index(index, tools, fasta, out_dir, genome, threads):
    """Phases that build `index` for genome `genome` in out_dir (created by the caller)."""
    b = lambda tool, exe: os.path.join(tools, "bin", tool, exe)
    p = os.path.join(out_dir, genome)
    if index == "batmis":
        return [("build", [b("batmis", "bwt_build_index"), "-o", out_dir, fasta], None)]
    if index == "bwa":
        return [("build", [b("bwa", "bwa"), "index", "-p", p, fasta], None)]
    if index == "bowtie":
        return [("build", [b("bowtie", "bowtie-build"), "--threads", str(threads), "-q", fasta, p], None)]
    if index == "bowtie2":
        return [("build", [b("bowtie2", "bowtie2-build"), "--threads", str(threads), "-q", fasta, p], None)]
    if index == "columba":
        return [("build", [b("columba", "columba_build"), "-r", p, "-f", fasta], None)]
    if index == "mrsfast":
        # mrsFAST writes <fasta>.index beside the FASTA it is given: give it a copy here.
        return [("copy", ["cp", fasta, p + ".fa"], None),
                ("build", [b("mrsfast", "mrsfast"), "--index", p + ".fa"], None)]
    raise ValueError(index)


def map_reads(tool, tools, idx_dir, genome, reads, length, k, mode, threads, out, extra=()):
    """Phases that map `reads` with `tool`; the SAM output goes to `out`."""
    b = lambda t, exe: os.path.join(tools, "bin", t, exe)
    p = os.path.join(idx_dir, genome)
    first = mode == "first"
    t = str(threads)
    if tool in ("batmis1", "batmis2"):
        return [("map", [b("batmis", tool), "-g", os.path.join(idx_dir, os.path.basename(genome_fasta_name(genome))),
                         "-q", reads, "-n", str(k), "-m", "1" if first else "0", "--sam", "-t", t, "-o", out,
                         "--timing", *extra], None)]
    if tool == "original":
        g = os.path.join(idx_dir, os.path.basename(genome_fasta_name(genome)))
        return [("align", [b("original", "batman"), "-g", g, "-q", reads, "-n", str(k),
                           *(["-m", "1"] if first else ["-mall"]), "-o", out + ".bin"], None),
                ("decode", [b("original", "batdecode"), "-i", out + ".bin", "-g", g, "-o", out], None)]
    if tool.startswith("bwa_aln"):
        build = {"bwa_aln": "bwa", "bwa_aln_stock": "bwa_stock"}[tool]
        aln = [b(build, "bwa"), "aln", "-t", t, "-n", str(k), "-o", "0", "-e", "0", "-l", "1000", "-k", str(k)]
        if not first:
            aln += ["-N", "-m", "20000000"]
        samse = [b(build, "bwa"), "samse"] + ([] if first else ["-n", "1000000"]) + [p, out + ".sai", reads]
        return [("align", aln + [p, reads], out + ".sai"), ("decode", samse, out)]
    if tool == "bwa_mem":
        return [("map", [b("bwa", "bwa"), "mem", "-t", t, *([] if first else ["-a"]), p, reads], out)]
    if tool == "bowtie":
        return [("map", [b("bowtie", "bowtie"), "-p", t, "-v", str(k), *(["-k", "1", "--best"] if first else ["-a"]),
                         "--sam", "-x", p, reads, out], None)]
    if tool == "bowtie2":
        return [("map", [b("bowtie2", "bowtie2"), "-p", t, "--end-to-end", "--very-sensitive",
                         *([] if first else ["-k", "100"]), "-x", p, "-U", reads, "-S", out], None)]
    if tool == "columba":
        if first:
            # The highest whole identity that still admits k mismatches; below 100 bp it may
            # also admit k + 1, whose hits the comparison ignores.
            ident = math.floor(100.0 * (length - k) / length)
            mode_args = ["-a", "best", "-I", str(ident)]
        else:
            mode_args = ["-a", "all", "-e", str(k)]
        return [("map", [b("columba", "columba"), "-r", p, "-f", reads, "-m", "hamming", *mode_args, "-t", t,
                         "-o", out, "-l", out + ".log"], None)]
    if tool == "mrsfast":
        return [("map", [b("mrsfast", "mrsfast"), "--search", p + ".fa", "--seq", reads, "-e", str(k),
                         "--threads", t, *(["--best"] if first else []), "-o", out, "--disable-nohits"], None)]
    raise ValueError(tool)


def genome_fasta_name(genome):
    return f"{genome}.fa"
