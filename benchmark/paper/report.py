#!/usr/bin/env python3
"""report.py RESULTS_DIR: the paper's tables (report.md) and figures (figures/*.pdf, .png)
from collect.py's TSVs.

Figures (one panel per read length, one line per tool; median of the timed repeats, the
range as a band):
  time_vs_k_<dataset group>_<mode>_t<threads>   search wall time against K
  memory_<genome>                               peak memory per tool, K = 2, all hits
  scaling_<dataset>_k<K>_<mode>                 wall time against threads
Tables: correctness per read set (missed / invalid / not-best hits), index build, and the
timing summary with the align/decode split.

Colours follow the tool, never its position: each tool keeps one colour and one marker in
every figure (the reference categorical palette, fixed order; markers as a second cue).
"""
import csv, glob, json, os, sys
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

OUT = sys.argv[1]
FIG = os.path.join(OUT, "figures")
os.makedirs(FIG, exist_ok=True)

# Fixed tool -> (colour, marker, label). Exact tools in the first slots; the two BatMis
# versions take the strongest pair.
STYLE = {
    "batmis2": ("#2a78d6", "o", "BatMis2 (bidirectional)"),
    "batmis1": ("#eb6834", "s", "BatMis1"),
    "original": ("#1baf7a", "^", "BatMis 3.00"),
    "bwa_aln": ("#eda100", "D", "BWA aln (exhaustive)"),
    "bowtie": ("#e87ba4", "v", "Bowtie"),
    "columba": ("#008300", "P", "Columba"),
    "mrsfast": ("#6250d6", "X", "mrsFAST-Ultra"),
    "bwa_mem": ("#e34948", "*", "BWA-MEM"),
    "bowtie2": ("#7a7a7a", "h", "Bowtie 2"),
    "bwa_aln_stock": ("#eda100", "1", "BWA aln, stock build"),
}
ORDER = list(STYLE)
MAIN = [t for t in ORDER if not t.startswith("bwa_aln_")]  # the supplementary BWA builds get their own figure
INK, MUTED, GRID = "#1f1f1f", "#5c5c5c", "#e4e4e4"

plt.rcParams.update({
    "font.size": 8, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED, "ytick.color": MUTED,
    "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.color": GRID,
    "grid.linewidth": 0.6, "lines.linewidth": 2, "lines.markersize": 5, "legend.frameon": False,
    "figure.dpi": 150, "savefig.bbox": "tight",
})


def read_tsv(name):
    p = os.path.join(OUT, name)
    if not os.path.exists(p):
        return []
    return list(csv.DictReader(open(p), delimiter="\t"))


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def plain_log(axis):
    """Log axis with plain numbers (0.2, 1, 20) instead of 2x10^-1."""
    from matplotlib.ticker import FuncFormatter, LogLocator
    axis.set_major_locator(LogLocator(base=10, subs=(1.0, 2.0, 5.0)))
    axis.set_minor_formatter(FuncFormatter(lambda v, _: ""))
    axis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))


def save(fig, name):
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(FIG, f"{name}.{ext}"))
    plt.close(fig)


timing = read_tsv("timing.tsv")
for r in timing:
    r["k"] = int(r["k"]); r["threads"] = int(r["threads"]); r["length"] = int(r["length"])

# ---------------------------------------------------------------- time against K
groups = defaultdict(list)
for r in timing:
    if r["kind"] in ("timing", "realdata") and r["k"] >= 0 and r["status"] == "ok":
        base = r["dataset"].rsplit("_", 2)[0] if r["kind"] == "timing" else r["dataset"].rsplit("_", 1)[0]
        groups[(base, r["mode"], r["threads"])].append(r)
for (base, mode, threads), rs in sorted(groups.items()):
    lengths = sorted({r["length"] for r in rs})
    fig, axes = plt.subplots(1, len(lengths), figsize=(2.6 * len(lengths), 2.6), sharey=False, squeeze=False)
    for ax, L in zip(axes[0], lengths):
        for tool in MAIN:
            pts = sorted((r["k"], num(r["wall_s_median"]), num(r["wall_s_min"]), num(r["wall_s_max"]))
                         for r in rs if r["tool"] == tool and r["length"] == L)
            if not pts:
                continue
            c, m, lab = STYLE[tool]
            ks = [p[0] for p in pts]
            ax.plot(ks, [p[1] for p in pts], color=c, marker=m, label=lab)
            ax.fill_between(ks, [p[2] for p in pts], [p[3] for p in pts], color=c, alpha=0.15, linewidth=0)
        ax.set_yscale("log")
        plain_log(ax.yaxis)
        ax.set_title(f"{L} bp", color=INK)
        ax.set_xlabel("mismatches K")
    axes[0][0].set_ylabel("wall time (s)")
    h, l = axes[0][0].get_legend_handles_labels()
    fig.legend(h, l, loc="upper center", ncol=4, bbox_to_anchor=(0.5, 1.12))
    fig.suptitle(f"{base}, {'first hit' if mode == 'first' else 'all hits'}, {threads} thread{'s' * (threads > 1)}",
                 y=1.18, color=INK)
    save(fig, f"time_vs_k_{base}_{mode}_t{threads}")

# ---------------------------------------------------------------- BWA builds (supplementary)
bwa = defaultdict(list)
for r in timing:
    if r["kind"] == "timing" and r["tool"].startswith("bwa_aln") and r["status"] == "ok":
        bwa[(r["dataset"], r["mode"], r["threads"])].append(r)
for (dataset, mode, threads), rs in sorted(bwa.items()):
    if len({r["tool"] for r in rs}) < 2:
        continue
    fig, ax = plt.subplots(figsize=(3.2, 2.4))
    for tool, color, marker in (("bwa_aln_stock", "#7a7a7a", "s"), ("bwa_aln", "#eda100", "D")):
        pts = sorted((r["k"], num(r["wall_s_median"])) for r in rs if r["tool"] == tool)
        if pts:
            ax.plot([p[0] for p in pts], [p[1] for p in pts], color=color, marker=marker, label=STYLE[tool][2])
    ax.set_yscale("log"); plain_log(ax.yaxis)
    ax.set_xlabel("mismatches K"); ax.set_ylabel("wall time (s)")
    ax.set_title(f"{dataset}, {'first hit' if mode == 'first' else 'all hits'}, {threads} thread{'s' * (threads > 1)}",
                 color=INK, fontsize=8)
    ax.legend(fontsize=6)
    save(fig, f"bwa_builds_{dataset}_{mode}_t{threads}")

# ---------------------------------------------------------------- memory
mem = defaultdict(dict)
for r in timing:
    if r["kind"] in ("timing", "realdata") and r["k"] in (2, -1) and r["mode"] == "all" and r["status"] == "ok":
        mem[(r["genome"], r["dataset"], r["threads"])][r["tool"]] = num(r["max_rss_kb_median"]) / 1024
for (genome, dataset, threads), by_tool in sorted(mem.items()):
    tools = [t for t in MAIN if t in by_tool]
    fig, ax = plt.subplots(figsize=(3.6, 0.28 * len(tools) + 0.6))
    ax.barh(range(len(tools)), [by_tool[t] for t in tools], color=[STYLE[t][0] for t in tools], height=0.7)
    ax.set_yticks(range(len(tools)), [STYLE[t][2] for t in tools])
    ax.invert_yaxis()
    ax.set_xlabel("peak memory (MB)")
    ax.grid(axis="y", visible=False)
    for i, t in enumerate(tools):
        ax.text(by_tool[t], i, f" {by_tool[t]:,.0f}", va="center", color=INK, fontsize=7)
    ax.set_title(f"{dataset}, K = 2 (approximate tools: default), all hits, {threads} threads", color=INK, fontsize=8)
    save(fig, f"memory_{dataset}_t{threads}")

# ---------------------------------------------------------------- thread scaling
scal = defaultdict(list)
for r in timing:
    if r["kind"] in ("scaling", "timing", "realdata") and r["status"] == "ok":
        scal[(r["dataset"], r["k"], r["mode"])].append(r)
for (dataset, k, mode), rs in sorted(scal.items()):
    if not any(r["kind"] == "scaling" for r in rs):
        continue
    fig, ax = plt.subplots(figsize=(3.2, 2.6))
    for tool in MAIN:
        pts = sorted((r["threads"], num(r["wall_s_median"])) for r in rs if r["tool"] == tool)
        if len(pts) < 2:
            continue
        c, m, lab = STYLE[tool]
        ax.plot([p[0] for p in pts], [p[1] for p in pts], color=c, marker=m, label=lab)
    ax.set_xscale("log", base=2); ax.set_yscale("log")
    plain_log(ax.yaxis)
    ax.set_xlabel("threads"); ax.set_ylabel("wall time (s)")
    ax.set_title(f"{dataset}, K = {k}, {'first hit' if mode == 'first' else 'all hits'}", color=INK)
    ax.legend(fontsize=6, loc="upper right")
    save(fig, f"scaling_{dataset}_k{k}_{mode}")

# ---------------------------------------------------------------- tables
lines = ["# Paper benchmark results", ""]
man = os.path.join(OUT, "manifest.json")
if os.path.exists(man):
    m = json.load(open(man))
    lines += [f"Machine: {m['cpu']}, {m['cores']} CPUs, {m['memory']}; {m['compiler']}; flags `{m['opt_flags']}`.", "",
              "| binary | flags | popcnt | prefetch |", "|---|---|---:|---:|"]
    lines += [f"| {k} | `{v.get('flags', m['opt_flags'])}` | {v['popcnt']} | {v['prefetch']} |" for k, v in sorted(m["binaries"].items())
              if k.split("/")[1] in ("batmis1", "batmis2", "batman", "batdecode", "bwa", "bowtie-align-s", "bowtie2-align-s",
                                     "columba", "mrsfast")]
    lines.append("")

idx = read_tsv("index.tsv")
if idx:
    lines += ["## Index build", "", "| genome | index | wall (s) | peak memory (MB) | size (MB) |", "|---|---|---:|---:|---:|"]
    for r in idx:
        lines.append(f"| {r['genome']} | {r['index']} | {num(r['wall_s']):.1f} | {num(r['max_rss_kb']) / 1024:,.0f} | "
                     f"{num(r['index_bytes']) / 2**20:,.0f} |")
    lines.append("")

cor = read_tsv("correctness.tsv")
if cor:
    lines += ["## Correctness", "",
              "All hits: hits missed and invalid hits reported against the truth (brute force, or the verified union of the "
              "exact tools). First hit: reads with a hit that got none, and reads whose hit is not at the read's best "
              "mismatch count.", "",
              "| read set | tool | K | mode | truth | missed | invalid | missed reads | not best | true origin (k <= K) |",
              "|---|---|---:|---|---:|---:|---:|---:|---:|---:|"]
    for r in cor:
        o = json.loads(r["origin_found_by_k"]) if r["origin_found_by_k"] not in ("", '""') else {}
        k = int(r["k"])
        f = sum(v[0] for kk, v in o.items() if k < 0 or int(kk) <= k)
        t = sum(v[1] for kk, v in o.items() if k < 0 or int(kk) <= k)
        lines.append(f"| {r['dataset']} | {r['tool']} | {k} | {r['mode']} | {r['truth'] or r['truth_reads']} | "
                     f"{r['missed']} | {r['invalid']} | {r['missed_reads']} | {r['not_best']} | "
                     + (f"{100 * f / t:.2f}%" if t else "") + " |")
    lines.append("")

disc = sorted(glob.glob(os.path.join(OUT, "discrepancies", "*", "summary.tsv")))
if disc:
    lines += ["## Discrepancies", "",
              "Every missed or extra hit of the exact tools, by cause (`discrepancies.py`; examples in "
              "`discrepancies/<read set>/examples.txt`). Read sets without any are left out.", "",
              "| read set | tool | K | mode | kind | cause | count |", "|---|---|---:|---|---|---|---:|"]
    for p in disc:
        name = os.path.basename(os.path.dirname(p))
        for r in csv.DictReader(open(p), delimiter="\t"):
            lines.append(f"| {name} | {r['tool']} | {r['k']} | {r['mode']} | {r['kind']} | {r['cause']} | {r['count']} |")
    lines.append("")

if timing:
    lines += ["## Timing", "", "Median wall seconds over the timed repeats (min-max); align / decode where the tool "
              "separates them; peak memory.", "",
              "| read set | tool | K | mode | threads | wall (s) | align (s) | decode (s) | memory (MB) | status |",
              "|---|---|---:|---|---:|---:|---:|---:|---:|---|"]
    for r in sorted(timing, key=lambda r: (r["kind"], r["dataset"], r["mode"], r["threads"], r["k"], ORDER.index(r["tool"]) if r["tool"] in ORDER else 99)):
        if r["kind"] not in ("timing", "realdata", "load"):
            continue
        w = num(r["wall_s_median"])
        lines.append(f"| {r['dataset']} | {r['tool']} | {r['k']} | {r['mode']} | {r['threads']} | "
                     + (f"{w:.2f} ({num(r['wall_s_min']):.2f}-{num(r['wall_s_max']):.2f})" if w is not None else "")
                     + f" | {r['align_wall_s_median'] and f'{num(r['align_wall_s_median']):.2f}'} | "
                     f"{r['decode_wall_s_median'] and f'{num(r['decode_wall_s_median']):.2f}'} | "
                     + (f"{num(r['max_rss_kb_median']) / 1024:,.0f}" if r['max_rss_kb_median'] else "") + f" | {r['status']} |")
open(os.path.join(OUT, "report.md"), "w").write("\n".join(lines) + "\n")
print(f"report: {OUT}/report.md, {len(os.listdir(FIG))} figure files", file=sys.stderr)
