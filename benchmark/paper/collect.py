#!/usr/bin/env python3
"""collect.py RUNS_DIR RESULTS_DIR MANIFEST: merges the benchmark's JSON records into TSVs.

  runs.tsv         one line per mapping run (every repeat), with the phase split
  timing.tsv       per series (dataset, tool, k, mode, threads): median and range over the
                   timed repeats (rep >= 1; rep 0 is the warm-up)
  index.tsv        index build time, memory and size
  correctness.tsv  the scores of compare_hits.py, one line per run (origin_found_by_k as JSON)
  manifest.json    a copy of the build manifest
"""
import csv, glob, json, os, shutil, statistics, sys

runs_dir, out_dir, manifest = sys.argv[1], sys.argv[2], sys.argv[3]
os.makedirs(out_dir, exist_ok=True)

RUN_COLS = ["kind", "dataset", "genome", "length", "tool", "k", "mode", "threads", "rep", "status", "reason",
            "wall_s", "user_s", "sys_s", "max_rss_kb", "output_bytes", "align_wall_s", "decode_wall_s",
            "align_max_rss_kb", "decode_max_rss_kb", "tool_index_load", "tool_search", "tool_output_format",
            "tool_total", "finished"]
recs = []
for p in glob.glob(os.path.join(runs_dir, "*", "**", "run.json"), recursive=True):
    r = json.load(open(p))
    if r.get("kind") == "index":
        continue
    recs.append(r)
recs.sort(key=lambda r: (r["kind"], r["dataset"], r["tool"], r["k"], r["mode"], r["threads"], r["rep"]))
with open(os.path.join(out_dir, "runs.tsv"), "w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(RUN_COLS)
    for r in recs:
        w.writerow([r.get(c, "") for c in RUN_COLS])

series = {}
for r in recs:
    key = (r["kind"], r["dataset"], r["genome"], r["length"], r["tool"], r["k"], r["mode"], r["threads"])
    series.setdefault(key, []).append(r)
with open(os.path.join(out_dir, "timing.tsv"), "w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    metrics = ["wall_s", "max_rss_kb", "align_wall_s", "decode_wall_s", "tool_index_load", "tool_search", "tool_output_format"]
    w.writerow(["kind", "dataset", "genome", "length", "tool", "k", "mode", "threads", "status", "n"] +
               [f"{m}_{s}" for m in metrics for s in ("median", "min", "max")])
    for key, rs in sorted(series.items(), key=lambda kv: tuple(map(str, kv[0]))):
        timed = [r for r in rs if r["rep"] >= 1 and r.get("status") == "ok"] or [r for r in rs if r.get("status") == "ok"]
        statuses = {r.get("status") for r in rs if r["rep"] >= 1} or {r.get("status") for r in rs}
        status = "ok" if timed else ("timeout" if "timeout" in statuses else "/".join(sorted(map(str, statuses))))
        row = list(key) + [status, len(timed)]
        for m in metrics:
            v = [r[m] for r in timed if r.get(m) not in (None, "")]
            row += [statistics.median(v), min(v), max(v)] if v else ["", "", ""]
        w.writerow(row)

with open(os.path.join(out_dir, "index.tsv"), "w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["genome", "index", "status", "threads", "wall_s", "max_rss_kb", "index_bytes"])
    for p in sorted(glob.glob(os.path.join(runs_dir, "index", "*", "*", "run.json"))):
        r = json.load(open(p))
        w.writerow([r["genome"], r["index"], r["status"], r["threads"], r["wall_s"], r["max_rss_kb"], r["index_bytes"]])

COR_COLS = ["dataset", "tool", "k", "mode", "truth_source", "reported", "valid", "invalid", "over_k", "over_n", "truth", "missed",
            "extra_not_in_truth", "truth_reads", "reported_reads", "missed_reads", "not_best", "reads_invalid_only", "reads_over_k_only",
            "origin_found_by_k"]
with open(os.path.join(out_dir, "correctness.tsv"), "w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(COR_COLS)
    for p in sorted(glob.glob(os.path.join(runs_dir, "correct", "*.json"))):
        d = os.path.basename(p)[:-5]
        for s in json.load(open(p))["scores"]:
            s["dataset"] = d
            s["origin_found_by_k"] = json.dumps(s.get("origin_found_by_k", ""))
            w.writerow([s.get(c, "") for c in COR_COLS])

if os.path.exists(manifest):
    shutil.copy(manifest, os.path.join(out_dir, "manifest.json"))
print(f"collect: {len(recs)} runs, {len(series)} series -> {out_dir}", file=sys.stderr)
