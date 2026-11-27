#!/usr/bin/env python3
"""bench.py: the paper benchmark's driver (README.md in this directory).

  bench.py [-c config.yaml] fetch                 genomes into <work>/genomes
  bench.py [-c ...] simulate [--genomes g,..]     simulated read sets into <work>/reads
  bench.py [-c ...] prepare [--real name,..]      real read sets into <work>/reads
  bench.py [-c ...] index [--genomes g,..] [--indexes i,..]   build and time every index
  bench.py [-c ...] calibrate                     time 5,000 reads per tool (K = 0 and 5) on the
                                                  calibrated genomes, and size their timing runs
  bench.py [-c ...] estimate                      node time of the sized runs still to do
  bench.py [-c ...] list [filters]                the runs of the matrix, and which are done
  bench.py [-c ...] status [filters]              done/timeout/error per tool, the newest finished
                                                  run, and the runs in progress with their age
  bench.py [-c ...] run [filters] [--dry-run]     run the matrix's missing runs, in order
  bench.py [-c ...] truth [--datasets d,..]       brute-force hits of the correctness sets
  bench.py [-c ...] correct [filters]             score the correctness runs
  bench.py [-c ...] discrepancies [filters] [--read NAME]   explain each missed or extra hit
  bench.py [-c ...] collect                       merge every result into <results>/*.tsv

Filters select runs by field: --set kind=timing,genome=ecoli,tool=batmis2,k=2 (values may be
lists split by '|', e.g. tool=batmis1|batmis2). Every run writes <work>/runs/<run id>/run.json;
a run whose run.json exists is skipped, so `run` can be stopped and started again.
"""
import argparse, glob, gzip, json, os, re, shutil, signal, subprocess, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import tools as T  # noqa: E402

# ---------------------------------------------------------------- configuration


def load_config(path):
    import yaml
    text = open(path).read()
    text = re.sub(r"\$\{(\w+)\}", lambda m: os.environ.get(m.group(1), m.group(0)), text)
    cfg = yaml.safe_load(text)
    if "${" in str(cfg["results"]):  # BENCH_RESULTS unset: the default of the job scripts
        cfg["results"] = os.path.join(HERE, "results")
    for key in ("work", "tools", "results"):
        if "${" in str(cfg[key]):
            sys.exit(f"config: {key} uses an unset environment variable: {cfg[key]}")
    if "${" in str(cfg.get("local_output") or ""):  # e.g. $USER unset in a batch job
        import tempfile
        cfg["local_output"] = os.path.join(tempfile.gettempdir(), f"batmis_bench_out_{os.getuid()}")
    return cfg


def wpath(cfg, *parts):
    return os.path.join(cfg["work"], *parts)


def genome_fa(cfg, g):
    return wpath(cfg, "genomes", f"{g}.fa")


def index_dir(cfg, g, index):
    return wpath(cfg, "index", g, index)


def count_reads(fq):
    """Reads in a FASTQ file; cached beside it (<fq>.nreads), since sets reach 10M reads."""
    cache = fq + ".nreads"
    if os.path.exists(cache) and os.path.getmtime(cache) >= os.path.getmtime(fq):
        return int(open(cache).read())
    n = 0
    with open(fq, "rb") as f:
        for n, _ in enumerate(f, 1):
            pass
    try:
        open(cache, "w").write(str(n // 4))
    except OSError:
        pass
    return n // 4


def subset_reads(fq, n):
    """A set of n reads drawn from fq at random (fixed seed, file order kept), so that it keeps
    the set's mix. Not a stride: the simulator cycles k = 0..max_k read by read, and every
    10th read of an 8-cycle is only odd k. Written once as <fq without .fq>.n<n>.fq; fq itself
    when n covers it."""
    import numpy as np
    import zlib
    total = count_reads(fq)
    if n >= total:
        return fq
    out = f"{fq[:-3]}.n{n}.fq"
    if os.path.exists(out):
        return out
    seed = zlib.crc32(f"{os.path.basename(fq)}:{n}".encode())
    keep = np.zeros(total, dtype=bool)
    keep[np.random.default_rng(seed).choice(total, size=n, replace=False)] = True
    with open(fq, "rb") as f, open(out + ".tmp", "wb") as o:
        for i in range(total):
            rec = [f.readline() for _ in range(4)]
            if keep[i]:
                o.writelines(rec)
    os.rename(out + ".tmp", out)
    return out


def calibration(cfg):
    return cfg.get("calibration") or {}


def calibrated(cfg, genome):
    return genome in calibration(cfg).get("genomes", [])


def correctness_reads(cfg, genome):
    s = cfg["simulated"]
    return s.get("correctness_reads_by_genome", {}).get(genome, s["correctness_reads"])


def correctness_fq(cfg, ds):
    """The reads of a correctness set (check, real, indel) as every run and score uses them:
    check sets at the configured size, cut by stride from a larger existing file; the realistic
    sets cut to correctness_reads_by_genome where a genome has one (human), else whole."""
    if ds["kind"] == "check":
        return subset_reads(ds["fq"], correctness_reads(cfg, ds["genome"]))
    by = cfg["simulated"].get("correctness_reads_by_genome", {})
    return subset_reads(ds["fq"], by[ds["genome"]]) if ds["genome"] in by else ds["fq"]


# ---------------------------------------------------------------- fetch, simulate, prepare


def cmd_fetch(cfg, a):
    os.makedirs(wpath(cfg, "genomes"), exist_ok=True)
    for g, spec in cfg["genomes"].items():
        out = genome_fa(cfg, g)
        if os.path.exists(out):
            continue
        if "url" in spec:
            print(f"fetch {g}: {spec['url']}", file=sys.stderr)
            tmp = out + ".download"
            subprocess.run(["curl", "-sSfL", "--retry", "5", "-o", tmp, spec["url"]], check=True)
            if spec["url"].endswith(".gz"):
                with gzip.open(tmp, "rb") as i, open(out + ".tmp", "wb") as o:
                    shutil.copyfileobj(i, o, 1 << 24)
                os.remove(tmp)
            else:
                os.rename(tmp, out + ".tmp")
            os.rename(out + ".tmp", out)
        else:  # contigs of another genome
            src = genome_fa(cfg, spec["from"])
            want = set(spec["contigs"])
            with open(src, "rb") as i, open(out + ".tmp", "wb") as o:
                on = False
                for line in i:
                    if line.startswith(b">"):
                        on = line[1:].split()[0].decode() in want
                    if on:
                        o.write(line)
            os.rename(out + ".tmp", out)


def sim_sets(cfg):
    """(name, genome, length, mode, reads, extra args) of every simulated read set."""
    s = cfg["simulated"]
    out = []
    for g in s["genomes"]:
        for L in s.get("lengths_by_genome", {}).get(g, s["lengths"]):
            out.append((f"sim_{g}_{L}_timing", g, L, "stratified", s["timing_reads"], []))
            out.append((f"sim_{g}_{L}_check", g, L, "stratified", correctness_reads(cfg, g), []))
            out.append((f"sim_{g}_{L}_real", g, L, "realistic", s["realistic_reads"], []))
            out.append((f"sim_{g}_{L}_indel", g, L, "realistic", s["realistic_reads"], ["--indel-rate", "0.001"]))
    for g in s.get("pilot_genomes", []):
        for L in s["lengths"]:
            out.append((f"sim_{g}_{L}_check", g, L, "stratified", s.get("pilot_reads", 10000), []))
    return out


def cmd_simulate(cfg, a):
    os.makedirs(wpath(cfg, "reads"), exist_ok=True)
    s = cfg["simulated"]
    want = set(a.genomes.split(",")) if a.genomes else None
    for i, (name, g, L, mode, n, extra) in enumerate(sim_sets(cfg)):
        out = wpath(cfg, "reads", name + ".fq")
        if (want and g not in want) or os.path.exists(out):
            continue
        print(f"simulate {name}", file=sys.stderr)
        subprocess.run([sys.executable, os.path.join(HERE, "simulate_reads.py"), genome_fa(cfg, g), out + ".tmp",
                        "--length", str(L), "--reads", str(n), "--mode", mode, "--max-k", str(s["max_k"]),
                        "--random-fraction", str(s["random_fraction"]), "--seed", str(s["seed"] + i), *extra], check=True)
        os.rename(out + ".tmp", out)


def cmd_prepare(cfg, a):
    os.makedirs(wpath(cfg, "reads"), exist_ok=True)
    want = set(a.real.split(",")) if a.real else None
    for name, spec in cfg["real"].items():
        if want and name not in want:
            continue
        prefix = wpath(cfg, "reads", f"real_{name}")
        if os.path.exists(prefix + ".json"):
            continue
        print(f"prepare {name}: {spec['url']}", file=sys.stderr)
        subprocess.run([sys.executable, os.path.join(HERE, "prepare_reads.py"), spec["url"], prefix,
                        "--lengths", ",".join(map(str, spec["lengths"])), "--fraction", str(spec.get("sample_fraction", 1.0)),
                        "--max-reads", str(spec.get("max_reads", 0))], check=True)
        # Single-thread runs use the first single_thread_reads reads of each set.
        n1 = cfg.get("single_thread_reads", 0)
        for L in spec["lengths"]:
            fq = f"{prefix}.{L}.fq"
            os.rename(fq, f"{prefix}_{L}.fq")
            if n1:
                subprocess.run(f"head -n {4 * n1} {prefix}_{L}.fq > {prefix}_{L}_t1.fq", shell=True, check=True)


# ---------------------------------------------------------------- running one measured command


TIME_FIELDS = {
    "Elapsed (wall clock) time (h:mm:ss or m:ss)": "wall_s",
    "User time (seconds)": "user_s",
    "System time (seconds)": "sys_s",
    "Maximum resident set size (kbytes)": "max_rss_kb",
    "Exit status": "exit",
}


def parse_time_v(path):
    out = {}
    for line in open(path, errors="replace"):
        line = line.strip()
        if line.startswith("Command terminated by signal"):
            out["signal"] = int(line.split()[-1])
        for k, v in TIME_FIELDS.items():
            if line.startswith(k + ":"):
                val = line.split(": ", 1)[1].strip()
                if v == "wall_s":
                    parts = [float(x) for x in val.split(":")]
                    sec = 0.0
                    for p in parts:
                        sec = sec * 60 + p
                    out[v] = sec
                else:
                    out[v] = float(val) if "." in val else int(val)
    return out


def parse_batmis_timing(path):
    for line in open(path, errors="replace"):
        if " timing: " in line:
            fields = dict(re.findall(r"(\w+)=([^\s]+)", line))
            return {k: float(v.rstrip("s")) for k, v in fields.items()
                    if k in ("index_load", "table_build", "input_parse", "input_read", "search", "output_format", "total")}
    return {}


def run_phases(phases, run_dir, time_limit):
    """Runs the phases in order, each under /usr/bin/time -v; stops at the first failure."""
    results, deadline = [], time.time() + time_limit
    for name, argv, stdout_path in phases:
        tv = os.path.join(run_dir, f"{name}.time")
        err = os.path.join(run_dir, f"{name}.stderr")
        out = open(stdout_path, "wb") if stdout_path else subprocess.DEVNULL
        start = time.time()
        # In the run's directory: original BatMis writes its logs to the working directory.
        p = subprocess.Popen(["/usr/bin/time", "-v", "-o", tv, *argv], stdout=out, stderr=open(err, "wb"),
                             start_new_session=True, cwd=run_dir)
        try:
            p.wait(timeout=max(1, deadline - time.time()))
            status = "ok" if p.returncode == 0 else "error"  # a crash shows as "signal" from time
        except subprocess.TimeoutExpired:
            os.killpg(p.pid, signal.SIGKILL)
            p.wait()
            status = "timeout"
        if stdout_path:
            out.close()
        r = {"phase": name, "status": status, "argv": argv, "outer_wall_s": time.time() - start}
        if os.path.exists(tv):
            r.update(parse_time_v(tv))
        r.update({f"tool_{k}": v for k, v in parse_batmis_timing(err).items()})
        results.append(r)
        if status != "ok":
            break
    return results


# ---------------------------------------------------------------- the matrix


def datasets(cfg):
    """Every read set the runs use: dict name -> {genome, length, kind, fq, fq_t1}."""
    d = {}
    for name, g, L, mode, n, extra in sim_sets(cfg):
        kind = name.rsplit("_", 1)[1]  # timing, check, real, indel
        d[name] = {"genome": g, "length": L, "kind": kind, "fq": wpath(cfg, "reads", name + ".fq"), "reads": n}
    for name, spec in cfg["real"].items():
        for L in spec["lengths"]:
            fq = wpath(cfg, "reads", f"real_{name}_{L}.fq")
            d[f"real_{name}_{L}"] = {"genome": spec["genome"], "length": L, "kind": "realdata", "fq": fq,
                                     "fq_t1": wpath(cfg, "reads", f"real_{name}_{L}_t1.fq") if cfg.get("single_thread_reads") else fq}
    return d


def matrix(cfg):
    """Every run, as a dict; 'id' is unique and names its directory."""
    runs = []
    ks = cfg["k_exact"]
    kx = cfg.get("k_extra", [])
    reps = cfg["repeats"]
    main_t = cfg["threads"]
    tmax = max(main_t)
    for dname, ds in datasets(cfg).items():
        kind = ds["kind"]
        if kind in ("timing", "realdata"):
            for tool in T.EXACT + T.NOREF + T.SUPPLEMENTARY:
                for k in ks + kx:
                    if not T.accepts(tool, k) or (k in kx and tool not in T.HIGH_K):
                        continue
                    if tool in T.SUPPLEMENTARY and (k not in (0, 2, 4) or kind != "timing"):
                        continue
                    for mode in ("first", "all"):
                        for t in main_t:
                            if t > 1 and not T.THREADED.get(tool, True):
                                continue
                            for rep in range(reps + 1):
                                runs.append(dict(kind=kind, dataset=dname, genome=ds["genome"], length=ds["length"],
                                                 tool=tool, k=k, mode=mode, threads=t, rep=rep))
            for tool in T.APPROXIMATE:
                for mode in ("first", "all"):
                    for t in main_t:
                        for rep in range(reps + 1):
                            runs.append(dict(kind=kind, dataset=dname, genome=ds["genome"], length=ds["length"],
                                             tool=tool, k=-1, mode=mode, threads=t, rep=rep))
        elif kind == "check":
            for tool in T.EXACT + T.NOREF:
                for k in ks + kx:
                    if not T.accepts(tool, k) or (k in kx and tool not in T.HIGH_K):
                        continue
                    for mode in ("first", "all"):
                        runs.append(dict(kind=kind, dataset=dname, genome=ds["genome"], length=ds["length"],
                                         tool=tool, k=k, mode=mode, threads=tmax if T.THREADED.get(tool, True) else 1, rep=1))
            for tool in T.APPROXIMATE:
                runs.append(dict(kind=kind, dataset=dname, genome=ds["genome"], length=ds["length"],
                                 tool=tool, k=-1, mode="first", threads=tmax, rep=1))
        elif kind in ("real", "indel"):  # true-origin recovery on realistic reads
            for tool in T.EXACT + T.APPROXIMATE:
                for k in ([-1] if tool in T.APPROXIMATE else [k for k in ks if T.accepts(tool, k)]):
                    runs.append(dict(kind=kind, dataset=dname, genome=ds["genome"], length=ds["length"],
                                     tool=tool, k=k, mode="first", threads=tmax if T.THREADED.get(tool, True) else 1, rep=1))
    # Thread scaling: human 100 bp, timing set, K = 2 and 4, first and all.
    sc = cfg.get("scaling", {})
    for dname in sc.get("datasets", []):
        ds = datasets(cfg)[dname]
        for tool in sc.get("tools", []):
            for k in sc.get("k", []):
                for mode in ("first", "all"):
                    for t in cfg["scaling_threads"]:
                        if t in main_t:
                            continue
                        runs.append(dict(kind="scaling", dataset=dname, genome=ds["genome"], length=ds["length"],
                                         tool=tool, k=k, mode=mode, threads=t, rep=1))
    # Calibration (kind calib): 5,000 reads, one thread, at the calibration's two K, for every
    # timing and real-data set of a calibrated genome (bench.py calibrate runs these).
    cal = calibration(cfg)
    for dname, ds in datasets(cfg).items():
        if ds["kind"] not in ("timing", "realdata") or not calibrated(cfg, ds["genome"]):
            continue
        for tool in T.EXACT + T.NOREF + T.SUPPLEMENTARY + T.APPROXIMATE:
            kk = [-1] if tool in T.APPROXIMATE else sorted({min(k, T.MAX_K.get(tool, 7)) for k in cal.get("k", [0, 5])})
            for k in kk:
                for mode in ("first", "all"):
                    runs.append(dict(kind="calib", dataset=dname, genome=ds["genome"], length=ds["length"],
                                     tool=tool, k=k, mode=mode, threads=1, rep=1))
    # Index load: one read, per genome and tool.
    for g in {ds["genome"] for ds in datasets(cfg).values()}:
        for tool in T.EXACT + T.NOREF + T.SUPPLEMENTARY + T.APPROXIMATE:
            runs.append(dict(kind="load", dataset=f"one_{g}", genome=g, length=0, tool=tool,
                             k=0 if tool in T.APPROXIMATE else min(2, T.MAX_K.get(tool, 7)), mode="first", threads=1, rep=1))
    for r in runs:
        r["id"] = f"{r['kind']}/{r['dataset']}/{r['tool']}/k{r['k']}/{r['mode']}/t{r['threads']}/r{r['rep']}"
    return runs


def select(runs, a):
    out = runs
    for cond in (a.set or "").split(","):
        if not cond:
            continue
        key, val = cond.split("=", 1)
        vals = set(val.split("|"))
        out = [r for r in out if str(r.get(key)) in vals]
    return out


def run_dir(cfg, r):
    return wpath(cfg, "runs", r["id"])


def done(cfg, r):
    return os.path.exists(os.path.join(run_dir(cfg, r), "run.json"))


def time_limit(cfg, r):
    """The run's time limit: time_limit_s, or the genome's entry in time_limit_s_genome."""
    return (cfg.get("time_limit_s_genome") or {}).get(r["genome"], cfg["time_limit_s"])


def blocked(cfg, r):
    """A reason not to run r: a smaller K of the same series timed out, its warm-up timed out
    (a timed repeat would only time out again), or its warm-up ran longer than
    single_repeat_above_s (then one timed repeat is enough)."""
    if r["rep"] >= 1:
        warm = os.path.join(wpath(cfg, "runs", r["id"].rsplit("/r", 1)[0] + "/r0"), "run.json")
        if os.path.exists(warm):
            w = json.load(open(warm))
            if w.get("status") == "timeout":
                return "warm-up timed out"
            if r["rep"] >= 2 and w.get("wall_s", 0) > cfg.get("single_repeat_above_s", 1800):
                return "long run: one timed repeat"
    if r["k"] > 0:
        for k in range(r["k"]):
            prev = os.path.join(wpath(cfg, "runs", re.sub(r"/k\d+/", f"/k{k}/", r["id"])), "run.json")
            if os.path.exists(prev) and json.load(open(prev)).get("status") == "timeout":
                return f"k={k} timed out"
    return None


def reads_for(cfg, r):
    if r["kind"] == "load":
        p = wpath(cfg, "reads", f"one_{r['genome']}.fq")
        if not os.path.exists(p):
            # One read of the genome's first ACGT window, 100 bp.
            seq = b""
            with open(genome_fa(cfg, r["genome"]), "rb") as f:
                for line in f:
                    if line.startswith(b">"):
                        continue
                    seq += line.strip().upper()
                    m = re.search(rb"[ACGT]{100}", seq)
                    if m:
                        seq = m.group(0)
                        break
                    seq = seq[-200:]
            open(p, "wb").write(b"@one\n" + seq + b"\n+\n" + b"I" * 100 + b"\n")
        return p, 100
    ds = datasets(cfg)[r["dataset"]]
    if r["kind"] == "calib":
        return subset_reads(ds["fq"], calibration(cfg).get("reads", 5000)), ds["length"]
    if r["kind"] in ("check", "real", "indel"):
        return correctness_fq(cfg, ds), ds["length"]
    if calibrated(cfg, ds["genome"]) and r["kind"] in ("timing", "realdata", "scaling"):
        n = planned_reads(cfg, r)
        return (subset_reads(ds["fq"], n), ds["length"]) if n else (None, ds["length"])
    fq = ds.get("fq_t1", ds["fq"]) if r["threads"] == 1 else ds["fq"]
    return fq, ds["length"]


# ---------------------------------------------------------------- calibration


def calib_plan_path(cfg):
    return wpath(cfg, "runs", "calib", "plan.json")


_PLAN = {}


def load_plan(cfg):
    p = calib_plan_path(cfg)
    if p not in _PLAN:
        _PLAN[p] = json.load(open(p)) if os.path.exists(p) else None
    return _PLAN[p]


def per_read_seconds(cfg, r):
    """The calibrated one-thread time per read of r's tool, set and mode at r's K: log-linear in
    K between (and beyond) the two calibration points. None if not calibrated."""
    plan = load_plan(cfg)
    if not plan:
        return None
    e = plan["fits"].get(f"{r['dataset']}|{r['tool']}|{r['mode']}")
    if not e:
        return None
    pts = sorted((int(k), t) for k, t in e["per_read_s"].items())
    if len(pts) == 1 or r["k"] < 0:
        return pts[0][1]
    (k0, t0), (k1, t1) = pts[0], pts[-1]
    return t0 * (t1 / t0) ** ((r["k"] - k0) / (k1 - k0))


def planned_reads(cfg, r):
    """Reads for a calibrated run: as many as fit the target time on its threads (linear in
    reads and threads), on the ladder, at least `reads`, at most the whole set. None before
    calibration."""
    t = per_read_seconds(cfg, r)
    if t is None:
        return None
    cal = calibration(cfg)
    speed = r["threads"] if T.THREADED.get(r["tool"], True) else 1
    want = cal.get("target_s", 600) * speed / t
    fq = datasets(cfg)[r["dataset"]]["fq"]
    if not os.path.exists(fq):
        return None
    total = count_reads(fq)
    n = cal.get("reads", 5000)
    for step in cal.get("ladder", []):
        if step <= want:
            n = max(n, step)
    return total if want >= total else min(n, total)


def startup_seconds(cfg, genome, tool):
    """A tool's start-up on a genome: the wall time of its one-read load run (0 if missing)."""
    for p in glob.glob(wpath(cfg, "runs", "load", f"one_{genome}", tool, "*", "*", "*", "*", "run.json")):
        rec = json.load(open(p))
        if rec.get("status") == "ok":
            return rec["wall_s"]
    return 0.0


def estimate_matrix(cfg, gens):
    """Prints the estimated node time of the calibrated genomes' runs still to do (one node,
    runs one after another), by kind and mode, so that target_s or repeats can be adjusted
    before the runs start. Runs over the time limit count with the limit."""
    plan = load_plan(cfg)
    if not plan:
        print(f"no calibration yet ({calib_plan_path(cfg)} missing): run `bench.py calibrate` first",
              file=sys.stderr)
        return
    # Calibration runs that did not give a fit, so that their runs cannot be sized.
    failed = []
    for r in matrix(cfg):
        if r["kind"] == "calib" and r["genome"] in gens:
            rp = os.path.join(run_dir(cfg, r), "run.json")
            st = json.load(open(rp)).get("status") if os.path.exists(rp) else "not run"
            if st not in ("ok", "timeout"):
                failed.append((r["id"], st))
    print(f"calibration: {len(plan.get('fits', {}))} fits ({calib_plan_path(cfg)}); "
          f"{len(failed)} calibration runs failed or missing", file=sys.stderr)
    for rid, st in failed[:10]:
        print(f"  {st:8s} {rid}  (see its *.stderr / *.time)", file=sys.stderr)
    if len(failed) > 10:
        print(f"  ... {len(failed) - 10} more", file=sys.stderr)
    total, parts, n = 0.0, {}, 0
    for r in matrix(cfg):
        if r["genome"] not in gens or r["kind"] not in ("timing", "realdata", "scaling") or done(cfg, r):
            continue
        t, reads = per_read_seconds(cfg, r), planned_reads(cfg, r)
        if t is None or reads is None:
            continue
        speed = r["threads"] if T.THREADED.get(r["tool"], True) else 1
        est = min(time_limit(cfg, r), startup_seconds(cfg, r["genome"], r["tool"]) + t * reads / speed)
        total += est
        n += 1
        key = f"{r['kind']} {r['mode']} t{r['threads']}"
        parts[key] = parts.get(key, 0.0) + est
    print(f"estimated time of the {n} sized runs still to do: {total / 3600:.1f} node-hours "
          f"({total / 86400:.1f} node-days); target_s {calibration(cfg).get('target_s')}, repeats {cfg['repeats']} + warm-up",
          file=sys.stderr)
    for key, v in sorted(parts.items()):
        print(f"  {key:24s} {v / 3600:8.1f} h", file=sys.stderr)


def cmd_estimate(cfg, a):
    estimate_matrix(cfg, calibration(cfg).get("genomes", []))


def cmd_calibrate(cfg, a):
    """Runs the load and calibration runs of the calibrated genomes, then fits each (set, tool,
    mode): per-read seconds at each calibration K, (wall - start-up) / reads. A calibration run
    over the time limit counts with the limit (a lower bound)."""
    gens = calibration(cfg).get("genomes", [])
    if not gens:
        sys.exit("config has no calibration genomes")
    runs = [r for r in matrix(cfg) if r["kind"] in ("load", "calib") and r["genome"] in gens]
    for r in runs:
        rp = os.path.join(run_dir(cfg, r), "run.json")
        # Failed runs are run again (e.g. after an index is rebuilt); timeouts are kept.
        if not done(cfg, r) or json.load(open(rp)).get("status") == "error":
            execute(cfg, r)
    fits = {}
    n = calibration(cfg).get("reads", 5000)
    for r in runs:
        if r["kind"] != "calib":
            continue
        rec = json.load(open(os.path.join(run_dir(cfg, r), "run.json")))
        if rec.get("status") not in ("ok", "timeout"):
            print(f"calibration failed ({rec.get('status')}): {r['id']}", file=sys.stderr)
            continue
        reads = rec.get("reads_n") or n
        t = max(1e-6, (rec["wall_s"] - startup_seconds(cfg, r["genome"], r["tool"])) / reads)
        key = f"{r['dataset']}|{r['tool']}|{r['mode']}"
        fits.setdefault(key, {"per_read_s": {}, "status": {}})
        fits[key]["per_read_s"][str(r["k"])] = t
        fits[key]["status"][str(r["k"])] = rec["status"]
    os.makedirs(os.path.dirname(calib_plan_path(cfg)), exist_ok=True)
    json.dump({"target_s": calibration(cfg).get("target_s", 600), "reads": n, "fits": fits},
              open(calib_plan_path(cfg), "w"), indent=1)
    _PLAN.clear()
    print(f"calibration: {len(fits)} fits -> {calib_plan_path(cfg)}", file=sys.stderr)
    estimate_matrix(cfg, gens)
    for r in matrix(cfg):
        if r["kind"] in ("timing", "realdata") and r["genome"] in gens and r["rep"] == 1 and r["threads"] == 1:
            t = per_read_seconds(cfg, r)
            if t is not None:
                print(f"  {r['dataset']:22s} {r['tool']:14s} K={r['k']:2d} {r['mode']:5s} {t * 1e6:10.1f} us/read "
                      f"-> {planned_reads(cfg, r):>9d} reads on 1 thread", file=sys.stderr)


def output_dir(cfg, r):
    """Where a run writes its SAM (and BWA's .sai, BatMis 3.00's .bin): local_output if set,
    so that output goes to node-local disk, not the shared filesystem; else the run's directory."""
    local = cfg.get("local_output")
    if not local:
        return run_dir(cfg, r)
    return os.path.join(local, r["id"].replace("/", "_"))


def execute(cfg, r, dry_run=False):
    d = run_dir(cfg, r)
    reads, L = reads_for(cfg, r)
    if reads is None:
        print(f"not calibrated yet (bench.py calibrate): {r['id']}", flush=True)
        return
    od = output_dir(cfg, r)
    out = os.path.join(od, "out.sam")
    extra = []
    opts = cfg.get("tool_options", {}).get(r["tool"].replace("_noref", ""), {})  # no-ref: batmis's options
    if "large" in opts and r["kind"] != "load" and os.path.exists(reads) and \
            os.path.getsize(reads) > cfg.get("large_reads", 500000) * (2 * L + 30):
        extra = opts["large"]
    phases = T.map_reads(r["tool"], cfg["tools"], index_dir(cfg, r["genome"], T.INDEX_OF[r["tool"]]), r["genome"],
                         reads, L, max(r["k"], 0), r["mode"], r["threads"], out, extra)
    if dry_run:
        for name, argv, so in phases:
            print(f"  {name}: {' '.join(argv)}" + (f" > {so}" if so else ""))
        return
    shutil.rmtree(d, ignore_errors=True)
    os.makedirs(d)
    if od != d:  # a stale local directory from an interrupted run
        shutil.rmtree(od, ignore_errors=True)
        os.makedirs(od)
    res = run_phases(phases, d, time_limit(cfg, r))
    status = "ok" if all(p["status"] == "ok" for p in res) and len(res) == len(phases) else \
        ("timeout" if any(p["status"] == "timeout" for p in res) else "error")
    rec = dict(r, status=status, reads_file=reads, phases=res,
               wall_s=sum(p.get("wall_s", p["outer_wall_s"]) for p in res),
               user_s=sum(p.get("user_s", 0) for p in res), sys_s=sum(p.get("sys_s", 0) for p in res),
               max_rss_kb=max((p.get("max_rss_kb", 0) for p in res), default=0),
               output_bytes=os.path.getsize(out) if os.path.exists(out) else 0,
               reads_n=count_reads(reads), startup_s=startup_seconds(cfg, r["genome"], r["tool"]),
               finished=time.strftime("%Y-%m-%dT%H:%M:%S"))
    t = per_read_seconds(cfg, r) if calibrated(cfg, r["genome"]) and r["kind"] in ("timing", "realdata", "scaling") else None
    if t is not None:  # the calibration's estimate, and how far the run is from it
        speed = r["threads"] if T.THREADED.get(r["tool"], True) else 1
        rec["estimated_s"] = rec["startup_s"] + t * rec["reads_n"] / speed
        # Flag a real miss only: short runs are dominated by fixed costs the one-read load run
        # does not see (threads starting, a lookup table built for large read sets).
        if rec["wall_s"] > 3 * rec["estimated_s"] and rec["wall_s"] - rec["estimated_s"] > 10:
            print(f"  note: {rec['wall_s']:.1f} s against an estimate of {rec['estimated_s']:.1f} s "
                  f"(not linear here): {r['id']}", flush=True)
    for p in res:  # tool-reported phase times (batmis --timing)
        for k, v in p.items():
            if k.startswith("tool_"):
                rec[k] = v
    if len(res) > 1:  # align / decode split
        for p in res:
            rec[f"{p['phase']}_wall_s"] = p.get("wall_s", p["outer_wall_s"])
            rec[f"{p['phase']}_max_rss_kb"] = p.get("max_rss_kb", 0)
    # Keep the hits of the correctness runs (normalised) and drop the SAM of everything else.
    if r["kind"] in ("check", "real", "indel") and status == "ok":
        subprocess.run([sys.executable, os.path.join(HERE, "normalize_hits.py"), out, os.path.join(d, "hits.tsv")], check=True)
    if cfg.get("keep_sam") and od != d:
        for f in glob.glob(os.path.join(od, "out.sam*")):
            shutil.move(f, d)
    if not cfg.get("keep_sam"):
        for f in glob.glob(os.path.join(od, "out.sam*")):
            os.remove(f)
    if od != d:
        shutil.rmtree(od, ignore_errors=True)
    json.dump(rec, open(os.path.join(d, "run.json.tmp"), "w"), indent=1)
    os.rename(os.path.join(d, "run.json.tmp"), os.path.join(d, "run.json"))
    print(f"{status:7s} {rec['wall_s']:9.2f}s {rec['max_rss_kb'] / 1024:8.0f}MB  {r['id']}", flush=True)


def cmd_list(cfg, a):
    runs = select(matrix(cfg), a)
    n_done = sum(done(cfg, r) for r in runs)
    for r in runs:
        n = ""
        if calibrated(cfg, r["genome"]) and r["kind"] in ("timing", "realdata", "scaling"):
            p = planned_reads(cfg, r)
            n = f"  [{p} reads]" if p else "  [not calibrated]"
        print(("done " if done(cfg, r) else "     ") + r["id"] + n)
    print(f"{len(runs)} runs, {n_done} done", file=sys.stderr)


def ago(s):
    s = int(s)
    return f"{s // 86400}d {s % 86400 // 3600}h" if s >= 86400 else f"{s // 3600}h {s % 3600 // 60}m"


def cmd_status(cfg, a):
    """Progress of the selected runs, read-only: done, timeouts, errors and skips per tool; the
    newest finished run; and the runs in progress (phase files but no run.json) with the age
    of their current phase. A run in progress longer than its time limit has a process the
    timeout could not kill, or was left by a job that is no longer running (it is redone)."""
    now, tools, newest, active = time.time(), {}, None, []
    for r in select(matrix(cfg), a):
        t = tools.setdefault(r["tool"], dict(runs=0, done=0, ok=0, timeout=0, error=0, skipped=0))
        t["runs"] += 1
        d = run_dir(cfg, r)
        rp = os.path.join(d, "run.json")
        if os.path.exists(rp):
            t["done"] += 1
            st = json.load(open(rp)).get("status", "error")
            t[st if st in t else "error"] += 1
            m = os.path.getmtime(rp)
            if newest is None or m > newest[0]:
                newest = (m, st, r["id"])
        else:
            if glob.glob(os.path.join(d, "*.stderr")):
                started = os.path.getmtime(d)  # the last file created there: the current phase's start
                last = max(os.path.getmtime(p) for p in glob.glob(os.path.join(d, "*")))
                active.append((now - started, now - last, time_limit(cfg, r), r["id"]))
    tot = {k: sum(t[k] for t in tools.values()) for k in ("runs", "done", "ok", "timeout", "error", "skipped")}
    print(f"{'tool':16s} {'runs':>6s} {'done':>6s} {'ok':>6s} {'timeout':>7s} {'error':>6s} {'skipped':>7s}")
    for name, t in sorted(tools.items()) + [("all", tot)]:
        print(f"{name:16s} {t['runs']:6d} {t['done']:6d} {t['ok']:6d} {t['timeout']:7d} {t['error']:6d} {t['skipped']:7d}")
    if newest:
        print(f"newest finished: {time.strftime('%Y-%m-%d %H:%M', time.localtime(newest[0]))} "
              f"({ago(now - newest[0])} ago) {newest[1]} {newest[2]}")
    for age, quiet, limit, rid in sorted(active):
        print(f"in progress: {rid}  phase started {ago(age)} ago, last write {ago(quiet)} ago"
              + (f"  OVER THE {ago(limit)} LIMIT, a hung process, or left by a stopped job" if age > limit + 600 else ""))
    if not active:
        print("in progress: none (no job is running these runs)")


def cmd_run(cfg, a):
    # Calibration runs belong to `calibrate`, which fits the plan after them.
    runs = [r for r in select(matrix(cfg), a) if not done(cfg, r) and r["kind"] != "calib"]
    waiting = [r for r in runs if calibrated(cfg, r["genome"]) and r["kind"] in ("timing", "realdata", "scaling")
               and planned_reads(cfg, r) is None]
    if waiting:
        print(f"{len(waiting)} runs wait for `bench.py calibrate` (genomes {', '.join(calibration(cfg)['genomes'])}); "
              "the others run now", file=sys.stderr)
        runs = [r for r in runs if r not in waiting]
    # Small genomes first, so that the pilot's results come early.
    order = {g: i for i, g in enumerate(["sarscov2", "ecoli", "chr21", "human"])}
    runs.sort(key=lambda r: (order.get(r["genome"], 9), r["dataset"], r["tool"], r["mode"], r["threads"], r["k"], r["rep"]))
    for r in runs:
        why = blocked(cfg, r)
        if why:
            d = run_dir(cfg, r)
            os.makedirs(d, exist_ok=True)
            json.dump(dict(r, status="skipped", reason=why), open(os.path.join(d, "run.json"), "w"), indent=1)
            print(f"skipped {why}: {r['id']}", flush=True)
            continue
        if a.dry_run:
            print(r["id"])
        execute(cfg, r, a.dry_run)


# ---------------------------------------------------------------- indexes


def sa_interval(sa_file):
    """The sampling interval stored in a BatMis .sa file (after inverseSa0 and 4 counts)."""
    import struct
    with open(sa_file, "rb") as f:
        return struct.unpack("<6I", f.read(24))[5]


def index_complete(cfg, g, i):
    """A finished, successful build of index i of genome g: its record says ok and, for BatMis,
    every file the mappers load is there (a job killed mid-build leaves some out)."""
    rp = wpath(cfg, "runs", "index", g, i, "run.json")
    if not os.path.exists(rp) or json.load(open(rp)).get("status") != "ok":
        return False
    if i == "batmis":
        d = index_dir(cfg, g, i)
        need = [f"{p}{g}.fa.{x}" for p in ("", "rev") for x in ("bwt", "fmv", "sa")] + [f"{g}.fa.pac"]
        return all(os.path.exists(os.path.join(d, f)) for f in need)
    return True


def cmd_index(cfg, a):
    gs = a.genomes.split(",") if a.genomes else sorted({ds["genome"] for ds in datasets(cfg).values()})
    idxs = a.indexes.split(",") if a.indexes else T.INDEXES
    threads = max(cfg["threads"])
    for g in gs:
        for i in idxs:
            d = index_dir(cfg, g, i)
            rd = wpath(cfg, "runs", "index", g, i)
            if index_complete(cfg, g, i):
                continue
            if os.path.exists(os.path.join(rd, "run.json")):
                print(f"index {g}/{i}: previous build failed or is incomplete; rebuilding", flush=True)
            shutil.rmtree(d, ignore_errors=True)
            shutil.rmtree(rd, ignore_errors=True)
            os.makedirs(d)
            os.makedirs(rd)
            fa = os.path.join(d, f"{g}.fa") if i == "batmis" else genome_fa(cfg, g)
            if i == "batmis":  # bwt_build_index names its files after the FASTA: <g>.fa.*
                os.symlink(genome_fa(cfg, g), fa)
            res = run_phases(T.build_index(i, cfg["tools"], fa, d, g, threads), rd, cfg.get("index_time_limit_s", 86400))
            # Every file the tool keeps in its index directory (mrsFAST's includes the FASTA
            # copy it searches with); not batmis's link to the genome.
            size = sum(os.path.getsize(f) for f in glob.glob(os.path.join(d, "*")) if os.path.isfile(f) and not os.path.islink(f))
            status = "ok" if all(p["status"] == "ok" for p in res) else "error"
            rec = dict(kind="index", genome=g, index=i, status=status, phases=res,
                       wall_s=sum(p.get("wall_s", 0) for p in res), max_rss_kb=max(p.get("max_rss_kb", 0) for p in res),
                       index_bytes=size, threads=threads)
            if i == "batmis" and status == "ok":  # the SA sampling the builder's budget chose
                rec["sa_interval"] = sa_interval(os.path.join(d, f"{g}.fa.sa"))
            json.dump(rec, open(os.path.join(rd, "run.json"), "w"), indent=1)
            print(f"index {g}/{i}: {status} {rec['wall_s']:.1f}s {rec['max_rss_kb'] / 1024:.0f}MB {size / 2**20:.0f}MB", flush=True)


# ---------------------------------------------------------------- correctness


def cmd_truth(cfg, a):
    """Brute-force hits (up to the largest K) of the correctness sets of brute-forceable genomes."""
    bf = os.path.join(cfg["tools"], "brute_force")
    if not os.path.exists(bf):
        subprocess.run(["g++", "-O3", "-march=native", "-std=c++17", "-pthread", "-o", bf,
                        os.path.join(HERE, "brute_force.cpp")], check=True)
    kmax = max(cfg["k_exact"] + cfg.get("k_extra", []))
    want = set(a.datasets.split(",")) if a.datasets else None
    for name, ds in datasets(cfg).items():
        if ds["kind"] != "check" or ds["genome"] not in cfg.get("brute_force_genomes", []) or (want and name not in want):
            continue
        out = wpath(cfg, "truth", name + ".tsv")
        if os.path.exists(out):
            continue
        os.makedirs(os.path.dirname(out), exist_ok=True)
        print(f"truth {name} (K <= {kmax})", file=sys.stderr, flush=True)
        t0 = time.time()
        with open(out + ".tmp", "wb") as o:
            p1 = subprocess.Popen([bf, genome_fa(cfg, ds["genome"]), correctness_fq(cfg, ds), str(kmax)],
                                  stdout=subprocess.PIPE)
            subprocess.run(["sort", "-u", "-S", "25%"], stdin=p1.stdout, stdout=o, check=True, env={"LC_ALL": "C", "PATH": "/usr/bin:/bin"})
            p1.wait()
        os.rename(out + ".tmp", out)
        print(f"  {time.time() - t0:.0f}s", file=sys.stderr)


def cmd_correct(cfg, a):
    runs = [r for r in select(matrix(cfg), a) if r["kind"] in ("check", "real", "indel") and done(cfg, r)]
    by_set = {}
    for r in runs:
        by_set.setdefault(r["dataset"], []).append(r)
    for name, rs in sorted(by_set.items()):
        ds = datasets(cfg)[name]
        truth = wpath(cfg, "truth", name + ".tsv")
        args = [sys.executable, os.path.join(HERE, "compare_hits.py"), "--genome", genome_fa(cfg, ds["genome"]),
                "--reads", correctness_fq(cfg, ds),
                "--out", wpath(cfg, "runs", "correct", name + ".json")]
        if os.path.exists(truth):
            args += ["--truth", truth]
        for r in rs:
            hits = os.path.join(run_dir(cfg, r), "hits.tsv")
            if os.path.exists(hits):
                args += ["--run", f"{r['tool']}:{r['k']}:{r['mode']}:{hits}"]
        os.makedirs(wpath(cfg, "runs", "correct"), exist_ok=True)
        print(f"correct {name}: {len(rs)} runs" + (" against brute force" if os.path.exists(truth) else " against the verified union"),
              file=sys.stderr, flush=True)
        subprocess.run(args, check=True)


def cmd_discrepancies(cfg, a):
    """discrepancies.py per correctness set (check kind: the sets with all-hits runs), into
    <results>/discrepancies/<set>/; with --read, one read of one set (--set dataset=...)."""
    runs = [r for r in select(matrix(cfg), a) if r["kind"] == "check" and done(cfg, r)]
    by_set = {}
    for r in runs:
        by_set.setdefault(r["dataset"], []).append(r)
    if a.read and len(by_set) != 1:
        sys.exit("--read needs --set dataset=<one read set>")
    for name, rs in sorted(by_set.items()):
        ds = datasets(cfg)[name]
        truth = wpath(cfg, "truth", name + ".tsv")
        args = [sys.executable, os.path.join(HERE, "discrepancies.py"), "--genome", genome_fa(cfg, ds["genome"]),
                "--reads", correctness_fq(cfg, ds)]
        args += ["--read", a.read] if a.read else ["--out", os.path.join(cfg["results"], "discrepancies", name)]
        if os.path.exists(truth):
            args += ["--truth", truth]
        for r in rs:
            hits = os.path.join(run_dir(cfg, r), "hits.tsv")
            if os.path.exists(hits):
                args += ["--run", f"{r['tool']}:{r['k']}:{r['mode']}:{hits}"]
        if not a.read:
            print(f"discrepancies {name}", file=sys.stderr, flush=True)
        subprocess.run(args, check=True)


def cmd_collect(cfg, a):
    os.makedirs(cfg["results"], exist_ok=True)
    # The calibration, for the report (which genomes' runs compare per read).
    plan = load_plan(cfg) or {}
    json.dump(dict(plan, genomes=calibration(cfg).get("genomes", [])),
              open(os.path.join(cfg["results"], "calibration.json"), "w"), indent=1)
    # The runs of the current matrix: collect leaves out records of runs no longer in it (e.g.
    # human 40/150 bp runs made before human was restricted to 100 bp).
    ids = os.path.join(cfg["results"], "matrix_ids.txt")
    with open(ids, "w") as f:
        f.writelines(r["id"] + "\n" for r in matrix(cfg))
    subprocess.run([sys.executable, os.path.join(HERE, "collect.py"), wpath(cfg, "runs"), cfg["results"],
                    os.path.join(cfg["tools"], "manifest.json"), ids], check=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-c", "--config", default=os.path.join(HERE, "config.yaml"))
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("fetch")
    p = sub.add_parser("simulate"); p.add_argument("--genomes")
    p = sub.add_parser("prepare"); p.add_argument("--real")
    p = sub.add_parser("index"); p.add_argument("--genomes"); p.add_argument("--indexes")
    sub.add_parser("calibrate")
    sub.add_parser("estimate")
    for name in ("list", "status", "run", "correct", "discrepancies"):
        p = sub.add_parser(name); p.add_argument("--set")
        if name == "run":
            p.add_argument("--dry-run", action="store_true")
        if name == "discrepancies":
            p.add_argument("--read", help="show one read of the one selected set")
    p = sub.add_parser("truth"); p.add_argument("--datasets")
    sub.add_parser("collect")
    a = ap.parse_args()
    cfg = load_config(a.config)
    globals()["cmd_" + a.cmd](cfg, a)


if __name__ == "__main__":
    main()
