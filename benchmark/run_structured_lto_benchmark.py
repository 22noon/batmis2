#!/usr/bin/env python3
"""Alternating minimized, structured, and LTO structured exact-search trials."""
import argparse
import hashlib
import json
import os
import re
import resource
import subprocess
import time


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def timings(stderr):
    line = next((s for s in stderr.splitlines() if " timing:" in s), None)
    if line is None:
        raise RuntimeError(f"No application timing line in stderr: {stderr}")
    result = {}
    for key, value, _unit in re.findall(r"([a-z_]+)=([0-9.]+)(s)?", line):
        result[key] = float(value)
    for key in ("records", "reported_hits"):
        match = re.search(rf"{key}=(\d+)", line)
        if not match:
            raise RuntimeError(f"Missing {key} in timing line: {line}")
        result[key] = int(match.group(1))
    return result


def run(binary, variant, block, order_position, reads, prefix, metadata):
    command = [binary, "-q", reads, "-g", prefix, "-m", "0", "--quiet", "--timing"]
    load_before = os.getloadavg()
    cpu_before = resource.getrusage(resource.RUSAGE_CHILDREN).ru_utime
    wall_start = time.monotonic()
    completed = subprocess.run(command, text=True, capture_output=True, check=True)
    wall = time.monotonic() - wall_start
    cpu = resource.getrusage(resource.RUSAGE_CHILDREN).ru_utime - cpu_before
    row = {
        **metadata,
        "variant": variant,
        "block": block,
        "order_position": order_position,
        "binary": os.path.realpath(binary),
        "binary_sha256": sha256(binary),
        "command": command,
        "external_wall_s": wall,
        "user_cpu_s": cpu,
        "load_before": load_before,
        "load_after": os.getloadavg(),
        **timings(completed.stderr),
    }
    return row


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--minimized", required=True)
    parser.add_argument("--structured", required=True)
    parser.add_argument("--lto", required=True)
    parser.add_argument("--reads", required=True)
    parser.add_argument("--genome-prefix", required=True)
    parser.add_argument("--trials", type=int, default=10)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    prefix = os.path.abspath(args.genome_prefix)
    directory, basename = os.path.split(prefix)
    reverse_prefix = os.path.join(directory, "rev" + basename)
    metadata = {
        "reads": os.path.abspath(args.reads),
        "reads_sha256": sha256(args.reads),
        "genome_prefix": prefix,
        "reverse_bwt_sha256": sha256(reverse_prefix + ".bwt"),
        "reverse_fmv_sha256": sha256(reverse_prefix + ".fmv"),
        "workload": "100000 100bp reads; m=0, quiet, no SA load/location/output",
    }
    binaries = {
        "minimized": args.minimized,
        "structured": args.structured,
        "lto": args.lto,
    }
    orders = [
        ("minimized", "structured", "lto"),
        ("lto", "minimized", "structured"),
        ("structured", "lto", "minimized"),
        ("structured", "minimized", "lto"),
        ("lto", "structured", "minimized"),
        ("minimized", "lto", "structured"),
    ]
    rows = []
    for block in range(1, args.trials + 1):
        for position, variant in enumerate(orders[(block - 1) % len(orders)], start=1):
            row = run(binaries[variant], variant, block, position, args.reads, prefix, metadata)
            rows.append(row)
            print(json.dumps(row, sort_keys=True), flush=True)

    counts = {(r["records"], r["reported_hits"]) for r in rows}
    if len(counts) != 1:
        raise SystemExit(f"Configuration count mismatch: {sorted(counts)}")
    with open(args.output, "w") as stream:
        for row in rows:
            stream.write(json.dumps(row, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
