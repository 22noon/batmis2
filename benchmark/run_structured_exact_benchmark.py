#!/usr/bin/env python3
"""Alternating baseline/structured BatMIS exact-search benchmark."""
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


def timing_fields(stderr):
    line = next((line for line in stderr.splitlines() if " timing:" in line), None)
    if line is None:
        raise RuntimeError(f"No application timing line in stderr: {stderr}")
    fields = {}
    for key, value, _unit in re.findall(r"([a-z_]+)=([0-9.]+)(s)?", line):
        fields[key] = float(value)
    for key in ("records", "reported_hits"):
        match = re.search(rf"{key}=(\d+)", line)
        if not match:
            raise RuntimeError(f"Missing {key} in timing line: {line}")
        fields[key] = int(match.group(1))
    return fields


def run(binary, variant, pair, position, reads, genome_prefix, metadata):
    command = [binary, "-q", reads, "-g", genome_prefix, "-m", "0", "--quiet", "--timing"]
    load_before = os.getloadavg()
    cpu_before = resource.getrusage(resource.RUSAGE_CHILDREN).ru_utime
    wall_start = time.monotonic()
    result = subprocess.run(command, text=True, capture_output=True, check=True)
    wall = time.monotonic() - wall_start
    cpu = resource.getrusage(resource.RUSAGE_CHILDREN).ru_utime - cpu_before
    load_after = os.getloadavg()
    row = {
        **metadata,
        "variant": variant,
        "pair": pair,
        "order_position": position,
        "binary": os.path.realpath(binary),
        "binary_sha256": sha256(binary),
        "command": command,
        "external_wall_s": wall,
        "user_cpu_s": cpu,
        "load_before": load_before,
        "load_after": load_after,
        **timing_fields(result.stderr),
    }
    return row


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--baseline", required=True)
    parser.add_argument("--structured", required=True)
    parser.add_argument("--reads", required=True)
    parser.add_argument("--genome-prefix", required=True)
    parser.add_argument("--trials", type=int, default=5)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    base = os.path.abspath(args.genome_prefix)
    directory, basename = os.path.split(base)
    reverse_prefix = os.path.join(directory, "rev" + basename)
    metadata = {
        "reads": os.path.abspath(args.reads),
        "reads_sha256": sha256(args.reads),
        "genome_prefix": base,
        "reverse_bwt_sha256": sha256(reverse_prefix + ".bwt"),
        "reverse_fmv_sha256": sha256(reverse_prefix + ".fmv"),
        "workload": "100000 100bp reads; quiet, no SA load/location/output; m=0",
    }

    rows = []
    for pair in range(1, args.trials + 1):
        order = (("baseline", args.baseline), ("structured", args.structured))
        if pair % 2 == 0:
            order = tuple(reversed(order))
        for position, (variant, binary) in enumerate(order, start=1):
            row = run(binary, variant, pair, position, args.reads, base, metadata)
            rows.append(row)
            print(json.dumps(row, sort_keys=True), flush=True)

    counts = {(row["records"], row["reported_hits"]) for row in rows}
    if len(counts) != 1:
        raise SystemExit(f"Configuration count mismatch: {sorted(counts)}")
    with open(args.output, "w") as stream:
        for row in rows:
            stream.write(json.dumps(row, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
