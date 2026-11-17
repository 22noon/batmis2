#!/usr/bin/env python3
"""manifest.py PREFIX OPT: prints a JSON record of the benchmark's build environment.

Records the machine (CPU, cores, memory, kernel, THP), the compiler, the flags, the source
versions, and for every binary under PREFIX/bin its size and how many popcnt and prefetch
instructions it holds (from objdump), so that a reader can check how each tool was built.
"""
import json, os, platform, re, socket, subprocess, sys


def run(cmd):
    try:
        return subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=600).stdout.strip()
    except Exception as e:  # noqa: BLE001 (a manifest field, never fatal)
        return f"error: {e}"


def count_insns(path):
    """popcnt and prefetch* instructions in an ELF binary (0, 0 for scripts)."""
    if run(f"file -b {path!r}").find("ELF") < 0:
        return 0, 0
    out = subprocess.run(["objdump", "-d", "--no-show-raw-insn", path], capture_output=True, text=True).stdout
    pop = len(re.findall(r"\tpopcnt\b", out))
    pre = len(re.findall(r"\tprefetch", out))
    return pop, pre


def main():
    prefix, opt = sys.argv[1], sys.argv[2]
    binaries = {}
    for tool in sorted(os.listdir(os.path.join(prefix, "bin"))):
        d = os.path.join(prefix, "bin", tool)
        for f in sorted(os.listdir(d)):
            p = os.path.join(d, f)
            if os.path.isfile(p) and os.access(p, os.X_OK):
                pop, pre = count_insns(p)
                binaries[f"{tool}/{f}"] = {"bytes": os.path.getsize(p), "popcnt": pop, "prefetch": pre}
    sources = {}
    src = os.path.join(prefix, "src")
    for d in sorted(os.listdir(src)) if os.path.isdir(src) else []:
        if os.path.isdir(os.path.join(src, d, ".git")):
            sources[d] = run(f"git -C {os.path.join(src, d)!r} describe --tags --always")
    repo = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
    sources["batmis"] = run(f"git -C {repo!r} rev-parse HEAD") + (" (modified)" if run(f"git -C {repo!r} status --porcelain --untracked-files=no") else "")
    thp = run("cat /sys/kernel/mm/transparent_hugepage/enabled")
    print(json.dumps({
        "host": socket.gethostname(),
        "cpu": run("lscpu | grep 'Model name' | sed 's/.*: *//'"),
        "lscpu": run("lscpu"),
        "cores": os.cpu_count(),
        "memory": run("free -g | awk '/Mem:/{print $2\" GB\"}'"),
        "kernel": platform.release(),
        "thp": thp,
        "compiler": run("gcc --version | head -1"),
        "opt_flags": opt,
        "sources": sources,
        "binaries": binaries,
    }, indent=1))


if __name__ == "__main__":
    main()
