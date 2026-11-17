#!/bin/bash
#SBATCH --job-name=bench_build
#SBATCH --exclusive
#SBATCH --time=02:00:00
# Builds every tool on the node itself (-march=native), and writes tools/manifest.json.
# Needs ORIG_TARBALL (BatMis-3.00.tar.gz) and network access to github.com.
# Submit from this directory: sbatch copies the script, so $0 is not its location.
source "${SLURM_SUBMIT_DIR:-$(dirname "$0")}/common.sh"
node_report
PREFIX=$BENCH_WORK/tools "$PAPER/build_tools.sh"
g++ -O3 -march=native -std=c++17 -pthread -o "$BENCH_WORK/tools/brute_force" "$PAPER/brute_force.cpp"
stage_out tools
