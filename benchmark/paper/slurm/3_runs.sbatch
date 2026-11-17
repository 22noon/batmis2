#!/bin/bash
#SBATCH --job-name=bench_runs
#SBATCH --exclusive
#SBATCH --time=7-00:00:00
# The run matrix, or a part of it: sbatch --export=ALL,SELECT='kind=timing,genome=human' 3_runs.sbatch
# Runs are resumable: resubmit with the same SELECT to continue after a time limit.
# Submit from this directory: sbatch copies the script, so $0 is not its location.
source "${SLURM_SUBMIT_DIR:-$(dirname "$0")}/common.sh"
stage_in
node_report
bench list ${SELECT:+--set "$SELECT"} | tail -1
bench run ${SELECT:+--set "$SELECT"}
stage_out runs
