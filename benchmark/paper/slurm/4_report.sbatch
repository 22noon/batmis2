#!/bin/bash
#SBATCH --job-name=bench_report
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
# Scores the correctness runs, merges everything and writes the report and figures.
# Submit from this directory: sbatch copies the script, so $0 is not its location.
source "${SLURM_SUBMIT_DIR:-$(dirname "$0")}/common.sh"
stage_in
bench correct
bench collect
python3 "$PAPER/report.py" "$BENCH_RESULTS"
stage_out runs
