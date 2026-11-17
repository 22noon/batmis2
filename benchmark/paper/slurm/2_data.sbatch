#!/bin/bash
#SBATCH --job-name=bench_data
#SBATCH --exclusive
#SBATCH --time=12:00:00
# Genomes, simulated and real read sets, every index (timed), and the brute-force truth.
# Submit from this directory: sbatch copies the script, so $0 is not its location.
source "${SLURM_SUBMIT_DIR:-$(dirname "$0")}/common.sh"
stage_in
node_report
bench fetch
bench simulate
bench prepare
bench index
bench truth
stage_out genomes reads index truth runs
