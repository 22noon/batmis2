#!/bin/bash
#SBATCH --job-name=bench_pilot
#SBATCH --exclusive
#SBATCH --time=12:00:00
# The pilot (config_pilot.yaml): SARS-CoV-2 and E. coli, small read sets, every step from
# download to report, to check the toolchain before the real runs. Submit from this
# directory after 1_build.sbatch:
#   sbatch --reservation=<name> --mem=0 0_pilot.sbatch
# It uses the tools in $BENCH_WORK/tools, works in $BENCH_WORK/pilot, and writes the report
# to benchmark/paper/results_pilot (or $BENCH_RESULTS_PILOT).
set -euo pipefail
HERE=${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")" && pwd)}
MAIN_WORK=${BENCH_WORK:-/tmp/$USER/batmis_bench}
if [ ! -d "$MAIN_WORK/tools/bin" ]; then
    echo "no tools in $MAIN_WORK/tools: run 1_build.sbatch first (with the same BENCH_WORK)" >&2
    exit 1
fi
export BENCH_WORK=$MAIN_WORK/pilot
export BENCH_RESULTS=${BENCH_RESULTS_PILOT:-$HERE/../results_pilot}
export CONFIG=$HERE/../config_pilot.yaml
mkdir -p "$BENCH_WORK"
ln -sfn "$MAIN_WORK/tools" "$BENCH_WORK/tools"
source "$HERE/common.sh"
node_report
bench fetch
bench simulate
bench prepare
bench index
bench truth
bench run
bench correct
bench discrepancies
bench collect
python3 "$PAPER/report.py" "$BENCH_RESULTS"
echo "pilot report: $BENCH_RESULTS/report.md" >&2
