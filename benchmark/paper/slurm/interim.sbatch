#!/bin/bash
#SBATCH --job-name=bench_interim
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
# An interim report while 3_runs jobs are still running: scores and merges the runs finished
# so far, into a separate results directory (the final results/ is not touched).
#   sbatch --export=ALL,BENCH_WORK=... interim.sbatch
#   SELECT (default genome=sarscov2|ecoli) limits the correctness scoring, the heavy step;
#   SELECT=none skips it. BENCH_RESULTS defaults to <work>/../interim_<date>.
# Only finished runs (run.json written) are read; the running jobs are not disturbed, and
# 3_runs.sbatch is --exclusive, so this job never shares their node.
# Submit from this directory: sbatch copies the script, so $0 is not its location.
export BENCH_RESULTS=${BENCH_RESULTS:-$(dirname "${BENCH_WORK:-/tmp/$USER/batmis_bench}")/interim_$(date +%m%d)}
source "${SLURM_SUBMIT_DIR:-$(dirname "$0")}/common.sh"
SELECT=${SELECT-genome=sarscov2|ecoli}
if [ "$SELECT" != none ]; then
    bench correct ${SELECT:+--set "$SELECT"}
    bench discrepancies ${SELECT:+--set "$SELECT"}
fi
bench collect
python3 "$PAPER/report.py" "$BENCH_RESULTS"
echo "interim report: $BENCH_RESULTS/report.md" >&2
