# Sourced by every sbatch script: the environment of the paper benchmark on the node.
#   BENCH_WORK     node-local scratch that holds tools, genomes, reads, indexes and runs
#                  (default /tmp/$USER/batmis_bench; it must survive between jobs on the
#                  reserved node, or be re-staged from BENCH_STORE)
#   BENCH_RESULTS  where the merged TSVs, figures and report go (shared filesystem)
#   BENCH_STORE    optional shared directory: stage_in copies tools/genomes/reads/index from
#                  it if they are missing, stage_out saves runs/ to it
#   GCC_MODULE     optional module to load for the builds (recorded in the manifest)
set -euo pipefail
REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)  # common.sh is sourced from its own path
PAPER=$REPO/benchmark/paper
export BENCH_WORK=${BENCH_WORK:-/tmp/$USER/batmis_bench}
export BENCH_RESULTS=${BENCH_RESULTS:-$PAPER/results}
CONFIG=${CONFIG:-$PAPER/config.yaml}
mkdir -p "$BENCH_WORK" "$BENCH_RESULTS"
if [ -n "${GCC_MODULE:-}" ] && command -v module >/dev/null; then module load "$GCC_MODULE"; fi
bench() { python3 "$PAPER/bench.py" -c "$CONFIG" "$@"; }

stage_in() {
    [ -n "${BENCH_STORE:-}" ] || return 0
    for d in tools genomes reads index truth; do
        if [ ! -e "$BENCH_WORK/$d" ] && [ -e "$BENCH_STORE/$d" ]; then
            echo "stage in $d" >&2
            cp -a "$BENCH_STORE/$d" "$BENCH_WORK/"
        fi
    done
}
stage_out() {
    [ -n "${BENCH_STORE:-}" ] || return 0
    mkdir -p "$BENCH_STORE"
    for d in "$@"; do rsync -a "$BENCH_WORK/$d" "$BENCH_STORE/"; done
}
node_report() {
    echo "node $(hostname), $(nproc) CPUs, $(lscpu | sed -n 's/^Model name: *//p'), load $(cut -d' ' -f1-3 /proc/loadavg)" >&2
}
