# Paper benchmark: batmis1 and batmis2 against established mappers

This directory reproduces the paper's comparison: correctness, time to the first hit and to
all hits within K mismatches, align/decode split, memory, index build, and thread scaling.

## Tools

**Exact** (every hit within K mismatches, Hamming distance):

| tool | version | first hit | all hits |
|---|---|---|---|
| batmis1, batmis2 | this repository | `-n K -m 1` | `-n K -m 0` |
| BatMis 3.00 (original) | `BatMis-3.00.tar.gz` | `batman -m 1` + `batdecode` | `batman -mall` + `batdecode` |
| BWA aln | 0.7.19 | `aln -n K -o 0 -e 0 -l 1000 -k K` + `samse` | adds `-N -m 20000000`, `samse -n 1000000` |
| Bowtie | 1.3.1 | `-v K -k 1 --best` (K ≤ 3) | `-v K -a` |
| Columba | 2.0 (Vanilla) | `-m hamming -a best -I <identity>` | `-m hamming -a all -e K` |
| mrsFAST-Ultra | 3.4.2 | `-e K --best` | `-e K` |

**Approximate** (heuristic, gapped; recovery and speed only): BWA-MEM 0.7.19 (default, `-a`)
and Bowtie 2 2.5.4 (`--end-to-end --very-sensitive`, default and `-k 100`).

Notes for the methods section:
- **BWA aln is made exhaustive.** Seeding is disabled (`-l` longer than the read), gaps are
  off (`-o 0 -e 0`) and `-N` turns off the iterative search, so it reports every hit within K.
  That is not how it is normally run; it is the only setting that makes it comparable.
- **BWA is the unmodified release**, built with `-O3 -march=native -mpopcnt`. It counts bases
  with a bit trick and a byte table rather than a popcount call, so the flag may add no
  popcount instruction; the manifest shows the count. `bwa_aln_stock` (BWA's own Makefile
  flags) is timed at K = 0, 2, 4 as a supplement.
- **Columba has no one-hit mode.** Its best mode reports every hit of the best stratum,
  with a whole-number identity threshold.
- **Original BatMis** is single-threaded. One guard is added to its `batdecode`
  (`build_tools.sh`): `Location_To_Genome` reads the key of `map::end()` for every hit in the
  reference's last sequence (all hits, in a one-sequence genome), which is undefined
  behaviour; optimised GCC 13 builds crashed on it. The guard removes the read and changes
  no result.
- **`bwa samse` is single-threaded**, which the align/decode split shows.
- **Overlapping occurrences in tandem repeats:** on chr21 one 40 bp read has 992 hits at
  K = 4, most of them overlapping (a (TA)n microsatellite). Columba reported 529 and exhaustive
  BWA aln 991; batmis1/2, the original, Bowtie and mrsFAST all 992. Columba appears to report
  only some of a set of overlapping occurrences. Expect such "missed" hits in repeats, and
  name them as such in the paper.
- **mrsFAST-Ultra** is exact only while the read holds K + 1 whole hash windows (12 bases by
  default): at 40 bp up to K = 2. Beyond that it misses hits by design.
- **Reference N:** BatMis stores it as C; BWA, Columba and mrsFAST put some base there. Hits
  whose window covers a reference N are left out of every comparison (`over_n`).

## Data

- **Genomes:** SARS-CoV-2 (NC_045512.2), E. coli K-12 MG1655 (NC_000913.3), human GRCh38
  no-alt analysis set; chr21 for the human pipeline's pilot.
- **Read lengths:** 40, 100 and 150 bp for SARS-CoV-2 and E. coli; **100 bp only for human**
  (every human set: timing, real data, correctness, realistic), which keeps the human part of
  the matrix, by far the most expensive, to a third (`lengths_by_genome`, `real.hg002.lengths`).
- **Simulated** (`simulate_reads.py`), both strands.
  - `timing`: 1M stratified reads, exactly k = 0..7 substitutions in equal shares, plus 10%
    random reads;
  - `check`: 100k of the same kind, scored against brute force (small genomes) or the verified
    union of the exact tools (human);
  - `real`: a position-dependent substitution profile (0.1% to 1.5% along the read);
  - `indel`: the same with 0.1% indels, for the approximate tools.
- **Real** (`prepare_reads.py`), single-end R1, trimmed to each length, reads with N dropped:
  - SARS-CoV-2 ARTIC amplicons, ERR4207974 (COG-UK, NextSeq 550);
  - GIAB HG002, SRR15395856 (NovaSeq 6000, 2×151), about 10M reads by a name-hash sample.
- 1-thread runs on real data use the first 1M reads of each set (`single_thread_reads`).

## Running

Needs: GCC with C++17, CMake, git, curl, zlib, groff (mrsFAST's help page), Python 3 with
numpy, PyYAML and matplotlib, and network access to github.com, NCBI and ENA for the build
and data jobs.

Validated on a workstation (`config_pilot.yaml`, and a chr21 version of `config.yaml`):
every exact tool's hits against brute force for SARS-CoV-2, E. coli and chr21 at 40, 100 and
150 bp, K = 0-7; and the human method (the verified union of the exact tools' hits) gave
exactly the brute-force hit sets on chr21 at every length and K.

Everything is driven by `bench.py` and `config.yaml` (`config_pilot.yaml` is a workstation-
sized version). On a SLURM node reserved with `--exclusive`, from `benchmark/paper/slurm/`:

```sh
export BENCH_WORK=/local/scratch/$USER/batmis_bench   # node-local disk
export BENCH_STORE=/shared/$USER/batmis_bench         # optional: survives the node
export ORIG_TARBALL=/path/to/BatMis-3.00.tar.gz
sbatch 1_build.sbatch                                 # tools, built on the node, manifest
sbatch 0_pilot.sbatch                                 # optional: the pilot on the small genomes (about 1 h)
sbatch 2_data.sbatch                                  # genomes, reads, indexes, truth, calibration
sbatch --export=ALL,SELECT='genome=sarscov2|ecoli' 3_runs.sbatch
sbatch --export=ALL,SELECT='genome=human,kind=check|load' 3_runs.sbatch
sbatch --export=ALL,SELECT='genome=human,kind=timing|realdata' 3_runs.sbatch
sbatch --export=ALL,SELECT='kind=scaling' 3_runs.sbatch
sbatch 4_report.sbatch                                # scores, TSVs, report.md, figures
```

Run by hand, `bench.py` needs `BENCH_WORK` set as for the jobs; `BENCH_RESULTS` defaults to
`benchmark/paper/results`. The same steps by hand: `bench.py fetch | simulate | prepare | index | truth | run | correct |
collect`, then `report.py <results>`. `bench.py list --set ...` shows the matrix and what is
done; `bench.py run --set ... --dry-run` prints the commands.

While jobs run: `bench.py status --set ...` gives done/timeout/error counts per tool, the newest
finished run, and the runs in progress with the age of their phase (one over its time limit is a
hung process, or was left by a job that stopped). A run writes `run.json` only when it finishes,
so a job can be cancelled and resubmitted with the same `SELECT`: only the run in progress is
redone. Never run two jobs whose `SELECT`s overlap; they delete each other's run directories.
`sbatch interim.sbatch` writes an interim report from the finished runs into
`interim_<date>/` beside the work directory (`SELECT` limits its scoring, default the small
genomes); `bench.py collect` plus `report.py` by hand does the same without the scoring.

Run behaviour:
- Each run is one directory under `<work>/runs/<kind>/<read set>/<tool>/k<K>/<mode>/t<threads>/r<rep>`
  with `run.json`; existing runs are skipped, so jobs can be resubmitted after a time limit.
- `rep` 0 is a warm-up (page cache); the tables use the median of `repeats` timed runs.
  A warm-up longer than 30 min gets one timed run.
- A run over `time_limit_s` (6 h) is a timeout, and the larger K of that series are skipped.
- Each phase runs under `/usr/bin/time -v` (wall, CPU, peak memory); batmis also reports its
  own phases (`--timing`).
- Output SAM goes to node-local disk (`local_output`, default `/tmp/$USER/batmis_bench_out`) and
  is deleted after the run; correctness runs keep their hits, normalised (`normalize_hits.py`),
  in the run directory. Inputs, indexes and results can sit on a shared filesystem.

## BatMis with and without the reference

batmis1 and batmis2 normally load the packed genome (`.pac`) and finish a single-row search
path against it ("reference loaded"). With `--no-ref` they walk the path through the index to
its end instead: less memory, more index steps. Every cross-tool comparison uses the
reference-loaded runs; `batmis1_noref` and `batmis2_noref` repeat batmis's whole timing
matrix (every timing and real-data set, all K, both modes, both thread counts) and its
correctness and load runs, with `--no-ref` added. Both write SAM, which needs the SA samples
(`.sa`) either way, so the memory difference is the `.pac` file.

`report.py` draws `batmis_ref_vs_noref_*` (time against K, reference loaded solid, without
dashed) and tabulates wall time, batmis's own search and index-load times, and peak memory
both ways, with the ratio.

## Human runs: sized by calibration

Human runs with many hits per read or high K take hours at the full read sets (1M simulated,
about 10M HG002 reads), so their read counts are set from a calibration (`calibration:` in
`config.yaml`):

1. **Calibrate** (`bench.py calibrate`, run by `2_data.sbatch`): every tool, both modes, every
   human timing and real-data set, timed on 5,000 reads (a fixed-seed random sample that keeps
   the set's k mix) on one thread, at K = 0 and K = 5 (Bowtie: 3). Start-up is taken from the
   one-read load runs and subtracted; the result is seconds per read at each K, written to
   `runs/calib/plan.json`.
2. **Size each run:** per-read time at its K (log-linear between and beyond the two points),
   divided by its threads (assumed linear), gives the reads that fit 10 minutes; rounded down
   to 5k, 10k, 20k, 50k, ... 10M, at least 5,000, at most the whole set. `bench.py list` shows
   each run's count; `bench.py run` holds human timing runs back until the calibration exists.
   `bench.py calibrate` (and `bench.py estimate` any time later) prints the estimated node time
   of the sized runs still to do, by kind, mode and threads: lower `target_s` or `repeats` in
   `config.yaml` if it is too much, then run `estimate` again (no recalibration needed).
3. **Compare per read:** every run records `reads_n` and its start-up, and `timing.tsv` gives
   `us_per_read` and `s_per_million` (seconds per million reads, start-up included); the report
   plots human results in seconds per million reads. A run more than 3x (and 10 s) over its
   estimate is noted in the run log.

The human correctness sets (and the realistic sets for true-origin recovery) are 5,000 reads
(`correctness_reads_by_genome`), the same reads for every tool, so the verified union of hits
still applies. Existing larger sets are sampled down; nothing is simulated again.

In the pilot (E. coli calibrated, with a tiny target to force small read counts) runs took a
median 1.1x their estimate; BatMis runs under a second on 4 threads took up to 15x, from fixed
costs that only matter at that scale.

## Analysing discrepancies

`bench.py correct` counts each tool's missed and extra hits; `bench.py discrepancies` says why.
For every correctness set it writes `<results>/discrepancies/<set>/`:

- `summary.tsv`: the discrepancies of each tool, K and mode, by cause;
- `pairwise.tsv`: for each K, the valid hits tool A found and tool B did not;
- `<tool>_k<K>_<mode>.tsv`: every discrepancy, with its read and the read's truth-hit count;
- `examples.txt`: up to 10 per cause, aligned against the genome.

Causes (`discrepancies.py` has the exact rules):
- extra: `over_n` (window over reference N), `clipped_or_gapped`, `off_contig`, `over_k`
  (more than K mismatches), `not_in_truth`;
- missed: `window_limit` (mrsFAST, L < (K+1)·12), `overlapping` (the tool reports a hit of the
  read within L bases: tandem repeats), `many_hits_read` (> 100 truth hits, some reported),
  `contig_edge`, `low_complexity`, `unexplained`;
- first-hit runs: `missed_read` (with its best hit's cause) and `not_best`.

`unexplained` rows are the ones to look into by hand. One read at a time:

```sh
bench.py discrepancies --set dataset=sim_chr21_40_check --read 'sim4579_k3_chr21_23204462_+'
```

This prints the read, its true origin, each run's hit count against the truth, and the truth
hits aligned to the genome with the runs that miss each. In the pilot every discrepancy had a
cause (none `unexplained`); `report.py` lists the non-zero causes in `report.md`.

## Before the main runs

1. Check `tools/manifest.json`: popcnt and prefetch counts per binary (batmis1/2 must pass
   `check_prefetch.sh`, which `build_tools.sh` runs).
2. Run batmis2 on the human timing set with `--batch auto` a few times, and check its
   `batch auto:` line: if 16 is not within a few percent of the fastest, set `BATMIS_BATCH`.
3. Estimate node-hours from the small-genome runs and trim `config.yaml` if needed
   (`k_exact`, `threads`, `repeats`, read counts).

## Files

| file | role |
|---|---|
| `config.yaml`, `config_pilot.yaml` | the matrix, data and paths |
| `build_tools.sh`, `manifest.py` | builds every tool with one compiler and `-O3 -march=native`; manifest |
| `tools.py` | each tool's index and mapping commands |
| `bench.py` | the driver |
| `simulate_reads.py`, `prepare_reads.py` | read sets |
| `normalize_hits.py` | SAM (and BWA's XA) to one hit list |
| `brute_force.cpp` | ground truth by scanning the genome |
| `compare_hits.py` | correctness scores |
| `discrepancies.py` | the cause of every missed or extra hit; one-read view |
| `collect.py`, `report.py` | TSVs, report.md, figures |
| `slurm/` | the jobs |
