#!/bin/bash
#SBATCH --job-name=ff_adjudicate_r150_a100_nu030_s10readj
#SBATCH --output=ff_adjudicate_r150_a100_nu030_s10readj_%j.out
#SBATCH --error=ff_adjudicate_r150_a100_nu030_s10readj_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=64
#SBATCH --time=03:00:00
#SBATCH --requeue
#SBATCH --partition=requeue

# Re-adjudication rerun of Sec 6.5 / Table S.4's r150/a100 nu=0.30 result
# (original job 2381555) under the CURRENT SOLVER_VERSION=2026-07-10.s10
# checkpoint, per FUTURE_CHAT_PROMPTS_2026-08-13.md item 4. Reuses the
# SAME generic probe_ff_adjudicate_candidates.py already at the package
# root (restored 2026-08-08 for r200_a100, geometry-agnostic, no
# per-geometry fork needed) -- just pointed at the freshly-merged r150
# candidates file.
#
# TWO differences from the archived original submit_ff_adjudicate_r150_a100_nu030.sh:
#   1. --partition=requeue in place of the now-billed --partition=general
#      (2026-08-04 policy, project memory cluster-requeue-partition-billing).
#   2. ANSYS_DIR points at Ansys/NewAnsys (where
#      geomsweep_{oop,ip}_r150_a100_nu030.txt actually live today and where
#      the r200 pipeline's own ANSYS_DIR points), NOT the original's
#      FutureWork/ansys_geomsweep_verification path -- per this task's own
#      instruction. Confirm both .txt files are present there before
#      submitting; they should already be (generated during the original
#      r150_a100_nu030 closure).
#
# 3 things must be done first:
#   1. submit_ff_r150_a100_capture_nu030_oop_s10readj.sh + _ip_s10readj.sh
#      (the fresh scan)
#   2. submit_merge_r150_a100_nu030_s10readj.sh (combines the two capture
#      JSONs into r150_a100_nu030_candidates.json)
#   3. Confirm Ansys/NewAnsys/geomsweep_{oop,ip}_r150_a100_nu030.txt exist
#      (they already should -- generated for the original closure, see
#      generate_r150_a100_nu030_ansys_decks.py in the archive folder).
#
# cpus-per-task=64: the fresh capture (jobs 2419785 OOP / 2419786 IP,
# 2026-08-14) reproduced the ORIGINAL job 2381555 candidate set exactly --
# 30 OOP + 64 IP = 94 candidates, bit-identical Omega values to the Aug 3
# run (job 2362730/2362721) -- a clean reproducibility confirmation for the
# capture/scan step on its own. The original one-candidate-per-worker
# sizing (cpus-per-task=94) hit this cluster's evident per-job core ceiling
# (~64), so this is capped at 64 instead, with the time budget widened
# accordingly (02:00:00 -> 03:00:00). Per the cost model
# (plate-solver-cluster-probes skill: ~150s/point serial, ~150-490s/point
# observed in the original run's own log), 94 candidates over 64 workers is
# at most 2 dispatch rounds (ProcessPoolExecutor backfills a slot as soon
# as it frees, not strict batches), so worst case is roughly 2x the
# slowest single candidate (~500s) plus overhead -- expect well under
# 30 min in practice; the original 94-worker run finished all candidates
# in 9.9 min. 03:00:00 is generous headroom, not an estimate of actual
# runtime.
#
# Run from the SAME folder as r150_a100_nu030_candidates.json (i.e.
# wherever submit_merge_r150_a100_nu030_s10readj.sh wrote it -- your main
# PythonMill/Plate_Solver_Package/ folder).

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export BLIS_NUM_THREADS=1

export PKG_PATH=/home/ghmkfh/PythonMill/Plate_Solver_Package
export DPS=40
export EXPECT_SOLVER_VERSION=2026-07-10.s10
export CANDIDATES_JSON=$(pwd)/r150_a100_nu030_candidates.json
export ANSYS_DIR=/home/ghmkfh/PythonMill/Plate_Solver_Package/Ansys/NewAnsys
export FE_MATCH_TOL_PCT=3.0
export N_WORKERS=${SLURM_CPUS_PER_TASK:-64}

module load python/3.12.1
source /home/ghmkfh/PythonMill/Plate_Solver_Package/venv/bin/activate

if [ ! -f probe_ff_adjudicate_candidates.py ]; then
  echo "FATAL: probe_ff_adjudicate_candidates.py not found in $(pwd)."
  exit 2
fi
if [ ! -f r150_a100_nu030_candidates.json ]; then
  echo "FATAL: r150_a100_nu030_candidates.json not found in $(pwd)."
  echo "  Run submit_merge_r150_a100_nu030_s10readj.sh first (after both"
  echo "  nu=0.30 capture jobs complete)."
  exit 2
fi

echo "Job $SLURM_JOB_ID on $(hostname): cpus-per-task=$SLURM_CPUS_PER_TASK, N_WORKERS=$N_WORKERS"
echo "Start: $(date)"

python3 -u probe_ff_adjudicate_candidates.py

echo "End: $(date)"
