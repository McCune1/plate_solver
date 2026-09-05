#!/bin/bash
#SBATCH --job-name=ff_adjudicate_r150_a100_nu035_s10
#SBATCH --output=ff_adjudicate_r150_a100_nu035_s10_%j.out
#SBATCH --error=ff_adjudicate_r150_a100_nu035_s10_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=64
#SBATCH --time=03:00:00
#SBATCH --requeue
#SBATCH --partition=requeue

# Genuinely NEW ν=0.35 one-to-one adjudication for r150/a100 under
# SOLVER_VERSION=2026-07-10.s10, at production n_dofs=20 -- the original
# job 2347329 used the coarse geometry-sweep's own smaller-basis
# (n_dofs=16/14) accepted lists, not a fresh production capture. Reuses
# the same generic probe_ff_adjudicate_candidates.py already at the
# package root.
#
# cpus-per-task=64 (not one-worker-per-candidate): this cluster's actual
# per-job core ceiling is ~64, learned from the nu=0.30 rerun
# (submit_ff_adjudicate_r150_a100_nu030_s10readj.sh originally requested
# 94 and had to be cut back). --time widened to 03:00:00 accordingly; the
# nu=0.30 rerun's own 64-worker adjudicate finished in 8.7 min against 94
# candidates, so this should land in a similar range even if the nu=0.35
# candidate count differs somewhat from the nu=0.30 count (30 OOP + 64 IP).
#
# ANSYS_DIR points at FutureWork/ansys_geomsweep_verification -- NOT
# Ansys/NewAnsys (that's where the nu=0.30 FE files live; nu=0.35's FE
# data was never regenerated there, it's still only in the original
# geometry-sweep location). Confirm geomsweep_{oop,ip}_r150_a100.txt
# (no "_nu030" suffix) are present there before submitting -- they should
# be, unchanged since the original nu=0.35 closure.
#
# 3 things must be done first:
#   1. submit_ff_r150_a100_capture_nu035_oop.sh + _ip.sh (the fresh scan)
#   2. submit_merge_r150_a100_nu035.sh (combines the two capture JSONs
#      into r150_a100_nu035_candidates.json)
#   3. Confirm FutureWork/ansys_geomsweep_verification/
#      geomsweep_{oop,ip}_r150_a100.txt exist (they already should).
#
# Run from the SAME folder as r150_a100_nu035_candidates.json.

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export BLIS_NUM_THREADS=1

export PKG_PATH=/home/ghmkfh/PythonMill/Plate_Solver_Package
export DPS=40
export EXPECT_SOLVER_VERSION=2026-07-10.s10
export CANDIDATES_JSON=$(pwd)/r150_a100_nu035_candidates.json
export ANSYS_DIR=/home/ghmkfh/PythonMill/Plate_Solver_Package/FutureWork/ansys_geomsweep_verification
export FE_MATCH_TOL_PCT=3.0
export N_WORKERS=${SLURM_CPUS_PER_TASK:-64}

module load python/3.12.1
source /home/ghmkfh/PythonMill/Plate_Solver_Package/venv/bin/activate

if [ ! -f probe_ff_adjudicate_candidates.py ]; then
  echo "FATAL: probe_ff_adjudicate_candidates.py not found in $(pwd)."
  exit 2
fi
if [ ! -f r150_a100_nu035_candidates.json ]; then
  echo "FATAL: r150_a100_nu035_candidates.json not found in $(pwd)."
  echo "  Run submit_merge_r150_a100_nu035.sh first (after both nu=0.35"
  echo "  capture jobs complete)."
  exit 2
fi

echo "Job $SLURM_JOB_ID on $(hostname): cpus-per-task=$SLURM_CPUS_PER_TASK, N_WORKERS=$N_WORKERS"
echo "Start: $(date)"

python3 -u probe_ff_adjudicate_candidates.py

echo "End: $(date)"
