#!/bin/bash
#SBATCH --job-name=ff_residual_screen_r167
#SBATCH --output=ff_residual_screen_r167_%j.out
#SBATCH --error=ff_residual_screen_r167_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=54
#SBATCH --time=01:00:00
#SBATCH --partition=general

# Runs probe_ff_residual_screen_geomsweep.py for RATIO_TAG=r167 (54
# cut-off-adjacent candidates: r0_2b=1.5, all 6 angles, both parts). See
# that file's docstring for the full context, target-list provenance, and
# the HONEST SCOPE CAVEAT about applying FF-P1-calibrated thresholds to a
# new geometry/material for the first time.
#
# Cost estimate: 54 points x ~150s/point serial (1 full_search + 1
# residual-screen build + 1 cheap sigma_min reuse, dps=40, mirroring
# probe_ff_oop_spurious_table.py's own per-point cost) = ~2.25 CPU-hours,
# parallelized 1 point/worker across 54 workers -- expect a few minutes
# wall time plus queue overhead; --time=01:00:00 is a generous margin.

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export BLIS_NUM_THREADS=1

export PKG_PATH=/home/ghmkfh/PythonMill/Plate_Solver_Package
export DPS=40
export EXPECT_SOLVER_VERSION=2026-07-10.s8
export RATIO_TAG=r167
export N_WORKERS=${SLURM_CPUS_PER_TASK:-53}

module load python/3.12.1
source /home/ghmkfh/PythonMill/Plate_Solver_Package/venv/bin/activate

# NOTE: no explicit `cd` here on purpose -- PKG_PATH (absolute) is what
# lets the probe find the plate_solver package regardless of where this
# script/its probe live; submit this from the SAME directory as
# probe_ff_residual_screen_geomsweep.py (cd there, then `sbatch
# submit_ff_residual_screen_r167.sh`).
if [ ! -f probe_ff_residual_screen_geomsweep.py ]; then
  echo "FATAL: probe_ff_residual_screen_geomsweep.py not found in $(pwd)."
  exit 2
fi

echo "Job $SLURM_JOB_ID on $(hostname): cpus-per-task=$SLURM_CPUS_PER_TASK, N_WORKERS=$N_WORKERS, RATIO_TAG=$RATIO_TAG"
echo "Start: $(date)"

python3 -u probe_ff_residual_screen_geomsweep.py

echo "End: $(date)"
