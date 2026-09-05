#!/bin/bash
#SBATCH --job-name=ff_recon_fix
#SBATCH --output=ff_residual_screen_reconstruction_fix_%j.out
#SBATCH --error=ff_residual_screen_reconstruction_fix_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=44
#SBATCH --time=03:00:00
#SBATCH --partition=general

# Runs probe_ff_residual_screen_reconstruction_fix.py -- retries the 44
# cut-off-adjacent candidates (12 OOP, 32 IP) that jobs 2327073-2327076
# could not reconstruct, using the already-validated _derive_branches_
# robust escalation ladder (probe_threshold_sensitivity.py, job 2325369)
# instead of the original bare cold-search call. See that probe's
# docstring for the full context, the concrete grounded hypothesis, the
# in-run fidelity gate, and the pre-registered interpretation. This
# directly targets an item PAPER1_FREEFREE_DRAFT.tex Sec 6.4 currently
# states as "not root-caused."
#
# MUST run from the same directory as probe_ff_residual_screen_
# geomsweep.py (FutureWork/geometry_sweep/) -- this probe imports that
# file's TARGETS list verbatim rather than retyping 44 rows of
# floating-point data by hand. RATIO_TAG is exported below only to
# satisfy that imported module's own CLI guard (SystemExit if unset) --
# this probe's own logic does not use it; it always processes all 44
# failures from all 4 original ratios together in one run.
#
# Cost estimate: up to ~1200s/point worst case (a point that needs the
# guarded xmax_cap=90 escalation stage) x 44 points = ~14.7 CPU-hours
# serial; parallelized 1 point/worker across up to 44 workers, expect a
# few minutes if most resolve early/cheap, up to ~20-25 min if many need
# the expensive final stage. --time=03:00:00 is a generous margin.

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export BLIS_NUM_THREADS=1

export PKG_PATH=/home/ghmkfh/PythonMill/Plate_Solver_Package
export DPS=40
export EXPECT_SOLVER_VERSION=2026-07-10.s8
export RATIO_TAG=r150
export N_WORKERS=${SLURM_CPUS_PER_TASK:-44}

module load python/3.12.1
source /home/ghmkfh/PythonMill/Plate_Solver_Package/venv/bin/activate

if [ ! -f probe_ff_residual_screen_reconstruction_fix.py ]; then
  echo "FATAL: probe_ff_residual_screen_reconstruction_fix.py not found in $(pwd)."
  echo "This probe must be placed in FutureWork/geometry_sweep/ alongside"
  echo "probe_ff_residual_screen_geomsweep.py -- cd there before sbatch."
  exit 2
fi
if [ ! -f probe_ff_residual_screen_geomsweep.py ]; then
  echo "FATAL: probe_ff_residual_screen_geomsweep.py (source of the TARGETS"
  echo "list this probe imports) not found in $(pwd)."
  exit 2
fi

echo "Job $SLURM_JOB_ID on $(hostname): cpus-per-task=$SLURM_CPUS_PER_TASK, N_WORKERS=$N_WORKERS"
echo "Start: $(date)"

python3 -u probe_ff_residual_screen_reconstruction_fix.py

echo "End: $(date)"
