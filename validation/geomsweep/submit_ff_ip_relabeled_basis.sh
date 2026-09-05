#!/bin/bash
#SBATCH --job-name=ip_relabel_basis
#SBATCH --output=ff_ip_relabeled_basis_%j.out
#SBATCH --error=ff_ip_relabeled_basis_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=25
#SBATCH --time=01:00:00
#SBATCH --partition=general

# Runs probe_ff_ip_relabeled_basis.py -- direct follow-up to job 2327101.
# Group A (24 pts): the IP candidates stuck at exactly n_dofs-2 branches
# even at xmax_cap=90 -- evaluates the residual screen at the achieved
# basis size instead of continuing to demand the original target (mirrors
# the paper's existing used@28/used@36 relabeling precedent). Group B (1
# pt, Om=0.547282): the one shallow-lock-on case -- re-converges via
# seed_polish_modes instead of more branch-search escalation. See the
# probe's docstring for full reasoning and pre-registered interpretation.
#
# MUST run from FutureWork/geometry_sweep/ (imports TARGETS from
# probe_ff_residual_screen_geomsweep.py). RATIO_TAG is exported only to
# satisfy that imported module's CLI guard, unused by this probe directly.
#
# Cost estimate: 24 Group-A points re-run job 2327101's own escalation
# ladder (deterministic, ~1150-1300s/point per that job's log) plus a
# residual evaluation; 1 Group-B point costs a seed_polish_modes
# convergence (~150-250s) plus escalated derivation. Total ~8.5 CPU-hours
# serial; parallelized across 25 workers, expect ~20-25 min wall time.
# --time=01:00:00 is generous headroom.

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
export N_WORKERS=${SLURM_CPUS_PER_TASK:-25}

module load python/3.12.1
source /home/ghmkfh/PythonMill/Plate_Solver_Package/venv/bin/activate

if [ ! -f probe_ff_ip_relabeled_basis.py ]; then
  echo "FATAL: probe_ff_ip_relabeled_basis.py not found in $(pwd)."
  echo "Must run from FutureWork/geometry_sweep/, alongside"
  echo "probe_ff_residual_screen_geomsweep.py."
  exit 2
fi
if [ ! -f probe_ff_residual_screen_geomsweep.py ]; then
  echo "FATAL: probe_ff_residual_screen_geomsweep.py (TARGETS source) not found."
  exit 2
fi

echo "Job $SLURM_JOB_ID on $(hostname): cpus-per-task=$SLURM_CPUS_PER_TASK, N_WORKERS=$N_WORKERS"
echo "Start: $(date)"

python3 -u probe_ff_ip_relabeled_basis.py

echo "End: $(date)"
