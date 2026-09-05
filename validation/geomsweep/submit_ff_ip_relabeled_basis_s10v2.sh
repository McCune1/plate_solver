#!/bin/bash
#SBATCH --job-name=ip_relabel_basis_s10v2
#SBATCH --output=ff_ip_relabeled_basis_s10v2_%j.out
#SBATCH --error=ff_ip_relabeled_basis_s10v2_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=25
#SBATCH --time=01:00:00
#SBATCH --partition=general

# Runs probe_ff_ip_relabeled_basis_s10v2.py -- direct follow-up to job
# 2417234 (probe_ff_residual_screen_geomsweep_s10v2_reconstruction_fix.py),
# which closed all 12 OOP cut-off-adjacent candidates but only 7/32 IP.
# The remaining 25 IP candidates split into the SAME Group A (24 pts,
# stuck at exactly n_dofs-2 even at xmax_cap=90) / Group B (1 pt,
# Om=0.547221, full fill but shallow smin=-2.681) structure the pre-s10
# investigation already solved via probe_ff_ip_relabeled_basis.py (job
# 2327405): Group A gets evaluated at its achieved basis size (mirrors
# the paper's existing used@28/used@36 relabeling precedent); Group B
# gets re-converged via seed_polish_modes instead of more branch search.
# See the probe's docstring for full reasoning and pre-registered
# interpretation.
#
# MUST run from FutureWork/geometry_sweep/ -- reads figures_ff_{r150,
# r167,r200,r250}_s10v2/freefree_manifest.json directly (same manifests
# job 2417234 used); does NOT import TARGETS from any other probe script,
# unlike the pre-s10 version.
#
# Cost estimate: 24 Group-A points re-run the escalation ladder
# (deterministic, up to ~1200s/point per job 2417234's own log) plus a
# residual evaluation; 1 Group-B point costs a seed_polish_modes
# convergence (~150-250s) plus escalated derivation. Total on the order
# of ~8 CPU-hours serial; parallelized across 25 workers, expect roughly
# 20-25 min wall time if most points don't need the full guarded stage,
# up to the full hour if several do. --time=01:00:00 mirrors the pre-s10
# job's own budget for this exact cost profile.

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export BLIS_NUM_THREADS=1

export PKG_PATH=/home/ghmkfh/PythonMill/Plate_Solver_Package
export DPS=40
export EXPECT_SOLVER_VERSION=2026-07-10.s10
export N_WORKERS=${SLURM_CPUS_PER_TASK:-25}

module load python/3.12.1
source /home/ghmkfh/PythonMill/Plate_Solver_Package/venv/bin/activate

if [ ! -f probe_ff_ip_relabeled_basis_s10v2.py ]; then
  echo "FATAL: probe_ff_ip_relabeled_basis_s10v2.py not found in $(pwd)."
  echo "cd to FutureWork/geometry_sweep/ before sbatch."
  exit 2
fi
MISSING=0
for TAG in r150 r167 r200 r250; do
  if [ ! -f "figures_ff_${TAG}_s10v2/freefree_manifest.json" ]; then
    echo "FATAL: figures_ff_${TAG}_s10v2/freefree_manifest.json not found. If"
    echo "figures_ff_${TAG}_s10v2.zip is present instead, unzip it first:"
    echo "  mkdir -p figures_ff_${TAG}_s10v2 && unzip -o -j figures_ff_${TAG}_s10v2.zip freefree_manifest.json -d figures_ff_${TAG}_s10v2"
    MISSING=1
  fi
done
if [ "$MISSING" -eq 1 ]; then
  exit 2
fi

echo "Job $SLURM_JOB_ID on $(hostname): cpus-per-task=$SLURM_CPUS_PER_TASK, N_WORKERS=$N_WORKERS"
echo "Start: $(date)"

python3 -u probe_ff_ip_relabeled_basis_s10v2.py

echo "End: $(date)"
