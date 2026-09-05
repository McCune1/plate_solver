#!/bin/bash
#SBATCH --job-name=ff_recon_fix_s10v2
#SBATCH --output=ff_res_screen_geomsweep_s10v2_reconstruction_fix_%j.out
#SBATCH --error=ff_res_screen_geomsweep_s10v2_reconstruction_fix_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=44
#SBATCH --time=03:00:00
#SBATCH --partition=general

# Runs probe_ff_residual_screen_geomsweep_s10v2_reconstruction_fix.py --
# retries the 44 cut-off-adjacent candidates (12 OOP, 32 IP) that jobs
# 2416866-2416869 could not reconstruct with a bare cold search, using
# the already-validated _derive_branches_robust escalation ladder
# (probe_threshold_sensitivity.py, job 2325369; previously used to close
# the identical 44-candidate failure under the pre-s10 checkpoint via
# probe_ff_residual_screen_reconstruction_fix.py). See that probe's
# docstring for the full context, the concrete grounded hypothesis, the
# in-run fidelity gate, and the pre-registered interpretation. This
# closes the "spot-check, escalation not re-run" caveat 2026-08-13 added
# to PAPER1_FREEFREE_DRAFT.tex's Table 13 caption and Appendix C.3
# (LESSONS_LEARNED.md Sec 18.10).
#
# MUST run from FutureWork/geometry_sweep/ -- reads figures_ff_{r150,
# r167,r200,r250}_s10v2/freefree_manifest.json as relative paths from
# cwd (the same manifests jobs 2416866-2416869 used). Processes all 44
# failures from all 4 ratios together in ONE run (unlike the original
# geomsweep_s10v2 probe, which was split per-ratio) -- matches the
# pre-s10 reconstruction-fix job's own convention.
#
# Cost estimate: up to ~1200s/point worst case (a point needing the
# guarded xmax_cap=90 escalation stage) x 44 points = ~14.7 CPU-hours
# serial; parallelized 1 point/worker across up to 44 workers, expect a
# few minutes if most resolve early/cheap, up to ~20-25 min if many need
# the expensive final stage. --time=03:00:00 is a generous margin,
# identical to the pre-s10 job's own budget for this exact cost profile.

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export BLIS_NUM_THREADS=1

export PKG_PATH=/home/ghmkfh/PythonMill/Plate_Solver_Package
export DPS=40
export EXPECT_SOLVER_VERSION=2026-07-10.s10
export N_WORKERS=${SLURM_CPUS_PER_TASK:-44}

module load python/3.12.1
source /home/ghmkfh/PythonMill/Plate_Solver_Package/venv/bin/activate

if [ ! -f probe_ff_residual_screen_geomsweep_s10v2_reconstruction_fix.py ]; then
  echo "FATAL: probe_ff_residual_screen_geomsweep_s10v2_reconstruction_fix.py not found in $(pwd)."
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

python3 -u probe_ff_residual_screen_geomsweep_s10v2_reconstruction_fix.py

echo "End: $(date)"
