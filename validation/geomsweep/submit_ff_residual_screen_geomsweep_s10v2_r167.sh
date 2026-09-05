#!/bin/bash
#SBATCH --job-name=ff_res_screen_geomsweep_s10v2_r167
#SBATCH --output=ff_res_screen_geomsweep_s10v2_r167_%j.out
#SBATCH --error=ff_res_screen_geomsweep_s10v2_r167_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=64
#SBATCH --time=01:00:00
#SBATCH --partition=requeue

# Runs probe_ff_residual_screen_geomsweep_s10v2.py for RATIO_TAG=r167
# (r0/2b=1.66667). Re-screens the FULL cut-off-adjacent candidate population
# (the same population PAPER1_FREEFREE_DRAFT.tex's Table 13/Appendix C.3
# call "119 flexural, 96 extensional candidates") against the s10v2
# checkpoint -- NOT the same job as
# submit_ff_residual_screen_s10v2_growth_r167.sh, which only screened the
# DELTA between the old and new candidate sets (new-only + cross-geometry
# suspects). See the probe's own docstring and LESSONS_LEARNED.md Sec 18.9
# for the full context: this is the one piece of the s10v2 geomsweep
# integration still open as of 2026-08-13, explicitly flagged in Table
# 13's own caption.
#
# Cost: dry-run target discovery (sandbox, 2026-08-13, against the actual
# figures_ff_r167_s10v2.zip manifest) found 54 targets (30 OOP +
# 24 IP) -- identical in count to the pre-s10 script's own r167 total,
# as expected (n_cutoffs is a geometry property, not a candidate-set
# property). ~150s/point serial (mirrors both prior geomsweep-residual
# probes' own estimate) = ~2.25 CPU-hours, parallelized 1 point/worker
# across 64 workers -- expect well under 10 min wall time plus queue
# overhead; --time=01:00:00 is a generous margin.

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export BLIS_NUM_THREADS=1

export PKG_PATH=/home/ghmkfh/PythonMill/Plate_Solver_Package
export DPS=40
export EXPECT_SOLVER_VERSION=2026-07-10.s10
export RATIO_TAG=r167
export N_WORKERS=${SLURM_CPUS_PER_TASK:-64}

module load python/3.12.1
source /home/ghmkfh/PythonMill/Plate_Solver_Package/venv/bin/activate

# NOTE: no explicit `cd` to PKG_PATH -- PKG_PATH (absolute) is what lets
# the probe find the plate_solver package. Submit this from the SAME
# directory as probe_ff_residual_screen_geomsweep_s10v2.py, i.e.
# FutureWork/geometry_sweep/ (cd there, then
# `sbatch submit_ff_residual_screen_geomsweep_s10v2_r167.sh`) -- the
# script reads figures_ff_r167_s10v2/freefree_manifest.json as a
# RELATIVE path from cwd.
if [ ! -f probe_ff_residual_screen_geomsweep_s10v2.py ]; then
  echo "FATAL: probe_ff_residual_screen_geomsweep_s10v2.py not found in $(pwd)."
  exit 2
fi
if [ ! -f figures_ff_r167_s10v2/freefree_manifest.json ]; then
  echo "FATAL: figures_ff_r167_s10v2/freefree_manifest.json not found. If"
  echo "figures_ff_r167_s10v2.zip is present instead, unzip it first:"
  echo "  mkdir -p figures_ff_r167_s10v2 && unzip -o -j figures_ff_r167_s10v2.zip freefree_manifest.json -d figures_ff_r167_s10v2"
  exit 2
fi

echo "Job $SLURM_JOB_ID on $(hostname): cpus-per-task=$SLURM_CPUS_PER_TASK, N_WORKERS=$N_WORKERS, RATIO_TAG=$RATIO_TAG"
echo "Start: $(date)"

python3 -u probe_ff_residual_screen_geomsweep_s10v2.py

echo "End: $(date)"
