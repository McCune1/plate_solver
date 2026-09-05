#!/bin/bash
#SBATCH --job-name=ff_oop_spurious_table
#SBATCH --output=ff_oop_spurious_table_%j.out
#SBATCH --error=ff_oop_spurious_table_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=13
#SBATCH --time=00:45:00
#SBATCH --partition=general

# Reconfirms, fresh against the currently deployed checkpoint, the
# residual-screen numbers (block/rV/rM/maxW/verdict) for all 13 raw Omega
# in the FF-P1 OOP free-free spectrum: the 8 physical modes currently in
# PAPER1_FREEFREE_DRAFT.tex's tab:oop, plus the 5 previously-catalogued
# spurious determinant zeros (LESSONS_LEARNED.md Sec 10.8 / midpoint
# report Sec 3). Requested 2026-07-19 in response to a SuperGrok review
# (ReviewGrok2.txt) asking for a supplementary table of the spurious
# zeros' raw Omega/residual data -- the stale 2026-07-09/07-13 numbers are
# NOT reused directly; this job regenerates them so the paper only carries
# freshly-verified values, per this project's own validate-by-probe
# standard. See the probe's own docstring for full pre-registered
# interpretation criteria and cost accounting.
#
# Diagnostic/data-generation only: writes a plain-text log to stdout only.
# No package changes, no SOLVER_VERSION implications.

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export BLIS_NUM_THREADS=1

export N_WORKERS=13
export DPS=40
export PKG_PATH=/home/ghmkfh/PythonMill/Plate_Solver_Package
export EXPECT_SOLVER_VERSION=2026-07-10.s10

module load python/3.12.1
source /home/ghmkfh/PythonMill/Plate_Solver_Package/venv/bin/activate

cd /home/ghmkfh/PythonMill/Plate_Solver_Package/

PROBE=""
if [ -f validation/probe_ff_oop_spurious_table.py ]; then
  PROBE=validation/probe_ff_oop_spurious_table.py
elif [ -f probe_ff_oop_spurious_table.py ]; then
  PROBE=probe_ff_oop_spurious_table.py
else
  echo "FATAL: probe_ff_oop_spurious_table.py not found in validation/ or $(pwd)."
  exit 2
fi
echo "probe: $PROBE"

echo "Job $SLURM_JOB_ID on $(hostname)"
echo "Start: $(date)"

python3 -u "$PROBE"

echo "End: $(date)"
