#!/bin/bash
#SBATCH --job-name=p5_thick
#SBATCH --output=piezo_p5_thickness_sweep_2026-09-20_%j.out
#SBATCH --error=piezo_p5_thickness_sweep_2026-09-20_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:20:00
#SBATCH --mem=4G
#SBATCH --partition=requeue
#SBATCH --requeue

# Paper 5 roadmap Sec 8 item 6 optional-cheap: thickness sweep (H/r_o),
# the homogeneous analogue of Paper 4's h1/2h curve. Pure probe-level
# scan on PiezoMonolithicOutOfPlaneSolver, F-F n=0 only -- no new
# solver code, no SOLVER_VERSION bump. See
# probe_piezo_p5_thickness_sweep_2026-09-20.py's own docstring for the
# full G_found/G_theorem/G_longwave pre-registration and why H drops
# out of the long-wave estimate (PAPER5_DERIVATION.md Sec 6), making
# this a genuinely falsifiable check, not just a curve for its own sake.
#
# In-sandbox timing this session (dps=30, inner_iters=35, scan_n=200,
# 9 H/r_o points from 0.005 to 0.18): ~101s end to end, SENTINEL
# PASS_ALL. Split stays flat (~2.35%) below H/r_o~0.025 as the long-
# wave estimate predicts, then drifts down to 2.25% by H/r_o=0.18
# (-3.98% relative) -- a real, small finite-thickness correction.
# 1 CPU is enough; 20 min is generous headroom.

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMBA_NUM_THREADS=1

echo "job=$SLURM_JOB_ID host=$(hostname) start=$(date)"

module load python/3.12.1
source /home/ghmkfh/PythonMill/Plate_Solver_Package/venv/bin/activate
cd /home/ghmkfh/PythonMill/Plate_Solver_Package/Paper5_Monolithic

export PKG_PATH=/home/ghmkfh/PythonMill/Plate_Solver_Package
export EXPECT_SOLVER_VERSION=2026-07-10.s10
export PIEZO_DPS=60
export PIEZO_INNER_ITERS=45
export PIEZO_SCAN_N=250

if [ ! -f probe_piezo_p5_thickness_sweep_2026-09-20.py ]; then
  echo "FATAL: probe_piezo_p5_thickness_sweep_2026-09-20.py not found in $(pwd)"
  exit 2
fi
if ! grep -q "def coupled_bisect" "$PKG_PATH/plate_solver/piezo_monolithic.py"; then
  echo "FATAL: $PKG_PATH/plate_solver/piezo_monolithic.py has no coupled_bisect"
  exit 2
fi

python3 -u probe_piezo_p5_thickness_sweep_2026-09-20.py
echo "end=$(date)"
