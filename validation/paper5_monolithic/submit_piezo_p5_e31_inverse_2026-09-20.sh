#!/bin/bash
#SBATCH --job-name=p5_e31_inv
#SBATCH --output=piezo_p5_e31_inverse_2026-09-20_%j.out
#SBATCH --error=piezo_p5_e31_inverse_2026-09-20_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:30:00
#SBATCH --mem=4G
#SBATCH --partition=requeue
#SBATCH --requeue

# Paper 5 Option A synthetic inverse: recover e31 from the SC-vs-elastic
# F-F n=0 flexural gap on the homogeneous thickness-poled PZT-4 Kirchhoff
# ring, c11_bar (C11E/C12E/C13E/C33E) held known. NOT an OC-vs-SC
# inversion -- Phase 0 theorem (PAPER5_DERIVATION.md Sec 4, LESSONS_
# LEARNED.md Sec 18.199/18.203, FE job 2520476) says fully-electroded OC
# flexure coincides with SC flexure on this stacking, for every
# mechanical BC, so there is no OC frequency solver to invert.
# See probe_piezo_p5_e31_inverse_2026-09-20.py's own docstring for the
# full G_worker/G0/G1/G2/G_reflection pre-registration and the
# newly-found e31 reflection-twin structural fact.
#
# Fully serial, same as probe_piezo_p4_sequential_e_e31_2026-09-16.py:
# the outer shooting loop needs each inner coupled_bisect call in
# sequence, so ProcessPoolExecutor buys nothing here. 1 CPU is enough.
# In-sandbox timing this session (dps=60, inner_iters=45, outer_iters=30,
# both F-F and C-C legs): a few minutes end to end; 30 min is generous
# headroom.

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
export PIEZO_OUTER_ITERS=30
export PIEZO_RUN_CC=1

if [ ! -f probe_piezo_p5_e31_inverse_2026-09-20.py ]; then
  echo "FATAL: probe_piezo_p5_e31_inverse_2026-09-20.py not found in $(pwd)"
  exit 2
fi
if ! grep -q "def coupled_bisect" "$PKG_PATH/plate_solver/piezo_monolithic.py"; then
  echo "FATAL: $PKG_PATH/plate_solver/piezo_monolithic.py has no coupled_bisect"
  echo "  Push plate_solver/piezo_monolithic.py first."
  exit 2
fi
if ! grep -q "_piezo_mono_coupled_ff_root_worker" "$PKG_PATH/plate_solver/workers.py"; then
  echo "FATAL: $PKG_PATH/plate_solver/workers.py has no"
  echo "  _piezo_mono_coupled_ff_root_worker. Push plate_solver/workers.py first."
  exit 2
fi

python3 -u probe_piezo_p5_e31_inverse_2026-09-20.py
echo "end=$(date)"
