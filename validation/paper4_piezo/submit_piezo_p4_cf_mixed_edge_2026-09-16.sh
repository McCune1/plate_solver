#!/bin/bash
#SBATCH --job-name=p4_piezo_cf
#SBATCH --output=piezo_p4_cf_mixed_edge_2026-09-16_%j.out
#SBATCH --error=piezo_p4_cf_mixed_edge_2026-09-16_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:30:00
#SBATCH --mem=4G
#SBATCH --partition=requeue
#SBATCH --requeue

# Paper 4 lever 10: C-F / F-C mixed mechanical BC (elastic, SC 6x6, OC).
# See probe_piezo_p4_cf_mixed_edge_2026-09-16.py and
# LESSONS_LEARNED.md Sec 18.192. 4x4 + a few 6x6; 30 min is generous.

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMBA_NUM_THREADS=1

echo "job=$SLURM_JOB_ID host=$(hostname) start=$(date)"

module load python/3.12.1
source /home/ghmkfh/PythonMill/Plate_Solver_Package/venv/bin/activate
cd /home/ghmkfh/PythonMill/Plate_Solver_Package/Paper4_Piezo

export PKG_PATH=/home/ghmkfh/PythonMill/Plate_Solver_Package
export EXPECT_SOLVER_VERSION=2026-07-10.s10
export PIEZO_DPS=60
export PIEZO_ITERS=40
export PIEZO_OUTER_ITERS=26

if [ ! -f probe_piezo_p4_cf_mixed_edge_2026-09-16.py ]; then
  echo "FATAL: probe_piezo_p4_cf_mixed_edge_2026-09-16.py not found in $(pwd)"
  exit 2
fi
if ! grep -q "def oc_cf_det" "$PKG_PATH/plate_solver/piezo_solver.py"; then
  echo "FATAL: $PKG_PATH/plate_solver/piezo_solver.py has no oc_cf_det"
  echo "  scp/push plate_solver/piezo_solver.py and workers.py first."
  exit 2
fi

python3 -u probe_piezo_p4_cf_mixed_edge_2026-09-16.py
echo "end=$(date)"
