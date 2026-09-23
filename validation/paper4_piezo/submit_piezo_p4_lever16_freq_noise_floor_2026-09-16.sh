#!/bin/bash
#SBATCH --job-name=p4_piezo_l16
#SBATCH --output=piezo_p4_lever16_freq_noise_floor_2026-09-16_%j.out
#SBATCH --error=piezo_p4_lever16_freq_noise_floor_2026-09-16_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:30:00
#SBATCH --mem=4G
#SBATCH --partition=requeue
#SBATCH --requeue

# Paper 4 lever 16: non-0.01% noise curve (F-F and C-F n=0 OC e31
# uncertainty vs assumed relative frequency error, plus the
# Boeringa-McCune impedance/SLDV lab floor). Same-operator sweep of
# oc_ff_bisect / oc_cf_bisect -- no new derivation, no package change.
# See probe_piezo_p4_lever16_freq_noise_floor_2026-09-16.py and
# LESSONS_LEARNED.md (lever 16 entry, appended alongside Sec 18.192).
# In-sandbox SENTINEL PASS_ALL, ~85s wall at dps=60 (2026-09-16).

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

if [ ! -f probe_piezo_p4_lever16_freq_noise_floor_2026-09-16.py ]; then
  echo "FATAL: probe_piezo_p4_lever16_freq_noise_floor_2026-09-16.py not found in $(pwd)"
  exit 2
fi
if ! grep -q "def oc_cf_bisect" "$PKG_PATH/plate_solver/piezo_solver.py"; then
  echo "FATAL: $PKG_PATH/plate_solver/piezo_solver.py has no oc_cf_bisect"
  echo "  scp/push plate_solver/piezo_solver.py and workers.py first."
  exit 2
fi

python3 -u probe_piezo_p4_lever16_freq_noise_floor_2026-09-16.py
echo "end=$(date)"
