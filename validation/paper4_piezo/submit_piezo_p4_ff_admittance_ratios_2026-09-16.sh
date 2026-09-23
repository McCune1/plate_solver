#!/bin/bash
#SBATCH --job-name=p4_piezo_y_ratios
#SBATCH --output=piezo_p4_ff_admittance_ratios_2026-09-16_%j.out
#SBATCH --error=piezo_p4_ff_admittance_ratios_2026-09-16_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:30:00
#SBATCH --mem=4G
#SBATCH --partition=requeue
#SBATCH --requeue

# Lever 5: Y(omega) Q/V sweep at h1/2h=1/8 and 1/5. Clone of job
# 2491920's python side. Push plate_solver/piezo_solver.py (must
# contain driven_ff_qv) AND this probe before sbatch.
# From Paper4_Piezo/: sbatch submit_piezo_p4_ff_admittance_ratios_2026-09-16.sh

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

echo "job=$SLURM_JOB_ID host=$(hostname) start=$(date)"
module load python/3.12.1
source /home/ghmkfh/PythonMill/Plate_Solver_Package/venv/bin/activate
cd /home/ghmkfh/PythonMill/Plate_Solver_Package/Paper4_Piezo

export PKG_PATH=/home/ghmkfh/PythonMill/Plate_Solver_Package
export EXPECT_SOLVER_VERSION=2026-07-10.s10
export PIEZO_DPS=60

if [ ! -f probe_piezo_p4_ff_admittance_ratios_2026-09-16.py ]; then
  echo "FATAL: probe not found in $(pwd)"
  exit 2
fi
if ! grep -q "def driven_ff_qv" "$PKG_PATH/plate_solver/piezo_solver.py"; then
  echo "FATAL: deployed piezo_solver.py has no driven_ff_qv -- push it first"
  exit 2
fi

python3 -u probe_piezo_p4_ff_admittance_ratios_2026-09-16.py
echo "end=$(date)"
