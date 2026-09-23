#!/bin/bash
#SBATCH --job-name=p4_piezo_n0ot
#SBATCH --output=piezo_p4_ff_n0_overtones_2026-09-16_%j.out
#SBATCH --error=piezo_p4_ff_n0_overtones_2026-09-16_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:20:00
#SBATCH --mem=4G
#SBATCH --partition=requeue
#SBATCH --requeue

# Lever 1: F-F n=0 OC radial overtones. From Paper4_Piezo/:
#   sbatch submit_piezo_p4_ff_n0_overtones_2026-09-16.sh
# Push plate_solver/piezo_solver.py (must have oc_ff_bisect).

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
echo "job=$SLURM_JOB_ID host=$(hostname) start=$(date)"
module load python/3.12.1
source /home/ghmkfh/PythonMill/Plate_Solver_Package/venv/bin/activate
cd /home/ghmkfh/PythonMill/Plate_Solver_Package/Paper4_Piezo
export PKG_PATH=/home/ghmkfh/PythonMill/Plate_Solver_Package
export EXPECT_SOLVER_VERSION=2026-07-10.s10
export PIEZO_DPS=40
if [ ! -f probe_piezo_p4_ff_n0_overtones_2026-09-16.py ]; then
  echo "FATAL: probe not found in $(pwd)"; exit 2
fi
if ! grep -q "def oc_ff_bisect" "$PKG_PATH/plate_solver/piezo_solver.py"; then
  echo "FATAL: no oc_ff_bisect -- push piezo_solver.py"; exit 2
fi
python3 -u probe_piezo_p4_ff_n0_overtones_2026-09-16.py
echo "end=$(date)"
