#!/bin/bash
#SBATCH --job-name=shi14_t3_n28
#SBATCH --output=probe_oop_shi2014_t3_larger_basis_v1_%j.out
#SBATCH --error=probe_oop_shi2014_t3_larger_basis_v1_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=03:00:00
#SBATCH --mem=8G
#SBATCH --partition=requeue
#SBATCH --requeue

# Rank 10: Shi 2014 target 38.428 at n=20/28/36.
# Job 2411852 already showed n=20 bottoms at -3.297. This only
# asks whether a larger achieved basis clears -3.5.
# Submit from /home/ghmkfh/PythonMill/Plate_Solver_Package.

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMBA_NUM_THREADS=1

echo "job=$SLURM_JOB_ID host=$(hostname) start=$(date)"

module load python/3.12.1
source /home/ghmkfh/PythonMill/Plate_Solver_Package/venv/bin/activate
cd /home/ghmkfh/PythonMill/Plate_Solver_Package

export DPS=40
export EXPECT_SOLVER_VERSION=2026-07-10.s10
export PKG_PATH=/home/ghmkfh/PythonMill/Plate_Solver_Package
export N_WORKERS=8

if [ ! -f probe_oop_shi2014_t3_larger_basis_v1.py ]; then
  echo "FATAL: probe_oop_shi2014_t3_larger_basis_v1.py not found in $(pwd)"
  exit 2
fi

echo "N_WORKERS=$N_WORKERS"
python3 -u probe_oop_shi2014_t3_larger_basis_v1.py
echo "end=$(date)"
