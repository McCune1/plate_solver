#!/bin/bash
#SBATCH --job-name=ip06_n36_thr
#SBATCH --output=probe_ip_ffp1_ip06_n36_branch_threshold_v1_%j.out
#SBATCH --error=probe_ip_ffp1_ip06_n36_branch_threshold_v1_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:45:00
#SBATCH --mem=4G

echo "job=$SLURM_JOB_ID host=$(hostname) start=$(date)"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMBA_NUM_THREADS=1

module load python/3.12.1
source /home/ghmkfh/PythonMill/Plate_Solver_Package/venv/bin/activate

cd /home/ghmkfh/PythonMill/Plate_Solver_Package

export DPS=40
export EXPECT_SOLVER_VERSION=2026-07-10.s10
python3 -u probe_ip_ffp1_ip06_n36_branch_threshold_v1.py

echo "end=$(date)"