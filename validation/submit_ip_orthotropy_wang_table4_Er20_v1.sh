#!/bin/bash
#SBATCH --job-name=ip_wang_Er20r
#SBATCH --output=probe_ip_orthotropy_wang_table4_Er20_v1_%j.out
#SBATCH --error=probe_ip_orthotropy_wang_table4_Er20_v1_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=12:00:00
#SBATCH --mem=32G

echo "job=$SLURM_JOB_ID host=$(hostname) start=$(date)"
echo "SLURM_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK"
echo "RESUME: targets 5-8 only (1-4 seeded from job 2408466)"

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
export N_WORKERS=$SLURM_CPUS_PER_TASK

python3 -u probe_ip_orthotropy_wang_table4_Er20_v1.py

echo "end=$(date)"