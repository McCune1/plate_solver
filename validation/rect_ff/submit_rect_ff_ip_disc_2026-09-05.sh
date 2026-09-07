#!/bin/bash
#SBATCH --job-name=rect_ff_ip_disc
#SBATCH --output=rect_ff_ip_disc_%j.out
#SBATCH --error=rect_ff_ip_disc_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=06:00:00
#SBATCH --requeue
#SBATCH --partition=requeue
#
# Rectangular FFFF IP discovery [0.02, 2.50] at five l/b, persist trial
# n_cpair=3 (NOT Screen B), pair to job 2455569 FE.
# COPY: plate_solver/core_solvers.py (must have RectIPAssembler bc=free_free),
#       p5_rect_ff_lib.py (must have IP_FE_LISTS), this probe, this script.
# CPU-only. Do not run as MAPDL. Independent of OOP leftover jobs.

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export BLIS_NUM_THREADS=1
export DPS=30
export EXPECT_SOLVER_VERSION=2026-07-10.s10
export PKG_PATH=/home/ghmkfh/PythonMill/Plate_Solver_Package
export CHECKPOINT_DIR=./p5_ip_disc_checkpoints

echo "Job $SLURM_JOB_ID on $(hostname)  start: $(date)"
module load python/3.12.1
source /home/ghmkfh/PythonMill/Plate_Solver_Package/venv/bin/activate
cd /home/ghmkfh/PythonMill/Plate_Solver_Package/
python3 -u probe_rect_ff_ip_disc_2026-09-05.py
echo "Job complete: $(date)"
