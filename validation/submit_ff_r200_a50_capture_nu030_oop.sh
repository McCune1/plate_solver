#!/bin/bash
#SBATCH --job-name=ff_r200a50_oop
#SBATCH --output=ff_r200a50_nu030_oop_%j.out
#SBATCH --error=ff_r200a50_nu030_oop_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=10:00:00
#SBATCH --requeue
#SBATCH --partition=requeue

# Rank 16 OOP capture: r0/2b=2.0, 2Theta/pi=0.5, nu=0.30.
# Crossed geometry (r200 ratio × FF-P1 angle). Plain production default.
# Submit from /home/ghmkfh/PythonMill/Plate_Solver_Package.

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export BLIS_NUM_THREADS=1

export PKG_PATH=/home/ghmkfh/PythonMill/Plate_Solver_Package
export DPS=40
export EXPECT_SOLVER_VERSION=2026-07-10.s10
export PART=1
export N_WORKERS=${SLURM_CPUS_PER_TASK:-16}

module load python/3.12.1
source /home/ghmkfh/PythonMill/Plate_Solver_Package/venv/bin/activate
cd /home/ghmkfh/PythonMill/Plate_Solver_Package

if [ ! -f probe_ff_r200_a50_capture_nu030.py ]; then
  echo "FATAL: probe_ff_r200_a50_capture_nu030.py not found in $(pwd)."
  exit 2
fi

echo "Job $SLURM_JOB_ID on $(hostname): PART=$PART (OOP), N_WORKERS=$N_WORKERS"
echo "Start: $(date)"

python3 -u probe_ff_r200_a50_capture_nu030.py

echo "End: $(date)"
