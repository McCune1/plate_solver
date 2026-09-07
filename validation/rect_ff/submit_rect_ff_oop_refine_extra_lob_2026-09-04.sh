#!/bin/bash
#SBATCH --job-name=rect_ff_ref_elob
#SBATCH --output=rect_ff_ref_elob_%j.out
#SBATCH --error=rect_ff_ref_elob_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=03:00:00
#SBATCH --requeue
#SBATCH --partition=requeue
#
# Extra-lob refine (2453876 candidates vs 2454654 FE). CPU-only.
# Independent of miss-windows and highlam-extra-lob -- may run in parallel.
# COPY: p5_rect_ff_lib.py (MUST be the copy with l/b=2.0/3.0 FE_LISTS),
#       probe_rect_ff_oop_refine_extra_lob_2026-09-04.py, this script

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export BLIS_NUM_THREADS=1
export DPS=30
export EXPECT_SOLVER_VERSION=2026-07-10.s10
export PKG_PATH=/home/ghmkfh/PythonMill/Plate_Solver_Package
export CHECKPOINT_DIR=./p5_refine_extra_lob_checkpoints

echo "Job $SLURM_JOB_ID on $(hostname)  start: $(date)"
module load python/3.12.1
source /home/ghmkfh/PythonMill/Plate_Solver_Package/venv/bin/activate
cd /home/ghmkfh/PythonMill/Plate_Solver_Package/
python3 -u probe_rect_ff_oop_refine_extra_lob_2026-09-04.py
echo "Job complete: $(date)"
