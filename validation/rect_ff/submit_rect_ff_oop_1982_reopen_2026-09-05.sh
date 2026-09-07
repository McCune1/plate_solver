#!/bin/bash
#SBATCH --job-name=rect_ff_1982_reopen
#SBATCH --output=rect_ff_1982_reopen_%j.out
#SBATCH --error=rect_ff_1982_reopen_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --requeue
#SBATCH --partition=requeue
#
# Fine-step (0.0002) reopen of l/b=1.0 SYM FE=1.98249 at the PERSIST
# basis (6,3), same basis job 2455568 used when it printed
# sig@FE=1.73e-06 but "dips: NONE" at step 0.002. ~201 pts at persist
# basis; expect well under 10 min (2455568's 61-pt persist window at
# this same basis took 60.6s, i.e. ~1.0s/pt).
# COPY: p5_rect_ff_lib.py, this probe, this script.
# Independent of all IP jobs. Does not touch the corner term or Screen B.

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export BLIS_NUM_THREADS=1
export DPS=30
export EXPECT_SOLVER_VERSION=2026-07-10.s10
export PKG_PATH=/home/ghmkfh/PythonMill/Plate_Solver_Package
export CHECKPOINT_DIR=./p5_1982_reopen_checkpoints

echo "Job $SLURM_JOB_ID on $(hostname)  start: $(date)"
module load python/3.12.1
source /home/ghmkfh/PythonMill/Plate_Solver_Package/venv/bin/activate
cd /home/ghmkfh/PythonMill/Plate_Solver_Package/
python3 -u probe_rect_ff_oop_1982_reopen_2026-09-05.py
echo "Job complete: $(date)"
