#!/bin/bash
#SBATCH --job-name=rect_ff_ip_anti_wide
#SBATCH --output=rect_ff_ip_anti_wide_%j.out
#SBATCH --error=rect_ff_ip_anti_wide_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
#SBATCH --requeue
#SBATCH --partition=requeue
#
# Widen the IP FFFF ANTI persist criterion (n_cpair 3->5, window
# +/-0.02->+/-0.05) on the 26 near-miss (<=6%) FAIL_P dips job 2455570
# found in the ANTI branch, to test whether ANTI's much lower
# persist-pass rate than SYM is a persist-basis/window limitation
# rather than a genuine absence of those modes.
# Cost estimate: 26 targets x 51 pts/window x ~1.0-1.3 s/pt (n_cpair=5,
# im_cap=7, extrapolated from the ncp0->ncp3 ~1.8x scaling seen in
# 2455570) ~ 26x51x1.2 ~ 1600s (~27 min). Budgeted generously.
# COPY: p5_rect_ff_lib.py, this probe, this script.
# Independent of the OOP reopen job and the discovery-gap job.

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export BLIS_NUM_THREADS=1
export DPS=30
export EXPECT_SOLVER_VERSION=2026-07-10.s10
export PKG_PATH=/home/ghmkfh/PythonMill/Plate_Solver_Package
export CHECKPOINT_DIR=./p5_ip_anti_wide_checkpoints

echo "Job $SLURM_JOB_ID on $(hostname)  start: $(date)"
module load python/3.12.1
source /home/ghmkfh/PythonMill/Plate_Solver_Package/venv/bin/activate
cd /home/ghmkfh/PythonMill/Plate_Solver_Package/
python3 -u probe_rect_ff_ip_anti_persist_enlarged_2026-09-05.py
echo "Job complete: $(date)"
