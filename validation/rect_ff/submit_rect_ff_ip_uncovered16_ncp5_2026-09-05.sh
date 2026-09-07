#!/bin/bash
#SBATCH --job-name=rect_ff_ip_uncov16
#SBATCH --output=rect_ff_ip_uncov16_%j.out
#SBATCH --error=rect_ff_ip_uncov16_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:30:00
#SBATCH --requeue
#SBATCH --partition=requeue
#
# Mop-up for the 16 FE targets job 2456002's n_cpair=3 full re-discovery
# still couldn't claim (out of the original ~130+ FE targets in
# [0.02,2.50], MATCH count went 38 -> 76 there). Two of the 16 are former
# OLD matches that stopped reproducing -- (2.0,ANTI,1.27148) and
# (3.0,SYM,1.66346) -- both immediately adjacent to a different FE root
# that WAS claimed instead, consistent with closely-spaced-root crowding,
# not a real disappearance. Scans a +/-0.05 window centered directly on
# each FE value (not seeded from a possibly-crowded discovery dip) at
# n_cpair=5 -- the strongest basis 2455677 already validated recovers
# ANTI targets n_cpair=3 alone missed, now applied to SYM too and to this
# specific residual list.
# Cost: 2455677 ran 26 targets (window +/-0.05, step 0.002, n_cpair=5) at
# 35-40s each (~0.7-0.8 s/pt over ~51 pts). 16 targets here -> ~10 min.
# COPY: p5_rect_ff_lib.py (must have IP_FE_LISTS), this probe, this
# script.

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export BLIS_NUM_THREADS=1
export DPS=30
export EXPECT_SOLVER_VERSION=2026-07-10.s10
export PKG_PATH=/home/ghmkfh/PythonMill/Plate_Solver_Package
export CHECKPOINT_DIR=./p5_ip_uncov16_checkpoints

echo "Job $SLURM_JOB_ID on $(hostname)  start: $(date)"
module load python/3.12.1
source /home/ghmkfh/PythonMill/Plate_Solver_Package/venv/bin/activate
cd /home/ghmkfh/PythonMill/Plate_Solver_Package/
python3 -u probe_rect_ff_ip_uncovered16_ncp5_2026-09-05.py
echo "Job complete: $(date)"
