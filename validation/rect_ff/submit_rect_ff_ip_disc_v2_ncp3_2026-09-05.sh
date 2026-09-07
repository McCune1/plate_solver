#!/bin/bash
#SBATCH --job-name=rect_ff_ip_discv2
#SBATCH --output=rect_ff_ip_discv2_%j.out
#SBATCH --error=rect_ff_ip_discv2_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=03:00:00
#SBATCH --requeue
#SBATCH --partition=requeue
#
# Production-candidate re-discovery: same full [0.02,2.50] sweep as
# 2455570, same step convention (0.01 SYM / 0.02 ANTI), but discovery
# basis is now n_cpair=3 (persist basis) directly instead of n_cpair=0.
# Motivated by 2455676 (14 basis-starved dips in [1.20,1.50], invisible
# at n_cpair=0) and 2455677 (20/26 ANTI FAIL_P recovered at n_cpair=5,
# most already recoverable at n_cpair=3 within a wider window). Adopts
# the 3% MATCH cut 2455678 validated (37/38 held), superseding 2455570's
# 5% trial cut. Every dip found is fine-refined (+/-0.02, step 0.002,
# same basis) via persist_one_ip before being scored, same as 2455678.
# Cross-references every MATCH against 2455570's original 38 (NEW vs
# OLD) so the paper's match list can be updated by diff, not rewritten.
# Cost: 2455570 (n_cpair=0 discovery + persist-on-hits only) ran 2164.8s.
# This job's discovery step alone is ~1.9x costlier per point (n_cpair=3
# ~0.74s/pt vs n_cpair=0 ~0.39s/pt per 2455676's own two passes) over the
# same ~1870 discovery points (~1380s), plus refine on however many dips
# turn up (likely more than 47, since n_cpair=0 undercounted) at ~15.5s
# each. Budget 3h; expect well under half that.
# COPY: p5_rect_ff_lib.py (must have IP_FE_LISTS), this probe, this
# script. Independent of the three 2026-09-05 follow-up jobs.

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export BLIS_NUM_THREADS=1
export DPS=30
export EXPECT_SOLVER_VERSION=2026-07-10.s10
export PKG_PATH=/home/ghmkfh/PythonMill/Plate_Solver_Package
export CHECKPOINT_DIR=./p5_ip_disc_v2_checkpoints

echo "Job $SLURM_JOB_ID on $(hostname)  start: $(date)"
module load python/3.12.1
source /home/ghmkfh/PythonMill/Plate_Solver_Package/venv/bin/activate
cd /home/ghmkfh/PythonMill/Plate_Solver_Package/
python3 -u probe_rect_ff_ip_disc_v2_ncp3_2026-09-05.py
echo "Job complete: $(date)"
