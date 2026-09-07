#!/bin/bash
#SBATCH --job-name=rect_ff_ip_match_refine
#SBATCH --output=rect_ff_ip_match_refine_%j.out
#SBATCH --error=rect_ff_ip_match_refine_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --requeue
#SBATCH --partition=requeue
#
# Refine 2455570's 38 persist-trial MATCHES to a 3% cut: production
# basis (n_real=5,n_cpair=0) rescan +/-0.02 step 0.001 around each,
# then reconfirm persist (n_cpair=3) at the refined point.
# Cost estimate: 38 targets x 41 pts x ~0.35 s/pt (production basis)
# ~ 545s, plus 38 persist confirmations x 21 pts x ~0.65 s/pt ~ 518s.
# Total ~1060s (~18 min).
# COPY: p5_rect_ff_lib.py (must have IP_FE_LISTS), this probe, this
# script. Independent of the other three P5 IP/OOP jobs.

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export BLIS_NUM_THREADS=1
export DPS=30
export EXPECT_SOLVER_VERSION=2026-07-10.s10
export PKG_PATH=/home/ghmkfh/PythonMill/Plate_Solver_Package
export CHECKPOINT_DIR=./p5_ip_match_refine_checkpoints

echo "Job $SLURM_JOB_ID on $(hostname)  start: $(date)"
module load python/3.12.1
source /home/ghmkfh/PythonMill/Plate_Solver_Package/venv/bin/activate
cd /home/ghmkfh/PythonMill/Plate_Solver_Package/
python3 -u probe_rect_ff_ip_match_refine_2026-09-05.py
echo "Job complete: $(date)"
