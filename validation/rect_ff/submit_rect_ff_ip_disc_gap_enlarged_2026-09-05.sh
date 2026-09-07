#!/bin/bash
#SBATCH --job-name=rect_ff_ip_gap
#SBATCH --output=rect_ff_ip_gap_%j.out
#SBATCH --error=rect_ff_ip_gap_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --requeue
#SBATCH --partition=requeue
#
# Enlarged-basis (n_cpair 0 vs 3) rescan of the aspect-ratio-independent
# Omega in [1.20,1.50] band where job 2455570 found NO dip at all (not
# just FAIL_P) at nearly every l/b/parity, despite an FE target there in
# almost every case. Tests the 2026-09-05 CLAUDE.md hypothesis directly
# ("likely basis n_real=5 n_cpair=0, not a missing corner").
# Cost estimate from 2455570's own timings: n_cpair=0 ~0.35-0.38 s/pt,
# n_cpair=3 (persist) ~0.65 s/pt; 31 pts x 10 (lob,sym) x 2 bases
# ~ 310x(0.35+0.65) ~ 310s, plus persist confirmations on any hits.
# Expect well under 20 minutes.
# COPY: p5_rect_ff_lib.py (must have IP_FE_LISTS), this probe, this
# script. Independent of the OOP reopen job and the ANTI-persist job.

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export BLIS_NUM_THREADS=1
export DPS=30
export EXPECT_SOLVER_VERSION=2026-07-10.s10
export PKG_PATH=/home/ghmkfh/PythonMill/Plate_Solver_Package
export CHECKPOINT_DIR=./p5_ip_gap_checkpoints

echo "Job $SLURM_JOB_ID on $(hostname)  start: $(date)"
module load python/3.12.1
source /home/ghmkfh/PythonMill/Plate_Solver_Package/venv/bin/activate
cd /home/ghmkfh/PythonMill/Plate_Solver_Package/
python3 -u probe_rect_ff_ip_disc_gap_enlarged_2026-09-05.py
echo "Job complete: $(date)"
