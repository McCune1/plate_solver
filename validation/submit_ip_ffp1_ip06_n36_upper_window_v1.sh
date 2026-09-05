#!/bin/bash
#SBATCH --job-name=ip06_upper
#SBATCH --output=probe_ip_ffp1_ip06_n36_upper_window_v1_%j.out
#SBATCH --error=probe_ip_ffp1_ip06_n36_upper_window_v1_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --mem=4G
#SBATCH --partition=requeue
#SBATCH --requeue

# Finishes the IP-06 n=36 basis-count investigation that job 2411851
# started in the wrong window.
#
# 2411851 established (keep): at Omega=1.19704984, n_dofs=36, five repeats
# give cnt=30 identically, and a scan over [1.196, 1.199] is constant at
# 30. No per-call randomness, no threshold effect THERE.
#
# The gap: IP-06's tabulated n=36 frequency is 1.203379, OUTSIDE that
# window. The region where the root actually sits was never tested, so the
# paper's "not reproducible run to run (30 vs 34); read as unresolved"
# wording in Appendix C.1 cannot yet be revised.
#
# This run covers [1.199, 1.205] -- butting directly onto 2411851's upper
# edge so the two together span [1.196, 1.205] with no gap -- plus the
# three tabulated IP-06 frequencies directly and a repeatability set at
# 1.203379 itself.
#
# Part 0 is a fidelity gate: it re-evaluates 2411851's own fixed point and
# exits 4 unless it returns cnt=30, n_raw=9. Without that the two runs are
# not comparable and nothing else in the log means anything.
#
# All four pre-registered outcomes (A deterministic / B threshold /
# C nondeterminism / D stable-34) are usable; only A licenses changing the
# paper. See the probe docstring.
#
# Single CPU, ~40 evaluations at ~34 s = ~25 min. Partition requeue
# (general/GPU became paid 2026-08-04); if preempted just resubmit, there
# is no checkpoint and none is needed.
#
# Diagnostic only: stdout log, no package changes, no SOLVER_VERSION impact.

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export BLIS_NUM_THREADS=1
export NUMBA_NUM_THREADS=1

export DPS=40
export PKG_PATH=/home/ghmkfh/PythonMill/Plate_Solver_Package
export EXPECT_SOLVER_VERSION=2026-07-10.s10

module load python/3.12.1
source /home/ghmkfh/PythonMill/Plate_Solver_Package/venv/bin/activate

cd /home/ghmkfh/PythonMill/Plate_Solver_Package/

if [ ! -f probe_ip_ffp1_ip06_n36_upper_window_v1.py ]; then
  echo "FATAL: probe_ip_ffp1_ip06_n36_upper_window_v1.py not found in $(pwd)."
  exit 2
fi

echo "job=$SLURM_JOB_ID host=$(hostname) start=$(date)"

python3 -u probe_ip_ffp1_ip06_n36_upper_window_v1.py

echo "end=$(date)"
