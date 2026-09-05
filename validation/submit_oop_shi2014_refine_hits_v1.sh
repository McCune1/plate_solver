#!/bin/bash
#SBATCH --job-name=shi14_refine
#SBATCH --output=probe_oop_shi2014_refine_hits_v1_%j.out
#SBATCH --error=probe_oop_shi2014_refine_hits_v1_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH --partition=requeue
#SBATCH --requeue

# Refines the five confirmed Shi2014 Table-3 HITs (job 2410389) from
# "a singularity exists inside the predicted window" to an independently
# located root with a quotable err% against the published value.
#
# WHY: job 2410389 reported rel% = 0.000 for all five, which looks perfect
# but only means its scan grid contained the predicted centre and arg-max
# picked it -- no refinement ever happened. That is why the paper's Sec. 6.9
# can currently claim only a coarse window confirmation, and why the primary
# flexural tables still have to be described as FE-validated only. This run
# is the single largest evidentiary upgrade available to the paper.
#
# Target 3 (38.428) is excluded: job 2411852's 33-point dense scan of the
# tight +-5% window confirmed its dip bottoms at -3.297 and never reaches
# -3.5, so it is a genuine near-miss, not an under-resolution artifact.
#
# Includes a fidelity gate (each predicted centre must reproduce job
# 2410389's depth to within 0.05 decades) and a bracket-shift stability
# replica per target -- a refined root that moves when the bracket is
# nudged is the golden-section lock-on failure this project has hit before,
# and is reported as PARTIAL rather than quoted.
#
# Partition: general/GPU became paid 2026-08-04, so this uses requeue.
# ~35 min of real work; if preempted, just resubmit (no checkpoint needed).
#
# Diagnostic only: stdout log, no package changes, no SOLVER_VERSION impact.

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export BLIS_NUM_THREADS=1
export NUMBA_NUM_THREADS=1

export N_WORKERS=5
export DPS=40
export PKG_PATH=/home/ghmkfh/PythonMill/Plate_Solver_Package
export EXPECT_SOLVER_VERSION=2026-07-10.s10

module load python/3.12.1
source /home/ghmkfh/PythonMill/Plate_Solver_Package/venv/bin/activate

cd /home/ghmkfh/PythonMill/Plate_Solver_Package/

if [ ! -f probe_oop_shi2014_refine_hits_v1.py ]; then
  echo "FATAL: probe_oop_shi2014_refine_hits_v1.py not found in $(pwd)."
  exit 2
fi

echo "job=$SLURM_JOB_ID host=$(hostname) start=$(date)"

python3 -u probe_oop_shi2014_refine_hits_v1.py

echo "end=$(date)"
