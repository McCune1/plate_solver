#!/bin/bash
#SBATCH --job-name=ff_oop_artifact_2p8176
#SBATCH --output=ff_oop_artifact_2p8176_%j.out
#SBATCH --error=ff_oop_artifact_2p8176_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3
#SBATCH --time=00:30:00
#SBATCH --partition=requeue
#SBATCH --requeue

# Computes the pointwise edge-residual screen (block/rV/rM/maxW/verdict)
# for the FF-P1 OOP spurious zero at raw Omega=2.817613 -- the EVEN zero
# shown in the paper's Fig. 1 and Fig. 2(d), which sits above the
# validated flexural window and so is not in Table 5.
#
# WHY NOW: the 2026-08-11 review pass found that Fig. 2's caption told the
# reader to "see Sec. 5 for the residual values" of this zero, but no
# rV/rM had ever been computed for it anywhere in the project. The caption
# was rewritten to stop promising them; this job produces them so the
# caption can state them instead.
#
# Runs the unknown alongside two published controls (0.383695 REAL-like,
# 2.121685 ARTIFACT-like) in the SAME run, so a surprising reading can be
# separated from screen drift between s8 (job 2325529, which produced
# Table 5) and the currently deployed s10. See the probe docstring for the
# full pre-registered interpretation -- controls are read FIRST.
#
# Partition note: general/GPU became paid 2026-08-04, so this uses the
# preemptible requeue partition. The job is ~2 min of real work; if it is
# preempted, just resubmit -- there is no checkpoint and none is needed.
#
# Diagnostic/data-generation only: stdout log, no package changes, no
# SOLVER_VERSION implications.

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export BLIS_NUM_THREADS=1

export N_WORKERS=3
export DPS=40
export PKG_PATH=/home/ghmkfh/PythonMill/Plate_Solver_Package
export EXPECT_SOLVER_VERSION=2026-07-10.s10

module load python/3.12.1
source /home/ghmkfh/PythonMill/Plate_Solver_Package/venv/bin/activate

cd /home/ghmkfh/PythonMill/Plate_Solver_Package/

if [ ! -f probe_ff_oop_artifact_2p8176_v1.py ]; then
  echo "FATAL: probe_ff_oop_artifact_2p8176_v1.py not found in $(pwd)."
  exit 2
fi

echo "Job $SLURM_JOB_ID on $(hostname)"
echo "Start: $(date)"

python3 -u probe_ff_oop_artifact_2p8176_v1.py

echo "End: $(date)"
