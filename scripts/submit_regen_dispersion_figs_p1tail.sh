#!/bin/bash
#SBATCH --job-name=regen_dispersion_figs_p1tail
#SBATCH --output=regen_dispersion_figs_p1tail_%j.out
#SBATCH --error=regen_dispersion_figs_p1tail_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=00:30:00
#SBATCH --partition=general

# Same as submit_regen_dispersion_figs.sh (see that file's header for the
# full rationale) but for the p1tail checkpoint/figdir pair: whichever of
# 2*Theta=1.0/1.25/1.5 pi (r0/(2b)=1.25), 1.0 pi (r0/(2b)=5/3), 1.0 pi
# (r0/(2b)=5/2) are checkpointed in research23_checkpoint_p1tail.json at
# the time this runs, figures_p1tail/.
#
# TIMING NOTE: only run this once submit_cantilever_table_regen_p1tail_
# resume.sh (job 2324437 or its successor) has fully finished, or at least
# once geometry 5 (r0/(2b)=2.5) has been checkpointed -- this script only
# reads the checkpoint once at start, so if it runs while that job is
# still mid-geometry-5 it will simply regenerate figures for whichever
# geometries are checkpointed so far (harmless, just incomplete) and need
# a cheap re-run afterward to pick up the last one. Not destructive either
# way, just wasted effort if run too early.
#
# PUSH BEFORE RUNNING: plate_solver/plotting.py (the 2026-07-16 fix) and
# this script + probe_regen_dispersion_figs.py.
# PULL AFTER RUNNING: the .out log only, to confirm every geometry OK with
# no FAILED lines -- PNGs are overwritten in place in figures_p1tail/.

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export BLIS_NUM_THREADS=1

export PKG_PATH=/home/ghmkfh/PythonMill/Plate_Solver_Package
export EXPECT_SOLVER_VERSION=2026-07-10.s8
export CHECKPOINT=research23_checkpoint_p1tail.json
export FIGDIR=figures_p1tail
export N_WORKERS=8

module load python/3.12.1
source /home/ghmkfh/PythonMill/Plate_Solver_Package/venv/bin/activate

cd /home/ghmkfh/PythonMill/Plate_Solver_Package/

if [ ! -f probe_regen_dispersion_figs.py ]; then
  echo "FATAL: probe_regen_dispersion_figs.py not found in $(pwd)."
  exit 2
fi
if [ ! -f "$CHECKPOINT" ]; then
  echo "FATAL: $CHECKPOINT not found in $(pwd)."
  exit 3
fi

echo "Job $SLURM_JOB_ID on $(hostname)"
echo "Start: $(date)"

python3 -u probe_regen_dispersion_figs.py

echo "End: $(date)"
