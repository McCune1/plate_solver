#!/bin/bash
#SBATCH --job-name=regen_modeshape_figs_p1tail
#SBATCH --output=regen_modeshape_figs_p1tail_%j.out
#SBATCH --error=regen_modeshape_figs_p1tail_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=01:00:00
#SBATCH --partition=general

# 2026-07-19: bumped 8->32 cpus-per-task (cluster reported open). The
# --time=00:15:00 first tried here was WRONG, same mistake as the main
# script: computed from a stale "21 modes" claim. Job 2324836 (32 workers)
# was CANCELLED AT THE TIME LIMIT at 915s with only 76/147 modes done
# (~52%) -- direct JSON read (2026-07-19) confirms this checkpoint actually
# has 147 modes across 5 geometries: r0/(2b)=1.25 at 2*Theta=1.0pi (21),
# 1.25pi (27), 1.5pi (33); r0/(2b)=1.667 at 1.0pi (26); r0/(2b)=2.5 at
# 1.0pi (40). (The 1.0pi/1.25pi-at-r0=1.25 entries DUPLICATE what
# submit_regen_modeshape_figs.sh's checkpoint already covers -- see that
# script's header; not fixed here, flagged for the user.)
# probe_regen_modeshape_figs.py has no resume logic, so a resubmission
# redoes all 147 modes from scratch. Extrapolating the observed 76
# modes/915s throughput to 147 modes gives ~1770s (~29.5 min); --time=
# 01:00:00 gives ~2x margin, not a tight re-estimate, per-mode cost may
# not be uniform across these 5 geometries (different mode counts,
# possibly different per-mode cost) and this project has a history of
# underestimating these costs (see submit_escalation_root_swap_audit.sh's
# own comment re: job 2324755's 15-20x miss).
#
# Same as submit_regen_modeshape_figs.sh (see that file's header for the
# full B18 rationale) but for the p1tail checkpoint/figdir pair,
# figures_p1tail/. Run alongside (or after) submit_regen_modeshape_figs.sh
# -- independent jobs, no ordering dependency between them.
#
# PUSH BEFORE RUNNING: core_solvers.py (B18 fix), this script, and
# probe_regen_modeshape_figs.py.
# PULL AFTER RUNNING: the .out log only, to confirm 21/21 modes
# regenerated with no FAILED lines -- PNGs are overwritten in place in
# figures_p1tail/.

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
export N_WORKERS=32

module load python/3.12.1
source /home/ghmkfh/PythonMill/Plate_Solver_Package/venv/bin/activate

cd /home/ghmkfh/PythonMill/Plate_Solver_Package/

if [ ! -f probe_regen_modeshape_figs.py ]; then
  echo "FATAL: probe_regen_modeshape_figs.py not found in $(pwd)."
  exit 2
fi
if [ ! -f "$CHECKPOINT" ]; then
  echo "FATAL: $CHECKPOINT not found in $(pwd)."
  exit 3
fi

echo "Job $SLURM_JOB_ID on $(hostname)"
echo "Start: $(date)"

python3 -u probe_regen_modeshape_figs.py

echo "End: $(date)"
