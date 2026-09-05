#!/bin/bash
#SBATCH --job-name=paper_figures_oop
#SBATCH --output=paper_figures_oop_%j.out
#SBATCH --error=paper_figures_oop_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=01:30:00
#SBATCH --partition=general

# Generates the two figures PAPER1_FREEFREE_DRAFT.tex Sec.5/6 currently
# carry as \todo placeholders (added this session per external feedback on
# the draft): (A) 4 OOP mode-shape PNGs (3 physical + 1 confirmed artifact,
# all at already-validated raw Omega -- no new mode search), and (B) a
# sigma_min(Omega) trace over [2.55,2.90] spanning the last validated real
# mode, a confirmed artifact, and the still-open OOP ceiling ambiguity
# (Sec 20.4 of LESSONS_LEARNED), so the figure directly illustrates the
# paper's own "one bounded ambiguity" discussion.
#
# SCOPE LIMIT: OOP (Part 1) only -- plate_solver.plotting.plot_mode_shape /
# mode_shape_grid is not implemented for the in-plane solver yet (checked
# 2026-07-14 against plotting.py's own docstring). An in-plane figure job
# is a separate, not-yet-buildable future task.
#
# 2026-07-17 RERUN NOTE: this is an unmodified re-run of a script that
# already exists/ran before (job 2323612, 2324040) -- being resubmitted to
# pick up plotting.py's MODE_SHAPE_Z_RATIO fix (0.55 -> 0.18, see
# CLAUDE.md item 6): Part A's plot_mode_shape calls will now render with
# the corrected Z-scaling automatically, since it's the same plotting.py
# import, not a code change to this script.
#
# Diagnostic/data-generation only: writes PNGs + a plain-text sigma_min
# log under paper_figures/, nothing else. No package changes, no
# SOLVER_VERSION implications. See the probe's own docstring for the full
# target list, cost accounting, and pre-registered scope.

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export BLIS_NUM_THREADS=1

export N_WORKERS=8
export DPS=40
export PKG_PATH=/home/ghmkfh/PythonMill/Plate_Solver_Package
export EXPECT_SOLVER_VERSION=2026-07-10.s10
# Mode-shape PDFs only. The sigma_min trace is already shipped as a vector
# PDF reconstructed from the 37-point .txt; do not spend ~1 hour redoing it.
export SKIP_TRACE=1
export OUTDIR=paper_figures

module load python/3.12.1
source /home/ghmkfh/PythonMill/Plate_Solver_Package/venv/bin/activate

cd /home/ghmkfh/PythonMill/Plate_Solver_Package/

if [ ! -f probe_paper_figures_oop.py ]; then
  echo "FATAL: probe_paper_figures_oop.py not found in $(pwd)."
  exit 2
fi

echo "Job $SLURM_JOB_ID on $(hostname)"
echo "Start: $(date)"

python3 -u probe_paper_figures_oop.py

echo "End: $(date)"
