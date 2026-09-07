#!/bin/bash
#SBATCH --job-name=persist_thresh_sens
#SBATCH --output=rect_ff_persist_threshold_sensitivity_%j.out
#SBATCH --error=rect_ff_persist_threshold_sensitivity_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#
# Sensitivity/margin check for Paper 2's frozen Screen-B persistence cut
# (p5_rect_ff_lib.PERSIST_DL = 0.01), the rectangular companion's analogue
# of Paper 1's frozen SUBDOM 0.70 threshold and geometry-sweep tolerance
# insensitivity checks. Paper 2 does not yet have one of its own for this
# threshold -- this closes that specific rigor-gap item.
#
# Calls the exact, unmodified, already-frozen production functions
# p5_rect_ff_lib.persist_one / .persist_one_ip (imported, never copied or
# edited) once for EVERY published production match in Paper 2 -- all 26
# rows of Table tab:oop-primary, all 30 rows of tab:oop-high, and all 92
# rows of tab:ip92 (148 candidates total) -- then reclassifies each one's
# already-computed dip distance against a dozen alternate thresholds in
# plain Python arithmetic (zero extra plate_solver calls past the initial
# 148). No monkeypatch, no retuning of PERSIST_DL, no paper edits.
#
# Reuses the same production im_cap/basis settings baked into
# persist_one/persist_one_ip themselves -- no new branch-resolution risk.
# Each candidate self-checkpoints to its own JSON file, so a resubmitted
# job resumes for free if it runs long.
#
# Needs NO ANSYS/FE data at all (purely internal to plate_solver) --
# unlike the MAC/symsweep/symvsanti jobs, this one does not check for the
# .QUEUE_OK sentinels.
#
# Copy/paste on the cluster:
#   squeue -u $USER
#   cd /home/ghmkfh/PythonMill/Plate_Solver_Package
#   sbatch submit_rect_ff_persist_threshold_sensitivity_2026-09-07.sh
#
# Fidelity-gate-only smoke test (seconds, 6 candidates only):
#   SMOKE=1 sbatch submit_rect_ff_persist_threshold_sensitivity_2026-09-07.sh
#
# Internal soft wall-clock budget the probe checkpoints against (minutes;
# default 100, leaving headroom under this script's 2h SBATCH limit --
# expect roughly 148 candidates x 21-point window; per earlier jobs'
# sigma-only cost this should land well under an hour, but the budget
# and per-candidate checkpointing make a slow node harmless):
#   MAX_MINUTES=80 sbatch submit_rect_ff_persist_threshold_sensitivity_2026-09-07.sh
#
set -euo pipefail
cd /home/ghmkfh/PythonMill/Plate_Solver_Package

module load python/3.12.1
source venv/bin/activate
export PKG_PATH="."
export DPS="${DPS:-30}"
export EXPECT_SOLVER_VERSION=2026-07-10.s10
export PYTHONUNBUFFERED=1
export SMOKE="${SMOKE:-0}"
export MAX_MINUTES="${MAX_MINUTES:-100}"
python3 -u probe_rect_ff_persist_threshold_sensitivity_2026-09-07.py
