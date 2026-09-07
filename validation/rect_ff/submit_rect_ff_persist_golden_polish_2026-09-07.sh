#!/bin/bash
#SBATCH --job-name=persist_golden
#SBATCH --output=rect_ff_persist_golden_polish_%j.out
#SBATCH --error=rect_ff_persist_golden_polish_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#
# Decisive follow-up to the fine-grid recheck (job 2458006), which found
# that 2 of the 3 delta=0.01000 borderline candidates from job 2457815
# got WORSE (0.0105, 0.0110 -- both > frozen PERSIST_DL=0.01) at 4x finer
# grid resolution, not better -- a uniform grid, however fine, cannot by
# itself settle whether the true continuous minimum is inside or outside
# the frozen cut. This job runs golden-section search (a continuous
# minimizer, not another grid) bracketed tightly around each candidate's
# own job-2458006 fine-grid dip location, to pin down the true minimum
# precisely.
#
# Reuses golden-section code already cluster-run elsewhere in this
# project (probe_rect_ff_oop_mac_anti_remain_2026-09-06.py's golden()),
# copied verbatim for the OOP family; the IP family gets an identical
# copy with only sigma_at swapped for sigma_at_ip -- the project's
# exact-copy / minimal-marked-change discipline, not a rewrite.
#
# Does NOT retune PERSIST_DL, does NOT edit plate_solver internals, does
# NOT touch a tabulated frequency, no SOLVER_VERSION bump, no paper edit.
# This is read-only diagnostic input to a human decision, not an
# automatic edit -- if any candidate resolves FAILS_SCREEN_B_AT_CONTINUUM,
# that does not retroactively invalidate the published frequency under
# Paper 2's stated (fixed 0.002-step) Screen B procedure, but it is a
# genuine discretization-sensitivity caveat to disclose.
#
# Copy/paste on the cluster:
#   cd /home/ghmkfh/PythonMill/Plate_Solver_Package
#   sbatch submit_rect_ff_persist_golden_polish_2026-09-07.sh
#
# Runtime: 3 candidates x 25 golden-section iterations, ~150 sigma calls
# total -- expect well under 2 minutes; 20-minute SBATCH limit is
# generous headroom.
#
set -euo pipefail
cd /home/ghmkfh/PythonMill/Plate_Solver_Package

module load python/3.12.1
source venv/bin/activate
export PKG_PATH="."
export DPS="${DPS:-30}"
export EXPECT_SOLVER_VERSION=2026-07-10.s10
export PYTHONUNBUFFERED=1
python3 -u probe_rect_ff_persist_golden_polish_2026-09-07.py
