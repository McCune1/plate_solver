#!/bin/bash
#SBATCH --job-name=rectff_relgate
#SBATCH --output=rect_ff_release_gate_%j.out
#SBATCH --error=rect_ff_release_gate_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#
# Rectangular FFFF release gate. Confirms, on the cluster's own deployed
# tree, the four things Paper 2 review pass #3 showed the .tex-only
# currency check cannot see:
#
#   G1  the free-free corner term is the closed-contour checkerboard
#       (TT - TB - WT + WB) and NOT the naive same-sign sum, which
#       vanishes identically under the parity identity and therefore
#       fails silently rather than raising;
#   G2  RectIPAssembler actually exposes bc='free_free';
#   G3  the two stored single-point sigma_min anchors reproduce exactly
#       (l/b=1.5 SYM Lambda*=0.964 and l/b=2.0 ANTI Omega*=0.5200), at
#       dps 26/30/40;
#   G4  in-plane free-free and clamped-free really do assemble
#       differently at the same point.
#
# The released github_repo/ copy had drifted a month behind on exactly
# one file (core_solvers.py) and shipped the pre-2026-09-03 assembler.
# Run this before cutting the Paper 2 manuscript-matching tag.
#
# Copy/paste on the cluster:
#   cd /home/ghmkfh/PythonMill/Plate_Solver_Package
#   sbatch submit_rect_ff_release_gate_2026-09-07.sh
#
# The G5 release-vs-deployed md5 comparison is skipped here (the release
# tree lives on the workstation, not the cluster). Run it locally with:
#   PKG_PATH=<...>/Plate_Solver_Package \
#   RELEASE_PKG=<...>/github_repo \
#   python3 probe_rect_ff_release_gate_2026-09-07.py
#
# Read-only. Does NOT retune PERSIST_DL or the 0.744 MAC bar, does NOT
# edit plate_solver, does NOT bump SOLVER_VERSION, does NOT touch a
# tabulated frequency. A FAIL_ANCHOR verdict means a code change altered
# a computed rectangular FFFF frequency -- bump SOLVER_VERSION and
# re-derive the tables; do not edit the stored anchors to make it pass.
#
# Runtime: six single-point sigma_min probes plus source checks. Expect
# well under a minute; the 15-minute limit is headroom.
#
set -euo pipefail
cd /home/ghmkfh/PythonMill/Plate_Solver_Package

module load python/3.12.1
source venv/bin/activate
export PKG_PATH="."
export EXPECT_SOLVER_VERSION=2026-07-10.s10
export DPS_LIST="26,30,40"
export RTOL=1e-12
export PYTHONUNBUFFERED=1
python3 -u probe_rect_ff_release_gate_2026-09-07.py
