#!/bin/bash
#SBATCH --job-name=persist_finegrid
#SBATCH --output=rect_ff_persist_finegrid_recheck_%j.out
#SBATCH --error=rect_ff_persist_finegrid_recheck_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#
# Tiny, surgical follow-up to job 2457815
# (rect_ff_persist_threshold_sensitivity_2026-09-07), which scored all 148
# of Paper 2's published production matches against the frozen Screen-B
# persistence cut (PERSIST_DL=0.01). All 148 persist, as required -- but
# exactly THREE sit at the outer edge (delta=0.01000, the max the coarse
# 0.002-step window scan can report before failing):
#
#   OOP l/b=1.5 SYM  Lambda*=5.380
#   OOP l/b=2.5 ANTI Lambda*=6.260
#   IP  l/b=3.0 ANTI Omega*=2.2400
#
# Because persist_one/persist_one_ip's window scan is quantized to 0.002
# steps, "delta=0.01000" for these three is ambiguous: it could be the
# true nearest dip, or a closer dip the coarse grid stepped over. This
# job rescans ONLY these three candidates' identical +/-0.02 window at a
# 4x finer 0.0005 step, using the same unmodified, already-public
# p5_rect_ff_lib building blocks (sigma_at/sigma_at_ip, local_minima,
# frange, n_basis/n_basis_ip) -- same production basis and im_cap as
# persist_one/persist_one_ip use internally. Not a monkeypatch; does not
# touch persist_one/persist_one_ip, does not retune PERSIST_DL, does not
# touch a tabulated frequency, no SOLVER_VERSION bump, no paper edit.
#
# A per-candidate fidelity gate (coarse rescan must reproduce job
# 2457815's own delta before the fine-grid result is trusted) runs first
# for each of the three.
#
# Copy/paste on the cluster:
#   cd /home/ghmkfh/PythonMill/Plate_Solver_Package
#   sbatch submit_rect_ff_persist_finegrid_recheck_2026-09-07.sh
#
# Runtime: 3 candidates x 81-point window = 243 sigma calls, expect well
# under a minute total; 15-minute SBATCH limit is generous headroom.
#
set -euo pipefail
cd /home/ghmkfh/PythonMill/Plate_Solver_Package

module load python/3.12.1
source venv/bin/activate
export PKG_PATH="."
export DPS="${DPS:-30}"
export EXPECT_SOLVER_VERSION=2026-07-10.s10
export PYTHONUNBUFFERED=1
python3 -u probe_rect_ff_persist_finegrid_recheck_2026-09-07.py
