#!/bin/bash
#SBATCH --job-name=rect_ff_oop_mac_anti
#SBATCH --output=rect_ff_oop_mac_anti_%j.out
#SBATCH --error=rect_ff_oop_mac_anti_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#
# OPTIONAL. Seven persist-basis golden-section polishes of remaining
# Table-1 ANTI MAC misses. Not required for the paper fold-in already
# in github_repo/paper/. Does not touch ANSYS (no extra license seat).
#
#   cd /home/ghmkfh/PythonMill/Plate_Solver_Package
#   sbatch submit_rect_ff_oop_mac_anti_remain_2026-09-06.sh
#
set -euo pipefail
cd /home/ghmkfh/PythonMill/Plate_Solver_Package
module load python/3.12.1
source venv/bin/activate
export PKG_PATH="."
export DPS=40
export EXPECT_SOLVER_VERSION=2026-07-10.s10
export PYTHONUNBUFFERED=1
python3 -u probe_rect_ff_oop_mac_anti_remain_2026-09-06.py
