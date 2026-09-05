#!/bin/bash
#SBATCH --job-name=regen_dispersion_figs
#SBATCH --output=regen_dispersion_figs_%j.out
#SBATCH --error=regen_dispersion_figs_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=00:30:00
#SBATCH --partition=general

# Regenerates every *_fig2.png dispersion-curve figure in figures/ from
# research23_checkpoint.json's geometry list (2*Theta=0.25/0.5/0.75/1.0/
# 1.25 pi at r0/(2b)=1.25), picking up the 2026-07-16 plotting.py fix:
# compute_dispersion_branches' root search was cheapened (ngrid=6/nax=140
# vs. full_search's own ngrid=10/nax=240 default) purely for speed, and at
# that density it sometimes failed to resolve a branch's root at one Omega
# sample while finding it fine on neighbouring samples -- the branch
# tracker then ended that branch right there, producing visible gaps in
# the plotted line not present in Seok & Tiersten's own Fig. 2. Fixed by
# restoring the search density to full_search's own default and doubling
# the Omega-sample count (24->48).
#
# This does NOT re-run the sigma_min discovery sweep that originally found
# each geometry's tabulated frequencies (hours per geometry) --
# plot_dispersion_curves never depends on any found frequency, only the
# geometry/material and the cut-off-derived Omega window, and uses the
# fast float64 engine throughout, not the mp.dps=40 discovery path. Cost
# is one dispersion sweep per geometry, a few minutes each -- the whole
# job should finish in well under the 30-minute budget.
#
# Independent of (can run in parallel with): everything else queued.
# Covers Part 1 only (this checkpoint has no Part 2 entries) -- see
# submit_regen_dispersion_figs_p2.sh for the in-plane geometries and
# submit_regen_dispersion_figs_p1tail.sh for the p1tail checkpoint.
#
# PUSH BEFORE RUNNING: plate_solver/plotting.py (this fix) and this script
# + probe_regen_dispersion_figs.py.
# PULL AFTER RUNNING: the .out log only, to confirm every geometry OK with
# no FAILED lines -- PNGs are overwritten in place in figures/.

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export BLIS_NUM_THREADS=1

export PKG_PATH=/home/ghmkfh/PythonMill/Plate_Solver_Package
export EXPECT_SOLVER_VERSION=2026-07-10.s8
export CHECKPOINT=research23_checkpoint.json
export FIGDIR=figures
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
