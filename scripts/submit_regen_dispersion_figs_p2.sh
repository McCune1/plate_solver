#!/bin/bash
#SBATCH --job-name=regen_dispersion_figs_p2
#SBATCH --output=regen_dispersion_figs_p2_%j.out
#SBATCH --error=regen_dispersion_figs_p2_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=00:30:00
#SBATCH --partition=general

# Same as submit_regen_dispersion_figs.sh (see that file's header for the
# full rationale) but for the p2 (in-plane) checkpoint/figdir pair:
# 2*Theta=0.25/0.5/1.0 pi at r0/(2b)=1.25, figures_p2/. This is the ONLY
# figure-regen job that touches Part 2 entries -- plot_mode_shape/
# mode_shape_grid (Fig. 3) is Part-1-only (see probe_regen_modeshape_figs.
# py's SCOPE note), so Part 2 geometries only ever had *_fig2.png to begin
# with; this script is their complete fix.
#
# PUSH BEFORE RUNNING: plate_solver/plotting.py (the 2026-07-16 fix) and
# this script + probe_regen_dispersion_figs.py.
# PULL AFTER RUNNING: the .out log only, to confirm every geometry OK with
# no FAILED lines -- PNGs are overwritten in place in figures_p2/.

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export BLIS_NUM_THREADS=1

export PKG_PATH=/home/ghmkfh/PythonMill/Plate_Solver_Package
export EXPECT_SOLVER_VERSION=2026-07-10.s8
export CHECKPOINT=research23_checkpoint_p2.json
export FIGDIR=figures_p2
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
