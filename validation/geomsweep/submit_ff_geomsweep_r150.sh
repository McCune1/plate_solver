#!/bin/bash
#SBATCH --job-name=ff_geomsweep_r150
#SBATCH --output=ff_geomsweep_r150_%j.out
#SBATCH --error=ff_geomsweep_r150_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=64
#SBATCH --time=32:00:00
#SBATCH --partition=general

# ──────────────────────────────────────────────────────────────────────────
#  Free-free GEOMETRIC BENCHMARKING sweep -- new radius ratio r0/(2b)=1.5
#
#  2026-07-19: extends the free-free validation (closed 2026-07-12/13 at
#  ONE geometry, r0/(2b)=1.25, "FF-P1") to a second radius ratio, mirroring
#  how Seok & Tiersten's own Table 2 covers MULTIPLE (r0/(2b), 2Theta/pi)
#  pairs, not a single geometry -- the goal being to check the free-free
#  boundary condition (FreeFreeOOP/FreeFreeIP, boundary.py) generalizes
#  across geometry rather than being an artifact of the one pinned case.
#
#  NO NEW PYTHON REQUIRED: run_overnight()/_solve_freefree() (validation.py)
#  already derive their scan windows PURELY from the xi=0/zeta=0 cut-off
#  frequencies of the geometry actually passed in (_scan_cfg_part1/2) -- no
#  paper-table lookups, so any (r0_2b, two_T_pi) pair is already "fully
#  general" by design. This script only widens the R40_FF_P1/R40_FF_P2
#  geometry-list env vars cli.py -> run_overnight() already reads.
#
#  Uses its own FIGDIR/checkpoint (figures_ff_r150/) so it can run
#  concurrently with the other new-ratio jobs in this folder without
#  clobbering their checkpoints, and independently of the original FF-P1
#  figures/ directory.
#
#  Angle set matches the existing FF-P1 P1 angle grid (0.25 .. 1.5 * pi)
#  exactly, run for BOTH Part 1 (OOP) and Part 2 (IP) -- the original
#  free-free work only ran 3 of these 6 angles for Part 2; this sweep runs
#  the full 6 for both parts at the new ratio for symmetric coverage.
#
#  COST: submit_plate_solver.sh (9 geometries: 6 P1 + 3 P2, r0_2b=1.25)
#  budgets 10h on 64 cores. This job runs 12 geometries (6 P1 + 6 P2) at a
#  DIFFERENT ratio -- no prior timing data for r0_2b=1.5 specifically, so
#  --time is set to 16h as a conservative ~1.3x scale-up. RESUME=1 (default)
#  means a re-submit after a timeout picks up from the checkpoint rather
#  than re-paying for completed geometries.
# ──────────────────────────────────────────────────────────────────────────

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export BLIS_NUM_THREADS=1

export N_WORKERS=${SLURM_CPUS_PER_TASK:-64}

export R40_OVERNIGHT=1
export R40_SPOTCHECK=1
export R40_FREEFREE=1
export R40_ORTHO_TEST=1
export R40_FF_MODES=6

# r0/(2b) = 1.5, all 6 angles, both parts.
export R40_FF_P1=1.5:0.25,1.5:0.5,1.5:0.75,1.5:1.0,1.5:1.25,1.5:1.5
export R40_FF_P2=1.5:0.25,1.5:0.5,1.5:0.75,1.5:1.0,1.5:1.25,1.5:1.5

export MAKE_FIGS=0
export FIGDIR=figures_ff_r150
export CHECKPOINT=research23_checkpoint_ff_r150.json
export RESUME=1

module load python/3.12.1
source /home/ghmkfh/PythonMill/Plate_Solver_Package/venv/bin/activate

cd /home/ghmkfh/PythonMill/Plate_Solver_Package/
mkdir -p "$FIGDIR"

echo "Job $SLURM_JOB_ID on $(hostname): cpus-per-task=$SLURM_CPUS_PER_TASK, N_WORKERS=$N_WORKERS"
echo "Start: $(date)"

python3 -u -m plate_solver.cli

echo "End: $(date)"
