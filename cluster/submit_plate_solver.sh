#!/bin/bash
#SBATCH --job-name=plate_solver
#SBATCH --output=plate_solver_%j.out
#SBATCH --error=plate_solver_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=64
#SBATCH --time=10:00:00
#SBATCH --partition=general

# ──────────────────────────────────────────────────────────────────────────────
#  Annular & Rectangular Sector Plate Vibration Solver — plate_solver package
#  The Mill, Missouri University of Science and Technology
#
#  2026-07-01: replaces submit_research42.sh. Research50.py/rect_int.py were
#  merged and split into the plate_solver/ package (see LESSONS_LEARNED.md
#  Sec. 22). Every env var below is UNCHANGED from submit_research42.sh --
#  the module split and the new argparse layer (plate_solver/cli.py) are both
#  purely additive; nothing here needed to change except the two lines that
#  actually invoke Python. If you diff this against submit_research42.sh, the
#  only differences should be the header comment and those two lines.
#
#  Pipeline (run_overnight):
#    (1) cantilever regression spot-check (Step-3/6 neutrality, ~3 min),
#    (2) FREE-FREE spectra sweep (both radial edges free, 64-core parallel),
#        tabulated in the native Seok-Tiersten parameter and the literature
#        parameter Omega_lit = w*R_o^2*sqrt(rho*H/D),
#    (3) self-checking free-free literature comparison, and
#    (4) orthotropic-material scaffold + isotropic-reduction sanity test.
#
#  Deploy: copy the plate_solver/ directory (the package, not just this
#  script) to the same directory as this submit script, or `pip install -e .`
#  it into the SLURM job's venv.
#  Per-geometry checkpoint figures/freefree_checkpoint.pkl => re-submit RESUMES.
# ──────────────────────────────────────────────────────────────────────────────

# Thread-pinning: one BLAS thread/process; parallelism is the mp ProcessPool
# (LESSONS §8 — threaded BLAS oversubscribes every worker and starves the run).
# NOTE: plate_solver/config.py now sets these with a DIRECT assignment (not
# setdefault), so it wins even if SLURM injected different values -- but
# exporting them here first is still good practice and costs nothing.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export BLIS_NUM_THREADS=1

# One worker per allocated core (LESSONS §8 — size from SLURM, not cpu_count()).
export N_WORKERS=${SLURM_CPUS_PER_TASK:-64}

# ── Overnight driver controls (names match the old Research42.py: R40_*) ──────
export R40_OVERNIGHT=1        # run the all-steps driver
export R40_SPOTCHECK=1        # (1) cantilever regression spot-check
export R40_FREEFREE=1         # (2) free-free spectra sweep
export R40_ORTHO_TEST=1       # (4) orthotropic scaffold sanity
export R40_FF_MODES=6         # free-free modes wanted per geometry

# Free-free geometry lists, "r0_2b:twoTheta_over_pi" comma-separated.
export R40_FF_P1=1.25:0.5,1.25:1.0,1.25:1.5,1.25:0.25,1.25:0.75,1.25:1.25
export R40_FF_P2=1.25:0.5,1.25:1.0,1.25:0.25

# ── Material (default = ν=0.35 validation material; reproduces all prior runs)─
# export MAT_E=210e9
# export MAT_NU=0.35
# export MAT_RHO=7800.0

# ------------------------------------------------------------------------------
#  STEP 1 (orthotropic) -- SHI/LV 2016 FFFF DIRECT ANNULAR BENCHMARK (R40_SHI=1)
#  Uses Research50's CORRECTED material defaults (G=3.51e9, rho=7850) and the
#  real 8-mode Ansys reference -- NOT rect_int.py's stale G=7.3e9/rho=7800 +
#  fabricated table, which was discarded (not kept) in the plate_solver merge.
#  See LESSONS_LEARNED Sec. 21.3/22. Single multi-hour OOP solve (~40-57 min
#  on 64 cores); resumes from figures/shi_checkpoint.pkl.
# ------------------------------------------------------------------------------
# export R40_SHI=1
# export R40_FF_MODES=11

# Figures off for the spectra sweep (values are written to figures/freefree_results.txt).
export MAKE_FIGS=0

# Quick shake-out (~30-40 min): halves every scan range/point count.
# export FAST=1

module load python/3.12.1
source /home/ghmkfh/PythonMill/Plate_Solver_Package/venv/bin/activate

echo "Job $SLURM_JOB_ID on $(hostname): cpus-per-task=$SLURM_CPUS_PER_TASK, N_WORKERS=$N_WORKERS, OMP=$OMP_NUM_THREADS"
echo "Start: $(date)"

python3 -u -m plate_solver.cli

echo "End: $(date)"
