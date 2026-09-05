#!/bin/bash
#SBATCH --job-name=shi_ffff
#SBATCH --output=shi_ffff_%j.out
#SBATCH --error=shi_ffff_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=64
#SBATCH --time=02:00:00
#SBATCH --partition=general
#
# ──────────────────────────────────────────────────────────────────────────────
#  Orthotropic annular sector FFFF — Shi, Liang, Wang & Teng (2016),
#  J. Vibroengineering 18(5), paper 2111, DOI 10.21595/jve.2016.17004.
#
#  Dedicated launcher for R40_SHI=1 (was bundled/commented-out inside
#  submit_research42.sh; broken out into its own script, same pattern as
#  submit_rect_ip.sh, per LESSONS_LEARNED Sec. 21.4/22 — a single multi-hour
#  free-free OOP solve doesn't need the full overnight driver's 10h walltime
#  or its RUN_P1/RUN_P2/free-free-sweep machinery around it).
#
#  Uses the plate_solver package's CORRECTED material defaults (G=3.51e9,
#  rho=7850, real 8-mode Ansys reference) via `--shi` / R40_SHI=1 --
#  Research50.py's fix, NOT rect_int.py's stale G=7.3e9/rho=7800 + fabricated
#  reference table (see LESSONS_LEARNED Sec. 21.3/22.1). Validates the
#  polar-orthotropic (T,R,mu_theta) derivation and the beta-term fix against
#  an independent Ansys SHELL281 benchmark; does NOT touch the isotropic
#  path (isotropic frequencies are byte-identical regardless of this run).
#
#  Single multi-hour free-free OOP solve (~40-57 min on 64 cores); resumes
#  from figures/shi_checkpoint.pkl on re-submit.
#
#  Deploy: copy the plate_solver/ package directory alongside this script.
# ──────────────────────────────────────────────────────────────────────────────

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export BLIS_NUM_THREADS=1

export N_WORKERS=${SLURM_CPUS_PER_TASK:-64}

# ── Shi/Lv validation controls (all optional -- defaults = paper section 3) ──
export R40_FF_MODES=11            # solve 11 modes; first 3 (~0) are rigid-body
export SHI_OM_FLOOR=1.0           # drop omega_lit below this (rigid-body cut)
# material overrides (defaults = CORRECTED: G=3.51e9, rho=7850):
# export SHI_E_R=40e9 SHI_E_THETA=70e9 SHI_NU_R=0.3 SHI_G_RTHETA=3.51e9 SHI_RHO=7850
# geometry overrides (defaults = Shi a/b=0.5, phi=90deg):
# export SHI_R0_2B=1.5 SHI_TWO_T=0.5
# override the built-in Ansys reference spectrum (comma-separated), if a
# mesh-converged re-run produces updated target values:
# export VAL_REF=12.40,15.57,36.13,43.65,66.76,87.46,89.82,95.35

export MAKE_FIGS=0

echo "=== DIAGNOSTIC ==="
echo "--- module command available? ---"
type module 2>&1

echo "--- module avail python (before load) ---"
module avail python 2>&1

echo "--- module load python/3.12.1 ---"
module load python/3.12.1
echo "exit code: $?"
module list 2>&1

echo "--- which python3 / version (after module load, before venv) ---"
which python3
python3 --version

echo "--- venv pyvenv.cfg ---"
cat /home/ghmkfh/PythonMill/Plate_Solver_Package/venv/pyvenv.cfg 2>&1

echo "--- venv bin contents ---"
ls -la /home/ghmkfh/PythonMill/Plate_Solver_Package/venv/bin/ 2>&1

echo "--- activating venv ---"
source /home/ghmkfh/PythonMill/Plate_Solver_Package/venv/bin/activate
echo "VIRTUAL_ENV=$VIRTUAL_ENV"
which python3
python3 --version
echo "=== END DIAGNOSTIC ==="


module load python/3.12.1
source /home/ghmkfh/PythonMill/Plate_Solver_Package/venv/bin/activate

echo "Job $SLURM_JOB_ID on $(hostname): cpus-per-task=$SLURM_CPUS_PER_TASK, N_WORKERS=$N_WORKERS"
echo "Start: $(date)"

python3 -u -m plate_solver.cli --shi

echo "End: $(date)"
