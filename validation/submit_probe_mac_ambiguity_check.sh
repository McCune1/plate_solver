#!/bin/bash
#SBATCH --job-name=mac_ambiguity_check
#SBATCH --output=mac_ambiguity_check_%j.out
#SBATCH --error=mac_ambiguity_check_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=00:30:00
#SBATCH --partition=general

# Second half of the MAC ambiguity-resolution pair (first half:
# ansys_ff_oop_mac_extract.inp / submit_ansys_ff_oop_mac_extract.sh, which
# MUST be run first, on an ANSYS-licensed allocation). As of 2026-07-15 that
# submit script auto-copies its mac_eigvec_mesh96.txt output up to this
# package root itself -- no manual copy step should be needed before this
# job is submitted (see LESSONS_LEARNED Sec 13.5). This script also checks
# the Ansys/ subfolder directly as a defense-in-depth fallback (older Ansys
# submit script, or a copy that failed) before giving up; if the file truly
# can't be found anywhere, it fails fast with a clear message instead of
# guessing.
#
# Computes the Modal Assurance Criterion between the FE eigenvector at
# Omega_lit=111.4627 and each of the two open OOP ceiling-ambiguity
# candidates (raw Omega=2.828013, 2.838515 -- LESSONS_LEARNED Sec 20.4),
# the sixth and last-identified lever after five prior techniques (basis
# bump, sensitivity map, null-vector gap, parity readout, ODD-block
# topology scan) failed to force the assignment. Pre-registered
# interpretation criteria are printed inline by the probe itself; see its
# docstring for the full pre-registration and cost accounting (cheap:
# ~90-150s per candidate reconstruction, single CPU sufficient, 2 requested
# only so both candidates run without contention).
#
# Diagnostic only -- no SOLVER_VERSION action, no package changes. Does
# NOT by itself reopen or re-decide the closed Sec 20.4 ambiguity record;
# a human must read the printed MAC numbers against the pre-registered
# criteria before updating LESSONS_LEARNED or either paper/report document.

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export BLIS_NUM_THREADS=1

export PKG_PATH=/home/ghmkfh/PythonMill/Plate_Solver_Package
export EXPECT_SOLVER_VERSION=2026-07-10.s10
# FE_EIGVEC_PATH intentionally NOT set here -- leaving it unset lets
# probe_mac_ambiguity_check.py's own candidate search (package root, then
# Ansys/ subfolder) run; setting it here would force one exact path with no
# fallback, defeating that search. Set it manually only if you need to point
# at a nonstandard location.
#
# MAC_N_DOFS / MAC_XI_MAX: uncomment BOTH to re-run at a larger basis
# (LESSONS_LEARNED Sec 20.4 next step, added 2026-07-15 after job 2324018's
# marginal 0.9913-vs-0.9998 MAC gap between the two ODD candidates) --
# defaults (unset) are 20/20.0, matching every prior MAC job. MUST move
# together (project's standard "n_dofs+8/xi_max+6" convention) -- job
# 2324023 set MAC_N_DOFS=28 alone and hard-failed ("only 24/28 branches
# filled") because the old xi_max=20.0 didn't contain enough branches to
# fill the larger basis. A widening gap at n_dofs=28 would support cand_B
# (raw 2.838515); a flat gap would support genuine near-degeneracy instead.
# export MAC_N_DOFS=28
# export MAC_XI_MAX=26.0

module load python/3.12.1
source /home/ghmkfh/PythonMill/Plate_Solver_Package/venv/bin/activate

cd /home/ghmkfh/PythonMill/Plate_Solver_Package/

PROBE=""
if [ -f validation/probe_mac_ambiguity_check.py ]; then
  PROBE=validation/probe_mac_ambiguity_check.py
elif [ -f probe_mac_ambiguity_check.py ]; then
  PROBE=probe_mac_ambiguity_check.py
else
  echo "FATAL: probe_mac_ambiguity_check.py not found in validation/ or $(pwd)."
  exit 2
fi
echo "probe: $PROBE"

# Defense-in-depth: if the Ansys submit script's auto-copy didn't run (older
# copy, or the copy step itself failed) but the file exists in Ansys/, stage
# it here before handing off to Python.
if [ ! -f mac_eigvec_mesh96.txt ] && [ -f Ansys/mac_eigvec_mesh96.txt ]; then
  echo "mac_eigvec_mesh96.txt not found at package root but found in Ansys/ -- copying up."
  cp -v Ansys/mac_eigvec_mesh96.txt .
fi

if [ ! -f mac_eigvec_mesh96.txt ] && [ ! -f Ansys/mac_eigvec_mesh96.txt ]; then
  echo "FATAL: mac_eigvec_mesh96.txt not found at the package root or in Ansys/."
  echo "  Run submit_ansys_ff_oop_mac_extract.sh (from Ansys/) first, on an"
  echo "  ANSYS-licensed allocation. Its submit script auto-copies the output"
  echo "  here, so if it's still missing after that job completed, check its"
  echo "  /COM lines and .out log for a mode-selection or UZ-dominance"
  echo "  failure -- the file may simply never have been produced."
  exit 3
fi

echo "Job $SLURM_JOB_ID on $(hostname)"
echo "Start: $(date)"

python3 -u "$PROBE"

echo "End: $(date)"
