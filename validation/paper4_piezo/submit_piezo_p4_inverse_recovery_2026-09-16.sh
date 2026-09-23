#!/bin/bash
#SBATCH --job-name=p4_piezo_inverse
#SBATCH --output=piezo_p4_inverse_recovery_2026-09-16_%j.out
#SBATCH --error=piezo_p4_inverse_recovery_2026-09-16_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=00:30:00
#SBATCH --mem=4G
#SBATCH --partition=requeue
#SBATCH --requeue

# Paper 4 (piezoelectric ring), PAPER4_PIEZO_ROADMAP.md Sec 6/Sec 8 step 4:
# the inverse/parameter-estimation demonstration, following Paper 3's
# Southwell inverse-b/a methodology (LESSONS_LEARNED.md Sec 18.119) --
# self-consistency round-trip (G0), root-quality gate with a genuine
# negative control (G1), on-branch sensitivity (G2, reported not gated) --
# on the now-complete, worker-path-verified C-C piezo-coupled forward
# model (LESSONS_LEARNED.md Sec 18.171-18.173). Two independent single-
# parameter recoveries: host elastic modulus E, and piezoelectric stress
# constant e31, each with the other material constants held known -- per
# the roadmap's own explicit recommendation to do these separately before
# any joint/multi-parameter fit.
#
# PRE-REGISTERED EXPECTATION (see the .py's own docstring for the full
# reasoning and this session's sandbox sweep that motivated it): E is
# expected to be well-conditioned (comparable relative uncertainty to the
# assumed frequency-measurement precision); e31 is expected to be
# essentially UNIDENTIFIABLE from this single short-circuit mode alone
# (LESSONS_LEARNED Sec 18.141 F / Duan2005 Fig. 3's own finding that
# short-circuit coupling is only a ~0.01-0.016% correction on top of
# ordinary added-layer stiffness) -- a large propagated e31 uncertainty is
# the correctly-anticipated finding, not a probe failure.
#
# PRE-REGISTERED PASS/FAIL: G0 (<1e-6 relative self-consistency, exact
# synthetic target, no injected measurement noise) AND G1 (genuine sign
# flip + negative control) must BOTH pass for BOTH parameters for
# SENTINEL PASS_ALL -- this is a pure numerical round-trip/root-quality
# check and is expected to pass tightly regardless of the very different
# physical conditioning between E and e31. G2's sensitivity/propagated-
# uncertainty numbers are reported either way, not gated -- read
# piezo_p4_inverse_recovery_results.json for the full per-parameter
# detail, and the SENTINEL line for the plain PASS_ALL/FAIL_ALL verdict.
#
# Cost model (this session's own sandbox timing at dps=100, reduced-
# iteration smoke test on the actual deployed package): ~0.085s per
# determinant evaluation at this mode; full defaults (INNER_ITERS=40,
# OUTER_ITERS=26) put each parameter at roughly 90 inner-bisect-equivalent
# calls, ~5-6 minutes; two parameters run concurrently
# (ProcessPoolExecutor, max_workers=2) so wall time should not exceed
# ~10 minutes. The 30-minute budget here is generous headroom, matching
# this project's own stated practice for small jobs, not a tight
# estimate -- some evaluation points elsewhere in this piezo-coupled
# family have run 5-6x slower than this mode's typical cost (Sec 18.148's
# table4 gate saw up to ~25s/point at some ratios), so the margin is
# deliberate, not padding.
#
# Submit from
# /home/ghmkfh/PythonMill/Plate_Solver_Package/Paper4_Piezo
# (the venv is one level up and still sourced from there -- only the cwd
# for running the .py and writing its outputs moves into this subfolder).
# Deliberately imports plate_solver.piezo_solver.PiezoOutOfPlaneSolver
# directly (unlike most of this project's standalone-mpmath probes) --
# see the .py's own "Standalone-vs-package note" for why that's safe here
# (PiezoOutOfPlaneSolver is not on the SOLVER_VERSION checkpoint-
# compatibility path at all).

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMBA_NUM_THREADS=1

echo "job=$SLURM_JOB_ID host=$(hostname) start=$(date)"

module load python/3.12.1
source /home/ghmkfh/PythonMill/Plate_Solver_Package/venv/bin/activate
cd /home/ghmkfh/PythonMill/Plate_Solver_Package/Paper4_Piezo

export PKG_PATH=/home/ghmkfh/PythonMill/Plate_Solver_Package
export EXPECT_SOLVER_VERSION=2026-07-10.s10
export N_WORKERS=2
export PIEZO_DPS=100
export PIEZO_INNER_ITERS=40
export PIEZO_OUTER_ITERS=26

if [ ! -f probe_piezo_p4_inverse_recovery_2026-09-16.py ]; then
  echo "FATAL: probe_piezo_p4_inverse_recovery_2026-09-16.py not found in $(pwd)"
  exit 2
fi

python3 -u probe_piezo_p4_inverse_recovery_2026-09-16.py
echo "end=$(date)"
