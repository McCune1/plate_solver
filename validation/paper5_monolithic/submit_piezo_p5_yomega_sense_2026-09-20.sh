#!/bin/bash
#SBATCH --job-name=p5_yomega_sense
#SBATCH --output=piezo_p5_yomega_sense_2026-09-20_%j.out
#SBATCH --error=piezo_p5_yomega_sense_2026-09-20_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:20:00
#SBATCH --mem=4G
#SBATCH --partition=requeue
#SBATCH --requeue

# Paper 5 roadmap Sec 8 item 6, last of the four optional-cheap extras.
# The voltage-driven Y(omega) was proven algebraically TRIVIAL this
# session (LESSONS_LEARNED.md Sec 18.210: Y=j*omega*C0, k_eff^2==0, a
# direct corollary of the closed Phase 0 OC==SC theorem). Redirected,
# by the user's explicit choice, to a force-driven sensing admittance:
# a harmonic ring load on the F-F SC ring, charge read off a SEGMENT
# of the still-grounded face (the full-face charge is ALSO identically
# zero -- the reciprocal finding, Sec 18.211). Derivation handed to
# Grok (PAPER5_YOMEGA_SENSE_DERIVATION.md); this probe re-derives and
# re-checks its three validation gates independently rather than
# trusting the derivation doc at face value, plus a fourth diagnostic
# gate (G_mu) that documents a real bug found in an earlier interactive
# probe this session (scaling e31 without e33 does not actually drive
# e31_bar -> 0 at this material point -- see the probe's own docstring
# and LESSONS_LEARNED.md Sec 18.21x for the full account).
#
# New solver code: plate_solver/piezo_monolithic.py gained
# _equilibrate_solve (row+column-equilibrated mp.lu_solve, needed
# because the driven 12x12's third chi-cubic branch produces Bessel
# arguments large enough that raw entries span >50 orders of
# magnitude) plus _driven_vecs_at/driven_ff_force_sc/_driven_eval/
# Q_segment/Y_sense. plate_solver/workers.py gained
# _piezo_mono_driven_force_sc_worker (a point evaluation, not a root
# bisection -- returns (Re, Im) as two floats). SOLVER_VERSION is NOT
# bumped (new unused-by-default entry points, same pattern as every
# other method in this module).
#
# See probe_piezo_p5_yomega_sense_2026-09-20.py's own docstring for
# the full G1/G2/G3/G_mu/G_worker pre-registration.
#
# In-sandbox timing this session (dps=60, dps_hi=80, inner_iters=60):
# ~107s end to end, SENTINEL PASS_ALL. G1 identity(G) rel_err~1e-61;
# G2 pole rel_err=0.0 (coupled_bisect and the driven-12x12 homogeneous
# determinant agree to all computed digits); G3 log-log slope of Q_in
# vs e31_bar -> 0.999993 (linear vanishing); G_mu mechanical branches
# match +-elastic-mu to rel_err=0.0 once e31 and e33 are scaled
# together correctly; G_worker bit-identical at a generic point and
# near the pole. Y_sense sweep shows the expected large-magnitude
# sign flip through the F-F SC fundamental (473.2496...).
# 1 CPU is enough; 20 min is generous headroom.

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMBA_NUM_THREADS=1

echo "job=$SLURM_JOB_ID host=$(hostname) start=$(date)"

module load python/3.12.1
source /home/ghmkfh/PythonMill/Plate_Solver_Package/venv/bin/activate
cd /home/ghmkfh/PythonMill/Plate_Solver_Package/Paper5_Monolithic

export PKG_PATH=/home/ghmkfh/PythonMill/Plate_Solver_Package
export EXPECT_SOLVER_VERSION=2026-07-10.s10
export PIEZO_DPS=60
export PIEZO_DPS_HI=80
export PIEZO_INNER_ITERS=60

if [ ! -f probe_piezo_p5_yomega_sense_2026-09-20.py ]; then
  echo "FATAL: probe_piezo_p5_yomega_sense_2026-09-20.py not found in $(pwd)"
  exit 2
fi
if ! grep -q "def Y_sense" "$PKG_PATH/plate_solver/piezo_monolithic.py"; then
  echo "FATAL: $PKG_PATH/plate_solver/piezo_monolithic.py has no Y_sense -- push the updated file first"
  exit 2
fi
if ! grep -q "_piezo_mono_driven_force_sc_worker" "$PKG_PATH/plate_solver/workers.py"; then
  echo "FATAL: $PKG_PATH/plate_solver/workers.py has no _piezo_mono_driven_force_sc_worker -- push the updated file first"
  exit 2
fi

python3 -u probe_piezo_p5_yomega_sense_2026-09-20.py
echo "end=$(date)"
