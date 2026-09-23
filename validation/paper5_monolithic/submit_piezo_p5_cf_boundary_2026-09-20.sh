#!/bin/bash
#SBATCH --job-name=p5_cf_bc
#SBATCH --output=piezo_p5_cf_boundary_2026-09-20_%j.out
#SBATCH --error=piezo_p5_cf_boundary_2026-09-20_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:20:00
#SBATCH --mem=4G
#SBATCH --partition=requeue
#SBATCH --requeue

# Paper 5 roadmap Sec 8 item 6: C-F boundary condition. This one adds
# real new solver code (plate_solver/piezo_monolithic.py gained
# _elastic_mixed_det/_coupled_mixed_det + elastic_cf/fc_det/bisect +
# cf/fc_coupled_det/bisect, plus four new worker functions in
# plate_solver/workers.py), unlike the pure-probe overtones/thickness-
# sweep extras (Sec 18.206/18.207). SOLVER_VERSION is NOT bumped (new
# unused-by-default entry points, same pattern as every other method
# in this module).
#
# See probe_piezo_p5_cf_boundary_2026-09-20.py's own docstring for the
# full G_reduce/G_found/G_theorem/G1/G_worker pre-registration.
# G_reduce is the load-bearing gate: the new mixed-edge row-selection
# code must reproduce the already-validated elastic_det/elastic_cc_det/
# coupled_det/cc_coupled_det determinant values BIT-IDENTICALLY when
# both edges get the same label, which is a much stronger correctness
# proof than a fresh derivation would give.
#
# In-sandbox timing this session (dps=35, inner_iters=35, scan_n=200):
# ~33s end to end, SENTINEL PASS_ALL. C-F: el=1605.44 sc=1644.84
# split=2.454%. F-C: el=556.88 sc=566.60 split=1.745%. Both stiffen
# (SC > elastic), generalizing the Sec 18.199 theorem to mixed edges.
# Genuine finding (not predicted in advance): the C-F elastic
# fundamental sits 90.4% of the way from the F-F fundamental to the
# C-C fundamental, while F-C sits only 7.5% of the way there -- for
# this r_i=0.1/r_o=0.6 geometry, the INNER edge condition dominates
# bending stiffness far more than the outer edge condition (plausible
# mechanism: the n=0 moment condition's -(2*A1/r)*dZ term is ~6x more
# sensitive at r_i than at r_o, since r_o/r_i=6 here -- reported as a
# hypothesis, not asserted as proven).
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
export PIEZO_INNER_ITERS=45
export PIEZO_SCAN_N=600

if [ ! -f probe_piezo_p5_cf_boundary_2026-09-20.py ]; then
  echo "FATAL: probe_piezo_p5_cf_boundary_2026-09-20.py not found in $(pwd)"
  exit 2
fi
if ! grep -q "def cf_coupled_bisect" "$PKG_PATH/plate_solver/piezo_monolithic.py"; then
  echo "FATAL: $PKG_PATH/plate_solver/piezo_monolithic.py has no cf_coupled_bisect -- push the updated file first"
  exit 2
fi
if ! grep -q "_piezo_mono_coupled_cf_root_worker" "$PKG_PATH/plate_solver/workers.py"; then
  echo "FATAL: $PKG_PATH/plate_solver/workers.py has no _piezo_mono_coupled_cf_root_worker -- push the updated file first"
  exit 2
fi

python3 -u probe_piezo_p5_cf_boundary_2026-09-20.py
echo "end=$(date)"
