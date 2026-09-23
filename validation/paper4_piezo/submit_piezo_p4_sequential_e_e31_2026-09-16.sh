#!/bin/bash
#SBATCH --job-name=p4_seq_e_e31
#SBATCH --output=piezo_p4_sequential_e_e31_2026-09-16_%j.out
#SBATCH --error=piezo_p4_sequential_e_e31_2026-09-16_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:30:00
#SBATCH --mem=4G
#SBATCH --partition=requeue
#SBATCH --requeue

# Paper 4 Lever 15: the inverse chapter's own prescribed experiment, run
# SEQUENTIALLY for the first time -- E from C-C SC (job 2491531's
# operator), then e31 from F-F n=0 OC holding E_hat (not E_true), then
# compound uncertainty propagation and a negative control on C-C SC/OC
# as the e31 observable. See probe_piezo_p4_sequential_e_e31_2026-09-16.py
# docstring and LESSONS_LEARNED.md Sec 18.182.
#
# SANDBOX RESULT this session (dps=100, full production precision, not a
# reduced smoke test): SENTINEL PASS_ALL. This cluster run is archival
# confirmation through the real deploy + worker path, same discipline as
# job 2491764 confirming Sec 18.175's own in-sandbox PASS_ALL. See the
# probe's own docstring "SANDBOX RESULT" paragraph for the full numbers.
#
# Fully serial (step 2 needs step 1's own E_hat, step 4 needs it too) --
# no ProcessPoolExecutor, unlike the sibling inverse-recovery probe's
# 2-worker E/e31 split (those two tasks were independent; these are not).
# 1 CPU is enough. 30 min is generous headroom (sandbox total was under
# 5 minutes end to end at these settings).

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
export PIEZO_DPS=100
export PIEZO_SC_INNER_ITERS=40
export PIEZO_OC_INNER_ITERS=50
export PIEZO_OUTER_ITERS=26

if [ ! -f probe_piezo_p4_sequential_e_e31_2026-09-16.py ]; then
  echo "FATAL: probe_piezo_p4_sequential_e_e31_2026-09-16.py not found in $(pwd)"
  exit 2
fi
if ! grep -q "def oc_ff_det" "$PKG_PATH/plate_solver/piezo_solver.py"; then
  echo "FATAL: $PKG_PATH/plate_solver/piezo_solver.py has no oc_ff_det"
  echo "  (job 2491661's exact failure mode: probe was on the mill,"
  echo "  package methods were not). Push plate_solver/piezo_solver.py"
  echo "  first."
  exit 2
fi
if ! grep -q "def cc_coupled_det" "$PKG_PATH/plate_solver/piezo_solver.py"; then
  echo "FATAL: $PKG_PATH/plate_solver/piezo_solver.py has no cc_coupled_det"
  exit 2
fi

python3 -u probe_piezo_p4_sequential_e_e31_2026-09-16.py
echo "end=$(date)"
