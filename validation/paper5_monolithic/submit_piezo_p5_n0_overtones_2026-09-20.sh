#!/bin/bash
#SBATCH --job-name=p5_n0_ot
#SBATCH --output=piezo_p5_n0_overtones_2026-09-20_%j.out
#SBATCH --error=piezo_p5_n0_overtones_2026-09-20_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:20:00
#SBATCH --mem=4G
#SBATCH --partition=requeue
#SBATCH --requeue

# Paper 5 roadmap Sec 8 item 6 optional-cheap: n=0 overtones. Pure
# probe-level scan on PiezoMonolithicOutOfPlaneSolver -- no new solver
# code, no SOLVER_VERSION bump. See
# probe_piezo_p5_n0_overtones_2026-09-20.py's own docstring for the
# full G_found/G_theorem/G_reflection/G2 pre-registration.
#
# In-sandbox timing this session (dps=35, inner_iters=35,
# ff_scan_n=350, cc_scan_n=300): ~71s end to end, SENTINEL PASS_ALL,
# F-F 3 pairs found, C-C 2 pairs found, stiffening theorem and e31
# reflection-twin identity (Sec 18.204) both hold at every pair.
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
export PIEZO_FF_SCAN_N=600
export PIEZO_CC_SCAN_N=500

if [ ! -f probe_piezo_p5_n0_overtones_2026-09-20.py ]; then
  echo "FATAL: probe_piezo_p5_n0_overtones_2026-09-20.py not found in $(pwd)"
  exit 2
fi
if ! grep -q "def coupled_bisect" "$PKG_PATH/plate_solver/piezo_monolithic.py"; then
  echo "FATAL: $PKG_PATH/plate_solver/piezo_monolithic.py has no coupled_bisect"
  exit 2
fi

python3 -u probe_piezo_p5_n0_overtones_2026-09-20.py
echo "end=$(date)"
