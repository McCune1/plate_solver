#!/bin/bash
#SBATCH --job-name=p4_piezo_cc_gate
#SBATCH --output=piezo_p4_cc_table4_gate_2026-09-15_%j.out
#SBATCH --error=piezo_p4_cc_table4_gate_2026-09-15_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=27
#SBATCH --time=00:20:00
#SBATCH --mem=8G
#SBATCH --partition=requeue
#SBATCH --requeue

# Paper 4 (piezoelectric ring), C-C forward model validation gate.
# LESSONS_LEARNED.md Sec 18.147: cubic_coeffs() had every A2*omega^2
# (inertial) term scaled by h1**2 instead of h1**4 -- found and fixed
# 2026-09-15, confirmed on 6 hand-picked Duan2005 Table 4 points in the
# sandbox at dps=60-300 (all <0.05% rel. err). This job reproduces ALL 27
# non-trivial CPT/r0-h=60 Table 4 rows (p=0,1,2 x n=0,1,2 x h1/2h in
# {1/12,1/8,1/5}) at production precision as the real validation pass
# before this model is trusted for the paper or extended to F-F/C-F.
#
# Standalone script -- does NOT import plate_solver (pure mpmath), so no
# PKG_PATH/EXPECT_SOLVER_VERSION needed. Sandbox timing at dps=60-100
# (extraprec=200-300) was 1-25s/point/iteration-batch depending on mode
# (worst case p=2/n=2/h1=1/5 ~22s for 30 iters at dps=60) -- 27 points on
# 27 workers should finish in a few minutes; 20 min budget is generous
# headroom, not a tight estimate.
#
# PRE-REGISTERED PASS/FAIL (see the .py's own docstring): all 27 points
# within 0.15% of their Table 4 target AND a genuine bracket sign flip on
# each -- read piezo_p4_cc_table4_gate_results.json's SENTINEL line
# (PASS_ALL vs FAIL_ALL) plainly; a single miss is a real finding to
# report back, not something to average away.
#
# Submit from /home/ghmkfh/PythonMill/Plate_Solver_Package/Paper4_Piezo
# (the venv is one level up and still sourced from there -- only the cwd
# for running the .py and writing its outputs moves into this subfolder).

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

export N_WORKERS=27
export PIEZO_DPS=100
export PIEZO_EXTRAPREC=300
export PIEZO_BISECT_ITERS=45

if [ ! -f probe_piezo_p4_cc_table4_gate_2026-09-15.py ]; then
  echo "FATAL: probe_piezo_p4_cc_table4_gate_2026-09-15.py not found in $(pwd)"
  exit 2
fi

python3 -u probe_piezo_p4_cc_table4_gate_2026-09-15.py
echo "end=$(date)"
