#!/bin/bash
#SBATCH --job-name=ip_pair_bias
#SBATCH --output=probe_ip_pairing_fallback_bias_v1_%j.out
#SBATCH --error=probe_ip_pairing_fallback_bias_v1_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=02:30:00
#SBATCH --partition=requeue
#SBATCH --requeue

# Cartesian-midside pairing fallback vs exact polar pairing, at the
# geometry where exact pairing works (r200/a100, PHI=180). Uses the
# existing Rank 18 dump -- no new ANSYS job.
#
# Copy onto Mill (package root):
#   probe_ip_pairing_fallback_bias_v1.py
#   submit_ip_pairing_fallback_bias_v1.sh
# Already on Mill:
#   ip_mac_bar.py
#   probe_ip_mac_adversarial_fullblock_v1.py
#   Ansys/NewAnsys/ip_mac_fullblock_r200a100_m{4..60}_mesh96.txt
# Do NOT copy plate_solver/. Do NOT bump SOLVER_VERSION.
# Cost: 8 candidates x 2 pairing tables ~= 50 min.

export PKG_PATH=/home/ghmkfh/PythonMill/Plate_Solver_Package
export EXPECT_SOLVER_VERSION=2026-07-10.s10
export DPS=40
export RULE_DELTA=0.001
export R40_BC_KIND=free_free

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMBA_NUM_THREADS=1

module load python/3.12.1
source /home/ghmkfh/PythonMill/Plate_Solver_Package/venv/bin/activate
cd /home/ghmkfh/PythonMill/Plate_Solver_Package/ || exit 2

if [ ! -f probe_ip_pairing_fallback_bias_v1.py ]; then
  echo "FATAL: probe_ip_pairing_fallback_bias_v1.py not in $(pwd)"
  exit 2
fi
N=$(ls -1 Ansys/NewAnsys/ip_mac_fullblock_r200a100_m*_mesh96.txt 2>/dev/null | wc -l)
echo "Rank 18 dumps: $N/57"
if [ "$N" -lt 57 ]; then
  echo "FATAL: Rank 18 full-block dump missing."
  exit 3
fi

echo "Job $SLURM_JOB_ID on $(hostname)"
echo "Start: $(date)"
python3 -u probe_ip_pairing_fallback_bias_v1.py
rc=$?
echo "probe exit $rc"
echo "End: $(date)"
exit $rc
