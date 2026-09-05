#!/bin/bash
#SBATCH --job-name=ip_1115_mac
#SBATCH --output=probe_ip_ffp1_1115_pair_mac_v1_%j.out
#SBATCH --error=probe_ip_ffp1_1115_pair_mac_v1_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=02:30:00
#SBATCH --partition=requeue
#SBATCH --requeue

# MAC on Table 7 1114.9/1115.3 Hz isolates (job 2428351).
# Submit AFTER submit_ansys_ip_ffp1_mac_fullblock.sh:
#   sbatch --dependency=afterok:$ANSYS_JOB submit_ip_ffp1_1115_pair_mac_v1.sh
#
# Copy onto Mill (package root):
#   probe_ip_ffp1_close_pair_mac_v1.py
#   submit_ip_ffp1_1115_pair_mac_v1.sh
# Already on Mill: ip_mac_bar.py, probe_ip_mac_adversarial_fullblock_v1.py
# Do NOT copy plate_solver/. Do NOT bump SOLVER_VERSION.
# Cost: 2 gates + 3 Omegas, ~5 min each ~= 25 min.

export PKG_PATH=/home/ghmkfh/PythonMill/Plate_Solver_Package
export EXPECT_SOLVER_VERSION=2026-07-10.s10
export DPS=40
export RULE_DELTA=0.001
export R40_BC_KIND=free_free
export PAIR=1115
export FB_PREFIX=ip_mac_ffp1_nu030

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMBA_NUM_THREADS=1

module load python/3.12.1
source /home/ghmkfh/PythonMill/Plate_Solver_Package/venv/bin/activate
cd /home/ghmkfh/PythonMill/Plate_Solver_Package/ || exit 2

if [ ! -f probe_ip_ffp1_close_pair_mac_v1.py ]; then
  echo "FATAL: probe_ip_ffp1_close_pair_mac_v1.py not in $(pwd)"
  exit 2
fi
if [ ! -f ip_mac_bar.py ]; then
  echo "FATAL: ip_mac_bar.py missing."
  exit 2
fi
N=$(ls -1 Ansys/NewAnsys/${FB_PREFIX}_m*_mesh96.txt 2>/dev/null | wc -l)
echo "dumps for ${FB_PREFIX}: $N/57"
if [ "$N" -lt 57 ]; then
  echo "FATAL: short dump. ANSYS job must finish first (afterok)."
  exit 3
fi

echo "Job $SLURM_JOB_ID on $(hostname)"
echo "PAIR=$PAIR PREFIX=$FB_PREFIX"
echo "Start: $(date)"
python3 -u probe_ip_ffp1_close_pair_mac_v1.py
rc=$?
echo "probe exit $rc"
echo "End: $(date)"
exit $rc
