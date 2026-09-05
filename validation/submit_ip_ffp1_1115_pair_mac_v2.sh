#!/bin/bash
#SBATCH --job-name=ip_1115_v2
#SBATCH --output=probe_ip_ffp1_1115_pair_mac_v2_%j.out
#SBATCH --error=probe_ip_ffp1_1115_pair_mac_v2_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=02:30:00
#SBATCH --partition=requeue
#SBATCH --requeue

# v2 of 2434286. Dump already on Mill (ip_mac_ffp1_nu030, job 2434285).
# GATE_ART is residual+pick+CALL, not 8-target MAC vs full-block max.
#
# Copy onto Mill (package root):
#   probe_ip_ffp1_close_pair_mac_v2.py
#   submit_ip_ffp1_1115_pair_mac_v2.sh
# Already on Mill: ip_mac_bar.py, probe_ip_mac_adversarial_fullblock_v1.py,
#   Ansys/NewAnsys/ip_mac_ffp1_nu030_m{4..60}_mesh96.txt
# Do NOT copy plate_solver/. Do NOT bump SOLVER_VERSION.

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

if [ ! -f probe_ip_ffp1_close_pair_mac_v2.py ]; then
  echo "FATAL: probe_ip_ffp1_close_pair_mac_v2.py not in $(pwd)"
  exit 2
fi
N=$(ls -1 Ansys/NewAnsys/${FB_PREFIX}_m*_mesh96.txt 2>/dev/null | wc -l)
echo "dumps for ${FB_PREFIX}: $N/57"
if [ "$N" -lt 57 ]; then
  echo "FATAL: dump missing. Job 2434285 should already be on this node."
  exit 3
fi

echo "Job $SLURM_JOB_ID on $(hostname)"
echo "PAIR=$PAIR PREFIX=$FB_PREFIX"
echo "Start: $(date)"
python3 -u probe_ip_ffp1_close_pair_mac_v2.py
rc=$?
echo "probe exit $rc"
echo "End: $(date)"
exit $rc
