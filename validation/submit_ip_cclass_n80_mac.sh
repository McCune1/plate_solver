#!/bin/bash
#SBATCH --job-name=ip_c80_mac
#SBATCH --output=probe_ip_cclass_n80_mac_v1_%j.out
#SBATCH --error=probe_ip_cclass_n80_mac_v1_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --partition=requeue
#SBATCH --requeue

# C-class n80 extract-completeness MAC bar -- ONE key (a125 or a150).
# 16-way candidate pool. plate_solver N_WORKERS=1 (no nested pools).
# Required: RATIO_TAG  ANGLE_TAG=a125|a150
# Submit AFTER the matching n80 ANSYS dump is QUEUE_OK.

export RATIO_TAG=${RATIO_TAG:?set RATIO_TAG=r150|r167|r200|r250}
export ANGLE_TAG=${ANGLE_TAG:?set ANGLE_TAG=a125|a150}
export GS_NU=${GS_NU:-0.35}
export PKG_PATH=/home/ghmkfh/PythonMill/Plate_Solver_Package
export EXPECT_SOLVER_VERSION=2026-07-10.s10
export DPS=40
export RULE_DELTA=0.001
export R40_BC_KIND=free_free
export FB_PREFIX="ip_mac_gs_${RATIO_TAG}_${ANGLE_TAG}_n80"
export FB_MODE_LO=4
export FB_MODE_HI=83
export N_WORKERS=1
export C80_WORKERS=${SLURM_CPUS_PER_TASK:-16}

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMBA_NUM_THREADS=1

module load python/3.12.1
source /home/ghmkfh/PythonMill/Plate_Solver_Package/venv/bin/activate
cd /home/ghmkfh/PythonMill/Plate_Solver_Package/ || exit 2

if [ ! -f probe_ip_cclass_n80_mac_v1.py ]; then
  echo "FATAL: probe_ip_cclass_n80_mac_v1.py not in $(pwd)"
  exit 2
fi
if [ ! -f ip_mac_bar.py ]; then
  echo "FATAL: ip_mac_bar.py missing."
  exit 2
fi

N=$(ls -1 Ansys/NewAnsys/${FB_PREFIX}_m*_mesh96.txt 2>/dev/null | wc -l)
echo "dumps for ${FB_PREFIX}: $N/80"
if [ "$N" -lt 80 ]; then
  echo "FATAL: short n80 dump ($N/80). Finish the ANSYS queue first."
  exit 2
fi
echo "C80_WORKERS=$C80_WORKERS  N_WORKERS=$N_WORKERS  SLURM_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK"

python3 -u probe_ip_cclass_n80_mac_v1.py
echo "End: $(date)"
