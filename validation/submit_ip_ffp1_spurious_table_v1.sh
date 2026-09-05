#!/bin/bash
#SBATCH --job-name=ip_spurious
#SBATCH --output=probe_ip_ffp1_spurious_table_v1_%j.out
#SBATCH --error=probe_ip_ffp1_spurious_table_v1_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --partition=requeue
#SBATCH --requeue

# IP analogue of Table 6: residual-screen the 8 FF-P1 MAC calibration
# candidates (4 Table-7 REALs + 4 catalogued artifacts). No ANSYS dump.
#
# Copy onto Mill (package root):
#   probe_ip_ffp1_spurious_table_v1.py
#   submit_ip_ffp1_spurious_table_v1.sh
# Already on Mill: ip_mac_bar.py, probe_ip_mac_adversarial_fullblock_v1.py
# Do NOT copy plate_solver/. Do NOT bump SOLVER_VERSION.
# Cost: 8 reconstructions ~= 40 min.

export PKG_PATH=/home/ghmkfh/PythonMill/Plate_Solver_Package
export EXPECT_SOLVER_VERSION=2026-07-10.s10
export DPS=40
export R40_BC_KIND=free_free
export MAC_N_DOFS=20

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMBA_NUM_THREADS=1

module load python/3.12.1
source /home/ghmkfh/PythonMill/Plate_Solver_Package/venv/bin/activate
cd /home/ghmkfh/PythonMill/Plate_Solver_Package/ || exit 2

if [ ! -f probe_ip_ffp1_spurious_table_v1.py ]; then
  echo "FATAL: probe_ip_ffp1_spurious_table_v1.py not in $(pwd)"
  exit 2
fi

echo "Job $SLURM_JOB_ID on $(hostname)"
echo "Start: $(date)"
python3 -u probe_ip_ffp1_spurious_table_v1.py
rc=$?
echo "probe exit $rc"
echo "End: $(date)"
exit $rc
