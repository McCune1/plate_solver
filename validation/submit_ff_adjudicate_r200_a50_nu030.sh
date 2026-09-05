#!/bin/bash
#SBATCH --job-name=ff_adj_r200a50
#SBATCH --output=ff_adjudicate_r200_a50_nu030_%j.out
#SBATCH --error=ff_adjudicate_r200_a50_nu030_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=64
#SBATCH --time=03:00:00
#SBATCH --requeue
#SBATCH --partition=requeue

# Rank 16 adjudicate. Run ONLY after:
#   1. both r200_a50_nu030_capture_{oop,ip}.json have complete=true
#   2. python3 merge_r200_a50_nu030_candidates.py succeeded
#   3. Ansys/NewAnsys/geomsweep_{oop,ip}_r200_a50_nu030.txt exist
#      and have the expected row counts (40 OOP / 80 IP)
#
# 64 cpus matches the nu=0.35 r150 adjudicate (cluster core ceiling).
# a50 should have fewer candidates than r200/a100's 108.

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export BLIS_NUM_THREADS=1

cd /home/ghmkfh/PythonMill/Plate_Solver_Package

export PKG_PATH=/home/ghmkfh/PythonMill/Plate_Solver_Package
export DPS=40
export EXPECT_SOLVER_VERSION=2026-07-10.s10
export CANDIDATES_JSON=/home/ghmkfh/PythonMill/Plate_Solver_Package/r200_a50_nu030_candidates.json
export ANSYS_DIR=/home/ghmkfh/PythonMill/Plate_Solver_Package/Ansys/NewAnsys
export FE_MATCH_TOL_PCT=3.0
export N_WORKERS=${SLURM_CPUS_PER_TASK:-64}

module load python/3.12.1
source /home/ghmkfh/PythonMill/Plate_Solver_Package/venv/bin/activate

if [ ! -f probe_ff_adjudicate_candidates.py ]; then
  echo "FATAL: probe_ff_adjudicate_candidates.py not found in $(pwd)."
  exit 2
fi
if [ ! -f r200_a50_nu030_candidates.json ]; then
  echo "FATAL: r200_a50_nu030_candidates.json missing. Merge first."
  exit 2
fi
if [ ! -f "$ANSYS_DIR/geomsweep_oop_r200_a50_nu030.txt" ] || \
   [ ! -f "$ANSYS_DIR/geomsweep_ip_r200_a50_nu030.txt" ]; then
  echo "FATAL: ANSYS txt files missing in $ANSYS_DIR"
  exit 2
fi

echo "Job $SLURM_JOB_ID on $(hostname): N_WORKERS=$N_WORKERS"
echo "Start: $(date)"
python3 -u probe_ff_adjudicate_candidates.py
echo "End: $(date)"
