#!/bin/bash
#SBATCH --job-name=tables_gate
#SBATCH --output=tables_gate_%j.out
#SBATCH --error=tables_gate_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=24:00:00
#SBATCH --partition=general

# FULL regression gates on the deployed s10 tree against the UPDATED pinned
# tables (test_validated_tables.py, updated 2026-07-11: FFP1_OOP_RAW now 13
# modes incl. 2.320741/2.602500; SHI_REAL_LIT now 8/8 incl. 67.2235).
# These heavy gates have NOT been run since those updates -- running them
# is itself a completion requirement for free-free validation (the gate
# class that caught the s7 regression, job 2316094).
#   GATE_FFP1: fully-default FF-P1 OOP discovery must reproduce all 13
#              pinned raw Omegas to |d|<0.003 (~1.5-3h on 16 cores).
#   GATE_SHI:  fully-default Shi orthotropic FFFF discovery must match all
#              8 pinned Omega_lit values <1% (~3-5h on 16 cores).
# Plus the always-run presence/version pins (seconds).
#
# Independent of (can run in parallel with):
#   submit_oop_ceiling_parity_adjudication.sh, submit_ip_pin_prep.sh
#
# Budget 16h (A.25: budget wide; both gates sequential in one job).
# Read: any FAILURE here is a REGRESSION -- suspect
# SIGMIN_GAP_REFINE_MAX_GAPS truncation FIRST (Sec 16), not detector logic.

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export BLIS_NUM_THREADS=1

export N_WORKERS=${SLURM_CPUS_PER_TASK:-16}
export DPS=40
export PKG_PATH=/home/ghmkfh/PythonMill/Plate_Solver_Package
export EXPECT_SOLVER_VERSION=2026-07-10.s10
export GATE_FFP1=1
export GATE_SHI=1
export GATE_FFP1_IP=1   # added 2026-07-12: 43-location IP scan pin (job 2317336)

module load python/3.12.1
source /home/ghmkfh/PythonMill/Plate_Solver_Package/venv/bin/activate

cd /home/ghmkfh/PythonMill/Plate_Solver_Package/

echo "Job $SLURM_JOB_ID on $(hostname); N_WORKERS=$N_WORKERS"
echo "Start: $(date)"

if [ -f validation/test_validated_tables.py ]; then
  MOD=validation.test_validated_tables
elif [ -f tests/test_validated_tables.py ]; then
  MOD=tests.test_validated_tables
elif [ -f test_validated_tables.py ]; then
  MOD=test_validated_tables
else
  echo "FATAL: test_validated_tables.py not found (validation/, tests/, or root)."
  exit 2
fi
echo "running module: $MOD"

python3 -u -m unittest "$MOD" -v

echo "End: $(date)"
