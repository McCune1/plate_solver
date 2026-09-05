#!/bin/bash
#SBATCH --job-name=ansys_ff_isoval
#SBATCH --output=ansys_ff_isoval_%j.out
#SBATCH --error=ansys_ff_isoval_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=00:45:00
#SBATCH --partition=general

# Two-pass headless Mechanical APDL modal run: calibration (PHI=120, pins the
# Omega convention) then the actual free-free isotropic target (PHI=90),
# matched geometry/material to plate_solver's submit_ff_isoval.sh
# (R40_FF_P1=1.5:0.5, MAT_NU=0.30). See ansys_ff_isoval.inp header for the
# full comparison plan (LESSONS_LEARNED.md §23.7 item 2 / §23.6).
#
# Adjust the module name + solver binary to your cluster's Ansys install
# (check: 'module avail ansys').
module load ansys/2024R2 2>/dev/null || module load ansys 2>/dev/null

mapdl -b -smp -np ${SLURM_CPUS_PER_TASK:-8} -i ansys_ff_isoval.inp -o ansys_ff_isoval_out.txt

echo "=== CALIBRATION pass (expect first elastic Omega ~10.19) ==="
cat sector_modes_calib.txt
echo
echo "=== TARGET pass (compare elastic modes to plate_solver's"
echo "    [12.6596, 15.1389, 24.3442, 33.3391, 37.9182, 39.5575, 54.0755, 63.2524]) ==="
cat sector_modes_target.txt
