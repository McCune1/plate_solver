#!/bin/bash
#SBATCH --job-name=ansys_ff_tagged
#SBATCH --output=ansys_ff_isoval_tagged_%j.out
#SBATCH --error=ansys_ff_isoval_tagged_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=00:45:00
#SBATCH --partition=general

# Same baseline mesh (NRAD=NCIRC=48) as the already-reported TARGET run,
# but Pass 2 now also tags each mode by UzFrac = SUM(Uz^2)/SUM(Ux^2+Uy^2+Uz^2)
# to check whether FE modes 7-14 are really the pure-bending sequence the
# positional comparison assumed, or have membrane/coupled modes interleaved
# (LESSONS §23.7 item 2, follow-up 1). See ansys_ff_isoval_tagged.inp header.
#
# Adjust the module name + solver binary to your cluster's Ansys install
# (check: 'module avail ansys').
module load ansys/2024R2 2>/dev/null || module load ansys 2>/dev/null

mapdl -b -smp -np ${SLURM_CPUS_PER_TASK:-8} -i ansys_ff_isoval_tagged.inp -o ansys_ff_isoval_tagged_out.txt

echo "=== CALIBRATION pass (expect first elastic Omega ~10.19) ==="
cat sector_modes_calib.txt
echo
echo "=== TARGET pass (compare elastic modes to plate_solver's"
echo "    [12.6596, 15.1389, 24.3442, 33.3391, 37.9182, 39.5575, 54.0755, 63.2524]) ==="
cat sector_modes_target.txt
echo
echo "=== TARGET pass Uz-dominance tags (>0.9 = bending-dominant, <0.5 = in-plane/coupled) ==="
cat sector_modes_target_tags.txt
