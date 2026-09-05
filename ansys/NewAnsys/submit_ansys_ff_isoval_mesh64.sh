#!/bin/bash
#SBATCH --job-name=ansys_ff_mesh64
#SBATCH --output=ansys_ff_isoval_mesh64_%j.out
#SBATCH --error=ansys_ff_isoval_mesh64_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=01:15:00
#SBATCH --partition=general

# Mesh-convergence check, NRAD=NCIRC=64 (up from the baseline 48), TARGET
# case only (PHI=90). LESSONS §23.7 item 2, follow-up 2 -- the original
# deck's own trailing comment flagged this as unrun. Includes the same
# Uz-dominance tagging as the tagged-baseline run. Walltime bumped from the
# baseline's 45 min since the finer mesh has ~1.8x the DOFs; adjust if your
# queue needs more.
module load ansys/2024R2 2>/dev/null || module load ansys 2>/dev/null

mapdl -b -smp -np ${SLURM_CPUS_PER_TASK:-8} -i ansys_ff_isoval_mesh64.inp -o ansys_ff_isoval_mesh64_out.txt

echo "=== TARGET pass, mesh 64x64 ==="
cat sector_modes_target_mesh64.txt
echo
echo "=== TARGET pass, mesh 64x64 -- Uz-dominance tags ==="
cat sector_modes_target_mesh64_tags.txt
