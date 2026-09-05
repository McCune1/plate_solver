#!/bin/bash
#SBATCH --job-name=ansys_ff_mesh80
#SBATCH --output=ansys_ff_isoval_mesh80_%j.out
#SBATCH --error=ansys_ff_isoval_mesh80_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=02:00:00
#SBATCH --partition=general

# Mesh-convergence check, NRAD=NCIRC=80 (up from the baseline 48), TARGET
# case only (PHI=90). LESSONS §23.7 item 2, follow-up 2 -- the finer of the
# two mesh-convergence points the original deck's comment named
# ("NRAD=NCIRC=64,80"). Includes the same Uz-dominance tagging. Walltime
# bumped further than mesh64's since 80x80 has ~2.8x the baseline's DOFs;
# adjust if your queue needs more -- this is the most expensive of the
# three new jobs.
module load ansys/2024R2 2>/dev/null || module load ansys 2>/dev/null

mapdl -b -smp -np ${SLURM_CPUS_PER_TASK:-8} -i ansys_ff_isoval_mesh80.inp -o ansys_ff_isoval_mesh80_out.txt

echo "=== TARGET pass, mesh 80x80 ==="
cat sector_modes_target_mesh80.txt
echo
echo "=== TARGET pass, mesh 80x80 -- Uz-dominance tags ==="
cat sector_modes_target_mesh80_tags.txt
