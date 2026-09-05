#!/bin/bash
#SBATCH --job-name=ansys_oop_r200_a50
#SBATCH --output=ansys_oop_r200_a50_%j.out
#SBATCH --error=ansys_oop_r200_a50_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=01:00:00
#SBATCH --requeue
#SBATCH --partition=requeue

# Rank 16 FE. Submit from Ansys/NewAnsys/.
# -j r200a50oop avoids leftover NewAnsys/file.lock (jobs 2421046/2421049).
module load ansys/2024R2 2>/dev/null || module load ansys 2>/dev/null

echo "Job $SLURM_JOB_ID on $(hostname)"
echo "Start: $(date)"

mapdl -b -smp -np ${SLURM_CPUS_PER_TASK:-8} -j r200a50oop -i ansys_oop_r200_a50_nu030.inp -o ansys_oop_r200_a50_nu030_out.txt

echo "End: $(date)"
for f in geomsweep_oop_r200_a50_nu030.txt; do
  echo "--- $f ---"; cat "$f" 2>/dev/null; echo
done
