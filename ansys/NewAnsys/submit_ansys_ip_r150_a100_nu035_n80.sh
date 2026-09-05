#!/bin/bash
#SBATCH --job-name=ansys_ip_r150_a100_nu035_n80
#SBATCH --output=ansys_ip_r150_a100_nu035_n80_%j.out
#SBATCH --error=ansys_ip_r150_a100_nu035_n80_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=01:00:00
#SBATCH --partition=general

# Raised-NMODES rematch of the a100 pass of ansys_geomsweep_ip_r150.inp
# (35 -> 80, nu=0.35). Submit from the SAME directory as
# ansys_ip_r150_a100_nu035_n80.inp (Ansys/NewAnsys/ on the cluster).
# Does not overwrite the published 35-mode file geomsweep_ip_r150_a100.txt.
#
# -j gives this run its own scratch prefix so it does not collide with
# the leftover file.lock / file.db sitting in NewAnsys/ from earlier
# MAPDL jobs (that is what killed jobs 2421046 / 2421049).
module load ansys/2024R2 2>/dev/null || module load ansys 2>/dev/null

echo "Job $SLURM_JOB_ID on $(hostname)"
echo "Start: $(date)"

mapdl -b -smp -np ${SLURM_CPUS_PER_TASK:-8} -j r150ip35n80 -i ansys_ip_r150_a100_nu035_n80.inp -o ansys_ip_r150_a100_nu035_n80_out.txt

echo "End: $(date)"
echo "--- geomsweep_ip_r150_a100_nu035_n80.txt ---"
cat geomsweep_ip_r150_a100_nu035_n80.txt 2>/dev/null
echo
echo "Expect 80 data rows after the header; last frequency should be"
echo "comfortably above ~972 Hz (highest unmatched nu=0.35 solver candidate)."
