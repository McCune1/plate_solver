#!/bin/bash
#SBATCH --job-name=ansys_hs_parity
#SBATCH --output=ansys_hs_parity_%j.out
#SBATCH --error=ansys_hs_parity_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=01:30:00
#SBATCH --partition=general

# Half-sector SYMM/ANTI parity classification of the FF-P1 free-free FE
# modes (addendum 2026-07-11 Sec I's "genuinely new lever" for the OOP
# ceiling item; LESSONS_LEARNED Sec 10.7 item 8). Four modal solves
# (SYMM/ANTI x mesh64/mesh96-equivalent) on a half sector, each smaller
# than the full-sector NMODES=60 runs that fit in 45 min -- 1.5h is ample.
#
# Independent of (can run in parallel with):
#   ../submit_ip_baseline_capture.sh
#   ../submit_ip_pairmember_finescan.sh
#
# See ansys_ff_oop_halfsector.inp header for the pre-registered
# convergence gate, calibration gate (four known-parity modes MUST match
# Sec 10.6-H before the 111.46 read), and the decisive check.
module load ansys/2024R2 2>/dev/null || module load ansys 2>/dev/null

mapdl -b -smp -np ${SLURM_CPUS_PER_TASK:-8} -i ansys_ff_oop_halfsector.inp -o ansys_ff_oop_halfsector_out.txt

echo "=== SYMM (EVEN) mesh64-equiv ==="
cat oop_hs_symm_mesh64.txt
echo
echo "=== SYMM (EVEN) mesh96-equiv ==="
cat oop_hs_symm_mesh96.txt
echo
echo "=== ANTI (ODD) mesh64-equiv ==="
cat oop_hs_anti_mesh64.txt
echo
echo "=== ANTI (ODD) mesh96-equiv ==="
cat oop_hs_anti_mesh96.txt
echo
echo "Calibration gate (Sec 10.6-H): 15.128->SYMM, 24.190->ANTI,"
echo "37.818->ANTI, 53.828->SYMM. Then: which pass holds Omega_lit~111.46?"
