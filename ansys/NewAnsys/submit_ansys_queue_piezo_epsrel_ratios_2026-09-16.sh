#!/bin/bash
#SBATCH --job-name=p4_epsrel_ratio
#SBATCH --output=ansys_queue_piezo_epsrel_ratios_%j.out
#SBATCH --error=ansys_queue_piezo_epsrel_ratios_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=00:20:00
#SBATCH --mem=8G
#SBATCH --requeue
#SBATCH --partition=requeue
#SBATCH --chdir=/home/ghmkfh/PythonMill/Plate_Solver_Package/Ansys/NewAnsys

# Paper 4 (piezo) F-F epsrel LANB deck, h1/2h=1/8 and 1/5 -- two decks
# sequentially in this one submission (one ANSYS seat). Direct follow-up
# to job 2490706 (LESSONS_LEARNED.md Sec 18.165); both decks are clones
# of the already-cluster-confirmed h1/2h=1/12 deck, no new hypothesis.
#
# PUSH these files to
#   /home/ghmkfh/PythonMill/Plate_Solver_Package/Ansys/NewAnsys/
# then from that directory:
#   sbatch submit_ansys_queue_piezo_epsrel_ratios_2026-09-16.sh
#
# Files:
#   ansys_p4_epsrel_ff_lanb_h18_2026-09-16.inp
#   ansys_p4_epsrel_ff_lanb_h15_2026-09-16.inp
#   ansys_queue_manifest_piezo_epsrel_ratios_2026-09-16.txt
#   submit_ansys_queue_piezo_epsrel_ratios_2026-09-16.sh
#   build_epsrel_ratio_decks_2026-09-16.py   (provenance only; not needed to run)
#
# After a clean run, for EACH deck read its own
#   ansys_p4_epsrel_ff_lanb_<tag>_modes.txt
#   ansys_p4_epsrel_ff_lanb_<tag>_discriminator.txt
#   ansys_p4_epsrel_ff_lanb_<tag>_2026-09-16_out.txt
# and grep the _out.txt for "solve:" (ill-conditioning residual, should
# be ABSENT per job 2490706's precedent) and "*** ERROR ***".
# Pre-registered criteria are in the manifest.

echo "Job $SLURM_JOB_ID on $(hostname)"
echo "Start: $(date)"
if [ ! -f run_ansys_queue.sh ]; then
  echo "FATAL: must be Ansys/NewAnsys/"
  exit 2
fi
module load ansys/2024R2 2>/dev/null || module load ansys 2>/dev/null
export CPUS=${SLURM_CPUS_PER_TASK:-8}
export MANIFEST=ansys_queue_manifest_piezo_epsrel_ratios_2026-09-16.txt
export STOP_ON_FAIL=0
bash run_ansys_queue.sh
echo "End: $(date)"
