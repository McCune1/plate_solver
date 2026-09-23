#!/bin/bash
#SBATCH --job-name=p4_lever7_oc_sc
#SBATCH --output=ansys_queue_piezo_ff_oc_scboth_ratios_%j.out
#SBATCH --error=ansys_queue_piezo_ff_oc_scboth_ratios_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=00:30:00
#SBATCH --mem=8G
#SBATCH --requeue
#SBATCH --partition=requeue
#SBATCH --chdir=/home/ghmkfh/PythonMill/Plate_Solver_Package/Ansys/NewAnsys

# Paper 4 (piezo) Lever 7 -- four decks, clones of the confirmed 1/12
# OC (job 2491761) and 1/12 scboth (job 2491831) decks with only H1
# and the closed-form target changed (LESSONS_LEARNED.md Sec 18.178,
# PAPER4_LEVERS_2026-09-16.md item 7).
#
# PUSH these files to
#   /home/ghmkfh/PythonMill/Plate_Solver_Package/Ansys/NewAnsys/
# then from that directory:
#   sbatch submit_ansys_queue_piezo_ff_oc_scboth_ratios_2026-09-16.sh
#
# Files:
#   ansys_p4_epsrel_ff_lanb_h18_oc_2026-09-16.inp
#   ansys_p4_epsrel_ff_lanb_h18_scboth_2026-09-16.inp
#   ansys_p4_epsrel_ff_lanb_h15_oc_2026-09-16.inp
#   ansys_p4_epsrel_ff_lanb_h15_scboth_2026-09-16.inp
#   ansys_queue_manifest_piezo_ff_oc_scboth_ratios_2026-09-16.txt
#   submit_ansys_queue_piezo_ff_oc_scboth_ratios_2026-09-16.sh
#   build_ff_oc_scboth_ratio_decks_2026-09-16.py   (provenance only; not needed to run)
#
# After a clean run, read each deck's own
#   ansys_p4_epsrel_ff_lanb_<tag>_modes.txt
#   ansys_p4_epsrel_ff_lanb_<tag>_discriminator.txt
#   ansys_p4_epsrel_ff_lanb_<tag>_2026-09-16_out.txt
# and grep the _out.txt for "solve:" and "*** ERROR ***".
# Pre-registered criteria are in the manifest.

echo "Job $SLURM_JOB_ID on $(hostname)"
echo "Start: $(date)"
if [ ! -f run_ansys_queue.sh ]; then
  echo "FATAL: must be Ansys/NewAnsys/"
  exit 2
fi
module load ansys/2024R2 2>/dev/null || module load ansys 2>/dev/null
export CPUS=${SLURM_CPUS_PER_TASK:-8}
export MANIFEST=ansys_queue_manifest_piezo_ff_oc_scboth_ratios_2026-09-16.txt
export STOP_ON_FAIL=0
bash run_ansys_queue.sh
echo "End: $(date)"
