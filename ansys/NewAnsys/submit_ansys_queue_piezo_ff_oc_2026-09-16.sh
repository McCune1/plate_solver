#!/bin/bash
#SBATCH --job-name=p4_epsrel_ffoc
#SBATCH --output=ansys_queue_piezo_ff_oc_%j.out
#SBATCH --error=ansys_queue_piezo_ff_oc_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=00:15:00
#SBATCH --mem=8G
#SBATCH --requeue
#SBATCH --partition=requeue
#SBATCH --chdir=/home/ghmkfh/PythonMill/Plate_Solver_Package/Ansys/NewAnsys

# Paper 4 (piezo) F-F OPEN-CIRCUIT -- one deck, clone of job 2490706
# h1/2h=1/12 SC LANB with only the electrical BC changed
# (LESSONS_LEARNED.md Sec 18.175).
#
# PUSH these files to
#   /home/ghmkfh/PythonMill/Plate_Solver_Package/Ansys/NewAnsys/
# then from that directory:
#   sbatch submit_ansys_queue_piezo_ff_oc_2026-09-16.sh
#
# Files:
#   ansys_p4_epsrel_ff_lanb_oc_2026-09-16.inp
#   ansys_queue_manifest_piezo_ff_oc_2026-09-16.txt
#   submit_ansys_queue_piezo_ff_oc_2026-09-16.sh
#   build_ff_oc_deck_2026-09-16.py   (provenance only; not needed to run)
#
# After a clean run, read
#   ansys_p4_epsrel_ff_lanb_oc_modes.txt
#   ansys_p4_epsrel_ff_lanb_oc_discriminator.txt
#   ansys_p4_epsrel_ff_lanb_oc_2026-09-16_out.txt
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
export MANIFEST=ansys_queue_manifest_piezo_ff_oc_2026-09-16.txt
export STOP_ON_FAIL=0
bash run_ansys_queue.sh
echo "End: $(date)"
