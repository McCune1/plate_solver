#!/bin/bash
#SBATCH --job-name=p4_lever8_zp4
#SBATCH --output=ansys_queue_piezo_ff_zp4_%j.out
#SBATCH --error=ansys_queue_piezo_ff_zp4_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=00:45:00
#SBATCH --mem=8G
#SBATCH --requeue
#SBATCH --partition=requeue
#SBATCH --chdir=/home/ghmkfh/PythonMill/Plate_Solver_Package/Ansys/NewAnsys

# Paper 4 (piezo) Lever 8 -- six decks, NDIV_ZP 2 -> 4 clones of the
# confirmed 1/12 OC (job 2491761) and scboth (job 2491831) pair at
# h1/2h=1/12, 1/8, 1/5. Mesh refinement only.
#
# PUSH the six .inp files, this script, and
#   ansys_queue_manifest_piezo_ff_zp4_2026-09-16.txt
# to Ansys/NewAnsys/, then:
#   sbatch submit_ansys_queue_piezo_ff_zp4_2026-09-16.sh
#
# After a clean run, read each deck's own _modes.txt /
# _discriminator.txt / _out.txt. Use NUMBER OF ERROR MESSAGES
# ENCOUNTERED and RUN COMPLETED. Do not trust the queue script's
# own PASS/FAIL label.

echo "Job $SLURM_JOB_ID on $(hostname)"
echo "Start: $(date)"
if [ ! -f run_ansys_queue.sh ]; then
  echo "FATAL: must be Ansys/NewAnsys/"
  exit 2
fi
module load ansys/2024R2 2>/dev/null || module load ansys 2>/dev/null
export CPUS=${SLURM_CPUS_PER_TASK:-8}
export MANIFEST=ansys_queue_manifest_piezo_ff_zp4_2026-09-16.txt
export STOP_ON_FAIL=0
bash run_ansys_queue.sh
echo "End: $(date)"
