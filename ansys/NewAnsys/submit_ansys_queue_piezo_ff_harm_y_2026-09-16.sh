#!/bin/bash
#SBATCH --job-name=p4_lever9_harmy
#SBATCH --output=ansys_queue_piezo_ff_harm_y_%j.out
#SBATCH --error=ansys_queue_piezo_ff_harm_y_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --requeue
#SBATCH --partition=requeue
#SBATCH --chdir=/home/ghmkfh/PythonMill/Plate_Solver_Package/Ansys/NewAnsys

# Paper 4 (piezo) Lever 9 -- one harmonic driven-admittance deck at
# h1/2h=1/12. New family (ANTYPE,HARMIC), not a modal clone.
#
# PUSH
#   ansys_p4_harm_ff_y_h112_2026-09-16.inp
#   ansys_queue_manifest_piezo_ff_harm_y_2026-09-16.txt
#   this script
# then:
#   sbatch submit_ansys_queue_piezo_ff_harm_y_2026-09-16.sh
#
# After a clean run, read ansys_p4_harm_ff_y_h112_sweep.csv and the
# _out.txt. Use NUMBER OF ERROR MESSAGES ENCOUNTERED and RUN
# COMPLETED. Do not trust the queue script's own PASS/FAIL label.

echo "Job $SLURM_JOB_ID on $(hostname)"
echo "Start: $(date)"
if [ ! -f run_ansys_queue.sh ]; then
  echo "FATAL: must be Ansys/NewAnsys/"
  exit 2
fi
module load ansys/2024R2 2>/dev/null || module load ansys 2>/dev/null
export CPUS=${SLURM_CPUS_PER_TASK:-8}
export MANIFEST=ansys_queue_manifest_piezo_ff_harm_y_2026-09-16.txt
export STOP_ON_FAIL=0
bash run_ansys_queue.sh
echo "End: $(date)"
