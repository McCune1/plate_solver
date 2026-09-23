#!/bin/bash
#SBATCH --job-name=p4_lever9_yv2
#SBATCH --output=ansys_queue_piezo_ff_harm_y_v2_%j.out
#SBATCH --error=ansys_queue_piezo_ff_harm_y_v2_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --requeue
#SBATCH --partition=requeue
#SBATCH --chdir=/home/ghmkfh/PythonMill/Plate_Solver_Package/Ansys/NewAnsys

# Lever 9 v2 -- CHRG dump fix of job 2494669. Same mesh/HARFRQ.
#   sbatch submit_ansys_queue_piezo_ff_harm_y_v2_2026-09-16.sh
# After: read ansys_p4_harm_ff_y_v2_h112_sweep.csv and _out.txt
# (PRRSOL,CHRG, NUMBER OF ERROR MESSAGES ENCOUNTERED, RUN COMPLETED).

echo "Job $SLURM_JOB_ID on $(hostname)"
echo "Start: $(date)"
if [ ! -f run_ansys_queue.sh ]; then
  echo "FATAL: must be Ansys/NewAnsys/"
  exit 2
fi
module load ansys/2024R2 2>/dev/null || module load ansys 2>/dev/null
export CPUS=${SLURM_CPUS_PER_TASK:-8}
export MANIFEST=ansys_queue_manifest_piezo_ff_harm_y_v2_2026-09-16.txt
export STOP_ON_FAIL=0
bash run_ansys_queue.sh
echo "End: $(date)"
