#!/bin/bash
#SBATCH --job-name=p4_epsrel_scboth
#SBATCH --output=ansys_queue_piezo_ff_scboth_%j.out
#SBATCH --error=ansys_queue_piezo_ff_scboth_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=00:15:00
#SBATCH --mem=8G
#SBATCH --requeue
#SBATCH --partition=requeue
#SBATCH --chdir=/home/ghmkfh/PythonMill/Plate_Solver_Package/Ansys/NewAnsys

# F-F both-faces SC control for job 2491761's OC FE (Sec 18.177).
# PUSH ansys_p4_epsrel_ff_lanb_scboth_2026-09-16.inp, this submit, and
# ansys_queue_manifest_piezo_ff_scboth_2026-09-16.txt, then:
#   sbatch submit_ansys_queue_piezo_ff_scboth_2026-09-16.sh
# Read _modes.txt / _discriminator.txt / _out.txt, never .QUEUE_OK alone.

echo "Job $SLURM_JOB_ID on $(hostname)"
echo "Start: $(date)"
if [ ! -f run_ansys_queue.sh ]; then
  echo "FATAL: must be Ansys/NewAnsys/"
  exit 2
fi
module load ansys/2024R2 2>/dev/null || module load ansys 2>/dev/null
export CPUS=${SLURM_CPUS_PER_TASK:-8}
export MANIFEST=ansys_queue_manifest_piezo_ff_scboth_2026-09-16.txt
export STOP_ON_FAIL=0
bash run_ansys_queue.sh
echo "End: $(date)"
