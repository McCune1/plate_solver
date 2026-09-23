#!/bin/bash
#SBATCH --job-name=p4_epsrel_cf
#SBATCH --output=ansys_queue_piezo_cf_oc_scboth_%j.out
#SBATCH --error=ansys_queue_piezo_cf_oc_scboth_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=00:15:00
#SBATCH --mem=8G
#SBATCH --requeue
#SBATCH --partition=requeue
#SBATCH --chdir=/home/ghmkfh/PythonMill/Plate_Solver_Package/Ansys/NewAnsys

# Optional C-F FE: clones of jobs 2491761/2491831, mechanical BC only.
# PUSH the two .inp files, this submit, and
#   ansys_queue_manifest_piezo_cf_oc_scboth_2026-09-16.txt
# then from Ansys/NewAnsys/:
#   sbatch submit_ansys_queue_piezo_cf_oc_scboth_2026-09-16.sh
#
# After a clean run, read (never the .QUEUE_OK sentinel alone):
#   ansys_p4_epsrel_cf_lanb_scboth_modes.txt
#   ansys_p4_epsrel_cf_lanb_scboth_discriminator.txt
#   ansys_p4_epsrel_cf_lanb_scboth_2026-09-16_out.txt
#   ansys_p4_epsrel_cf_lanb_oc_modes.txt
#   ansys_p4_epsrel_cf_lanb_oc_discriminator.txt
#   ansys_p4_epsrel_cf_lanb_oc_2026-09-16_out.txt
# Judge MAPDL from NUMBER OF ERROR MESSAGES ENCOUNTERED and RUN
# COMPLETED. Pre-registered bars are in the manifest.

echo "Job $SLURM_JOB_ID on $(hostname)"
echo "Start: $(date)"
if [ ! -f run_ansys_queue.sh ]; then
  echo "FATAL: must be Ansys/NewAnsys/"
  exit 2
fi
module load ansys/2024R2 2>/dev/null || module load ansys 2>/dev/null
export CPUS=${SLURM_CPUS_PER_TASK:-8}
export MANIFEST=ansys_queue_manifest_piezo_cf_oc_scboth_2026-09-16.txt
export STOP_ON_FAIL=0
bash run_ansys_queue.sh
echo "End: $(date)"
