#!/bin/bash
#SBATCH --job-name=p4_epsrel
#SBATCH --output=ansys_queue_piezo_epsrel_%j.out
#SBATCH --error=ansys_queue_piezo_epsrel_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=00:30:00
#SBATCH --mem=8G
#SBATCH --requeue
#SBATCH --partition=requeue
#SBATCH --chdir=/home/ghmkfh/PythonMill/Plate_Solver_Package/Ansys/NewAnsys

# Paper 4 (piezo) PLANE223 relative-permittivity diagnostic -- four decks
# sequentially in this one submission (one ANSYS seat).
#
# PUSH these files to
#   /home/ghmkfh/PythonMill/Plate_Solver_Package/Ansys/NewAnsys/
# then from that directory:
#   sbatch submit_ansys_queue_piezo_epsrel_2026-09-16.sh
#
# Files:
#   ansys_p4_epsrel_ff_e0_2026-09-16.inp
#   ansys_p4_epsrel_ff_lanb_2026-09-16.inp
#   ansys_p4_epsrel_ff_subsp_2026-09-16.inp
#   ansys_p4_epsrel_cc_lanb_2026-09-16.inp
#   ansys_queue_manifest_piezo_epsrel_2026-09-16.txt
#   submit_ansys_queue_piezo_epsrel_2026-09-16.sh
#   build_epsrel_decks_2026-09-16.py   (provenance only; not needed to run)
#
# After a clean run, for EACH deck read its own
#   ansys_p4_epsrel_<tag>_modes.txt
#   ansys_p4_epsrel_<tag>_discriminator.txt
#   ansys_p4_epsrel_<tag>_2026-09-16_out.txt
# and grep the _out.txt for "solve:" (ill-conditioning residual) and
# "*** ERROR ***". Pre-registered criteria are in the manifest.

echo "Job $SLURM_JOB_ID on $(hostname)"
echo "Start: $(date)"
if [ ! -f run_ansys_queue.sh ]; then
  echo "FATAL: must be Ansys/NewAnsys/"
  exit 2
fi
module load ansys/2024R2 2>/dev/null || module load ansys 2>/dev/null
export CPUS=${SLURM_CPUS_PER_TASK:-8}
export MANIFEST=ansys_queue_manifest_piezo_epsrel_2026-09-16.txt
export STOP_ON_FAIL=0
bash run_ansys_queue.sh
echo "End: $(date)"
