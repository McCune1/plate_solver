#!/bin/bash
#SBATCH --job-name=p5_mono_ff
#SBATCH --output=ansys_queue_piezo_p5_mono_ff_%j.out
#SBATCH --error=ansys_queue_piezo_p5_mono_ff_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=00:20:00
#SBATCH --mem=8G
#SBATCH --requeue
#SBATCH --partition=requeue
#SBATCH --chdir=/home/ghmkfh/PythonMill/Plate_Solver_Package/Ansys/NewAnsys

# Paper 5 Option A -- F-F monolithic PZT-4: e0, scboth, OC in that order.
# Clone of Paper 4 F-F epsrel LANB (jobs 2491761 / 2491831).
#
# PUSH these files to
#   /home/ghmkfh/PythonMill/Plate_Solver_Package/Ansys/NewAnsys/
# then from that directory:
#   sbatch submit_ansys_queue_piezo_p5_mono_ff_2026-09-19.sh
#
# Files:
#   ansys_p5_epsrel_ff_lanb_e0_2026-09-19.inp
#   ansys_p5_epsrel_ff_lanb_scboth_2026-09-19.inp
#   ansys_p5_epsrel_ff_lanb_oc_2026-09-19.inp
#   ansys_queue_manifest_piezo_p5_mono_ff_2026-09-19.txt
#   submit_ansys_queue_piezo_p5_mono_ff_2026-09-19.sh
#   build_p5_mono_ff_decks_2026-09-19.py   (provenance only; not needed to run)
#
# After a clean run, read (never the .QUEUE_OK sentinel alone):
#   ansys_p5_epsrel_ff_lanb_e0_modes.txt
#   ansys_p5_epsrel_ff_lanb_e0_discriminator.txt
#   ansys_p5_epsrel_ff_lanb_e0_2026-09-19_out.txt
#   ansys_p5_epsrel_ff_lanb_scboth_modes.txt
#   ansys_p5_epsrel_ff_lanb_scboth_discriminator.txt
#   ansys_p5_epsrel_ff_lanb_scboth_2026-09-19_out.txt
#   ansys_p5_epsrel_ff_lanb_oc_modes.txt
#   ansys_p5_epsrel_ff_lanb_oc_discriminator.txt
#   ansys_p5_epsrel_ff_lanb_oc_2026-09-19_out.txt
# Judge MAPDL from NUMBER OF ERROR MESSAGES ENCOUNTERED and RUN
# COMPLETED. Pre-registered bars are in the manifest.
# Do not retune mesh/units/TBDATA/LANB.

echo "Job $SLURM_JOB_ID on $(hostname)"
echo "Start: $(date)"
if [ ! -f run_ansys_queue.sh ]; then
  echo "FATAL: must be Ansys/NewAnsys/"
  exit 2
fi
module load ansys/2024R2 2>/dev/null || module load ansys 2>/dev/null
export CPUS=${SLURM_CPUS_PER_TASK:-8}
export MANIFEST=ansys_queue_manifest_piezo_p5_mono_ff_2026-09-19.txt
export STOP_ON_FAIL=0
bash run_ansys_queue.sh
echo "End: $(date)"
