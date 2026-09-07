#!/bin/bash
#SBATCH --job-name=ansys_rect_ff_oop_mac
#SBATCH --output=ansys_rect_ff_oop_mac_%j.out
#SBATCH --error=ansys_rect_ff_oop_mac_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=04:00:00
#SBATCH --mem=16G
#SBATCH --requeue
#SBATCH --partition=requeue
#
# Isolated wrapper around run_ansys_queue.sh for
# ansys_queue_manifest_rect_ff_oop_mac.txt only -- five SHELL281 half-model
# FFFF OOP (flexural) MAC-extraction decks (l/b=1.0/1.5/2.0/2.5/3.0), 2
# MODAL passes (SYMM/ASYM) each, NMODES=20/pass, plus per-mode nodal UZ
# eigenvector dumps for an off-cluster MAC computation. Same one-seat,
# .QUEUE_OK resume as submit_ansys_queue.sh.
#
# MUST be submitted from Ansys/NewAnsys/ (not Ansys/).
# DO NOT run concurrent with any other MAPDL job. Check squeue first.
#
# After the run completes, the five *.QUEUE_OK sentinels confirm success;
# the 250 output files (rect_ff_oop_mac_{symm,asym}_freqs_lob<tag>.txt and
# rect_ff_oop_mac_{symm,asym}_m<1..20>_lob<tag>.txt) are the input to the
# (not yet written) off-cluster MAC computation against this project's
# analytical OOP mode shapes.

module load ansys/2024R2 2>/dev/null || module load ansys 2>/dev/null

echo "Job $SLURM_JOB_ID on $(hostname)"
echo "Start: $(date)"

if [ ! -f run_ansys_queue.sh ]; then
  echo "FATAL: run_ansys_queue.sh not found in $(pwd) (must be Ansys/NewAnsys/)."
  exit 2
fi

export CPUS=${SLURM_CPUS_PER_TASK:-8}
export MANIFEST=ansys_queue_manifest_rect_ff_oop_mac.txt
bash run_ansys_queue.sh

echo "End: $(date)"
