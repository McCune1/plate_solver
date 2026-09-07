#!/bin/bash
#SBATCH --job-name=ansys_rect_ff_ip
#SBATCH --output=ansys_rect_ff_ip_%j.out
#SBATCH --error=ansys_rect_ff_ip_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=04:00:00
#SBATCH --mem=16G
#SBATCH --requeue
#SBATCH --partition=requeue
#
# Isolated wrapper around run_ansys_queue.sh for
# ansys_queue_manifest_rect_ff_ip.txt only -- five PLANE183 half-model
# FFFF IP decks (l/b=1.0/1.5/2.0/2.5/3.0), 4 MODAL passes each.
# Same one-seat, .QUEUE_OK resume as submit_ansys_queue.sh.
#
# MUST be submitted from Ansys/NewAnsys/ (not Ansys/).
# DO NOT run concurrent with any other MAPDL job. Check squeue first.
# After: Omega = f[Hz] / 804.481 from
#   rect_ff_ip_{symm,asym}_mesh{1,2}_lob{100,150,200,250,300}.txt

module load ansys/2024R2 2>/dev/null || module load ansys 2>/dev/null

echo "Job $SLURM_JOB_ID on $(hostname)"
echo "Start: $(date)"

if [ ! -f run_ansys_queue.sh ]; then
  echo "FATAL: run_ansys_queue.sh not found in $(pwd) (must be Ansys/NewAnsys/)."
  exit 2
fi

export CPUS=${SLURM_CPUS_PER_TASK:-8}
export MANIFEST=ansys_queue_manifest_rect_ff_ip.txt
bash run_ansys_queue.sh

echo "End: $(date)"
