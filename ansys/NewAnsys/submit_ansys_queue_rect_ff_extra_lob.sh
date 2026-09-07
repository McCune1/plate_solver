#!/bin/bash
#SBATCH --job-name=ansys_rect_ff_xtra
#SBATCH --output=ansys_rect_ff_xtra_%j.out
#SBATCH --error=ansys_rect_ff_xtra_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=03:00:00
#SBATCH --mem=16G
#SBATCH --requeue
#SBATCH --partition=requeue
#
# Isolated wrapper around this directory's run_ansys_queue.sh, scoped to
# ONLY ansys_queue_manifest_rect_ff_extra_lob.txt (l/b=2.0 and l/b=3.0
# half-model FFFF OOP decks) via the MANIFEST env var -- deliberately
# does NOT touch ansys_queue_manifest.txt or the already-run p5phaseA
# manifest (lob100 / lob250). Same one-seat-license, run-one-at-a-time
# discipline as submit_ansys_queue.sh, same .QUEUE_OK resume-skip sentinel.
#
# --time=03:00:00 is headroom: each deck is 4 MODAL passes (SHELL281,
# NMODES=30). l/b=2.0 is NDIVX=80/120; l/b=3.0 is NDIVX=120/180 (larger
# than the already-run lob250 100/150 pair). Expect ~45-90 min of MAPDL;
# the rest is license-queue headroom.
#
# MUST be submitted from Ansys/NewAnsys/ (not Ansys/) -- see
# run_ansys_queue.sh's own header for why (MAPDL only finds its input file
# from this exact directory on this cluster).
#
# DO NOT run this at the same time as submit_ansys_queue.sh, the p5phaseA
# queue, or any other MAPDL job -- one ANSYS license seat, project-wide.
# Check squeue first. After the job: convert Hz with KFAC=24.6644 from
# the first mode-1 *VWRITE block in
#   rect_ff_{symm,asym}_mesh{1,2}_lob{200,300}.txt

module load ansys/2024R2 2>/dev/null || module load ansys 2>/dev/null

echo "Job $SLURM_JOB_ID on $(hostname)"
echo "Start: $(date)"

if [ ! -f run_ansys_queue.sh ]; then
  echo "FATAL: run_ansys_queue.sh not found in $(pwd) (must be Ansys/NewAnsys/)."
  exit 2
fi

export CPUS=${SLURM_CPUS_PER_TASK:-8}
export MANIFEST=ansys_queue_manifest_rect_ff_extra_lob.txt
bash run_ansys_queue.sh

echo "End: $(date)"
