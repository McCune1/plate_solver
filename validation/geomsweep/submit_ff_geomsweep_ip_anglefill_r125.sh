#!/bin/bash
#SBATCH --job-name=ff_ip_anglefill_r125
#SBATCH --output=ff_ip_anglefill_r125_%j.out
#SBATCH --error=ff_ip_anglefill_r125_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=64
#SBATCH --time=16:00:00
#SBATCH --partition=general

# Cheap fill-in: the ORIGINAL FF-P1 free-free work only ran Part 2 (IP) at
# 3 of the 6 angles used for Part 1 (0.25, 0.5, 1.0 pi -- see
# submit_plate_solver.sh's R40_FF_P2). This job adds the 3 missing angles
# (0.75, 1.25, 1.5 pi) at the SAME r0_2b=1.25 ratio so the validated ratio
# has full 6-angle IP coverage matching Part 1, before spending compute on
# the new ratios' IP sweeps in the other 4 scripts in this folder.
# R40_FF_P1 is left EMPTY (no Part-1 work here -- already exhaustively done
# and closed, see LESSONS_LEARNED.md item 2).
# Small job: 3 geometries, own FIGDIR/checkpoint so it can run alongside
# the other geometry-sweep jobs without collision.

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export BLIS_NUM_THREADS=1

export N_WORKERS=${SLURM_CPUS_PER_TASK:-64}

export R40_OVERNIGHT=1
export R40_SPOTCHECK=1
export R40_FREEFREE=1
export R40_ORTHO_TEST=0
export R40_FF_MODES=6

export R40_FF_P1=
export R40_FF_P2=1.25:0.75,1.25:1.25,1.25:1.5

export MAKE_FIGS=0
export FIGDIR=figures_ff_ip_anglefill_r125
export CHECKPOINT=research23_checkpoint_ff_ip_anglefill_r125.json
export RESUME=1

module load python/3.12.1
source /home/ghmkfh/PythonMill/Plate_Solver_Package/venv/bin/activate

cd /home/ghmkfh/PythonMill/Plate_Solver_Package/
mkdir -p "$FIGDIR"

echo "Job $SLURM_JOB_ID on $(hostname): cpus-per-task=$SLURM_CPUS_PER_TASK, N_WORKERS=$N_WORKERS"
echo "Start: $(date)"

python3 -u -m plate_solver.cli

echo "End: $(date)"
