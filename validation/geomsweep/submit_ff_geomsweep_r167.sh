#!/bin/bash
#SBATCH --job-name=ff_geomsweep_r167
#SBATCH --output=ff_geomsweep_r167_%j.out
#SBATCH --error=ff_geomsweep_r167_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=64
#SBATCH --time=32:00:00
#SBATCH --partition=general

# Free-free geometric-benchmarking sweep -- r0/(2b) = 5/3 ~ 1.66667.
# Same rationale/cost notes as submit_ff_geomsweep_r150.sh (read that file's
# header first). This ratio already has ONE tabulated cantilever point in
# PAPER_PART1 (r0_2b=5\3, mode 1, 2Theta=pi -> 0.01805), so it is a
# geometry regime the paper-faithful detector is already known to handle
# correctly in the clamped-free case; this job checks the SAME ratio under
# free-free across the full 6-angle grid, both parts.

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
export R40_ORTHO_TEST=1
export R40_FF_MODES=6

export R40_FF_P1=1.66667:0.25,1.66667:0.5,1.66667:0.75,1.66667:1.0,1.66667:1.25,1.66667:1.5
export R40_FF_P2=1.66667:0.25,1.66667:0.5,1.66667:0.75,1.66667:1.0,1.66667:1.25,1.66667:1.5

export MAKE_FIGS=0
export FIGDIR=figures_ff_r167
export CHECKPOINT=research23_checkpoint_ff_r167.json
export RESUME=1

module load python/3.12.1
source /home/ghmkfh/PythonMill/Plate_Solver_Package/venv/bin/activate

cd /home/ghmkfh/PythonMill/Plate_Solver_Package/
mkdir -p "$FIGDIR"

echo "Job $SLURM_JOB_ID on $(hostname): cpus-per-task=$SLURM_CPUS_PER_TASK, N_WORKERS=$N_WORKERS"
echo "Start: $(date)"

python3 -u -m plate_solver.cli

echo "End: $(date)"
