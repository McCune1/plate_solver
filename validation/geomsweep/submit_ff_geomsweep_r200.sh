#!/bin/bash
#SBATCH --job-name=ff_geomsweep_r200
#SBATCH --output=ff_geomsweep_r200_%j.out
#SBATCH --error=ff_geomsweep_r200_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=64
#SBATCH --time=32:00:00
#SBATCH --partition=general

# Free-free geometric-benchmarking sweep -- r0/(2b) = 2.0.
# Same rationale/cost notes as submit_ff_geomsweep_r150.sh. r0_2b>=2.0 is
# the regime _scan_cfg_part1/2 already widens xi_max/ze_max for (large-
# annulus branch．count escalation, same code path used by the validated
# r0_2b=2.5 cantilever point) -- not a new numerical regime, just not yet
# checked under free-free.

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

export R40_FF_P1=2.0:0.25,2.0:0.5,2.0:0.75,2.0:1.0,2.0:1.25,2.0:1.5
export R40_FF_P2=2.0:0.25,2.0:0.5,2.0:0.75,2.0:1.0,2.0:1.25,2.0:1.5

export MAKE_FIGS=0
export FIGDIR=figures_ff_r200
export CHECKPOINT=research23_checkpoint_ff_r200.json
export RESUME=1

module load python/3.12.1
source /home/ghmkfh/PythonMill/Plate_Solver_Package/venv/bin/activate

cd /home/ghmkfh/PythonMill/Plate_Solver_Package/
mkdir -p "$FIGDIR"

echo "Job $SLURM_JOB_ID on $(hostname): cpus-per-task=$SLURM_CPUS_PER_TASK, N_WORKERS=$N_WORKERS"
echo "Start: $(date)"

python3 -u -m plate_solver.cli

echo "End: $(date)"
