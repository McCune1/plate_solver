#!/bin/bash
#SBATCH --job-name=rect_ff_1982_wide
#SBATCH --output=rect_ff_1982_wide_%j.out
#SBATCH --error=rect_ff_1982_wide_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:30:00
#SBATCH --requeue
#SBATCH --partition=requeue
#
# Follow-up to 2455675: that job's +/-0.02 fine reopen window sat on a
# flat ~1.7e-6 sigma plateau with no 10x rise on either flank -- neither
# a clean isolated notch nor a clean falling tail within that narrow a
# window. This job widens the scan to [1.73,2.23] (step 0.002, 250 pts)
# at the same persist basis (6,3) to see whether the plateau ever walls
# off (real but broad feature) or keeps falling toward the next SYM FE
# root at 2.45439 (tail-of-neighbor). Purely descriptive; does not
# change the Table 3 MISS verdict for l/b=1.0 SYM 1.98249 either way.
# Cost: 2455675 ran 201 pts (0.04 wide, step 0.0002) in 215s at this
# basis, ~1.07 s/pt. 250 pts here -> ~270s. Expect well under 10 min.
# COPY: p5_rect_ff_lib.py, this probe, this script.

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export BLIS_NUM_THREADS=1
export DPS=30
export EXPECT_SOLVER_VERSION=2026-07-10.s10
export PKG_PATH=/home/ghmkfh/PythonMill/Plate_Solver_Package
export CHECKPOINT_DIR=./p5_1982_wide_gap_checkpoints

echo "Job $SLURM_JOB_ID on $(hostname)  start: $(date)"
module load python/3.12.1
source /home/ghmkfh/PythonMill/Plate_Solver_Package/venv/bin/activate
cd /home/ghmkfh/PythonMill/Plate_Solver_Package/
python3 -u probe_rect_ff_oop_1982_wide_gap_2026-09-05.py
echo "Job complete: $(date)"
