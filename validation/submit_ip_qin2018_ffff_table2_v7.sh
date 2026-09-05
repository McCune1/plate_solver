# submit_ip_qin2018_ffff_table2_v7.sh
#!/bin/bash
#SBATCH --job-name=qin_ff7
#SBATCH --output=probe_ip_qin2018_ffff_table2_v7_%j.out
#SBATCH --error=probe_ip_qin2018_ffff_table2_v7_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=06:00:00
#SBATCH --mem=8G

echo "job=$SLURM_JOB_ID host=$(hostname) start=$(date)"
echo "SLURM_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMBA_NUM_THREADS=1

module load python/3.12.1
source /home/ghmkfh/PythonMill/Plate_Solver_Package/venv/bin/activate

cd /home/ghmkfh/PythonMill/Plate_Solver_Package

export DPS=40
export EXPECT_SOLVER_VERSION=2026-07-10.s10
export N_WORKERS=$SLURM_CPUS_PER_TASK

echo "N_WORKERS=$N_WORKERS (should be 32)"

python3 -u probe_ip_qin2018_ffff_table2_v7.py

echo "end=$(date)"