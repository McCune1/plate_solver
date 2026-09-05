#!/bin/bash
#SBATCH --job-name=geomsweep_fe_match
#SBATCH --output=geomsweep_fe_match_%j.out
#SBATCH --error=geomsweep_fe_match_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:10:00
#SBATCH --partition=general

# Runs probe_geomsweep_fe_match.py -- closes the reviewer-flagged Sec 6.4
# reproducibility gap ("the script that produced the 84.0%/1232/1466
# geometry-generalization statistic was never committed"). This is a
# from-scratch RECONSTRUCTION (2026-07-23), not a re-run of the lost
# original -- see the probe's own module docstring for the full account,
# including the one methodological choice (extensional-contamination
# removal tolerance) that could not be recovered and had to be re-decided.
#
# Pure Python stdlib (json/csv/re/glob) -- no plate_solver import, no
# mpmath, no SOLVER_VERSION dependency, no multiprocessing. Single CPU,
# runs in well under a minute; the SLURM wrapper exists only to give this
# evidence a job number consistent with the rest of this project's
# provenance trail, not because the compute needs a cluster.
#
# MUST run from FutureWork/geometry_sweep/ -- the manifests
# (figures_ff_r*/freefree_manifest.json) live here, and
# GEOMSWEEP_FE_DIR points at the sibling ansys_geomsweep_verification/
# directory where the independent FE per-angle results live. Do NOT point
# this at validation/geomsweep/ or ansys/geomsweep/ -- those are the
# curated github_repo/ copies; this run should be against the canonical
# working-tree data.

export GEOMSWEEP_MANIFEST_DIR=.
export GEOMSWEEP_FE_DIR=../ansys_geomsweep_verification
export GEOMSWEEP_MATCH_CSV=geomsweep_fe_match_consolidated.csv
export GEOMSWEEP_CONTAM_TOL=1.0

module load python/3.12.1
source /home/ghmkfh/PythonMill/Plate_Solver_Package/venv/bin/activate

if [ ! -f probe_geomsweep_fe_match.py ]; then
  echo "FATAL: probe_geomsweep_fe_match.py not found in $(pwd)."
  echo "Must run from FutureWork/geometry_sweep/."
  exit 2
fi
if [ ! -d "$GEOMSWEEP_FE_DIR" ]; then
  echo "FATAL: $GEOMSWEEP_FE_DIR not found relative to $(pwd)."
  echo "Expected FutureWork/ansys_geomsweep_verification/ as a sibling of"
  echo "FutureWork/geometry_sweep/."
  exit 2
fi

echo "Job $SLURM_JOB_ID on $(hostname)"
echo "Start: $(date)"

python3 -u probe_geomsweep_fe_match.py

echo "End: $(date)"
