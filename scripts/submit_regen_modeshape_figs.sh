#!/bin/bash
#SBATCH --job-name=regen_modeshape_figs
#SBATCH --output=regen_modeshape_figs_%j.out
#SBATCH --error=regen_modeshape_figs_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=00:40:00
#SBATCH --partition=general

# 2026-07-19: bumped 8->32 cpus-per-task (cluster reported open). The
# --time=00:20:00 first tried here was WRONG: it was computed from this
# comment block's stale "33 modes" claim. Job 2324835 (32 workers) was
# CANCELLED AT THE TIME LIMIT at 1220s with 80/81 modes done -- the real
# checkpoint now has 81 modes (5 geometries: 2*Theta=0.25pi/0.5pi/0.75pi/
# 1.0pi/1.25pi at r0/(2b)=1.25 with 5/12/16/21/27 modes respectively), not
# 33 -- the geometry set grew after this comment was written and was never
# updated. probe_regen_modeshape_figs.py has NO resume/skip-if-exists
# logic, so a resubmission redoes all 81 modes from scratch (not just the
# missing one) at the observed ~1220s/81 modes throughput -- --time=00:40:00
# gives ~2x margin over that (~20 min) rather than a tight re-estimate,
# consistent with this project's history of underestimating these costs.
#
# Regenerates every Part-1 (out-of-plane) mode-shape *_fig3.png figure in
# figures/ from research23_checkpoint.json's already-recorded Omega values,
# picking up the B18 fix (LESSONS_LEARNED bug ledger, Sec 3, and Sec 13.5 /
# 20.4 addenda, all dated 2026-07-15) to `OutOfPlaneSolver.mode_shape_grid`
# -- cos(xi*theta+ph) was used where the validated, frequency-determining
# _build_K_real assembly uses sin(...); a 90-degree theta-parity error that
# affected every mode-shape figure this project has ever produced, but NO
# frequency/table value (mode_shape_grid is not on the full_search/
# select_fill/_build_K_real path that finds Omega) -- so this is a pure
# re-render, not a re-discovery: it does NOT re-run the original
# sigma_min discovery sweep (hours per geometry on this checkpoint's 3
# geometries), only full_search + mode_shape_grid reconstruction per mode
# (~60-90s each), the same cost driver as probe_paper_figures_oop.py and
# the MAC probe.
#
# Covers: r0/(2b)=1.25 at 2*Theta=0.25pi (5), 0.5pi (12), 0.75pi (16),
# 1.0pi (21), 1.25pi (27) -- 81 modes total (corrected count, see above;
# the old "33 modes / 3 geometries" comment was stale). CONFIRMED (2026-07-
# 19, direct JSON read) that research23_checkpoint_p1tail.json DUPLICATES
# the 1.0pi (21) and 1.25pi (27) entries from THIS checkpoint verbatim (48
# modes computed twice) -- harmless (different FIGDIR: figures/ vs
# figures_p1tail/, no file collision) but wasted compute if both jobs run;
# not changed here since it's unclear which checkpoint is authoritative --
# flagged for the user to decide, not fixed unilaterally.
#
# paper_figures/'s 4 mode-shape PNGs are a SEPARATE deliverable (produced
# by probe_paper_figures_oop.py directly, not from this checkpoint) --
# that one just needs a plain resubmission of its existing, unchanged
# submit_paper_figures_oop.sh once this fixed core_solvers.py is on the
# cluster; no new script needed for it.
#
# Part 2 (in-plane) entries in the checkpoint have no mode-shape figure
# counterpart (plot_mode_shape/mode_shape_grid is OOP-only) and are
# silently skipped by the probe -- figures_p2/ is unaffected by B18.
#
# PUSH BEFORE RUNNING: core_solvers.py (the B18 fix itself) and this
# script + probe_regen_modeshape_figs.py must all be on the cluster first.
# PULL AFTER RUNNING: the .out log (to confirm all 33 modes regenerated
# OK/suspect as expected, no FAILED lines) -- the PNGs themselves are
# overwritten in place in figures/ and don't need pulling back into this
# chat, just verify the job log looks clean.

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export BLIS_NUM_THREADS=1

export PKG_PATH=/home/ghmkfh/PythonMill/Plate_Solver_Package
export EXPECT_SOLVER_VERSION=2026-07-10.s8
export CHECKPOINT=research23_checkpoint.json
export FIGDIR=figures
export N_WORKERS=32

module load python/3.12.1
source /home/ghmkfh/PythonMill/Plate_Solver_Package/venv/bin/activate

cd /home/ghmkfh/PythonMill/Plate_Solver_Package/

if [ ! -f probe_regen_modeshape_figs.py ]; then
  echo "FATAL: probe_regen_modeshape_figs.py not found in $(pwd)."
  exit 2
fi
if [ ! -f "$CHECKPOINT" ]; then
  echo "FATAL: $CHECKPOINT not found in $(pwd)."
  exit 3
fi

echo "Job $SLURM_JOB_ID on $(hostname)"
echo "Start: $(date)"

python3 -u probe_regen_modeshape_figs.py

echo "End: $(date)"
