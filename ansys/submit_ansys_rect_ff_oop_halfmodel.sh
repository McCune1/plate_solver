#!/bin/bash
#SBATCH --job-name=ansys_rect_ff
#SBATCH --output=ansys_rect_ff_%j.out
#SBATCH --error=ansys_rect_ff_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=01:30:00
#SBATCH --partition=general

# First FE cross-check for the rectangular free-free (FFFF) OOP candidates
# (job 2322815 discovery + job 2322845 basis-bump confirmation, both corner
# variants agreeing on Lambda* for 6 of 7 candidates -- LESSONS_LEARNED
# item 16). Four modal solves (SYMM/ASYM half-model x mesh1/mesh2) on a
# flat rectangular half-plate, mirroring the annular half-sector parity
# technique but in plain global Cartesian coordinates (no local CS needed
# for a flat-plate mirror plane). EXPLORATORY -- first FE pass on one
# geometry instance, not yet a validated gate. No SOLVER_VERSION action.
#
# See ansys_rect_ff_oop_halfmodel.inp header for the full Lambda->Hz
# derivation (independently re-derived from Seok/Tiersten/Scarton 2004
# rect Part 1 Eqs. 24-28/37-38/41-42 and cross-checked against
# RectOOPAssembler's isotropic T=R=1, mu=1.7 reduction), the solver-side
# reference table, and the pre-registered reading order (convergence gate
# first, parity-matched comparison second, common-factor-miss check third).

module load ansys/2024R2 2>/dev/null || module load ansys 2>/dev/null

mapdl -b -smp -np ${SLURM_CPUS_PER_TASK:-8} -i ansys_rect_ff_oop_halfmodel.inp -o ansys_rect_ff_oop_halfmodel_out.txt

echo "=== SYMM (EVEN) mesh1 (x60/y20) ==="
cat rect_ff_symm_mesh1.txt
echo
echo "=== SYMM (EVEN) mesh2 (x90/y30, convergence) ==="
cat rect_ff_symm_mesh2.txt
echo
echo "=== ASYM (ODD) mesh1 (x60/y20) ==="
cat rect_ff_asym_mesh1.txt
echo
echo "=== ASYM (ODD) mesh2 (x90/y30, convergence) ==="
cat rect_ff_asym_mesh2.txt
echo
echo "Solver reference (job 2322845, Lambda*): SYM 0.056/0.412/0.431/0.864/~0.96,"
echo "ANTI 0.475 (ANTI ~0.97 excluded, not yet trustworthy). Read against"
echo "the pre-registered criteria in ansys_rect_ff_oop_halfmodel.inp's header:"
echo "convergence gate first, parity-matched comparison second, common-factor"
echo "miss check third."
