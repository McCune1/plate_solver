#!/bin/bash
#SBATCH --job-name=ansys_ortho_ff
#SBATCH --output=ansys_ortho_ff_%j.out
#SBATCH --error=ansys_ortho_ff_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=01:30:00
#SBATCH --partition=general

# Independent regeneration of the orthotropic free-free (FFFF) flexural
# benchmark (2026-07-19 investigation): PAPER1_FREEFREE_DRAFT.tex claims
# tab:shi matches "the published orthotropic free-free tables of Shi,
# Liang and Wang" 8/8 to <=1.1%, against 8 targets LESSONS_LEARNED.md Sec 0
# labels "Ansys ortho-FFFF truth". Re-checking the actual cited paper
# (jve18517004.pdf) this session found none of those 8 numbers anywhere in
# it -- its own same-geometry FFFF row is 2-3 orders of magnitude smaller.
# This job does not try to settle that by more paper-reading -- it reruns
# the FE model fresh (SHELL281, orthotropic Er/Etheta/Grtheta/nu_r from
# Shi's paper, same R_i=4/R_o=8/PHI=90/H=0.08 physical plate as every other
# FF-P1 deck in this project) and reports the OOP-dominant frequency
# sequence directly, so the target list can be checked against a live run
# instead of a several-week-old note.
#
# Unlike the isotropic ansys_ff_oop_ext.inp (which uses a companion
# PLANE183 run to discard interleaved in-plane modes by frequency
# matching), this deck computes a per-mode UzFrac energy fraction
# (SUM(UZ^2)/SUM(UX^2+UY^2+UZ^2)) directly from the SHELL281 eigenvectors
# -- no companion IP model needed. NMODES=100 (vs the isotropic decks'
# 60) since this orthotropic material's in-plane stiffness (Er=40/
# Etheta=70 GPa) is much softer than the isotropic FF-P1 material
# (E=210 GPa), so more IP modes are expected to interleave below the same
# Omega_lit~100 ceiling -- generous margin so the 8th OOP target
# (~95.35) isn't cut off by running out of requested modes.
#
# NOT dry-run by the authoring session (no licensed ANSYS available there)
# -- built line-for-line against this project's already-validated
# ansys_ff_oop_ext.inp (geometry/mesh/two-pass structure) and
# ansys_ff_oop_mac_extract.inp (the *VOPER/*VSCFUN UzFrac idiom, already
# proven working on the cluster after its 2026-07-15 ABMAX fix) -- this is
# its first actual run.

module load ansys/2024R2 2>/dev/null || module load ansys 2>/dev/null

cd /home/ghmkfh/PythonMill/Plate_Solver_Package/Ansys/

if [ ! -f ansys_ff_ortho_oop.inp ]; then
  echo "FATAL: ansys_ff_ortho_oop.inp not found in $(pwd)."
  exit 2
fi

echo "Job $SLURM_JOB_ID on $(hostname)"
echo "Start: $(date)"

mapdl -b -smp -np ${SLURM_CPUS_PER_TASK:-8} -i ansys_ff_ortho_oop.inp -o ansys_ff_ortho_oop_out.txt

echo "--- ortho_ff_mesh64.txt / mesh96.txt existence check ---"
ls -la ortho_ff_mesh64.txt ortho_ff_mesh96.txt 2>&1

echo "--- NU_theta/C11BAR/DD/KLIT sanity lines (grep from the log) ---"
grep -E "MESH64:|MESH96:" ansys_ff_ortho_oop_out.txt

# --- error/warning scan -----------------------------------------------------
# Same discipline as submit_ansys_ff_oop_mac_extract.sh: a clean-looking
# banner above does not by itself prove the *DO/*VGET/*VSCFUN loop (100
# modes x 2 passes) ran to completion -- surface real APDL errors here
# rather than making a human dig through the full log by hand.
echo "--- error/warning scan (full log) ---"
grep -n -i -E '\*\*\* ERROR \*\*\*|\*\*\*ERROR\*\*\*|NUMBER OF ERROR MESSAGES|NUMBER OF WARNING MESSAGES|FATAL' \
  ansys_ff_ortho_oop_out.txt \
  || echo "(no error markers matched -- if the .txt files are still missing, inspect ansys_ff_ortho_oop_out.txt directly, and also check ansys_ortho_ff_${SLURM_JOB_ID}.err)"

# --- quick OOP-candidate preview (UzFrac > 0.9, Omega_lit not ~0) -----------
echo "--- mesh96 rows with UzFrac > 0.9 and Omega_lit > 1.0 (candidate OOP-flexural modes) ---"
if [ -f ortho_ff_mesh96.txt ]; then
  awk 'NR>2 { if ($4+0 > 0.9 && $3+0 > 1.0) print }' ortho_ff_mesh96.txt
else
  echo "ortho_ff_mesh96.txt missing -- see error scan above."
fi

# --- auto-copy to the package root -----------------------------------------
# Same reasoning as submit_ansys_ff_oop_mac_extract.sh (LESSONS_LEARNED Sec
# 13.5): whatever reads these results next runs from the package root, one
# folder up from here -- copy automatically rather than relying on a manual
# step between two separate SLURM submissions in two different folders.
PKG_ROOT="$(dirname "$(pwd)")"
echo "--- copying output files up to package root ($PKG_ROOT) ---"
for fname in ortho_ff_mesh64.txt ortho_ff_mesh96.txt; do
  if [ -f "$fname" ]; then
    cp -v "$fname" "$PKG_ROOT/$fname"
  else
    echo "WARNING: $fname was not produced by this run -- nothing to copy."
    echo "  Check the error scan above before assuming this is just a"
    echo "  missing-copy problem -- it may be a genuine run failure."
  fi
done

echo "End: $(date)"
