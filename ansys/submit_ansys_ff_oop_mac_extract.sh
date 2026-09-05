#!/bin/bash
#SBATCH --job-name=ansys_mac_extract
#SBATCH --output=ansys_mac_extract_%j.out
#SBATCH --error=ansys_mac_extract_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=00:45:00
#SBATCH --partition=general

# Extracts the FE eigenvector (nodal UX/UY/UZ) for the OOP mode nearest
# Omega_lit=111.4627 -- the still-open ceiling ambiguity (Sec 20.4) between
# two solver candidates (raw 2.828013 / 2.838515). Same model as the
# validated ansys_ff_oop_ext.inp (job 2316746); this deck only adds the
# dynamic mode-location + nodal eigenvector dump. Modal solves on this
# mesh size are fast (minutes), not the bottleneck -- 45 min is generous.
#
# Output: mac_eigvec_mesh64.txt / mac_eigvec_mesh96.txt (r, theta_deg, UX,
# UY, UZ per node) plus /COM lines in the .out log reporting which mode
# was selected and its UZ-dominance -- READ THOSE FIRST, before trusting
# the eigenvector (see the deck's own header/footer for what to check).
#
# NOT dry-run by the authoring session (no licensed ANSYS available there)
# -- built line-for-line against ansys_ff_oop_ext.inp's already-validated
# model block and ansys_ff_isoval_tagged.inp's already-used nodal-readout
# idiom, but this is its first actual run.

module load ansys/2024R2 2>/dev/null || module load ansys 2>/dev/null

cd /home/ghmkfh/PythonMill/Plate_Solver_Package/Ansys/

if [ ! -f ansys_ff_oop_mac_extract.inp ]; then
  echo "FATAL: ansys_ff_oop_mac_extract.inp not found in $(pwd)."
  exit 2
fi

echo "Job $SLURM_JOB_ID on $(hostname)"
echo "Start: $(date)"

mapdl -b -smp -np ${SLURM_CPUS_PER_TASK:-8} -i ansys_ff_oop_mac_extract.inp -o ansys_ff_oop_mac_extract_out.txt

echo "--- mac_eigvec_mesh64.txt / mesh96.txt existence check ---"
ls -la mac_eigvec_mesh64.txt mac_eigvec_mesh96.txt 2>&1

echo "--- mode-selection / UZ-dominance summary (grep from the log) ---"
grep -E "MESH64|MESH96" ansys_ff_oop_mac_extract_out.txt

# --- error/warning scan -----------------------------------------------------
# Added 2026-07-15 after job 2323987: the MESH64/MESH96 grep above only
# catches this deck's OWN /COM banners, so a real APDL fatal error elsewhere
# in the deck (e.g. in the *GET/*VGET/*VSCFUN/*CFOPEN eigenvector-dump block)
# silently disappears from this summary even though it's what actually
# aborted the run and left the .txt files missing. Surface it here instead of
# making a human dig through the full log by hand.
echo "--- error/warning scan (full log -- this is what actually killed a run" \
     "that has a clean MESH64/MESH96 banner above but no output files) ---"
grep -n -i -E '\*\*\* ERROR \*\*\*|\*\*\*ERROR\*\*\*|NUMBER OF ERROR MESSAGES|NUMBER OF WARNING MESSAGES|FATAL' \
  ansys_ff_oop_mac_extract_out.txt \
  || echo "(no error markers matched -- if the eigenvector files are still missing, inspect ansys_ff_oop_mac_extract_out.txt directly, and also check ansys_mac_extract_${SLURM_JOB_ID}.err)"

# --- auto-copy to the package root -----------------------------------------
# probe_mac_ambiguity_check.py runs from the PACKAGE ROOT (one folder up from
# this Ansys/ folder), not from here -- this project always runs Ansys decks
# from Ansys/ and Python probes from the root above it. Previously the deck's
# own footer comment told a human to copy mac_eigvec_mesh96.txt up by hand;
# that manual step is exactly the kind of thing that's easy to forget between
# two separate SLURM submissions in two different folders, and was the
# reason probe_mac_ambiguity_check.py had nothing to read on the first
# attempt (2026-07-15). Do the copy here instead, automatically, right after
# the file is known to exist -- see LESSONS_LEARNED Sec 13.5.
PKG_ROOT="$(dirname "$(pwd)")"
echo "--- copying eigenvector files up to package root ($PKG_ROOT) ---"
for fname in mac_eigvec_mesh64.txt mac_eigvec_mesh96.txt; do
  if [ -f "$fname" ]; then
    cp -v "$fname" "$PKG_ROOT/$fname"
  else
    echo "WARNING: $fname was not produced by this run -- nothing to copy."
    echo "  Check the MESH64/MESH96 /COM lines above and the .out log for a"
    echo "  mode-selection or UZ-dominance failure before assuming this is"
    echo "  just a missing-copy problem -- it may be a genuine extraction bug."
  fi
done

echo "End: $(date)"
