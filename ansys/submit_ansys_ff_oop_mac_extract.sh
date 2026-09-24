#!/bin/bash
#SBATCH --job-name=ansys_mac_extract
#SBATCH --output=ansys_mac_extract_%j.out
#SBATCH --error=ansys_mac_extract_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=00:45:00
#SBATCH --partition=requeue
# (partition switched from general -- billed since 2026-08-04)
#SBATCH --requeue
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

# ---- ANSYS environment (2026-09-23) -----------------------------------
# Replaces `module load ansys/2024R2 2>/dev/null || module load ansys
# 2>/dev/null`, which fails silently now: Mill dropped ansys/2024R2, and the
# default `ansys` (2026R1) is a container with no `mapdl` on PATH (job
# 2530336: every deck exit 127). Same logic as the block in
# Ansys/NewAnsys/submit_ansys_queue_piezo_p4_3d_nge1_2026-09-23.sh, which
# ran job 2530423 on ansys/2025R1. Tries ansys/2025R1, ansys/2024R1 (native
# launcher), then ansys/2026R1 through an UNVERIFIED container wrapper.
# Overrides: MAPDL_EXE=/path/to/launcher, ANSYS_MODULES="ansys/2024R1".
# Stops (exit 3) before running anything if no MAPDL launcher is found.
_ansys_e=0; _ansys_u=0; _ansys_pf=0
case $- in *e*) _ansys_e=1;; esac
case $- in *u*) _ansys_u=1;; esac
if shopt -qo pipefail; then _ansys_pf=1; fi
set +e +u +o pipefail
if ! type module >/dev/null 2>&1; then
  for _i in /etc/profile.d/lmod.sh /etc/profile.d/modules.sh /usr/share/lmod/lmod/init/bash; do
    if [ -f "$_i" ]; then . "$_i"; break; fi
  done
fi
_ansys_find() {
  if [ -n "${MAPDL_EXE:-}" ] && [ -x "$MAPDL_EXE" ]; then echo "$MAPDL_EXE"; return 0; fi
  if command -v mapdl 2>/dev/null; then return 0; fi
  local d f
  for d in $(printf '%s' "$PATH" | tr ':' ' '); do
    for f in "$d"/ansys[0-9][0-9][0-9]; do
      if [ -x "$f" ]; then echo "$f"; return 0; fi
    done
  done
  return 0
}
_ansys_shim="${TMPDIR:-/tmp}/mapdl_shim_${SLURM_JOB_ID:-$$}"
_ansys_mods="${ANSYS_MODULES:-ansys/2025R1 ansys/2024R1 ansys/2026R1}"
_ansys_exe="$(_ansys_find | head -1)"
if [ -z "$_ansys_exe" ] && type module >/dev/null 2>&1; then
  for _m in $_ansys_mods; do
    module unload ansys >/dev/null 2>&1
    module load "$_m" || continue
    _ansys_exe="$(_ansys_find | head -1)"
    if [ -n "$_ansys_exe" ]; then echo "ANSYS module: $_m"; break; fi
    _l="$(type ansys 2>/dev/null | grep -o '/[^" ]*ansys_launcher_logic\.sh' | head -1)"
    if [ -n "$_l" ] && [ -x "$_l" ]; then
      mkdir -p "$_ansys_shim"
      printf '#!/bin/bash\nexec "%s" mapdl "$@"\n' "$_l" > "$_ansys_shim/mapdl_container"
      chmod +x "$_ansys_shim/mapdl_container"
      _ansys_exe="$_ansys_shim/mapdl_container"
      echo "ANSYS module: $_m (CONTAINER wrapper, unverified: read the first transcript)"
      break
    fi
  done
fi
if [ -n "$_ansys_exe" ] && [ "$_ansys_exe" != "$(command -v mapdl 2>/dev/null)" ]; then
  mkdir -p "$_ansys_shim" && ln -sf "$_ansys_exe" "$_ansys_shim/mapdl"
  export PATH="$_ansys_shim:$PATH"
fi
if ! command -v mapdl >/dev/null 2>&1; then
  echo "FATAL: no MAPDL launcher found (tried: $_ansys_mods). Nothing was run."
  type module >/dev/null 2>&1 && module spider ansys 2>&1 | sed 's/^/  /'
  exit 3
fi
echo "MAPDL launcher: $_ansys_exe"
if [ "$_ansys_e" = 1 ]; then set -e; fi
if [ "$_ansys_u" = 1 ]; then set -u; fi
if [ "$_ansys_pf" = 1 ]; then set -o pipefail; fi
# ---- end ANSYS environment ---------------------------------------------

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
