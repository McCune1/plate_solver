#!/bin/bash
#SBATCH --job-name=ansys_ip_ext
#SBATCH --output=ansys_ip_ext_%j.out
#SBATCH --error=ansys_ip_ext_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=00:45:00
#SBATCH --partition=requeue
# (partition switched from general -- billed since 2026-08-04)
#SBATCH --requeue
# Recreates the isotropic free-free extensional (in-plane) PLANE183
# benchmark as a genuine two-mesh (64/96) pair (2026-07-19 investigation):
# PAPER1_FREEFREE_DRAFT.tex Sec 6.3 claims this model was "meshed to the
# same two-mesh <0.3% standard" as the flexural SHELL281 benchmark, but no
# PLANE183 deck was found anywhere in the local project tree -- only its
# output filename (ip_modes_mesh80.txt, job 2316743) is referenced, always
# at a single mesh density. Likely lost in one of this project's several
# documented cleanup passes, not evidence the original work didn't happen.
# This job reruns the model fresh at both mesh densities so the two-mesh
# claim is backed by an actual matching run.
#
# Small 2D plane-stress model (PLANE183, no shell/orthotropic complexity) --
# 45 min is generous for two passes at NMODES=35.
#
# NOT dry-run by the authoring session (no licensed ANSYS available there)
# -- built off this project's already-validated ansys_ff_oop_ext.inp
# geometry/two-pass structure, swapping SHELL281 for PLANE183 in plane
# stress and dropping the ESYS/orthotropic machinery (isotropic material
# needs no directional orientation) -- this is its first actual run.

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

if [ ! -f ansys_ff_ip_ext.inp ]; then
  echo "FATAL: ansys_ff_ip_ext.inp not found in $(pwd)."
  exit 2
fi

echo "Job $SLURM_JOB_ID on $(hostname)"
echo "Start: $(date)"

mapdl -b -smp -np ${SLURM_CPUS_PER_TASK:-8} -i ansys_ff_ip_ext.inp -o ansys_ff_ip_ext_out.txt

echo "--- ip_ext_mesh64.txt / mesh96.txt existence check ---"
ls -la ip_ext_mesh64.txt ip_ext_mesh96.txt 2>&1

# --- error/warning scan -----------------------------------------------------
echo "--- error/warning scan (full log) ---"
grep -n -i -E '\*\*\* ERROR \*\*\*|\*\*\*ERROR\*\*\*|NUMBER OF ERROR MESSAGES|NUMBER OF WARNING MESSAGES|FATAL' \
  ansys_ff_ip_ext_out.txt \
  || echo "(no error markers matched -- if the .txt files are still missing, inspect ansys_ff_ip_ext_out.txt directly, and also check ansys_ip_ext_${SLURM_JOB_ID}.err)"

# --- quick side-by-side preview (skip the 3 rigid modes) --------------------
echo "--- mesh64 vs mesh96, modes 4 up (side by side, first 30 rows each) ---"
if [ -f ip_ext_mesh64.txt ] && [ -f ip_ext_mesh96.txt ]; then
  paste <(awk 'NR>4 && NR<=34 {printf "%-8s %s\n",$1,$2}' ip_ext_mesh64.txt) \
        <(awk 'NR>4 && NR<=34 {printf "%-8s %s\n",$1,$2}' ip_ext_mesh96.txt)
else
  echo "one or both output files missing -- see error scan above."
fi

# --- auto-copy to the package root -----------------------------------------
PKG_ROOT="$(dirname "$(pwd)")"
echo "--- copying output files up to package root ($PKG_ROOT) ---"
for fname in ip_ext_mesh64.txt ip_ext_mesh96.txt; do
  if [ -f "$fname" ]; then
    cp -v "$fname" "$PKG_ROOT/$fname"
  else
    echo "WARNING: $fname was not produced by this run -- nothing to copy."
    echo "  Check the error scan above before assuming this is just a"
    echo "  missing-copy problem -- it may be a genuine run failure."
  fi
done

echo "End: $(date)"
