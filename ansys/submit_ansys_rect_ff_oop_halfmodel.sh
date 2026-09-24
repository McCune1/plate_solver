#!/bin/bash
#SBATCH --job-name=ansys_rect_ff
#SBATCH --output=ansys_rect_ff_%j.out
#SBATCH --error=ansys_rect_ff_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=01:30:00
#SBATCH --partition=requeue
# (partition switched from general -- billed since 2026-08-04)
#SBATCH --requeue
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
