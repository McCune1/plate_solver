#!/bin/bash
#SBATCH --job-name=ansys_hs_parity
#SBATCH --output=ansys_hs_parity_%j.out
#SBATCH --error=ansys_hs_parity_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=01:30:00
#SBATCH --partition=requeue
# (partition switched from general -- billed since 2026-08-04)
#SBATCH --requeue
# Half-sector SYMM/ANTI parity classification of the FF-P1 free-free FE
# modes (addendum 2026-07-11 Sec I's "genuinely new lever" for the OOP
# ceiling item; LESSONS_LEARNED Sec 10.7 item 8). Four modal solves
# (SYMM/ANTI x mesh64/mesh96-equivalent) on a half sector, each smaller
# than the full-sector NMODES=60 runs that fit in 45 min -- 1.5h is ample.
#
# Independent of (can run in parallel with):
#   ../submit_ip_baseline_capture.sh
#   ../submit_ip_pairmember_finescan.sh
#
# See ansys_ff_oop_halfsector.inp header for the pre-registered
# convergence gate, calibration gate (four known-parity modes MUST match
# Sec 10.6-H before the 111.46 read), and the decisive check.
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

mapdl -b -smp -np ${SLURM_CPUS_PER_TASK:-8} -i ansys_ff_oop_halfsector.inp -o ansys_ff_oop_halfsector_out.txt

echo "=== SYMM (EVEN) mesh64-equiv ==="
cat oop_hs_symm_mesh64.txt
echo
echo "=== SYMM (EVEN) mesh96-equiv ==="
cat oop_hs_symm_mesh96.txt
echo
echo "=== ANTI (ODD) mesh64-equiv ==="
cat oop_hs_anti_mesh64.txt
echo
echo "=== ANTI (ODD) mesh96-equiv ==="
cat oop_hs_anti_mesh96.txt
echo
echo "Calibration gate (Sec 10.6-H): 15.128->SYMM, 24.190->ANTI,"
echo "37.818->ANTI, 53.828->SYMM. Then: which pass holds Omega_lit~111.46?"
