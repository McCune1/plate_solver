#!/bin/bash
#SBATCH --job-name=p5_seg
#SBATCH --output=ansys_queue_piezo_p5seg_%j.out
#SBATCH --error=ansys_queue_piezo_p5seg_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=01:00:00
#SBATCH --mem=16G
#SBATCH --requeue
#SBATCH --partition=requeue
#SBATCH --chdir=/home/ghmkfh/PythonMill/Plate_Solver_Package/Ansys/NewAnsys

# Paper 5 segmented-electrode sensing charge -- harmonic FE (2026-09-24).
# Decides the FUTURE_WORK "sensing admittance vs e15*gamma_rz" item, option (a).
# Pre-registration: Project Knowledge/PAPER5_E15_SEGMENT_CHARGE_DERIVATION_2026-09-24.md
# Six axisymmetric PLANE223 decks: C44E = 73 / 26 GPa x mesh 40x6, 80x12, 160x12.
# ONE ANSYS seat: this script calls run_ansys_queue.sh. Do NOT also sbatch
# the decks separately or submit concurrent ANSYS jobs.
#
# PUSH to /home/ghmkfh/PythonMill/Plate_Solver_Package/Ansys/NewAnsys/:
#   ansys_p5seg_c73_m40_2026-09-24.inp   ansys_p5seg_c26_m40_2026-09-24.inp
#   ansys_p5seg_c73_m80_2026-09-24.inp   ansys_p5seg_c26_m80_2026-09-24.inp
#   ansys_p5seg_c73_m160_2026-09-24.inp  ansys_p5seg_c26_m160_2026-09-24.inp
#   ansys_queue_manifest_piezo_p5seg_2026-09-24.txt
#   submit_ansys_queue_piezo_p5seg_2026-09-24.sh
#   targets_p5seg_2026-09-24.json
#   score_p5seg_fe_2026-09-24.py
#   build_p5seg_decks_2026-09-24.py   (optional; needs the predictions JSON)
# then from that directory:
#   sbatch submit_ansys_queue_piezo_p5seg_2026-09-24.sh
# after it ends (the script also runs it):
#   python3 score_p5seg_fe_2026-09-24.py | tee score_p5seg_fe_2026-09-24.txt
# Read the score file and each deck's _out.txt error count, never .QUEUE_OK alone.

echo "Job $SLURM_JOB_ID on $(hostname)"
echo "Start: $(date)"
if [ ! -f run_ansys_queue.sh ]; then
  echo "FATAL: must be Ansys/NewAnsys/"
  exit 2
fi
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
export CPUS=${SLURM_CPUS_PER_TASK:-8}
export MANIFEST=ansys_queue_manifest_piezo_p5seg_2026-09-24.txt
export STOP_ON_FAIL=0
bash run_ansys_queue.sh
echo "End: $(date)"
if [ -f score_p5seg_fe_2026-09-24.py ]; then
  python3 score_p5seg_fe_2026-09-24.py | tee score_p5seg_fe_2026-09-24.txt || true
fi
