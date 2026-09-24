#!/bin/bash
#SBATCH --job-name=p4_liu_disk_sc
#SBATCH --output=ansys_queue_piezo_liu_disk_%j.out
#SBATCH --error=ansys_queue_piezo_liu_disk_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH --requeue
#SBATCH --partition=requeue
#SBATCH --chdir=/home/ghmkfh/PythonMill/Plate_Solver_Package/Ansys/NewAnsys

# Paper 4 lever 11 -- Liu 2002 solid-disk SC FE (Gate 2 CPT Tables 2 & 6).
# Mesh fix after job 2530974 (RI_MESH=3 mm annular LESIZE). Makise primary:
#   github_repo/ansys/NewAnsys (lowercase); this chdir is Mill Ansys/NewAnsys.
# ONE ANSYS seat: this script calls run_ansys_queue.sh. Do NOT also sbatch
# the decks separately or submit concurrent ANSYS jobs.
#
# PUSH to /home/ghmkfh/PythonMill/Plate_Solver_Package/Ansys/NewAnsys/:
#   ansys_p4_liu_disk_sc_C_2026-09-24.inp
#   ansys_p4_liu_disk_sc_S_2026-09-24.inp
#   ansys_queue_manifest_piezo_liu_disk_2026-09-24.txt
#   submit_ansys_queue_piezo_liu_disk_2026-09-24.sh
#   targets_liu_disk_2026-09-24.json
#   LIU_DISK_FE_README.md
#   build_liu_disk_sc_decks_2026-09-24.py
# then from that directory:
#   sbatch submit_ansys_queue_piezo_liu_disk_2026-09-24.sh
#
# Two SOLID226 moderate decks (~6912 elems, O(2e5) DOF). Read each
# deck's _modes.txt / _samples.txt / _out.txt. Never trust .QUEUE_OK alone.
# Recompute % externally from raw Hz vs targets_liu_disk_2026-09-24.json.
#
# Env block matches sibling piezo submits (2026-09-23). Resources sized
# for SOLID226 (not the 15min/8G PLANE223 scboth defaults).

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
export MANIFEST=ansys_queue_manifest_piezo_liu_disk_2026-09-24.txt
export STOP_ON_FAIL=0
bash run_ansys_queue.sh
echo "End: $(date)"
