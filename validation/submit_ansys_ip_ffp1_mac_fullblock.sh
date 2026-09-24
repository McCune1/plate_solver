#!/bin/bash
#SBATCH --job-name=ffp1_ip_mac
#SBATCH --output=ansys_ip_ffp1_mac_fullblock_%j.out
#SBATCH --error=ansys_ip_ffp1_mac_fullblock_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=03:00:00
#SBATCH --partition=requeue
#SBATCH --requeue

# FF-P1 (r0/2b=1.5, 2Theta/pi=0.5, nu=0.30) IP eigenvector dump.
# Distinct prefix ip_mac_ffp1_nu030 so it cannot collide with the
# Rank 14 nu=0.35 dump ip_mac_gs_r150_a050_*.
# Prerequisite for submit_ip_ffp1_1115_pair_mac_v1.sh and
# submit_ip_qin_pair_mac_v1.sh (--dependency=afterok:$this).

export RATIO_TAG=r150
export ANGLE_TAG=a050
export GS_NU=0.30
export FB_PREFIX=ip_mac_ffp1_nu030
export FB_DECK=ansys_ip_ffp1_nu030_mac_fullblock.inp
export FB_JOB=ffp1ipmac
export PKG_PATH=/home/ghmkfh/PythonMill/Plate_Solver_Package

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
module load python/3.12.1
source /home/ghmkfh/PythonMill/Plate_Solver_Package/venv/bin/activate

cd /home/ghmkfh/PythonMill/Plate_Solver_Package/ || {
  echo "FATAL: cannot cd to package root."
  exit 2
}

if [ ! -f generate_ansys_ip_geomsweep_mac_fullblock.py ]; then
  echo "FATAL: generator not in $(pwd)."
  exit 2
fi

python3 -u generate_ansys_ip_geomsweep_mac_fullblock.py || {
  echo "FATAL: generator failed."
  exit 2
}

cd /home/ghmkfh/PythonMill/Plate_Solver_Package/Ansys/NewAnsys || {
  echo "FATAL: cannot cd to Ansys/NewAnsys."
  exit 2
}

DECK="$FB_DECK"
PREFIX="$FB_PREFIX"
if [ ! -f "$DECK" ]; then
  echo "FATAL: $DECK not in $(pwd) after generate."
  exit 2
fi

echo "Job $SLURM_JOB_ID on $(hostname)"
echo "Working directory: $(pwd)"
echo "GS_NU=$GS_NU PREFIX=$PREFIX JOB=$FB_JOB"
echo "Start: $(date)"

mapdl -b -smp -np ${SLURM_CPUS_PER_TASK:-8} -j "$FB_JOB" -i "$DECK" \
      -o "ansys_ip_ffp1_nu030_mac_fullblock_out.txt" || {
  echo "MAPDL_FAILED (exit $?) -- stopping before the summary greps."
  exit 1
}

OUTTXT="ansys_ip_ffp1_nu030_mac_fullblock_out.txt"
if head -40 "$OUTTXT" 2>/dev/null | grep -q -i "no input file specified"; then
  echo "FATAL: MAPDL never read the deck."
  exit 1
fi

echo "--- frequency table ---"
ls -la "${PREFIX}_freqs.txt" 2>&1
head -8 "${PREFIX}_freqs.txt" 2>/dev/null

echo "--- eigenvector file count (57 expected) ---"
N=$(ls -1 ${PREFIX}_m*_mesh96.txt 2>/dev/null | wc -l)
echo "found $N"
if [ "$N" -lt 57 ]; then
  echo "WARNING: only $N/57 mode files."
fi

echo "--- error/warning scan ---"
grep -n -i -E '\*\*\* ERROR \*\*\*|NUMBER OF ERROR MESSAGES|NUMBER OF WARNING MESSAGES|FATAL' \
  "$OUTTXT" || echo "(no error markers matched)"

echo "End: $(date)"
