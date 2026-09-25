#!/bin/bash
#SBATCH --job-name=p5x_pcm_e31m
#SBATCH --output=ansys_queue_combined_2026-09-24b_%j.out
#SBATCH --error=ansys_queue_combined_2026-09-24b_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=04:00:00
#SBATCH --mem=16G
#SBATCH --requeue
#SBATCH --partition=requeue
#SBATCH --chdir=/home/ghmkfh/PythonMill/Plate_Solver_Package/Ansys/NewAnsys

# ONE job, ONE ANSYS seat, 80 decks (LESSONS Sec 18.241-18.244):
#   ansys_p5mix_e31m_*   Paper 5 mixed-edge sensing, e31 = -4.1        (24)
#   ansys_pcm_*          contour-mode n = 0, e31 = +-4.1               (12)
#   ansys_parashar_*     Parashar 2013 scout (PIC 255, radial poling)  (6)
#   *_e31m_2026-09-24    e31 = -4.1 clones of every coupled P4/P5/P6 deck (32)
#   ansys_p5seg_e31m_*   Paper 5 segment charge (e15), e31 = -4.1      (6)
# At the end it runs the four scorers and bundle_results_2026-09-24b.py, which
# writes results_2026-09-24b_<jobid>.zip with EVERY output needed -- download
# that one file.
#
# PUSH to /home/ghmkfh/PythonMill/Plate_Solver_Package/Ansys/NewAnsys/:
#   the contents of _push_NewAnsys_2026-09-24b.zip (decks, manifest, this
#   script, targets JSONs, scorers, bundler). Parent decks of the clones are
#   already on the cluster and are not needed to run.
# then:  sbatch submit_ansys_queue_combined_2026-09-24b.sh
# Do NOT also submit any per-family script while this runs (one seat).

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
export MANIFEST=ansys_queue_manifest_combined_2026-09-24b.txt
export STOP_ON_FAIL=0
bash run_ansys_queue.sh
echo "End: $(date)"
# ---- scoring (never fatal) ----
python3 score_p5mix_fe_2026-09-24.py targets_p5mix_e31m_2026-09-24.json > score_p5mix_e31m_fe_2026-09-24.txt 2>&1 || true
python3 score_pcm_fe_2026-09-24.py > score_pcm_fe_2026-09-24.txt 2>&1 || true
python3 score_parashar_fe_2026-09-24.py > score_parashar_fe_2026-09-24.txt 2>&1 || true
python3 score_p5seg_fe_2026-09-24.py targets_p5seg_e31m_2026-09-24.json > score_p5seg_e31m_fe_2026-09-24.txt 2>&1 || true
tail -n 3 score_p5mix_e31m_fe_2026-09-24.txt score_pcm_fe_2026-09-24.txt score_parashar_fe_2026-09-24.txt score_p5seg_e31m_fe_2026-09-24.txt || true
# ---- one-file results bundle ----
python3 bundle_results_2026-09-24b.py "${SLURM_JOB_ID:-nojob}" || true
echo "Download: $(pwd)/results_2026-09-24b_${SLURM_JOB_ID:-nojob}.zip"
