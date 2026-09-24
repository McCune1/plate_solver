#!/bin/bash
#SBATCH --job-name=p4_3d_nge1
#SBATCH --output=ansys_queue_piezo_p4_3d_nge1_%j.out
#SBATCH --error=ansys_queue_piezo_p4_3d_nge1_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=06:00:00
#SBATCH --mem=64G
#SBATCH --requeue
#SBATCH --partition=requeue
#SBATCH --chdir=/home/ghmkfh/PythonMill/Plate_Solver_Package/Ansys/NewAnsys

# Paper 4 3-D FE check of coupled n >= 1 roots. ONE ANSYS seat: this
# script calls run_ansys_queue.sh. Do not also sbatch the decks separately.
#
# Cloned from submit_ansys_queue_piezo_p6_sector_2026-09-22.sh (job
# 2528994). Absolute --chdir. Per-deck -j jobname, -np, the failure-token
# scan, the no-input-file guard and the _out.txt copy all live in
# run_ansys_queue.sh. CPUS is ${SLURM_CPUS_PER_TASK:-8}.
#
# Budget: ten SOLID226 decks, the largest ~8.3e5 DOF (Paper 6's 1.6e6-DOF
# fine deck took 60 s wall and 12 GB scratch). --time=06:00:00 is generous
# on purpose (requeue preemption restarts only unfinished decks: finished
# ones have a .QUEUE_OK newer than the .inp and are skipped).
#
# PUSH the ten .inp files, this script, and
#   ansys_queue_manifest_piezo_p4_3d_nge1_2026-09-23.txt
#   score_p4_3d_nge1_2026-09-23.py
#   targets_p4_3d_nge1_2026-09-23.json
# to /home/ghmkfh/PythonMill/Plate_Solver_Package/Ansys/NewAnsys/
# then from that directory:
#   sbatch submit_ansys_queue_piezo_p4_3d_nge1_2026-09-23.sh
#
# After the run, score from the transcripts (python3 score_p4_3d_nge1_2026-09-23.py).
# Do not treat .QUEUE_OK as the result.

echo "Job $SLURM_JOB_ID on $(hostname)"
echo "Start: $(date)"
if [ ! -f run_ansys_queue.sh ]; then
  echo "FATAL: must be Ansys/NewAnsys/"
  exit 2
fi
# ---- ANSYS environment (rewritten after job 2530336) --------------------
# Job 2530336 died in 1 s: every deck got "mapdl: command not found"
# (exit 127). The old one-liner
#   module load ansys/2024R2 2>/dev/null || module load ansys 2>/dev/null
# hid the reason: `module` itself was not defined in the batch shell, or the
# module failed to load, and both errors went to /dev/null. Job 2528994
# ran the same line on the same node the night before and worked. The
# likeliest cause is the submitting shell: sbatch copies its environment,
# including the exported `module` function. Now we initialise modules
# ourselves, print every error, fall back to a versioned launcher, and
# stop BEFORE the queue if there is still no mapdl. (Login-node check the
# same evening: `module load ansys/2024R2` -> Lmod "unknown module";
# `module load ansys` -> no error, still no mapdl on PATH.)
if ! type module >/dev/null 2>&1; then
  for init in /etc/profile.d/lmod.sh /etc/profile.d/modules.sh \
              /usr/share/lmod/lmod/init/bash /usr/share/Modules/init/bash; do
    if [ -f "$init" ]; then
      echo "module not defined in batch shell; sourcing $init"
      . "$init"
      break
    fi
  done
fi
# 2026-09-23 (Mill module change): ansys/2024R2 is gone. `module spider
# ansys` lists 2023R2, 2024_EDT, 2024R1, 2025R1, 2026R1. The default,
# 2026R1, is a CONTAINER: its modulefile only defines a shell function
# `ansys` -> /share/apps/common/ansys/2026R1/bin/ansys_launcher_logic.sh,
# with no mapdl on PATH and no AWP_ROOT variables.
# Order tried here (override: ANSYS_MODULES="ansys/2024R1" sbatch ...):
#   ansys/2025R1, ansys/2024R1  -- nearest to the 2024R2 that ran every
#                                  earlier Paper 4-6 deck; used only if the
#                                  module puts a real launcher on PATH
#   ansys/2026R1                -- container, last resort, via a wrapper
#                                  that runs `<launcher> mapdl <args>`
# MAPDL_EXE=/full/path/to/launcher still overrides everything.
ANSYS_MODULES="${ANSYS_MODULES:-ansys/2025R1 ansys/2024R1 ansys/2026R1}"
find_mapdl() {
  if [ -n "${MAPDL_EXE:-}" ] && [ -x "$MAPDL_EXE" ]; then echo "$MAPDL_EXE"; return; fi
  command -v mapdl 2>/dev/null && return
  local d f v
  for d in $(printf '%s' "$PATH" | tr ':' ' '); do
    for f in "$d"/ansys[0-9][0-9][0-9]; do
      [ -x "$f" ] && { echo "$f"; return; }
    done
  done
  for v in $(env | sed -n 's/^\(AWP_ROOT[0-9]*\|ANSYS[0-9]*_DIR\)=.*/\1/p' | sort -r); do
    d="${!v}"
    for f in "$d"/ansys/bin/mapdl "$d"/ansys/bin/ansys[0-9][0-9][0-9] \
             "$d"/bin/mapdl "$d"/bin/ansys[0-9][0-9][0-9]; do
      [ -x "$f" ] && { echo "$f"; return; }
    done
  done
}
shim="$PWD/.mapdl_shim_${SLURM_JOB_ID:-local}"
exe="" ; used_module="" ; via=""
exe="$(find_mapdl | head -1)"
if [ -z "$exe" ] && type module >/dev/null 2>&1; then
  for m in $ANSYS_MODULES; do
    module unload ansys >/dev/null 2>&1
    echo "trying module $m"
    module load "$m" || continue
    exe="$(find_mapdl | head -1)"
    if [ -n "$exe" ]; then used_module="$m"; via="native launcher"; break; fi
    launcher="$(type ansys 2>/dev/null | grep -o '/[^" ]*ansys_launcher_logic\.sh' | head -1)"
    if [ -n "$launcher" ] && [ -x "$launcher" ]; then
      mkdir -p "$shim"
      printf '#!/bin/bash\nexec "%s" mapdl "$@"\n' "$launcher" > "$shim/mapdl_container"
      chmod +x "$shim/mapdl_container"
      exe="$shim/mapdl_container"; used_module="$m"; via="container wrapper ($launcher mapdl ...)"
      break
    fi
    echo "  $m gives no MAPDL launcher"
  done
elif [ -z "$exe" ]; then
  echo "WARNING: no module command available even after init scripts"
fi
if [ -n "$exe" ] && [ "$exe" != "$(command -v mapdl 2>/dev/null)" ]; then
  mkdir -p "$shim" && ln -sf "$exe" "$shim/mapdl"
  export PATH="$shim:$PATH"
fi
if ! command -v mapdl >/dev/null 2>&1; then
  echo "FATAL: no MAPDL launcher found. No deck was started."
  echo "  tried modules: $ANSYS_MODULES"
  echo "  PATH=$PATH"
  if type module >/dev/null 2>&1; then
    for m in $ANSYS_MODULES; do
      echo "  module show $m:"
      module show "$m" 2>&1 | sed 's/^/    /'
    done
  fi
  echo "  Fix: resubmit with MAPDL_EXE=/path/to/launcher or ANSYS_MODULES=..."
  exit 3
fi
type module >/dev/null 2>&1 && module list 2>&1 | sed 's/^/  [module list] /'
echo "ANSYS module: ${used_module:-(none; launcher already on PATH or MAPDL_EXE)}"
echo "MAPDL launcher: $exe ${via:+[$via]}"
if [ "$used_module" = "ansys/2026R1" ]; then
  echo "NOTE: container route. Read deck 1's _out.txt first: if it shows"
  echo "  the MAPDL banner and RUN COMPLETED the wrapper works; if it is empty"
  echo "  or shows a launcher usage message, the container's MAPDL command"
  echo "  is not 'mapdl' -- see /share/apps/common/ansys/2026R1/bin/ansys_launcher_logic.sh."
fi
export CPUS=${SLURM_CPUS_PER_TASK:-8}
export MANIFEST=ansys_queue_manifest_piezo_p4_3d_nge1_2026-09-23.txt
export STOP_ON_FAIL=0
bash run_ansys_queue.sh
rc=$?
rm -rf ".mapdl_shim_${SLURM_JOB_ID:-local}" 2>/dev/null
if command -v python3 >/dev/null 2>&1 && [ -f score_p4_3d_nge1_2026-09-23.py ]; then
  echo ""
  echo "==== scorer (convenience; re-run it yourself on the synced files) ===="
  python3 score_p4_3d_nge1_2026-09-23.py || true
fi
echo "End: $(date)"
exit $rc
