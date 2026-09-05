#!/bin/bash
# Generate 8 C-class n80 IP MAC decks into Ansys/NewAnsys/.
# Does NOT run MAPDL. Login-node python3 on Mill is 3.6 -- load 3.12 first.
#
# Prefix ip_mac_gs_*_n80 does NOT overwrite Rank 14's 57-mode dumps.
# FB_MODE_LO=4 FB_MODE_HI=83 -> 80 elastic modes. FB_NMODES=90 for the solve.

set -u
export PKG_PATH="${PKG_PATH:-/home/ghmkfh/PythonMill/Plate_Solver_Package}"
cd "$PKG_PATH" || exit 2
if [ ! -f generate_ansys_ip_geomsweep_mac_fullblock.py ]; then
  echo "FATAL: generate_ansys_ip_geomsweep_mac_fullblock.py not in $(pwd)"
  exit 2
fi

# Mill login python3 is 3.6. The generator is 3.6-safe, but prefer the
# package venv when it exists so this matches every other probe.
if [ -f /usr/share/Modules/init/bash ] || [ -f /etc/profile.d/modules.sh ]; then
  # shellcheck disable=SC1091
  source /usr/share/Modules/init/bash 2>/dev/null || true
  module load python/3.12.1 2>/dev/null || true
fi
if [ -f "$PKG_PATH/venv/bin/activate" ]; then
  # shellcheck disable=SC1091
  source "$PKG_PATH/venv/bin/activate"
fi
echo "python: $(command -v python3)  $($(command -v python3) --version 2>&1)"

export FB_MODE_LO=4
export FB_MODE_HI=83
export FB_NMODES=90
export GS_NU=0.35

ANSYS_DIR="$PKG_PATH/Ansys/NewAnsys"
mkdir -p "$ANSYS_DIR"
MANIFEST="$ANSYS_DIR/ansys_queue_manifest_cclass_n80.txt"
TMP_MANIFEST="${MANIFEST}.tmp"
{
  echo "# C-class n80 IP MAC extracts. Generated $(date -Iseconds 2>/dev/null || date)."
  echo "# 80 elastic modes (4..83). Do not mix with the default Rank 14 dumps."
} > "$TMP_MANIFEST"

n_ok=0
for RATIO_TAG in r150 r167 r200 r250; do
  for ANGLE_TAG in a125 a150; do
    export RATIO_TAG ANGLE_TAG
    export FB_PREFIX="ip_mac_gs_${RATIO_TAG}_${ANGLE_TAG}_n80"
    export FB_DECK="ansys_ip_gs_${RATIO_TAG}_${ANGLE_TAG}_n80_mac_fullblock.inp"
    export FB_OUT="$FB_DECK"
    export FB_JOB="gs${RATIO_TAG}${ANGLE_TAG}n80"
    python3 -u generate_ansys_ip_geomsweep_mac_fullblock.py || {
      echo "FATAL: generator failed for $RATIO_TAG/$ANGLE_TAG"
      rm -f "$TMP_MANIFEST"
      exit 2
    }
    if [ ! -f "$ANSYS_DIR/$FB_DECK" ]; then
      echo "FATAL: $ANSYS_DIR/$FB_DECK was not written"
      rm -f "$TMP_MANIFEST"
      exit 2
    fi
    echo "$FB_DECK" >> "$TMP_MANIFEST"
    n_ok=$((n_ok + 1))
  done
done

if [ "$n_ok" -ne 8 ]; then
  echo "FATAL: expected 8 decks, wrote $n_ok. Not replacing manifest."
  rm -f "$TMP_MANIFEST"
  exit 2
fi
mv "$TMP_MANIFEST" "$MANIFEST"
echo "wrote 8 decks and $MANIFEST"
echo "ls:"
ls -1 "$ANSYS_DIR"/ansys_ip_gs_*_n80_mac_fullblock.inp
echo
echo "Next (one licence seat -- do NOT sbatch eight ANSYS jobs):"
echo "  cd $PKG_PATH/Ansys/NewAnsys"
echo "  MANIFEST=ansys_queue_manifest_cclass_n80.txt sbatch submit_ansys_queue.sh"
echo "Do NOT submit submit_ip_cclass_n80_mac.sh until each deck has a"
echo ".QUEUE_OK sentinel AND ls Ansys/NewAnsys/ip_mac_gs_<key>_n80_m*_mesh96.txt"
echo "shows 80 files. Then, CPU jobs (can overlap):"
echo "  RATIO_TAG=r150 ANGLE_TAG=a125 sbatch submit_ip_cclass_n80_mac.sh"
echo "  (repeat for r150/r167/r200/r250 x a125/a150)."
