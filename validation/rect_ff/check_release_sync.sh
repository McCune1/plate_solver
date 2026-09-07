#!/bin/bash
# check_release_sync.sh -- 2026-09-07
#
# Standalone, dependency-free guard against the defect Paper 2 review pass #3
# found: github_repo/plate_solver/ silently falling behind the deployed
# Plate_Solver_Package/plate_solver/ tree. It needs nothing but coreutils, so
# it runs on the workstation VM (which has no scipy and cannot import
# plate_solver) as well as anywhere else.
#
# Run from the project root:
#     bash Plate_Solver_Package/check_release_sync.sh
# or with explicit roots:
#     bash check_release_sync.sh <deployed_pkg_root> <release_root>
#
# Exit 0 = in sync. Exit 1 = drift; sync release <- deployed (never the other
# way), then rebuild PLATE_SOLVER_REPO_DIGEST.md before doing anything else.
set -uo pipefail

DEP="${1:-}"
REL="${2:-}"
if [ -z "$DEP" ]; then
  here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
  DEP="$here"
  REL="$(dirname "$here")/github_repo"
fi

dep_pkg="$DEP/plate_solver"
rel_pkg="$REL/plate_solver"
for d in "$dep_pkg" "$rel_pkg"; do
  [ -d "$d" ] || { echo "ERROR: no such directory: $d"; exit 2; }
done

echo "deployed: $dep_pkg"
echo "release : $rel_pkg"
echo

drift=0
n=0
for f in "$dep_pkg"/*.py; do
  b="$(basename "$f")"
  n=$((n+1))
  if [ ! -f "$rel_pkg/$b" ]; then
    echo "MISSING  $b  (not in release)"; drift=$((drift+1)); continue
  fi
  a="$(md5sum < "$f" | cut -d' ' -f1)"
  c="$(md5sum < "$rel_pkg/$b" | cut -d' ' -f1)"
  if [ "$a" != "$c" ]; then
    echo "DRIFT    $b  deployed ${a:0:12}  release ${c:0:12}"; drift=$((drift+1))
  fi
done
for f in "$rel_pkg"/*.py; do
  b="$(basename "$f")"
  [ -f "$dep_pkg/$b" ] || { echo "EXTRA    $b  (in release only)"; drift=$((drift+1)); }
done

echo
echo "compared $n deployed modules"
if [ "$drift" -eq 0 ]; then
  echo "RELEASE_SYNC: OK"
  exit 0
fi
echo "RELEASE_SYNC: DRIFT ($drift file(s))"
echo "Fix:  cp \$DEP/plate_solver/<file> \$REL/plate_solver/<file>"
echo "Then: python github_repo/scripts/build_repo_digest.py"
exit 1
