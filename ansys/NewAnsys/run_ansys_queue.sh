#!/bin/bash
# run_ansys_queue.sh -- runs a LIST of .inp decks through MAPDL one at a
# time, in a single held allocation/session, so a single-seat ANSYS
# license is used efficiently without round-tripping through the SLURM
# queue once per deck. Meant to be invoked by submit_ansys_queue.sh (or
# directly on a login/interactive node with an ANSYS module already
# loaded), FROM Ansys/NewAnsys/ -- per this project's own established
# lesson, MAPDL jobs on this cluster only find their input file when run
# from THIS directory, not Ansys/ itself (a .inp sitting one level up
# fails with "ERROR: no input file specified" as its very first line, and
# a naive `cat` of old output can then display STALE results instead of
# a real failure -- see project memory ansys-newansys-execution-dir).
#
# WHY A MANIFEST, NOT AUTO-GLOB: Ansys/NewAnsys/ accumulates .inp decks
# from many unrelated threads over time (isoval calibration meshes,
# geometry-sweep verification, this-week's r200 adjudication, ...).
# Blindly globbing ansys_*.inp would re-run old, already-closed decks and
# burn license time for no reason. A manifest is an explicit, auditable
# list of what THIS queue run is for -- add new deck names as new
# geometries come up, don't touch old lines.
#
# RESUME-SKIP: same discipline as this project's own python capture
# probes (a "complete": true JSON sentinel). Here the sentinel is a
# "<deck>.QUEUE_OK" marker file, written ONLY after mapdl exits 0, and
# only trusted if it is NEWER than the .inp itself -- editing a deck and
# forgetting to clear the sentinel will not silently skip a stale result.
#
# FAILURE HANDLING: a failed deck is logged clearly and the queue MOVES ON
# to the next deck by default (one bad mesh/deck shouldn't cost the
# license-seat time of every deck behind it in the queue). Set
# STOP_ON_FAIL=1 to abort the whole queue on the first failure instead.
#
# Env vars:
#   MANIFEST      -- path to the deck list (default: ansys_queue_manifest.txt,
#                     one .inp filename per line, blank lines and #-comments
#                     ignored)
#   CPUS          -- mapdl -np worker count per deck (default: 8, matching
#                     this project's existing per-deck submit scripts)
#   STOP_ON_FAIL  -- if "1", abort the queue on the first mapdl failure
#                     instead of continuing to the next deck (default: 0)

set -u

MANIFEST="${MANIFEST:-ansys_queue_manifest.txt}"
CPUS="${CPUS:-8}"
STOP_ON_FAIL="${STOP_ON_FAIL:-0}"

if [ ! -f "$MANIFEST" ]; then
  echo "FATAL: manifest '$MANIFEST' not found in $(pwd)."
  echo "  Create it (one .inp filename per line) or set MANIFEST=<path>."
  exit 2
fi

echo "================================================================"
echo "  run_ansys_queue.sh -- $(date)"
echo "  cwd=$(pwd)  manifest=$MANIFEST  cpus/deck=$CPUS  stop_on_fail=$STOP_ON_FAIL"
echo "================================================================"

n_total=0
n_ok=0
n_skip=0
n_fail=0
declare -a SUMMARY

t_queue_start=$(date +%s)

while IFS= read -r line || [ -n "$line" ]; do
  # strip comments, CR (Windows manifests), and surrounding whitespace
  line="${line%$'\r'}"
  deck="$(printf '%s' "$line" | sed 's/#.*$//' | xargs)"
  [ -z "$deck" ] && continue
  n_total=$((n_total + 1))

  if [ ! -f "$deck" ]; then
    echo ""
    echo "---- [$n_total] $deck : MISSING -- not found in $(pwd), skipping ----"
    n_fail=$((n_fail + 1))
    SUMMARY+=("$deck : MISSING")
    [ "$STOP_ON_FAIL" = "1" ] && { echo "STOP_ON_FAIL=1, aborting queue."; break; }
    continue
  fi

  outlog="${deck%.inp}_out.txt"
  sentinel="${deck}.QUEUE_OK"

  if [ -f "$sentinel" ] && [ "$sentinel" -nt "$deck" ]; then
    echo ""
    echo "---- [$n_total] $deck : SKIP (already done -- $sentinel is newer than the deck) ----"
    n_skip=$((n_skip + 1))
    SUMMARY+=("$deck : SKIP (already done)")
    continue
  fi

  # PER-DECK ANSYS JOBNAME + STALE-LOCK CLEANUP (added 2026-08-26).
  # Previously every deck ran under MAPDL's default jobname "file", so all
  # decks shared file.db/file.rst/file.lock. If any deck terminates
  # abnormally -- a licence-queue timeout killed by SLURM is the live example,
  # jobs 2444733/2444734 -- the leftover lock blocks the NEXT deck with
  #   "*** ERROR *** Another ANSYS job with the same job name (file) is
  #    already running in this directory or the file.lock file has not been
  #    deleted from an abnormally terminated ANSYS run."
  # and MAPDL exits 100 having read nothing (jobs 2444927/2444928). A
  # per-deck jobname isolates the scratch, and we clear only THAT jobname's
  # lock -- never a blanket rm *.lock, which could clobber a concurrent job.
  # NOTE: a stale file.lock from 2026-08-14 is still in this directory; it
  # belongs to the old shared-jobname era and can simply be deleted.
  ansjob="$(basename "${deck%.inp}")"
  ansjob="${ansjob#ansys_}"
  ansjob="$(printf '%s' "$ansjob" | tr -c 'A-Za-z0-9' '_')"
  ansjob="${ansjob:0:24}"

  if [ -f "${ansjob}.lock" ]; then
    echo "     removing stale ${ansjob}.lock (previous run ended abnormally)"
    rm -f "${ansjob}.lock"
  fi

  echo ""
  echo "---- [$n_total] $deck : RUNNING (mapdl -np $CPUS -j $ansjob) -- $(date) ----"
  t0=$(date +%s)
  mapdl -b -smp -np "$CPUS" -j "$ansjob" -i "$deck" -o "$outlog"
  rc=$?
  dt=$(( $(date +%s) - t0 ))

  # EXIT 0 IS NOT ENOUGH (added 2026-08-26). MAPDL routinely exits 0 having
  # printed "*** ERROR ***" and skipped work -- an out-of-memory solve on a
  # late pass of a multi-pass deck is the live example. Without this check the
  # deck gets a .QUEUE_OK sentinel, and every later resubmit SKIPS it, so a
  # silently truncated result becomes permanent. Set IGNORE_DECK_ERRORS=1 to
  # override for a deck whose errors are known-benign.
  # grep -c prints 0 AND exits 1 when nothing matches, so the old
  # `|| echo 0` produced "0\n0" and an "integer expression expected"
  # warning (harmless, but noise). Fixed 2026-09-23.
  n_err=$(grep -c "\*\*\* ERROR \*\*\*" "$outlog" 2>/dev/null)
  n_err=${n_err:-0}
  deck_errors_ok=1
  if [ "$n_err" -gt 0 ] && [ "${IGNORE_DECK_ERRORS:-0}" != "1" ]; then
    deck_errors_ok=0
  fi

  if [ $rc -eq 0 ] && [ "$deck_errors_ok" = "1" ] \
     && ! grep -q "ERROR: no input file specified" "$outlog" 2>/dev/null; then
    touch "$sentinel"
    echo "---- [$n_total] $deck : OK (${dt}s), wrote $outlog, sentinel=$sentinel ----"
    n_ok=$((n_ok + 1))
    SUMMARY+=("$deck : OK (${dt}s)")
  else
    if [ "$deck_errors_ok" = "0" ]; then
      echo "---- [$n_total] $deck : FAILED -- MAPDL exited $rc but printed"
      echo "     $n_err '*** ERROR ***' line(s). NO sentinel written, so this"
      echo "     deck will re-run on the next submit rather than being skipped"
      echo "     with a truncated result. First errors:"
      grep -n -m5 "\*\*\* ERROR \*\*\*" "$outlog" 2>/dev/null | sed 's/^/       /'
    fi
    echo "---- [$n_total] $deck : FAILED (exit=$rc, ${dt}s) -- see $outlog ----"
    echo "     last 15 lines of $outlog:"
    tail -n 15 "$outlog" 2>/dev/null | sed 's/^/     /'
    n_fail=$((n_fail + 1))
    SUMMARY+=("$deck : FAILED (exit=$rc)")
    [ "$STOP_ON_FAIL" = "1" ] && { echo "STOP_ON_FAIL=1, aborting queue."; break; }
  fi
done < "$MANIFEST"

dt_queue=$(( $(date +%s) - t_queue_start ))

echo ""
echo "================================================================"
echo "  QUEUE SUMMARY -- $(date)  (total wall: $((dt_queue / 60)) min)"
echo "================================================================"
for line in "${SUMMARY[@]}"; do
  echo "  $line"
done
echo ""
echo "  total=$n_total  ok=$n_ok  skip=$n_skip  fail=$n_fail"

if [ "$n_fail" -gt 0 ]; then
  exit 1
fi
exit 0
