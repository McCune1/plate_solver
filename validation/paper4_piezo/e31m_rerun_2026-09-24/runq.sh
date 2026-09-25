#!/bin/bash
# usage: runq.sh name probe1 probe2 ...
q=$1; shift
for p in "$@"; do
  b=${p%.py}
  echo "=== $(date +%T) start $p" >> queue_$q.log
  PKG_PATH=.. timeout 7200 python3 -u $p > ${b}.log 2>&1
  echo "=== $(date +%T) end $p rc=$? $(grep -h SENTINEL ${b}.log | tail -1)" >> queue_$q.log
done
