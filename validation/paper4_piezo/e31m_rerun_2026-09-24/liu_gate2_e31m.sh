#!/bin/bash
# Liu disk gate 2 at the new default bundle (LIU_DISK_KWARGS e31 = -4.1), output to its own dir
mkdir -p liu_e31m && cp probe_piezo_p4_liu_disk_gate2_2026-09-24.py liu_e31m/ && cd liu_e31m && PKG_PATH=../.. python3 -u probe_piezo_p4_liu_disk_gate2_2026-09-24.py > liu_gate2_e31m.log 2>&1
