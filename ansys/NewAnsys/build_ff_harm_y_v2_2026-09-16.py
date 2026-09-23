#!/usr/bin/env python3
"""Lever 9 v2: same harmonic family as job 2494669, dump fixed.

Job 2494669: MAPDL RUN COMPLETED, 0 error messages, but I=1.28e-27
at every frequency. 11992 warnings `Unknown parameter name= QI`.
Cause: this deck's PLANE223 KEYOPT(1)=1001 (structural-piezoelectric)
uses force label CHRG for the VOLT DOF, not AMPS. AMPS is KEYOPT(1)=101
(current-based). Coupled-Field Guide: "KEYOPT(1)=1001, the force label
for the VOLT degree of freedom is CHRG." Also write reactions to the
rst (OUTRES,RSOL) and extract from the CP-master node only, so there
is one reaction to read, not a node loop.

Do not retune mesh/units/TBDATA/HARFRQ. Isolated dump+OUTRES+CP-master
change. No MAPDL error-token string in comments.
"""
from __future__ import annotations

import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
PARENT = os.path.join(HERE, "ansys_p4_harm_ff_y_h112_2026-09-16.inp")
OUT = os.path.join(HERE, "ansys_p4_harm_ff_y_v2_h112_2026-09-16.inp")
FORBIDDEN = "*** ERROR ***"

EXISTING = [
    "ansys_p4_epsrel_ff_e0_2026-09-16.inp",
    "ansys_p4_epsrel_ff_lanb_2026-09-16.inp",
    "ansys_p4_epsrel_ff_subsp_2026-09-16.inp",
    "ansys_p4_epsrel_cc_lanb_2026-09-16.inp",
    "ansys_p4_epsrel_ff_lanb_h18_2026-09-16.inp",
    "ansys_p4_epsrel_ff_lanb_h15_2026-09-16.inp",
    "ansys_p4_epsrel_cc_lanb_rmid_2026-09-16.inp",
    "ansys_p4_epsrel_ff_lanb_oc_2026-09-16.inp",
    "ansys_p4_epsrel_ff_lanb_scboth_2026-09-16.inp",
    "ansys_p4_epsrel_ff_lanb_h18_oc_2026-09-16.inp",
    "ansys_p4_epsrel_ff_lanb_h18_scboth_2026-09-16.inp",
    "ansys_p4_epsrel_ff_lanb_h15_oc_2026-09-16.inp",
    "ansys_p4_epsrel_ff_lanb_h15_scboth_2026-09-16.inp",
    "ansys_p4_zp4_ff_oc_h112_2026-09-16.inp",
    "ansys_p4_zp4_ff_sc_h112_2026-09-16.inp",
    "ansys_p4_zp4_ff_oc_h18_2026-09-16.inp",
    "ansys_p4_zp4_ff_sc_h18_2026-09-16.inp",
    "ansys_p4_zp4_ff_oc_h15_2026-09-16.inp",
    "ansys_p4_zp4_ff_sc_h15_2026-09-16.inp",
    "ansys_p4_harm_ff_y_h112_2026-09-16.inp",
]

SOLU = r"""
/SOLU
! Lever 9 v2. Same HARMIC / mesh / HARFRQ as job 2494669. Isolated
! dump fix: KEYOPT(1)=1001 force label is CHRG not AMPS; OUTRES,RSOL
! so reactions are on the rst; CP outer bus and D the master only so
! one node carries the total electrode charge. Do not retune mesh,
! units, TB,PIEZ, DMPRAT, or the frequency window.
ANTYPE,HARMIC
HROPT,FULL
HROUT,OFF
OUTRES,NSOL,ALL
OUTRES,RSOL,ALL
DMPRAT,1.0E-4
KBC,1
EQSLV,SPARSE

NSEL,S,LOC,Y,H
NSEL,A,LOC,Y,-H
D,ALL,VOLT,0
ALLSEL,ALL
NSEL,S,LOC,Y,H+H1
NSEL,A,LOC,Y,-H-H1
CP,1,VOLT,ALL
*GET,NMAST,NODE,0,NUM,MIN
D,NMAST,VOLT,1.0
ALLSEL,ALL

HARFRQ,118.50,120.20
NSUBST,35
LSWRITE,1

HARFRQ,80.0,80.0
NSUBST,1
LSWRITE,2

HARFRQ,140.0,140.0
NSUBST,1
LSWRITE,3

LSSOLVE,1,3
FINISH

/POST1
! KEYOPT(1)=1001 -> VOLT reaction label is CHRG (charge). For V=1,
! Ceff=Q and |Y|=omega*|Q|. Pole/zero frequencies of Q and Y coincide
! for omega!=0. Read the CSV, never the queue script's own label.
! One PRRSOL at the first dense substep so a human can see the label
! if *GET ever fails again.
SET,1,1,,0
NSEL,S,NODE,,NMAST
PRRSOL,CHRG
ALLSEL,ALL

NROW=37
*DIM,LSARR,ARRAY,NROW
*DIM,SBARR,ARRAY,NROW
*DIM,FRARR,ARRAY,NROW
*DIM,QREARR,ARRAY,NROW
*DIM,QIMARR,ARRAY,NROW
*DIM,YARR,ARRAY,NROW
*DIM,LSLIST,ARRAY,3
*DIM,NSBLIST,ARRAY,3
LSLIST(1)=1
LSLIST(2)=2
LSLIST(3)=3
NSBLIST(1)=35
NSBLIST(2)=1
NSBLIST(3)=1
K=0
*DO,IL,1,3
  LSNOW=LSLIST(IL)
  NSB=NSBLIST(IL)
  *DO,ISB,1,NSB
    K=K+1
    SET,LSNOW,ISB,,0
    *GET,FREQV,ACTIVE,0,SET,FREQ
    *GET,QREAC,NODE,NMAST,RF,CHRG
    SET,LSNOW,ISB,,1
    *GET,QIMAG,NODE,NMAST,RF,CHRG
    OMEGA=2*3.141592653589793*FREQV
    QABS=SQRT(QREAC*QREAC+QIMAG*QIMAG)
    LSARR(K)=LSNOW
    SBARR(K)=ISB
    FRARR(K)=FREQV
    QREARR(K)=QREAC
    QIMARR(K)=QIMAG
    YARR(K)=OMEGA*QABS
  *ENDDO
*ENDDO
*CFOPEN,ansys_p4_harm_ff_y_v2_h112_sweep,csv
*VWRITE
('ls,sb,freq_Hz,Qre,Qim,Yabs')
*VWRITE,LSARR(1),SBARR(1),FRARR(1),QREARR(1),QIMARR(1),YARR(1)
(F4.0,',',F4.0,',',E16.8,',',E16.8,',',E16.8,',',E16.8)
*CFCLOS
FINISH
"""


def jobname(fname: str) -> str:
    stem = os.path.basename(fname)
    if stem.endswith(".inp"):
        stem = stem[:-4]
    if stem.startswith("ansys_"):
        stem = stem[6:]
    stem = re.sub(r"[^A-Za-z0-9]", "_", stem)
    return stem[:24]


def main() -> None:
    with open(PARENT, "r", encoding="utf-8", newline="") as f:
        parent = f.read().replace("\r\n", "\n").replace("\r", "\n")
    marker = "\n/SOLU\n"
    if marker not in parent:
        raise SystemExit("FATAL: parent has no /SOLU")
    prep = parent.split(marker, 1)[0]
    prep = prep.replace(
        "/TITLE, Paper 4 F-F piezo ring EPSREL HARMIC driven Y, p=0, h1/2h=1/12, mm-N-s-tonne (2026-09-16)",
        "/TITLE, Paper 4 F-F piezo ring EPSREL HARMIC driven Y v2 CHRG dump, p=0, h1/2h=1/12 (2026-09-16)",
        1,
    )
    banner = (
        "!\n"
        "! ============================================================\n"
        "! LEVER 9 v2 (2026-09-16): dump fix of job 2494669. Same mesh,\n"
        "! HARFRQ, DMPRAT, KEYOPT(1)=1001. Isolated changes: force label\n"
        "! CHRG (not AMPS), OUTRES,RSOL, CP-master D,VOLT,1. PASS: RUN\n"
        "! COMPLETED, 0 MAPDL error lines, |Y| peak then dip in\n"
        "! 118.5-120.2 Hz with dip above peak, 80 Hz and 140 Hz finite\n"
        "! and not ~1e-27. Read the CSV and PRRSOL block, never the\n"
        "! queue script's own PASS/FAIL label.\n"
        "! ============================================================\n"
    )
    slot = "FINISH\n/CLEAR,NOSTART\n/TITLE,"
    if slot not in prep:
        raise SystemExit("FATAL: missing FINISH/CLEAR/TITLE")
    # keep the existing lever-9 banner; prepend the v2 note at the slot
    prep = prep.replace(slot, banner + slot, 1)
    text = prep + SOLU
    if not text.endswith("\n"):
        text += "\n"
    if "\r" in text:
        raise SystemExit("FATAL: CR leaked")
    if FORBIDDEN in text:
        raise SystemExit("FATAL: MAPDL error-token string present")
    if "RF,AMPS" in text:
        raise SystemExit("FATAL: leftover RF,AMPS")
    if "RF,CHRG" not in text:
        raise SystemExit("FATAL: missing RF,CHRG")
    if "ANTYPE,HARMIC" not in text:
        raise SystemExit("FATAL: missing ANTYPE,HARMIC")

    all_names = EXISTING + [os.path.basename(OUT)]
    seen = {}
    print("jobname uniqueness:")
    for fname in all_names:
        jn = jobname(fname)
        if jn in seen:
            raise SystemExit("FATAL jobname collision %r: %s vs %s"
                             % (jn, seen[jn], fname))
        seen[jn] = fname
        print("  jobname %r <- %s" % (jn, fname))

    with open(OUT, "w", encoding="utf-8", newline="\n") as f:
        f.write(text)
    print("wrote %s (%d lines)" % (OUT, text.count("\n")))


if __name__ == "__main__":
    main()
