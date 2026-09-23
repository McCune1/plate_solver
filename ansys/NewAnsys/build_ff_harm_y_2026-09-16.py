#!/usr/bin/env python3
"""Lever 9 (PAPER4_LEVERS_2026-09-16.md sec B.9): first harmonic
(ANTYPE,HARMIC) PLANE223 deck for the driven F-F n=0 admittance.

Clones PREP7 (mesh, mm-N-s-tonne, relative permittivity, TB,PIEZ,
NDIV_R/ZH/ZP) from the confirmed 1/12 OC modal parent (job 2491761).
Does NOT clone the modal /SOLU -- this is a new deck family.

Electrical BC for the DRIVEN problem (not OC, not scboth):
  inner electrodes (z=+-h): VOLT=0
  outer electrodes (z=+-(h+h1)): D,VOLT,1.0 on every outer node
    (fully electroded 1 V bus; CP is redundant when V is prescribed)

h1/2h=1/12 only. Do not attempt 1/8 or 1/5 until this family is
validated.

Closed-form (driven_ff_qv, job 2491920): pole at 747.986 rad/s =
119.047 Hz (elastic F-F), zero at 750.239 rad/s = 119.405 Hz (OC).
Fair FE modal pair on this mesh: scboth 118.811 Hz, OC 119.143 Hz.
The harmonic pole/zero are expected to sit near the FE modal pair
(mesh error), not bit-identical to closed-form.

PRE-REGISTERED (read the CSV, the NUMBER OF ERROR MESSAGES
ENCOUNTERED line, and RUN COMPLETED -- never the queue script's
own PASS/FAIL label):
  G_run    RUN COMPLETED, 0 MAPDL error lines.
  G_pole   |Y| peak in the dense window near the FE SC frequency
           (~118.8 Hz) / analytic 119.05 Hz, within ~0.5 Hz.
  G_zero   |Y| dip ABOVE that peak, near the FE OC frequency
           (~119.14 Hz) / analytic 119.40 Hz, within ~0.5 Hz.
  G_order  zero frequency > pole frequency (Butterworth-van Dyke).
  G_off    80 Hz and 140 Hz responses finite (no rigid-body blow-up
           at 80 Hz; 140 Hz is above the first pair).

Tiny DMPRAT=1e-4 keeps the pole finite. Compare frequencies, not
peak heights, to the undamped closed-form.

Do not put the MAPDL error-token string in this deck's comments.
"""
from __future__ import annotations

import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
OC_PARENT = os.path.join(HERE, "ansys_p4_epsrel_ff_lanb_oc_2026-09-16.inp")
OUT = os.path.join(HERE, "ansys_p4_harm_ff_y_h112_2026-09-16.inp")

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
]

FORBIDDEN = "*** ERROR ***"
SPLIT = "\n/SOLU\nANTYPE,MODAL\n"

SOLU = r"""
/SOLU
! Lever 9: harmonic (not modal). Driven 1 V on the outer bus, inners
! grounded. Tiny damping keeps the pole finite; compare pole/zero
! FREQUENCIES to driven_ff_qv, not peak heights.
ANTYPE,HARMIC
HROPT,FULL
HROUT,OFF
DMPRAT,1.0E-4
KBC,1
EQSLV,SPARSE

! Electrical BC -- DRIVEN F-F n=0 (imposed V, not OC, not scboth).
! Inners VOLT=0; outers every-node D,VOLT,1 (prescribed 1 V bus).
! No mechanical constraint (free-free). Do not start HARFRQ at 0
! (rigid-body singularity of F-F).
NSEL,S,LOC,Y,H
NSEL,A,LOC,Y,-H
D,ALL,VOLT,0
ALLSEL,ALL
NSEL,S,LOC,Y,H+H1
NSEL,A,LOC,Y,-H-H1
D,ALL,VOLT,1.0
ALLSEL,ALL

! LS1: dense window around the first n=0 pair.
! Analytic pole 119.047 Hz, zero 119.405 Hz; FE modal pair
! 118.811 / 119.143 Hz. 0.05 Hz steps.
HARFRQ,118.50,120.20
NSUBST,35
LSWRITE,1

! LS2: off-resonance below the pair (BvD: Ceff/C0 > 1).
HARFRQ,80.0,80.0
NSUBST,1
LSWRITE,2

! LS3: off-resonance above the pair (BvD: Ceff/C0 in (0,1)).
HARFRQ,140.0,140.0
NSUBST,1
LSWRITE,3

LSSOLVE,1,3
FINISH

/POST1
! Outer-bus current (harmonic reaction on VOLT = electrical current,
! Coupled-Field Guide). V=1 so Y = I. Real then imag via KIMG.
! Fill arrays then *VWRITE once -- same pattern as the modal
! _modes.txt dumps; scalar *VWRITE-in-a-loop is unreliable.
NROW=37
*DIM,LSARR,ARRAY,NROW
*DIM,SBARR,ARRAY,NROW
*DIM,FRARR,ARRAY,NROW
*DIM,IREARR,ARRAY,NROW
*DIM,IIMARR,ARRAY,NROW
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
    NSEL,S,LOC,Y,H+H1
    NSEL,A,LOC,Y,-H-H1
    IRE=0
    *GET,NCOUNT,NODE,0,COUNT
    *GET,NID,NODE,0,NUM,MIN
    *DO,II,1,NCOUNT
      *GET,QI,NODE,NID,RF,AMPS
      IRE=IRE+QI
      NID=NDNEXT(NID)
    *ENDDO
    SET,LSNOW,ISB,,1
    NSEL,S,LOC,Y,H+H1
    NSEL,A,LOC,Y,-H-H1
    IIM=0
    *GET,NCOUNT,NODE,0,COUNT
    *GET,NID,NODE,0,NUM,MIN
    *DO,II,1,NCOUNT
      *GET,QI,NODE,NID,RF,AMPS
      IIM=IIM+QI
      NID=NDNEXT(NID)
    *ENDDO
    ALLSEL,ALL
    LSARR(K)=LSNOW
    SBARR(K)=ISB
    FRARR(K)=FREQV
    IREARR(K)=IRE
    IIMARR(K)=IIM
    YARR(K)=SQRT(IRE*IRE+IIM*IIM)
  *ENDDO
*ENDDO
*CFOPEN,ansys_p4_harm_ff_y_h112_sweep,csv
*VWRITE
('ls,sb,freq_Hz,Ire,Iim,Yabs')
*VWRITE,LSARR(1),SBARR(1),FRARR(1),IREARR(1),IIMARR(1),YARR(1)
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
    with open(OC_PARENT, "r", encoding="utf-8", newline="") as f:
        parent = f.read().replace("\r\n", "\n").replace("\r", "\n")
    if SPLIT not in parent:
        raise SystemExit("FATAL: could not split OC parent on modal /SOLU")
    prep = parent.split(SPLIT, 1)[0]
    prep = prep.replace(
        "/TITLE, Paper 4 F-F piezo ring EPSREL LANB OPEN-CIRCUIT, p=0, h1/2h=1/12, mm-N-s-tonne (2026-09-16)",
        "/TITLE, Paper 4 F-F piezo ring EPSREL HARMIC driven Y, p=0, h1/2h=1/12, mm-N-s-tonne (2026-09-16)",
        1,
    )
    prep = prep.replace(
        "! ELECTRICAL BC: OPEN-CIRCUIT, Sec 18.175. VOLT=0 on the host/piezo\n"
        "! INTERFACE (z=+-h, inner electrodes grounded). Outer faces z=+-(h+h1)\n"
        "! are a conducting electrode (CP,VOLT, one parallel bus) left\n"
        "! unconstrained so Q=0 is the natural electrical condition. This is\n"
        "! the opposite of the parent SC deck. Do not re-open mesh/units/TB,PIEZ.\n",
        "! ELECTRICAL BC: DRIVEN 1 V (lever 9), not open-circuit and not\n"
        "! both-faces SC. VOLT=0 on the host/piezo INTERFACE (z=+-h). Outer\n"
        "! faces z=+-(h+h1) have D,VOLT,1.0 on every node (prescribed 1 V\n"
        "! bus). Do not re-open mesh/units/TB,PIEZ.\n",
        1,
    )
    prep = prep.replace(
        "! OPEN-CIRCUIT CLONE (2026-09-16): cloned from the cluster-confirmed\n"
        "! h1/2h=1/12 epsrel-fix F-F LANB deck (job 2490706, LESSONS_LEARNED.md\n"
        "! Sec 18.165). Mesh, mm-N-s-tonne, relative permittivity, FIX v6 TB,PIEZ,\n"
        "! LANB, free-free mechanical BC: UNCHANGED. Do not re-open those.\n"
        "!\n"
        "! THIS DECK changes ONLY the electrical BC, to match the analytic F-F\n"
        "! open-circuit model (Sec 18.175):\n"
        "!   inner electrodes (host/piezo interface, z=+-h): VOLT=0\n"
        "!   outer electrodes (z=+-(h+h1)): CP,VOLT as one parallel bus, unconstrained\n"
        "! Parent SC deck did the opposite (VOLT=0 on the OUTERS, interface free).\n"
        "!\n"
        "! Closed-form target: PiezoOutOfPlaneSolver.oc_ff_bisect n=0 h1/2h=1/12\n"
        "!   omega = 750.2394210586 rad/s = 119.4053 Hz\n"
        "! Parent SC FE mode 2 = 119.242 Hz (job 2490706). PASS: a bending mode\n"
        "! near 119.405 Hz, ABOVE 119.242 Hz (stiffening vs the same mesh), 0 Hz\n"
        "! rigid-body present, real O(0.1-1) bending discriminator, no\n"
        "! ill-conditioning warning. A result that reproduces 119.242 Hz is FAIL\n"
        "! -- the OC electrical BC did not take. Read this deck's own _modes.txt\n"
        "! / _discriminator.txt / _out.txt, never the .QUEUE_OK sentinel alone.\n",
        "! HARMONIC DRIVEN-Y CLONE (2026-09-16): PREP7 cloned from the\n"
        "! cluster-confirmed h1/2h=1/12 OC modal parent (job 2491761).\n"
        "! Mesh, mm-N-s-tonne, relative permittivity, FIX v6 TB,PIEZ,\n"
        "! free-free mechanical BC: UNCHANGED. Solution type is NEW\n"
        "! (ANTYPE,HARMIC, not LANB).\n"
        "!\n"
        "! THIS DECK's electrical BC is the DRIVEN problem:\n"
        "!   inner electrodes (host/piezo interface, z=+-h): VOLT=0\n"
        "!   outer electrodes (z=+-(h+h1)): D,VOLT,1.0 on every outer node\n"
        "! Analytic pole 119.047 Hz, zero 119.405 Hz (job 2491920). FE modal\n"
        "! pair 118.811 / 119.143 Hz. Read the sweep CSV, never .QUEUE_OK.\n",
        1,
    )
    banner = (
        "!\n"
        "! ============================================================\n"
        "! LEVER 9 HARMONIC ADMITTANCE (2026-09-16): PREP7 cloned from\n"
        "! the confirmed h1/2h=1/12 OC modal parent (job 2491761).\n"
        "! Mesh/units/TB,PIEZ unchanged. NEW solution type: ANTYPE,\n"
        "! HARMIC, driven 1 V on the outer bus. h1/2h=1/12 only.\n"
        "! Analytic pole 119.047 Hz, zero 119.405 Hz (job 2491920).\n"
        "! FE modal pair 118.811 / 119.143 Hz. Tiny DMPRAT=1e-4.\n"
        "! PASS: RUN COMPLETED, 0 MAPDL error lines, |Y| peak then\n"
        "! dip in 118.5-120.2 Hz with dip above peak. Read the CSV,\n"
        "! never the queue script's own PASS/FAIL label.\n"
        "! ============================================================\n"
    )
    # insert banner at the same FINISH /CLEAR slot the modal clones use
    marker = "FINISH\n/CLEAR,NOSTART\n/TITLE,"
    if marker not in prep:
        raise SystemExit("FATAL: missing FINISH/CLEAR/TITLE marker in PREP7")
    prep = prep.replace(marker, banner + marker, 1)
    text = prep + SOLU
    if not text.endswith("\n"):
        text += "\n"
    if "\r" in text:
        raise SystemExit("FATAL: CR leaked")
    if FORBIDDEN in text:
        raise SystemExit("FATAL: MAPDL error-token string present")
    if "ANTYPE,MODAL" in text:
        raise SystemExit("FATAL: leftover ANTYPE,MODAL")
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
