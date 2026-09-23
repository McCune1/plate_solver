#!/usr/bin/env python3
"""Build Paper 5 Option A F-F FE decks: elastic/e0, both-faces SC, OC.

Clone of the confirmed Paper 4 F-F epsrel LANB family (jobs 2491761 /
2491831): mm-N-s-tonne, PLANE223 KEYOPT(1)=1001 KEYOPT(3)=1, NDIV_R=40,
NDIV_Z=6, LANB, TB,PIEZ FIX v6 Y-poled layout, relative permittivity.
NOT a retune of that family. Geometry change only: one PZT-4 rectangle
|z|<=H=10 mm, no steel host, no piezo skins.

Closed-form (PiezoMonolithicOutOfPlaneSolver, dps=60, LESSONS Sec 18.201):
  elastic F-F n=0  462.3976497355817 rad/s = 73.577 Hz
  SC      F-F n=0  473.2496088764243 rad/s = 75.304 Hz
  OC flexure = SC at this ansatz (PAPER5_DERIVATION.md Sec 4).

Output is LF. Jobnames after run_ansys_queue.sh 24-char truncation are
asserted unique against existing ansys_*.inp stems.
"""
from __future__ import annotations

import glob
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
DATE = "2026-09-19"

EL_HZ = 73.577
SC_HZ = 75.304


def jobname(fname: str) -> str:
    stem = os.path.basename(fname)
    if stem.endswith(".inp"):
        stem = stem[:-4]
    if stem.startswith("ansys_"):
        stem = stem[6:]
    stem = re.sub(r"[^A-Za-z0-9]", "_", stem)
    return stem[:24]


def existing_jobnames():
    out = {}
    for path in glob.glob(os.path.join(HERE, "ansys_*.inp")):
        out[jobname(path)] = os.path.basename(path)
    return out


HEADER = r"""! ansys_p5_epsrel_ff_lanb_{tag}_2026-09-19.inp
! Paper 5 Option A -- homogeneous thickness-poled PZT-4 annulus, F-F, p=0.
! Clone of Paper 4 F-F epsrel LANB (jobs 2491761 / 2491831): same units,
! element, mesh density, TB,PIEZ layout, LANB. Do not retune those.
! Geometry: ONE ceramic |z|<=H, electrodes on the two faces only.
! Not a flag on PiezoOutOfPlaneSolver. Not Parashar d15.
!
! UNITS: mm-N-s-tonne, real volts. Same conversion as
! ansys_p4_epsrel_ff_lanb_scboth_2026-09-16.inp (do not re-derive).
! Frequencies come out in Hz.
!
! Closed-form (plate_solver.piezo_monolithic, LESSONS Sec 18.201):
!   elastic 73.577 Hz   SC 75.304 Hz   OC = SC (Phase 0 theorem)
!
! F-F: unconstrained. Expect a rigid-body mode near 0 Hz, then bending.
! Discriminator at r=RO, z=+-H: bending has UY same sign, UX opposite.
!
! Do not put the MAPDL failure token in this file (queue grep false-fail).
! ============================================================
FINISH
/CLEAR,NOSTART
/TITLE, Paper 5 F-F monolithic PZT-4 {title}, p=0, mm-N-s-tonne (2026-09-19)

RI     = 100.0
RO     = 600.0
H      = 10.0
NDIV_R = 40
NDIV_Z = 6
NMODES = 20

C11E = 1.320E5
C12E = 7.100E4
C13E = 7.300E4
C33E = 1.150E5
C44E = 7.300E4
C66E = (C11E-C12E)/2
E31  = 4.1E-3
E33  = 1.41E-2
E15  = 1.05E-2
EPS0   = 8.854E-12
X11ABS = 7.124E-9
X33ABS = 5.841E-9
X11    = X11ABS/EPS0
X33    = X33ABS/EPS0
RHO_PZT = 7.500E-9

/PREP7
ET,1,PLANE223
KEYOPT,1,1,1001
KEYOPT,1,3,1

TB,ANEL,1,1,,0
TBDATA,1,C11E,C13E,C12E,0,0,0
TBDATA,7,C33E,C13E,0,0,0
TBDATA,12,C11E,0,0,0
TBDATA,16,C44E,0,0
TBDATA,19,C66E,0
TBDATA,21,C44E
TB,PIEZ,1
{piez}
MP,PERX,1,X11
MP,PERY,1,X33
MP,PERZ,1,X11
MP,DENS,1,RHO_PZT

RECTNG,RI,RO,-H,H
ASEL,S,AREA,,1
AATT,1,,1
LSLA,S
LSEL,R,LOC,X,RI
LESIZE,ALL,,,NDIV_Z,,,,,1
LSLA,S
LSEL,R,LOC,X,RO
LESIZE,ALL,,,NDIV_Z,,,,,1
LSLA,S
LSEL,R,LOC,Y,-H
LESIZE,ALL,,,NDIV_R,,,,,1
LSLA,S
LSEL,R,LOC,Y,H
LESIZE,ALL,,,NDIV_R,,,,,1
ALLSEL,ALL
ASEL,S,AREA,,1
MSHAPE,0,2D
MSHKEY,1
AMESH,ALL
ALLSEL,ALL
FINISH

/SOLU
ANTYPE,MODAL
MODOPT,LANB,NMODES,-1.0,,,ON
MXPAND,NMODES,,,NO
{ebc}
SOLVE
FINISH

/POST1
*DIM,FR1,ARRAY,NMODES
*DIM,IDXARR,ARRAY,NMODES
*DO,I,1,NMODES
  *GET,FR1(I),MODE,I,FREQ
  IDXARR(I)=I
*ENDDO
*CFOPEN,ansys_p5_epsrel_ff_lanb_{tag}_modes,txt
*VWRITE
('mode f[Hz] tag=p5_ff_{tag} target={target}Hz')
*VWRITE,IDXARR(1),FR1(1)
(F4.0, 3X, E16.8)
*CFCLOS
NSEL,S,LOC,X,RO
NSEL,R,LOC,Y,H
*GET,NTOP,NODE,0,NUM,MIN
NSEL,S,LOC,X,RO
NSEL,R,LOC,Y,-H
*GET,NBOT,NODE,0,NUM,MIN
ALLSEL,ALL
*DIM,UXT,ARRAY,NMODES
*DIM,UYT,ARRAY,NMODES
*DIM,UXB,ARRAY,NMODES
*DIM,UYB,ARRAY,NMODES
*DO,I,1,NMODES
  SET,1,I
  *GET,UXT(I),NODE,NTOP,U,X
  *GET,UYT(I),NODE,NTOP,U,Y
  *GET,UXB(I),NODE,NBOT,U,X
  *GET,UYB(I),NODE,NBOT,U,Y
*ENDDO
*CFOPEN,ansys_p5_epsrel_ff_lanb_{tag}_discriminator,txt
*VWRITE
('mode   UX_top        UY_top        UX_bot        UY_bot      (bending: UY_top~UY_bot, UX_top~-UX_bot)')
*VWRITE,IDXARR(1),UXT(1),UYT(1),UXB(1),UYB(1)
(F4.0,3X,E12.4,3X,E12.4,3X,E12.4,3X,E12.4)
*CFCLOS
FINISH
"""

PIEZ_COUPLED = """TBDATA,1,0,E31,0,0,E33,0
TBDATA,7,0,E31,0,E15,0,0
TBDATA,13,0,0,E15,0,0,0"""

PIEZ_E0 = """TBDATA,1,0,0,0,0,0,0
TBDATA,7,0,0,0,0,0,0
TBDATA,13,0,0,0,0,0,0"""

EBC_SCBOTH = """! Electrical BC -- both faces VOLT=0 (analytic short-circuit).
NSEL,S,LOC,Y,H
NSEL,A,LOC,Y,-H
D,ALL,VOLT,0
ALLSEL,ALL"""

EBC_OC = """! Electrical BC -- bottom grounded, top floating bus (analytic OC).
! Do not also D,VOLT,0 the top face -- that restores short-circuit.
NSEL,S,LOC,Y,-H
D,ALL,VOLT,0
ALLSEL,ALL
NSEL,S,LOC,Y,H
CP,1,VOLT,ALL
ALLSEL,ALL"""

EBC_E0 = """! Coupling-off control: TB,PIEZ is zero. Both faces VOLT=0 so VOLT
! stays well-posed. Mechanical spectrum must match elastic 4x4.
NSEL,S,LOC,Y,H
NSEL,A,LOC,Y,-H
D,ALL,VOLT,0
ALLSEL,ALL"""

DECKS = [
    dict(tag="e0", title="coupling-off e=0", piez=PIEZ_E0, ebc=EBC_E0,
         target="%.3f" % EL_HZ),
    dict(tag="scboth", title="both-faces SC", piez=PIEZ_COUPLED, ebc=EBC_SCBOTH,
         target="%.3f" % SC_HZ),
    dict(tag="oc", title="open-circuit", piez=PIEZ_COUPLED, ebc=EBC_OC,
         target="%.3f" % SC_HZ),
]


def main():
    taken = existing_jobnames()
    written = []
    for spec in DECKS:
        fname = "ansys_p5_epsrel_ff_lanb_%s_%s.inp" % (spec["tag"], DATE)
        path = os.path.join(HERE, fname)
        jn = jobname(fname)
        if jn in taken and os.path.basename(taken[jn]) != fname:
            raise SystemExit("FATAL jobname collision %r vs %s"
                             % (jn, taken[jn]))
        body = HEADER.format(
            tag=spec["tag"], title=spec["title"], piez=spec["piez"],
            ebc=spec["ebc"], target=spec["target"])
        bad = "***" + " ERROR " + "***"
        if bad in body:
            raise SystemExit("FATAL: deck body contains MAPDL failure token")
        text = body.replace("\r\n", "\n")
        if not text.endswith("\n"):
            text += "\n"
        with open(path, "w", newline="\n", encoding="ascii") as f:
            f.write(text)
        taken[jn] = fname
        written.append((fname, jn))
        print("wrote %s  jobname=%s" % (fname, jn))
    print("OK %d decks" % len(written))


if __name__ == "__main__":
    main()
