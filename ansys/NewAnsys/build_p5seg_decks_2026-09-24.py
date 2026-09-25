# -*- coding: utf-8 -*-
"""
build_p5seg_decks_2026-09-24.py -- Paper 5 segmented-electrode harmonic FE.

Pre-registration: Project Knowledge/PAPER5_E15_SEGMENT_CHARGE_DERIVATION_2026-09-24.md
Predictions:      Plate_Solver_Package/Paper5_Monolithic/
                  piezo_p5_e15_segment_charge_predictions.json
                  (probe_piezo_p5_e15_segment_charge_2026-09-24.py, PASS_ALL)

Writes six decks (c44 = 73, 26 GPa x NDIV_R/NDIV_Z = 40/6, 80/12,
160/12), the queue manifest and targets_p5seg_2026-09-24.json.
PREP7 material/element block is the Sec 18.201 Paper 5 F-F deck
(ansys_p5_epsrel_ff_lanb_scboth_2026-09-19.inp, job 2520476) with C44E
made a parameter. New: ANTYPE,HARMIC, ring load, per-node CHRG dump.
Every *VWRITE is a single call outside any *DO loop (ansys-kfac-runtime-bug).
No derived quantity is computed in APDL; the scorer does all arithmetic.

Run from Ansys/NewAnsys/:  python3 build_p5seg_decks_2026-09-24.py
"""
import json
import math
import os

import sys

HERE = os.path.dirname(os.path.abspath(__file__))
# e31 sign migration (LESSONS Sec 18.244): "--e31m" builds the e31 = -4.1
# (Liu 2002 / Duan 2005 Table 1 as printed) set; all names carry "_e31m".
E31M = "--e31m" in sys.argv
TAG = "_e31m" if E31M else ""
E31_DECK = "-4.1E-3" if E31M else "4.1E-3"
_PRED_NAME = "piezo_p5_e15_segment_charge_predictions%s.json" % TAG
_PRED_DIRS = [os.path.join(HERE, "..", "..", "Paper5_Monolithic"),          # cluster mirror
              os.path.join(HERE, "..", "..", "validation", "paper5_monolithic")]  # github_repo
PRED = next((os.path.join(d, _PRED_NAME) for d in _PRED_DIRS
             if os.path.exists(os.path.join(d, _PRED_NAME))),
            os.path.join(_PRED_DIRS[0], _PRED_NAME))
DATE = "2026-09-24"
MESHES = [(40, 6), (80, 12), (160, 12)]
C44S = [("c73", 7.3e4, "c44_73"), ("c26", 2.6e4, "c44_26")]
OMEGA_CAL = 100.0

HEAD = """! {deck}
! Paper 5 segmented-electrode sensing charge, harmonic FE (2026-09-24).
! Built by build_p5seg_decks_2026-09-24.py. Do not hand-edit; rebuild.
! Pre-registration: PAPER5_E15_SEGMENT_CHARGE_DERIVATION_2026-09-24.md.
! Question: does the inner-segment face charge follow Kirchhoff
! (e15*gamma_rz dropped, Paper 5 Eq. Ysense) or the consistent
! leading-order model with the e15 shear term (kappa = {kap:+.4f} x Kirchhoff
! at C44E = {c44gpa:.0f} GPa)?
! Geometry/material/element: the Sec 18.201 Paper 5 F-F deck
! (ansys_p5_epsrel_ff_lanb_scboth_2026-09-19.inp), C44E parameterised.
! UNITS: mm-N-s-tonne, real volts. CHRG comes out in mC (1e-3 C).
! LS1: sign calibration. Top inner segment (r <= RSTAR) at 1 V, all
!      other face nodes 0 V, no force. Its CHRG sum must be > 0.
! LS2..LS6: all face nodes 0 V (both faces fully electroded, SC),
!      FY = 1 N total ring load at r = RLOAD (350 mm), spread over the through-
!      thickness node line (symmetric in z: flexure only).
! Output: one file per load step, raw per-node r and CHRG (Re, Im) on
! the top and bottom faces. The scorer does every sum and ratio.
! Do not put the MAPDL failure token in this file (queue grep false-fail).
! ============================================================
FINISH
/CLEAR,NOSTART
/TITLE, Paper 5 segmented-electrode harmonic charge, C44E={c44gpa:.0f} GPa, mesh {nr}x{nz} ({date})

RI     = 100.0
RO     = 600.0
H      = 10.0
RLOAD  = 350.0
RSTAR  = 300.0
TOL    = 1.0E-3
NDIV_R = {nr}
NDIV_Z = {nz}
FTOT   = 1.0

C11E = 1.320E5
C12E = 7.100E4
C13E = 7.300E4
C33E = 1.150E5
C44E = {c44:.4E}
C66E = (C11E-C12E)/2
E31  = {e31deck}
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
TBDATA,1,0,E31,0,0,E33,0
TBDATA,7,0,E31,0,E15,0,0
TBDATA,13,0,0,E15,0,0,0
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
! Mesh guards: RSTAR and RLOAD must be corner-node radii, and the element
! count must be NDIV_R*NDIV_Z. *MSG,FATAL stops the deck before any solve
! (same guard pattern as build_p4_3d_nge1_decks_2026-09-23.py).
*GET,NELEM,ELEM,0,COUNT
*IF,NELEM,NE,NDIV_R*NDIV_Z,THEN
  *MSG,FATAL
Mapped mesh element count does not equal NDIV_R*NDIV_Z.
*ENDIF
NSEL,S,LOC,Y,H-TOL,H+TOL
NSEL,R,LOC,X,RSTAR-TOL,RSTAR+TOL
*GET,NCHK1,NODE,0,COUNT
NSEL,S,LOC,X,RLOAD-TOL,RLOAD+TOL
*GET,NFN,NODE,0,COUNT
ALLSEL,ALL
*IF,NCHK1,NE,1,THEN
  *MSG,FATAL
No single top-face node at r = RSTAR.
*ENDIF
*IF,NFN,NE,2*NDIV_Z+1,THEN
  *MSG,FATAL
Force node line at r = RLOAD is not 2*NDIV_Z+1 nodes.
*ENDIF
FINISH

/SOLU
ANTYPE,HARMIC
HROPT,FULL
HROUT,OFF
OUTRES,NSOL,ALL
OUTRES,RSOL,ALL
KBC,1
EQSLV,SPARSE

! Both faces grounded everywhere.
NSEL,S,LOC,Y,H-TOL,H+TOL
NSEL,A,LOC,Y,-H-TOL,-H+TOL
D,ALL,VOLT,0
ALLSEL,ALL

! LS1 -- sign calibration: top inner segment at 1 V, no force.
NSEL,S,LOC,Y,H-TOL,H+TOL
NSEL,R,LOC,X,RI-TOL,RSTAR+TOL
D,ALL,VOLT,1.0
ALLSEL,ALL
HARFRQ,{fcal:.9f},{fcal:.9f}
NSUBST,1
LSWRITE,1

! LS2..LS6 -- SC, ring load at RLOAD.
NSEL,S,LOC,Y,H-TOL,H+TOL
NSEL,R,LOC,X,RI-TOL,RSTAR+TOL
D,ALL,VOLT,0
NSEL,S,LOC,X,RLOAD-TOL,RLOAD+TOL
F,ALL,FY,FTOT/NFN
ALLSEL,ALL
{lsblocks}
LSSOLVE,1,{nls}
FINISH

/POST1
! Human-readable label check (one PRRSOL, LS2 real part).
SET,2,1,,0
NSEL,S,LOC,Y,H-TOL,H+TOL
NSEL,R,LOC,X,RSTAR-TOL,RSTAR+TOL
PRRSOL,CHRG
ALLSEL,ALL

NSEL,S,LOC,Y,H-TOL,H+TOL
*GET,NTOP,NODE,0,COUNT
ALLSEL,ALL
*DIM,IDX,ARRAY,NTOP
*DIM,RT,ARRAY,NTOP
*DIM,RB,ARRAY,NTOP
*DIM,FRQ,ARRAY,{nls}
*DIM,LSN,ARRAY,{nls}
{postblocks}
*CFOPEN,{stem}_freqs,txt
*VWRITE
('ls   freq_Hz')
*VWRITE,LSN(1),FRQ(1)
(F4.0,3X,E20.12)
*CFCLOS
FINISH
"""

LS_BLOCK = """HARFRQ,{f:.9f},{f:.9f}
NSUBST,1
LSWRITE,{ls}
"""

POST_BLOCK = """! ---- load step {ls} ----
*DIM,QTR{ls},ARRAY,NTOP
*DIM,QTI{ls},ARRAY,NTOP
*DIM,QBR{ls},ARRAY,NTOP
*DIM,QBI{ls},ARRAY,NTOP
SET,{ls},1,,0
*GET,FRQ({ls}),ACTIVE,0,SET,FREQ
LSN({ls})={ls}
NSEL,S,LOC,Y,H-TOL,H+TOL
ND=0
*DO,I,1,NTOP
  ND=NDNEXT(ND)
  IDX(I)=I
  RT(I)=NX(ND)
  *GET,QTR{ls}(I),NODE,ND,RF,CHRG
*ENDDO
NSEL,S,LOC,Y,-H-TOL,-H+TOL
ND=0
*DO,I,1,NTOP
  ND=NDNEXT(ND)
  RB(I)=NX(ND)
  *GET,QBR{ls}(I),NODE,ND,RF,CHRG
*ENDDO
SET,{ls},1,,1
NSEL,S,LOC,Y,H-TOL,H+TOL
ND=0
*DO,I,1,NTOP
  ND=NDNEXT(ND)
  *GET,QTI{ls}(I),NODE,ND,RF,CHRG
*ENDDO
NSEL,S,LOC,Y,-H-TOL,-H+TOL
ND=0
*DO,I,1,NTOP
  ND=NDNEXT(ND)
  *GET,QBI{ls}(I),NODE,ND,RF,CHRG
*ENDDO
ALLSEL,ALL
*CFOPEN,{stem}_ls{ls},txt
*VWRITE
('i  r_top  QtopRe  QtopIm  r_bot  QbotRe  QbotIm   [mm, mC]')
*VWRITE,IDX(1),RT(1),QTR{ls}(1),QTI{ls}(1),RB(1),QBR{ls}(1),QBI{ls}(1)
(F6.0,1X,E20.12,1X,E20.12,1X,E20.12,1X,E20.12,1X,E20.12,1X,E20.12)
*CFCLOS
"""


def main():
    with open(PRED) as fh:
        pred = json.load(fh)
    if pred.get("verdict") != "PASS_ALL":
        raise SystemExit("predictions JSON is not PASS_ALL; refusing to build")
    omegas = [row["omega"] for row in pred["rows"]]
    fcal = OMEGA_CAL / (2 * math.pi)
    freqs = [om / (2 * math.pi) for om in omegas]
    nls = 1 + len(freqs)
    lsblocks = "".join(LS_BLOCK.format(f=f, ls=i + 2) for i, f in enumerate(freqs))
    decks = []
    for tag, c44, key in C44S:
        for nr, nz in MESHES:
            deck = "ansys_p5seg%s_%s_m%d_%s.inp" % (TAG, tag, nr, DATE)
            stem = "ansys_p5seg%s_%s_m%d" % (TAG, tag, nr)
            post = "".join(POST_BLOCK.format(ls=ls, stem=stem)
                           for ls in range(1, nls + 1))
            txt = HEAD.format(deck=deck, kap=pred["kappa_closed"][key],
                              c44gpa=c44 / 1e3, c44=c44, nr=nr, nz=nz,
                              date=DATE, fcal=fcal, lsblocks=lsblocks,
                              nls=nls, postblocks=post, stem=stem, e31deck=E31_DECK)
            assert "*** ERROR" not in txt
            with open(os.path.join(HERE, deck), "w", newline="\n") as fh:
                fh.write(txt)
            decks.append(dict(deck=deck, stem=stem, c44_key=key,
                              c44_MPa=c44, ndiv_r=nr, ndiv_z=nz))
    man = os.path.join(HERE, "ansys_queue_manifest_piezo_p5seg%s_%s.txt" % (TAG, DATE))
    with open(man, "w", newline="\n") as fh:
        fh.write("# Paper 5 segmented-electrode harmonic charge FE (2026-09-24).\n"
                 "# Built by build_p5seg_decks_2026-09-24.py. Score with\n"
                 "# python3 score_p5seg_fe_2026-09-24.py (bars fixed before the run\n"
                 "# in targets_p5seg_2026-09-24.json). ONE licence seat.\n#\n")
        for d in decks:
            fh.write(d["deck"] + "\n")
    targets = dict(
        built_by="build_p5seg_decks_2026-09-24.py",
        predictions_source="Paper5_Monolithic/" + _PRED_NAME, e31_deck=E31_DECK,
        # bars relative to max(|kappa|, kappa_floor); kappa_floor = 1 for the
        # e31 = -4.1 set (kappa(c44=73) = -0.055), 0 = original behaviour.
        kappa_floor=1.0 if E31M else 0.0,
        r_star_mm=300.0, r_i_mm=100.0, charge_unit_to_C=1e-3,
        omega_cal=OMEGA_CAL,
        load_steps=[dict(ls=i + 2, omega=om, f_Hz=freqs[i],
                         Q_K_face_C=pred["rows"][i]["Q_K_face"],
                         kappa={k: pred["rows"][i][k]["kappa"]
                                for k in ("c44_73", "c44_26")})
                    for i, om in enumerate(omegas)],
        kappa_closed=pred["kappa_closed"],
        bars=dict(G_sym_rel=0.01, G_mesh_rel=0.01, corrected_rel=0.05,
                  kirchhoff_rel=0.10,
                  cal_C_range_F=[5e-8, 3e-7]),
        decks=decks)
    with open(os.path.join(HERE, "targets_p5seg%s_%s.json" % (TAG, DATE)), "w",
              newline="\n") as fh:
        json.dump(targets, fh, indent=1)
    for d in decks:
        print("wrote", d["deck"])
    print("wrote", os.path.basename(man), "and targets_p5seg%s_%s.json" % (TAG, DATE))


if __name__ == "__main__":
    main()
