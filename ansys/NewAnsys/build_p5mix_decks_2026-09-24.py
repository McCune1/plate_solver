# -*- coding: utf-8 -*-
"""
build_p5mix_decks_2026-09-24.py -- Paper 5 segmented-electrode sensing
charge for the clamped-clamped and mixed-edge rings (C-C, C-F, F-C), with
an F-F control.  Harmonic, axisymmetric PLANE223, same deck family as
build_p5seg_decks_2026-09-24.py (job 2531984, CORRECTED).

Pre-registration: Project Knowledge/PAPER5_MIXED_EDGE_SENSING_2026-09-24.md
Predictions:      Plate_Solver_Package/Paper5_Monolithic/
                  piezo_p5_mixed_sense_predictions.json
                  (probe_piezo_p5_mixed_sense_2026-09-24.py)

Writes 24 decks: BC {ff, cc, cf, fc} x C44E {73, 26 GPa} x mesh
NDIV_R x NDIV_Z {80x12, 160x12, 320x24}, the queue manifest and
targets_p5mix_2026-09-24.json.  What is new relative to p5seg:
  * clamped rim(s): UX = UY = 0 on every node of r = RI and/or r = RO
    (the Paper 4 C-C deck convention, job 2490706); rims stay
    electrically insulated (no VOLT constraint on the rim);
  * per-BC test frequencies from the predictions JSON (5 per BC, each
    kept >= 8% from any zero of the Kirchhoff segment charge);
  * mesh guards for every scored cut radius (112.5, 125, 300, 450, 575 mm).
The per-node CHRG dump is unchanged, so every cut is scored post hoc.
Every *VWRITE is a single call outside any *DO loop (ansys-kfac-runtime-bug).
No derived quantity is computed in APDL; the scorer does all arithmetic.

Run from Ansys/NewAnsys/:  python3 build_p5mix_decks_2026-09-24.py
"""
import json
import math
import os

import sys

HERE = os.path.dirname(os.path.abspath(__file__))
# e31 sign migration (LESSONS Sec 18.244): "--e31m" builds the e31 = -4.1
# (Liu 2002 / Duan 2005 Table 1 as printed) set from the _e31m predictions.
E31M = "--e31m" in sys.argv
TAG = "_e31m" if E31M else ""
E31_DECK = "-4.1E-3" if E31M else "4.1E-3"
_PRED_NAME = "piezo_p5_mixed_sense_predictions%s.json" % TAG
_PRED_DIRS = [os.path.join(HERE, "..", "..", "Paper5_Monolithic"),
              os.path.join(HERE, "..", "..", "validation", "paper5_monolithic")]
PRED = next((os.path.join(d, _PRED_NAME) for d in _PRED_DIRS
             if os.path.exists(os.path.join(d, _PRED_NAME))),
            os.path.join(_PRED_DIRS[0], _PRED_NAME))
DATE = "2026-09-24"
MESHES = [(80, 12), (160, 12), (320, 24)]
C44S = [("c73", 7.3e4, "c44_73"), ("c26", 2.6e4, "c44_26")]
BCS = [("ff", "FF", "free inner, free outer"),
       ("cc", "CC", "clamped inner, clamped outer"),
       ("cf", "CF", "clamped inner, free outer"),
       ("fc", "FC", "free inner, clamped outer")]
CUTS_MM = [112.5, 125.0, 300.0, 450.0, 575.0]
OMEGA_CAL = 100.0

HEAD = """! {deck}
! Paper 5 segmented-electrode sensing charge, mixed edges, harmonic FE ({date}).
! Built by build_p5mix_decks_2026-09-24.py. Do not hand-edit; rebuild.
! Pre-registration: PAPER5_MIXED_EDGE_SENSING_2026-09-24.md.
! Mechanical BC: {bcname} ({bclong}). Rims electrically insulated.
! Question 1 (interior cuts r* = 300, 450 mm): does the e15-consistent
!   factor kappa = {kap:+.4f} (C44E = {c44gpa:.0f} GPa) carry over from F-F
!   (job 2531984) to this edge pairing?
! Question 2 (cuts within 2.5H of a clamped rim): which rim model, if any,
!   matches -- M1 phibar'=0 or M2 insulated rim with the e15 shear term?
! Geometry/material/element identical to ansys_p5seg_* (C44E parameterised).
! UNITS: mm-N-s-tonne, real volts. CHRG comes out in mC (1e-3 C).
!   ANSYS 2025R1 reports CHRG reaction = -(electrode charge) (job 2531984);
!   LS1 re-measures that sign in every deck.
! LS1: sign calibration. Top inner segment (r <= RSTAR) at 1 V, others 0 V,
!      no force.
! LS2..LS{nls}: all face nodes 0 V (SC), FY = 1 N total ring load at r = RLOAD
!      spread over the through-thickness node line (flexure only).
! Output: one file per load step, raw per-node r and CHRG (Re, Im) on the
! top and bottom faces. The scorer does every sum and ratio.
! Do not put the MAPDL failure token in this file (queue grep false-fail).
! ============================================================
FINISH
/CLEAR,NOSTART
/TITLE, Paper 5 mixed-edge segment charge {bcname}, C44E={c44gpa:.0f} GPa, mesh {nr}x{nz} ({date})

RI     = 100.0
RO     = 600.0
H      = 10.0
RLOAD  = 350.0
RSTAR  = 300.0
TOL    = 1.0E-3
NDIV_R = {nr}
NDIV_Z = {nz}
FTOT   = 1.0
CLAMPI = {clampi}
CLAMPO = {clampo}

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
! Mesh guards: every scored cut radius and RLOAD must be a corner-node
! radius, and the element count must be NDIV_R*NDIV_Z. *MSG,FATAL stops
! the deck before any solve.
*GET,NELEM,ELEM,0,COUNT
*IF,NELEM,NE,NDIV_R*NDIV_Z,THEN
  *MSG,FATAL
Mapped mesh element count does not equal NDIV_R*NDIV_Z.
*ENDIF
{cutguards}NSEL,S,LOC,X,RLOAD-TOL,RLOAD+TOL
*GET,NFN,NODE,0,COUNT
ALLSEL,ALL
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

! Mechanical clamps (all nodes of the rim line, UX = UY = 0).
*IF,CLAMPI,EQ,1,THEN
  NSEL,S,LOC,X,RI-TOL,RI+TOL
  D,ALL,UX,0
  D,ALL,UY,0
  ALLSEL,ALL
*ENDIF
*IF,CLAMPO,EQ,1,THEN
  NSEL,S,LOC,X,RO-TOL,RO+TOL
  D,ALL,UX,0
  D,ALL,UY,0
  ALLSEL,ALL
*ENDIF

! LS1 -- sign calibration: top inner segment at 1 V, no force.
NSEL,S,LOC,Y,H-TOL,H+TOL
NSEL,R,LOC,X,RI-TOL,RSTAR+TOL
D,ALL,VOLT,1.0
ALLSEL,ALL
HARFRQ,{fcal:.9f},{fcal:.9f}
NSUBST,1
LSWRITE,1

! LS2.. -- SC, ring load at RLOAD.
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

CUT_GUARD = """NSEL,S,LOC,Y,H-TOL,H+TOL
NSEL,R,LOC,X,{r:.4f}-TOL,{r:.4f}+TOL
*GET,NCK{i},NODE,0,COUNT
ALLSEL,ALL
*IF,NCK{i},NE,1,THEN
  *MSG,FATAL
No single top-face node at a scored cut radius (guard {i}).
*ENDIF
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
    algebra = ("G_reduce", "G_kernel", "G_force", "G_recip")
    if not all(pred.get(g) for g in algebra):
        raise SystemExit("predictions JSON fails an algebra gate %s; refusing to build"
                         % ({g: pred.get(g) for g in algebra},))
    if sorted(pred["bcs"]) != ["CC", "CF", "FC", "FF"]:
        raise SystemExit("predictions JSON does not cover all four BCs")
    fcal = OMEGA_CAL / (2 * math.pi)
    guards = "".join(CUT_GUARD.format(r=r, i=i + 1) for i, r in enumerate(CUTS_MM))
    decks = []
    load_steps = {}
    for tag, key, bclong in BCS:
        bcp = pred["bcs"][key]
        omegas = [row["omega"] for row in bcp["rows"]]
        freqs = [om / (2 * math.pi) for om in omegas]
        nls = 1 + len(freqs)
        load_steps[key] = [dict(ls=i + 2, omega=row["omega"],
                                f_Hz=freqs[i], cuts=row["cuts"])
                           for i, row in enumerate(bcp["rows"])]
        lsblocks = "".join(LS_BLOCK.format(f=f, ls=i + 2) for i, f in enumerate(freqs))
        for ctag, c44, ckey in C44S:
            for nr, nz in MESHES:
                deck = "ansys_p5mix%s_%s_%s_m%d_%s.inp" % (TAG, tag, ctag, nr, DATE)
                stem = "ansys_p5mix%s_%s_%s_m%d" % (TAG, tag, ctag, nr)
                post = "".join(POST_BLOCK.format(ls=ls, stem=stem)
                               for ls in range(1, nls + 1))
                txt = HEAD.format(deck=deck, date=DATE, bcname=key, bclong=bclong,
                                  kap=pred["kappa_closed"][ckey],
                                  c44gpa=c44 / 1e3, c44=c44, nr=nr, nz=nz,
                                  clampi=1 if key[0] == "C" else 0,
                                  clampo=1 if key[1] == "C" else 0,
                                  cutguards=guards, fcal=fcal,
                                  lsblocks=lsblocks, nls=nls,
                                  postblocks=post, stem=stem, e31deck=E31_DECK)
                assert "*** ERROR" not in txt
                with open(os.path.join(HERE, deck), "w", newline="\n") as fh:
                    fh.write(txt)
                decks.append(dict(deck=deck, stem=stem, bc=key, c44_key=ckey,
                                  c44_MPa=c44, ndiv_r=nr, ndiv_z=nz))
    man = os.path.join(HERE, "ansys_queue_manifest_piezo_p5mix%s_%s.txt" % (TAG, DATE))
    with open(man, "w", newline="\n") as fh:
        fh.write("# Paper 5 mixed-edge segmented-electrode harmonic charge FE (%s).\n"
                 "# Built by build_p5mix_decks_2026-09-24.py. Score with\n"
                 "# python3 score_p5mix_fe_2026-09-24.py (bars fixed before the run\n"
                 "# in targets_p5mix_2026-09-24.json). ONE licence seat.\n#\n" % DATE)
        for d in decks:
            fh.write(d["deck"] + "\n")
    targets = dict(
        built_by="build_p5mix_decks_2026-09-24.py",
        predictions_source="Paper5_Monolithic/" + _PRED_NAME, e31_deck=E31_DECK,
        predictions_verdict=pred["verdict"],
        predictions_G_kappa=pred["G_kappa"],
        predictions_G_kappa_fail_points=pred.get("G_kappa_fail_points_all", []),
        r_i_mm=100.0, r_o_mm=600.0, H_mm=10.0, r_star_cal_mm=300.0,
        cuts_mm={"in_1125": 112.5, "in_125": 125.0, "int_300": 300.0,
                 "int_450": 450.0, "out_575": 575.0},
        interior_cuts=["int_300", "int_450"],
        clamped_rim_cuts={"CC": ["in_1125", "in_125", "out_575"],
                          "CF": ["in_1125", "in_125"],
                          "FC": ["out_575"], "FF": []},
        free_rim_cuts={"CC": [], "CF": ["out_575"],
                       "FC": ["in_1125", "in_125"],
                       "FF": ["in_1125", "in_125", "out_575"]},
        charge_unit_to_C=1e-3, omega_cal=OMEGA_CAL,
        kappa_closed=pred["kappa_closed"],
        load_steps=load_steps,
        # Bars are relative to max(|kappa|, 1): at e31 = -4.1, c44 = 73 GPa,
        # kappa = -0.055 (near-cancellation); a pure relative bar there would
        # test the dropped O((Hk)^2) terms. Identical to the relative bar
        # whenever |kappa| >= 1 (every other case). Fixed before any run.
        kappa_floor=1.0,
        bars=dict(G_sym_rel=0.01, G_mesh_rel=0.01, corrected_rel=0.05,
                  kirchhoff_rel=0.10, rim_model_rel=0.03, rim_mesh_rel=0.02,
                  cal_C_range_F=[5e-8, 3e-7]),
        decks=decks)
    with open(os.path.join(HERE, "targets_p5mix%s_%s.json" % (TAG, DATE)), "w",
              newline="\n") as fh:
        json.dump(targets, fh, indent=1)
    print("wrote %d decks, %s, targets_p5mix%s_%s.json"
          % (len(decks), os.path.basename(man), TAG, DATE))


if __name__ == "__main__":
    main()
