# -*- coding: utf-8 -*-
"""
build_pcm_decks_2026-09-24.py -- piezoelectric contour-mode (in-plane,
n = 0 radial) FE for the Paper 5 monolithic PZT-4 ring. First FE of the
in-plane piezoelectric paper (Paper 4 lever 13 / Paper 5 claim C).

Pre-registration: Project Knowledge/CONTOUR_MODE_N0_DERIVATION_2026-09-24.md Sec 5
Predictions:      Plate_Solver_Package/Contour_Mode/piezo_contour_n0_predictions.json
                  (probe_piezo_contour_n0_2026-09-24.py, PASS_ALL)

12 decks: e31 {+4.1 (project convention), -4.1 (Duan/Liu print, physical
sign)} x {sc modal, oc modal, harm} x mesh NDIV_R x NDIV_Z {80x12, 160x24}.
  sc   : both faces VOLT = 0, LANB, 30 modes, discriminator at r = RO.
  oc   : bottom VOLT = 0, top face CP-coupled floating (as job 2520476).
  harm : top face VOLT = 1, bottom 0, HARMIC at 1 Hz (static calibration:
         free capacitance C_T and the CHRG sign) and at the model's
         off-resonance test points; per-node top-face CHRG dump.
Geometry, material block and element: identical to the Paper 5 decks
(ansys_p5_epsrel_ff_lanb_*_2026-09-19.inp), except E31 per case.
Every *VWRITE is a single call outside any *DO loop (ansys-kfac-runtime-bug).

Run from Ansys/NewAnsys/:  python3 build_pcm_decks_2026-09-24.py
"""
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
_PRED_NAME = "piezo_contour_n0_predictions.json"
_PRED_DIRS = [os.path.join(HERE, "..", "..", "Contour_Mode"),
              os.path.join(HERE, "..", "..", "validation", "contour_mode")]
PRED = next((os.path.join(d, _PRED_NAME) for d in _PRED_DIRS
             if os.path.exists(os.path.join(d, _PRED_NAME))),
            os.path.join(_PRED_DIRS[0], _PRED_NAME))
DATE = "2026-09-24"
MESHES = [(80, 12), (160, 24)]
CASES = [("p41", "e31_p41", "4.1E-3"), ("m41", "e31_m41", "-4.1E-3")]
NMODES = 30
F_STATIC = 1.0

COMMON = """FINISH
/CLEAR,NOSTART
/TITLE, Contour-mode PZT-4 ring {kind} e31={e31tag} mesh {nr}x{nz} ({date})
RI     = 100.0
RO     = 600.0
H      = 10.0
TOL    = 1.0E-3
NDIV_R = {nr}
NDIV_Z = {nz}
NMODES = {nmodes}
C11E = 1.320E5
C12E = 7.100E4
C13E = 7.300E4
C33E = 1.150E5
C44E = 7.300E4
C66E = (C11E-C12E)/2
E31  = {e31}
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
*GET,NELEM,ELEM,0,COUNT
*IF,NELEM,NE,NDIV_R*NDIV_Z,THEN
  *MSG,FATAL
Mapped mesh element count does not equal NDIV_R*NDIV_Z.
*ENDIF
FINISH
"""

HDR = """! {deck}
! Piezoelectric contour-mode (n = 0 radial) FE, Paper 5 monolithic ring ({date}).
! Built by build_pcm_decks_2026-09-24.py. Do not hand-edit; rebuild.
! Pre-registration: CONTOUR_MODE_N0_DERIVATION_2026-09-24.md Sec 5.
! e31 = {e31tag} C/m^2 ({e31note}). {kindnote}
! Thin-plate model (probe_piezo_contour_n0_2026-09-24.py), radial modes:
{modelines}
! UNITS: mm-N-s-tonne, real volts; Hz; CHRG in mC.
! Do not put the MAPDL failure token in this file (queue grep false-fail).
! ============================================================
"""

MODAL = """/SOLU
ANTYPE,MODAL
MODOPT,LANB,NMODES,-1.0,,,ON
MXPAND,NMODES,,,NO
{elec}SOLVE
FINISH
/POST1
*DIM,FR1,ARRAY,NMODES
*DIM,IDXARR,ARRAY,NMODES
*DO,I,1,NMODES
  *GET,FR1(I),MODE,I,FREQ
  IDXARR(I)=I
*ENDDO
NSEL,S,LOC,X,RO-TOL,RO+TOL
NSEL,R,LOC,Y,H-TOL,H+TOL
*GET,NTOP,NODE,0,NUM,MIN
NSEL,S,LOC,X,RO-TOL,RO+TOL
NSEL,R,LOC,Y,-H-TOL,-H+TOL
*GET,NBOT,NODE,0,NUM,MIN
NSEL,S,LOC,X,RO-TOL,RO+TOL
NSEL,R,LOC,Y,-TOL,TOL
*GET,NMID,NODE,0,NUM,MIN
ALLSEL,ALL
*DIM,UXT,ARRAY,NMODES
*DIM,UYT,ARRAY,NMODES
*DIM,UXB,ARRAY,NMODES
*DIM,UYB,ARRAY,NMODES
*DIM,UXM,ARRAY,NMODES
*DO,I,1,NMODES
  SET,1,I
  *GET,UXT(I),NODE,NTOP,U,X
  *GET,UYT(I),NODE,NTOP,U,Y
  *GET,UXB(I),NODE,NBOT,U,X
  *GET,UYB(I),NODE,NBOT,U,Y
  *GET,UXM(I),NODE,NMID,U,X
*ENDDO
*CFOPEN,{stem}_modes,txt
*VWRITE
('mode  f_Hz  UX_top  UY_top  UX_bot  UY_bot  UX_mid   [r = RO; extensional: UX_top ~ +UX_bot]')
*VWRITE,IDXARR(1),FR1(1),UXT(1),UYT(1),UXB(1),UYB(1),UXM(1)
(F4.0,1X,E20.12,1X,E13.5,1X,E13.5,1X,E13.5,1X,E13.5,1X,E13.5)
*CFCLOS
FINISH
"""

ELEC_SC = """NSEL,S,LOC,Y,H-TOL,H+TOL
NSEL,A,LOC,Y,-H-TOL,-H+TOL
D,ALL,VOLT,0
ALLSEL,ALL
"""
ELEC_OC = """! Bottom grounded, top face one floating equipotential (open circuit).
NSEL,S,LOC,Y,-H-TOL,-H+TOL
D,ALL,VOLT,0
NSEL,S,LOC,Y,H-TOL,H+TOL
CP,1,VOLT,ALL
ALLSEL,ALL
"""

HARM = """/SOLU
ANTYPE,HARMIC
HROPT,FULL
HROUT,OFF
OUTRES,NSOL,ALL
OUTRES,RSOL,ALL
KBC,1
EQSLV,SPARSE
NSEL,S,LOC,Y,-H-TOL,-H+TOL
D,ALL,VOLT,0
NSEL,S,LOC,Y,H-TOL,H+TOL
D,ALL,VOLT,1.0
ALLSEL,ALL
{lsblocks}
LSSOLVE,1,{nls}
FINISH
/POST1
NSEL,S,LOC,Y,H-TOL,H+TOL
*GET,NTOP,NODE,0,COUNT
ALLSEL,ALL
*DIM,IDX,ARRAY,NTOP
*DIM,RT,ARRAY,NTOP
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
*DIM,QR{ls},ARRAY,NTOP
*DIM,QI{ls},ARRAY,NTOP
SET,{ls},1,,0
*GET,FRQ({ls}),ACTIVE,0,SET,FREQ
LSN({ls})={ls}
NSEL,S,LOC,Y,H-TOL,H+TOL
ND=0
*DO,I,1,NTOP
  ND=NDNEXT(ND)
  IDX(I)=I
  RT(I)=NX(ND)
  *GET,QR{ls}(I),NODE,ND,RF,CHRG
*ENDDO
SET,{ls},1,,1
ND=0
*DO,I,1,NTOP
  ND=NDNEXT(ND)
  *GET,QI{ls}(I),NODE,ND,RF,CHRG
*ENDDO
ALLSEL,ALL
*CFOPEN,{stem}_ls{ls},txt
*VWRITE
('i  r_top  QtopRe  QtopIm   [mm, mC]')
*VWRITE,IDX(1),RT(1),QR{ls}(1),QI{ls}(1)
(F6.0,1X,E20.12,1X,E20.12,1X,E20.12)
*CFCLOS
"""


def main():
    pred = json.load(open(PRED))
    if pred.get("verdict") != "PASS_ALL":
        raise SystemExit("predictions JSON is not PASS_ALL; refusing to build")
    decks = []
    for etag, ckey, e31 in CASES:
        case = pred["cases"][ckey]
        modelines = "\n".join(
            "!   m%d  fr %.4f Hz  fa %.4f Hz  split %+.4f%%  keff^2 %.5f"
            % (m["m"], m["fr_Hz"], m["fa_Hz"], m["split_pct"], m["keff2"])
            for m in case["modes"])
        note = ("project convention, Papers 4-6" if etag == "p41"
                else "sign printed by Duan 2005 / Liu 2002; physical PZT-4 e31 < 0")
        freqs = [F_STATIC] + [p["f_Hz"] for p in case["harm_points"]]
        for kind in ("sc", "oc", "harm"):
            for nr, nz in MESHES:
                deck = "ansys_pcm_%s_%s_m%d_%s.inp" % (etag, kind, nr, DATE)
                stem = "ansys_pcm_%s_%s_m%d" % (etag, kind, nr)
                kindnote = {"sc": "Modal, short circuit (both faces grounded).",
                            "oc": "Modal, open circuit (top face floating, CP).",
                            "harm": "Harmonic, top face at 1 V, bottom grounded."}[kind]
                txt = HDR.format(deck=deck, date=DATE,
                                 e31tag="%+.1f (deck value %s in mm-N-s units)"
                                 % (float(e31) * 1e3, e31),
                                 e31note=note, kindnote=kindnote, modelines=modelines)
                txt += COMMON.format(kind=kind, e31tag=etag, nr=nr, nz=nz, date=DATE,
                                     nmodes=NMODES, e31=e31)
                if kind in ("sc", "oc"):
                    txt += MODAL.format(elec=ELEC_SC if kind == "sc" else ELEC_OC, stem=stem)
                else:
                    nls = len(freqs)
                    lsb = "".join(LS_BLOCK.format(f=f, ls=i + 1) for i, f in enumerate(freqs))
                    post = "".join(POST_BLOCK.format(ls=i + 1, stem=stem) for i in range(nls))
                    txt += HARM.format(lsblocks=lsb, nls=nls, postblocks=post, stem=stem)
                assert "*** ERROR" not in txt
                with open(os.path.join(HERE, deck), "w", newline="\n") as fh:
                    fh.write(txt)
                decks.append(dict(deck=deck, stem=stem, e31_key=ckey, kind=kind,
                                  ndiv_r=nr, ndiv_z=nz,
                                  harm_freqs_Hz=freqs if kind == "harm" else None))
    man = os.path.join(HERE, "ansys_queue_manifest_piezo_pcm_%s.txt" % DATE)
    with open(man, "w", newline="\n") as fh:
        fh.write("# Piezoelectric contour-mode (n=0 radial) FE, Paper 5 ring (%s).\n"
                 "# Built by build_pcm_decks_2026-09-24.py. Score with\n"
                 "# python3 score_pcm_fe_2026-09-24.py. ONE licence seat.\n#\n" % DATE)
        for d in decks:
            fh.write(d["deck"] + "\n")
    targets = dict(
        built_by="build_pcm_decks_2026-09-24.py",
        predictions_source="Contour_Mode/piezo_contour_n0_predictions.json",
        cases=dict((k, dict(modes=v["modes"], C_T_F=v["C_T_F"], C_S_F=v["C_S_F"],
                            harm_points=v["harm_points"]))
                   for k, v in pred["cases"].items()),
        charge_unit_to_C=1e-3, chrg_sign_expected=-1,
        bars=dict(fr_rel_m123=0.005, split_rel_m123=0.05, keff2_min=1e-3,
                  mesh_fr_rel=5e-4, static_rel=0.01, harm_rel_below_fr2=0.01),
        decks=decks)
    with open(os.path.join(HERE, "targets_pcm_%s.json" % DATE), "w", newline="\n") as fh:
        json.dump(targets, fh, indent=1)
    print("wrote %d decks, %s, targets_pcm_%s.json" % (len(decks), os.path.basename(man), DATE))


if __name__ == "__main__":
    main()
