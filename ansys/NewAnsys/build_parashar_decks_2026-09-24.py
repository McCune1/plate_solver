# -*- coding: utf-8 -*-
"""
build_parashar_decks_2026-09-24.py -- Phase B scouting: independent
axisymmetric FE of Parashar 2013 (J Intell Mater Syst Struct 24:1572),
radially poled PIC 255 annulus, face electrodes, n = 0 column of Table 2
(free-free, D 25 / d 10 / h 4 mm) and Table 3 (fixed-free, inner rim
clamped, D 24 / d 4 / h 3 mm).  Short circuit (both faces grounded), as
in Parashar's free-vibration model (U0 = 0).

Purpose: before any Mindlin/d15 derivation, find out whether Parashar's
Rayleigh-Ritz numbers or his measurements are the better external
benchmark (the project's usual "FE the published table first" step).
Pre-registration: Project Knowledge/PARASHAR_PHASEB_SCOUT_2026-09-24.md.

Material (Table 1, read from the page image; e31 = -5.6 printed with
its minus sign): rho 7800; cE11 1.108e11, cE12 6.326e10, cE13 6.896e10,
cE33 1.108e11, cE44 1.909e10 N/m^2; e31 -5.6, e33 12.8, e15 10.3 C/m^2;
eps33S/eps0 1161, eps11S/eps0 1023.  cE66 = (cE11 - cE12)/2 (not in the
table; transversely isotropic).
POLING IS RADIAL = ANSYS X in the axisymmetric PLANE223 model, so the
3-axis constants sit on X (not Y as in every other project deck):
  TB,ANEL  x: c33 c13 c13 | y: c11 c12 | z: c11 | xy: c44 | yz: c66 | xz: c44
  TB,PIEZ  x: (e33,0,0)  y: (e31,0,0)  z: (e31,0,0)  xy: (0,e15,0)
           yz: 0  xz: (0,0,e15)
  PERX = eps33S/eps0, PERY = PERZ = eps11S/eps0.
Decks: {ff, cf} x mesh {m60: 60x16, m120: 120x32} SC, plus {ff, cf} m120
with e = 0 (elastic control, shows the size of the piezo stiffening).
"""
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
DATE = "2026-09-24"
GEOM = {"ff": dict(RI=5.0, RO=12.5, H=2.0, clamp=0, table="Table 2 (free-free)",
                   model=[28.291, 131.949, 191.186], exp=[27.919, 131.275, None]),
        "cf": dict(RI=2.0, RO=12.0, H=1.5, clamp=1, table="Table 3 (fixed-free)",
                   model=[12.624, 54.878, 125.815, 197.124],
                   exp=[12.224, 54.752, None, None])}
MESHES = [(60, 16), (120, 32)]
NMODES = 30

DECK = """! {deck}
! Parashar 2013 radially poled PIC 255 annulus, {table}, n = 0, SC ({date}).
! Built by build_parashar_decks_2026-09-24.py. Do not hand-edit; rebuild.
! Pre-registration: PARASHAR_PHASEB_SCOUT_2026-09-24.md.
! Parashar n = 0 column, model (experiment) [kHz]: {tgt}
! Poling = radial = ANSYS X. {piezo_note}
! UNITS: mm-N-s-tonne, real volts; Hz.
! Do not put the MAPDL failure token in this file (queue grep false-fail).
! ============================================================
FINISH
/CLEAR,NOSTART
/TITLE, Parashar 2013 {bc} n=0 SC {piezo_tag} mesh {nr}x{nz} ({date})
RI     = {RI}
RO     = {RO}
H      = {H}
TOL    = 1.0E-4
NDIV_R = {nr}
NDIV_Z = {nz}
NMODES = {nmodes}
CLAMPI = {clamp}
C11E = 1.108E5
C12E = 6.326E4
C13E = 6.896E4
C33E = 1.108E5
C44E = 1.909E4
C66E = (C11E-C12E)/2
E31  = {e31}
E33  = {e33}
E15  = {e15}
P33  = 1161
P11  = 1023
RHO_PIC = 7.800E-9
/PREP7
ET,1,PLANE223
KEYOPT,1,1,1001
KEYOPT,1,3,1
TB,ANEL,1,1,,0
TBDATA,1,C33E,C13E,C13E,0,0,0
TBDATA,7,C11E,C12E,0,0,0
TBDATA,12,C11E,0,0,0
TBDATA,16,C44E,0,0
TBDATA,19,C66E,0
TBDATA,21,C44E
TB,PIEZ,1
TBDATA,1,E33,0,0,E31,0,0
TBDATA,7,E31,0,0,0,E15,0
TBDATA,13,0,0,0,0,0,E15
MP,PERX,1,P33
MP,PERY,1,P11
MP,PERZ,1,P11
MP,DENS,1,RHO_PIC
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
/SOLU
ANTYPE,MODAL
MODOPT,LANB,NMODES,-1.0,,,ON
MXPAND,NMODES,,,NO
! Short circuit: both electroded faces grounded. Rims insulated.
NSEL,S,LOC,Y,H-TOL,H+TOL
NSEL,A,LOC,Y,-H-TOL,-H+TOL
D,ALL,VOLT,0
ALLSEL,ALL
*IF,CLAMPI,EQ,1,THEN
  NSEL,S,LOC,X,RI-TOL,RI+TOL
  D,ALL,UX,0
  D,ALL,UY,0
  ALLSEL,ALL
*ENDIF
SOLVE
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
*DIM,UYM,ARRAY,NMODES
*DO,I,1,NMODES
  SET,1,I
  *GET,UXT(I),NODE,NTOP,U,X
  *GET,UYT(I),NODE,NTOP,U,Y
  *GET,UXB(I),NODE,NBOT,U,X
  *GET,UYB(I),NODE,NBOT,U,Y
  *GET,UYM(I),NODE,NMID,U,Y
*ENDDO
*CFOPEN,{stem}_modes,txt
*VWRITE
('mode  f_Hz  UX_top  UY_top  UX_bot  UY_bot  UY_mid   [r = RO; flexural: UY_top ~ +UY_bot, UX_top ~ -UX_bot]')
*VWRITE,IDXARR(1),FR1(1),UXT(1),UYT(1),UXB(1),UYB(1),UYM(1)
(F4.0,1X,E20.12,1X,E13.5,1X,E13.5,1X,E13.5,1X,E13.5,1X,E13.5)
*CFCLOS
FINISH
"""


def main():
    decks = []
    for bc, g in GEOM.items():
        tgt = ", ".join("%.3f (%s)" % (m, "%.3f" % e if e else "-")
                        for m, e in zip(g["model"], g["exp"]))
        for piezo in ("sc", "e0"):
            meshes = MESHES if piezo == "sc" else MESHES[1:]
            for nr, nz in meshes:
                deck = "ansys_parashar_%s_%s_m%d_%s.inp" % (bc, piezo, nr, DATE)
                stem = "ansys_parashar_%s_%s_m%d" % (bc, piezo, nr)
                on = piezo == "sc"
                txt = DECK.format(
                    deck=deck, table=g["table"], date=DATE, tgt=tgt, bc=bc,
                    piezo_note=("e31 = -5.6, e33 = 12.8, e15 = 10.3 C/m^2 (Table 1)." if on
                                else "ELASTIC CONTROL: all e = 0 (permittivity kept)."),
                    piezo_tag="piezo" if on else "e0",
                    RI=g["RI"], RO=g["RO"], H=g["H"], nr=nr, nz=nz, nmodes=NMODES,
                    clamp=g["clamp"],
                    e31="-5.6E-3" if on else "0.0", e33="12.8E-3" if on else "0.0",
                    e15="10.3E-3" if on else "0.0", stem=stem)
                assert "*** ERROR" not in txt
                with open(os.path.join(HERE, deck), "w", newline="\n") as fh:
                    fh.write(txt)
                decks.append(dict(deck=deck, stem=stem, bc=bc, piezo=piezo, ndiv_r=nr))
    man = os.path.join(HERE, "ansys_queue_manifest_parashar_%s.txt" % DATE)
    with open(man, "w", newline="\n") as fh:
        fh.write("# Parashar 2013 n=0 FE scout (%s). Built by build_parashar_decks_2026-09-24.py.\n"
                 "# Score with python3 score_parashar_fe_2026-09-24.py. ONE licence seat.\n#\n" % DATE)
        for d in decks:
            fh.write(d["deck"] + "\n")
    tg = dict(built_by="build_parashar_decks_2026-09-24.py",
              geometry=GEOM, decks=decks,
              bars=dict(model_rel_s12=0.03, mesh_rel=0.002))
    with open(os.path.join(HERE, "targets_parashar_%s.json" % DATE), "w", newline="\n") as fh:
        json.dump(tg, fh, indent=1)
    print("wrote %d decks" % len(decks))


if __name__ == "__main__":
    main()
