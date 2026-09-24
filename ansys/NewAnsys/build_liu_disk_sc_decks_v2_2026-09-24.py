#!/usr/bin/env python3
"""Build Liu 2002 solid-disk SC FE decks (Paper 4 lever 11) -- v2.

v2 (2026-09-24, Claude review of job 2530987): the v1 decks put global
UX=0 on every node of the r=RI_MESH hole surface through the full
thickness. An n=1 (cos theta) bending mode tilts the plate about the
y-axis at the centre, so u_x = -z dw/dx is NONZERO there; the constraint
acted as a partial clamp on that tilt and inflated every n=1 row
(+5.6/+6.7% for m=1, +2.2/+2.9% for m=2). n=0 and n>=2 have w'(0)=0 and
were barely touched. v2 leaves the hole FREE (a small free hole is a
regular perturbation; a constrained one is the singular tiny-hole trap).
Removing it frees rigid x-translation in the S deck (UZ ring + UY sym
plane do not block it), so S pins UX=0 at ONE midplane node on the outer
edge at theta=90 deg (x=0): pure flexure of a symmetric laminate has zero
midplane in-plane displacement, so this pin does not touch bending.
Stems/jobnames changed (liudisk2) so v1 .QUEUE_OK sentinels and the v1
jobname .err/.log files cannot collide.


Cloned from build_p4_3d_nge1_decks_2026-09-23.py (SOLID226 half-model,
mm-N-s-tonne, relative permittivity, stress-major TB,PIEZ, both-faces SC).

MESH FIX (2026-09-24, after failed job 2530974)
----------------------------------------------
True RI=0 half-disk volumes are 5-area wedges: each constant-z face is a
3-sided half-disk (2 radials + 1 outer arc). ANSYS mapped-hex rules require
the three sides of a triangular face to have *equal, even* divisions
(help: "If the area is bounded by three lines, the number of element
divisions must be even and equal on all sides"). The first drop applied the
annular opposite-edge pattern (radials NDIV_R=24, arcs NDIV_TH=72), so
Volume 1 failed with:

  Line 1 has 24 element divisions, but line 3 has 72 divisions.
  Volume 1 cannot be meshed with hexahedra.

Equalising all three triangle edges would satisfy the rule but destroys the
structured polar nr x nt node lattice needed for sample-node asserts and
(n,m) lobe ID. Instead we restore the parent *annular* 4-sided topology with
a tiny mesh hole:

  RI_MESH = 3.0 mm = 0.5% of RO   (Liu physical RI = 0)

Inner arcs get NDIV_TH (match outer arcs); radials get NDIV_R; 28 lines after
VGLUE — identical LESIZE pairing to build_p4_3d_nge1. v2: the hole is
FREE (no inner constraint of any kind; see v2 note above).

Geometry / physics (unchanged intent)
-------------------------------------
  * Liu Table 1: r0=0.6 m, h=0.01 m, h1=0.002 m
  * Liu Table 1 PZT-4 C55E = 26 GPa (NOT Duan C44E=73 GPa)
  * e31 = +4.1 (project convention)
  * Outer C: UX=UY=UZ=0 on the outer cylindrical surface
  * Outer S: UZ=0 on the outer circumference at midplane z=0 only
  * Electrical SC: VOLT=0 on both faces of each piezo (scboth)
  * HALF MODEL theta=0..180, UY=0 on y=0 (cos(n theta) family)

Mesh: moderate = parent coarse (24 x 72 x (2+1+1)) = 6912 elements.

Deck text is LF. The MAPDL failure token is never written into a deck.
Do not bump SOLVER_VERSION. Do not edit piezo_disk.py from this builder.
"""
# (no __future__: cluster python may be 3.6)

import math
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
DATE = "2026-09-24"

# Liu physical solid disk RI=0; FE uses a tiny hole so mapped hex works.
RI_LIU = 0.0
RI, RO, H, H1 = 3.0, 600.0, 10.0, 2.0  # mm; RI_MESH = 0.5% of RO
RMID = 0.5 * (RI + RO)
NMODES = 30  # need room for n=0..2, m=1..2 bending + extensional interlopers

# Moderate mesh = confirmed 3-D coarse parent density (annular LESIZE)
NDIV_R, NDIV_TH, NDIV_ZH, NDIV_ZP = 24, 72, 2, 1

DECKS = [
    # (tag, bc)  bc in {C, S}
    ("C", "C"),
    ("S", "S"),
]

# Sample at corner nodes of this mesh (exact radial stations of the tiny-hole grid)
_dr = (RO - RI) / float(NDIV_R)
SAMPLE_RADII = tuple(RI + _dr * k for k in (6, 12, 18, 24))
SAMPLE_THETAS = tuple(15.0 * k for k in range(13))  # 0..180

PIEZ_ON = """TBDATA,1,0,0,E31,0,0,E31
TBDATA,7,0,0,E33,0,0,0
TBDATA,13,0,E15,0,E15,0,0"""

BC_TEXT = {
    "C": "outer clamped: UX=UY=UZ=0 on r=RO cylindrical surface",
    "S": "outer SS Kirchhoff analogue: UZ=0 on r=RO midplane ring z=0 only",
}


def stem(bc):
    return "ansys_p4_liudisk2_sc_%s_%s" % (bc, DATE)


def jobname(fname):
    s = os.path.basename(fname)
    if s.endswith(".inp"):
        s = s[:-4]
    if s.startswith("ansys_"):
        s = s[6:]
    s = re.sub(r"[^A-Za-z0-9]", "_", s)
    return s[:24]


def face_nodes(n1, n2):
    """Nodes on one mapped face of 20-node bricks, n1 x n2 elements."""
    return (n1 + 1) * (n2 + 1) + n1 * (n2 + 1) + (n1 + 1) * n2


def samples():
    zt = H + H1
    out = []
    for r in SAMPLE_RADII:
        for th in SAMPLE_THETAS:
            out.append((r, th, zt))
    for th in SAMPLE_THETAS:
        out.append((RO, th, -zt))
    return out


def line_block(label, xlo, xhi, zlo, zhi, expect, ndiv):
    return (
        "! {label}\n"
        "LSEL,S,LOC,X,{xlo:.6f},{xhi:.6f}\n"
        "LSEL,R,LOC,Z,{zlo:.6f},{zhi:.6f}\n"
        "*GET,NL,LINE,0,COUNT\n"
        "*IF,NL,NE,{expect},THEN\n"
        "  *MSG,FATAL\n"
        "{label}: line count was not {expect}.\n"
        "*ENDIF\n"
        "LESIZE,ALL,,,{ndiv}\n"
        "ALLSEL,ALL\n"
    ).format(label=label, xlo=xlo, xhi=xhi, zlo=zlo, zhi=zhi,
             expect=expect, ndiv=ndiv)


def build(bc):
    nr, nt, nzh, nzp = NDIV_R, NDIV_TH, NDIV_ZH, NDIV_ZP
    nz = nzh + 2 * nzp
    h1 = H1
    nel = nr * nt * nz
    # Annular half-ring electrode face (parent formula; tiny hole included)
    n_elec_face = face_nodes(nr, nt)
    # Symmetry plane y=0: two cut faces (parent formula)
    n_sym = 2 * face_nodes(nr, nz)
    # Inner cylindrical surface (v2: FREE; count asserted only): nt x nz face
    n_inner = face_nodes(nt, nz)
    s = stem(bc)
    smp = samples()
    nsamp = len(smp)
    ztol = min(0.01, h1 / 20.0)
    ltol = h1 / 10.0

    dr = (RO - RI) / nr
    for r, th, _ in smp:
        kr = (r - RI) / (dr / 2.0)
        kt = th / (180.0 / nt / 2.0)
        assert abs(kr - round(kr)) < 1e-9, (r,)
        assert abs(kt - round(kt)) < 1e-9, (th,)
        assert int(round(kt)) % 2 == 0, (th,)

    zl = [-H - h1, -H, H, H + h1]
    layers = [("bottom skin", -H - h1 / 2.0, nzp),
              ("host", 0.0, nzh),
              ("top skin", H + h1 / 2.0, nzp)]

    L = []
    a = L.append
    a("! %s.inp" % s)
    a("! Paper 4 lever 11 -- Liu 2002 solid-disk SC FE (Gate 2 CPT targets).")
    a("! Built by build_liu_disk_sc_decks_v2_%s.py. Do not hand-edit; rebuild." % DATE)
    a("! MESH FIX after job 2530974: true RI=0 is a 3-sided wedge; ANSYS mapped")
    a("! hex needs equal even divisions on all three triangle edges. Equalising")
    a("! would break polar sample lattice. Restored parent annular LESIZE with")
    a("! RI_MESH=%.1f mm (Liu RI=0; hole = %.2f%% of RO). Inner arcs=NDIV_TH," % (
        RI, 100.0 * RI / RO))
    a("! radials=NDIV_R, 28 lines. v2: hole r=RI_MESH is FREE (v1 UX=0 there")
    a("! stiffened every n=1 mode: it blocks the centre tilt u_x=-z dw/dx).")
    a("! S deck: UX=0 at ONE outer midplane node (theta=90) kills rigid x-shift.")
    a("! SOLID226 KEYOPT(1)=1001, theta=0..180, UY=0 on y=0 (cos(n theta)).")
    a("! Both-faces SC. Outer BC: %s." % BC_TEXT[bc])
    a("! Mesh moderate: %d radial x %d circ x (%d host + %d + %d skins) = %d elements." % (
        nr, nt, nzh, nzp, nzp, nel))
    a("! Approx DOF O(2e5). Units mm-N-s-tonne, real volts, RELATIVE permittivity.")
    a("! Liu Table 1: C55E=26 GPa (TB C44E), e31=+4.1 project sign, h1=2 mm.")
    a("! Prints raw Hz and raw UX/UY/UZ/VOLT at %d sample nodes only." % nsamp)
    a("! Score externally from raw Hz; never trust deck percentage columns.")
    a("! ============================================================")
    a("FINISH")
    a("/CLEAR,NOSTART")
    a("/TITLE, Paper 4 Liu solid-disk SC outer-%s SOLID226 (%s)" % (bc, DATE))
    a("")
    a("RI     = %.1f            ! RI_MESH (Liu physical RI=0); tiny hole for mapped hex" % RI)
    a("RO     = %.1f            ! Liu r0 = 0.6 m" % RO)
    a("H      = %.1f            ! host HALF-thickness [mm]" % H)
    a("H1     = %.1f             ! each piezo skin [mm] (Liu h1=2 mm)" % H1)
    a("TH1    = 0.0")
    a("TH2    = 180.0")
    a("NDIV_R = %d" % nr)
    a("NDIV_TH = %d" % nt)
    a("NDIV_ZH = %d" % nzh)
    a("NDIV_ZP = %d" % nzp)
    a("NMODES = %d" % NMODES)
    a("")
    a("! ---- steel host (same TB,ANEL numbers as PLANE223 epsrel / 3-D nge1) ----")
    a("E_STEEL   = 2.00E5")
    a("NU_STEEL  = 0.3")
    a("RHO_STEEL = 7.800E-9")
    a("DS11 = E_STEEL*(1-NU_STEEL)/((1+NU_STEEL)*(1-2*NU_STEEL))")
    a("DS12 = E_STEEL*NU_STEEL/((1+NU_STEEL)*(1-2*NU_STEEL))")
    a("GS   = E_STEEL/(2*(1+NU_STEEL))")
    a("X_INERT = 1.0")
    a("")
    a("! ---- PZT-4, Liu 2002 Table 1 (NOT Duan C44E=73 GPa) ----")
    a("C11E = 1.320E5")
    a("C12E = 7.100E4")
    a("C13E = 7.300E4")
    a("C33E = 1.150E5")
    a("C44E = 2.600E4          ! Liu C55E = 26 GPa (yz/xz shear)")
    a("C66E = (C11E-C12E)/2")
    a("E31  = 4.1E-3           ! +4.1 project convention (Liu print -4.1)")
    a("E33  = 1.41E-2")
    a("E15  = 1.05E-2")
    a("EPS0   = 8.854E-12")
    a("X11ABS = 7.124E-9")
    a("X33ABS = 5.841E-9")
    a("X11    = X11ABS/EPS0     ! MP,PERx RELATIVE (epsrel pitfall)")
    a("X33    = X33ABS/EPS0")
    a("RHO_PZT = 7.500E-9")
    a("")
    a("/PREP7")
    a("ET,1,SOLID226")
    a("KEYOPT,1,1,1001")
    a("")
    a("! Material 1: steel. ANSYS order x, y, z, xy, yz, xz.")
    a("TB,ANEL,1,1,,0")
    a("TBDATA,1,DS11,DS12,DS12,0,0,0")
    a("TBDATA,7,DS11,DS12,0,0,0")
    a("TBDATA,12,DS11,0,0,0")
    a("TBDATA,16,GS,0,0")
    a("TBDATA,19,GS,0")
    a("TBDATA,21,GS")
    a("TB,PIEZ,1")
    a("TBDATA,1,0,0,0,0,0,0")
    a("TBDATA,7,0,0,0,0,0,0")
    a("TBDATA,13,0,0,0,0,0,0")
    a("MP,PERX,1,X_INERT")
    a("MP,PERY,1,X_INERT")
    a("MP,PERZ,1,X_INERT")
    a("MP,DENS,1,RHO_STEEL")
    a("")
    a("! Material 2: PZT-4, poled +Z. TB,PIEZ stress-major 6x3 (not field-major).")
    a("TB,ANEL,2,1,,0")
    a("TBDATA,1,C11E,C12E,C13E,0,0,0")
    a("TBDATA,7,C11E,C13E,0,0,0")
    a("TBDATA,12,C33E,0,0,0")
    a("TBDATA,16,C66E,0,0")
    a("TBDATA,19,C44E,0")
    a("TBDATA,21,C44E")
    a("TB,PIEZ,2")
    a(PIEZ_ON)
    a("MP,PERX,2,X11")
    a("MP,PERY,2,X11")
    a("MP,PERZ,2,X33")
    a("MP,DENS,2,RHO_PZT")
    a("")
    a("! ---- geometry: three stacked half-annuli (RI_MESH), glued ----")
    a("CSYS,0")
    a("CYLIND,RI,RO,-H,H,TH1,TH2")
    a("CYLIND,RI,RO,H,H+H1,TH1,TH2")
    a("CYLIND,RI,RO,-H-H1,-H,TH1,TH2")
    a("VGLUE,ALL")
    a("*GET,NV,VOLU,0,COUNT")
    a("*IF,NV,NE,3,THEN")
    a("  *MSG,FATAL")
    a("VGLUE did not leave exactly three volumes.")
    a("*ENDIF")
    a("*GET,NLT,LINE,0,COUNT")
    a("*IF,NLT,NE,28,THEN")
    a("  *MSG,FATAL")
    a("Glued half-annulus model does not have 28 lines.")
    a("*ENDIF")
    for lab, zc, mat in (("host", 0.0, 1),
                         ("top skin", H + h1 / 2.0, 2),
                         ("bottom skin", -H - h1 / 2.0, 2)):
        a("VSEL,S,LOC,Z,%.6f,%.6f" % (zc - ltol, zc + ltol))
        a("*GET,NV,VOLU,0,COUNT")
        a("*IF,NV,NE,1,THEN")
        a("  *MSG,FATAL")
        a("%s volume select did not return one volume." % lab)
        a("*ENDIF")
        a("VATT,%d,,1" % mat)
        a("ALLSEL,ALL")
    a("")
    a("! ---- line divisions (CSYS,1: LOC X is r). Parent annular pairing. ----")
    a("! Opposite edges of each 4-sided face: radials=NDIV_R, arcs=NDIV_TH.")
    a("CSYS,1")
    for k, z in enumerate(zl, start=1):
        a(line_block("radial lines z%d" % k, RMID - 1.0, RMID + 1.0,
                     z - ltol, z + ltol, 2, nr).rstrip("\n"))
        a(line_block("inner arc z%d" % k, RI - 0.3, RI + 0.3,
                     z - ltol, z + ltol, 1, nt).rstrip("\n"))
        a(line_block("outer arc z%d" % k, RO - 0.3, RO + 0.3,
                     z - ltol, z + ltol, 1, nt).rstrip("\n"))
    for lab, zc, nd in layers:
        a(line_block("inner thickness lines, %s" % lab, RI - 0.3, RI + 0.3,
                     zc - ltol, zc + ltol, 2, nd).rstrip("\n"))
        a(line_block("outer thickness lines, %s" % lab, RO - 0.3, RO + 0.3,
                     zc - ltol, zc + ltol, 2, nd).rstrip("\n"))
    a("CSYS,0")
    a("ALLSEL,ALL")
    a("MSHAPE,0,3D")
    a("MSHKEY,1")
    a("VMESH,ALL")
    a("ALLSEL,ALL")
    a("NUMMRG,NODES,1.0E-5")
    a("*GET,NEL,ELEM,0,COUNT")
    a("*IF,NEL,NE,%d,THEN" % nel)
    a("  *MSG,FATAL")
    a("Mapped mesh element count does not equal the division product.")
    a("*ENDIF")
    a("ESEL,S,MAT,,2")
    a("*GET,NELP,ELEM,0,COUNT")
    a("*IF,NELP,NE,%d,THEN" % (nr * nt * 2 * nzp))
    a("  *MSG,FATAL")
    a("Piezo-skin element count is wrong.")
    a("*ENDIF")
    a("ALLSEL,ALL")
    a("*GET,NNODE,NODE,0,COUNT")
    a("FINISH")
    a("")
    a("/SOLU")
    a("ANTYPE,MODAL")
    a("MODOPT,LANB,NMODES,-1.0,,,ON")
    a("MXPAND,NMODES,,,NO")
    a("CSYS,0")
    a("! symmetry plane y=0 (theta=0 and theta=180 faces): UY=0 only")
    a("NSEL,S,LOC,Y,-0.001,0.001")
    a("*GET,NSYM,NODE,0,COUNT")
    a("*IF,NSYM,NE,%d,THEN" % n_sym)
    a("  *MSG,FATAL")
    a("Symmetry-plane node count is wrong.")
    a("*ENDIF")
    a("D,ALL,UY,0")
    a("ALLSEL,ALL")
    a("! inner hole r=RI_MESH: FREE in v2 (mesh count asserted, NO constraint)")
    a("CSYS,1")
    a("NSEL,S,LOC,X,%.6f,%.6f" % (RI - 0.01, RI + 0.01))
    a("CSYS,0")
    a("*GET,NIN,NODE,0,COUNT")
    a("*IF,NIN,NE,%d,THEN" % n_inner)
    a("  *MSG,FATAL")
    a("Inner-hole node count is wrong.")
    a("*ENDIF")
    a("ALLSEL,ALL")
    a("! ---- outer mechanical BC ----")
    if bc == "C":
        a("! Outer clamped: all U = 0 on the outer cylindrical surface (r=RO).")
        a("CSYS,1")
        a("NSEL,S,LOC,X,%.6f,%.6f" % (RO - 0.01, RO + 0.01))
        a("CSYS,0")
        a("*GET,NCL,NODE,0,COUNT")
        a("*IF,NCL,LT,10,THEN")
        a("  *MSG,FATAL")
        a("Outer clamp node select returned fewer than 10 nodes.")
        a("*ENDIF")
        a("D,ALL,UX,0")
        a("D,ALL,UY,0")
        a("D,ALL,UZ,0")
        a("ALLSEL,ALL")
    else:
        a("! Outer simply-supported Kirchhoff analogue (classical w=Mr=0):")
        a("! UZ=0 on the outer circumference at midplane z=0 only.")
        a("! Mr is left natural (no moment/rotation constraint). Full-thickness")
        a("! UZ=0 on the outer cylinder would be a stiffer 3-D edge fixity;")
        a("! midplane-only is the closer thin-plate SS analogue used here.")
        a("! No project disk/SS SOLID226 precedent was found; this choice is")
        a("! documented explicitly rather than invented silently.")
        a("CSYS,1")
        a("NSEL,S,LOC,X,%.6f,%.6f" % (RO - 0.01, RO + 0.01))
        a("CSYS,0")
        a("NSEL,R,LOC,Z,%.6f,%.6f" % (-ztol, ztol))
        a("*GET,NSS,NODE,0,COUNT")
        a("*IF,NSS,LT,5,THEN")
        a("  *MSG,FATAL")
        a("Outer SS midplane ring select returned fewer than 5 nodes.")
        a("*ENDIF")
        a("D,ALL,UZ,0")
        a("ALLSEL,ALL")
        a("! v2: rigid x-translation pin -- UX=0 at the single outer midplane")
        a("! node at theta=90 deg (x=0, y=RO, z=0). Zero midplane in-plane")
        a("! displacement in symmetric-laminate flexure: bending untouched.")
        a("NSEL,S,LOC,X,-0.01,0.01")
        a("NSEL,R,LOC,Y,%.6f,%.6f" % (RO - 0.01, RO + 0.01))
        a("NSEL,R,LOC,Z,%.6f,%.6f" % (-ztol, ztol))
        a("*GET,NPIN,NODE,0,COUNT")
        a("*IF,NPIN,NE,1,THEN")
        a("  *MSG,FATAL")
        a("Rigid-x pin select is not exactly one node.")
        a("*ENDIF")
        a("D,ALL,UX,0")
        a("ALLSEL,ALL")
    a("! ---- electrical BC: both-faces SC (scboth / nge1 sc pattern) ----")
    a("! VOLT=0 on host/piezo interfaces AND on both outer faces. No CP.")
    a("NSEL,S,LOC,Z,%.6f,%.6f" % (H - ztol, H + ztol))
    a("NSEL,A,LOC,Z,%.6f,%.6f" % (-H - ztol, -H + ztol))
    a("*GET,NINT,NODE,0,COUNT")
    a("*IF,NINT,NE,%d,THEN" % (2 * n_elec_face))
    a("  *MSG,FATAL")
    a("Interface electrode node count is wrong.")
    a("*ENDIF")
    a("D,ALL,VOLT,0")
    a("ALLSEL,ALL")
    a("NSEL,S,LOC,Z,%.6f,%.6f" % (H + h1 - ztol, H + h1 + ztol))
    a("NSEL,A,LOC,Z,%.6f,%.6f" % (-H - h1 - ztol, -H - h1 + ztol))
    a("*GET,NOUT,NODE,0,COUNT")
    a("*IF,NOUT,NE,%d,THEN" % (2 * n_elec_face))
    a("  *MSG,FATAL")
    a("Outer electrode node count is wrong.")
    a("*ENDIF")
    a("D,ALL,VOLT,0")
    a("ALLSEL,ALL")
    a("SOLVE")
    a("FINISH")
    a("")
    a("/POST1")
    a("*DIM,FR1,ARRAY,NMODES")
    a("*DIM,IDXARR,ARRAY,NMODES")
    a("*DO,I,1,NMODES")
    a("  *GET,FR1(I),MODE,I,FREQ")
    a("  IDXARR(I)=I")
    a("*ENDDO")
    # CFOPEN stem without date suffix (match nge1 pattern)
    modes_tag = s[: -len("_" + DATE)]
    a("*CFOPEN,%s_modes,txt" % modes_tag)
    a("*VWRITE")
    a("('tag liu_disk_sc_%s')" % bc)
    a("*VWRITE,NEL")
    a("('elements ', F12.0)")
    a("*VWRITE,NNODE")
    a("('nodes ', F12.0)")
    a("*VWRITE,NSYM")
    a("('symmetry_nodes ', F12.0)")
    a("*VWRITE,NINT")
    a("('interface_nodes ', F12.0)")
    a("*VWRITE,NOUT")
    a("('outer_nodes ', F12.0)")
    a("*VWRITE")
    a("('mode f_Hz')")
    a("*VWRITE,IDXARR(1),FR1(1)")
    a("(F4.0, 3X, E20.12)")
    a("*CFCLOS")
    a("")
    a("! ---- sample nodes: located once, before any mode loop ----")
    a("NSAMP = %d" % nsamp)
    a("*DIM,SX,ARRAY,NSAMP")
    a("*DIM,SY,ARRAY,NSAMP")
    a("*DIM,SZ,ARRAY,NSAMP")
    a("*DIM,SR,ARRAY,NSAMP")
    a("*DIM,ST,ARRAY,NSAMP")
    a("*DIM,SKI,ARRAY,NSAMP")
    a("*DIM,NS,ARRAY,NSAMP")
    for k, (r, th, z) in enumerate(smp, start=1):
        x = r * math.cos(math.radians(th))
        y = r * math.sin(math.radians(th))
        if abs(y) < 1e-9:
            y = 0.0
        if abs(x) < 1e-9:
            x = 0.0
        a("SX(%d) = %.9f" % (k, x))
        a("SY(%d) = %.9f" % (k, y))
        a("SZ(%d) = %.9f" % (k, z))
        a("SR(%d) = %.6f" % (k, r))
        a("ST(%d) = %.3f" % (k, th))
        a("SKI(%d) = %d" % (k, k))
    a("STOL = 0.01")
    a("CSYS,0")
    a("*DO,K,1,NSAMP")
    a("  X1 = SX(K)-STOL")
    a("  X2 = SX(K)+STOL")
    a("  Y1 = SY(K)-STOL")
    a("  Y2 = SY(K)+STOL")
    a("  Z1 = SZ(K)-STOL")
    a("  Z2 = SZ(K)+STOL")
    a("  NSEL,S,LOC,X,X1,X2")
    a("  NSEL,R,LOC,Y,Y1,Y2")
    a("  NSEL,R,LOC,Z,Z1,Z2")
    a("  *GET,NC,NODE,0,COUNT")
    a("  *IF,NC,NE,1,THEN")
    a("    *MSG,FATAL")
    a("Sample node select did not return exactly one node.")
    a("  *ENDIF")
    a("  *GET,NN,NODE,0,NUM,MIN")
    a("  NS(K) = NN")
    a("*ENDDO")
    a("ALLSEL,ALL")
    a("*CFOPEN,%s_samplecoords,txt" % modes_tag)
    a("*VWRITE")
    a("('k r_mm theta_deg z_mm node')")
    a("*VWRITE,SKI(1),SR(1),ST(1),SZ(1),NS(1)")
    a("(F5.0, 3X, F10.3, 3X, F10.3, 3X, F10.4, 3X, F10.0)")
    a("*CFCLOS")
    a("")
    a("! ---- raw nodal values, every mode x every sample, flat arrays ----")
    a("NTOT = NMODES*NSAMP")
    a("*DIM,OMI,ARRAY,NTOT")
    a("*DIM,OKI,ARRAY,NTOT")
    a("*DIM,OUX,ARRAY,NTOT")
    a("*DIM,OUY,ARRAY,NTOT")
    a("*DIM,OUZ,ARRAY,NTOT")
    a("*DIM,OVV,ARRAY,NTOT")
    a("*DO,I,1,NMODES")
    a("  SET,1,I")
    a("  *DO,K,1,NSAMP")
    a("    J = (I-1)*NSAMP+K")
    a("    NN = NS(K)")
    a("    *GET,GVX,NODE,NN,U,X")
    a("    *GET,GVY,NODE,NN,U,Y")
    a("    *GET,GVZ,NODE,NN,U,Z")
    a("    *GET,GVV,NODE,NN,VOLT")
    a("    OMI(J) = I")
    a("    OKI(J) = K")
    a("    OUX(J) = GVX")
    a("    OUY(J) = GVY")
    a("    OUZ(J) = GVZ")
    a("    OVV(J) = GVV")
    a("  *ENDDO")
    a("*ENDDO")
    a("*CFOPEN,%s_samples,txt" % modes_tag)
    a("*VWRITE")
    a("('mode k UX UY UZ VOLT')")
    a("*VWRITE,OMI(1),OKI(1),OUX(1),OUY(1),OUZ(1),OVV(1)")
    a("(F4.0, 1X, F5.0, 1X, E16.8, 1X, E16.8, 1X, E16.8, 1X, E16.8)")
    a("*CFCLOS")
    a("FINISH")
    txt = "\n".join(L) + "\n"
    assert "***" not in txt, "MAPDL failure-token or *** leaked into deck"
    path = os.path.join(HERE, s + ".inp")
    with open(path, "w", newline="\n") as f:
        f.write(txt)
    return path, dict(
        elements=nel, sym=n_sym, elec_face=n_elec_face,
        inner=n_inner, nsamp=nsamp, job=jobname(path),
        RI_mesh=RI, RI_liu=RI_LIU,
    )


def main():
    names = []
    jobs = {}
    # collide against sibling inp in this folder if any
    for fn in os.listdir(HERE):
        if fn.startswith("ansys_") and fn.endswith(".inp"):
            jobs[jobname(fn)] = fn
    for _tag, bc in DECKS:
        p, info = build(bc)
        j = info["job"]
        if j in jobs and jobs[j] != os.path.basename(p):
            raise SystemExit("FATAL jobname collision %r with %s" % (j, jobs[j]))
        jobs[j] = os.path.basename(p)
        names.append(os.path.basename(p))
        print("%-52s job=%-24s %s" % (os.path.basename(p), j, info))
    return names


if __name__ == "__main__":
    main()
