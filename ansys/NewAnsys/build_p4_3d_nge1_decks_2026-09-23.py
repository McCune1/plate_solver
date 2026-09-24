#!/usr/bin/env python3
"""Build the Paper 4 3-D FE check of coupled n >= 1 roots (SOLID226 half ring).

FUTURE_WORK.md, Paper 4: "(Optional) 3-D FE check of a coupled n >= 1
root. None is published, so no paper claim depends on it."

Every Paper 4 FE deck so far is axisymmetric PLANE223, so it can only
see n = 0. This builder writes 3-D SOLID226 decks of the same layered
F-F ring (steel host |z| <= h, PZT-4 skins, same poling +Z, Duan 2005
geometry and constants, mm-N-s-tonne, relative permittivity) so that the
first flexural roots of n = 0..5 all appear in one modal solve.

HALF MODEL. theta = 0..180 deg with UY = 0 on the y = 0 plane (both cut
faces). That is the symmetry condition for the cos(n theta) family of
w = W(r) cos(n theta): u_theta ~ sin(n theta) vanishes on the cut. It
keeps one member of every degenerate pair, for every n, including n = 0.
VOLT is left natural on the cut (symmetric). Rigid modes left: X and Z
translation and rotation about Y (three near-zero modes).

MATERIAL BLOCKS. The Z-poled SOLID226 TB,ANEL / TB,PIEZ / MP,PER* layout
is copied from build_p6_sector_fe_decks_2026-09-22.py (FE-confirmed:
Paper 6 Gate E S3b PASS, job 2528994). The steel block and the inert
host permittivity (relative 1.0) are the PLANE223 epsrel decks'
(ansys_p4_epsrel_ff_lanb_*_2026-09-16.inp), with the same TB,ANEL
numbers. Units block is those decks', not re-derived.

ELECTRICAL CASES (same names as the axisymmetric decks):
  e0  TB,PIEZ,2 all zero; VOLT = 0 on interfaces and outers (as sc).
  sc  both-faces short circuit: VOLT = 0 on z = +-h and z = +-(h+h1).
  oc  VOLT = 0 on z = +-h; both outer faces one CP,VOLT bus, unconstrained.

MESHES. coarse 24 x 72 x (2 + 1 + 1), fine 40 x 144 x (4 + 2 + 2)
(radial x circumferential x through-thickness host + each skin). The fine
mesh is about 2.1e5 nodes, 8.3e5 DOF (Paper 6's fine SOLID226 deck was
1.6e6 DOF and solved in 49 s with 12 GB scratch).

DECKS (10). h1/2h = 1/12: e0/sc coarse, e0/sc/oc fine. h1/2h = 1/5: the
same five. The coarse pairs exist only to show that the fine-mesh
frequencies AND the tiny sc/e0 splits are mesh-stable.

OUTPUT. Raw Hz and raw nodal UX, UY, UZ, VOLT only. Exactly one *VWRITE
per output file, outside every *DO loop, on whole arrays (project memory
ansys-kfac-runtime-bug; the 2026-09-22 *VWRITE-in-*DO garble). The
scorer does every ratio.

Deck text is LF. The MAPDL failure token is never written into a deck,
so the queue's error grep of the echoed input cannot false-fail.
"""
# (no __future__ import: the compute nodes run python3 3.6, job 2530336)

import math
import os

HERE = os.path.dirname(os.path.abspath(__file__))
DATE = "2026-09-23"

RI, RO, H = 100.0, 600.0, 10.0
RMID = 0.5 * (RI + RO)
NMODES = 24

MESHES = {
    "coarse": dict(NDIV_R=24, NDIV_TH=72, NDIV_ZH=2, NDIV_ZP=1),
    "fine": dict(NDIV_R=40, NDIV_TH=144, NDIV_ZH=4, NDIV_ZP=2),
}

# (ratio tag, h1/2h denominator, case, mesh)
DECKS = [
    ("h112", 12, "e0", "coarse"),
    ("h112", 12, "sc", "coarse"),
    ("h112", 12, "e0", "fine"),
    ("h112", 12, "sc", "fine"),
    ("h112", 12, "oc", "fine"),
    ("h15", 5, "e0", "coarse"),
    ("h15", 5, "sc", "coarse"),
    ("h15", 5, "e0", "fine"),
    ("h15", 5, "sc", "fine"),
    ("h15", 5, "oc", "fine"),
]

SAMPLE_RADII = (RI, 225.0, RMID, 475.0, RO)
SAMPLE_THETAS = tuple(15.0 * k for k in range(13))  # 0..180

PIEZ_ON = """TBDATA,1,0,0,E31,0,0,E31
TBDATA,7,0,0,E33,0,0,0
TBDATA,13,0,E15,0,E15,0,0"""
PIEZ_OFF = """TBDATA,1,0,0,0,0,0,0
TBDATA,7,0,0,0,0,0,0
TBDATA,13,0,0,0,0,0,0"""

CASE_TEXT = {
    "e0": "coupling off (TB,PIEZ,2 zero), VOLT=0 on interfaces and outers",
    "sc": "both-faces short circuit, VOLT=0 on interfaces and outers",
    "oc": "open circuit, VOLT=0 on interfaces, outers one CP,VOLT bus",
}


def stem(tag, case, mesh):
    return "ansys_p4_3d_%s_%s_%s_%s" % (tag, case, mesh, DATE)


def jobname(fname):
    """Same rule as run_ansys_queue.sh: strip ansys_, non-alnum -> _, 24 chars."""
    s = os.path.basename(fname)
    if s.endswith(".inp"):
        s = s[:-4]
    if s.startswith("ansys_"):
        s = s[6:]
    s = "".join(c if c.isalnum() else "_" for c in s)
    return s[:24]


def face_nodes(n1, n2):
    """Nodes on one mapped face of 8-node quads, n1 x n2 elements."""
    return (n1 + 1) * (n2 + 1) + n1 * (n2 + 1) + (n1 + 1) * n2


def samples(h1):
    zt = H + h1
    out = []
    for r in SAMPLE_RADII:
        for th in SAMPLE_THETAS:
            out.append((r, th, zt))
    for th in SAMPLE_THETAS:
        out.append((RO, th, -zt))
    return out


def line_block(label, xlo, xhi, zlo, zhi, expect, ndiv):
    return """! {label}
LSEL,S,LOC,X,{xlo:.6f},{xhi:.6f}
LSEL,R,LOC,Z,{zlo:.6f},{zhi:.6f}
*GET,NL,LINE,0,COUNT
*IF,NL,NE,{expect},THEN
  *MSG,FATAL
{label}: line count was not {expect}.
*ENDIF
LESIZE,ALL,,,{ndiv}
ALLSEL,ALL
""".format(label=label, xlo=xlo, xhi=xhi, zlo=zlo, zhi=zhi,
           expect=expect, ndiv=ndiv)


def build(tag, den, case, mesh):
    m = MESHES[mesh]
    nr, nt, nzh, nzp = m["NDIV_R"], m["NDIV_TH"], m["NDIV_ZH"], m["NDIV_ZP"]
    nz = nzh + 2 * nzp
    h1 = 2.0 * H / den
    nel = nr * nt * nz
    n_elec_face = face_nodes(nr, nt)
    n_sym = 2 * face_nodes(nr, nz)
    s = stem(tag, case, mesh)
    smp = samples(h1)
    nsamp = len(smp)
    ztol = min(0.01, h1 / 20.0)
    ltol = h1 / 10.0

    # sanity: every sample must be a node of this mesh
    dr = (RO - RI) / nr
    for r, th, _ in smp:
        kr = (r - RI) / (dr / 2.0)
        kt = th / (180.0 / nt / 2.0)
        assert abs(kr - round(kr)) < 1e-9, (mesh, r)
        assert abs(kt - round(kt)) < 1e-9, (mesh, th)
        # theta must be a corner line: a (midside r, midside theta) point
        # would be a face centre, which a 20-node brick has no node at.
        assert int(round(kt)) % 2 == 0, (mesh, th)

    zl = [-H - h1, -H, H, H + h1]
    layers = [("bottom skin", -H - h1 / 2.0, nzp),
              ("host", 0.0, nzh),
              ("top skin", H + h1 / 2.0, nzp)]

    L = []
    a = L.append
    a("! %s.inp" % s)
    a("! Paper 4 -- 3-D FE check of coupled n >= 1 F-F roots (FUTURE_WORK.md).")
    a("! Built by build_p4_3d_nge1_decks_%s.py. Do not hand-edit; rebuild." % DATE)
    a("! Layered F-F annulus, SOLID226 (KEYOPT(1)=1001), HALF model theta=0..180,")
    a("! UY=0 on y=0 (cos(n theta) family). h1/2h = 1/%d. Case: %s." % (den, CASE_TEXT[case]))
    a("! Mesh %s: %d radial x %d circ x (%d host + %d + %d skins) = %d elements." % (
        mesh, nr, nt, nzh, nzp, nzp, nel))
    a("! Units mm-N-s-tonne, real volts, RELATIVE permittivity (epsrel decks).")
    a("! Poling +Z in both skins (same-poling parallel bimorph, paper Eq. Mp).")
    a("! Prints raw Hz and raw UX/UY/UZ/VOLT at %d sample nodes only." % nsamp)
    a("! Score with score_p4_3d_nge1_%s.py. Never trust .QUEUE_OK alone." % DATE)
    a("! ============================================================")
    a("FINISH")
    a("/CLEAR,NOSTART")
    a("/TITLE, Paper 4 3-D half ring SOLID226 %s %s h1/2h=1/%d (%s)" % (case, mesh, den, DATE))
    a("")
    a("RI     = %.1f" % RI)
    a("RO     = %.1f" % RO)
    a("H      = %.1f           ! host HALF-thickness [mm]" % H)
    a("H1     = 2*H/%d          ! skin thickness [mm]" % den)
    a("TH1    = 0.0")
    a("TH2    = 180.0")
    a("NDIV_R = %d" % nr)
    a("NDIV_TH = %d" % nt)
    a("NDIV_ZH = %d" % nzh)
    a("NDIV_ZP = %d" % nzp)
    a("NMODES = %d" % NMODES)
    a("")
    a("! ---- steel host (TB,ANEL numbers as the PLANE223 epsrel decks) ----")
    a("E_STEEL   = 2.00E5")
    a("NU_STEEL  = 0.3")
    a("RHO_STEEL = 7.800E-9")
    a("DS11 = E_STEEL*(1-NU_STEEL)/((1+NU_STEEL)*(1-2*NU_STEEL))")
    a("DS12 = E_STEEL*NU_STEEL/((1+NU_STEEL)*(1-2*NU_STEEL))")
    a("GS   = E_STEEL/(2*(1+NU_STEEL))")
    a("X_INERT = 1.0")
    a("")
    a("! ---- PZT-4, Duan 2005 Table 1 (as every Paper 4 deck) ----")
    a("C11E = 1.320E5")
    a("C12E = 7.100E4")
    a("C13E = 7.300E4")
    a("C33E = 1.150E5")
    a("C44E = 7.300E4")
    a("C66E = (C11E-C12E)/2")
    a("E31  = 4.1E-3")
    a("E33  = 1.41E-2")
    a("E15  = 1.05E-2")
    a("EPS0   = 8.854E-12")
    a("X11ABS = 7.124E-9")
    a("X33ABS = 5.841E-9")
    a("X11    = X11ABS/EPS0")
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
    a(PIEZ_OFF)
    a("MP,PERX,1,X_INERT")
    a("MP,PERY,1,X_INERT")
    a("MP,PERZ,1,X_INERT")
    a("MP,DENS,1,RHO_STEEL")
    a("")
    a("! Material 2: PZT-4, poled +Z (Paper 6 SOLID226 block).")
    a("TB,ANEL,2,1,,0")
    a("TBDATA,1,C11E,C12E,C13E,0,0,0")
    a("TBDATA,7,C11E,C13E,0,0,0")
    a("TBDATA,12,C33E,0,0,0")
    a("TBDATA,16,C66E,0,0")
    a("TBDATA,19,C44E,0")
    a("TBDATA,21,C44E")
    a("TB,PIEZ,2")
    a(PIEZ_OFF if case == "e0" else PIEZ_ON)
    a("MP,PERX,2,X11")
    a("MP,PERY,2,X11")
    a("MP,PERZ,2,X33")
    a("MP,DENS,2,RHO_PZT")
    a("")
    a("! ---- geometry: three stacked half-annulus volumes, glued ----")
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
    a("Glued model does not have 28 lines.")
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
    a("! ---- line divisions (CSYS,1: LOC X is r) ----")
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
    a("! inner electrodes (host/skin interfaces) grounded in every case")
    a("NSEL,S,LOC,Z,%.6f,%.6f" % (H - ztol, H + ztol))
    a("NSEL,A,LOC,Z,%.6f,%.6f" % (-H - ztol, -H + ztol))
    a("*GET,NINT,NODE,0,COUNT")
    a("*IF,NINT,NE,%d,THEN" % (2 * n_elec_face))
    a("  *MSG,FATAL")
    a("Interface electrode node count is wrong.")
    a("*ENDIF")
    a("D,ALL,VOLT,0")
    a("ALLSEL,ALL")
    a("! outer electrodes z = +-(H+H1)")
    a("NSEL,S,LOC,Z,%.6f,%.6f" % (H + h1 - ztol, H + h1 + ztol))
    a("NSEL,A,LOC,Z,%.6f,%.6f" % (-H - h1 - ztol, -H - h1 + ztol))
    a("*GET,NOUT,NODE,0,COUNT")
    a("*IF,NOUT,NE,%d,THEN" % (2 * n_elec_face))
    a("  *MSG,FATAL")
    a("Outer electrode node count is wrong.")
    a("*ENDIF")
    if case == "oc":
        a("CP,1,VOLT,ALL        ! one bus, both outer faces, unconstrained (Q=0)")
    else:
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
    a("*CFOPEN,%s_modes,txt" % s[:-len("_" + DATE)])
    a("*VWRITE")
    a("('tag %s_%s_%s')" % (tag, case, mesh))
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
        a("SR(%d) = %.3f" % (k, r))
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
    a("*CFOPEN,%s_samplecoords,txt" % s[:-len("_" + DATE)])
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
    a("*CFOPEN,%s_samples,txt" % s[:-len("_" + DATE)])
    a("*VWRITE")
    a("('mode k UX UY UZ VOLT')")
    a("*VWRITE,OMI(1),OKI(1),OUX(1),OUY(1),OUZ(1),OVV(1)")
    a("(F4.0, 1X, F5.0, 1X, E16.8, 1X, E16.8, 1X, E16.8, 1X, E16.8)")
    a("*CFCLOS")
    a("FINISH")
    txt = "\n".join(L) + "\n"
    assert "***" not in txt
    path = os.path.join(HERE, s + ".inp")
    with open(path, "w", newline="\n") as f:
        f.write(txt)
    return path, dict(elements=nel, sym=n_sym, elec_face=n_elec_face,
                      nsamp=nsamp)


def main():
    names = []
    jobs = set()
    for tag, den, case, mesh in DECKS:
        p, info = build(tag, den, case, mesh)
        j = jobname(p)
        assert j not in jobs, j
        jobs.add(j)
        names.append(os.path.basename(p))
        print("%-52s job=%-24s %s" % (os.path.basename(p), j, info))
    return names


if __name__ == "__main__":
    main()
