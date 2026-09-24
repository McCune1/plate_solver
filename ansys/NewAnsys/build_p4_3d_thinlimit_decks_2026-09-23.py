#!/usr/bin/env python3
"""Build the Paper 4 3-D thin-limit decks (bare steel half ring, t = 20/10/5 mm).

Purpose: settle the G3 FAIL of job 2530423 (LESSONS 18.232). There the
elastic n >= 1 FE roots of the layered ring sat 0.5-1.2% below the
Kirchhoff elastic-bilayer 4x4, beyond the 1% bar at three points. Two
readings:
  (a) thick-plate physics (transverse shear / rotary inertia in the
      interior, plus the free-edge boundary layer that a Kirchhoff
      effective-shear edge omits). Then FE/Kirchhoff - 1 -> 0 as t -> 0
      at every n, the n = 0 part ~t^2 and the n >= 1 edge part ~t.
  (b) an n-dependent defect in the thin-plate operator. Then the
      relative offset survives as t -> 0.
A bare steel plate at three thicknesses separates them: Kirchhoff roots
scale exactly with t, so any t-independent residue is an operator error.

Model: same half-ring topology as build_p4_3d_nge1_decks_2026-09-23.py
(theta = 0..180, UY = 0 on y = 0, cos(n theta) family), one steel volume,
SOLID186 with KEYOPT(2) = 1 (full integration, matching SOLID226's
structural integration) and MP,EX / NUXY / DENS in mm-N-s-tonne.
No VOLT dof exists, so the samples file carries a VOLT column of zeros
(the nge1 scorer's parser expects six columns).

Meshes: fine 48 x 144 x 4 at t = 20, 10, 5; coarse 24 x 72 x 2 at t = 20
and 5 (job 2530553). Round 2 adds coarse at t = 10 and xfine 96 x 288 x 8
at all three t (a ratio-2 sequence at every thickness). Samples, count checks, one *VWRITE per file outside every *DO: as
the nge1 decks. Output LF; the MAPDL failure token never appears.
"""
import importlib.util
import math
import os

HERE = os.path.dirname(os.path.abspath(__file__))
DATE = "2026-09-23"


def _nge1():
    p = os.path.join(HERE, "build_p4_3d_nge1_decks_%s.py" % DATE)
    spec = importlib.util.spec_from_file_location("nge1_builder", p)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


NG = _nge1()
RI, RO, RMID = NG.RI, NG.RO, NG.RMID
NMODES = 24
MESHES = {
    "coarse": dict(NDIV_R=24, NDIV_TH=72, NDIV_Z=2),
    "fine": dict(NDIV_R=48, NDIV_TH=144, NDIV_Z=4),
    # round 2 (2026-09-23, after job 2530553 left T3 UNRESOLVED on its
    # t = 5 mesh bar): one more uniform halving in every direction, so
    # coarse / fine / xfine is a ratio-2 sequence for Richardson.
    "xfine": dict(NDIV_R=96, NDIV_TH=288, NDIV_Z=8),
}
DECKS = [(20, "coarse"), (20, "fine"), (10, "fine"), (5, "fine"), (5, "coarse")]
# Round 2: the missing t = 10 coarse deck plus xfine at all three t.
DECKS_R2 = [(10, "coarse"), (20, "xfine"), (10, "xfine"), (5, "xfine")]


def stem(t, mesh):
    return "ansys_p4_3d_bare_t%d_%s_%s" % (t, mesh, DATE)


def samples(t):
    zt = t / 2.0
    out = [(r, th, zt) for r in NG.SAMPLE_RADII for th in NG.SAMPLE_THETAS]
    out += [(RO, th, -zt) for th in NG.SAMPLE_THETAS]
    return out


def build(t, mesh):
    m = MESHES[mesh]
    nr, nt, nz = m["NDIV_R"], m["NDIV_TH"], m["NDIV_Z"]
    nel = nr * nt * nz
    n_sym = 2 * NG.face_nodes(nr, nz)
    s = stem(t, mesh)
    short = s[:-len("_" + DATE)]
    smp = samples(t)
    nsamp = len(smp)
    ltol = t / 10.0
    dr = (RO - RI) / nr
    for r, th, _ in smp:
        kr = (r - RI) / (dr / 2.0)
        kt = th / (180.0 / nt / 2.0)
        assert abs(kr - round(kr)) < 1e-9 and abs(kt - round(kt)) < 1e-9
        assert int(round(kt)) % 2 == 0

    L = []
    a = L.append
    a("! %s.inp" % s)
    a("! Paper 4 -- 3-D thin-limit check for the G3 FAIL of job 2530423.")
    a("! Built by build_p4_3d_thinlimit_decks_%s.py. Do not hand-edit; rebuild." % DATE)
    a("! Bare steel F-F half annulus, SOLID186 full integration, theta=0..180,")
    a("! UY=0 on y=0 (cos(n theta) family). Thickness t = %d mm." % t)
    a("! Mesh %s: %d radial x %d circ x %d thick = %d elements." % (mesh, nr, nt, nz, nel))
    a("! Units mm-N-s-tonne. Prints raw Hz and raw UX/UY/UZ at %d sample nodes." % nsamp)
    a("! Score with score_p4_3d_thinlimit_%s.py. Never trust .QUEUE_OK alone." % DATE)
    a("! ============================================================")
    a("FINISH")
    a("/CLEAR,NOSTART")
    a("/TITLE, Paper 4 thin-limit bare steel half ring t=%dmm %s (%s)" % (t, mesh, DATE))
    a("")
    a("RI     = %.1f" % RI)
    a("RO     = %.1f" % RO)
    a("T      = %.1f           ! full thickness [mm]" % t)
    a("TH1    = 0.0")
    a("TH2    = 180.0")
    a("NMODES = %d" % NMODES)
    a("E_STEEL   = 2.00E5")
    a("NU_STEEL  = 0.3")
    a("RHO_STEEL = 7.800E-9")
    a("")
    a("/PREP7")
    a("ET,1,SOLID186")
    a("KEYOPT,1,2,1        ! full integration")
    a("MP,EX,1,E_STEEL")
    a("MP,NUXY,1,NU_STEEL")
    a("MP,DENS,1,RHO_STEEL")
    a("CSYS,0")
    a("CYLIND,RI,RO,-T/2,T/2,TH1,TH2")
    a("*GET,NV,VOLU,0,COUNT")
    a("*IF,NV,NE,1,THEN")
    a("  *MSG,FATAL")
    a("CYLIND did not create exactly one volume.")
    a("*ENDIF")
    a("*GET,NLT,LINE,0,COUNT")
    a("*IF,NLT,NE,12,THEN")
    a("  *MSG,FATAL")
    a("Volume does not have 12 lines.")
    a("*ENDIF")
    a("VSEL,ALL")
    a("VATT,1,,1")
    a("ALLSEL,ALL")
    a("CSYS,1")
    for k, z in enumerate((-t / 2.0, t / 2.0), start=1):
        a(NG.line_block("radial lines z%d" % k, RMID - 1.0, RMID + 1.0,
                        z - ltol, z + ltol, 2, nr).rstrip("\n"))
        a(NG.line_block("inner arc z%d" % k, RI - 0.3, RI + 0.3,
                        z - ltol, z + ltol, 1, nt).rstrip("\n"))
        a(NG.line_block("outer arc z%d" % k, RO - 0.3, RO + 0.3,
                        z - ltol, z + ltol, 1, nt).rstrip("\n"))
    a(NG.line_block("inner thickness lines", RI - 0.3, RI + 0.3,
                    -ltol, ltol, 2, nz).rstrip("\n"))
    a(NG.line_block("outer thickness lines", RO - 0.3, RO + 0.3,
                    -ltol, ltol, 2, nz).rstrip("\n"))
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
    a("*GET,NNODE,NODE,0,COUNT")
    a("FINISH")
    a("")
    a("/SOLU")
    a("ANTYPE,MODAL")
    a("MODOPT,LANB,NMODES,-1.0,,,ON")
    a("MXPAND,NMODES,,,NO")
    a("CSYS,0")
    a("NSEL,S,LOC,Y,-0.001,0.001")
    a("*GET,NSYM,NODE,0,COUNT")
    a("*IF,NSYM,NE,%d,THEN" % n_sym)
    a("  *MSG,FATAL")
    a("Symmetry-plane node count is wrong.")
    a("*ENDIF")
    a("D,ALL,UY,0")
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
    a("*CFOPEN,%s_modes,txt" % short)
    a("*VWRITE")
    a("('tag bare_t%d_%s')" % (t, mesh))
    a("*VWRITE,NEL")
    a("('elements ', F12.0)")
    a("*VWRITE,NNODE")
    a("('nodes ', F12.0)")
    a("*VWRITE,NSYM")
    a("('symmetry_nodes ', F12.0)")
    a("*VWRITE")
    a("('mode f_Hz')")
    a("*VWRITE,IDXARR(1),FR1(1)")
    a("(F4.0, 3X, E20.12)")
    a("*CFCLOS")
    a("")
    a("NSAMP = %d" % nsamp)
    for arr in ("SX", "SY", "SZ", "SR", "ST", "SKI", "NS"):
        a("*DIM,%s,ARRAY,NSAMP" % arr)
    for k, (r, th, z) in enumerate(smp, start=1):
        x = r * math.cos(math.radians(th))
        y = r * math.sin(math.radians(th))
        x = 0.0 if abs(x) < 1e-9 else x
        y = 0.0 if abs(y) < 1e-9 else y
        a("SX(%d) = %.9f" % (k, x))
        a("SY(%d) = %.9f" % (k, y))
        a("SZ(%d) = %.9f" % (k, z))
        a("SR(%d) = %.3f" % (k, r))
        a("ST(%d) = %.3f" % (k, th))
        a("SKI(%d) = %d" % (k, k))
    a("STOL = %.4f" % min(0.01, t / 20.0))
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
    a("*CFOPEN,%s_samplecoords,txt" % short)
    a("*VWRITE")
    a("('k r_mm theta_deg z_mm node')")
    a("*VWRITE,SKI(1),SR(1),ST(1),SZ(1),NS(1)")
    a("(F5.0, 3X, F10.3, 3X, F10.3, 3X, F10.4, 3X, F10.0)")
    a("*CFCLOS")
    a("")
    a("NTOT = NMODES*NSAMP")
    for arr in ("OMI", "OKI", "OUX", "OUY", "OUZ", "OVV"):
        a("*DIM,%s,ARRAY,NTOT" % arr)
    a("*DO,I,1,NMODES")
    a("  SET,1,I")
    a("  *DO,K,1,NSAMP")
    a("    J = (I-1)*NSAMP+K")
    a("    NN = NS(K)")
    a("    *GET,GVX,NODE,NN,U,X")
    a("    *GET,GVY,NODE,NN,U,Y")
    a("    *GET,GVZ,NODE,NN,U,Z")
    a("    OMI(J) = I")
    a("    OKI(J) = K")
    a("    OUX(J) = GVX")
    a("    OUY(J) = GVY")
    a("    OUZ(J) = GVZ")
    a("    OVV(J) = 0")
    a("  *ENDDO")
    a("*ENDDO")
    a("*CFOPEN,%s_samples,txt" % short)
    a("*VWRITE")
    a("('mode k UX UY UZ VOLT(zero: no VOLT dof)')")
    a("*VWRITE,OMI(1),OKI(1),OUX(1),OUY(1),OUZ(1),OVV(1)")
    a("(F4.0, 1X, F5.0, 1X, E16.8, 1X, E16.8, 1X, E16.8, 1X, E16.8)")
    a("*CFCLOS")
    a("FINISH")
    txt = "\n".join(L) + "\n"
    assert "***" not in txt
    path = os.path.join(HERE, s + ".inp")
    with open(path, "w", newline="\n") as f:
        f.write(txt)
    return path, dict(elements=nel, sym=n_sym, nsamp=nsamp)


def main():
    jobs = set()
    for t, mesh in DECKS + DECKS_R2:
        p, info = build(t, mesh)
        j = NG.jobname(p)
        assert j not in jobs, j
        jobs.add(j)
        print("%-44s job=%-24s %s" % (os.path.basename(p), j, info))


if __name__ == "__main__":
    main()
