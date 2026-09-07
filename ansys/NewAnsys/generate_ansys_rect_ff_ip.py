# Generate rectangular FFFF IP half-model PLANE183 decks for l/b=1.0, 1.5,
# 2.0, 2.5, 3.0 (same aspect-ratio set as the OOP FE). Hz-only *VWRITE
# AFTER the *DO (annular ansys_ff_ip_ext.inp pattern) -- never *VWRITE
# inside a *DO, never a deck Omega column (APDL SQRT is untrusted).
from __future__ import annotations

from pathlib import Path

HERE = Path(__file__).resolve().parent
OFAC = 804.481  # Hz per unit Omega_bar; B=1, E=210e9, nu=0.30, rho=7800
# Omega = 4*B*f / sqrt(c66/rho), c66=E/(2(1+nu)); f[Hz] / OFAC.


def ndiv(lob):
    return int(round(2 * lob * 20)), int(round(2 * lob * 30))


def one_pass(tag, lob, title, ndivx, ndivy, par, mesh, arr):
    """par is SYMM or ASYM; mesh is 1 or 2."""
    if par == "SYMM":
        bc = "NSEL, S, LOC, Y, 0\nD, ALL, UY, 0"
        note = "SYMM=EVEN; UY=0 on Y=0 (solver SYM, u2 odd)"
    else:
        bc = "NSEL, S, LOC, Y, 0\nD, ALL, UX, 0"
        note = "ASYM=ODD; UX=0 on Y=0 (solver ANTI, u1 odd)"
    cf = f"rect_ff_ip_{par.lower()}_mesh{mesh}_{tag}"
    conv = (
        "convergence pass; must agree with mesh1 <0.3%"
        if mesh == 2
        else "discard ~0 Hz rigid"
    )
    return f"""! ======================= PASS {par}, mesh{mesh} =======================
FINISH
/CLEAR,NOSTART
/TITLE, {title}

B      = 1.0
LOB    = {lob:.1f}
L      = LOB*B
H      = 0.04
E_iso  = 210.0e9
NU_iso = 0.30
RHO    = 7800.0
NDIVX  = {ndivx}
NDIVY  = {ndivy}
NMODES = 35
PI     = 3.14159265358979

/PREP7
ET, 1, PLANE183
KEYOPT, 1, 3, 3
R, 1, H
MP, EX,1,E_iso
MP, PRXY,1,NU_iso
MP, DENS,1,RHO
RECTNG, -L, L, 0, B
TYPE,1 $ MAT,1 $ REAL,1
LSEL, S, LOC, Y, 0  $ LSEL, A, LOC, Y, B
LESIZE, ALL, , , NDIVX
LSEL, S, LOC, X, -L $ LSEL, A, LOC, X, L
LESIZE, ALL, , , NDIVY
ALLSEL, ALL
MSHKEY, 1 $ AMESH, ALL
! mirror cut at Y=0: {note}
{bc}
ALLSEL, ALL
FINISH

/SOLU
ANTYPE, MODAL
MODOPT, LANB, NMODES, -1.0, , , ON
MXPAND, NMODES, , , YES
SOLVE
FINISH

/POST1
*DIM, MIDX{arr}, ARRAY, NMODES
*VFILL, MIDX{arr}, RAMP, 1, 1
*DIM, FR{arr}, ARRAY, NMODES
*DO, I, 1, NMODES
  *GET, FR{arr}(I), MODE, I, FREQ
*ENDDO
*CFOPEN, {cf}, txt
*VWRITE
('mode      f[Hz]           (PLANE183 {par} mesh{mesh} l/b={lob:.1f} x{ndivx}/y{ndivy}; {conv})')
*VWRITE, MIDX{arr}(1), FR{arr}(1)
(F4.0, 3X, E14.6)
*CFCLOS
FINISH

"""


def header(tag, lob, full_len, ndiv1, ndiv2):
    return f"""! ===========================================================================
!  ansys_rect_ff_ip_halfmodel_{tag}.inp  --  Mechanical APDL batch deck
!  Rectangular free-free (FFFF) IN-PLANE FE ground truth, l/b={lob:.1f}
!  half-model, PLANE183 plane-stress.
!
!  GENERATED 2026-09-05. Geometry/material match the OOP half-model family
!  (B=1 m, H=0.04 m, E=210 GPa, nu=0.30, rho=7800). Element is PLANE183
!  KEYOPT(3)=3 (plane stress), thickness via R,1,H -- same as the validated
!  annular FFFF IP deck ansys_ff_ip_ext.inp (Paper 1 §6.3).
!
!  *VWRITE DISCIPLINE: fill the result array inside the *DO, then *VWRITE
!  once AFTER the loop over the full array from index 1 (ansys_ff_ip_ext.inp).
!  Never *VWRITE an indexed array inside a *DO. Hz only -- do not write an
!  Omega column from APDL (runtime SQRT is untrusted; OOP decks hit this).
!
!  Omega_bar CONVERSION (Seok/Tiersten/Scarton 2004 rect Part 2 Eqs.22-23):
!    omega_bar = (pi/(2B)) * sqrt(c66/rho),  c66 = E/(2(1+nu))
!    Omega = omega/omega_bar = 4*B*f[Hz] / sqrt(c66/rho)
!  For this material and B=1 m: OFAC = 804.481 Hz per unit Omega.
!  Convert the first mode-1 *VWRITE Hz block: Omega = f / 804.481.
!  Independent of LOB.
!
!  HALF-MODEL mirror BC (matches RectIPAssembler sph):
!    SYMM: UY=0 on Y=0  -- solver SYM (u1 even, u2 odd)
!    ASYM: UX=0 on Y=0  -- solver ANTI (u1 odd, u2 even)
!  Mesh: NDIVX = round(2*LOB*20) = {ndiv1} (mesh1), round(2*LOB*30) = {ndiv2}
!  (mesh2); NDIVY = 20 / 30. Same element size as the OOP decks.
!
!  READING:
!    0. Convergence: mesh1 vs mesh2 Hz, same parity, <0.3%.
!    1. Discard near-zero rigid (half-model keeps one in-plane rigid
!       translation per parity). Convert survivors with OFAC=804.481.
!    2. Pair SYMM vs solver SYM, ASYM vs ANTI only.
!    Formula-agnostic FE. No SOLVER_VERSION action from this deck.
!
!  MUST run from Ansys/NewAnsys/ via submit_ansys_queue_rect_ff_ip.sh
!  (one ANSYS license seat). Do NOT submit concurrent with any other MAPDL.
!    mapdl -b -smp -np 8 -i ansys_rect_ff_ip_halfmodel_{tag}.inp -o ansys_rect_ff_ip_halfmodel_{tag}_out.txt
!  Outputs: rect_ff_ip_{{symm,asym}}_mesh{{1,2}}_{tag}.txt
! ===========================================================================

"""


def footer(lob):
    return (
        "! ---------------------------------------------------------------------------\n"
        "!  AFTER RUNNING:\n"
        "!   0. Convergence gate per pass (<0.3%, mesh1 vs mesh2).\n"
        "!   1. Convert the first mode-1 Hz block: Omega = f[Hz] / 804.481\n"
        f"!      (l/b={lob:.1f}; OFAC independent of LOB).\n"
        "!   2. SYMM -> SYM, ASYM -> ANTI. Discard rigid (~0 Hz).\n"
        "!   3. If every candidate misses by a common factor, suspect OFAC\n"
        "!      / Hz extraction before the physics.\n"
        "! ===========================================================================\n"
    )


SPECS = [
    dict(tag="lob100", lob=1.0, full_len=2.0),
    dict(tag="lob150", lob=1.5, full_len=3.0),
    dict(tag="lob200", lob=2.0, full_len=4.0),
    dict(tag="lob250", lob=2.5, full_len=5.0),
    dict(tag="lob300", lob=3.0, full_len=6.0),
]


def deck_body(tag, lob):
    n1, n2 = ndiv(lob)
    parts = [
        one_pass(tag, lob, f"rect IP FFFF l/b={lob:.1f} SYMM mesh1", n1, 20, "SYMM", 1, "S1"),
        one_pass(tag, lob, f"rect IP FFFF l/b={lob:.1f} SYMM mesh2", n2, 30, "SYMM", 2, "S2"),
        one_pass(tag, lob, f"rect IP FFFF l/b={lob:.1f} ASYM mesh1", n1, 20, "ASYM", 1, "A1"),
        one_pass(tag, lob, f"rect IP FFFF l/b={lob:.1f} ASYM mesh2", n2, 30, "ASYM", 2, "A2"),
    ]
    return "".join(parts)


def main():
    names = []
    for s in SPECS:
        n1, n2 = ndiv(s["lob"])
        text = (
            header(s["tag"], s["lob"], s["full_len"], n1, n2)
            + deck_body(s["tag"], s["lob"])
            + footer(s["lob"])
        )
        assert "\r" not in text
        path = HERE / f"ansys_rect_ff_ip_halfmodel_{s['tag']}.inp"
        path.write_bytes(text.encode("utf-8"))
        names.append(path.name)
        assert f"LOB    = {s['lob']:.1f}" in text
        assert text.count("*CFOPEN,") == 4
        assert "*VWRITE, MIDX" in text
        assert "PLANE183" in text
        assert "SHELL281" not in text
        print(f"wrote {path.name} bytes={path.stat().st_size} lf={text.count(chr(10))}")
    return names


if __name__ == "__main__":
    main()
