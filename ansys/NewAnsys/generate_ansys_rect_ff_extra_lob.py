# Generate l/b=2.0 and l/b=3.0 half-model decks from lob100.
# APDL body is a scripted copy (LOB / NDIVX / *CFOPEN / /TITLE only) so the
# known *VWRITE-inside-*DO trap is not retyped. Header/footer are rewritten
# for the extra-lob FE job (2453876 PAPER_CANDIDATES, KFAC=24.6644).
from __future__ import annotations

from pathlib import Path

HERE = Path(__file__).resolve().parent
SRC = HERE / "ansys_rect_ff_oop_halfmodel_lob100.inp"
BODY_START = "! ======================= PASS SYMM, mesh1 =======================\n"
BODY_END = "! ---------------------------------------------------------------------------\n"


def header(tag, lob, full_len, ndiv1, ndiv2):
    return f"""! ===========================================================================
!  ansys_rect_ff_oop_halfmodel_{tag}.inp  --  Mechanical APDL batch deck
!  Rectangular free-free (FFFF) OOP FE ground truth, l/b={lob:.1f} half-model.
!
!  GENERATED 2026-09-04 by copying ansys_rect_ff_oop_halfmodel_lob100.inp
!  (that file itself copied 2026-09-02 from the validated l/b=1.5 deck
!  ansys_rect_ff_oop_halfmodel.inp, job 2398617) and changing ONLY: LOB,
!  NDIVX, *CFOPEN basenames, and /TITLE. Mesh scaling keeps the SAME
!  element size as the l/b=1.5 deck -- elem size 1/20 of B for mesh1,
!  1/30 of B for mesh2, i.e. NDIVX = round(2*LOB*20) = {ndiv1} (mesh1)
!  and round(2*LOB*30) = {ndiv2} (mesh2). Tagged _{tag} so outputs do
!  not collide with lob100 / lob250 / l/b=1.5 when run from NewAnsys/.
!  The *VWRITE/*DO structure, SYMM/ASYM mirror BCs, material constants,
!  and KFAC formula are copied, not retyped.
!
!  PURPOSE: job 2453876 produced Screen-B PAPER_CANDIDATES at this l/b
!  (listed in the footer) but has no FE at l/b=2.0 or 3.0. This deck is
!  formula-agnostic FE ground truth. Convert Hz with KFAC=24.6644 from
!  the first mode-1 *VWRITE block; do NOT use the deck Lambda_FE column
!  (runtime SQRT missing). Pair SYMM vs SYM and ASYM vs ANTI only.
!  No SOLVER_VERSION action from this deck alone.
!
!  GEOMETRY/MATERIAL (same convention as every other deck in this family):
!    half-breadth B=1.0 m, aspect ratio l/b={lob:.1f} (full length {full_len:.1f} m,
!    full breadth 2.0 m), full thickness H=0.04 m (H/(2B)=0.02, thin-plate
!    valid), E=210 GPa, nu=0.30, rho=7800 kg/m^3.
!
!  Lambda -> Hz CONVERSION: identical KFAC derivation to the l/b=1.5 deck
!  (Seok/Tiersten/Scarton 2004 rect Part 1 Eqs.24-28,37-38,41-42) --
!  KFAC = (pi*(H/2))/(4*B^2) * sqrt(G/(6*rho*(1-nu))), independent of LOB.
!  At runtime the APDL SQRT is missing, so ignore Lambda_FE and use
!  Lambda = f[Hz] / 24.6644 from the first mode-1 Hz block.
!
!  WHY A HALF-MODEL / mirror BC pattern: unchanged from the l/b=1.5 deck --
!  SYMM (UY=0,ROTX=0,ROTZ=0) matches solver phin=pi/2 (SYM, even);
!  ASYM (UX=0,UZ=0,ROTY=0) matches solver phin=0 (ANTI, odd).
!
!  READING THIS DECK'S OUTPUT:
!    0. CONVERGENCE GATE FIRST: mesh1 vs mesh2 Hz lists (same parity)
!       must agree <0.3% before reading anything else.
!    1. Discard near-zero (rigid body) entries and anything with
!       f>30 Hz (membrane-mode window).
!    2. Convert surviving Hz with KFAC=24.6644. Compare SYMM-pass against
!       SYM candidates and ASYM-pass against ANTI candidates ONLY.
!    EXPLORATORY. No SOLVER_VERSION action, no table promotion from this
!    deck alone regardless of outcome.
!
!  MUST run from Ansys/NewAnsys/ (one ANSYS license seat, isolated queue
!  submit_ansys_queue_rect_ff_extra_lob.sh). Do NOT submit concurrently
!  with any other MAPDL job.
!    mapdl -b -smp -np 8 -i ansys_rect_ff_oop_halfmodel_{tag}.inp -o ansys_rect_ff_oop_halfmodel_{tag}_out.txt
!  Outputs: rect_ff_symm_mesh1_{tag}.txt, rect_ff_symm_mesh2_{tag}.txt,
!           rect_ff_asym_mesh1_{tag}.txt, rect_ff_asym_mesh2_{tag}.txt
! ===========================================================================

"""


def footer(lob, sym_cands, anti_cands, extra_notes):
    sym = ", ".join(f"{x:.3f}" for x in sym_cands)
    anti = ", ".join(f"{x:.3f}" for x in anti_cands)
    notes = "".join("!      " + n + "\n" for n in extra_notes)
    return (
        "! ---------------------------------------------------------------------------\n"
        "!  AFTER RUNNING (order matters):\n"
        "!   0. Convergence gate per pass (<0.3%, mesh1 vs mesh2).\n"
        "!   1. Do NOT use the deck Lambda_FE column (runtime SQRT missing).\n"
        "!      Convert the first mode-1 *VWRITE Hz block with KFAC=24.6644:\n"
        "!      Lambda = f[Hz] / 24.6644. Discard near-zero (rigid-body) Hz and\n"
        "!      anything implying f>30 Hz.\n"
        f"!   2. Compare SYMM-pass Lambdas to l/b={lob:.1f} SYM B-pass candidates\n"
        f"!      ({sym}); ASYM-pass to ANTI candidates ({anti}).\n"
        "!      Cross-parity matches do not count. Candidates are job 2453876\n"
        "!      PAPER_CANDIDATES (Screen B persist, n_cpair=3, |dL|<=0.01) -- NOT\n"
        "!      yet FE-validated. Exclude the known SYM 0.41 leak (and at l/b=3.0\n"
        "!      the ANTI 0.440 member of the 0.46 leak family).\n"
        + notes
        + "!   3. If every candidate misses by a common factor, suspect KFAC / Hz\n"
        "!      extraction before the physics.\n"
        "! ===========================================================================\n"
    )


def transform_body(body, *, tag, lob, ndiv1, ndiv2):
    t = body
    t = t.replace("l/b=1.0", f"l/b={lob:.1f}")
    t = t.replace("LOB    = 1.0", f"LOB    = {lob:.1f}")
    t = t.replace("NDIVX  = 40", f"NDIVX  = {ndiv1}")
    t = t.replace("NDIVX  = 60", f"NDIVX  = {ndiv2}")
    t = t.replace("_lob100", f"_{tag}")
    return t


SPECS = [
    dict(
        tag="lob200",
        lob=2.0,
        ndiv1=80,
        ndiv2=120,
        full_len=4.0,
        footer=footer(
            2.0,
            [1.510, 1.530, 1.930, 2.230, 2.270],
            [0.680, 1.480, 1.640, 2.080, 2.120],
            [
                "SYM 1.930 is shallow (sig~4e-2); ANTI 1.640 is the recurring",
                "unmatched extra seen at l/b=1.0/1.5/2.5 -- still compare.",
            ],
        ),
    ),
    dict(
        tag="lob300",
        lob=3.0,
        ndiv1=120,
        ndiv2=180,
        full_len=6.0,
        footer=footer(
            3.0,
            [0.240, 1.320, 1.920, 2.160, 2.260],
            [0.540, 0.940, 1.160, 1.520, 1.640, 1.920, 2.120, 2.260],
            [
                "Do NOT treat ANTI 0.440 as a paper candidate (0.46-family leak",
                "that Screen B passed only because 0.440 is just outside +/-0.015",
                "of 0.460). SYM 0.41 is the known leak. ANTI 1.640 is the recurring extra.",
            ],
        ),
    ),
]


def main():
    src = SRC.read_text(encoding="utf-8", newline="")
    if "\r" in src:
        raise SystemExit("source deck is not LF")
    i0 = src.find(BODY_START)
    i1 = src.find(BODY_END)
    if i0 < 0 or i1 < 0 or i1 <= i0:
        raise SystemExit("could not split lob100 header/body/footer")
    body = src[i0:i1]
    for s in SPECS:
        out_body = transform_body(
            body, tag=s["tag"], lob=s["lob"], ndiv1=s["ndiv1"], ndiv2=s["ndiv2"]
        )
        out = (
            header(s["tag"], s["lob"], s["full_len"], s["ndiv1"], s["ndiv2"])
            + out_body
            + s["footer"]
        )
        path = HERE / f"ansys_rect_ff_oop_halfmodel_{s['tag']}.inp"
        path.write_bytes(out.encode("utf-8"))
        assert f"LOB    = {s['lob']:.1f}" in out
        assert "LOB    = 1.0" not in out
        assert f"NDIVX  = {s['ndiv1']}" in out
        assert f"NDIVX  = {s['ndiv2']}" in out
        assert "NDIVX  = 40" not in out
        assert "_lob100" not in out_body
        assert "copying ansys_rect_ff_oop_halfmodel_lob100.inp" in out
        assert "B      = 1.0" in out
        assert "NDIVY  = 20" in out and "NDIVY  = 30" in out
        assert out.count(f"LOB    = {s['lob']:.1f}") == 4
        assert out.count(f"NDIVX  = {s['ndiv1']}") == 2
        assert out.count(f"NDIVX  = {s['ndiv2']}") == 2
        assert out.count(f"*CFOPEN, rect_ff_symm_mesh1_{s['tag']}") == 1
        assert out.count(f"*CFOPEN, rect_ff_symm_mesh2_{s['tag']}") == 1
        assert out.count(f"*CFOPEN, rect_ff_asym_mesh1_{s['tag']}") == 1
        assert out.count(f"*CFOPEN, rect_ff_asym_mesh2_{s['tag']}") == 1
        assert f"/TITLE, rect half-model FFFF l/b={s['lob']:.1f}" in out
        assert "24.6644" in out
        assert "\r" not in out
        print(f"wrote {path.name} bytes={path.stat().st_size} lf={out.count(chr(10))}")


if __name__ == "__main__":
    main()
