# -*- coding: utf-8 -*-
"""
probe_ip_ffp1_spurious_table_v1.py -- DIAGNOSTIC/DATA-GENERATION ONLY.
No SOLVER_VERSION bump.

Regenerates residual-screen numbers for the 8 FF-P1 in-plane MAC
calibration candidates (4 REAL from Table 7, 4 catalogued artifacts of
the same determinant) so an IP analogue of Table 6 can be printed.

n_dofs=20, same production residual-adjudication size as §5. Geometry
FF-P1: r0/2b=1.5, 2Theta/pi=0.5, nu=0.30.

Reconstruction imported from probe_ip_mac_adversarial_fullblock_v1.py.
No ANSYS dump required.

-----------------------------------------------------------------------------
PRE-REGISTERED INTERPRETATION
-----------------------------------------------------------------------------
ST-0 GATE_FAIL
    SOLVER_VERSION mismatch or a reconstruction UNRELIABLE. VOID.
ST-A TABLE_READY
    All 4 REAL-like have score=max(rTyy,rTyr) < 0.3, all 4 ARTIFACT-like
    have score > 4, and no label flips vs the published REAL/ARTIFACT
    assignment. -> print these 8 rows as the IP Table-6 analogue.
ST-B BAND_MOVED
    Labels hold but a residual has left its published band (<0.3 or >4).
    Report the new numbers; do not paste into the paper without looking.
ST-C LABEL_FLIP
    Any REAL came back ARTIFACT-like (score>4) or vice versa (score<0.3).
    Do not put a new table in the paper. Re-examine the seed list.

Honour the printed letter. Residual ratios in Table 7 stay as they are;
this only adds the missing artifact table.
"""
from __future__ import annotations

import importlib.util
import os
import sys
import time

sys.path.insert(0, os.environ.get("PKG_PATH", "."))
os.environ.setdefault("DPS", "40")

EXPECT_VER = os.environ.get("EXPECT_SOLVER_VERSION", "2026-07-10.s10")
PKG = os.environ.get("PKG_PATH", ".")
THETA = float(__import__("numpy").pi) / 4.0

# Same 8 as probe_ip_block_selection_rule_v1.py GEOM=ffp1
CANDS = [
    ("real_a", 0.401154, "REAL", 1),
    ("real_b", 1.142287, "REAL", 4),
    ("real_c", 1.598985, "REAL", 10),
    ("real_d", 2.521489, "REAL", 23),
    ("art_a",  0.542105, "ARTIFACT", None),
    ("art_b",  1.021144, "ARTIFACT", None),
    ("art_c",  1.454616, "ARTIFACT", None),
    ("art_d",  1.603259, "ARTIFACT", None),
]
REAL_HI = 0.3
ART_LO = 4.0


def hdr(s):
    print("\n" + "=" * 78 + f"\n  {s}\n" + "=" * 78, flush=True)


def find_file(*names):
    roots = [PKG, os.path.join(PKG, "Ansys", "NewAnsys"), ".", 
             os.path.join("Ansys", "NewAnsys")]
    for n in names:
        for d in roots:
            p = os.path.join(d, n)
            if os.path.isfile(p):
                return os.path.abspath(p)
    return None


def load_adv():
    path = find_file("probe_ip_mac_adversarial_fullblock_v1.py")
    if path is None:
        print("READING: ST-0 GATE_FAIL  (adv probe missing)", flush=True)
        raise SystemExit(2)
    spec = importlib.util.spec_from_file_location("adv_st", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    mod.THETA = float(THETA)
    return mod


def main():
    t_all = time.time()
    hdr("probe_ip_ffp1_spurious_table_v1")
    import ip_mac_bar as bar
    bar.ensure_free_free_worker_bc()
    import plate_solver as ps
    print(f"  SOLVER_VERSION={ps.SOLVER_VERSION!r}  expect {EXPECT_VER!r}",
          flush=True)
    if ps.SOLVER_VERSION != EXPECT_VER:
        print("READING: ST-0 GATE_FAIL  (SOLVER_VERSION)", flush=True)
        raise SystemExit(3)

    adv = load_adv()
    mat = ps.IsotropicMaterial(E=210e9, nu=0.30, rho=7800.0)
    geom = ps.make_geometry(1.5, 0.5)
    solver = ps.InPlaneSolver(geom, mat, M=80, n_quad=30,
                              boundary=ps.FreeFreeIP())

    rows = []
    for name, Om, lab, tab7 in CANDS:
        hdr(f"{name}  Omega={Om}  published {lab}")
        n_dofs = adv.N_DOFS
        brs, _x, st, cnt, _n = adv.get_filled_brs(solver, Om, n_dofs)
        print(f"  fill state={st} cnt={cnt}/{n_dofs}", flush=True)
        if st == "UNRELIABLE":
            print("READING: ST-0 GATE_FAIL  (UNRELIABLE)", flush=True)
            raise SystemExit(4)
        sh = adv.reconstruct_shared(solver, Om, brs, n_dofs=n_dofs)
        if not sh.get("ok"):
            print("READING: ST-0 GATE_FAIL  (reconstruct)", flush=True)
            raise SystemExit(4)
        blocks = {}
        for w in ("q0", "q1"):
            r = adv.block_reconstruction(solver, sh, w)
            if r.get("ok"):
                blocks[w] = r
                print(f"  {w}: rTyy={r['rTyy']:.6g}  rTyr={r['rTyr']:.6g}",
                      flush=True)
        s0, s1, _c = adv.sigma_pair_at(solver, Om)
        print(f"  log10 sig  q0={s0}  q1={s1}", flush=True)
        # residual screen as published: max over the reconstructed blocks
        scores = []
        for w, r in blocks.items():
            scores.append((max(r["rTyy"], r["rTyr"]), w, r))
        scores.sort()
        best = scores[0]  # smaller residual owns the "physical-looking" read
        # Published IP split uses the residual on the evaluation that the
        # paper quotes. Quote BOTH blocks; classify on the SMALLER score
        # (a REAL should have at least one block in-band).
        score_lo, w_lo, r_lo = best
        score_hi = max(max(r["rTyy"], r["rTyr"]) for r in blocks.values())
        if score_lo < REAL_HI:
            got = "REAL-like"
        elif score_lo > ART_LO:
            got = "ARTIFACT-like"
        else:
            got = "MIDDLE"
        print(f"  min score={score_lo:.6g} on {w_lo}  max-block={score_hi:.6g}  "
              f"got={got}  expect={lab}", flush=True)
        rows.append(dict(
            name=name, Om=Om, lab=lab, tab7=tab7, pick=w_lo,
            rTyy=r_lo["rTyy"], rTyr=r_lo["rTyr"], score=score_lo,
            got=got, sig0=s0, sig1=s1,
        ))

    hdr("TABLE (copy into SM if ST-A)")
    print(f"  {'name':8s} {'Om':10s} {'blk':3s} {'rTyy':10s} {'rTyr':10s} "
          f"{'score':10s} got          pub  T7", flush=True)
    for r in rows:
        t7 = f"#{r['tab7']}" if r["tab7"] else "—"
        print(f"  {r['name']:8s} {r['Om']:10.6f} {r['pick']:3s} "
              f"{r['rTyy']:10.5g} {r['rTyr']:10.5g} {r['score']:10.5g} "
              f"{r['got']:12s} {r['lab']:9s} {t7}", flush=True)

    flips = [r for r in rows
             if (r["lab"] == "REAL" and r["got"] == "ARTIFACT-like")
             or (r["lab"] == "ARTIFACT" and r["got"] == "REAL-like")]
    band = []
    for r in rows:
        if r["lab"] == "REAL" and not (r["score"] < REAL_HI):
            band.append(r["name"])
        if r["lab"] == "ARTIFACT" and not (r["score"] > ART_LO):
            band.append(r["name"])

    hdr("READING (pre-registered -- honour the letter)")
    if flips:
        letter = "ST-C LABEL_FLIP"
        note = f"flips={ [r['name'] for r in flips] }. Do not add a table."
    elif band:
        letter = "ST-B BAND_MOVED"
        note = f"out of band: {band}. Report; do not paste blindly."
    else:
        letter = "ST-A TABLE_READY"
        note = ("4 REAL score<0.3, 4 ART score>4. These 8 rows are the IP "
                "analogue of Table 6.")
    print(f"\nREADING: {letter}", flush=True)
    print(f"  {note}", flush=True)
    print(f"elapsed_s={time.time()-t_all:.1f}", flush=True)


if __name__ == "__main__":
    main()
