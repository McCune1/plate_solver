# -*- coding: utf-8 -*-
"""
probe_ip_ffp1_spurious_table_v2.py -- DIAGNOSTIC/DATA-GENERATION ONLY.
No SOLVER_VERSION bump.

v2 of job 2434288 (ST-B BAND_MOVED). That letter fired because v1 scored
the SMALLER residual of the two parity blocks. Artifacts' quiet block is
not the root-owning one; on the deeper block the four artifacts were
6.91 / 8.29 / 7.03 / 5.20 and matched GATE bit-for-bit. real_d's 0.686
on q1 is also the published GATE number -- it was never inside the ~0.3
bar.

v2 scores RULE_LOCALMIN_RESID's owning block, gates residuals against
GATE (job 2431872), and classifies ARTIFACT as owning-block
score=max(rTyy,rTyr) > 4. REAL is "not ARTIFACT" (score is not >4);
real_d at ~0.686 is expected and is not a band failure.

No ANSYS dump required.

-----------------------------------------------------------------------------
PRE-REGISTERED INTERPRETATION
-----------------------------------------------------------------------------
ST-0 GATE_FAIL
    SOLVER_VERSION, UNRELIABLE fill, or any of the 8 owning-block
    (rTyy, rTyr) disagrees with GATE[name, pick] at rel 1e-4. VOID.
ST-A TABLE_READY
    Gates pass, 4 ART owning-score > 4, 4 REAL owning-score is not > 4,
    and in-run pick matches the deeper-sigma block. -> these 8 rows are
    the IP analogue of Table 6 (owning-block residuals).
ST-C LABEL_FLIP
    An owning-block REAL scores > 4, or an ART scores <= 4. Do not add
    a table.

There is no ST-B in v2: the ~0.3 REAL ceiling is not a close condition
(real_d is the known exception, published in GATE).
"""
from __future__ import annotations

import importlib.util
import os
import sys
import time

sys.path.insert(0, os.environ.get("PKG_PATH", "."))
os.environ.setdefault("DPS", "40")

import numpy as np

EXPECT_VER = os.environ.get("EXPECT_SOLVER_VERSION", "2026-07-10.s10")
PKG = os.environ.get("PKG_PATH", ".")
THETA = float(np.pi) / 4.0
RULE_DELTA = float(os.environ.get("RULE_DELTA", "0.001"))
RES_RTOL = 1e-4
ART_LO = 4.0

CANDS = [
    ("real_a", 0.401154, "REAL", 1, "t1"),
    ("real_b", 1.142287, "REAL", 4, "t2"),
    ("real_c", 1.598985, "REAL", 10, "t3"),
    ("real_d", 2.521489, "REAL", 23, "t4"),
    ("art_a",  0.542105, "ARTIFACT", None, "t5"),
    ("art_b",  1.021144, "ARTIFACT", None, "t6"),
    ("art_c",  1.454616, "ARTIFACT", None, "t7"),
    ("art_d",  1.603259, "ARTIFACT", None, "t8"),
]
GATE = {
    ("t1", "q0"): (3.48852, 1.57244),
    ("t1", "q1"): (0.0305187, 0.0180822),
    ("t2", "q0"): (0.0498155, 0.0353367),
    ("t2", "q1"): (0.740035, 0.35957),
    ("t3", "q0"): (1.9065, 1.00399),
    ("t3", "q1"): (0.0485492, 0.0216588),
    ("t4", "q0"): (2.66673, 2.13439),
    ("t4", "q1"): (0.685768, 0.320466),
    ("t5", "q0"): (6.90722, 4.80055),
    ("t5", "q1"): (3.28749, 1.68239),
    ("t6", "q0"): (3.17953, 1.695),
    ("t6", "q1"): (6.05769, 8.28675),
    ("t7", "q0"): (7.03462, 3.28804),
    ("t7", "q1"): (1.94664, 1.06654),
    ("t8", "q0"): (5.19862, 0.916939),
    ("t8", "q1"): (0.735771, 0.434507),
}


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
    spec = importlib.util.spec_from_file_location("adv_st2", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    mod.THETA = float(THETA)
    return mod


def rel_err(got, pub):
    pub = float(pub)
    if abs(pub) < 1e-15:
        return abs(float(got) - pub)
    return abs(float(got) - pub) / abs(pub)


def main():
    t_all = time.time()
    hdr("probe_ip_ffp1_spurious_table_v2  (owning-block residuals)")
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
    gate_fail = []
    for name, Om, lab, tab7, tag in CANDS:
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
            if not r.get("ok"):
                print(f"  {w}: FAIL {r.get('reason')}", flush=True)
                continue
            blocks[w] = r
            print(f"  {w}: rTyy={r['rTyy']:.6g}  rTyr={r['rTyr']:.6g}",
                  flush=True)
        s0, s1, _c = adv.sigma_pair_at(solver, Om)
        m0, m1, _c0 = adv.sigma_pair_at(solver, Om - RULE_DELTA)
        p0, p1, _c1 = adv.sigma_pair_at(solver, Om + RULE_DELTA)
        print(f"  log10 sig  q0={s0}  q1={s1}", flush=True)
        pick, how = bar.rule_localmin_resid(
            s0, s1, m0, p0, m1, p1, blocks["q0"]["rTyy"], blocks["q1"]["rTyy"],
        )
        print(f"  RULE_LOCALMIN_RESID -> {pick}  [{how}]", flush=True)
        r_own = blocks[pick]
        score = max(r_own["rTyy"], r_own["rTyr"])
        pub = GATE[(tag, pick)]
        e0 = rel_err(r_own["rTyy"], pub[0])
        e1 = rel_err(r_own["rTyr"], pub[1])
        print(f"  owning {pick}: score={score:.6g}  GATE rel "
              f"{e0:.2e}/{e1:.2e}", flush=True)
        if e0 > RES_RTOL or e1 > RES_RTOL:
            gate_fail.append(f"{name} residual vs GATE")
        if score > ART_LO:
            got = "ARTIFACT-like"
        else:
            got = "REAL-like"
        print(f"  got={got}  expect={lab}", flush=True)
        rows.append(dict(
            name=name, Om=Om, lab=lab, tab7=tab7, tag=tag, pick=pick,
            how=how, rTyy=r_own["rTyy"], rTyr=r_own["rTyr"], score=score,
            got=got, sig0=s0, sig1=s1,
        ))

    if gate_fail:
        print(f"\nREADING: ST-0 GATE_FAIL  {gate_fail}", flush=True)
        raise SystemExit(4)

    hdr("TABLE (copy into SM if ST-A; owning-block residuals)")
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

    hdr("READING (pre-registered -- honour the letter)")
    if flips:
        letter = "ST-C LABEL_FLIP"
        note = f"flips={[r['name'] for r in flips]}. Do not add a table."
    else:
        letter = "ST-A TABLE_READY"
        note = ("Owning-block residuals match GATE; 4 ART score>4; 4 REAL "
                "are not ARTIFACT. These 8 rows are the IP analogue of "
                "Table 6. real_d score~0.686 is expected (GATE t4 q1).")
    print(f"\nREADING: {letter}", flush=True)
    print(f"  {note}", flush=True)
    print(f"elapsed_s={time.time()-t_all:.1f}", flush=True)


if __name__ == "__main__":
    main()
