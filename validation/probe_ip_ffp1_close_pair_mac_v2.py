# -*- coding: utf-8 -*-
"""
probe_ip_ffp1_close_pair_mac_v2.py -- DIAGNOSTIC ONLY.
No SOLVER_VERSION bump. No package detector change.

v2 of job 2434286/2434287. Those printed CP-0 because GATE_ART's
full-block max MAC (0.6281 at FE#9) was gated against the published
8-target MAC 0.0067 -- a different statistic (same selection effect as
the r200/a100 adversarial full-block test). Reconstruction was healthy:
residuals bit-matched GATE, pick q0, CALL ARTIFACT-like at 0.744.
The isolates were never scored.

v2 gates:
  GATE_REAL  Omega=0.401154  pick q1, residuals match GATE, max-same
             vs published 1.0000 to 5e-4, CALL REAL-like.
  GATE_ART   Omega=1.454616  pick q0, residuals match GATE (rel 1e-4),
             CALL ARTIFACT-like at 0.744. Full-block max-same is
             PRINTED, not gated (job 2434286: 0.6281).
Then scores the isolates. Same dump ip_mac_ffp1_nu030. No new ANSYS.

PAIR=1115  Table 7 1114.9/1115.3 Hz (job 2428351 isolates)
PAIR=qin   Qin 4.4293/4.4754 (1.14274 and 1.14758)

-----------------------------------------------------------------------------
PRE-REGISTERED INTERPRETATION (honour the printed letter)
-----------------------------------------------------------------------------
CP-0 GATE_FAIL
    SOLVER_VERSION, missing dumps, mac() identity, THETA, residual
    mismatch vs GATE, wrong pick, GATE_REAL not REAL-like, or GATE_ART
    not ARTIFACT-like at 0.744. VOID. Do not edit the paper.
CP-A TWO_SHAPES
    Both isolates MAC-REAL, distinct best same-class FE modes.
CP-B ONE_SHAPE
    Both isolates MAC-REAL against the SAME best same-class FE mode.
CP-C SPLIT
    One MAC-REAL, one MAC-ART. Detector still lists one root.
CP-D BOTH_ART
    Neither isolate is MAC-REAL.
CP-E UNRELIABLE
    Gates pass but an isolate failed to reconstruct.

Do not move 0.744. Do not rewrite Table 7 or Table 9 from this run.
The footnote currently says "MAC was not run" -- that is what this
replaces, if the letter is A/B/C/D.
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
PAIR = os.environ.get("PAIR", "1115").strip().lower()
PREFIX = os.environ.get("FB_PREFIX", "ip_mac_ffp1_nu030")
RULE_DELTA = float(os.environ.get("RULE_DELTA", "0.001"))
MAC_ATOL = 5e-4
RES_RTOL = 1e-4
PHI_DEG = 90.0
THETA = float(np.pi) / 4.0
WBAR = 2527.351

PAIRS = {
    "1115": {
        "label": "Table 7 1114.9/1115.3 Hz",
        "cands": [
            ("PROD", 2.771955, "production Table 7 #26"),
            ("ISO_A", 2.772422, "isolate -> 1115.07 Hz"),
            ("ISO_B", 2.773399, "isolate -> 1115.46 Hz"),
        ],
        "fe_hz": (1114.9, 1115.3),
    },
    "qin": {
        "label": "Qin 4.4293/4.4754",
        "cands": [
            ("PROD", 1.14274, "production n=20/28, Qin 4.4529"),
            ("ISO_B", 1.14758, "isolated dip, Qin 4.4717"),
        ],
        "fe_hz": (459.5, 460.6),
    },
}

# job 2431872 / block-selection GATE residuals (both blocks)
GATE_RES = {
    ("GATE_REAL", "q0"): (3.48852, 1.57244),
    ("GATE_REAL", "q1"): (0.0305187, 0.0180822),
    ("GATE_ART", "q0"): (7.03462, 3.28804),
    ("GATE_ART", "q1"): (1.94664, 1.06654),
}
GATES = [
    ("GATE_REAL", 0.401154, "q1", "REAL-like", 1.0000),
    ("GATE_ART", 1.454616, "q0", "ARTIFACT-like", None),
]


def hdr(s):
    print("\n" + "=" * 78 + f"\n  {s}\n" + "=" * 78, flush=True)


def find_file(*names):
    roots = [PKG, os.path.join(PKG, "Ansys", "NewAnsys"),
             ".", os.path.join("Ansys", "NewAnsys")]
    for n in names:
        for d in roots:
            p = os.path.join(d, n)
            if os.path.isfile(p):
                return os.path.abspath(p)
    return None


def load_adv():
    path = find_file("probe_ip_mac_adversarial_fullblock_v1.py")
    if path is None:
        print("READING: CP-0 GATE_FAIL  (adv probe missing)", flush=True)
        raise SystemExit(2)
    spec = importlib.util.spec_from_file_location("adv_ffp1_cp2", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    mod.THETA = float(THETA)
    return mod


def f_hz(om):
    return float(om) * WBAR / (2.0 * np.pi)


def rel_err(got, pub):
    pub = float(pub)
    if abs(pub) < 1e-15:
        return abs(float(got) - pub)
    return abs(float(got) - pub) / abs(pub)


def score_one(name, Om, adv, bar, solver, FEB, mode_cls, modes):
    t0 = time.time()
    n_dofs = adv.N_DOFS
    brs, _x, st, cnt, _n = adv.get_filled_brs(solver, Om, n_dofs)
    print(f"  fill state={st} cnt={cnt}/{n_dofs}  f={f_hz(Om):.4f} Hz",
          flush=True)
    if st == "UNRELIABLE":
        return dict(name=name, Om=Om, ok=False, why="UNRELIABLE")
    sh = adv.reconstruct_shared(solver, Om, brs, n_dofs=n_dofs)
    if not sh.get("ok"):
        return dict(name=name, Om=Om, ok=False, why=str(sh.get("reason")))
    blocks = {}
    for w in ("q0", "q1"):
        r = adv.block_reconstruction(solver, sh, w)
        if not r.get("ok"):
            print(f"  block {w}: FAIL {r.get('reason')}", flush=True)
            continue
        blocks[w] = r
        print(f"  block {w}: rTyy={r['rTyy']:.6g}  rTyr={r['rTyr']:.6g}",
              flush=True)
    s0, s1, _c = adv.sigma_pair_at(solver, Om)
    m0, m1, _c0 = adv.sigma_pair_at(solver, Om - RULE_DELTA)
    p0, p1, _c1 = adv.sigma_pair_at(solver, Om + RULE_DELTA)
    if None in (s0, s1, m0, m1, p0, p1) or "q0" not in blocks or "q1" not in blocks:
        return dict(name=name, Om=Om, ok=False, why="rule inputs")
    pick, how = bar.rule_localmin_resid(
        s0, s1, m0, p0, m1, p1, blocks["q0"]["rTyy"], blocks["q1"]["rTyy"],
    )
    print(f"  RULE_LOCALMIN_RESID -> {pick}  [{how}]", flush=True)
    want_cls = bar.block_class(pick)
    rp, tp, Ur, Uth = adv.build_grid_from_cfull(
        solver, sh, blocks[pick]["c_full"])
    ux_s, uy_s, bad = adv.solver_field_at_fe(rp, tp, Ur, Uth, FEB[modes[0]])
    keep = ~bad
    bvec = np.concatenate([ux_s[keep], uy_s[keep]])
    Ablk = np.vstack([np.concatenate([FEB[m][2][keep], FEB[m][3][keep]])
                      for m in modes])
    macs = adv.mac_matrix(Ablk, bvec)
    mac_by = {m: float(macs[i]) for i, m in enumerate(modes)}
    same = [m for m in modes if mode_cls.get(m) == want_cls]
    if not same:
        return dict(name=name, Om=Om, ok=False, why="no same-class FE")
    mac_best, best = max((mac_by[m], m) for m in same)
    call = bar.classify_mac(mac_best)
    verb = adv.mac_against(rp, tp, Ur, Uth, FEB[best])
    print(f"  class {want_cls}  max-same={mac_best:.4f} at #{best}  "
          f"verbatim={verb:.6f}  CALL {call}", flush=True)
    print(f"  [{time.time()-t0:.0f}s]", flush=True)
    return dict(
        name=name, Om=Om, ok=True, pick=pick, how=how, want_cls=want_cls,
        mac_best=mac_best, best=best, call=call, blocks=blocks,
        rTyy=blocks[pick]["rTyy"], rTyr=blocks[pick]["rTyr"],
        f_hz=f_hz(Om),
    )


def main():
    t_all = time.time()
    if PAIR not in PAIRS:
        print(f"FATAL: PAIR={PAIR!r} not in {sorted(PAIRS)}", flush=True)
        raise SystemExit(2)
    spec = PAIRS[PAIR]
    hdr(f"probe_ip_ffp1_close_pair_mac_v2  PAIR={PAIR}  {spec['label']}")
    print("  v2: GATE_ART is residual+pick+CALL, not 8-target MAC vs "
          "full-block max", flush=True)

    import ip_mac_bar as bar
    bar.ensure_free_free_worker_bc()
    import plate_solver as ps
    print(f"  SOLVER_VERSION={ps.SOLVER_VERSION!r}  expect {EXPECT_VER!r}",
          flush=True)
    print(f"  PREFIX={PREFIX}  R40_BC_KIND={os.environ.get('R40_BC_KIND')!r}",
          flush=True)
    if ps.SOLVER_VERSION != EXPECT_VER:
        print("READING: CP-0 GATE_FAIL  (SOLVER_VERSION)", flush=True)
        raise SystemExit(3)

    adv = load_adv()
    print(f"  adv.THETA={float(adv.THETA):.6f}  expect {THETA:.6f}", flush=True)
    if abs(float(adv.THETA) - THETA) > 1e-12:
        print("READING: CP-0 GATE_FAIL  (THETA)", flush=True)
        raise SystemExit(4)

    rng = np.random.default_rng(21)
    a, b = rng.normal(size=200), rng.normal(size=200)
    dmac = abs(bar.mac(a, b) - adv.mac(a, b))
    print(f"  mac() identity |d|={dmac:.2e}", flush=True)
    if dmac > 1e-15:
        print("READING: CP-0 GATE_FAIL  (mac)", flush=True)
        raise SystemExit(4)

    freq_path = find_file(f"{PREFIX}_freqs.txt")
    if freq_path is None:
        print(f"READING: CP-0 GATE_FAIL  (missing {PREFIX}_freqs.txt)",
              flush=True)
        raise SystemExit(2)
    FREQ = adv.load_freqs(freq_path)
    FEB = {}
    for m in range(4, 61):
        p = find_file(f"{PREFIX}_m{m}_mesh96.txt")
        if p is None:
            print(f"READING: CP-0 GATE_FAIL  (missing mode {m})", flush=True)
            raise SystemExit(2)
        FEB[m] = adv.load_fe_eigvec(p)
    modes = list(range(4, 61))
    print(f"  dumps 57/57  n_nodes={len(FEB[4][0])}  {freq_path}", flush=True)

    mode_cls = {}
    for m in modes:
        r, th, ux, uy = FEB[m]
        cls, sr, ar, miss = bar.classify_fe_mode(
            r, th, ux, uy, phi_deg=PHI_DEG)
        mode_cls[m] = cls
    n_s = sum(1 for c in mode_cls.values() if c == "SYMM")
    n_a = sum(1 for c in mode_cls.values() if c == "ANTI")
    n_u = sum(1 for c in mode_cls.values() if c == "AMBIGUOUS")
    print(f"  live parity SYMM={n_s} ANTI={n_a} AMBIGUOUS={n_u}", flush=True)

    mat = ps.IsotropicMaterial(E=210e9, nu=0.30, rho=7800.0)
    geom = ps.make_geometry(1.5, 0.5)
    solver = ps.InPlaneSolver(geom, mat, M=80, n_quad=30,
                              boundary=ps.FreeFreeIP())

    gate_fail = []
    for gname, Om, blk_pub, call_pub, mac_pub in GATES:
        hdr(f"{gname}  Omega={Om}  expect pick {blk_pub} CALL {call_pub}")
        row = score_one(gname, Om, adv, bar, solver, FEB, mode_cls, modes)
        if not row.get("ok"):
            gate_fail.append(f"{gname} {row.get('why')}")
            continue
        for w in ("q0", "q1"):
            pub = GATE_RES[(gname, w)]
            got = (row["blocks"][w]["rTyy"], row["blocks"][w]["rTyr"])
            e0, e1 = rel_err(got[0], pub[0]), rel_err(got[1], pub[1])
            print(f"  residual {w}: got {got[0]:.6g}/{got[1]:.6g}  "
                  f"pub {pub[0]:.6g}/{pub[1]:.6g}  rel {e0:.2e}/{e1:.2e}",
                  flush=True)
            if e0 > RES_RTOL or e1 > RES_RTOL:
                gate_fail.append(f"{gname} {w} residual")
        if row["pick"] != blk_pub:
            gate_fail.append(f"{gname} pick {row['pick']}!={blk_pub}")
        if row["call"] != call_pub:
            gate_fail.append(f"{gname} CALL {row['call']}!={call_pub}")
        if mac_pub is not None:
            d = abs(row["mac_best"] - mac_pub)
            print(f"  GATE_REAL MAC vs {mac_pub:.4f}: {row['mac_best']:.4f} "
                  f"|d|={d:.2e}", flush=True)
            if d > MAC_ATOL:
                gate_fail.append(f"{gname} MAC {row['mac_best']:.4f}")
        else:
            print(f"  GATE_ART full-block max-same={row['mac_best']:.4f} at "
                  f"#{row['best']}  (informational; 2434286 was 0.6281; "
                  f"NOT gated against 0.0067)", flush=True)

    if gate_fail:
        print(f"\nREADING: CP-0 GATE_FAIL  {gate_fail}", flush=True)
        raise SystemExit(4)

    rows = []
    for name, Om, note in spec["cands"]:
        hdr(f"{name}  Omega={Om}  {note}")
        row = score_one(name, Om, adv, bar, solver, FEB, mode_cls, modes)
        row["note"] = note
        rows.append(row)
        if not row.get("ok"):
            print(f"  failed: {row.get('why')}", flush=True)

    hdr("SUMMARY")
    print(f"  FE pair Hz (paper): {spec['fe_hz']}", flush=True)
    for r in rows:
        if not r.get("ok"):
            print(f"  {r['name']:8s} FAIL {r.get('why')}", flush=True)
            continue
        print(f"  {r['name']:8s} Om={r['Om']:.6f} f={r['f_hz']:.3f}  "
              f"{r['pick']:3s} {r['want_cls']:5s}  "
              f"MAC={r['mac_best']:.4f} @{r['best']}  {r['call']}",
              flush=True)

    if PAIR == "1115":
        pair_rows = [r for r in rows if r["name"] in ("ISO_A", "ISO_B")]
    else:
        pair_rows = [r for r in rows if r["name"] in ("PROD", "ISO_B")]

    hdr("READING (pre-registered -- honour the letter)")
    if any(not r.get("ok") for r in pair_rows) or len(pair_rows) < 2:
        letter = "CP-E UNRELIABLE"
        note = "An isolate failed to reconstruct. No table edit."
    else:
        calls = [r["call"] for r in pair_rows]
        n_real = sum(c == "REAL-like" for c in calls)
        if n_real == 0:
            letter = "CP-D BOTH_ART"
            note = ("Neither isolate is MAC-REAL. Footnote: MAC run; both "
                    "ARTIFACT-like at 0.744. Do not split the table.")
        elif n_real == 1:
            letter = "CP-C SPLIT"
            note = ("One REAL-like, one ARTIFACT-like. Detector still lists "
                    "one root. Do not split the published table.")
        else:
            b0, b1 = pair_rows[0]["best"], pair_rows[1]["best"]
            if b0 != b1:
                letter = "CP-A TWO_SHAPES"
                note = ("Both MAC-REAL, distinct best FE modes. Two shapes; "
                        "production still one root. Footnote can say so.")
            else:
                letter = "CP-B ONE_SHAPE"
                note = ("Both MAC-REAL against the same FE mode. One shape, "
                        "two frequencies. Table stays one listing.")
    print(f"\nREADING: {letter}", flush=True)
    print(f"  {note}", flush=True)
    print("Diagnostic only -- no SOLVER_VERSION, no silent table rewrite.",
          flush=True)
    print(f"elapsed_s={time.time()-t_all:.1f}", flush=True)


if __name__ == "__main__":
    main()
