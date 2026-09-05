# -*- coding: utf-8 -*-
"""
probe_ip_pairing_fallback_bias_v1.py -- DIAGNOSTIC ONLY.
No SOLVER_VERSION bump. Does not edit ip_mac_bar.py.

Question: does the Cartesian-midside pairing fallback change a MAC call
relative to exact polar pairing?

Control geometry: r200/a100, PHI=180, where exact pairing is known to
succeed (ip_mac_bar._partner_index returns missing=0). Compare
  (exact)   production _partner_index
  (forced)  Cartesian nearest-neighbour only, skipping exact keys
on (a) FE SYMM/ANTI labels of all 57 elastic modes, (b) the 8 calibration
candidates' class-safe MAC calls at threshold 0.744.

Uses the existing Rank 18 dump ip_mac_fullblock_r200a100_m{4..60}_mesh96.txt
-- no new ANSYS job.

-----------------------------------------------------------------------------
PRE-REGISTERED INTERPRETATION
-----------------------------------------------------------------------------
PB-0 GATE_FAIL
    Missing dump, SOLVER_VERSION, mac() identity, or GATE_REAL t1 /
    GATE_ART t5 MAC disagrees with published 0.9993 / 0.1969 by >5e-4
    under production (exact) pairing. VOID.
PB-A NO_BIAS
    Zero FE-class flips among the 57, and zero of the 8 calibration
    MAC calls flip across 0.744. -> one sentence: fallback was checked
    for verdict bias at the exact-pairing geometry and found none.
PB-B CLASS_ONLY
    Some FE SYMM/ANTI labels move, but none of the 8 MAC calls flip.
    Report the class churn; the bar's call is stable. Optional sentence,
    not a paper defect.
PB-C CALL_FLIP
    At least one of the 8 calibration calls flips REAL-like <-> ARTIFACT-like
    when fallback is forced. -> do not claim the fallback is verdict-neutral.
    Do not retune 0.744. Investigate before any SM wording that the
    primary 90-deg calibration "used the fallback" is read as a threat.

Honour the printed letter. Do not move 0.744.
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
RULE_DELTA = float(os.environ.get("RULE_DELTA", "0.001"))
MAC_ATOL = 5e-4
PHI_DEG = 180.0
THETA = float(np.pi) / 2.0
PREFIX = "ip_mac_fullblock_r200a100"

CANDS = [
    ("t1", 0.304284, "REAL", 0.9993, "q1"),
    ("t2", 0.915784, "REAL", None, None),
    ("t3", 1.420335, "REAL", None, None),
    ("t4", 2.286070, "REAL", None, None),
    ("t5", 0.117785, "ARTIFACT", 0.1969, "q0"),
    ("t6", 0.443610, "ARTIFACT", None, None),
    ("t7", 0.507596, "ARTIFACT", None, None),
    ("t8", 0.590794, "ARTIFACT", None, None),
]


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
        print("READING: PB-0 GATE_FAIL  (adv missing)", flush=True)
        raise SystemExit(2)
    spec = importlib.util.spec_from_file_location("adv_pb", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    mod.THETA = float(THETA)
    return mod


def partner_cartesian_only(r, th, phi_deg, bar):
    """Force the Cartesian fallback, ignoring exact polar keys."""
    r = np.asarray(r, dtype=float)
    th = np.asarray(th, dtype=float)
    n = len(r)
    th_rad = np.deg2rad(th)
    xy = np.column_stack((r * np.cos(th_rad), r * np.sin(th_rad)))
    th_img = np.deg2rad(float(phi_deg) - th)
    xy_img = np.column_stack((r * np.cos(th_img), r * np.sin(th_img)))
    span = float(np.ptp(r)) if n else 0.0
    max_dist = max(0.5 * span / 96.0, 1.0e-4)
    try:
        from scipy.spatial import cKDTree
        dist, idx = cKDTree(xy).query(xy_img, k=1)
        dist = np.asarray(dist, dtype=float)
        idx = np.asarray(idx, dtype=np.int64)
    except Exception:
        dist = np.full(n, np.inf)
        idx = np.zeros(n, dtype=np.int64)
        for i in range(n):
            d2 = ((xy[:, 0] - xy_img[i, 0]) ** 2
                  + (xy[:, 1] - xy_img[i, 1]) ** 2)
            j = int(np.argmin(d2))
            dist[i] = float(np.sqrt(d2[j]))
            idx[i] = j
    partner = np.empty(n, dtype=np.int64)
    missing = 0
    for i in range(n):
        if dist[i] <= max_dist:
            partner[i] = int(idx[i])
        else:
            partner[i] = -1
            missing += 1
    return partner, missing


def classify_with(r, th, ux, uy, partner, missing, bar, miss_frac=0.02):
    n = len(r)
    if n == 0 or (missing / float(n)) > miss_frac:
        return "AMBIGUOUS", missing
    ur, ut = bar.polar(ux, uy, th)
    ratios = bar.parity_ratios(ur, ut, partner)
    cls, _s, _a = bar.classify_parity(*ratios)
    return cls, missing


def score_one(name, Om, adv, bar, solver, FEB, mode_hz, mode_cls, modes):
    n_dofs = adv.N_DOFS
    brs, _x, st, cnt, _n = adv.get_filled_brs(solver, Om, n_dofs)
    if st == "UNRELIABLE":
        return None
    sh = adv.reconstruct_shared(solver, Om, brs, n_dofs=n_dofs)
    if not sh.get("ok"):
        return None
    blocks = {}
    for w in ("q0", "q1"):
        r = adv.block_reconstruction(solver, sh, w)
        if r.get("ok"):
            blocks[w] = r
    s0, s1, _c = adv.sigma_pair_at(solver, Om)
    m0, m1, _c0 = adv.sigma_pair_at(solver, Om - RULE_DELTA)
    p0, p1, _c1 = adv.sigma_pair_at(solver, Om + RULE_DELTA)
    if None in (s0, s1, m0, m1, p0, p1) or "q0" not in blocks or "q1" not in blocks:
        return None
    pick, how = bar.rule_localmin_resid(
        s0, s1, m0, p0, m1, p1, blocks["q0"]["rTyy"], blocks["q1"]["rTyy"],
    )
    want = bar.block_class(pick)
    rp, tp, Ur, Uth = adv.build_grid_from_cfull(
        solver, sh, blocks[pick]["c_full"])
    ux_s, uy_s, bad = adv.solver_field_at_fe(rp, tp, Ur, Uth, FEB[modes[0]])
    keep = ~bad
    bvec = np.concatenate([ux_s[keep], uy_s[keep]])
    Ablk = np.vstack([np.concatenate([FEB[m][2][keep], FEB[m][3][keep]])
                      for m in modes])
    macs = adv.mac_matrix(Ablk, bvec)
    mac_by = {m: float(macs[i]) for i, m in enumerate(modes)}
    same = [m for m in modes if mode_cls.get(m) == want]
    if not same:
        return dict(name=name, pick=pick, want=want, mac=float("nan"),
                    best=None, call="INVALID")
    mac_best, best = max((mac_by[m], m) for m in same)
    return dict(name=name, pick=pick, want=want, mac=mac_best, best=best,
                call=bar.classify_mac(mac_best), how=how)


def main():
    t_all = time.time()
    hdr("probe_ip_pairing_fallback_bias_v1  r200/a100 PHI=180")
    import ip_mac_bar as bar
    bar.ensure_free_free_worker_bc()
    import plate_solver as ps
    print(f"  SOLVER_VERSION={ps.SOLVER_VERSION!r}", flush=True)
    if ps.SOLVER_VERSION != EXPECT_VER:
        print("READING: PB-0 GATE_FAIL  (SOLVER_VERSION)", flush=True)
        raise SystemExit(3)

    adv = load_adv()
    a = np.random.default_rng(18).normal(size=200)
    b = np.random.default_rng(19).normal(size=200)
    if abs(bar.mac(a, b) - adv.mac(a, b)) > 1e-15:
        print("READING: PB-0 GATE_FAIL  (mac)", flush=True)
        raise SystemExit(4)

    FEB = {}
    for m in range(4, 61):
        p = find_file(f"{PREFIX}_m{m}_mesh96.txt")
        if p is None:
            print(f"READING: PB-0 GATE_FAIL  (missing {PREFIX}_m{m})",
                  flush=True)
            raise SystemExit(2)
        FEB[m] = adv.load_fe_eigvec(p)
    freq_path = find_file(f"{PREFIX}_freqs.txt")
    FREQ = adv.load_freqs(freq_path) if freq_path else {}
    modes = list(range(4, 61))
    mode_hz = {m: float(FREQ[m]) for m in modes if m in FREQ}
    print(f"  dumps 57/57  n_nodes={len(FEB[4][0])}", flush=True)

    cls_exact = {}
    cls_cart = {}
    n_exact_miss = 0
    flips = []
    for m in modes:
        r, th, ux, uy = FEB[m]
        p_ex, miss_ex = bar._partner_index(r, th, phi_deg=PHI_DEG)
        p_ca, miss_ca = partner_cartesian_only(r, th, PHI_DEG, bar)
        c_ex, _ = classify_with(r, th, ux, uy, p_ex, miss_ex, bar)
        c_ca, _ = classify_with(r, th, ux, uy, p_ca, miss_ca, bar)
        cls_exact[m] = c_ex
        cls_cart[m] = c_ca
        n_exact_miss += miss_ex
        if c_ex != c_ca:
            flips.append((m, c_ex, c_ca, miss_ex, miss_ca))
    print(f"  exact-pairing total missing nodes over 57 modes: {n_exact_miss} "
          f"(expect 0 at PHI=180)", flush=True)
    print(f"  FE class flips exact vs forced-cartesian: {len(flips)}",
          flush=True)
    for row in flips[:12]:
        print(f"    FE#{row[0]}  exact={row[1]}  cart={row[2]}  "
              f"miss_ex={row[3]} miss_ca={row[4]}", flush=True)

    mat = ps.IsotropicMaterial(E=210e9, nu=0.30, rho=7800.0)
    geom = ps.make_geometry(2.0, 1.0)
    solver = ps.InPlaneSolver(geom, mat, M=80, n_quad=30,
                              boundary=ps.FreeFreeIP())

    gate_fail = []
    calls_ex, calls_ca = [], []
    for name, Om, lab, mac_pub, blk_pub in CANDS:
        hdr(f"{name}  Omega={Om}  {lab}")
        rex = score_one(name, Om, adv, bar, solver, FEB, mode_hz,
                        cls_exact, modes)
        rca = score_one(name, Om, adv, bar, solver, FEB, mode_hz,
                        cls_cart, modes)
        if rex is None or rca is None:
            gate_fail.append(f"{name} reconstruct")
            continue
        print(f"  exact  pick={rex['pick']} {rex['want']} "
              f"MAC={rex['mac']:.4f} @{rex['best']} {rex['call']}", flush=True)
        print(f"  cart   pick={rca['pick']} {rca['want']} "
              f"MAC={rca['mac']:.4f} @{rca['best']} {rca['call']}", flush=True)
        calls_ex.append(rex)
        calls_ca.append(rca)
        if mac_pub is not None:
            d = abs(rex["mac"] - mac_pub)
            print(f"  GATE vs published {mac_pub:.4f}: |d|={d:.2e}", flush=True)
            if d > MAC_ATOL:
                gate_fail.append(f"{name} MAC {rex['mac']:.4f}!={mac_pub:.4f}")

    if gate_fail:
        print(f"\nREADING: PB-0 GATE_FAIL  {gate_fail}", flush=True)
        raise SystemExit(4)

    call_flips = []
    for a, b in zip(calls_ex, calls_ca):
        if a["call"] != b["call"]:
            call_flips.append((a["name"], a["call"], b["call"]))

    hdr("READING (pre-registered -- honour the letter)")
    if call_flips:
        letter = "PB-C CALL_FLIP"
        note = f"MAC call flips: {call_flips}. Fallback is not verdict-neutral."
    elif flips:
        letter = "PB-B CLASS_ONLY"
        note = (f"{len(flips)} FE class labels moved; 0 of 8 MAC calls flipped. "
                "Bar is stable.")
    else:
        letter = "PB-A NO_BIAS"
        note = ("0 FE class flips, 0 MAC-call flips. Fallback checked for "
                "verdict bias at the exact-pairing geometry: none.")
    print(f"\nREADING: {letter}", flush=True)
    print(f"  {note}", flush=True)
    print(f"elapsed_s={time.time()-t_all:.1f}", flush=True)


if __name__ == "__main__":
    main()
