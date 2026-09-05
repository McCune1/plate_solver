# -*- coding: utf-8 -*-
"""
probe_oop_shi2014_table3_t3_local_v1.py -- DIAGNOSTIC ONLY (no SOLVER_VERSION impact).

Hypothesis (one sentence)
-------------------------
A finer local scan around Shi2014 Table 3 target 3's predicted native
location (Ω=1.40168536) finds a deeper σ_min dip than the coarse 13-point
±15% grid did, clearing the -3.5 PASS bar -- the same mechanism as job
2315631, where a 6x-finer rescan recovered a genuine dip the default coarse
step had stepped over.

Pre-registered criteria
-----------------------
  - PASS if: a scan denser than the existing 0.035-Ω coarse step, confined
    to a TIGHT window around the predicted center (Ω in [1.3316, 1.4718],
    i.e. ±5% -- deliberately narrower than the original ±15% so it cannot
    wander into the unrelated neighboring feature at Ω=1.22647 found in the
    original coarse run), finds log_A <= -3.5 at some point.
  - FAIL/kill if: the densest feasible scan in that window still bottoms
    out shallower than -3.5 everywhere (e.g. worse than -3.4) -- in that
    case target 3 is confirmed a genuine sub-threshold near-miss, not a
    resolution artifact, and no paper text changes.

Geometry / BC / Ω (from hand-off; differs from skill baseline)
-------------------------------------------------------------
  R0_2B = 1.166666667 (Ri/Ro = 0.4)
  TWO_THETA_OVER_PI = 0.5
  IsotropicMaterial(E=210e9, nu=0.30, rho=7800.0)
  FreeFreeOOP, n_dofs=20, xmax=18.0, M=80, n_quad=30, dps=40
  Om_lit = Om * 27.4155677808 (conversion factor, already confirmed)
  lit target 3 = 38.428

Existing coarse-grid data around target 3 (already sampled -- do NOT
re-evaluate these exact points; refine BETWEEN them):
  Om=1.36664323 log_A=-1.329 (2.5% low side)
  Om=1.40168536 log_A=-3.264 (predicted center, gap=1.90 decades over
                                 local background at coarse resolution)
  Om=1.43672749 log_A=-1.922 (2.5% high side)
  Om=1.22647469 log_A=-3.309 (12.5% away -- UNRELATED; OUT OF SCOPE)

Evaluation function (verbatim from job 2410389 / hand-off; do not redesign)
---------------------------------------------------------------------------
  log_A = float(solver.sigma_min(Om, sel))  # already log10(sigma_min)
"""
from __future__ import annotations
import os
import sys
import time
import math

sys.path.insert(0, os.environ.get("PKG_PATH", "."))
os.environ.setdefault("DPS", "40")

from mpmath import mp
mp.dps = int(os.environ["DPS"])

import numpy as np
import plate_solver as ps
from plate_solver import (
    OutOfPlaneSolver, IsotropicMaterial, make_geometry, FreeFreeOOP,
)

try:
    from plate_solver.detectors import full_search, select_fill
except ImportError:
    from plate_solver import full_search, select_fill

EXPECT_VER = os.environ.get("EXPECT_SOLVER_VERSION", "2026-07-10.s10")

# Geometry for Ri/Ro = 0.4 (identical to probe_oop_shi2014_table3_r04_v*)
EPS = 0.4
R0_2B = (1.0 + EPS) / (2.0 * (1.0 - EPS))  # 1.166666...
TWO_THETA = 0.5
N_DOFS = 20
XMAX = 18.0
M = 80
N_QUAD = 30

FACTOR = 27.4155677808
LIT = 38.428
NAT_PRED = 1.40168536
WINDOW_LO = 1.3316   # ±5% of NAT_PRED
WINDOW_HI = 1.4718
DEPTH_THR = -3.5

# Known coarse points (reference only; never re-evaluated)
COARSE_KNOWN = [
    (1.36664323, -1.329),
    (1.40168536, -3.264),
    (1.43672749, -1.922),
]


def hdr(s: str) -> None:
    print("\n" + "=" * 78 + f"\n  {s}\n" + "=" * 78, flush=True)


def eval_point(solver, Om, N_DOFS=20, XMAX=18.0):
    """
    Verbatim working evaluation from job 2410389 / hand-off.
    log_A is ALREADY log10(sigma_min); do not re-log.
    """
    raw = full_search(solver.fast, Om, xmax=XMAX)
    sel, n_ret = select_fill(raw, N_DOFS)
    if n_ret < max(4, N_DOFS // 2):
        return dict(log_A=float("nan"), n_ret=n_ret, n_br=len(raw), ok=False)
    try:
        log_A = float(solver.sigma_min(Om, sel))
    except Exception as exc:
        print(f"    sigma_min failed at Om={Om:.8f}: {exc}", flush=True)
        return dict(log_A=float("nan"), n_ret=n_ret, n_br=len(raw), ok=False)
    return dict(log_A=log_A, n_ret=n_ret, n_br=len(raw), ok=True)


def main():
    t0 = time.time()
    print(f"probe_oop_shi2014_table3_t3_local_v1  start={time.strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"EXPECT_SOLVER_VERSION={EXPECT_VER}  DPS={mp.dps}")
    print(f"SOLVER_VERSION={getattr(ps, 'SOLVER_VERSION', '?')}")
    if getattr(ps, "SOLVER_VERSION", None) != EXPECT_VER:
        print("PREFLIGHT WARN: version mismatch (continuing)")

    mat = IsotropicMaterial(E=210e9, nu=0.30, rho=7800.0)
    geom = make_geometry(R0_2B, TWO_THETA)
    Ri = float(geom.R_i)
    Ro = float(geom.R_o)
    print(f"  geometry R0_2B={R0_2B:.12g}  Ri={Ri:.10g}  Ro={Ro:.10g}  Ri/Ro={Ri/Ro:.10g}")
    if abs(Ri / Ro - 0.4) > 1e-6:
        print("FATAL: Ri/Ro != 0.4")
        raise SystemExit(3)

    solver = OutOfPlaneSolver(geom, mat, M=M, n_quad=N_QUAD, boundary=FreeFreeOOP())

    # Sanity: Omega=0.50 should not be sentinel (post double-log fix)
    hdr("SANITY (Omega=0.50)")
    res0 = eval_point(solver, 0.50)
    print(f"  log_A={res0['log_A']:+.6f}  n_br={res0['n_br']}  n_ret={res0['n_ret']}")
    if res0["log_A"] == -300.0 or (not math.isnan(res0["log_A"]) and res0["log_A"] < -100):
        print("FATAL: official sigma_min still returns sentinel at Omega=0.50")
        raise SystemExit(3)
    print("  sanity PASS")

    hdr("REFERENCE: known coarse points (NOT re-evaluated)")
    for om, la in COARSE_KNOWN:
        om_lit = om * FACTOR
        rel = abs(om_lit - LIT) / LIT * 100.0
        print(f"  Om={om:.8f}  log_A={la:+.3f}  Om_lit={om_lit:.5f}  rel%={rel:.3f}")

    # ------------------------------------------------------------------
    # Stage 1: denser linspace over the tight ±5% window
    # ------------------------------------------------------------------
    N_STAGE1 = 15
    oms1 = np.linspace(WINDOW_LO, WINDOW_HI, N_STAGE1)
    # Guard against accidental re-evaluation of the exact coarse centers
    # by a tiny jitter if any grid point lands exactly on a known Om.
    known_set = {round(o, 8) for o, _ in COARSE_KNOWN}
    oms1 = [float(o) if round(float(o), 8) not in known_set else float(o) + 1e-7
            for o in oms1]

    hdr(f"STAGE 1: {N_STAGE1}-pt linspace [{WINDOW_LO:.4f}, {WINDOW_HI:.4f}]")
    stage1 = []
    best_log = 0.0
    best_om = None
    for j, om in enumerate(oms1):
        t1 = time.time()
        res = eval_point(solver, om)
        dt = time.time() - t1
        om_lit = om * FACTOR
        rel = abs(om_lit - LIT) / LIT * 100.0
        la = res["log_A"]
        print(f"  [{j+1:2d}/{N_STAGE1}] Om={om:.8f}  log_A={la:+8.3f}  "
              f"n_br={res['n_br']:3d} n_ret={res['n_ret']:2d}  "
              f"Om_lit={om_lit:10.5f}  rel%={rel:6.3f}  dt={dt:.1f}s", flush=True)
        stage1.append(dict(om=om, log_A=la, n_ret=res["n_ret"], n_br=res["n_br"],
                           om_lit=om_lit, rel=rel, ok=res["ok"]))
        if res["ok"] and not math.isnan(la) and la < best_log:
            best_log = la
            best_om = om

    print(f"\n  Stage-1 deepest: Om={best_om}  log_A={best_log:+.4f}")

    # ------------------------------------------------------------------
    # Stage 2: denser local grid around the Stage-1 minimum
    # Half-width ~0.012 (covers ~1/3 of the coarse step) so we cannot
    # wander out of the ±5% window.
    # ------------------------------------------------------------------
    HALF_W = 0.012
    lo2 = max(WINDOW_LO, best_om - HALF_W)
    hi2 = min(WINDOW_HI, best_om + HALF_W)
    N_STAGE2 = 11
    oms2 = np.linspace(lo2, hi2, N_STAGE2)
    # Again avoid exact known coarse points
    oms2 = [float(o) if round(float(o), 8) not in known_set else float(o) + 1e-7
            for o in oms2]

    hdr(f"STAGE 2: {N_STAGE2}-pt denser grid around Stage-1 min  [{lo2:.6f}, {hi2:.6f}]")
    stage2 = []
    best2_log = 0.0
    best2_om = None
    for j, om in enumerate(oms2):
        t1 = time.time()
        res = eval_point(solver, om)
        dt = time.time() - t1
        om_lit = om * FACTOR
        rel = abs(om_lit - LIT) / LIT * 100.0
        la = res["log_A"]
        print(f"  [{j+1:2d}/{N_STAGE2}] Om={om:.8f}  log_A={la:+8.3f}  "
              f"n_br={res['n_br']:3d} n_ret={res['n_ret']:2d}  "
              f"Om_lit={om_lit:10.5f}  rel%={rel:6.3f}  dt={dt:.1f}s", flush=True)
        stage2.append(dict(om=om, log_A=la, n_ret=res["n_ret"], n_br=res["n_br"],
                           om_lit=om_lit, rel=rel, ok=res["ok"]))
        if res["ok"] and not math.isnan(la) and la < best2_log:
            best2_log = la
            best2_om = om

    print(f"\n  Stage-2 deepest: Om={best2_om}  log_A={best2_log:+.4f}")

    # Optional Stage 3: one more tight 7-pt if Stage-2 improved past -3.2
    # (budget still < 35 total)
    best_final_log = best2_log
    best_final_om = best2_om
    if best2_log < -3.2 and best2_om is not None:
        HALF_W3 = 0.004
        lo3 = max(WINDOW_LO, best2_om - HALF_W3)
        hi3 = min(WINDOW_HI, best2_om + HALF_W3)
        N_STAGE3 = 7
        oms3 = np.linspace(lo3, hi3, N_STAGE3)
        oms3 = [float(o) if round(float(o), 8) not in known_set else float(o) + 1e-7
                for o in oms3]
        hdr(f"STAGE 3: {N_STAGE3}-pt ultra-local around Stage-2 min  [{lo3:.6f}, {hi3:.6f}]")
        for j, om in enumerate(oms3):
            t1 = time.time()
            res = eval_point(solver, om)
            dt = time.time() - t1
            om_lit = om * FACTOR
            rel = abs(om_lit - LIT) / LIT * 100.0
            la = res["log_A"]
            print(f"  [{j+1:2d}/{N_STAGE3}] Om={om:.8f}  log_A={la:+8.3f}  "
                  f"n_br={res['n_br']:3d} n_ret={res['n_ret']:2d}  "
                  f"Om_lit={om_lit:10.5f}  rel%={rel:6.3f}  dt={dt:.1f}s", flush=True)
            if res["ok"] and not math.isnan(la) and la < best_final_log:
                best_final_log = la
                best_final_om = om
        print(f"\n  Stage-3 deepest: Om={best_final_om}  log_A={best_final_log:+.4f}")

    # ------------------------------------------------------------------
    # READING (pre-registered)
    # ------------------------------------------------------------------
    hdr("READING (pre-registered)")
    print(f"  Tight window: [{WINDOW_LO:.4f}, {WINDOW_HI:.4f}]  (±5% of {NAT_PRED})")
    print(f"  Depth threshold for PASS: log_A <= {DEPTH_THR}")
    print(f"  Deepest found: Om={best_final_om}  log_A={best_final_log:+.4f}")
    if best_final_om is not None:
        om_lit = best_final_om * FACTOR
        rel = abs(om_lit - LIT) / LIT * 100.0
        print(f"  Mapped Om_lit={om_lit:.5f}  rel% to lit={LIT} = {rel:.3f}%")

    if best_final_log <= DEPTH_THR:
        verdict = "PASS"
        note = ("Finer local scan recovered a dip deeper than -3.5 inside the "
                "tight ±5% window. Target 3 is a resolution artifact of the "
                "original coarse 13-pt grid (same mechanism as job 2315631).")
    else:
        verdict = "FAIL / KILL"
        note = ("Densest feasible scan still bottoms out shallower than -3.5 "
                "everywhere in the tight window. Target 3 is confirmed a "
                "genuine sub-threshold near-miss; no paper text changes.")

    print(f"\n  VERDICT: {verdict}")
    print(f"  {note}")
    print(f"\n  wall_time={time.time()-t0:.1f}s")
    print("end of probe_oop_shi2014_table3_t3_local_v1")


if __name__ == "__main__":
    main()