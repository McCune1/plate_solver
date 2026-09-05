# -*- coding: utf-8 -*-
"""
probe_oop_shi2014_table3_r04_v3.py -- DIAGNOSTIC ONLY (no SOLVER_VERSION impact).

Hypothesis
----------
With the double-log bug fixed (treat oop.sigma_min's return value as already
log10(sigma_min), no second conversion), the project's own standard,
previously-validated sigma_min detector reproduces job 2410239's
manual-equilibrated-SVD 5/6 recovery of Shi 2014 Table 3 (isotropic OOP,
a/b=0.4) directly -- giving a clean, citable result with no novel metric needed.

Pre-registered criteria (identical numerical bars to job 2410239, corrected detector)
------------------------------------------------------------------------------------
  - PASS:    >=5/6 targets recovered within 1%
             (log_A <= -3.5 inside the +-15% window whose mapped Om_lit is
              within 1% of the tabulated value), using the CORRECTED reading:
             log_A = float(oop.sigma_min(om, sel))  -- no second wrapper.
  - PARTIAL: 3-4/6 within 1%, or all recovered only within 1-3%.
  - FAIL:    <=2/6 within 3% or no dips near any target.

Also report the manual-SVD metric (job 2410239 harness) side-by-side at every
point so the two detectors can be compared directly.

Geometry / BC / material (identical to v1/v2)
---------------------------------------------
  R0_2B = 1.166666667 (Ri/Ro=0.4), TWO_THETA_OVER_PI=0.5
  IsotropicMaterial(E=210e9, nu=0.30, rho=7800.0)
  FreeFreeOOP, n_dofs=20, M=80, n_quad=30, xmax=18, dps=40

Targets / windows (identical to job 2410239)
--------------------------------------------
  lit=15.647  native=0.57073412
  lit=23.576  native=0.85994936
  lit=38.428  native=1.40168536
  lit=53.649  native=1.95688087
  lit=63.390  native=2.31218994
  lit=70.739  native=2.58024932
  +-15% window, 13-point grid.
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

try:
    from plate_solver.geometry import _kbar_wbar, _omega_lit
except ImportError:
    try:
        from plate_solver.core_solvers import _kbar_wbar, _omega_lit
    except ImportError:
        from plate_solver import _kbar_wbar, _omega_lit

EXPECT_VER = os.environ.get("EXPECT_SOLVER_VERSION", "2026-07-10.s10")

EPS = 0.4
R0_2B = (1.0 + EPS) / (2.0 * (1.0 - EPS))
TWO_THETA = 0.5
N_DOFS = 20
XMAX = 18.0
M = 80
N_QUAD = 30
N_GRID = 13
HALF_FRAC = 0.15
DEPTH_THR = -3.5

# Predicted natives from job 2410239 (factor = 27.4155677808)
TARGETS = [
    # (lit, native_pred)
    (15.647, 0.57073412),
    (23.576, 0.85994936),
    (38.428, 1.40168536),
    (53.649, 1.95688087),
    (63.390, 2.31218994),
    (70.739, 2.58024932),
]
FACTOR = 27.4155677808


def hdr(s: str) -> None:
    print("\n" + "=" * 78 + f"\n  {s}\n" + "=" * 78, flush=True)


def equilibrate_mp(K_mp, n, passes=2):
    """Verbatim from job 2410239."""
    K = K_mp.copy()
    for _ in range(passes):
        for i in range(n):
            s = max(abs(K[i, j]) for j in range(n))
            if s > 0:
                inv = mp.mpf(1) / s
                for j in range(n):
                    K[i, j] *= inv
        for j in range(n):
            s = max(abs(K[i, j]) for i in range(n))
            if s > 0:
                inv = mp.mpf(1) / s
                for i in range(n):
                    K[i, j] *= inv
    return K


def to_numpy(K_mp, n):
    Knp = np.empty((n, n), dtype=complex)
    for i in range(n):
        for j in range(n):
            Knp[i, j] = complex(K_mp[i, j])
    return Knp


def log10_svd_smin(K_mp, n):
    """Manual equilibrated SVD (job 2410239 metric)."""
    try:
        Keq = equilibrate_mp(K_mp, n)
        Knp = to_numpy(Keq, n)
        s = np.linalg.svd(Knp, compute_uv=False)
        smin = float(s[-1]) if len(s) else 0.0
        if smin <= 0.0:
            return -300.0, 0.0
        s2 = float(s[-2]) if len(s) >= 2 else smin
        return math.log10(smin), math.log10(s2) - math.log10(smin)
    except Exception as exc:
        print(f"    SVD failed: {exc}", flush=True)
        return float("nan"), float("nan")


def eval_both(solver, Om):
    """
    Returns:
      log_A   = corrected official (direct return of oop.sigma_min)
      log_SVD = manual equilibrated SVD
      gap     = SVD gap
      n_br, n_ret
    """
    try:
        raw = full_search(solver.fast, Om, xmax=XMAX)
        n_br = len(raw)
        sel, n_ret = select_fill(raw, N_DOFS)
        if n_ret < max(4, N_DOFS // 2):
            return dict(log_A=float("nan"), log_SVD=float("nan"), gap=float("nan"),
                        n_br=n_br, n_ret=n_ret, ok=False)

        # Corrected official reading: return value IS already log10(sigma_min)
        try:
            log_A = float(solver.sigma_min(Om, sel))
        except Exception as exc:
            print(f"    official failed: {exc}", flush=True)
            log_A = float("nan")

        # Manual SVD on the same (Om, sel)
        try:
            K_mp, n = solver._build_K_real(Om, sel, fast_scan=False)
            log_SVD, gap = log10_svd_smin(K_mp, n)
        except Exception as exc:
            print(f"    build/SVD failed: {exc}", flush=True)
            log_SVD, gap = float("nan"), float("nan")

        return dict(log_A=log_A, log_SVD=log_SVD, gap=gap,
                    n_br=n_br, n_ret=n_ret, ok=True)
    except Exception as exc:
        print(f"    eval_both failed: {exc}", flush=True)
        return dict(log_A=float("nan"), log_SVD=float("nan"), gap=float("nan"),
                    n_br=0, n_ret=0, ok=False)


def main():
    t0 = time.time()
    print(f"probe_oop_shi2014_table3_r04_v3  start={time.strftime('%Y-%m-%d %H:%M:%S')}")
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

    # Quick sanity at Omega=0.50 (expect log_A ~ -2, not -300)
    hdr("SANITY (Omega=0.50, corrected official)")
    res0 = eval_both(solver, 0.50)
    print(f"  log_A (corrected) = {res0['log_A']:+.6f}")
    print(f"  log_SVD           = {res0['log_SVD']:+.6f}  gap={res0['gap']:.2f}")
    print(f"  n_br={res0['n_br']}  n_ret={res0['n_ret']}")
    if res0["log_A"] == -300.0 or (not math.isnan(res0["log_A"]) and res0["log_A"] < -100):
        print("FATAL: corrected official still looks like a sentinel at Omega=0.50")
        raise SystemExit(3)
    print("  sanity PASS (official is no longer forced to -300)")

    hdr("LOCAL GRIDS (+-15% about each predicted native)")
    results = []
    all_points = []  # for the full side-by-side dump

    for idx, (lit, nat) in enumerate(TARGETS):
        lo = nat * (1.0 - HALF_FRAC)
        hi = nat * (1.0 + HALF_FRAC)
        print(f"\n--- target {idx+1}/6  lit={lit:.5f}  native_pred={nat:.8f} ---", flush=True)
        oms = np.linspace(lo, hi, N_GRID)
        deep_pts = []  # (om, log_A, log_SVD, gap, rel_pct)
        best_A = 0.0
        best_om_A = None

        for j, om in enumerate(oms):
            t1 = time.time()
            res = eval_both(solver, float(om))
            dt = time.time() - t1
            om_lit = float(om) * FACTOR
            rel = abs(om_lit - lit) / lit * 100.0
            print(f"  [{j+1:2d}/{N_GRID}] Om={om:.8f}  "
                  f"log_A={res['log_A']:+8.3f}  log_SVD={res['log_SVD']:+8.3f}  "
                  f"gap={res['gap']:5.2f}  n_br={res['n_br']:3d}  "
                  f"Om_lit={om_lit:10.5f}  rel%={rel:6.3f}  dt={dt:.1f}s", flush=True)
            all_points.append(dict(
                target=idx+1, lit=lit, om=float(om), log_A=res["log_A"],
                log_SVD=res["log_SVD"], gap=res["gap"], rel=rel, n_br=res["n_br"],
            ))
            if res["ok"] and not math.isnan(res["log_A"]):
                if res["log_A"] < best_A:
                    best_A = res["log_A"]
                    best_om_A = float(om)
                if res["log_A"] <= DEPTH_THR:
                    deep_pts.append((float(om), res["log_A"], res["log_SVD"],
                                     res["gap"], rel))

        recovered_1 = recovered_3 = False
        chosen_om = chosen_A = chosen_SVD = chosen_gap = chosen_rel = None
        if deep_pts:
            deep_pts.sort(key=lambda x: x[4])  # closest mapped lit
            chosen_om, chosen_A, chosen_SVD, chosen_gap, chosen_rel = deep_pts[0]
            if chosen_rel <= 1.0:
                recovered_1 = True
            if chosen_rel <= 3.0:
                recovered_3 = True
        else:
            chosen_om = best_om_A
            chosen_A = best_A
            if best_om_A is not None:
                chosen_rel = abs(best_om_A * FACTOR - lit) / lit * 100.0
            else:
                chosen_rel = float("nan")

        results.append(dict(
            lit=lit, nat=nat, chosen_om=chosen_om, chosen_A=chosen_A,
            chosen_SVD=chosen_SVD, chosen_gap=chosen_gap, chosen_rel=chosen_rel,
            recovered_1=recovered_1, recovered_3=recovered_3, n_deep=len(deep_pts),
        ))
        st = "HIT_1%" if recovered_1 else ("HIT_3%" if recovered_3 else "MISS")
        print(f"  >> chosen Om={chosen_om}  log_A={chosen_A}  rel%={chosen_rel}  [{st}]")

    # ------------------------------------------------------------------
    # Full side-by-side dump
    # ------------------------------------------------------------------
    hdr("SIDE-BY-SIDE TABLE (all 78 points)")
    print(f"{'#':>2}  {'lit':>8}  {'Om':>10}  {'log_A':>9}  {'log_SVD':>9}  "
          f"{'gap':>6}  {'rel%':>7}  {'n_br':>4}")
    for p in all_points:
        print(f"{p['target']:2d}  {p['lit']:8.3f}  {p['om']:10.6f}  "
              f"{p['log_A']:+9.3f}  {p['log_SVD']:+9.3f}  {p['gap']:6.2f}  "
              f"{p['rel']:7.3f}  {p['n_br']:4d}")

    # ------------------------------------------------------------------
    # Pre-registered reading (official corrected only)
    # ------------------------------------------------------------------
    hdr("PRE-REGISTERED READING (corrected official detector)")
    n1 = sum(1 for r in results if r["recovered_1"])
    n3 = sum(1 for r in results if r["recovered_3"])
    print(f"  recovered within 1%: {n1}/6")
    print(f"  recovered within 3%: {n3}/6")
    print()
    print(f"  {'#':>2}  {'lit':>10}  {'nat_pred':>12}  {'found_Om':>12}  "
          f"{'log_A':>9}  {'log_SVD':>9}  {'rel%':>8}  status")
    for i, r in enumerate(results):
        st = "HIT_1%" if r["recovered_1"] else ("HIT_3%" if r["recovered_3"] else "MISS")
        fo = r["chosen_om"] if r["chosen_om"] is not None else float("nan")
        la = r["chosen_A"] if r["chosen_A"] is not None else float("nan")
        ls = r["chosen_SVD"] if r["chosen_SVD"] is not None else float("nan")
        re = r["chosen_rel"] if r["chosen_rel"] is not None else float("nan")
        print(f"  {i+1:2d}  {r['lit']:10.5f}  {r['nat']:12.8f}  {fo:12.8f}  "
              f"{la:9.3f}  {ls:9.3f}  {re:8.3f}  {st}")

    print()
    if n1 >= 5:
        verdict = "PASS"
        note = (">=5/6 within 1% with corrected official detector -- "
                "reproduces the 2410239 recovery; citable with the standard metric.")
    elif n1 >= 3 or n3 >= 5:
        verdict = "PARTIAL"
        note = "3-4 within 1% or >=5 within 3% -- conversion partially successful."
    else:
        verdict = "FAIL"
        note = "<=2/6 within 3% -- corrected official does not recover the table."
    print(f"  VERDICT: {verdict}")
    print(f"  {note}")
    print(f"\n  elapsed = {time.time()-t0:.1f} s")
    print(f"probe_oop_shi2014_table3_r04_v3  end={time.strftime('%Y-%m-%d %H:%M:%S')}")


if __name__ == "__main__":
    main()