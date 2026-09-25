# e31 = -4.1 COPY (LESSONS Sec 18.246) of probe_piezo_p4_ff_n0_overtones_wide_2026-09-16.py.
# Written by migrate_e31m.py. Changes, all fixed before running:
#   * E31_TRUE = -4.1 (Liu/Duan Table 1 as printed, Sec 18.244);
#   * e31 inverse bracket [1, 10] -> [-10, -1] (e31bar = e31 - 8.95 is
#     monotone there; the reflection twin sits at 22.0);
#   * relative uncertainties divide by |E31_TRUE|;
#   * coupling-magnitude sanity windows (quantities ~ e31bar^2) scaled
#     by (13.05/4.85)^2 = 7.239;
#   * output files carry the suffix _e31m.
# Anchor gates that pin +4.1 manuscript numbers are reported against
# those numbers and are NOT expected to pass at -4.1.
# -*- coding: utf-8 -*-
"""
probe_piezo_p4_ff_n0_overtones_wide_2026-09-16.py

Paper 4 lever 1 FOLLOW-UP: wider n=0 F-F OC overtone scan. Does NOT
re-score job 2494865 (that job's G_found bar of >=2 roots in
(800, 4000) rad/s stays FAIL). This probe asks a different question:
are there further n=0 overtones above 4000?

Job 2494865 found one overtone in (800, 4000):
  elastic 3504.952818 rad/s, OC 3509.763208 rad/s, split +0.137%.

PRE-REGISTERED (this follow-up only):
  G_n1      n=1 oc_ff_det == elastic_det at a sample omega (identity).
  G_anchor  the 3505/3510 pair reappears (each within 1 rad/s of
            job 2494865). Confirms the wide scan still sees the
            already-measured overtone.
  G_more    at least one additional n=0 pair with BOTH roots > 4000.
            If none: that is a result (no second overtone in this
            window, or Bessel breakdown) -- do not retune SCAN_HI.
  G_pair    every paired overtone has OC > elastic (stiffening).
  G2        reported, not gated, on every pair.

Scan: (800, 15000) rad/s, 400 steps (~35 rad/s). n!=0 is still not
this lever. No SOLVER_VERSION bump. Linear-only production path.
Writes piezo_p4_ff_n0_overtones_wide_results_e31m.json -- does not
overwrite the 2494865 JSON.
"""
from __future__ import annotations

import json
import math
import os
import sys
import time

sys.path.insert(0, os.environ.get(
    "PKG_PATH", os.path.join(os.path.dirname(__file__), "..")))
import plate_solver as ps  # noqa: E402
from plate_solver.piezo_solver import PiezoOutOfPlaneSolver  # noqa: E402

EXPECT_SOLVER_VERSION = os.environ.get("EXPECT_SOLVER_VERSION", "2026-07-10.s10")
DPS = int(os.environ.get("PIEZO_DPS", "60"))
ITERS = int(os.environ.get("PIEZO_ITERS", "40"))
SCAN_LO = 800.0
SCAN_HI = 15000.0
SCAN_N = int(os.environ.get("PIEZO_SCAN_N", "400"))

R_I, R_O, H = 0.1, 0.6, 0.01
E_TRUE, NU, RHO = 200e9, 0.3, 7800.0
C11E, C12E, C13E, C33E = 132e9, 71e9, 73e9, 115e9
E31_TRUE, E33 = -4.1, 14.1
X11, X33 = 7.124e-9, 5.841e-9
RHO_PZT = 7500.0
H1 = (2.0 / 12) * H
FUND_OC = 750.2394210586
ANCHOR_EL = 3504.952818
# e31=-4.1 re-anchor: OC overtone of the independent narrow scan
# probe_piezo_p4_ff_n0_overtones_e31m (3536.762179); +4.1 anchor was 3509.763208.
ANCHOR_OC = 3536.762179


def _kwargs(e31=E31_TRUE):
    return dict(
        r_i=R_I, r_o=R_O, h=H, E=E_TRUE, nu=NU, rho=RHO, h1=H1,
        C11E=C11E, C12E=C12E, C13E=C13E, C33E=C33E,
        e31=e31, e33=E33, X11=X11, X33=X33, rho_pzt=RHO_PZT, dps=DPS,
    )


def _sign_flips(detfn, lo, hi, nstep):
    step = (hi - lo) / nstep
    xs = [lo + i * step for i in range(nstep + 1)]
    prev_x = xs[0]
    try:
        prev_s = detfn(prev_x).real
    except Exception as e:
        print("scan start failed at %.3f: %s" % (prev_x, e), flush=True)
        return []
    flips = []
    n_fail = 0
    for x in xs[1:]:
        try:
            s = detfn(x).real
        except Exception:
            n_fail += 1
            continue
        if prev_s != 0 and s != 0 and (prev_s > 0) != (s > 0):
            flips.append((prev_x, x))
        prev_x, prev_s = x, s
    if n_fail:
        print("scan: %d det evals failed (Bessel/overflow); flips kept"
              % n_fail, flush=True)
    return flips


def main():
    t_all = time.time()
    print("job start: dps=%s iters=%s scan=[%.0f,%.0f]/%d"
          % (DPS, ITERS, SCAN_LO, SCAN_HI, SCAN_N), flush=True)
    print("follow-up to 2494865; does not re-score that job's G_found",
          flush=True)
    if ps.SOLVER_VERSION != EXPECT_SOLVER_VERSION:
        print("PREFLIGHT FAIL: SOLVER_VERSION=%r != %r"
              % (ps.SOLVER_VERSION, EXPECT_SOLVER_VERSION), flush=True)
        sys.exit(2)
    if not hasattr(PiezoOutOfPlaneSolver, "oc_ff_bisect"):
        print("PREFLIGHT FAIL: no oc_ff_bisect -- push piezo_solver.py",
              flush=True)
        sys.exit(2)
    print("preflight OK: SOLVER_VERSION=%s" % ps.SOLVER_VERSION, flush=True)

    s = PiezoOutOfPlaneSolver(**_kwargs())
    g_n1 = (s.oc_ff_det(ANCHOR_EL, 1) == s.elastic_det(ANCHOR_EL, 1))
    print("G_n1: n=1 oc==elastic %s" % g_n1, flush=True)

    t0 = time.time()
    el_flips = _sign_flips(lambda w: s.elastic_det(w, 0), SCAN_LO, SCAN_HI, SCAN_N)
    oc_flips = _sign_flips(lambda w: s.oc_ff_det(w, 0), SCAN_LO, SCAN_HI, SCAN_N)
    print("scan: %d elastic flips, %d OC flips (%.1fs)"
          % (len(el_flips), len(oc_flips), time.time() - t0), flush=True)

    el_roots, oc_roots = [], []
    for lo, hi in el_flips:
        try:
            el_roots.append(s.elastic_bisect(lo, hi, 0, iters=ITERS))
        except Exception as e:
            print("elastic bisect failed on [%.3f, %.3f]: %s" % (lo, hi, e),
                  flush=True)
    for lo, hi in oc_flips:
        try:
            oc_roots.append(s.oc_ff_bisect(lo, hi, 0, iters=ITERS))
        except Exception as e:
            print("OC bisect failed on [%.3f, %.3f]: %s" % (lo, hi, e),
                  flush=True)
    print("elastic roots: %s" % ["%.6f" % x for x in el_roots], flush=True)
    print("OC roots:      %s" % ["%.6f" % x for x in oc_roots], flush=True)

    n_pair = min(len(el_roots), len(oc_roots))
    g_anchor = False
    g_more = False
    g_pair = n_pair > 0
    rows = []
    for i in range(n_pair):
        w_el, w_oc = el_roots[i], oc_roots[i]
        if not (w_oc > w_el > FUND_OC):
            g_pair = False
        if abs(w_el - ANCHOR_EL) < 1.0 and abs(w_oc - ANCHOR_OC) < 1.0:
            g_anchor = True
        if w_el > 4000.0 and w_oc > 4000.0:
            g_more = True
        split = (w_oc - w_el) / w_el
        k2 = 1.0 - (w_el / w_oc) ** 2
        de = 0.01 * E31_TRUE
        s_plus = PiezoOutOfPlaneSolver(**_kwargs(E31_TRUE + de))
        s_minus = PiezoOutOfPlaneSolver(**_kwargs(E31_TRUE - de))
        lo, hi = 0.97 * w_oc, 1.05 * w_oc
        try:
            w_plus = s_plus.oc_ff_bisect(lo, hi, 0, iters=ITERS)
            w_minus = s_minus.oc_ff_bisect(lo, hi, 0, iters=ITERS)
            dw_de = (w_plus - w_minus) / (2.0 * de)
            g2 = ((1e-4 * w_oc / abs(dw_de)) / abs(E31_TRUE) * 100.0
                  if dw_de else float("inf"))
        except Exception as e:
            dw_de, g2 = float("nan"), float("nan")
            print("G2 failed on overtone %d: %s" % (i + 1, e), flush=True)
        row = dict(
            index=i + 1, omega_el=w_el, omega_oc=w_oc,
            f_el_hz=w_el / (2 * math.pi), f_oc_hz=w_oc / (2 * math.pi),
            split=split, split_pct=100.0 * split, k2=k2,
            above_4000=(w_el > 4000.0 and w_oc > 4000.0),
            domega_de31=dw_de, g2_pct_001=g2,
        )
        rows.append(row)
        print("overtone %d: el=%.6f oc=%.6f split=%.4f%% k2=%.6f "
              "G2@0.01%%=%.3f%% above4000=%s"
              % (i + 1, w_el, w_oc, 100.0 * split, k2, g2,
                 row["above_4000"]), flush=True)

    all_pass = g_n1 and g_anchor and g_more and g_pair
    out_dir = os.path.dirname(__file__) or "."
    results = dict(
        dps=DPS, solver_version=ps.SOLVER_VERSION,
        scan_lo=SCAN_LO, scan_hi=SCAN_HI, scan_n=SCAN_N,
        parent_job="2494865",
        g_n1=g_n1, g_anchor=g_anchor, g_more=g_more, g_pair=g_pair,
        n_overtones=n_pair, all_pass=all_pass, rows=rows,
        elapsed_s=time.time() - t_all,
    )
    path = os.path.join(out_dir, "piezo_p4_ff_n0_overtones_wide_results_e31m.json")
    with open(path, "w") as f:
        json.dump(results, f, indent=2)
    print("SENTINEL " + ("PASS_ALL" if all_pass else "FAIL_ALL")
          + " n1=%s anchor=%s more=%s pair=%s n=%d"
          % (g_n1, g_anchor, g_more, g_pair, n_pair), flush=True)


if __name__ == "__main__":
    main()
