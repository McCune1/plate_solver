# -*- coding: utf-8 -*-
"""
probe_piezo_p4_ff_n0_overtones_2026-09-16.py

Paper 4 lever 1: F-F n=0 open-circuit RADIAL OVERTONES of the same
linear-potential operator as the fundamental (oc_ff_bisect /
elastic_bisect). n!=0 is NOT this lever -- it is bit-identical to
the elastic bilayer and does not identify e31.

PRE-REGISTERED:
  G_fund   n=0 elastic and OC each have at least 2 roots above the
                  fundamental in (800, 4000) rad/s (sign flips of the
                  4x4 det). The fundamental ~748-750 is excluded.
  G_pair          for each overtone, an OC root sits above its elastic
                  partner (stiffening). Pair by order, not by seeding
                  published numbers.
  G_n1            oc_ff_det(n=1) == elastic_det(n=1) at a sample omega
                  (identity; n!=0 is not an e31 observable).
  G2              reported, not gated: e31 uncertainty at 0.01% freq
                  on each overtone. Expect the split and G2 to survive
                  at higher radial wavenumber, not necessarily at the
                  same 1.98%.

No SOLVER_VERSION bump. Linear-only production path.
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
DPS = int(os.environ.get("PIEZO_DPS", "40"))
ITERS = int(os.environ.get("PIEZO_ITERS", "40"))
SCAN_LO = 800.0
SCAN_HI = 4000.0
SCAN_N = int(os.environ.get("PIEZO_SCAN_N", "160"))

R_I, R_O, H = 0.1, 0.6, 0.01
E_TRUE, NU, RHO = 200e9, 0.3, 7800.0
C11E, C12E, C13E, C33E = 132e9, 71e9, 73e9, 115e9
E31_TRUE, E33 = 4.1, 14.1
X11, X33 = 7.124e-9, 5.841e-9
RHO_PZT = 7500.0
H1 = (2.0 / 12) * H
FUND_EL = 747.9863950953
FUND_OC = 750.2394210586


def _kwargs(e31=E31_TRUE):
    return dict(
        r_i=R_I, r_o=R_O, h=H, E=E_TRUE, nu=NU, rho=RHO, h1=H1,
        C11E=C11E, C12E=C12E, C13E=C13E, C33E=C33E,
        e31=e31, e33=E33, X11=X11, X33=X33, rho_pzt=RHO_PZT, dps=DPS,
    )


def _sign_flips(detfn, lo, hi, nstep):
    step = (hi - lo) / nstep
    xs = [lo + i * step for i in range(nstep + 1)]
    prev_x, prev_s = xs[0], detfn(xs[0]).real
    flips = []
    for x in xs[1:]:
        s = detfn(x).real
        if prev_s != 0 and s != 0 and (prev_s > 0) != (s > 0):
            flips.append((prev_x, x))
        prev_x, prev_s = x, s
    return flips


def main():
    t_all = time.time()
    print("job start: dps=%s iters=%s scan=[%.0f,%.0f]/%d"
          % (DPS, ITERS, SCAN_LO, SCAN_HI, SCAN_N), flush=True)
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
    det_n1_oc = s.oc_ff_det(FUND_EL, 1)
    det_n1_el = s.elastic_det(FUND_EL, 1)
    g_n1 = (det_n1_oc == det_n1_el)
    print("G_n1: n=1 oc==elastic %s" % g_n1, flush=True)

    t0 = time.time()
    el_flips = _sign_flips(lambda w: s.elastic_det(w, 0), SCAN_LO, SCAN_HI, SCAN_N)
    oc_flips = _sign_flips(lambda w: s.oc_ff_det(w, 0), SCAN_LO, SCAN_HI, SCAN_N)
    print("scan: %d elastic flips, %d OC flips (%.1fs)"
          % (len(el_flips), len(oc_flips), time.time() - t0), flush=True)

    el_roots, oc_roots = [], []
    for lo, hi in el_flips:
        el_roots.append(s.elastic_bisect(lo, hi, 0, iters=ITERS))
    for lo, hi in oc_flips:
        oc_roots.append(s.oc_ff_bisect(lo, hi, 0, iters=ITERS))
    print("elastic roots: %s" % ["%.6f" % x for x in el_roots], flush=True)
    print("OC roots:      %s" % ["%.6f" % x for x in oc_roots], flush=True)

    n_pair = min(len(el_roots), len(oc_roots))
    g_fund = n_pair >= 2
    rows = []
    g_pair = True
    for i in range(n_pair):
        w_el, w_oc = el_roots[i], oc_roots[i]
        if not (w_oc > w_el > FUND_OC):
            g_pair = False
        split = (w_oc - w_el) / w_el
        k2 = 1.0 - (w_el / w_oc) ** 2
        de = 0.01 * E31_TRUE
        s_plus = PiezoOutOfPlaneSolver(**_kwargs(E31_TRUE + de))
        s_minus = PiezoOutOfPlaneSolver(**_kwargs(E31_TRUE - de))
        lo, hi = 0.97 * w_oc, 1.05 * w_oc
        w_plus = s_plus.oc_ff_bisect(lo, hi, 0, iters=ITERS)
        w_minus = s_minus.oc_ff_bisect(lo, hi, 0, iters=ITERS)
        dw_de = (w_plus - w_minus) / (2.0 * de)
        g2 = (1e-4 * w_oc / abs(dw_de)) / E31_TRUE * 100.0 if dw_de else float("inf")
        row = dict(
            index=i + 1, omega_el=w_el, omega_oc=w_oc,
            f_el_hz=w_el / (2 * math.pi), f_oc_hz=w_oc / (2 * math.pi),
            split=split, split_pct=100.0 * split, k2=k2,
            domega_de31=dw_de, g2_pct_001=g2,
        )
        rows.append(row)
        print("overtone %d: el=%.6f oc=%.6f split=%.4f%% k2=%.6f G2@0.01%%=%.3f%%"
              % (i + 1, w_el, w_oc, 100.0 * split, k2, g2), flush=True)

    all_pass = g_n1 and g_fund and g_pair and n_pair >= 2
    out_dir = os.path.dirname(__file__) or "."
    results = dict(
        dps=DPS, solver_version=ps.SOLVER_VERSION,
        g_n1=g_n1, g_fund=g_fund, g_pair=g_pair,
        n_overtones=n_pair, all_pass=all_pass,
        fund_el=FUND_EL, fund_oc=FUND_OC, rows=rows,
        elapsed_s=time.time() - t_all,
    )
    path = os.path.join(out_dir, "piezo_p4_ff_n0_overtones_results.json")
    with open(path, "w") as f:
        json.dump(results, f, indent=2)
    print("SENTINEL " + ("PASS_ALL" if all_pass else "FAIL_ALL")
          + " n1=%s found=%s pair=%s n=%d" % (g_n1, g_fund, g_pair, n_pair),
          flush=True)


if __name__ == "__main__":
    main()
