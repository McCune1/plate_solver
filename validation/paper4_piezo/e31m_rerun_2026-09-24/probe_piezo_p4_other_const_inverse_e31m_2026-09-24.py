# e31 = -4.1 COPY (LESSONS Sec 18.246) of probe_piezo_p4_other_const_inverse_2026-09-16.py.
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
probe_piezo_p4_other_const_inverse_2026-09-16.py

Paper 4 lever 3: Southwell G0/G1/G2 on constants OTHER than E and e31,
same F-F n=0 OC observable (oc_ff_bisect) unless noted.

  X33_C0   dielectric from clamped C0 (algebraic). Expect EXCELLENT.
  X33_freq dielectric from F-F n=0 OC frequency. Expect POOR --
           alpha ~ e31_bar^2 / Xi33_bar, frequency is a weak
           thermometer for the dielectric.
  nu       host Poisson ratio from F-F n=0 OC. Expect GOOD (elastic).
  C11E     piezo-layer raw C11^E from F-F n=0 OC. Expect MODERATE
           (added-layer stiffness).

Do not pretend a constant is identified if G2 says otherwise.
G2 is reported, not gated. G0/G1 gated for SENTINEL.

No SOLVER_VERSION bump. Do not joint-fit. Do not invert C-C OC.
"""
from __future__ import annotations

import json
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
OUTER = int(os.environ.get("PIEZO_OUTER_ITERS", "26"))
G0_TOL = 1e-6

R_I, R_O, H = 0.1, 0.6, 0.01
E_TRUE, NU_TRUE, RHO = 200e9, 0.3, 7800.0
C11E_TRUE, C12E, C13E, C33E = 132e9, 71e9, 73e9, 115e9
E31_TRUE, E33 = -4.1, 14.1
X11, X33_TRUE = 7.124e-9, 5.841e-9
RHO_PZT = 7500.0
H1 = (2.0 / 12) * H
OC_LO, OC_HI = 720.0, 800.0


def _base(**overrides):
    kw = dict(
        r_i=R_I, r_o=R_O, h=H, E=E_TRUE, nu=NU_TRUE, rho=RHO, h1=H1,
        C11E=C11E_TRUE, C12E=C12E, C13E=C13E, C33E=C33E,
        e31=E31_TRUE, e33=E33, X11=X11, X33=X33_TRUE, rho_pzt=RHO_PZT,
        dps=DPS,
    )
    kw.update(overrides)
    return PiezoOutOfPlaneSolver(**kw)


def _omega_oc(**overrides):
    return _base(**overrides).oc_ff_bisect(OC_LO, OC_HI, 0, iters=ITERS)


def _invert(fwd, target, lo, hi):
    g_lo = fwd(lo) - target
    g_hi = fwd(hi) - target
    sign_flip = (g_lo > 0) != (g_hi > 0)
    a, b, ga = lo, hi, g_lo
    for _ in range(OUTER):
        mid = 0.5 * (a + b)
        gm = fwd(mid) - target
        if (gm > 0) == (ga > 0):
            a, ga = mid, gm
        else:
            b = mid
    return 0.5 * (a + b), sign_flip


def _g2(fwd, p_true, omega_target, lo, hi):
    dp = 0.01 * p_true
    domega_dp = (fwd(p_true + dp) - fwd(p_true - dp)) / (2 * dp)
    g2 = {}
    for rel in (1e-4, 1e-3):
        dp_prop = (rel * omega_target) / abs(domega_dp) if domega_dp else float("inf")
        g2["%.4f" % (100 * rel)] = dp_prop / abs(p_true) * 100.0
    return domega_dp, g2


def _one_freq(name, p_true, lo, hi, fwd, expect):
    t0 = time.time()
    omega = fwd(p_true)
    p_hat, flip = _invert(fwd, omega, lo, hi)
    g0_rel = abs(p_hat - p_true) / abs(p_true)
    g0 = flip and g0_rel < G0_TOL
    mid = 0.5 * (p_hat + hi)
    g_mid = fwd(mid) - omega
    g_hi = fwd(hi) - omega
    g1 = flip and ((g_mid > 0) == (g_hi > 0))
    domega_dp, g2 = _g2(fwd, p_true, omega, lo, hi)
    rec = dict(
        name=name, p_true=p_true, omega=omega, p_hat=p_hat,
        g0_rel=g0_rel, g0=g0, g1=g1, expect=expect,
        domega_dp=domega_dp, g2_pct_001=g2["0.0100"],
        g2_pct_01=g2["0.1000"], elapsed_s=time.time() - t0,
    )
    print("%s: omega=%.6f p_hat=%.10g g0_rel=%.3e g0=%s g1=%s "
          "G2@0.01%%=%.4f%% expect=%s (%.1fs)"
          % (name, omega, p_hat, g0_rel, g0, g1, rec["g2_pct_001"],
             expect, rec["elapsed_s"]), flush=True)
    return rec


def main():
    t_all = time.time()
    print("job start: dps=%s iters=%s outer=%s" % (DPS, ITERS, OUTER), flush=True)
    if ps.SOLVER_VERSION != EXPECT_SOLVER_VERSION:
        print("PREFLIGHT FAIL: SOLVER_VERSION=%r != %r"
              % (ps.SOLVER_VERSION, EXPECT_SOLVER_VERSION), flush=True)
        sys.exit(2)
    if not hasattr(PiezoOutOfPlaneSolver, "oc_ff_bisect"):
        print("PREFLIGHT FAIL: no oc_ff_bisect", flush=True)
        sys.exit(2)
    print("preflight OK", flush=True)

    recs = []

    # --- X33 from C0 (algebraic) ---
    t0 = time.time()
    s = _base()
    c0 = s.clamped_C0()
    A = math_pi = __import__("math").pi * (R_O ** 2 - R_I ** 2)
    xi33b_hat = c0 * H1 / (2.0 * A)
    x33_hat = xi33b_hat - (E33 ** 2) / C33E
    g0_c0 = abs(x33_hat - X33_TRUE) / X33_TRUE < G0_TOL
    recs.append(dict(
        name="X33_C0", p_true=X33_TRUE, c0=c0, p_hat=x33_hat,
        g0_rel=abs(x33_hat - X33_TRUE) / X33_TRUE, g0=g0_c0, g1=True,
        expect="excellent", g2_pct_001=0.0, elapsed_s=time.time() - t0,
    ))
    print("X33_C0: C0=%.6e X33_hat=%.10g g0_rel=%.3e g0=%s (algebraic)"
          % (c0, x33_hat, recs[-1]["g0_rel"], g0_c0), flush=True)

    recs.append(_one_freq(
        "X33_freq", X33_TRUE, 3.0e-9, 9.0e-9,
        lambda x: _omega_oc(X33=x), "poor"))
    recs.append(_one_freq(
        "nu", NU_TRUE, 0.20, 0.40,
        lambda n: _omega_oc(nu=n), "good"))
    recs.append(_one_freq(
        "C11E", C11E_TRUE, 100e9, 160e9,
        lambda c: _omega_oc(C11E=c), "moderate"))

    g0_all = all(r["g0"] for r in recs)
    g1_all = all(r["g1"] for r in recs)
    all_pass = g0_all and g1_all
    out_dir = os.path.dirname(__file__) or "."
    results = dict(
        dps=DPS, solver_version=ps.SOLVER_VERSION,
        all_pass=all_pass, recs=recs, elapsed_s=time.time() - t_all,
    )
    path = os.path.join(out_dir, "piezo_p4_other_const_inverse_results_e31m.json")
    with open(path, "w") as f:
        json.dump(results, f, indent=2)
    print("SENTINEL " + ("PASS_ALL" if all_pass else "FAIL_ALL")
          + " g0=%s g1=%s" % (g0_all, g1_all), flush=True)


if __name__ == "__main__":
    main()
