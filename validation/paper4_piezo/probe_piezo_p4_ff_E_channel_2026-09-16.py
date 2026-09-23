# -*- coding: utf-8 -*-
"""
probe_piezo_p4_ff_E_channel_2026-09-16.py

Paper 4 lever 4: host E from F-F as a second channel, vs job 2491531's
C-C SC result (0.0246% of E at 0.01% freq).

  E_el   from F-F n=0 elastic bilayer (elastic_bisect).
  E_oc   from F-F n=0 open-circuit (oc_ff_bisect). The OC shift is
         small, so this is almost the bilayer-elastic inverse.

PRE-REGISTERED:
  G0  round-trip <1e-6 relative on exact data, both channels.
  G1  sign flip + hi-side negative control, both channels.
  G2  reported: expect ~0.02% of E at 0.01% freq (omega ~ sqrt(E)),
      same order as the C-C SC 0.0246%. A second channel that
      disagrees by orders of magnitude is a FAIL of the channel,
      not a new physics claim.

No SOLVER_VERSION bump. Do not invert C-C OC. Do not joint-fit.
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
CC_G2_ANCHOR = 0.0246  # percent of E at 0.01% freq, job 2491531

R_I, R_O, H = 0.1, 0.6, 0.01
E_TRUE, NU, RHO = 200e9, 0.3, 7800.0
C11E, C12E, C13E, C33E = 132e9, 71e9, 73e9, 115e9
E31_TRUE, E33 = 4.1, 14.1
X11, X33 = 7.124e-9, 5.841e-9
RHO_PZT = 7500.0
H1 = (2.0 / 12) * H
E_LO, E_HI = 140e9, 260e9


def _solver(E):
    return PiezoOutOfPlaneSolver(
        r_i=R_I, r_o=R_O, h=H, E=E, nu=NU, rho=RHO, h1=H1,
        C11E=C11E, C12E=C12E, C13E=C13E, C33E=C33E,
        e31=E31_TRUE, e33=E33, X11=X11, X33=X33, rho_pzt=RHO_PZT, dps=DPS,
    )


def _window(E):
    center = 747.9863950953 * (E / E_TRUE) ** 0.5
    return 0.80 * center, 1.20 * center


def _omega_el(E):
    lo, hi = _window(E)
    return _solver(E).elastic_bisect(lo, hi, 0, iters=ITERS)


def _omega_oc(E):
    lo, hi = _window(E)
    return _solver(E).oc_ff_bisect(lo, hi, 0, iters=ITERS)


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


def _channel(name, fwd):
    t0 = time.time()
    omega = fwd(E_TRUE)
    e_hat, flip = _invert(fwd, omega, E_LO, E_HI)
    g0_rel = abs(e_hat - E_TRUE) / E_TRUE
    g0 = flip and g0_rel < G0_TOL
    mid = 0.5 * (e_hat + E_HI)
    g_mid = fwd(mid) - omega
    g_hi = fwd(E_HI) - omega
    g1 = flip and ((g_mid > 0) == (g_hi > 0))
    dE = 0.01 * E_TRUE
    domega_dE = (fwd(E_TRUE + dE) - fwd(E_TRUE - dE)) / (2 * dE)
    g2_001 = (1e-4 * omega / abs(domega_dE)) / E_TRUE * 100.0
    g2_01 = (1e-3 * omega / abs(domega_dE)) / E_TRUE * 100.0
    rec = dict(
        name=name, omega=omega, e_hat=e_hat, g0_rel=g0_rel,
        g0=g0, g1=g1, domega_dE=domega_dE,
        g2_pct_001=g2_001, g2_pct_01=g2_01,
        ratio_to_cc=g2_001 / CC_G2_ANCHOR,
        elapsed_s=time.time() - t0,
    )
    print("%s: omega=%.6f E_hat=%.10g g0_rel=%.3e g0=%s g1=%s "
          "G2@0.01%%=%.4f%% vs C-C 0.0246%% (ratio=%.3f) (%.1fs)"
          % (name, omega, e_hat, g0_rel, g0, g1, g2_001,
             rec["ratio_to_cc"], rec["elapsed_s"]), flush=True)
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

    recs = [_channel("E_el", _omega_el), _channel("E_oc", _omega_oc)]
    g0_all = all(r["g0"] for r in recs)
    g1_all = all(r["g1"] for r in recs)
    # G2 same order as C-C (within 5x) -- a consistency check, gated
    g2_order = all(0.2 < r["ratio_to_cc"] < 5.0 for r in recs)
    all_pass = g0_all and g1_all and g2_order
    out_dir = os.path.dirname(__file__) or "."
    results = dict(
        dps=DPS, solver_version=ps.SOLVER_VERSION,
        cc_g2_anchor_pct=CC_G2_ANCHOR, all_pass=all_pass,
        recs=recs, elapsed_s=time.time() - t_all,
    )
    path = os.path.join(out_dir, "piezo_p4_ff_E_channel_results.json")
    with open(path, "w") as f:
        json.dump(results, f, indent=2)
    print("SENTINEL " + ("PASS_ALL" if all_pass else "FAIL_ALL")
          + " g0=%s g1=%s g2_order=%s" % (g0_all, g1_all, g2_order),
          flush=True)


if __name__ == "__main__":
    main()
