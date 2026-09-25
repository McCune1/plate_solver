# e31 = -4.1 COPY (LESSONS Sec 18.246) of probe_piezo_p4_h1_ratio_sweep_2026-09-16.py.
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
probe_piezo_p4_h1_ratio_sweep_2026-09-16.py

Paper 4 lever 2 (PAPER4_LEVERS_2026-09-16.md): denser h1/2h sweep of the
already-gated F-F n=0 linear-potential open-circuit operator
(oc_ff_bisect / elastic_bisect). Turns Table 2 / Fig. 2 from three
Duan2005 Table-4 ratios into a curve, and characterizes the thin-layer
asymptote as h1/2h -> 0.

Does NOT invert sequentially (lever 15 is closed). Does NOT use the
sine-enriched OC operator (lever 12 keeps linear-only as production).
Does NOT bump SOLVER_VERSION.

PRE-REGISTERED:
  G_anchor_split  at 1/12, 1/8, 1/5 the OC/elastic split recovers
                  +0.301%, +0.450%, +0.700% to <0.005 percentage points
                  (dps=40 vs the dps=100 anchors in LESSONS Sec 18.175/
                  18.182; 4x4 Bessel, this slack is numerical not physics).
  G_anchor_k2     k_eff^2 = 1-(w_el/w_oc)^2 recovers 0.00600, 0.00893,
                  0.01385 to <1e-4 absolute (manuscript Table 2).
  G_anchor_g2     e31 G2 at 0.01% freq recovers 1.98%, 1.33%, 0.86% of
                  e31_true to <0.20 percentage points.
  G_mono          split, k_eff^2 strictly increasing in h1/2h; G2
                  strictly decreasing (thicker layer, better e31).
  G_thin          as h1/2h -> 0, split -> 0 and split/(h1/2h) approaches
                  a finite intercept (reported, not gated to a number).
                  Rank-1 alpha ~ h1, so leading thin-layer scaling is
                  linear in h1/2h, not quadratic.

G0 at omega->0 is singular (rigid body) -- not computed here.
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

R_I, R_O, H = 0.1, 0.6, 0.01
E_TRUE, NU, RHO = 200e9, 0.3, 7800.0
C11E, C12E, C13E, C33E = 132e9, 71e9, 73e9, 115e9
E31_TRUE, E33 = -4.1, 14.1
X11, X33 = 7.124e-9, 5.841e-9
RHO_PZT = 7500.0

# (label, h1/2h, expected_split_pct or None, expected_k2 or None, expected_g2_pct or None)
# e31=-4.1 re-anchor (2026-09-24): anchors are cross-probe values from the
# independent -4.1 open_circuit and admittance (Y-zero) probes; +4.1 anchors were
# 0.280/0.00557/2.13, 0.405/0.00805/1.48, 0.601/0.01191/1.00 (LESSONS Sec 18.175).
SWEEP = [
    ("1/40", 1.0 / 40, None, None, None),
    ("1/30", 1.0 / 30, None, None, None),
    ("1/24", 1.0 / 24, None, None, None),
    ("1/20", 1.0 / 20, None, None, None),
    ("1/16", 1.0 / 16, None, None, None),
    ("1/12", 1.0 / 12, 1.9496, 0.03788, 0.870),  # e31=-4.1: open_circuit_e31m (drel, g3_k2, G2)
    ("1/10", 1.0 / 10, None, None, None),
    ("1/8",  1.0 / 8,  2.7763, 0.05330, None),  # e31=-4.1: admittance_ratios_e31m Y-zero route
    ("1/6",  1.0 / 6,  None, None, None),
    ("1/5",  1.0 / 5,  4.0155, 0.07572, None),  # e31=-4.1: admittance_ratios_e31m Y-zero route
    ("1/4",  1.0 / 4,  None, None, None),
]


def _kwargs(h1, e31=E31_TRUE):
    return dict(
        r_i=R_I, r_o=R_O, h=H, E=E_TRUE, nu=NU, rho=RHO, h1=h1,
        C11E=C11E, C12E=C12E, C13E=C13E, C33E=C33E,
        e31=e31, e33=E33, X11=X11, X33=X33, rho_pzt=RHO_PZT, dps=DPS,
    )


def _bisect_pair(h1, lo, hi, e31=E31_TRUE):
    s = PiezoOutOfPlaneSolver(**_kwargs(h1, e31))
    w_el = s.elastic_bisect(lo, hi, 0, iters=ITERS)
    w_oc = s.oc_ff_bisect(lo, hi, 0, iters=ITERS)
    return w_el, w_oc


def main():
    t_all = time.time()
    print("job start: dps=%s iters=%s" % (DPS, ITERS), flush=True)
    if ps.SOLVER_VERSION != EXPECT_SOLVER_VERSION:
        print("PREFLIGHT FAIL: SOLVER_VERSION=%r != %r"
              % (ps.SOLVER_VERSION, EXPECT_SOLVER_VERSION), flush=True)
        sys.exit(2)
    if not hasattr(PiezoOutOfPlaneSolver, "oc_ff_bisect"):
        print("PREFLIGHT FAIL: no oc_ff_bisect -- push piezo_solver.py",
              flush=True)
        sys.exit(2)
    print("preflight OK: SOLVER_VERSION=%s oc_ff_bisect present"
          % ps.SOLVER_VERSION, flush=True)

    rows = []
    lo, hi = 700.0, 780.0
    g_anchor_split = True
    g_anchor_k2 = True
    g_anchor_g2 = True
    for label, ratio, exp_split, exp_k2, exp_g2 in SWEEP:
        t0 = time.time()
        h1 = ratio * (2.0 * H)
        w_el, w_oc = _bisect_pair(h1, lo, hi)
        if not (w_oc > w_el > 0):
            print("FATAL: no stiffening at %s (el=%.6f oc=%.6f lo=%.1f hi=%.1f)"
                  % (label, w_el, w_oc, lo, hi), flush=True)
            sys.exit(3)
        split = (w_oc - w_el) / w_el
        k2 = 1.0 - (w_el / w_oc) ** 2
        de = 0.01 * E31_TRUE
        _, w_plus = _bisect_pair(h1, lo, hi, e31=E31_TRUE + de)
        _, w_minus = _bisect_pair(h1, lo, hi, e31=E31_TRUE - de)
        dw_de = (w_plus - w_minus) / (2.0 * de)
        g2_001 = (1e-4 * w_oc / abs(dw_de)) / abs(E31_TRUE) * 100.0
        g2_01 = (1e-3 * w_oc / abs(dw_de)) / abs(E31_TRUE) * 100.0
        split_over_ratio = split / ratio
        row = dict(
            label=label, ratio=ratio, h1=h1,
            omega_el=w_el, omega_oc=w_oc,
            f_el_hz=w_el / (2.0 * math.pi),
            f_oc_hz=w_oc / (2.0 * math.pi),
            split=split, split_pct=100.0 * split,
            k2=k2, split_over_ratio=split_over_ratio,
            domega_de31=dw_de, g2_pct_001=g2_001, g2_pct_01=g2_01,
            elapsed_s=time.time() - t0,
        )
        rows.append(row)
        a_split = a_k2 = a_g2 = ""
        if exp_split is not None:
            ok_s = abs(100.0 * split - exp_split) < 0.005
            g_anchor_split = g_anchor_split and ok_s
            a_split = " anchor_split=%s (exp %.3f)" % (ok_s, exp_split)
        if exp_k2 is not None:
            ok_k = abs(k2 - exp_k2) < 1e-4
            g_anchor_k2 = g_anchor_k2 and ok_k
            a_k2 = " anchor_k2=%s" % ok_k
        if exp_g2 is not None:
            ok_g = abs(g2_001 - exp_g2) < 0.20
            g_anchor_g2 = g_anchor_g2 and ok_g
            a_g2 = " anchor_g2=%s (exp %.2f)" % (ok_g, exp_g2)
        print(
            "ratio=%-5s h1=%.6f  el=%.6f oc=%.6f  split=%.4f%%  k2=%.6f  "
            "G2@0.01%%=%.3f%%  split/ratio=%.5f  (%.1fs)%s%s%s"
            % (label, h1, w_el, w_oc, 100.0 * split, k2, g2_001,
               split_over_ratio, row["elapsed_s"], a_split, a_k2, a_g2),
            flush=True,
        )
        lo = w_el * 0.97
        hi = w_oc * 1.10

    splits = [r["split"] for r in rows]
    k2s = [r["k2"] for r in rows]
    g2s = [r["g2_pct_001"] for r in rows]
    g_mono = (
        all(splits[i] < splits[i + 1] for i in range(len(splits) - 1))
        and all(k2s[i] < k2s[i + 1] for i in range(len(k2s) - 1))
        and all(g2s[i] > g2s[i + 1] for i in range(len(g2s) - 1))
    )
    thin = rows[:3]
    thin_splits = [r["split"] for r in thin]
    g_thin = (thin_splits[0] < thin_splits[1] < thin_splits[2]
              and thin[0]["split"] < 0.002 * 7.24)  # e31=-4.1: +4.1 bound x (13.05/4.85)^2
    intercepts = [r["split_over_ratio"] for r in thin]

    all_pass = (g_anchor_split and g_anchor_k2 and g_anchor_g2
                and g_mono and g_thin)
    out_dir = os.path.dirname(__file__) or "."
    results = dict(
        dps=DPS, iters=ITERS, solver_version=ps.SOLVER_VERSION,
        g_anchor_split=g_anchor_split, g_anchor_k2=g_anchor_k2,
        g_anchor_g2=g_anchor_g2, g_mono=g_mono, g_thin=g_thin,
        thin_split_over_ratio=intercepts,
        all_pass=all_pass, elapsed_s=time.time() - t_all, rows=rows,
    )
    json_path = os.path.join(out_dir, "piezo_p4_h1_ratio_sweep_results_e31m.json")
    with open(json_path, "w") as f:
        json.dump(results, f, indent=2)
    csv_path = os.path.join(out_dir, "piezo_p4_h1_ratio_sweep_e31m.csv")
    with open(csv_path, "w") as f:
        f.write("label,ratio,h1,omega_el,omega_oc,f_el_hz,f_oc_hz,"
                "split_pct,k2,g2_pct_001,split_over_ratio\n")
        for r in rows:
            f.write("%s,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g\n"
                    % (r["label"], r["ratio"], r["h1"], r["omega_el"],
                       r["omega_oc"], r["f_el_hz"], r["f_oc_hz"],
                       r["split_pct"], r["k2"], r["g2_pct_001"],
                       r["split_over_ratio"]))
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        xs = [r["ratio"] for r in rows]
        fig, ax = plt.subplots(3, 1, figsize=(7.2, 8.4), sharex=True)
        ax[0].plot(xs, [r["split_pct"] for r in rows], "ko-", lw=1.2, ms=4)
        for r in rows:
            if r["label"] in ("1/12", "1/8", "1/5"):
                ax[0].plot(r["ratio"], r["split_pct"], "C3o", ms=7)
        ax[0].set_ylabel(r"OC split (\%)")
        ax[1].plot(xs, [r["k2"] for r in rows], "ko-", lw=1.2, ms=4)
        ax[1].set_ylabel(r"$k_{\mathrm{eff}}^2$")
        ax[2].plot(xs, [r["g2_pct_001"] for r in rows], "ko-", lw=1.2, ms=4)
        ax[2].set_ylabel(r"$e_{31}$ G2 at 0.01\% freq (\%)")
        ax[2].set_xlabel(r"$h_1/2h$")
        for a in ax:
            a.axvline(1.0 / 12, color="0.7", ls=":", lw=0.8)
            a.axvline(1.0 / 8, color="0.7", ls=":", lw=0.8)
            a.axvline(1.0 / 5, color="0.7", ls=":", lw=0.8)
        fig.tight_layout()
        png = os.path.join(out_dir, "piezo_p4_h1_ratio_sweep_e31m.png")
        fig.savefig(png, dpi=140)
        plt.close()
        print("wrote %s" % png, flush=True)
    except Exception as e:
        print("plot skipped: %s" % e, flush=True)

    print("thin-layer split/(h1/2h) at 1/40,1/30,1/24 = %s" % intercepts,
          flush=True)
    print("SENTINEL " + ("PASS_ALL" if all_pass else "FAIL_ALL")
          + " split=%s k2=%s g2=%s mono=%s thin=%s"
          % (g_anchor_split, g_anchor_k2, g_anchor_g2, g_mono, g_thin),
          flush=True)


if __name__ == "__main__":
    main()
