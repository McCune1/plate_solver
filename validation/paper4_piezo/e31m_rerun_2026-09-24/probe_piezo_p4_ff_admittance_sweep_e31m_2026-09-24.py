# e31 = -4.1 COPY (LESSONS Sec 18.246) of probe_piezo_p4_ff_admittance_sweep_2026-09-16.py.
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
probe_piezo_p4_ff_admittance_sweep_2026-09-16.py

Paper 4, remaining piece of item 3 (impedance / antiresonance / k_eff):
the driven F-F n=0 admittance Y(omega) = j*omega*Q/V on the same
linear-potential model as oc_ff_det (LESSONS_LEARNED.md Sec 18.175/18.179).

Imposed V on the outer-electrode bus, F-F mechanical BC, inners grounded.
Q from the Sec 18.175 charge formula. Undamped, so Y is purely imaginary
and C_eff(omega) = Q/V is real.

PRE-REGISTERED:
  G_n     n!=0: Q/V == C0 (no motional charge on a continuous electrode).
  G_zero  C_eff crosses 0 at omega_OC (independent bisection of Q/V
          recovers oc_ff_bisect to <1e-6 relative).
  G_pole  C_eff changes sign across omega_el AND |C_eff| >> C0 on both
          sides (the SC / elastic F-F root is a pole of Y).
  G_shape Butterworth-van Dyke: C_eff/C0 > 1 below the pole, C_eff/C0
          in (0, 1) above the zero (approaches C0 from below).
  G_keff  1 - (omega_el/omega_zero)^2 matches the Sec 18.175 k_eff^2
          = 0.00600 to <1e-4 absolute (same two frequencies, now
          identified as pole and zero of the driven problem).

A coarse 600-900 rad/s sweep plus a dense window around the first
n=0 pair is written to JSON/CSV. Optional matplotlib PNG if available.

Cost: ~200 4x4 lu_solves at dps=60, sandbox-seconds, not a 6x6 search.
"""
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

R_I, R_O, H = 0.1, 0.6, 0.01
E_TRUE, NU, RHO = 200e9, 0.3, 7800.0
C11E, C12E, C13E, C33E = 132e9, 71e9, 73e9, 115e9
E31_TRUE, E33 = -4.1, 14.1
X11, X33 = 7.124e-9, 5.841e-9
RHO_PZT = 7500.0
H1 = (2.0 / 12) * H
ELASTIC_N0 = 747.9863950952586
# 2026-09-23 LESSONS Sec 18.231: anchors moved to the consistent OC charge arm (h+h1/2);
# old (h+h1) duan anchors in the trailing comments.
# e31 = -4.1 anchors from probe_piezo_p4_h1_ratio_sweep_e31m_2026-09-24 (1/12 row);
# +4.1 anchors were 750.0794987774267 and 0.005573231299129566.
OC_N0 = 762.5693785053465
K2_ANCHOR = 0.03788126062047603


def _solver():
    return PiezoOutOfPlaneSolver(
        r_i=R_I, r_o=R_O, h=H, E=E_TRUE, nu=NU, rho=RHO, h1=H1,
        C11E=C11E, C12E=C12E, C13E=C13E, C33E=C33E,
        e31=E31_TRUE, e33=E33, X11=X11, X33=X33, rho_pzt=RHO_PZT, dps=DPS)


def _bisect_qv_zero(s, lo, hi, iters=40):
    glo = s.driven_ff_qv(lo, 0)["q_over_v"]
    ghi = s.driven_ff_qv(hi, 0)["q_over_v"]
    if (glo > 0) == (ghi > 0):
        return float("nan"), False, glo, ghi
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        gm = s.driven_ff_qv(mid, 0)["q_over_v"]
        if (gm > 0) == (glo > 0):
            lo, glo = mid, gm
        else:
            hi = mid
    return 0.5 * (lo + hi), True, glo, ghi


def _frange(lo, hi, n):
    if n <= 1:
        return [lo]
    step = (hi - lo) / (n - 1)
    return [lo + i * step for i in range(n)]


def main():
    t0 = time.time()
    print("job start: dps=%s" % DPS, flush=True)
    if ps.SOLVER_VERSION != EXPECT_SOLVER_VERSION:
        print("PREFLIGHT FAIL: SOLVER_VERSION=%r != %r"
              % (ps.SOLVER_VERSION, EXPECT_SOLVER_VERSION), flush=True)
        sys.exit(2)
    if not hasattr(PiezoOutOfPlaneSolver, "driven_ff_qv"):
        print("PREFLIGHT FAIL: no driven_ff_qv -- push piezo_solver.py",
              flush=True)
        sys.exit(2)
    print("preflight OK: SOLVER_VERSION=%s driven_ff_qv present"
          % ps.SOLVER_VERSION, flush=True)

    s = _solver()
    c0 = s.clamped_C0()
    print("C0=%.6e F" % c0, flush=True)

    r_n = s.driven_ff_qv(750.0, 1)
    g_n = abs(r_n["q_over_v"] - c0) / c0 < 1e-12
    print("G_n: n!=0 q/v==C0 pass=%s" % g_n, flush=True)

    r_el_lo = s.driven_ff_qv(ELASTIC_N0 - 0.5, 0)["q_over_v"]
    r_el_hi = s.driven_ff_qv(ELASTIC_N0 + 0.5, 0)["q_over_v"]
    g_pole = (r_el_lo * r_el_hi < 0.0) and (abs(r_el_lo) > 2 * c0) and (abs(r_el_hi) > 2 * c0)
    print("G_pole: sign_flip=%s |qv|/C0 lo=%.3f hi=%.3f pass=%s"
          % (r_el_lo * r_el_hi < 0.0, abs(r_el_lo) / c0, abs(r_el_hi) / c0, g_pole),
          flush=True)

    zhat, zflip, _, _ = _bisect_qv_zero(s, ELASTIC_N0 + 0.2, OC_N0 + 1.0)
    zrel = abs(zhat - OC_N0) / OC_N0 if zhat == zhat else float("inf")
    g_zero = zflip and zrel < 1e-6
    print("G_zero: qv_bisect=%.10f oc_ff=%.10f rel=%s pass=%s"
          % (zhat, OC_N0, zrel, g_zero), flush=True)

    r400 = s.driven_ff_qv(400.0, 0)["q_over_v"] / c0
    r800 = s.driven_ff_qv(800.0, 0)["q_over_v"] / c0
    g_shape = (r400 > 1.0) and (0.5 < r800 < 1.0)
    print("G_shape: Ceff/C0 @400=%.4f @800=%.4f pass=%s"
          % (r400, r800, g_shape), flush=True)

    k2 = 1.0 - (ELASTIC_N0 / zhat) ** 2 if zhat == zhat else float("nan")
    g_keff = abs(k2 - K2_ANCHOR) < 1e-4
    print("G_keff: k2_from_Y=%.6e k2_anchor=%.6e pass=%s" % (k2, K2_ANCHOR, g_keff),
          flush=True)

    omegas = _frange(600.0, 900.0, 61) + _frange(745.0, 753.0, 81)
    omegas = sorted(set(round(w, 8) for w in omegas))
    sweep = []
    for w in omegas:
        r = s.driven_ff_qv(w, 0)
        qv = r["q_over_v"]
        yabs = abs(w * qv) if qv == qv and abs(qv) < 1e30 else float("inf")
        sweep.append({
            "omega": w,
            "f_hz": w / (2.0 * math.pi),
            "q_over_v": qv,
            "ceff_over_c0": qv / c0 if c0 else float("nan"),
            "y_abs": yabs,
        })
    print("sweep %d points" % len(sweep), flush=True)

    all_pass = g_n and g_pole and g_zero and g_shape and g_keff
    out_dir = os.path.dirname(__file__) or "."
    results = {
        "c0": c0,
        "omega_el": ELASTIC_N0,
        "omega_oc": OC_N0,
        "omega_zero_from_qv": zhat,
        "k2_from_Y": k2,
        "k2_anchor": K2_ANCHOR,
        "g_n": g_n,
        "g_pole": g_pole,
        "g_zero": g_zero,
        "g_shape": g_shape,
        "g_keff": g_keff,
        "ceff_over_c0_400": r400,
        "ceff_over_c0_800": r800,
        "elapsed_s": time.time() - t0,
        "sweep": sweep,
    }
    json_path = os.path.join(out_dir, "piezo_p4_ff_admittance_sweep_results_e31m.json")
    with open(json_path, "w") as f:
        json.dump(results, f, indent=2)
    csv_path = os.path.join(out_dir, "piezo_p4_ff_admittance_sweep_e31m.csv")
    with open(csv_path, "w") as f:
        f.write("omega,f_hz,q_over_v,ceff_over_c0,y_abs\n")
        for row in sweep:
            f.write("%s,%s,%s,%s,%s\n" % (
                row["omega"], row["f_hz"], row["q_over_v"],
                row["ceff_over_c0"], row["y_abs"]))
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fs = [row["f_hz"] for row in sweep]
        rat = [row["ceff_over_c0"] for row in sweep]
        yab = [row["y_abs"] if row["y_abs"] < 1e6 else float("nan") for row in sweep]
        fig, ax = plt.subplots(2, 1, figsize=(7.2, 6.2), sharex=True)
        ax[0].plot(fs, rat, "k-", lw=1.2)
        ax[0].axhline(1.0, color="0.6", ls="--", lw=0.8)
        ax[0].axvline(ELASTIC_N0 / (2 * math.pi), color="C0", ls=":", lw=1,
                      label="f_r (SC)")
        ax[0].axvline(OC_N0 / (2 * math.pi), color="C3", ls=":", lw=1,
                      label="f_a (OC)")
        ax[0].set_ylabel(r"$C_{\mathrm{eff}}/C_0$")
        ax[0].set_ylim(-8, 8)
        ax[0].legend(frameon=False, fontsize=8)
        ax[1].semilogy(fs, [abs(v) if v == v and v != 0 else 1e-18 for v in yab],
                       "k-", lw=1.2)
        ax[1].axvline(ELASTIC_N0 / (2 * math.pi), color="C0", ls=":", lw=1)
        ax[1].axvline(OC_N0 / (2 * math.pi), color="C3", ls=":", lw=1)
        ax[1].set_xlabel("frequency (Hz)")
        ax[1].set_ylabel(r"$|Y|=|\omega\,Q/V|$ (S)")
        fig.tight_layout()
        png = os.path.join(out_dir, "piezo_p4_ff_admittance_sweep_e31m.png")
        fig.savefig(png, dpi=140)
        plt.close()
        print("wrote %s" % png, flush=True)
    except Exception as e:
        print("plot skipped: %s" % e, flush=True)

    print("SENTINEL " + ("PASS_ALL" if all_pass else "FAIL_ALL")
          + " all_pass=%s k2=%.6e" % (all_pass, k2), flush=True)


if __name__ == "__main__":
    main()
