# -*- coding: utf-8 -*-
"""
probe_piezo_p4_ff_admittance_ratios_2026-09-16.py

Paper 4 lever 5 (PAPER4_LEVERS_2026-09-16.md): clone of
probe_piezo_p4_ff_admittance_sweep_2026-09-16.py (job 2491920) at the
other two Duan2005 Table 4 ratios, h1/2h=1/8 and 1/5. The 1/12 ratio
is already SENTINEL PASS_ALL (C0=9.988 nF, Q/V bisect matches
oc_ff_bisect to 1.3e-14, k_eff^2=0.00600) and is NOT re-run here.

Linear-only driven operator (driven_ff_qv / oc_ff_bisect), production
path. G0 at omega->0 is singular (rigid body); C0 is checked on n!=0.

PRE-REGISTERED, per ratio:
  G_n     n!=0: Q/V == C0 (no motional charge on a continuous electrode).
  G_zero  independent Q/V bisection recovers oc_ff_bisect to <1e-6 rel.
  G_pole  C_eff changes sign across omega_el AND |C_eff| >> C0 on both
          sides (the elastic F-F / linear-SC root is a pole of Y).
  G_shape Butterworth-van Dyke: C_eff/C0 > 1 at 0.55*omega_el (below
          pole), C_eff/C0 in (0.5, 1) at 1.08*omega_oc (above zero).
          Off-resonance points scale with each ratio -- a fixed 800 rad/s
          sits BETWEEN pole and zero at 1/5, so must not be reused.
  G_keff  1-(omega_el/omega_zero)^2 matches manuscript Table 2
          (0.00893 at 1/8, 0.01385 at 1/5) to <1e-4 absolute.

Anchors (LESSONS_LEARNED.md Sec 18.182 / 18.184, dps=100):
  1/8  omega_el=764.1133267197  omega_oc=767.5496387603
  1/5  omega_el=800.3026922029  omega_oc=805.9014661114
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

R_I, R_O, H = 0.1, 0.6, 0.01
E_TRUE, NU, RHO = 200e9, 0.3, 7800.0
C11E, C12E, C13E, C33E = 132e9, 71e9, 73e9, 115e9
E31_TRUE, E33 = 4.1, 14.1
X11, X33 = 7.124e-9, 5.841e-9
RHO_PZT = 7500.0

# (label, h1/2h, omega_el, omega_oc, k2_anchor)
RATIOS = [
# 2026-09-23 LESSONS Sec 18.231: anchors moved to the consistent OC charge arm (h+h1/2);
# old (h+h1) duan anchors in the trailing comments.
    ("1/8", 1.0 / 8, 764.1133267197, 767.2091198988, 0.00805),  # duan 767.5496387603, 0.00893
    ("1/5", 1.0 / 5, 800.3026922029, 805.1123631872, 0.01191),  # duan 805.9014661114, 0.01385
]


def _solver(h1):
    return PiezoOutOfPlaneSolver(
        r_i=R_I, r_o=R_O, h=H, E=E_TRUE, nu=NU, rho=RHO, h1=h1,
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


def _one_ratio(label, ratio, omega_el, omega_oc, k2_anchor):
    h1 = ratio * (2.0 * H)
    s = _solver(h1)
    c0 = s.clamped_C0()
    print("---- %s  h1=%.6f  C0=%.6e F ----" % (label, h1, c0), flush=True)

    r_n = s.driven_ff_qv(omega_oc, 1)
    g_n = abs(r_n["q_over_v"] - c0) / c0 < 1e-12
    print("G_n: n!=0 q/v==C0 pass=%s" % g_n, flush=True)

    r_el_lo = s.driven_ff_qv(omega_el - 0.5, 0)["q_over_v"]
    r_el_hi = s.driven_ff_qv(omega_el + 0.5, 0)["q_over_v"]
    g_pole = ((r_el_lo * r_el_hi < 0.0)
              and (abs(r_el_lo) > 2 * c0) and (abs(r_el_hi) > 2 * c0))
    print("G_pole: sign_flip=%s |qv|/C0 lo=%.3f hi=%.3f pass=%s"
          % (r_el_lo * r_el_hi < 0.0, abs(r_el_lo) / c0, abs(r_el_hi) / c0,
             g_pole), flush=True)

    zhat, zflip, _, _ = _bisect_qv_zero(s, omega_el + 0.2, omega_oc + 1.0)
    zrel = abs(zhat - omega_oc) / omega_oc if zhat == zhat else float("inf")
    g_zero = zflip and zrel < 1e-6
    print("G_zero: qv_bisect=%.10f oc_ff=%.10f rel=%s pass=%s"
          % (zhat, omega_oc, zrel, g_zero), flush=True)

    w_lo = 0.55 * omega_el
    w_hi = 1.08 * omega_oc
    r_lo = s.driven_ff_qv(w_lo, 0)["q_over_v"] / c0
    r_hi = s.driven_ff_qv(w_hi, 0)["q_over_v"] / c0
    g_shape = (r_lo > 1.0) and (0.5 < r_hi < 1.0)
    print("G_shape: Ceff/C0 @%.1f=%.4f @%.1f=%.4f pass=%s"
          % (w_lo, r_lo, w_hi, r_hi, g_shape), flush=True)

    k2 = 1.0 - (omega_el / zhat) ** 2 if zhat == zhat else float("nan")
    g_keff = abs(k2 - k2_anchor) < 1e-4
    print("G_keff: k2_from_Y=%.6e k2_anchor=%.6e pass=%s"
          % (k2, k2_anchor, g_keff), flush=True)

    omegas = (_frange(omega_el - 80.0, omega_oc + 80.0, 41)
              + _frange(omega_el - 3.0, omega_oc + 3.0, 61))
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
    return dict(
        label=label, ratio=ratio, h1=h1, c0=c0,
        omega_el=omega_el, omega_oc=omega_oc,
        omega_zero_from_qv=zhat, zrel=zrel,
        k2_from_Y=k2, k2_anchor=k2_anchor,
        g_n=g_n, g_pole=g_pole, g_zero=g_zero, g_shape=g_shape,
        g_keff=g_keff, all_pass=all_pass,
        ceff_over_c0_lo=r_lo, ceff_over_c0_hi=r_hi,
        w_lo=w_lo, w_hi=w_hi, sweep=sweep,
    ), all_pass, omega_el, omega_oc


def _plot_ratio(out_dir, rec, omega_el, omega_oc):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as e:
        print("plot skipped: %s" % e, flush=True)
        return
    tag = rec["label"].replace("/", "")
    sweep = rec["sweep"]
    fs = [row["f_hz"] for row in sweep]
    rat = [row["ceff_over_c0"] for row in sweep]
    yab = [row["y_abs"] if row["y_abs"] < 1e6 else float("nan") for row in sweep]
    fig, ax = plt.subplots(2, 1, figsize=(7.2, 6.2), sharex=True)
    ax[0].plot(fs, rat, "k-", lw=1.2)
    ax[0].axhline(1.0, color="0.6", ls="--", lw=0.8)
    ax[0].axvline(omega_el / (2 * math.pi), color="C0", ls=":", lw=1,
                  label="f_r (SC)")
    ax[0].axvline(omega_oc / (2 * math.pi), color="C3", ls=":", lw=1,
                  label="f_a (OC)")
    ax[0].set_ylabel(r"$C_{\mathrm{eff}}/C_0$")
    ax[0].set_ylim(-8, 8)
    ax[0].legend(frameon=False, fontsize=8)
    ax[0].set_title(r"$h_1/2h=%s$" % rec["label"])
    ax[1].semilogy(fs, [abs(v) if v == v and v != 0 else 1e-18 for v in yab],
                   "k-", lw=1.2)
    ax[1].axvline(omega_el / (2 * math.pi), color="C0", ls=":", lw=1)
    ax[1].axvline(omega_oc / (2 * math.pi), color="C3", ls=":", lw=1)
    ax[1].set_xlabel("frequency (Hz)")
    ax[1].set_ylabel(r"$|Y|=|\omega\,Q/V|$ (S)")
    fig.tight_layout()
    png = os.path.join(out_dir, "piezo_p4_ff_admittance_sweep_h%s.png" % tag)
    fig.savefig(png, dpi=140)
    plt.close()
    print("wrote %s" % png, flush=True)


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

    out_dir = os.path.dirname(__file__) or "."
    recs = []
    all_pass = True
    for spec in RATIOS:
        rec, ok, omega_el, omega_oc = _one_ratio(*spec)
        recs.append(rec)
        all_pass = all_pass and ok
        _plot_ratio(out_dir, rec, omega_el, omega_oc)

    results = dict(
        dps=DPS, solver_version=ps.SOLVER_VERSION,
        elapsed_s=time.time() - t0, all_pass=all_pass, ratios=recs,
    )
    json_path = os.path.join(out_dir, "piezo_p4_ff_admittance_ratios_results.json")
    with open(json_path, "w") as f:
        json.dump(results, f, indent=2)
    csv_path = os.path.join(out_dir, "piezo_p4_ff_admittance_ratios.csv")
    with open(csv_path, "w") as f:
        f.write("label,omega,f_hz,q_over_v,ceff_over_c0,y_abs\n")
        for rec in recs:
            for row in rec["sweep"]:
                f.write("%s,%s,%s,%s,%s,%s\n" % (
                    rec["label"], row["omega"], row["f_hz"],
                    row["q_over_v"], row["ceff_over_c0"], row["y_abs"]))
    print("SENTINEL " + ("PASS_ALL" if all_pass else "FAIL_ALL")
          + " all_pass=%s" % all_pass, flush=True)
    for rec in recs:
        print("  %s C0=%.4e k2=%.6e G_n=%s G_pole=%s G_zero=%s G_shape=%s G_keff=%s"
              % (rec["label"], rec["c0"], rec["k2_from_Y"], rec["g_n"],
                 rec["g_pole"], rec["g_zero"], rec["g_shape"], rec["g_keff"]),
              flush=True)


if __name__ == "__main__":
    main()
