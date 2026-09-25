# e31 = -4.1 COPY (LESSONS Sec 18.246) of probe_piezo_p5_thickness_sweep_2026-09-20.py.
# Written by migrate_e31m_p5.py. Only e31 changes, plus the windows that
# bracket the F-F SC root (471.36 -> 519.84 rad/s at -4.1, i.e. every
# window near that root shifts by +48.47), the e31 inverse bracket
# [1,8] -> [-8,-1] (wrong-sign control [-8,-1] -> [1,8]; the twin
# 2*8.95 - e31 moves to 22.0), relative errors divide by |E31_TRUE|,
# and output files carry the suffix _e31m. Fixed before running.
# -*- coding: utf-8 -*-
"""
probe_piezo_p5_thickness_sweep_2026-09-20.py

Paper 5 roadmap Sec 8 item 6, "optional cheap: a thickness sweep (the
homogeneous analogue of Paper 4's h1/2h curve is now H/r_o)". Pure
probe-level scan on the existing PiezoMonolithicOutOfPlaneSolver, no
new solver code.

WHY THIS IS A REAL QUESTION, NOT JUST A CURVE FOR ITS OWN SAKE:
PAPER5_DERIVATION.md Sec 6's long-wave reduction gives
    D_SC/d - 1 = 12*e31_bar^2 / (pi^2 * Xi33_bar * c11_bar),
and H cancels out of that ratio entirely (d ~ H^3 and the coupling
term ~ H^3 the same way) -- so to LEADING long-wave order, the SC-vs-
elastic percentage split should be independent of thickness H at fixed
r_i, r_o. That is a genuinely falsifiable prediction of the closed-form
estimate, and the full coupled_bisect (chi-cubic + Bessel, not just
the long-wave truncation) is the right tool to check it against, the
same way Sec 18.204/18.205 checked the long-wave estimate's VALUE at
one H. This probe checks whether that H-independence actually holds
as H grows away from the already-validated H/r_o=1/60 point, or
whether finite-thickness (higher-order, non-long-wave) corrections
appear -- and if so, how large and in which direction.

F-F n=0 only (roadmap item 6 does not specify a BC; C-C is a cheap
future extension if wanted, not done here to keep this probe scoped).

Sweep points (H/r_o), found interactively this session via a linear-
scaling bracket guess (omega_fundamental scales close to linearly in H
at fixed r_i/r_o -- verified, not assumed, before finalizing this
list): 0.005, 0.01, 1/60 (=0.016667, the already-validated Sec 18.201
point), 0.025, 0.04, 0.06, 0.09, 0.13, 0.18. This spans a factor of 36
in thickness, from thinner than the validated point to H/r_o=0.18 --
well outside where Kirchhoff kinematics is usually considered
trustworthy, included deliberately to see where (not whether) the
long-wave estimate visibly departs, not because the thick end is
proposed as a physical design point.

Root-finding per H: rather than hardcoding a bracket per H (fragile if
this sweep is ever rerun at different geometry/material), each point
scans a window derived from the fundamental's own near-linear H-
scaling (guess = OMEGA_ANCHOR_H001 * H / 0.01, window
[max(20, 0.1*guess), 3.5*guess], 250 steps) for BOTH elastic_det and
coupled_det, takes the FIRST sign flip (the fundamental), and bisects
it. This was verified interactively across the whole sweep before
being written into this probe -- every point returned exactly one
flip in that window for both dets.

PRE-REGISTERED PASS/FAIL:
  G_found: every H in the sweep returns exactly one elastic root and
     one coupled root in its scan window (no missed/extra fundamental).
  G_theorem: SC > elastic at every H (stiffening theorem holds across
     the whole thickness range, not just the validated point).
  G_longwave: at the THINNEST H in the sweep (H/r_o=0.005), the
     computed frequency split is within a FACTOR OF 2 of the closed-
     form long-wave estimate 12*e31_bar^2/(pi^2*Xi33_bar*c11_bar)/2
     (the /2 converts the D-ratio estimate to a frequency-ratio
     estimate, omega ~ sqrt(D)). This is an orientation check, not a
     tight numerical gate -- PAPER5_DERIVATION.md Sec 6 itself frames
     the ~2.2% number as "orientation ... not a bar to tune to", and
     the annulus's finite r_i/r_o=1/6 means even the thinnest point is
     not the infinite-long-wave limit.
  H-independence trend (reported, not gated -- this IS the finding):
     the split percentage across the whole sweep, and whether/where it
     departs measurably from the thinnest point's value.

Standing bans respected: no SOLVER_VERSION bump, no OC solver, no
joint-fit, no .tex, no new plate_solver/*.py code (pure probe).

Cost model (this session, dps~30-40 interactively): one 250-step scan
of elastic_det or coupled_det over the guessed window costs roughly
7-25s (grows mildly with H since the Bessel arguments involved grow);
one bisect (40 iters) is a few seconds. 9 H points x (2 scans + 2
bisects) is comfortably a few minutes serial at interactive dps;
production dps=60 cluster archival will cost more per call but the
same total call count.
"""
import json
import os
import sys
import time

sys.path.insert(0, os.environ.get("PKG_PATH",
                                  os.path.join(os.path.dirname(__file__), "..")))
import plate_solver as ps  # noqa: E402
from plate_solver.piezo_monolithic import PiezoMonolithicOutOfPlaneSolver  # noqa: E402

EXPECT_SOLVER_VERSION = os.environ.get("EXPECT_SOLVER_VERSION", "2026-07-10.s10")
DPS = int(os.environ.get("PIEZO_DPS", "60"))
INNER_ITERS = int(os.environ.get("PIEZO_INNER_ITERS", "45"))
SCAN_N = int(os.environ.get("PIEZO_SCAN_N", "250"))

# Duan Table 1 / Paper 4 BASE_KWARGS (Sec 18.200); r_i, r_o fixed,
# only H varies. e31/e33/X11/X33 unchanged from the closed inverse work.
R_I, R_O = 0.1, 0.6
C11E, C12E, C13E, C33E = 132e9, 71e9, 73e9, 115e9
RHO = 7500.0
E31_TRUE, E33 = -4.1, 14.1
X11, X33 = 7.124e-9, 5.841e-9

# Sec 18.201's own validated point: F-F n=0 elastic root at H=0.01.
OMEGA_ANCHOR_H001 = 462.3976497355817

H_OVER_RO_LIST = [0.005, 0.01, 1.0 / 60.0, 0.025, 0.04, 0.06, 0.09, 0.13, 0.18]

# Closed-form long-wave orientation estimate (PAPER5_DERIVATION.md
# Sec 6): D_SC/d - 1 = 12*e31_bar^2/(pi^2*Xi33_bar*c11_bar), independent
# of H at this order; frequency-ratio estimate is half that (omega ~
# sqrt(D)). Computed once, from raw constants, not copied from the
# derivation doc's rounded PZT-4 example.
import math  # noqa: E402

_e31_bar = E31_TRUE - (C13E / C33E) * E33
_c11_bar = C11E - C13E ** 2 / C33E
_Xi33_bar = X33 + E33 ** 2 / C33E
# 2026-09-23 LESSONS Sec 18.231: the package's default (consistent,
# quadratic Galerkin) projection has the EXACT thin-plate long-wave
# increment e31_bar^2/(Xi33_bar*c11_bar); the legacy 'duan' projection's
# was 12/pi^2 times that (the formula this line used before).
LONGWAVE_D_SPLIT = _e31_bar ** 2 / (_Xi33_bar * _c11_bar)
LONGWAVE_FREQ_SPLIT = 0.5 * LONGWAVE_D_SPLIT


def _solver(H, dps=DPS):
    return PiezoMonolithicOutOfPlaneSolver(
        r_i=R_I, r_o=R_O, H=H,
        C11E=C11E, C12E=C12E, C13E=C13E, C33E=C33E, rho=RHO,
        e31=E31_TRUE, e33=E33, X11=X11, X33=X33, dps=dps)


def _sign_flips(detfn, lo, hi, nstep):
    step = (hi - lo) / nstep
    xs = [lo + i * step for i in range(nstep + 1)]
    prev_x = xs[0]
    try:
        prev_s = float(detfn(prev_x, 0).real)
    except Exception as e:
        print("scan start failed at %.3f: %s" % (prev_x, e), flush=True)
        return []
    flips = []
    for x in xs[1:]:
        try:
            s = float(detfn(x, 0).real)
        except Exception:
            continue
        if prev_s != 0 and s != 0 and (prev_s > 0) != (s > 0):
            flips.append((prev_x, x))
        prev_x, prev_s = x, s
    return flips


def main():
    t_all = time.time()
    print("job start: dps=%s inner_iters=%s scan_n=%s"
          % (DPS, INNER_ITERS, SCAN_N), flush=True)
    if ps.SOLVER_VERSION != EXPECT_SOLVER_VERSION:
        print("PREFLIGHT FAIL: SOLVER_VERSION=%r != %r"
              % (ps.SOLVER_VERSION, EXPECT_SOLVER_VERSION), flush=True)
        sys.exit(2)
    print("preflight OK: SOLVER_VERSION=%s" % ps.SOLVER_VERSION, flush=True)
    if not hasattr(PiezoMonolithicOutOfPlaneSolver, "coupled_bisect"):
        print("PREFLIGHT FAIL: deployed PiezoMonolithicOutOfPlaneSolver "
              "has no coupled_bisect. Push plate_solver/piezo_monolithic.py "
              "before re-running.", flush=True)
        sys.exit(2)
    print("preflight OK: coupled_bisect present", flush=True)
    print("closed-form long-wave orientation: D-split=%.6f freq-split=%.6f "
          "(%.4f%%)" % (LONGWAVE_D_SPLIT, LONGWAVE_FREQ_SPLIT,
                         100.0 * LONGWAVE_FREQ_SPLIT), flush=True)

    log = []

    def note(msg):
        log.append(msg)
        print(msg, flush=True)

    rows = []
    g_found = True
    g_theorem = True
    for h_over_ro in H_OVER_RO_LIST:
        H = h_over_ro * R_O
        s = _solver(H)
        guess = OMEGA_ANCHOR_H001 * (H / 0.01)
        lo, hi = max(20.0, 0.1 * guess), 3.5 * guess

        t0 = time.time()
        el_flips = _sign_flips(s.elastic_det, lo, hi, SCAN_N)
        sc_flips = _sign_flips(s.coupled_det, lo, hi, SCAN_N)
        found_ok = (len(el_flips) == 1 and len(sc_flips) == 1)
        if not found_ok:
            g_found = False
            note("H/r_o=%.5f H=%.6f: UNEXPECTED flip count el=%d sc=%d "
                 "in window [%.2f,%.2f] -- not scoring this point"
                 % (h_over_ro, H, len(el_flips), len(sc_flips), lo, hi))
            rows.append(dict(h_over_ro=h_over_ro, H=H, found_ok=False))
            continue

        el = s.elastic_bisect(el_flips[0][0], el_flips[0][1], 0,
                               iters=INNER_ITERS)
        sc = s.coupled_bisect(sc_flips[0][0], sc_flips[0][1], 0,
                               iters=INNER_ITERS)
        split = (sc - el) / el
        if not (sc > el):
            g_theorem = False
        note("H/r_o=%.5f H=%.6f el=%.6f sc=%.6f split=%.5f%% (%.1fs)"
             % (h_over_ro, H, el, sc, 100.0 * split, time.time() - t0))
        rows.append(dict(h_over_ro=h_over_ro, H=H, found_ok=True,
                          omega_el=el, omega_sc=sc, split=split,
                          split_pct=100.0 * split,
                          elapsed_s=time.time() - t0))

    thin_rows = [r for r in rows if r.get("found_ok") and r["h_over_ro"] == 0.005]
    g_longwave = False
    thin_split = None
    if thin_rows:
        thin_split = thin_rows[0]["split"]
        ratio = thin_split / LONGWAVE_FREQ_SPLIT
        g_longwave = (0.5 <= ratio <= 2.0)
        note("G_longwave: thinnest split=%.5f%% vs closed-form long-wave "
             "%.5f%% ratio=%.4f pass(0.5-2.0x)=%s"
             % (100.0 * thin_split, 100.0 * LONGWAVE_FREQ_SPLIT, ratio,
                g_longwave))
    else:
        note("G_longwave: thinnest point was not scored (found_ok=False) "
             "-- cannot evaluate")

    scored = [r for r in rows if r.get("found_ok")]
    if len(scored) >= 2:
        thickest = scored[-1]
        thinnest = scored[0]
        drift_rel = ((thickest["split"] - thinnest["split"])
                     / thinnest["split"])
        note("H-independence trend: split at H/r_o=%.5f is %.5f%%; at "
             "H/r_o=%.5f is %.5f%% (relative drift %.4f%%)"
             % (thinnest["h_over_ro"], thinnest["split_pct"],
                thickest["h_over_ro"], thickest["split_pct"],
                100.0 * drift_rel))

    all_pass = g_found and g_theorem and g_longwave

    results = dict(
        longwave_D_split=LONGWAVE_D_SPLIT,
        longwave_freq_split=LONGWAVE_FREQ_SPLIT,
        rows=rows, g_found=g_found, g_theorem=g_theorem,
        g_longwave=g_longwave, thin_split=thin_split,
        all_pass=all_pass, elapsed_s=time.time() - t_all, log=log,
    )
    out_path = os.path.join(os.path.dirname(__file__) or ".",
                             "piezo_p5_thickness_sweep_results_e31m.json")
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)

    print("--- summary ---", flush=True)
    print("  g_found=%s g_theorem=%s g_longwave=%s"
          % (g_found, g_theorem, g_longwave), flush=True)
    print("SENTINEL " + ("PASS_ALL" if all_pass else "FAIL_ALL")
          + " all_pass=%s (%.1fs total)" % (all_pass, time.time() - t_all),
          flush=True)


if __name__ == "__main__":
    main()
