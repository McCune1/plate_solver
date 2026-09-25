# e31 = -4.1 COPY (LESSONS Sec 18.246) of probe_piezo_p5_n0_overtones_2026-09-20.py.
# Written by migrate_e31m_p5.py. Only e31 changes, plus the windows that
# bracket the F-F SC root (471.36 -> 519.84 rad/s at -4.1, i.e. every
# window near that root shifts by +48.47), the e31 inverse bracket
# [1,8] -> [-8,-1] (wrong-sign control [-8,-1] -> [1,8]; the twin
# 2*8.95 - e31 moves to 22.0), relative errors divide by |E31_TRUE|,
# and output files carry the suffix _e31m. Fixed before running.
# -*- coding: utf-8 -*-
"""
probe_piezo_p5_n0_overtones_2026-09-20.py

Paper 5 roadmap Sec 8 item 6, "optional cheap: n=0 overtones" -- a pure
probe-level scan on the existing PiezoMonolithicOutOfPlaneSolver, no
new solver code. Direct follow-up to the closed inverse work
(LESSONS_LEARNED.md Sec 18.204/18.205): this asks whether the SC-vs-
elastic F-F/C-C n=0 fundamental's stiffening theorem and reflection-
twin identity (Sec 18.204) generalize to higher n=0 modes (more radial
nodal circles), and whether the SC-vs-elastic split percentage stays
near the fundamental's ~2.3%-2.35%(F-F)/2.16%(C-C) or drifts.

Same style as Paper 4's probe_piezo_p4_ff_n0_overtones_wide_2026-09-16.py
(sign-flip scan -> bisect each bracket -> pair elastic with coupled by
index), re-pointed at SC-vs-elastic (this stacking has no OC channel,
Phase 0 theorem, Sec 18.199/18.203) and covering BOTH F-F and C-C
(the fixed geometry/material of Sec 18.200-18.201: PZT-4 Duan Table 1
/ Paper 4 BASE_KWARGS, r_i=0.1, r_o=0.6, H=0.01).

Scan windows (found interactively this session, not tuned to a target
count -- report whatever the scan actually finds):
  F-F n=0: (455, 6000) rad/s -- expected to contain the fundamental
    (~462.4/473.2) plus two further pairs (~2170/2211, ~5271/5374).
  C-C n=0: (1700, 6000) rad/s -- expected to contain the fundamental
    (~1726.9/1764.3) plus one further pair (~4777/4879).
If the actual scan finds a different count, that is itself reported
(G_found's bar is ">=1 pair beyond the fundamental" per leg, not a
specific number) -- do not retune the window to force a specific
count.

PRE-REGISTERED PASS/FAIL:
  G_found: at least 2 pairs found on F-F, at least 1 pair (i.e. the
     fundamental plus >=1 more) found on C-C, in the windows above.
  G_theorem: SC > elastic for EVERY found pair, both legs -- the
     Sec 18.199/PAPER5_DERIVATION.md Sec 6 stiffening-sign theorem
     (D_SC > d) generalizes beyond the n=0 fundamental, not just
     asserted there.
  G_reflection: the e31 reflection twin (Sec 18.204:
     twin = 2*(C13E/C33E)*e33 - e31, = 13.800869565217392 at e31=4.1)
     reproduces EVERY found SC root to <1e-9 relative, both legs --
     the reflection-symmetry proof in Sec 18.204 was derived for a
     generic omega, so it should hold at every root, not just the
     fundamental; this is a genuine test of that generality, not a
     re-run of an already-known number.
  G2 (reported, not gated, every pair): forward FD domega/de31;
     propagated e31 uncertainty at 0.01%/0.1% assumed frequency
     precision. Pre-registered expectation: informational only -- the
     fundamental's own sensitivity (Sec 18.204: FF 0.26%/2.65%, CC
     0.28%/2.82%) is not asserted to hold at higher modes; report
     whatever is found.
  Split-percentage trend (reported, not gated): whether the SC-vs-
     elastic split shrinks, grows, or stays flat with mode index. No
     a priori theoretical prediction exists for this beyond the
     fundamental's own long-wave estimate (PAPER5_DERIVATION.md Sec 6),
     so this is genuinely open until scored.

Standing bans respected: no SOLVER_VERSION bump, no OC solver, no
joint-fit, no .tex, no new plate_solver/*.py code (pure probe).

Cost model (this session, dps=40, INNER_ITERS=40, scan steps below):
one 250-step scan of a smooth 4x4/6x6 det over a ~1 octave-ish window
costs roughly 5-20s depending on omega magnitude (larger omega ->
larger Bessel arguments -> somewhat more mpmath work per eval); one
bisect call (40 iters) is a few seconds. Total this probe: ~2 scans
(elastic, coupled) x 2 legs plus ~6-8 bisects, comfortably a few
minutes serial. Backgrounding is not reliable in this environment
(each shell is fresh), so this was run in pieces interactively and
the numbers below the cluster archival section are the production
(dps=60) cluster run, not the interactive smoke values.
"""
import json
import math
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
FF_SCAN_LO, FF_SCAN_HI = 455.0, 6000.0
FF_SCAN_N = int(os.environ.get("PIEZO_FF_SCAN_N", "600"))
CC_SCAN_LO, CC_SCAN_HI = 1700.0, 6000.0
CC_SCAN_N = int(os.environ.get("PIEZO_CC_SCAN_N", "500"))
REFLECTION_TOL_REL = 1e-9

# Duan Table 1 / Paper 4 BASE_KWARGS (Sec 18.200), Duan radii (Sec
# 18.201), full ceramic thickness 2H = 0.02 m -- unchanged from the
# closed inverse work.
R_I, R_O, H = 0.1, 0.6, 0.01
C11E, C12E, C13E, C33E = 132e9, 71e9, 73e9, 115e9
RHO = 7500.0
E31_TRUE, E33 = -4.1, 14.1
X11, X33 = 7.124e-9, 5.841e-9
REFLECTION_TWIN = 2.0 * (C13E / C33E) * E33 - E31_TRUE

BASE_KWARGS = dict(
    r_i=R_I, r_o=R_O, H=H,
    C11E=C11E, C12E=C12E, C13E=C13E, C33E=C33E, rho=RHO,
    e33=E33, X11=X11, X33=X33, dps=DPS,
)


def _solver(e31=E31_TRUE):
    kw = dict(BASE_KWARGS)
    kw["e31"] = e31
    return PiezoMonolithicOutOfPlaneSolver(**kw)


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
    n_fail = 0
    for x in xs[1:]:
        try:
            s = float(detfn(x, 0).real)
        except Exception:
            n_fail += 1
            continue
        if prev_s != 0 and s != 0 and (prev_s > 0) != (s > 0):
            flips.append((prev_x, x))
        prev_x, prev_s = x, s
    if n_fail:
        print("scan: %d det evals failed; flips kept anyway" % n_fail,
              flush=True)
    return flips


def _run_leg(name, s, elastic_det, elastic_bisect, coupled_det,
             coupled_bisect, scan_lo, scan_hi, scan_n, note):
    t0 = time.time()
    el_flips = _sign_flips(elastic_det, scan_lo, scan_hi, scan_n)
    sc_flips = _sign_flips(coupled_det, scan_lo, scan_hi, scan_n)
    note("%s scan: %d elastic flips, %d coupled flips in [%.1f, %.1f] (%.1fs)"
         % (name, len(el_flips), len(sc_flips), scan_lo, scan_hi,
            time.time() - t0))

    el_roots, sc_roots = [], []
    for lo, hi in el_flips:
        try:
            el_roots.append(elastic_bisect(lo, hi, 0, iters=INNER_ITERS))
        except Exception as e:
            note("%s elastic bisect failed on [%.3f,%.3f]: %s"
                 % (name, lo, hi, e))
    for lo, hi in sc_flips:
        try:
            sc_roots.append(coupled_bisect(lo, hi, 0, iters=INNER_ITERS))
        except Exception as e:
            note("%s coupled bisect failed on [%.3f,%.3f]: %s"
                 % (name, lo, hi, e))
    note("%s elastic roots: %s" % (name, ["%.6f" % x for x in el_roots]))
    note("%s coupled roots: %s" % (name, ["%.6f" % x for x in sc_roots]))

    n_pair = min(len(el_roots), len(sc_roots))
    rows = []
    theorem_pass = True
    reflection_pass = True
    for i in range(n_pair):
        w_el, w_sc = el_roots[i], sc_roots[i]
        split = (w_sc - w_el) / w_el
        if not (w_sc > w_el):
            theorem_pass = False

        t0 = time.time()
        # reflection check: rebuild solver at REFLECTION_TWIN and bisect
        # in a bracket centered on this root.
        s_twin = _solver(REFLECTION_TWIN)
        cb_twin = getattr(s_twin, coupled_bisect.__name__)
        lo_t, hi_t = 0.98 * w_sc, 1.02 * w_sc
        try:
            w_twin_root = cb_twin(lo_t, hi_t, 0, iters=INNER_ITERS)
            refl_rel = abs(w_twin_root - w_sc) / w_sc
        except Exception as e:
            w_twin_root, refl_rel = float("nan"), float("inf")
            note("%s reflection check failed on pair %d: %s" % (name, i + 1, e))
        if not (refl_rel < REFLECTION_TOL_REL):
            reflection_pass = False

        dp = 0.01 * E31_TRUE
        s_plus, s_minus = _solver(E31_TRUE + dp), _solver(E31_TRUE - dp)
        cb_plus = getattr(s_plus, coupled_bisect.__name__)
        cb_minus = getattr(s_minus, coupled_bisect.__name__)
        lo_g, hi_g = 0.95 * w_sc, 1.05 * w_sc
        try:
            w_plus = cb_plus(lo_g, hi_g, 0, iters=INNER_ITERS)
            w_minus = cb_minus(lo_g, hi_g, 0, iters=INNER_ITERS)
            domega_de31 = (w_plus - w_minus) / (2 * dp)
        except Exception as e:
            domega_de31 = float("nan")
            note("%s G2 failed on pair %d: %s" % (name, i + 1, e))

        g2 = {}
        for rel_prec in (1e-4, 1e-3):
            d_om = rel_prec * w_sc
            dp_prop = (d_om / abs(domega_de31)
                       if domega_de31 == domega_de31 and domega_de31 != 0
                       else float("nan"))
            g2["%g" % rel_prec] = dict(
                propagated_de31=dp_prop,
                propagated_de31_rel=(dp_prop / abs(E31_TRUE)
                                      if dp_prop == dp_prop else float("nan")),
            )

        row = dict(
            index=i + 1, omega_el=w_el, omega_sc=w_sc,
            f_el_hz=w_el / (2 * math.pi), f_sc_hz=w_sc / (2 * math.pi),
            split=split, split_pct=100.0 * split,
            reflection_omega_at_twin=w_twin_root, reflection_rel=refl_rel,
            domega_de31=domega_de31, g2=g2, elapsed_s=time.time() - t0,
        )
        rows.append(row)
        note("%s pair %d: el=%.6f sc=%.6f split=%.5f%% reflection_rel=%.3e "
             "domega/de31=%.6e g2@0.01%%=%.4f%% (%.1fs)"
             % (name, i + 1, w_el, w_sc, 100.0 * split, refl_rel,
                domega_de31,
                row["g2"]["0.0001"]["propagated_de31_rel"] * 100
                if row["g2"]["0.0001"]["propagated_de31_rel"] == row["g2"]["0.0001"]["propagated_de31_rel"]
                else float("nan"),
                time.time() - t0))

    return dict(n_pair=n_pair, theorem_pass=theorem_pass,
                reflection_pass=reflection_pass, rows=rows)


def main():
    t_all = time.time()
    print("job start: dps=%s inner_iters=%s ff_scan_n=%s cc_scan_n=%s"
          % (DPS, INNER_ITERS, FF_SCAN_N, CC_SCAN_N), flush=True)
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

    log = []

    def note(msg):
        log.append(msg)
        print(msg, flush=True)

    s = _solver(E31_TRUE)

    ff = _run_leg("FF", s, s.elastic_det, s.elastic_bisect,
                  s.coupled_det, s.coupled_bisect,
                  FF_SCAN_LO, FF_SCAN_HI, FF_SCAN_N, note)
    cc = _run_leg("CC", s, s.elastic_cc_det, s.elastic_cc_bisect,
                  s.cc_coupled_det, s.cc_coupled_bisect,
                  CC_SCAN_LO, CC_SCAN_HI, CC_SCAN_N, note)

    g_found_ff = ff["n_pair"] >= 2
    # "at least 1 pair beyond the fundamental" on CC means n_pair >= 2;
    # the fundamental alone (n_pair==1) would not satisfy the roadmap's
    # own framing of an "overtone" scan finding something new.
    g_found_cc = cc["n_pair"] >= 2
    g_found = g_found_ff and g_found_cc
    all_pass = (g_found and ff["theorem_pass"] and cc["theorem_pass"]
                and ff["reflection_pass"] and cc["reflection_pass"])

    results = dict(
        ff=ff, cc=cc, g_found_ff=g_found_ff, g_found_cc=g_found_cc,
        g_found=g_found, all_pass=all_pass,
        elapsed_s=time.time() - t_all, log=log,
    )
    out_path = os.path.join(os.path.dirname(__file__) or ".",
                             "piezo_p5_n0_overtones_results_e31m.json")
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)

    print("--- summary ---", flush=True)
    print("  FF: n_pair=%d theorem=%s reflection=%s"
          % (ff["n_pair"], ff["theorem_pass"], ff["reflection_pass"]),
          flush=True)
    print("  CC: n_pair=%d theorem=%s reflection=%s"
          % (cc["n_pair"], cc["theorem_pass"], cc["reflection_pass"]),
          flush=True)
    print("SENTINEL " + ("PASS_ALL" if all_pass else "FAIL_ALL")
          + " g_found=%s (%.1fs total)" % (g_found, time.time() - t_all),
          flush=True)


if __name__ == "__main__":
    main()
