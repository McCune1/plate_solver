# e31 = -4.1 COPY (LESSONS Sec 18.246) of probe_piezo_p5_cf_boundary_2026-09-20.py.
# Written by migrate_e31m_p5.py. Only e31 changes, plus the windows that
# bracket the F-F SC root (471.36 -> 519.84 rad/s at -4.1, i.e. every
# window near that root shifts by +48.47), the e31 inverse bracket
# [1,8] -> [-8,-1] (wrong-sign control [-8,-1] -> [1,8]; the twin
# 2*8.95 - e31 moves to 22.0), relative errors divide by |E31_TRUE|,
# and output files carry the suffix _e31m. Fixed before running.
# -*- coding: utf-8 -*-
"""
probe_piezo_p5_cf_boundary_2026-09-20.py

Paper 5 roadmap Sec 8 item 6, "C-F boundary condition" -- the first of
the two remaining optional-cheap extras (n=0 overtones and thickness
sweep closed, Sec 18.206/18.207), picked up separately per the user's
own scoping ("we will tackle the other two one at a time").

UNLIKE the overtones/thickness-sweep probes, this one required real
new solver code: PiezoMonolithicOutOfPlaneSolver had only F-F (both
edges free) and C-C (both edges clamped) determinants. C-F/F-C mixed
mechanical edges did not exist. This probe validates the new methods
added to plate_solver/piezo_monolithic.py this session:

  _elastic_mixed_det(omega, n, inner, outer)   4x4, {'C','F'} per edge
  elastic_cf_det/bisect, elastic_fc_det/bisect thin wrappers
  _coupled_mixed_det(omega, n, inner, outer)   6x6 SC, {'C','F'} per edge
  cf_coupled_det/bisect, fc_coupled_det/bisect thin wrappers
  _signflip_bisect                              shared bisection loop

C-F = inner clamped (w=w'=0), outer free (M_rr=Q_r=0); F-C is the
swap. Same convention as PiezoOutOfPlaneSolver's existing C-F/F-C
family (Sec 18.192) and Paper 3's ring_disk.py. No open-circuit
variant is added or needed: Phase 0's OC=SC theorem (Sec 18.199) is a
property of the fully-electroded potential ansatz, independent of the
mechanical edge condition, so it covers C-F/F-C automatically -- this
respects the standing ban on a new OC frequency solver.

Design note on the row-selection code: rather than hand-write two new
4x4/6x6 determinants, the implementation generalizes the existing F-F/
C-C row-building into one _mixed_det per family that picks, at each
edge independently, either the (w, w') pair or the (M_rr, Q_r) pair
before building the same matrix/equilibrate/det pipeline already
validated for F-F and C-C. This means F-F and C-C are literally
special cases of the new code (inner=outer='F' or inner=outer='C'),
which gives a much stronger correctness gate than a fresh derivation
would: the new code must reproduce the OLD code's determinant values
BIT-IDENTICALLY under matching labels, not just "close."

PRE-REGISTERED PASS/FAIL:
  G_reduce (hard gate): _elastic_mixed_det(omega, n, 'F','F') ==
     elastic_det(omega, n) and (.., 'C','C') == elastic_cc_det(omega, n)
     EXACTLY (mpmath value equality, not a tolerance), at both n=0 and
     n=2; same for _coupled_mixed_det vs coupled_det/cc_coupled_det.
     This is the real proof that the new row-selection logic didn't
     silently change anything for the already-validated F-F/C-C paths.
  G_found: exactly one elastic and one coupled (SC) sign-flip found in
     each pre-registered window (CF: [1000, 2000]; FC: [400, 900] --
     found interactively this session, not tuned to a target value).
  G_theorem: SC > elastic for both CF and FC -- the Sec 18.199
     stiffening-sign theorem generalizes to mixed mechanical edges too,
     not just F-F/C-C.
  G1 (negative control): e31=0 (no e31/e33/X11/X33 supplied) raises
     ValueError from cf_coupled_det/fc_coupled_det, inherited from
     _require_coupled_consts() -- confirms the new coupled methods
     didn't bypass the existing singular-limit guard (Sec 18.149-style).
  G_worker: the four new worker functions added to plate_solver/
     workers.py (_piezo_mono_{elastic,coupled}_{cf,fc}_root_worker)
     reconstruct the solver from plain scalars and reproduce the
     in-process bisect result EXACTLY -- required before any of this
     can run under the cluster's ProcessPoolExecutor.

Also reported (not gated): where the CF/FC elastic fundamentals sit
relative to the FF/CC fundamentals. Pre-registered expectation before
running: NONE -- there is no a priori reason to expect the inner and
outer edges to contribute symmetrically to bending stiffness at n=0,
so whatever fraction is found is a genuine finding, not a confirmation
of a prediction.

No SOLVER_VERSION bump (new unused-by-default entry points, same
pattern as every other method in this module). No .tex, no FE, no
overtones/thickness-sweep re-litigation, no Y(omega) admittance (next,
separately, per the user's own scoping).
"""
import json
import os
import sys
import time

sys.path.insert(0, os.environ.get("PKG_PATH",
                                  os.path.join(os.path.dirname(__file__), "..")))
import plate_solver as ps  # noqa: E402
from plate_solver.piezo_monolithic import PiezoMonolithicOutOfPlaneSolver  # noqa: E402
from plate_solver.workers import (  # noqa: E402
    _piezo_mono_elastic_cf_root_worker,
    _piezo_mono_elastic_fc_root_worker,
    _piezo_mono_coupled_cf_root_worker,
    _piezo_mono_coupled_fc_root_worker,
)

EXPECT_SOLVER_VERSION = os.environ.get("EXPECT_SOLVER_VERSION", "2026-07-10.s10")
DPS = int(os.environ.get("PIEZO_DPS", "60"))
INNER_ITERS = int(os.environ.get("PIEZO_INNER_ITERS", "45"))
SCAN_N = int(os.environ.get("PIEZO_SCAN_N", "300"))

# Duan Table 1 / Paper 4 BASE_KWARGS (Sec 18.200), Duan radii (Sec
# 18.201), full ceramic thickness 2H = 0.02 m -- same fixed geometry/
# material as every other Paper 5 probe this session.
R_I, R_O, H = 0.1, 0.6, 0.01
C11E, C12E, C13E, C33E = 132e9, 71e9, 73e9, 115e9
RHO = 7500.0
E31, E33 = -4.1, 14.1
X11, X33 = 7.124e-9, 5.841e-9

# Reference F-F/C-C anchors (Sec 18.201/18.204), for the "where does
# CF/FC sit" report only -- not gated on these.
EL_FF, EL_CC = 462.3976497355817, 1726.889795840675

CF_SCAN_LO, CF_SCAN_HI = 1000.0, 2000.0
FC_SCAN_LO, FC_SCAN_HI = 400.0, 900.0

BASE_KWARGS = dict(
    r_i=R_I, r_o=R_O, H=H,
    C11E=C11E, C12E=C12E, C13E=C13E, C33E=C33E, rho=RHO,
    e31=E31, e33=E33, X11=X11, X33=X33, dps=DPS,
)


def _solver(coupled=True, dps=None):
    kw = dict(BASE_KWARGS)
    if dps is not None:
        kw["dps"] = dps
    if not coupled:
        for k in ("e31", "e33", "X11", "X33"):
            kw.pop(k)
    return PiezoMonolithicOutOfPlaneSolver(**kw)


def _sign_flips(detfn, lo, hi, nstep, n=0):
    step = (hi - lo) / nstep
    xs = [lo + i * step for i in range(nstep + 1)]
    prev_x = xs[0]
    prev_s = float(detfn(prev_x, n).real)
    flips = []
    for x in xs[1:]:
        s = float(detfn(x, n).real)
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
    if not hasattr(PiezoMonolithicOutOfPlaneSolver, "cf_coupled_bisect"):
        print("PREFLIGHT FAIL: deployed PiezoMonolithicOutOfPlaneSolver "
              "has no cf_coupled_bisect. Push plate_solver/"
              "piezo_monolithic.py and plate_solver/workers.py before "
              "re-running.", flush=True)
        sys.exit(2)
    print("preflight OK: cf_coupled_bisect present", flush=True)

    log = []

    def note(msg):
        log.append(msg)
        print(msg, flush=True)

    s = _solver(coupled=True)
    s_el = _solver(coupled=False)

    # ---------- G_reduce: bit-identical reduction to F-F / C-C ----------
    t0 = time.time()
    reduce_pass = True
    reduce_rows = []
    for n, omega in ((0, 500.0), (2, 1500.0)):
        checks = [
            ("elastic FF", s_el.elastic_det(omega, n),
             s_el._elastic_mixed_det(omega, n, "F", "F")),
            ("elastic CC", s_el.elastic_cc_det(omega, n),
             s_el._elastic_mixed_det(omega, n, "C", "C")),
            ("coupled FF", s.coupled_det(omega, n),
             s._coupled_mixed_det(omega, n, "F", "F")),
            ("coupled CC", s.cc_coupled_det(omega, n),
             s._coupled_mixed_det(omega, n, "C", "C")),
        ]
        for label, d_old, d_new in checks:
            ok = (d_old == d_new)
            reduce_pass = reduce_pass and ok
            reduce_rows.append(dict(n=n, omega=omega, label=label, equal=ok))
            note("G_reduce n=%d %s: old==new -> %s" % (n, label, ok))
    note("G_reduce: all_equal=%s (%.1fs)" % (reduce_pass, time.time() - t0))

    # ---------- root-finding: CF and FC, elastic and coupled ----------
    def leg(name, elastic_det, elastic_bisect, coupled_det, coupled_bisect,
            lo, hi):
        t0 = time.time()
        el_flips = _sign_flips(elastic_det, lo, hi, SCAN_N)
        sc_flips = _sign_flips(coupled_det, lo, hi, SCAN_N)
        note("%s scan: %d elastic flips, %d coupled flips in [%.1f, %.1f] (%.1fs)"
             % (name, len(el_flips), len(sc_flips), lo, hi, time.time() - t0))
        el_root = (elastic_bisect(*el_flips[0], 0, iters=INNER_ITERS)
                   if el_flips else float("nan"))
        sc_root = (coupled_bisect(*sc_flips[0], 0, iters=INNER_ITERS)
                   if sc_flips else float("nan"))
        found = (len(el_flips) == 1 and len(sc_flips) == 1)
        theorem = found and (sc_root > el_root)
        split = ((sc_root - el_root) / el_root
                 if found else float("nan"))
        note("%s: el=%.6f sc=%.6f split=%.5f%% found=%s theorem=%s (%.1fs)"
             % (name, el_root, sc_root, 100.0 * split, found, theorem,
                time.time() - t0))
        return dict(name=name, el_root=el_root, sc_root=sc_root,
                    split=split, found=found, theorem=theorem,
                    n_el_flips=len(el_flips), n_sc_flips=len(sc_flips))

    cf = leg("CF", s_el.elastic_cf_det, s_el.elastic_cf_bisect,
              s.cf_coupled_det, s.cf_coupled_bisect,
              CF_SCAN_LO, CF_SCAN_HI)
    fc = leg("FC", s_el.elastic_fc_det, s_el.elastic_fc_bisect,
              s.fc_coupled_det, s.fc_coupled_bisect,
              FC_SCAN_LO, FC_SCAN_HI)

    g_found = cf["found"] and fc["found"]
    g_theorem = cf["theorem"] and fc["theorem"]

    # ---------- G1: negative control, e31=0 ----------
    try:
        s_el.cf_coupled_det(cf["sc_root"] if cf["found"] else 1644.0, 0)
        g1_pass = False
        note("G1 NEG CONTROL FAIL: cf_coupled_det did not raise at e31=0")
    except ValueError as e:
        g1_pass = True
        note("G1 neg control (CF) PASS: %s" % (str(e)[:70],))
    try:
        s_el.fc_coupled_det(fc["sc_root"] if fc["found"] else 566.0, 0)
        g1_pass = False
        note("G1 NEG CONTROL FAIL: fc_coupled_det did not raise at e31=0")
    except ValueError as e:
        note("G1 neg control (FC) PASS: %s" % (str(e)[:70],))

    # ---------- G_worker: worker-path scalar reconstruction ----------
    t0 = time.time()
    worker_pass = True
    if cf["found"]:
        lo_el, hi_el = cf["el_root"] - 5.0, cf["el_root"] + 5.0
        lo_sc, hi_sc = cf["sc_root"] - 2.0, cf["sc_root"] + 2.0
        in_el = s_el.elastic_cf_bisect(lo_el, hi_el, 0, iters=INNER_ITERS)
        via_el = _piezo_mono_elastic_cf_root_worker(
            (R_I, R_O, H, C11E, C12E, C13E, C33E, RHO, DPS,
             0, lo_el, hi_el, INNER_ITERS))
        in_cp = s.cf_coupled_bisect(lo_sc, hi_sc, 0, iters=INNER_ITERS)
        via_cp = _piezo_mono_coupled_cf_root_worker(
            (R_I, R_O, H, C11E, C12E, C13E, C33E, RHO,
             E31, E33, X11, X33, DPS, 0, lo_sc, hi_sc, INNER_ITERS))
        worker_pass = worker_pass and (in_el == via_el) and (in_cp == via_cp)
        note("G_worker CF: elastic in==worker %s; coupled in==worker %s"
             % (in_el == via_el, in_cp == via_cp))
    if fc["found"]:
        lo_el, hi_el = fc["el_root"] - 5.0, fc["el_root"] + 5.0
        lo_sc, hi_sc = fc["sc_root"] - 2.0, fc["sc_root"] + 2.0
        in_el = s_el.elastic_fc_bisect(lo_el, hi_el, 0, iters=INNER_ITERS)
        via_el = _piezo_mono_elastic_fc_root_worker(
            (R_I, R_O, H, C11E, C12E, C13E, C33E, RHO, DPS,
             0, lo_el, hi_el, INNER_ITERS))
        in_cp = s.fc_coupled_bisect(lo_sc, hi_sc, 0, iters=INNER_ITERS)
        via_cp = _piezo_mono_coupled_fc_root_worker(
            (R_I, R_O, H, C11E, C12E, C13E, C33E, RHO,
             E31, E33, X11, X33, DPS, 0, lo_sc, hi_sc, INNER_ITERS))
        worker_pass = worker_pass and (in_el == via_el) and (in_cp == via_cp)
        note("G_worker FC: elastic in==worker %s; coupled in==worker %s"
             % (in_el == via_el, in_cp == via_cp))
    note("G_worker: all_pass=%s (%.1fs)" % (worker_pass, time.time() - t0))

    # ---------- where do CF/FC sit relative to FF/CC (reported only) ----------
    span = EL_CC - EL_FF
    if cf["found"]:
        frac_cf = (cf["el_root"] - EL_FF) / span
        note("CF elastic frequency is %.2f%% of the way from FF (%.4f) "
             "to CC (%.4f)" % (100.0 * frac_cf, EL_FF, EL_CC))
    if fc["found"]:
        frac_fc = (fc["el_root"] - EL_FF) / span
        note("FC elastic frequency is %.2f%% of the way from FF (%.4f) "
             "to CC (%.4f)" % (100.0 * frac_fc, EL_FF, EL_CC))

    all_pass = reduce_pass and g_found and g_theorem and g1_pass and worker_pass

    results = dict(
        cf=cf, fc=fc, reduce_pass=reduce_pass, reduce_rows=reduce_rows,
        g_found=g_found, g_theorem=g_theorem, g1_pass=g1_pass,
        worker_pass=worker_pass, all_pass=all_pass,
        elapsed_s=time.time() - t_all, log=log,
    )
    out_path = os.path.join(os.path.dirname(__file__) or ".",
                             "piezo_p5_cf_boundary_results_e31m.json")
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)

    print("--- summary ---", flush=True)
    print("  G_reduce=%s G_found=%s G_theorem=%s G1=%s G_worker=%s"
          % (reduce_pass, g_found, g_theorem, g1_pass, worker_pass),
          flush=True)
    print("SENTINEL " + ("PASS_ALL" if all_pass else "FAIL_ALL")
          + " (%.1fs total)" % (time.time() - t_all), flush=True)


if __name__ == "__main__":
    main()
