# -*- coding: utf-8 -*-
"""
probe_oop_shi2014_refine_hits_v1.py -- DIAGNOSTIC ONLY. No package changes,
no SOLVER_VERSION impact. Writes a plain-text log to stdout.

PURPOSE: upgrade Sec. 6.9 of PAPER1 from "a singularity exists inside the
predicted window" to "an independently refined root at Omega*, err% vs the
published value" for the five Shi2014 Table-3 targets already confirmed as
HITs by job 2410389.

WHY THIS IS THE HIGHEST-VALUE REMAINING RUN: job 2410389 reported rel% =
0.000 for all five HITs, which looks like a perfect match but is an
artifact of how that job worked -- its scan grid CONTAINED the predicted
centre as a grid point, and arg-max simply selected it. No refinement ever
happened. So the paper (correctly) can only claim a coarse window
confirmation, and must keep saying its primary flexural tables are
FE-validated only. Refining each dip to its actual minimum, and reporting
the distance from there to the published value, is what converts Sec. 6.9
into a genuine literature benchmark comparable to Sec. 6.7/6.8.

Target 3 (38.428) is DELIBERATELY EXCLUDED. Job 2411852 established with a
33-point 3-stage dense scan of the tight +-5% window that its dip bottoms
at log_A = -3.297 and never reaches -3.5 anywhere nearby. It is a confirmed
genuine near-miss, not an under-resolution artifact; refining it would only
polish a sub-threshold feature.

METHOD, per target: golden-section minimisation of log10(sigma_min) inside
a tight bracket around the predicted centre, then a bracket-shift stability
test. The stability test is the point of the exercise -- this project has
been burned before by golden-section locking onto a sharper neighbouring
artifact instead of the true zero (see the fixed-offset sensitivity-map
mitigation in Sec. 5 of the paper), so a refined root that MOVES when the
bracket is nudged is not a refined root.

PRE-REGISTERED INTERPRETATION (written before the job runs, honour it after):
  For each target independently --
  - PASS if golden-section converges to Omega* with log_A(Omega*) <= -3.5,
    AND |Omega_lit(Omega*) - published| / published < 2%, AND both
    bracket-shift replicas land within 0.1% of Omega*.
    -> quote Omega* and its err% in Sec. 6.9's table.
  - PARTIAL if it converges and is deep enough but a bracket-shift replica
    moves by more than 0.1% -> the minimum is not isolated at this
    resolution. Report as window-confirmed only; do NOT quote an err%.
    This is the lock-on failure mode, and reporting a number here would be
    exactly the "good enough error% is not a confirmed match" mistake the
    project's own methodology warns against.
  - FAIL / kill if log_A at the converged point drives to a floor
    (log_A <= -12) or returns the -300 sentinel. At dps=40 a genuine
    boundary-determinant root on this geometry reads about -3.3 to -4.9
    (job 2410389); -12 or deeper is a degenerate/zero matrix, not a better
    root. Report and stop -- do not treat a deeper number as a better one.
  - If a target's refined Omega* lands OUTSIDE its own bracket, that is a
    bug in this probe, not a result. Report it as such.

  Aggregate: 5/5 PASS -> Sec. 6.9 can be rewritten as a refined-root
  benchmark. 3-4 PASS -> report per-target, mixed. <=2 PASS -> the coarse
  window-confirmation framing stays as-is and this avenue is closed.

FIDELITY GATE: before any refinement, re-evaluate log_A at each of the five
predicted centres and confirm each reproduces job 2410389's tabulated depth
to within 0.05 decades. If any centre disagrees, the detector has moved
since 2410389 and the whole run is void.

COST: 5 targets x (1 gate + ~12 golden evaluations + 2 x ~10 shift
replicas) ~= 5 x 33 = 165 evaluations at ~57 s = ~2.6 CPU-hours. At 5
workers (1 target each) that is ~35 min wall. Budgeted 2 h below.
"""
from __future__ import annotations
import os
import sys
import time
import math
from concurrent.futures import ProcessPoolExecutor, as_completed

sys.path.insert(0, os.environ.get("PKG_PATH", "."))
os.environ.setdefault("DPS", "40")

EXPECT_VER = os.environ.get("EXPECT_SOLVER_VERSION", "2026-07-10.s10")

# Geometry for Ri/Ro = 0.4 -- identical to probe_oop_shi2014_table3_r04_v3
# and probe_oop_shi2014_table3_t3_local_v1.
EPS = 0.4
R0_2B = (1.0 + EPS) / (2.0 * (1.0 - EPS))     # 1.1666666...
TWO_THETA = 0.5
N_DOFS = 20
XMAX = 18.0
M = 80
N_QUAD = 30

# Omega_lit = FACTOR * Omega_native, established and scale-checked before
# job 2410389 ran.
FACTOR = 27.4155677808

DEPTH_THR = -3.5          # PASS bar on log10 sigma_min
FLOOR_THR = -12.0         # deeper than this is degenerate, not better
ERR_BAR = 2.0             # percent, Omega_lit vs published
SHIFT_TOL = 0.1           # percent, bracket-shift stability
BRACKET_HALF = 0.025      # +-2.5% of centre; job 2410389's own step size

# (label, published Omega_lit, predicted native Omega, job-2410389 depth)
# Target 3 (38.428, -3.264) excluded -- see docstring.
TARGETS = [
    ("T1", 15.647, 0.57073412, -4.857),
    ("T2", 23.576, 0.85994936, -4.205),
    ("T4", 53.649, 1.95688087, -3.983),
    ("T5", 63.390, 2.31218994, -4.164),
    ("T6", 70.739, 2.58024932, -3.952),
]

GOLDEN_ITERS = 12
INV_PHI = (math.sqrt(5.0) - 1.0) / 2.0


def _build():
    import plate_solver as ps
    from plate_solver import (OutOfPlaneSolver, IsotropicMaterial,
                              make_geometry, FreeFreeOOP)
    mat = IsotropicMaterial(E=210e9, nu=0.30, rho=7800.0)
    geom = make_geometry(R0_2B, TWO_THETA)
    if abs(float(geom.R_i) / float(geom.R_o) - 0.4) > 1e-6:
        raise RuntimeError("geometry Ri/Ro != 0.4")
    return OutOfPlaneSolver(geom, mat, M=M, n_quad=N_QUAD,
                            boundary=FreeFreeOOP())


def eval_point(solver, Om):
    """Verbatim evaluation from job 2410389 / 2411852. log_A is ALREADY
    log10(sigma_min) -- do not re-log (that double-log wrapper is what
    produced job 2410227's false 6/6 PASS with -300 everywhere)."""
    from plate_solver.detectors import full_search, select_fill
    raw = full_search(solver.fast, Om, xmax=XMAX)
    sel, n_ret = select_fill(raw, N_DOFS)
    if n_ret < max(4, N_DOFS // 2):
        return float("nan"), n_ret, len(raw)
    try:
        return float(solver.sigma_min(Om, sel)), n_ret, len(raw)
    except Exception:
        return float("nan"), n_ret, len(raw)


def golden_min(solver, lo, hi, iters, tag, log):
    """Golden-section minimisation of log_A on [lo, hi]."""
    a, b = lo, hi
    c = b - INV_PHI * (b - a)
    d = a + INV_PHI * (b - a)
    fc, _, _ = eval_point(solver, c)
    fd, _, _ = eval_point(solver, d)
    for _ in range(iters):
        if math.isnan(fc) or math.isnan(fd):
            return float("nan"), float("nan")
        if fc < fd:
            b, d, fd = d, c, fc
            c = b - INV_PHI * (b - a)
            fc, _, _ = eval_point(solver, c)
            log.append((tag, c, fc))
        else:
            a, c, fc = c, d, fd
            d = a + INV_PHI * (b - a)
            fd, _, _ = eval_point(solver, d)
            log.append((tag, d, fd))
    return (c, fc) if fc < fd else (d, fd)


def _worker(t):
    label, lit, nat, depth_2410389 = t
    t0 = time.time()
    out = dict(label=label, lit=lit, nat=nat, depth_ref=depth_2410389)
    try:
        solver = _build()

        # --- fidelity gate at the predicted centre -------------------
        g, n_ret, n_br = eval_point(solver, nat)
        out.update(gate_log=g, gate_n_ret=n_ret, gate_n_br=n_br)
        if math.isnan(g) or abs(g - depth_2410389) > 0.05:
            out.update(ok=False, verdict="GATE_FAIL",
                       note=f"centre reads {g:.4f}, job 2410389 had "
                            f"{depth_2410389:.4f}; detector has moved")
            out["dt"] = time.time() - t0
            return out

        # --- primary refinement -------------------------------------
        log = []
        lo, hi = nat * (1 - BRACKET_HALF), nat * (1 + BRACKET_HALF)
        om_star, f_star = golden_min(solver, lo, hi, GOLDEN_ITERS,
                                     "main", log)
        out.update(om_star=om_star, log_A=f_star, trace=log)

        if math.isnan(f_star):
            out.update(ok=False, verdict="EVAL_FAIL")
            out["dt"] = time.time() - t0
            return out
        if f_star <= FLOOR_THR:
            out.update(ok=True, verdict="FAIL_FLOOR",
                       note="converged into a numerical floor, not a root")
            out["dt"] = time.time() - t0
            return out
        if not (lo <= om_star <= hi):
            out.update(ok=False, verdict="PROBE_BUG",
                       note="refined point outside its own bracket")
            out["dt"] = time.time() - t0
            return out

        # --- bracket-shift stability replicas ------------------------
        reps = []
        for shift in (+0.25 * BRACKET_HALF, -0.25 * BRACKET_HALF):
            c2 = nat * (1 + shift)
            l2, h2 = c2 * (1 - BRACKET_HALF), c2 * (1 + BRACKET_HALF)
            o2, f2 = golden_min(solver, l2, h2, GOLDEN_ITERS, "shift", log)
            reps.append((o2, f2))
        out["replicas"] = reps
        moves = [abs(o2 - om_star) / om_star * 100.0
                 for o2, f2 in reps if not math.isnan(o2)]
        out["max_move_pct"] = max(moves) if moves else float("nan")

        om_lit = FACTOR * om_star
        err = abs(om_lit - lit) / lit * 100.0
        out.update(om_lit=om_lit, err_pct=err, ok=True)

        deep_enough = f_star <= DEPTH_THR
        close_enough = err < ERR_BAR
        stable = (not math.isnan(out["max_move_pct"])
                  and out["max_move_pct"] <= SHIFT_TOL)
        if deep_enough and close_enough and stable:
            out["verdict"] = "PASS"
        elif deep_enough and close_enough and not stable:
            out["verdict"] = "PARTIAL_UNSTABLE"
        elif deep_enough and not close_enough:
            out["verdict"] = "PARTIAL_FAR"
        else:
            out["verdict"] = "PARTIAL_SHALLOW"
    except Exception as exc:
        out.update(ok=False, verdict="EXC", note=f"{exc}")
    out["dt"] = time.time() - t0
    return out


def main():
    print("=" * 78)
    print("  probe_oop_shi2014_refine_hits_v1 -- start "
          + time.strftime("%Y-%m-%d %H:%M:%S"))
    print("=" * 78, flush=True)

    import plate_solver as ps
    print(f"SOLVER_VERSION={getattr(ps,'SOLVER_VERSION','?')} "
          f"(expect {EXPECT_VER})  DPS={os.environ['DPS']}", flush=True)
    if getattr(ps, "SOLVER_VERSION", None) != EXPECT_VER:
        print("PREFLIGHT FAIL: SOLVER_VERSION mismatch")
        raise SystemExit(2)
    print(f"geometry: R0_2B={R0_2B:.12g}  2Theta/pi={TWO_THETA}  "
          f"Ri/Ro={EPS}  N_DOFS={N_DOFS}  XMAX={XMAX}")
    print(f"bars: depth<={DEPTH_THR}  err<{ERR_BAR}%  "
          f"bracket-shift<={SHIFT_TOL}%  floor={FLOOR_THR}")
    print(f"{len(TARGETS)} targets (38.428 excluded, job 2411852 "
          f"confirmed genuine near-miss)\n", flush=True)

    n_workers = int(os.environ.get("N_WORKERS", "5"))
    results = []
    with ProcessPoolExecutor(max_workers=min(n_workers, len(TARGETS))) as ex:
        futs = {ex.submit(_worker, t): t for t in TARGETS}
        for fut in as_completed(futs):
            r = fut.result()
            results.append(r)
            if r.get("verdict") in ("GATE_FAIL", "EVAL_FAIL", "PROBE_BUG",
                                    "EXC"):
                print(f"  {r['label']}: {r['verdict']} -- "
                      f"{r.get('note','')} [{r['dt']:.0f}s]", flush=True)
            else:
                print(f"  {r['label']}: Om*={r.get('om_star',float('nan')):.8f} "
                      f"log_A={r.get('log_A',float('nan')):+.4f} "
                      f"Om_lit={r.get('om_lit',float('nan')):.4f} "
                      f"(pub {r['lit']}) err={r.get('err_pct',float('nan')):.3f}% "
                      f"shift={r.get('max_move_pct',float('nan')):.4f}% "
                      f"-> {r['verdict']} [{r['dt']:.0f}s]", flush=True)

    results.sort(key=lambda r: r["lit"])
    print("\n" + "=" * 78)
    print("  PRE-REGISTERED READING")
    print("=" * 78)
    print(f"{'lab':>4} {'published':>10} {'Om*':>12} {'Om_lit*':>10} "
          f"{'err%':>7} {'log_A':>8} {'shift%':>8}  verdict")
    for r in results:
        print(f"{r['label']:>4} {r['lit']:10.3f} "
              f"{r.get('om_star',float('nan')):12.8f} "
              f"{r.get('om_lit',float('nan')):10.4f} "
              f"{r.get('err_pct',float('nan')):7.3f} "
              f"{r.get('log_A',float('nan')):8.4f} "
              f"{r.get('max_move_pct',float('nan')):8.4f}  {r['verdict']}")

    n_pass = sum(1 for r in results if r["verdict"] == "PASS")
    print(f"\n  n_PASS = {n_pass}/{len(TARGETS)}")
    if n_pass == len(TARGETS):
        print("  -> ALL PASS. Sec. 6.9 may be rewritten as a refined-root")
        print("     benchmark quoting Om_lit* and err% per target.")
    elif n_pass >= 3:
        print("  -> MIXED. Quote err% only for the PASS rows; the rest stay")
        print("     window-confirmed. Do not average over the two classes.")
    else:
        print("  -> NOT SUPPORTED. Sec. 6.9's coarse window-confirmation")
        print("     framing stays exactly as written. This avenue is closed;")
        print("     do not re-refine the same windows.")
    print("\n  Reminder: a deeper log_A is not automatically a better root.")
    print("  Anything at or past", FLOOR_THR, "is a degenerate matrix.")


if __name__ == "__main__":
    main()
