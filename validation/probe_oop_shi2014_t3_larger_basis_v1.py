# -*- coding: utf-8 -*-
"""
probe_oop_shi2014_t3_larger_basis_v1.py -- DIAGNOSTIC ONLY.
No package changes, no SOLVER_VERSION impact. Plain-text log on stdout.

Rank 10 of PAPER1_CAVEATS / paper §6.9: Shi 2014 Table 3 target 38.428
is a genuine near-miss. Job 2411852's 33-point scan of the tight ±5%
window found a dip centered on the predicted frequency
(Ω=1.40170, log_A=-3.297) that never reached the -3.5 acceptance
depth. The remaining lever is a larger achieved basis at that same
point, not a finer scan of the same n=20 surface.

WHAT THIS IS NOT
----------------
It does not reopen the 12.5%-away feature at Ω=1.22647 (out of the
±5% window; already ruled unrelated). It does not re-run the n=20
33-point scan. The n=20 centre is only a fidelity gate.

METHOD
------
Same eval_point as jobs 2410389 / 2411852 / 2414020 (do not re-log;
log_A is already log10 σ_min).
  0. Fidelity: Ω=1.40170, n=20, xmax=18 must reproduce 2411852's
     -3.297 within 0.05 decades.
  1. Centre Ω=1.40170 at n=20, 28, 36. xmax starts at 18 and grows
     by 8 up to 90 if select_fill is short (same ladder idea as the
     adjudicator; do not force a corrupting fill).
  2. If no n>=28 centre is <= -3.5: 9-point grid on
     [1.3817, 1.4217] at n=28 and n=36 (stays inside the ±5% window
     and well away from 1.226).
  3. If a point at n=28 or 36 reaches -3.5: isolated golden-section
     + two half-shift replicas in a ±0.4% bracket (same isolation
     used for T4 / job 2421061).

PRE-REGISTERED
--------------
  (A) DEPTH_CLEARED: some in-window point at n=28 or n=36 has
      log_A <= -3.5, and the isolated refine is shift-stable
      (<=0.1%) with |Ω_lit-38.428|/38.428 < 2%.
      -> upgrade the 38.428 row from near-miss to refined root.
  (B) STILL_SHALLOW: deepest value at every n stays above -3.5
      (around -3.3 is the expected shape).
      -> leave the near-miss. Optional one-line "confirmed at n=28/36".
  (C) DIP_MOVED: the min walks more than 2% off 38.428, or the n=20
      centre is no longer a dip.
      -> do not edit the paper; show the raw .out.
  (D) GATE_FAIL / FLOOR (log_A <= -12): void / do not treat a floor
      as a better root.

GEOMETRY: identical to 2411852.
  Ri/Ro=0.4, R0_2B=(1+0.4)/(2*(1-0.4)), 2Θ/π=0.5,
  IsotropicMaterial(E=210e9, nu=0.30, rho=7800), FreeFreeOOP,
  M=80, n_quad=30, Ω_lit = 27.4155677808 * Ω_native.
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
N_WORKERS = int(os.environ.get("N_WORKERS", "8"))

EPS = 0.4
R0_2B = (1.0 + EPS) / (2.0 * (1.0 - EPS))
TWO_THETA = 0.5
M = 80
N_QUAD = 30
XMAX0 = 18.0
XMAX_STEP = 8.0
XMAX_CAP = 90.0
FACTOR = 27.4155677808

LIT = 38.428
# Job 2411852 deepest point (not the 2410389 grid-centre 1.40168536).
OM_GATE = 1.40170
DEPTH_REF = -3.297
GATE_TOL = 0.05

WINDOW_LO = 1.3316   # ±5% of predicted centre; excludes Ω=1.22647
WINDOW_HI = 1.4718
GRID_LO = 1.3817
GRID_HI = 1.4217
GRID_N = 9
N_LIST = (20, 28, 36)

DEPTH_THR = -3.5
FLOOR_THR = -12.0
ERR_BAR = 2.0
SHIFT_TOL = 0.1
ISOLATE_HALF = 0.004
GOLDEN_ITERS = 12
INV_PHI = (math.sqrt(5.0) - 1.0) / 2.0


def hdr(s):
    print("\n" + "=" * 78 + f"\n  {s}\n" + "=" * 78, flush=True)


def _build():
    from plate_solver import (OutOfPlaneSolver, IsotropicMaterial,
                              make_geometry, FreeFreeOOP)
    mat = IsotropicMaterial(E=210e9, nu=0.30, rho=7800.0)
    geom = make_geometry(R0_2B, TWO_THETA)
    if abs(float(geom.R_i) / float(geom.R_o) - 0.4) > 1e-6:
        raise RuntimeError("geometry Ri/Ro != 0.4")
    return OutOfPlaneSolver(geom, mat, M=M, n_quad=N_QUAD,
                            boundary=FreeFreeOOP())


def _fill(solver, Om, n_dofs):
    from plate_solver.detectors import full_search, select_fill
    cur = XMAX0
    raw = full_search(solver.fast, Om, xmax=cur)
    sel, cnt = select_fill(raw, n_dofs)
    while cnt < n_dofs and cur < XMAX_CAP:
        cur += XMAX_STEP
        raw = full_search(solver.fast, Om, xmax=cur)
        sel, cnt = select_fill(raw, n_dofs)
    return sel, cnt, cur, len(raw)


def eval_point(solver, Om, n_dofs):
    """Same sigma_min path as 2411852. log_A is already log10."""
    sel, cnt, xmax, n_br = _fill(solver, Om, n_dofs)
    if cnt < max(4, n_dofs // 2):
        return dict(ok=False, log_A=float("nan"), n_ret=cnt,
                    n_br=n_br, xmax=xmax)
    try:
        log_A = float(solver.sigma_min(Om, sel))
    except Exception:
        return dict(ok=False, log_A=float("nan"), n_ret=cnt,
                    n_br=n_br, xmax=xmax)
    return dict(ok=True, log_A=log_A, n_ret=cnt, n_br=n_br, xmax=xmax)


def _worker(args):
    Om, n_dofs = args
    os.environ.setdefault("DPS", "40")
    t0 = time.time()
    out = dict(Om=Om, n_dofs=n_dofs)
    try:
        solver = _build()
        res = eval_point(solver, Om, n_dofs)
        out.update(res)
        out["tag"] = "OK" if res.get("ok") else "EVAL_FAIL"
        if res.get("ok") and res["log_A"] <= FLOOR_THR:
            out["tag"] = "FLOOR"
    except Exception as exc:
        out.update(ok=False, tag="EXC", note=str(exc),
                   log_A=float("nan"))
    out["dt"] = time.time() - t0
    return out


def golden_min(solver, lo, hi, n_dofs, iters):
    a, b = lo, hi
    c = b - INV_PHI * (b - a)
    d = a + INV_PHI * (b - a)
    fc = eval_point(solver, c, n_dofs)["log_A"]
    fd = eval_point(solver, d, n_dofs)["log_A"]
    for _ in range(iters):
        if math.isnan(fc) or math.isnan(fd):
            return float("nan"), float("nan")
        if fc < fd:
            b, d, fd = d, c, fc
            c = b - INV_PHI * (b - a)
            fc = eval_point(solver, c, n_dofs)["log_A"]
        else:
            a, c, fc = c, d, fd
            d = a + INV_PHI * (b - a)
            fd = eval_point(solver, d, n_dofs)["log_A"]
    return (c, fc) if fc < fd else (d, fd)


def main():
    t0 = time.time()
    hdr("probe_oop_shi2014_t3_larger_basis_v1  start "
        + time.strftime("%Y-%m-%d %H:%M:%S"))
    import plate_solver as ps
    print(f"SOLVER_VERSION={getattr(ps,'SOLVER_VERSION','?')} "
          f"(expect {EXPECT_VER})  DPS={os.environ.get('DPS')}  "
          f"N_WORKERS={N_WORKERS}", flush=True)
    if getattr(ps, "SOLVER_VERSION", None) != EXPECT_VER:
        print("PREFLIGHT FAIL: SOLVER_VERSION mismatch")
        raise SystemExit(2)
    print(f"geometry R0_2B={R0_2B:.12g}  Ri/Ro=0.4  2Θ/π={TWO_THETA}",
          flush=True)
    print(f"gate Ω={OM_GATE}  expect log_A={DEPTH_REF} ± {GATE_TOL}  "
          f"(job 2411852)", flush=True)
    print("window stays inside [1.3316, 1.4718]; Ω=1.22647 is out of scope",
          flush=True)

    hdr("FIDELITY  n=20 centre vs job 2411852")
    gate = _worker((OM_GATE, 20))
    print(f"  Ω={OM_GATE:.8f}  log_A={gate.get('log_A', float('nan')):+.4f}  "
          f"n_ret={gate.get('n_ret')}  n_br={gate.get('n_br')}  "
          f"[{gate.get('dt', 0):.0f}s]", flush=True)
    if (not gate.get("ok")
            or abs(gate.get("log_A", 0) - DEPTH_REF) > GATE_TOL):
        print("  GATE_FAIL -- detector has moved; void the run", flush=True)
        raise SystemExit(3)
    print("  GATE PASS", flush=True)

    jobs = [(OM_GATE, n) for n in N_LIST]
    hdr("STAGE 1  centre at n=20/28/36")
    centres = {}
    with ProcessPoolExecutor(max_workers=min(N_WORKERS, 3)) as ex:
        futs = {ex.submit(_worker, j): j for j in jobs}
        for fut in as_completed(futs):
            r = fut.result()
            centres[r["n_dofs"]] = r
            print(f"  n={r['n_dofs']:<2d}  log_A={r.get('log_A', float('nan')):+.4f}  "
                  f"n_ret={r.get('n_ret')}  xmax={r.get('xmax')}  "
                  f"{r.get('tag')}  [{r.get('dt', 0):.0f}s]", flush=True)

    deep_hits = []
    for n, r in centres.items():
        if r.get("ok") and r.get("log_A", 0) <= DEPTH_THR and r.get("tag") != "FLOOR":
            deep_hits.append((OM_GATE, n, r["log_A"]))

    if not any(n >= 28 and r.get("ok") and r.get("log_A", 0) <= DEPTH_THR
               for n, r in centres.items()):
        hdr(f"STAGE 2  {GRID_N}-pt grid [{GRID_LO}, {GRID_HI}] at n=28 and 36")
        grid = [GRID_LO + i * (GRID_HI - GRID_LO) / (GRID_N - 1)
                for i in range(GRID_N)]
        jobs2 = [(om, n) for n in (28, 36) for om in grid]
        rows = []
        nw = max(1, min(N_WORKERS, len(jobs2)))
        with ProcessPoolExecutor(max_workers=nw) as ex:
            futs = {ex.submit(_worker, j): j for j in jobs2}
            done = 0
            for fut in as_completed(futs):
                r = fut.result()
                rows.append(r)
                done += 1
                print(f"  [{done:2d}/{len(jobs2)}] n={r['n_dofs']:<2d}  "
                      f"Ω={r['Om']:.8f}  Ω_lit={r['Om']*FACTOR:8.4f}  "
                      f"log_A={r.get('log_A', float('nan')):+.4f}  "
                      f"n_ret={r.get('n_ret')}  [{r.get('dt', 0):.0f}s]",
                      flush=True)
        for r in rows:
            if (r.get("ok") and r.get("log_A", 0) <= DEPTH_THR
                    and r.get("tag") != "FLOOR"):
                deep_hits.append((r["Om"], r["n_dofs"], r["log_A"]))

    refined = []
    if deep_hits:
        hdr("STAGE 3  isolated refine of points that cleared -3.5")
        solver = _build()
        seen = set()
        for om_c, n, la in sorted(deep_hits, key=lambda t: t[2]):
            key = (round(om_c, 6), n)
            if key in seen:
                continue
            seen.add(key)
            lo = om_c * (1.0 - ISOLATE_HALF)
            hi = om_c * (1.0 + ISOLATE_HALF)
            lo = max(lo, WINDOW_LO)
            hi = min(hi, WINDOW_HI)
            om_star, f_star = golden_min(solver, lo, hi, n, GOLDEN_ITERS)
            reps = []
            for shift in (+0.25 * ISOLATE_HALF, -0.25 * ISOLATE_HALF):
                c2 = om_c * (1.0 + shift)
                l2 = max(c2 * (1.0 - ISOLATE_HALF), WINDOW_LO)
                h2 = min(c2 * (1.0 + ISOLATE_HALF), WINDOW_HI)
                o2, f2 = golden_min(solver, l2, h2, n, GOLDEN_ITERS)
                reps.append((o2, f2))
            moves = [abs(o2 - om_star) / om_star * 100.0
                     for o2, _ in reps if not math.isnan(o2)]
            max_move = max(moves) if moves else float("nan")
            om_lit = om_star * FACTOR if not math.isnan(om_star) else float("nan")
            err = (abs(om_lit - LIT) / LIT * 100.0
                   if not math.isnan(om_lit) else float("nan"))
            print(f"  n={n} seed={om_c:.8f} -> Ω*={om_star:.8f}  "
                  f"Ω_lit={om_lit:.4f}  err={err:.3f}%  "
                  f"log_A={f_star:+.4f}  shift={max_move:.4f}%",
                  flush=True)
            refined.append(dict(n=n, om_star=om_star, log_A=f_star,
                                om_lit=om_lit, err=err, shift=max_move))

    hdr("PRE-REGISTERED READING")
    print("  centres:", flush=True)
    for n in N_LIST:
        r = centres.get(n, {})
        print(f"    n={n}  log_A={r.get('log_A', float('nan')):+.4f}  "
              f"n_ret={r.get('n_ret')}", flush=True)

    deepest = None
    for n, r in centres.items():
        if r.get("ok") and (deepest is None or r["log_A"] < deepest[2]):
            deepest = (OM_GATE, n, r["log_A"])
    if deep_hits:
        for h in deep_hits:
            if deepest is None or h[2] < deepest[2]:
                deepest = h
    print(f"  deepest seen: {deepest}", flush=True)

    n20 = centres.get(20, {})
    centre_ok = (n20.get("ok") and
                 abs(n20.get("log_A", 0) - DEPTH_REF) <= 0.3)
    cleared = [r for r in refined
               if r.get("log_A", 0) <= DEPTH_THR
               and r.get("err", 99) < ERR_BAR
               and (not math.isnan(r.get("shift", float("nan")))
                    and r["shift"] <= SHIFT_TOL)]
    floors = [r for r in centres.values() if r.get("tag") == "FLOOR"]
    floors += [r for r in refined if r.get("log_A", 0) <= FLOOR_THR]

    if floors and not cleared:
        verdict = "D"
        note = "Floor encountered and no stable -3.5 root. Do not edit."
    elif cleared:
        verdict = "A"
        note = ("Depth cleared at larger basis with a stable refine. "
                "Upgrade the 38.428 row to a refined root.")
    elif not centre_ok:
        verdict = "C"
        note = "n=20 centre is no longer the 2411852 dip. Do not edit."
    else:
        verdict = "B"
        note = ("Still shallower than -3.5 at n=20/28/36. Leave the "
                "near-miss; optional clause that the deficit survives "
                "a larger basis.")

    print(f"\n  VERDICT: {verdict}", flush=True)
    print(f"  {note}", flush=True)
    print(f"  wall {time.time()-t0:.1f}s", flush=True)
    hdr("probe_oop_shi2014_t3_larger_basis_v1  end")


if __name__ == "__main__":
    main()
