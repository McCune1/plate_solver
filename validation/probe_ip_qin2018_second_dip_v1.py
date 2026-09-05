# -*- coding: utf-8 -*-
"""
probe_ip_qin2018_second_dip_v1.py -- DIAGNOSTIC ONLY.
No package changes, no SOLVER_VERSION impact.

Rank 5 leftover. Jobs 2421062/2421063 closed B STILL SHARED: the
production detector returned one root (Ω=1.14274, Qin 4.4529). The
dense scan also showed a second local min at Ω=1.147551 (Qin 4.4716,
n=20 log_A=-3.653) that find_modes_sigmin did not accept. That second
dip was outside pre-registered A (detector split). This job only asks
whether that dip is an isolated, shift-stable root.

WHAT THIS IS NOT
----------------
It does not re-run find_modes_sigmin on the pair window. Even a STABLE
second dip does NOT change tab:qin to two roots — the production
detector still reports one. The most the paper can gain is a clause
that a second dense-scan dip exists and is/isn't isolated.

METHOD
------
Same Qin geometry as v7 / 2421062: r0/2b=1.5, 2Θ/π=0.5, h=0.001, b=0.25,
E=70e9, ν=0.3, ρ=2700, FreeFreeIP, xmax=60, n=20 then n=28.
  0. Fidelity: Ω=1.14274 n=20 must be log_A <= -3.0 (2421062: -3.531
     at the grid point; detector root 1.14274).
  1. Isolated ±0.15% golden-section of 1.147551 at n=20 and n=28,
     plus two half-shift replicas. Window cannot reach 1.14274
     (they are 0.43% apart).
  2. Same isolation of the 1.14274 control, to confirm the method
     recovers the known root.

PRE-REGISTERED
--------------
  (A) STABLE_SECOND: 1.147551 refine is shift-stable (<=0.1%),
      log_A <= -3.5, and stays >0.3% from 1.14274 at both n=20 and
      n=28. -> optional §6.8 clause. Do NOT change the shared-root
      table row.
  (B) SHOULDER: refine walks toward 1.14274, is unstable, or stays
      shallower than -3.5. -> leave the paper as written.
  (C) GATE_FAIL / FLOOR: void.

COST: ~40 evals. Budget 3 h / 8 cpus.
"""
from __future__ import annotations
import os
import sys
import time
import math

sys.path.insert(0, os.environ.get("PKG_PATH", "."))
os.environ.setdefault("DPS", "40")

from mpmath import mp
mp.dps = int(os.environ["DPS"])

EXPECT_VER = os.environ.get("EXPECT_SOLVER_VERSION", "2026-07-10.s10")
N_WORKERS = int(os.environ.get("N_WORKERS",
                               os.environ.get("SLURM_CPUS_PER_TASK", "8")))

M = 80
N_QUAD = 30
XMAX = 60.0
OM_DET = 1.14274           # production root, job 2421062
OM_DIP = 1.147551          # dense-scan second min
DEPTH_DET_REF = -3.531
GATE_LOOSE = -3.0
DEPTH_THR = -3.5
FLOOR_THR = -12.0
SHIFT_TOL = 0.1
SEP_BAR = 0.003
ISOLATE_HALF = 0.0015      # ±0.15%; cannot reach the 0.43% neighbor
GOLDEN_ITERS = 12
INV_PHI = (math.sqrt(5.0) - 1.0) / 2.0
FAC_PLAIN = 3.896666


def hdr(s):
    print("\n" + "=" * 78 + f"\n  {s}\n" + "=" * 78, flush=True)


def _build():
    from plate_solver import (InPlaneSolver, IsotropicMaterial,
                              make_geometry, FreeFreeIP)
    mat = IsotropicMaterial(E=70e9, nu=0.3, rho=2700.0)
    try:
        geom = make_geometry(r0_over_2b=1.5, two_Theta_over_pi=0.5,
                             h=0.001, b=0.25)
    except TypeError:
        geom = make_geometry(1.5, 0.5)
    return InPlaneSolver(geom, mat, M=M, n_quad=N_QUAD, boundary=FreeFreeIP())


def eval_point(solver, Om, n_dofs):
    from plate_solver.detectors import full_search, select_fill
    raw = full_search(solver.fast, Om, xmax=XMAX)
    sel, cnt = select_fill(raw, n_dofs)
    if cnt < max(4, n_dofs // 2):
        return dict(ok=False, log_A=float("nan"), n_ret=cnt, n_br=len(raw))
    try:
        log_A = float(solver.sigma_min(Om, sel))
    except Exception:
        return dict(ok=False, log_A=float("nan"), n_ret=cnt, n_br=len(raw))
    return dict(ok=True, log_A=log_A, n_ret=cnt, n_br=len(raw))


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


def isolate(solver, om_c, n_dofs):
    lo = om_c * (1.0 - ISOLATE_HALF)
    hi = om_c * (1.0 + ISOLATE_HALF)
    om_star, f_star = golden_min(solver, lo, hi, n_dofs, GOLDEN_ITERS)
    reps = []
    for shift in (+0.25 * ISOLATE_HALF, -0.25 * ISOLATE_HALF):
        c2 = om_c * (1.0 + shift)
        l2 = c2 * (1.0 - ISOLATE_HALF)
        h2 = c2 * (1.0 + ISOLATE_HALF)
        o2, f2 = golden_min(solver, l2, h2, n_dofs, GOLDEN_ITERS)
        reps.append((o2, f2))
    moves = [abs(o2 - om_star) / om_star * 100.0
             for o2, _ in reps if not math.isnan(o2)]
    max_move = max(moves) if moves else float("nan")
    return om_star, f_star, max_move


def main():
    t0 = time.time()
    hdr("probe_ip_qin2018_second_dip_v1  start "
        + time.strftime("%Y-%m-%d %H:%M:%S"))
    import plate_solver as ps
    print(f"SOLVER_VERSION={getattr(ps,'SOLVER_VERSION','?')} "
          f"(expect {EXPECT_VER})", flush=True)
    if getattr(ps, "SOLVER_VERSION", None) != EXPECT_VER:
        print("PREFLIGHT FAIL")
        raise SystemExit(2)

    solver = _build()
    b = float(getattr(solver.geom, "b", -1.0)) if hasattr(solver, "geom") else -1.0
    print(f"  Qin pair leftover dip {OM_DIP} vs detector {OM_DET}", flush=True)

    hdr("FIDELITY  detector root 1.14274 n=20")
    g = eval_point(solver, OM_DET, 20)
    print(f"  log_A={g.get('log_A', float('nan')):+.4f}  n_ret={g.get('n_ret')}  "
          f"(2421062 grid -3.531)", flush=True)
    if (not g.get("ok")) or g.get("log_A", 0) > GATE_LOOSE:
        print("  GATE_FAIL", flush=True)
        raise SystemExit(3)
    print("  GATE PASS", flush=True)

    hdr("ISOLATE both dips at n=20 and n=28")
    results = {}
    for n in (20, 28):
        for name, seed in (("DET", OM_DET), ("DIP", OM_DIP)):
            om, la, mv = isolate(solver, seed, n)
            sep = abs(om - OM_DET) / OM_DET * 100.0
            qin = om * FAC_PLAIN
            print(f"  {name} n={n} seed={seed:.8f} -> Ω*={om:.8f}  "
                  f"Qin={qin:.4f}  log_A={la:+.4f}  shift={mv:.4f}%  "
                  f"sep_det={sep:.3f}%", flush=True)
            results[(name, n)] = dict(om=om, log_A=la, shift=mv, sep=sep)

    hdr("PRE-REGISTERED READING")
    dips = [results[("DIP", n)] for n in (20, 28)]
    floors = [r for r in dips if r.get("log_A", 0) <= FLOOR_THR]
    stable = [r for r in dips
              if r.get("log_A", 0) <= DEPTH_THR
              and (not math.isnan(r.get("shift", float("nan")))
                   and r["shift"] <= SHIFT_TOL)
              and r.get("sep", 0) > SEP_BAR * 100.0]

    if floors:
        verdict = "C"
        note = "Floor on the second dip. Void."
    elif len(stable) == 2:
        verdict = "A"
        note = ("Second dip is isolated and shift-stable at n=20 and n=28. "
                "Optional §6.8 clause only. Do NOT change tab:qin to two "
                "roots — find_modes_sigmin still returned one.")
    else:
        verdict = "B"
        note = ("Second dip is a shoulder, unstable, or shallow. Leave "
                "§6.8 as written.")

    print(f"\n  VERDICT: {verdict}", flush=True)
    print(f"  {note}", flush=True)
    print(f"  wall {time.time()-t0:.1f}s", flush=True)
    hdr("probe_ip_qin2018_second_dip_v1  end")


if __name__ == "__main__":
    main()
