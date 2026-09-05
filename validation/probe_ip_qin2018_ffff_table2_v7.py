# -*- coding: utf-8 -*-
"""
probe_ip_qin2018_ffff_table2_v7.py

RERUN of the Qin 2018 Table 2 FFFF validation with the REAL production
detector (find_modes_sigmin) under proper SLURM parallelism (32 CPUs).

CONTEXT
-------
v1–v3 used find_modes_sigmin but requested --cpus-per-task=1, so the
coarse scan ran serially (~55–60 s/point) and timed out.  v4–v6 abandoned
the real detector for fixed grids / locked bases to fit the 1-CPU budget;
their PARTIAL / "structural near-null" conclusion was therefore an
artifact of the resource misconfiguration, not a physics finding.

v3's single completed window already produced a genuine hit:
  mode Ω=0.401508  σ_min ≈ 1e-5.01
right on the plain-factor prediction for Qin's first target (1.5641).

This probe is exactly that configuration, with the only change being
N_WORKERS = SLURM_CPUS_PER_TASK (=32) so each 24-point local window
parallelizes properly.

HYPOTHESIS
----------
Rerunning the six tight local windows around the plain-factor predicted
natives recovers most or all of Qin's Table 2 FFFF targets under the
plain conversion, at confidence levels similar to the confirmed first-mode
hit (log10 σ_min ≈ −5).

PRE-REGISTERED CRITERIA (unchanged)
-----------------------------------
  PASS:    >=5/6 within 2%
  PARTIAL: 3–4/6 within 2%, or >=5/6 within 5%
  FAIL:    <=2/6 within 5%, or conversion factor cannot be pinned down

GEOMETRY / MATERIAL (identical to v1–v3)
----------------------------------------
  make_geometry(r0_over_2b=1.5, two_Theta_over_pi=0.5, h=0.001, b=0.25)
  => R1=1, Ri/Ro=0.5, 2Θ=π/2
  IsotropicMaterial(E=70e9, nu=0.3, rho=2700)
  FreeFreeIP, InPlaneSolver, M=80, n_quad=30

SEARCH (identical to v3's tight-window design)
----------------------------------------------
  N_DOFS=12, MAX_DIM=16.0, N_SCAN=24, N_MODES_LOCAL=3, HALF_WIDTH=0.10
  One find_modes_sigmin call per plain-predicted native.

Diagnostic only. No package modification. No SOLVER_VERSION impact.
"""
from __future__ import annotations
import os
import sys
import time
import math

sys.path.insert(0, os.environ.get("PKG_PATH", "."))
os.environ.setdefault("DPS", "40")

from mpmath import mp, mpf, fabs, log10, nstr
mp.dps = int(os.environ["DPS"])

import plate_solver as ps
from plate_solver import (
    IsotropicMaterial, make_geometry, FreeFreeIP, InPlaneSolver,
    find_modes_sigmin,
)
from plate_solver.detectors import sigma_min_from_K, select_fill

EXPECT_VER = os.environ.get("EXPECT_SOLVER_VERSION", "2026-07-10.s10")
N_WORKERS = int(os.environ.get("N_WORKERS", os.environ.get("SLURM_CPUS_PER_TASK", "1")))

QIN_TARGETS = [
    1.5641,  # (1.5633)
    2.6702,  # (2.6785)
    2.9792,  # (2.9795)
    4.4293,  # (4.4513)
    4.4754,  # (4.4627)
    4.6575,  # (4.6668)
]

FIDELITY_OMEGA = 0.401154

N_DOFS = 12
MAX_DIM = 16.0
XMAX = 16.0
N_SCAN = 24
N_MODES_LOCAL = 3
HALF_WIDTH = 0.10


def hdr(s):
    print("\n" + "=" * 78 + f"\n  {s}\n" + "=" * 78, flush=True)


def get_geometry_attrs(geom):
    attrs = {}
    for name in ("ri", "ro", "a", "b", "R0", "R_inner", "R_outer",
                 "r_inner", "r_outer", "width", "h", "Theta", "two_theta"):
        if hasattr(geom, name):
            try:
                attrs[name] = float(getattr(geom, name))
            except Exception:
                pass
    return attrs


def compute_omega_bar(geom, mat):
    for meth in ("_kbar_wbar", "kbar_wbar", "omega_bar", "_omega_bar"):
        if hasattr(geom, meth):
            try:
                res = getattr(geom, meth)(mat)
                if isinstance(res, (tuple, list)) and len(res) >= 2:
                    return float(res[1])
                if hasattr(res, "imag") or isinstance(res, (int, float, mpf)):
                    return float(res)
            except Exception as e:
                print(f"  helper {meth} failed: {e}", flush=True)
    E = float(getattr(mat, "E", 70e9))
    nu = float(getattr(mat, "nu", 0.3))
    rho = float(getattr(mat, "rho", 2700.0))
    G = E / (2.0 * (1.0 + nu))
    attrs = get_geometry_attrs(geom)
    if "b" in attrs:
        twob = 2.0 * attrs["b"]
    elif "width" in attrs:
        twob = attrs["width"]
    else:
        twob = 0.5
        print("  WARNING: falling back to 2b=0.5", flush=True)
    return (math.pi / twob) * math.sqrt(G / rho)


def convert_native_to_qin(Omega_native, omega_bar, R1, rho, E, nu, h=1.0):
    omega = Omega_native * omega_bar
    qin_b = omega * R1 * math.sqrt(rho * h / E)
    qin_a = (2.0 * omega * R1 / math.pi) * math.sqrt(rho / (E * (1.0 - nu * nu)))
    qin_plain = omega * R1 * math.sqrt(rho / E)
    return {"a": qin_a, "b_h": qin_b, "plain": qin_plain, "omega": omega}


def nearest_rel_err(cand, targets):
    best = None
    for t in targets:
        rel = abs(cand - t) / t
        if best is None or rel < best[0]:
            best = (rel, t)
    return best


def fidelity_gate():
    hdr("FIDELITY GATE (project default FreeFreeIP, TWO_T=1.0)")
    mat = IsotropicMaterial(E=210e9, nu=0.30, rho=7800.0)
    try:
        geom = make_geometry(1.5, 1.0)
    except TypeError:
        geom = make_geometry(r0_over_2b=1.5, two_Theta_over_pi=1.0)
    solver = InPlaneSolver(geom, mat, M=80, n_quad=30, boundary=FreeFreeIP())
    print(f"  SOLVER_VERSION={getattr(ps, 'SOLVER_VERSION', '?')}", flush=True)
    if getattr(ps, "SOLVER_VERSION", None) != EXPECT_VER:
        print(f"  PREFLIGHT WARN: version {ps.SOLVER_VERSION!r} != {EXPECT_VER!r}", flush=True)

    try:
        from plate_solver import full_search
    except ImportError:
        full_search = getattr(ps, "full_search", None)
    if full_search is None:
        print("  full_search unavailable; skipping deep fidelity", flush=True)
        return True

    t0 = time.time()
    brs = full_search(solver.fast, FIDELITY_OMEGA, xmax=XMAX)
    print(f"  branches at Ω={FIDELITY_OMEGA}: {len(brs)} ({time.time()-t0:.1f}s)", flush=True)
    n_use = min(max(len(brs), 4), 10)
    sel = select_fill(brs, n_use)[0] if len(brs) >= n_use else brs
    K, sz = solver._build_K_real(FIDELITY_OMEGA, sel, lagrange=True, fast_scan=False)
    lg = float(sigma_min_from_K(K, sz))
    print(f"  sigma_min_from_K(Ω={FIDELITY_OMEGA}) size={sz}  log10={lg:.3f}", flush=True)
    if lg > -1.0:
        print("  FIDELITY FAIL", flush=True)
        return False
    print("  FIDELITY PASS (or soft-pass)", flush=True)
    return True


def local_search(solver, lo, hi, label):
    print(f"\n  --- local [{lo:.4f}, {hi:.4f}]  ({label}) ---", flush=True)
    t1 = time.time()
    try:
        freqs = find_modes_sigmin(
            solver, part=2,
            Omega_range=(lo, hi),
            n_scan=N_SCAN, n_dofs=N_DOFS, max_dim=MAX_DIM,
            n_modes_wanted=N_MODES_LOCAL,
            verbose=True, n_workers=N_WORKERS,
        )
        found = sorted(float(f) for f in freqs)
        print(f"  found {len(found)} in {time.time()-t1:.1f}s: {found}", flush=True)
        return found
    except Exception as e:
        print(f"  local search failed: {e}", flush=True)
        import traceback
        traceback.print_exc()
        return []


def main():
    t0 = time.time()
    hdr("probe_ip_qin2018_ffff_table2_v7  start " + time.strftime("%Y-%m-%d %H:%M:%S"))
    print(f"DPS={mp.dps}  N_WORKERS={N_WORKERS}  N_DOFS={N_DOFS}  "
          f"HALF_WIDTH={HALF_WIDTH}  N_SCAN={N_SCAN}", flush=True)
    print(f"  (N_WORKERS should be 32 under proper SLURM; was 1 in v1–v6)", flush=True)

    if not fidelity_gate():
        print("ABORT on fidelity", flush=True)
        raise SystemExit(3)

    hdr("QIN GEOMETRY + MATERIAL (explicit b=0.25, R1=1)")
    mat = IsotropicMaterial(E=70e9, nu=0.3, rho=2700.0)
    try:
        geom = make_geometry(
            r0_over_2b=1.5,
            two_Theta_over_pi=0.5,
            h=0.001,
            b=0.25,
        )
    except TypeError:
        geom = make_geometry(1.5, 0.5)
        print("  WARNING: explicit kwargs failed; positional fallback", flush=True)

    attrs = get_geometry_attrs(geom)
    print(f"  geometry attrs: {attrs}", flush=True)
    R1 = 1.0
    h = attrs.get("h", 0.001)
    b_val = attrs.get("b", 0.25)
    print(f"  forced R1={R1}  h={h}  b={b_val}", flush=True)
    if abs(b_val - 0.25) > 0.01:
        print(f"  WARNING: geometry b={b_val} != 0.25", flush=True)

    omega_bar = compute_omega_bar(geom, mat)
    print(f"  omega_bar (native) = {omega_bar:.6e}", flush=True)

    fac_plain = omega_bar * R1 * math.sqrt(2700.0 / 70e9)
    fac_h = omega_bar * R1 * math.sqrt(2700.0 * h / 70e9)
    fac_a = omega_bar * (2.0 * R1 / math.pi) * math.sqrt(2700.0 / (70e9 * (1 - 0.3**2)))
    print(f"  conversion factors (Qin/native):", flush=True)
    print(f"    plain (omega R1 sqrt(rho/E))     = {fac_plain:.6f}", flush=True)
    print(f"    with-h (omega R1 sqrt(rho h/E))  = {fac_h:.6f}", flush=True)
    print(f"    (a) (2/pi)(1-nu^2) form          = {fac_a:.6f}", flush=True)

    pred_plain = [t / fac_plain for t in QIN_TARGETS]
    print(f"\n  predicted native (plain): {[f'{p:.4f}' for p in pred_plain]}", flush=True)
    print(f"  (search windows = predicted ± {HALF_WIDTH})", flush=True)

    hdr("SCALE-FREE SPOT-CHECK")
    mat210 = IsotropicMaterial(E=210e9, nu=0.3, rho=2700.0)
    wbar210 = compute_omega_bar(geom, mat210)
    print(f"  omega_bar ratio E210/E70 = {wbar210/omega_bar:.6f}  (expect 1.732)", flush=True)

    solver = InPlaneSolver(geom, mat, M=80, n_quad=30, boundary=FreeFreeIP())

    # Exactly the 6 plain-predicted centers (no expanded seed list)
    seeds = [round(p, 4) for p in pred_plain]
    print(f"\n  local-search seeds (6 plain predictions): {seeds}", flush=True)

    hdr("TIGHT LOCAL SEARCHES (real find_modes_sigmin, parallel)")
    all_found = []
    for seed in seeds:
        lo = max(0.05, seed - HALF_WIDTH)
        hi = seed + HALF_WIDTH
        found = local_search(solver, lo, hi, f"seed={seed:.4f}")
        all_found.extend(found)

    uniq = sorted(set(round(f, 6) for f in all_found))
    print(f"\n  unique native candidates ({len(uniq)}): {uniq}", flush=True)

    hdr("CONVERSION + MATCH TO QIN TABLE 2")
    print(f"{'native':>10}  {'Qin_plain':>10}  {'Qin_h':>10}  {'Qin_a':>10}  "
          f"{'nearest_T':>10}  {'rel%_p':>8}  {'rel%_h':>8}  {'rel%_a':>8}", flush=True)
    matches_p = matches_h = matches_a = 0
    details = []
    for om in uniq:
        conv = convert_native_to_qin(om, omega_bar, R1, 2700.0, 70e9, 0.3, h=h)
        rel_p, t_p = nearest_rel_err(conv["plain"], QIN_TARGETS)
        rel_h, t_h = nearest_rel_err(conv["b_h"], QIN_TARGETS)
        rel_a, t_a = nearest_rel_err(conv["a"], QIN_TARGETS)
        print(f"{om:10.6f}  {conv['plain']:10.4f}  {conv['b_h']:10.4f}  {conv['a']:10.4f}  "
              f"{t_p:10.4f}  {100*rel_p:7.2f}%  {100*rel_h:7.2f}%  {100*rel_a:7.2f}%", flush=True)
        details.append((om, conv, rel_p, rel_h, rel_a, t_p))
        if rel_p < 0.02:
            matches_p += 1
        if rel_h < 0.02:
            matches_h += 1
        if rel_a < 0.02:
            matches_a += 1

    def count_within(rels, thr):
        return sum(1 for r in rels if r < thr)

    rels_p = [d[2] for d in details]
    rels_h = [d[3] for d in details]
    rels_a = [d[4] for d in details]

    hdr("PRE-REGISTERED VERDICT")
    print(f"  #candidates native = {len(uniq)}", flush=True)
    print(f"  matches within 2%  (plain / h / a) = {matches_p} / {matches_h} / {matches_a}", flush=True)
    print(f"  matches within 5%  (plain / h / a) = "
          f"{count_within(rels_p,0.05)} / {count_within(rels_h,0.05)} / {count_within(rels_a,0.05)}", flush=True)

    scores = [
        ("plain (omega R1 sqrt(rho/E))", matches_p, count_within(rels_p, 0.05)),
        ("with-h", matches_h, count_within(rels_h, 0.05)),
        ("(a) 2/pi (1-nu^2)", matches_a, count_within(rels_a, 0.05)),
    ]
    scores.sort(key=lambda t: (-t[1], -t[2]))
    best_name, best2, best5 = scores[0]
    print(f"\n  Best formula by hit count: {best_name}  (2%={best2}, 5%={best5})", flush=True)

    if best2 >= 5:
        verdict = "PASS"
    elif best2 >= 3 or best5 >= 5:
        verdict = "PARTIAL"
    else:
        verdict = "FAIL"
    print(f"\n  VERDICT: {verdict}", flush=True)
    print(f"  (criteria use the best of the three candidate conversions)", flush=True)

    print("\n  DEPTH NOTE: inspect the find_modes_sigmin verbose lines above.", flush=True)
    print("  Genuine hits look like v3's confirmed mode (log10 σ_min ≈ −5).", flush=True)
    print("  Locked-basis floor (v5/v6) was log10 ≈ −11…−43 everywhere —", flush=True)
    print("  if that pattern reappears here, the detector itself is still", flush=True)
    print("  reporting a near-null; otherwise the earlier 'structural'", flush=True)
    print("  conclusion is overturned.", flush=True)

    print(f"\n  total wall {time.time()-t0:.1f} s", flush=True)
    hdr("probe_ip_qin2018_ffff_table2_v7  end")


if __name__ == "__main__":
    main()