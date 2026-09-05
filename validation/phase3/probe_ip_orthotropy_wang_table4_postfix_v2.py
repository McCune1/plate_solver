# -*- coding: utf-8 -*-
"""
probe_ip_orthotropy_wang_table4_postfix_v2.py -- RESUME of
probe_ip_orthotropy_wang_table4_postfix_v1.py (job 2334099), which hit its
8h SLURM time limit partway through target 7's search (confirmed: the sum
of targets 1-6's own printed wall times is ~26502s = 7.36h, consistent with
target 7 starting and then being cut off by the 28800s/8h limit -- not a
computational failure).

Targets 1-6 are ALREADY CONFIRMED HIT under the post-fix (mu_theta, worker-
path-fixed) code and are hardcoded below rather than re-run (each cost
~1-1.7h; re-running them would just burn cluster time for no new
information):

  1 (Wang=1.0311): found=2.3070200248963486  rel%=0.061  HIT
  2 (Wang=1.7367): found=3.8916626377044565  rel%=0.213  HIT
  3 (Wang=2.0502): found=4.550587725467025   rel%=0.737  HIT
  4 (Wang=3.0618): found=6.809885435424392   rel%=0.533  HIT
  5 (Wang=3.1811): found=7.108489609037644   rel%=0.066  HIT
  6 (Wang=3.4093): found=7.624496535393413   rel%=0.014  HIT

This run only searches targets 7-8 (both fix-presence assertions re-checked
first, cheaply, in case this is a different checkout than job 2334099's).

Pre-registered reading: unchanged from v1 -- HIT if nearest accepted
singularity within |Omega-target|<0.15; overall PASS if hits (6 already
banked + hits here) >= 6 of 8 total (already guaranteed regardless of 7/8's
outcome, but report the full 8-target table for the record).
"""
import os, sys, time
sys.path.insert(0, os.environ.get("PKG_PATH", "."))
os.environ.setdefault("DPS", "40")

from mpmath import mp
mp.dps = int(os.environ["DPS"])

import plate_solver as ps
from plate_solver import (
    OrthotropicMaterial, make_geometry, FreeFreeIP, InPlaneSolver,
    find_modes_sigmin,
)

N_WORKERS = int(os.environ.get("N_WORKERS", os.environ.get("SLURM_CPUS_PER_TASK", "32")))

WANG_FFFF = [1.0311, 1.7367, 2.0502, 3.0618, 3.1811, 3.4093, 4.3023, 4.5752]
K_CONVERT = 1.0 / float(mp.sqrt(5))
TARGETS = [w / K_CONVERT for w in WANG_FFFF]

# Banked from job 2334099 (targets 1-6, all HIT)
BANKED = {
    0: 2.3070200248963486,
    1: 3.8916626377044565,
    2: 4.550587725467025,
    3: 6.809885435424392,
    4: 7.108489609037644,
    5: 7.624496535393413,
}

HALF_WIDTH = 0.20
N_DOFS = 14
MAX_DIM = 18.0
N_SCAN = 50
N_MODES_LOCAL = 4


def hdr(s):
    print("\n" + "=" * 78 + f"\n  {s}\n" + "=" * 78, flush=True)


def make_wang_material():
    Et = mp.mpf("70e9"); Er = 40 * Et; nu_r = mp.mpf("0.3"); rho = mp.mpf("7850")
    nu_theta = nu_r * Et / Er
    G = Et / (2 * (1 - nu_r * nu_theta))
    return OrthotropicMaterial(float(Er), float(Et), float(nu_r), float(G), float(rho))


def local_search(solver, lo, hi):
    freqs = find_modes_sigmin(
        solver, part=2,
        Omega_range=(lo, hi),
        n_scan=N_SCAN, n_dofs=N_DOFS, max_dim=MAX_DIM,
        n_modes_wanted=N_MODES_LOCAL,
        verbose=True, n_workers=N_WORKERS,
    )
    return [float(f) for f in freqs]


def main():
    t0 = time.time()
    print(f"probe_ip_orthotropy_wang_table4_postfix_v2 (resume of job 2334099)  "
          f"N_WORKERS={N_WORKERS}", flush=True)
    print(f"DPS={mp.dps}  SOLVER_VERSION={getattr(ps,'SOLVER_VERSION','?')}", flush=True)
    if not str(getattr(ps, 'SOLVER_VERSION', '')).endswith('s9'):
        print("  WARNING: SOLVER_VERSION does not end in s9 -- check this "
              "checkout has both IP orthotropy fixes before trusting results.",
              flush=True)

    hdr("Material + solver setup (NO monkeypatch -- fix is in the package)")
    geo = make_geometry(r0_over_2b=1.5, two_Theta_over_pi=1.0, h=0.001, b=0.25)
    mat = make_wang_material()
    ip_c = mat.ip_constants()
    print(f"  ip_constants() = {ip_c}", flush=True)
    print(f"  mat.mu_theta   = {mat.mu_theta}", flush=True)
    print(f"  mat.nu_bar     = {mat.nu_bar}", flush=True)
    assert abs(float(ip_c[2]) - float(mat.mu_theta)) < 1e-12, \
        "ip_constants()[2] != mu_theta -- fix 1 (geometry.py) not present, STOP"
    assert abs(float(ip_c[2]) - float(mat.nu_bar)) > 1e-3, \
        "ip_constants()[2] == nu_bar -- fix did not change anything, STOP"

    solver = InPlaneSolver(geo, mat, M=80, n_quad=30, boundary=FreeFreeIP())
    g_sc, m_sc = solver._p2_worker_params()
    print(f"  _p2_worker_params() m_sc[0] (worker's nu_bar slot) = {m_sc[0]}", flush=True)
    assert abs(m_sc[0] - float(mat.mu_theta)) < 1e-9, \
        "worker tuple's nu_bar slot != mu_theta -- fix 2 (core_solvers.py) not present, STOP"
    print("  Both fixes confirmed present and wired through the worker "
          "tuple. Proceeding to targets 7-8 only.", flush=True)

    results = []
    for i in range(6):
        tgt = TARGETS[i]
        near = BANKED[i]
        ae = abs(near - tgt)
        re = ae / tgt * 100
        results.append((i + 1, WANG_FFFF[i], tgt, near, ae, re, True, [near]))
        print(f"  BANKED target {i+1} (job 2334099): found={near:.6f}  "
              f"rel%={re:.3f}  HIT", flush=True)

    for i in range(6, 8):
        tgt = TARGETS[i]
        lo = max(0.5, tgt - HALF_WIDTH)
        hi = tgt + HALF_WIDTH
        hdr(f"TARGET {i+1}: Wang={WANG_FFFF[i]:.4f}  ->  Omega={tgt:.4f}  "
            f"window=[{lo:.3f},{hi:.3f}]")
        t1 = time.time()
        try:
            found = local_search(solver, lo, hi)
        except Exception as e:
            print(f"  SEARCH FAILED: {e}", flush=True)
            import traceback; traceback.print_exc()
            found = []
        elapsed = time.time() - t1
        print(f"  found: {found}  ({elapsed:.1f} s)", flush=True)

        if found:
            nearest = min(found, key=lambda f: abs(f - tgt))
            abs_e = abs(nearest - tgt)
            rel_e = abs_e / tgt * 100
            hit = abs_e < 0.15
        else:
            nearest, abs_e, rel_e, hit = None, float("nan"), float("nan"), False
        results.append((i + 1, WANG_FFFF[i], tgt, nearest, abs_e, rel_e, hit, found))
        print(f"  nearest={nearest}  abs_err={abs_e:.4e}  rel%={rel_e:.3f}  "
              f"{'HIT' if hit else 'MISS'}", flush=True)

    hdr("SUMMARY TABLE (all 8, targets 1-6 banked from job 2334099)")
    print(f"{'#':>3}  {'Wang':>8}  {'target':>10}  {'found':>10}  "
          f"{'abs_err':>10}  {'rel%':>8}  {'hit':>4}", flush=True)
    print("-" * 64, flush=True)
    n_hit = 0
    for i, wang, tgt, near, ae, re, hit, _ in results:
        near_s = f"{near:.6f}" if near is not None else "  ---"
        print(f"{i:3d}  {wang:8.4f}  {tgt:10.4f}  {near_s:>10}  "
              f"{ae:10.4e}  {re:8.3f}  {'YES' if hit else ' no'}", flush=True)
        if hit:
            n_hit += 1
    print("-" * 64, flush=True)
    print(f"  hits: {n_hit} / 8", flush=True)

    hdr("READING")
    ok = n_hit >= 6
    print(f"  Literature match (post worker-path fix): {'PASS' if ok else 'FAIL'}  ({n_hit}/8)", flush=True)
    print(f"  Wall time (this run, targets 7-8 only): {time.time()-t0:.1f} s", flush=True)
    if not ok:
        print("  DO NOT re-tune K_CONVERT or the tolerance to force a match.", flush=True)
        print("  Report the raw per-target errors for further investigation.", flush=True)
        sys.exit(1)
    print("\nALL CHECKS PASSED -- mu_theta fix genuinely validated through the worker path", flush=True)
    print("(6/6 banked from job 2334099 + this run's targets 7-8)", flush=True)


if __name__ == "__main__":
    main()
