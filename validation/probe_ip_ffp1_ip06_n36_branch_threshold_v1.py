# -*- coding: utf-8 -*-
"""
probe_ip_ffp1_ip06_n36_branch_threshold_v1.py

HYPOTHESIS
----------
IP-06's n_dofs=36 branch-fill discrepancy (30 branches retained in job
2324822 vs 34 in job 2324817, same nominal target) is caused by
run-to-run sensitivity of the converged Ω landing near a branch-count
threshold boundary, NOT true nondeterminism in select_fill / full_search
at a fixed Ω.

PRE-REGISTERED CRITERIA
-----------------------
PASS (confirms threshold-sensitivity):
  - At the FIXED, already-confirmed root Ω=1.19704984,
    full_search + select_fill(n_dofs=36) returns the IDENTICAL cnt
    across >=5 repeated calls in the same process
    (rules out true per-call randomness)
  - AND a dense scan across Ω in [1.1960, 1.1990] (step 0.0002)
    shows cnt taking at least two different values across that
    narrow window -- i.e. a real boundary effect in Ω.

FAIL / kill:
  - cnt varies across repeated calls AT THE SAME Ω
    (more serious nondeterminism; different follow-up needed)
  - OR cnt is constant across the whole scanned window
    (rules out simple threshold sensitivity; the 30-vs-34 gap
     between the two jobs would then need a different explanation,
     e.g. environment / package-version difference this probe
     cannot see).

GEOMETRY / BC / BASIS (exact match to hand-off)
-----------------------------------------------
  make_geometry(1.5, 0.5)          # FF-P1 primary, NOT geomsweep
  IsotropicMaterial(E=210e9, nu=0.30, rho=7800.0)
  InPlaneSolver + FreeFreeIP
  n_dofs=36, xmax=32.0 (confirmed max_dim pairing), M=80, n_quad=30, dps=40

CONTEXT
-------
IP-06's own converged root at n_dofs=28 is Ω=1.19704984, log10σ=-7.78
(jobs 2410035/2410044). That question is CLOSED. This probe only
asks why the branch COUNT at n_dofs=36 differed between two runs of
the same target.

Diagnostic only. No package modification. No SOLVER_VERSION impact.
"""
from __future__ import annotations
import os
import sys
import time

sys.path.insert(0, os.environ.get("PKG_PATH", "."))
os.environ.setdefault("DPS", "40")

from mpmath import mp
mp.dps = int(os.environ["DPS"])

import plate_solver as ps
from plate_solver import detectors as det
from plate_solver.detectors import select_fill

try:
    from plate_solver import full_search
except ImportError:
    full_search = getattr(ps, "full_search", None) or getattr(det, "full_search", None)

EXPECT_VER = os.environ.get("EXPECT_SOLVER_VERSION", "2026-07-10.s10")

R0_2B = 1.5
TWO_T_PI = 0.5
N_DOFS = 36
XMAX = 32.0
FIXED_OM = 1.19704984
OM_LO = 1.1960
OM_HI = 1.1990
OM_STEP = 0.0002
N_REPEATS = 5


def hdr(s):
    print("\n" + "=" * 78 + f"\n  {s}\n" + "=" * 78, flush=True)


def eval_point(r0_2b, two_T_pi, n_dofs, xmax, Om):
    """Exact formula supplied in the hand-off."""
    mat = ps.IsotropicMaterial(E=210e9, nu=0.30, rho=7800.0)
    geom = ps.make_geometry(r0_2b, two_T_pi)
    solver = ps.InPlaneSolver(geom, mat, M=80, n_quad=30,
                              boundary=ps.FreeFreeIP())
    raw = full_search(solver.fast, Om, xmax=xmax)
    sel, cnt = select_fill(raw, n_dofs)
    return cnt, sel, len(raw)


def main():
    t0 = time.time()
    print(f"job start  host={os.uname().nodename}  dps={mp.dps}", flush=True)
    print(f"EXPECT_SOLVER_VERSION={EXPECT_VER}", flush=True)

    # ---- version gate ----
    ver = getattr(ps, "SOLVER_VERSION", None) or getattr(ps, "__version__", "?")
    print(f"package SOLVER_VERSION = {ver!r}", flush=True)
    if ver != EXPECT_VER:
        print(f"FATAL: version mismatch (got {ver!r}, expected {EXPECT_VER!r})",
              flush=True)
        sys.exit(3)
    print("version gate PASS", flush=True)

    if full_search is None:
        print("FATAL: full_search not importable", flush=True)
        sys.exit(3)

    # ================================================================
    # PART 1 — repeatability at FIXED confirmed root
    # ================================================================
    hdr(f"PART 1  repeatability at fixed Ω={FIXED_OM}  (n_dofs={N_DOFS}, xmax={XMAX})")
    cnts = []
    n_raws = []
    for i in range(N_REPEATS):
        t1 = time.time()
        cnt, sel, n_raw = eval_point(R0_2B, TWO_T_PI, N_DOFS, XMAX, FIXED_OM)
        dt = time.time() - t1
        cnts.append(cnt)
        n_raws.append(n_raw)
        print(f"  rep {i+1}/{N_REPEATS}:  cnt={cnt:3d}  n_raw={n_raw:3d}  "
              f"dt={dt:.1f}s", flush=True)

    unique_cnt = sorted(set(cnts))
    print(f"\n  cnt sequence : {cnts}", flush=True)
    print(f"  unique cnt   : {unique_cnt}", flush=True)
    print(f"  n_raw sequence: {n_raws}", flush=True)

    if len(unique_cnt) == 1:
        print("  REPEATABILITY: IDENTICAL across all repeats  →  no per-call randomness",
              flush=True)
        fixed_cnt = unique_cnt[0]
        part1_ok = True
    else:
        print("  REPEATABILITY: cnt VARIES at fixed Ω  →  true nondeterminism",
              flush=True)
        fixed_cnt = None
        part1_ok = False

    # ================================================================
    # PART 2 — dense Ω scan across the narrow window
    # ================================================================
    hdr(f"PART 2  dense scan Ω ∈ [{OM_LO}, {OM_HI}]  step={OM_STEP}")
    oms = []
    scan_cnts = []
    scan_nraw = []
    om = OM_LO
    while om <= OM_HI + 1e-12:
        t1 = time.time()
        cnt, sel, n_raw = eval_point(R0_2B, TWO_T_PI, N_DOFS, XMAX, om)
        dt = time.time() - t1
        oms.append(om)
        scan_cnts.append(cnt)
        scan_nraw.append(n_raw)
        print(f"  Ω={om:.6f}  cnt={cnt:3d}  n_raw={n_raw:3d}  dt={dt:.1f}s",
              flush=True)
        om = round(om + OM_STEP, 10)   # avoid float drift

    unique_scan = sorted(set(scan_cnts))
    print(f"\n  scan cnt sequence : {scan_cnts}", flush=True)
    print(f"  unique scan cnt   : {unique_scan}", flush=True)
    print(f"  n_raw range       : {min(scan_nraw)} … {max(scan_nraw)}", flush=True)

    if len(unique_scan) >= 2:
        print("  SCAN: cnt takes >=2 distinct values  →  real Ω-threshold boundary present",
              flush=True)
        part2_ok = True
    else:
        print("  SCAN: cnt is CONSTANT across the whole window  →  no simple threshold effect",
              flush=True)
        part2_ok = False

    # ================================================================
    # PRE-REGISTERED READING
    # ================================================================
    hdr("READING (pre-registered)")
    print(f"  Part-1 (fixed-Ω identity) : {'PASS' if part1_ok else 'FAIL'}", flush=True)
    print(f"  Part-2 (scan variation)   : {'PASS' if part2_ok else 'FAIL'}", flush=True)

    if part1_ok and part2_ok:
        print("\n  VERDICT: PASS — confirms threshold-sensitivity hypothesis.", flush=True)
        print("  The 30-vs-34 discrepancy is explained by the converged Ω", flush=True)
        print("  landing on different sides of a select_fill branch-count", flush=True)
        print("  threshold inside a narrow Ω window; the machinery itself", flush=True)
        print("  is deterministic at fixed Ω.", flush=True)
    elif not part1_ok:
        print("\n  VERDICT: FAIL / kill — true nondeterminism at fixed Ω.", flush=True)
        print("  Follow-up must investigate randomness inside full_search", flush=True)
        print("  or select_fill (or an environmental source of non-reproducibility).", flush=True)
    else:
        print("\n  VERDICT: FAIL / kill — cnt constant across the scanned window.", flush=True)
        print("  Simple Ω-threshold sensitivity is ruled out. The original", flush=True)
        print("  30-vs-34 gap needs a different explanation (e.g. package", flush=True)
        print("  version / environment difference between the two jobs).", flush=True)

    print(f"\n  wall time: {time.time()-t0:.1f}s", flush=True)
    print("end", flush=True)


if __name__ == "__main__":
    main()