# -*- coding: utf-8 -*-
"""
probe_ip_ffp1_ip06_n36_upper_window_v1.py -- DIAGNOSTIC ONLY. No package
changes, no SOLVER_VERSION impact. Plain-text log to stdout.

PURPOSE: finish the job 2411851 investigation, which ruled out the wrong
window.

WHAT 2411851 ESTABLISHED (keep, do not redo): at Omega = 1.19704984 with
n_dofs = 36, xmax = 32.0, five repeats returned cnt = 30 identically
(n_raw = 9 every time), and a scan over Omega in [1.196, 1.199] step 0.0002
returned cnt = 30 at all 16 points. So there is no per-call randomness and
no threshold effect in THAT window. Its own verdict was FAIL/kill on the
threshold hypothesis.

THE GAP: IP-06's tabulated n = 36 frequency in the paper's Appendix C.1 is
**1.203379**, which is OUTSIDE [1.196, 1.199]. The region where the n = 36
root actually sits was never tested. The paper currently says IP-06's
used@36 "was not reproducible run to run (30 in job 2324822 vs 34 in job
2324817); read as unresolved", and that wording cannot be revised on
2411851's evidence alone, because the disagreement may simply be an
Omega-dependence that 2411851's window did not reach.

This probe covers the missing region, at the same n_dofs/xmax/geometry, so
its numbers are directly comparable to 2411851's.

  Part 1: cnt at each of the THREE tabulated IP-06 frequencies --
          1.197860 (n=20 label), 1.197592 (n=28 label), 1.203379 (n=36
          label). These are the points the paper's own table is built on.
  Part 2: scan Omega in [1.1990, 1.2050] step 0.0002 (31 points), which
          butts directly onto 2411851's upper edge, so the two runs
          together cover [1.196, 1.205] with no gap.
  Part 3: five repeats at Omega = 1.203379 specifically -- the n = 36
          label's own frequency, where 2411851 never tested repeatability.

PRE-REGISTERED INTERPRETATION (written before the job runs, honour it):

  A) If cnt == 30 at every point in Parts 1-3 and the Part-3 repeats are
     identical:
     -> The count is DETERMINISTIC and CONSTANT across the whole relevant
        region [1.196, 1.205]. The historical 34 (job 2324817) is not
        reproducible under s10 and is an artifact of that earlier run's
        environment/version, not of the routine.
     -> Appendix C.1 may then be revised: replace "not reproducible run to
        run ... read as unresolved" with a statement that the two 2026-07
        jobs disagreed and a dedicated re-test under the current version
        returns 30 deterministically. The conv-ip table's "30--34" cell
        becomes 30. THIS IS THE ONLY OUTCOME THAT LICENSES CHANGING THE
        PAPER.

  B) If cnt is 30 in some sub-range and 34 (or any other value) in another,
     with a clean boundary:
     -> There IS a threshold effect, sitting above 2411851's window. The
        30-vs-34 split is then real Omega-dependence, not irreproducibility.
     -> Report the boundary Omega. Appendix C.1's numbers stay, but the
        explanation changes from "unresolved" to "Omega-dependent basis
        fill", which is a better and more honest account. Do NOT collapse
        the table cell to a single number in this case.

  C) If the Part-3 repeats at 1.203379 DISAGREE with each other:
     -> Genuine nondeterminism, localised to a frequency 2411851 did not
        test. This is the alarming case and the only one that is a real
        defect rather than a bookkeeping question. Report it loudly; do not
        touch the paper; the next step would be to find what in the branch
        census is order- or state-dependent there.

  D) If cnt == 34 at 1.203379 but is stable:
     -> The original 34 is reproducible and simply belongs to that
        frequency. The paper's "30--34" range is correct as printed and the
        right fix is to say which Omega each belongs to, not to call it
        unresolved.

  Any of A-D is a usable result. There is no outcome here that wastes the
  run, which is the point of covering the window rather than arguing about
  it.

FIDELITY GATE: Part 0 re-evaluates 2411851's own fixed point
(Omega = 1.19704984) and must return cnt = 30, n_raw = 9. If it does not,
this run is not comparable to 2411851 and nothing below is interpretable.

COST: 3 + 31 + 5 + 1 = 40 evaluations at ~34 s each (2411851's measured
rate) = ~23 CPU-minutes serial. Run serially in Omega order so the scan
reads as a curve; ~25 min wall. Budgeted 1 h.
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

# Identical to job 2411851 -- do not change, or the runs stop being comparable.
R0_2B = 1.5
TWO_T_PI = 0.5
N_DOFS = 36
XMAX = 32.0

GATE_OM = 1.19704984      # 2411851's fixed point
GATE_CNT = 30             # what it returned, 5/5
GATE_NRAW = 9

# The paper's own three IP-06 labels (Appendix C.1, tab:conv-ip)
TABULATED = [
    ("IP-06 n=20 label", 1.197860),
    ("IP-06 n=28 label", 1.197592),
    ("IP-06 n=36 label", 1.203379),   # <-- the one 2411851 never reached
]

SCAN_LO, SCAN_HI, SCAN_STEP = 1.1990, 1.2050, 0.0002
REPEAT_OM = 1.203379
N_REPEATS = 5


def hdr(s):
    print("\n" + "=" * 78 + f"\n  {s}\n" + "=" * 78, flush=True)


def eval_point(Om):
    """Verbatim from job 2411851's eval_point -- same solver, same call."""
    mat = ps.IsotropicMaterial(E=210e9, nu=0.30, rho=7800.0)
    geom = ps.make_geometry(R0_2B, TWO_T_PI)
    solver = ps.InPlaneSolver(geom, mat, M=80, n_quad=30,
                              boundary=ps.FreeFreeIP())
    raw = full_search(solver.fast, Om, xmax=XMAX)
    sel, cnt = select_fill(raw, N_DOFS)
    return cnt, len(raw)


def main():
    t0 = time.time()
    print(f"job start  host={os.uname().nodename}  dps={mp.dps}", flush=True)
    print(f"EXPECT_SOLVER_VERSION={EXPECT_VER}", flush=True)
    ver = getattr(ps, "SOLVER_VERSION", None) or "?"
    print(f"package SOLVER_VERSION = {ver!r}", flush=True)
    if ver != EXPECT_VER:
        print(f"FATAL: version mismatch (got {ver!r})", flush=True)
        sys.exit(3)
    print("version gate PASS", flush=True)
    if full_search is None:
        print("FATAL: full_search not importable", flush=True)
        sys.exit(3)
    print(f"  geometry r0/2b={R0_2B}  2T/pi={TWO_T_PI}  "
          f"n_dofs={N_DOFS}  xmax={XMAX}", flush=True)

    # ---------------- Part 0: fidelity gate vs job 2411851 -----------
    hdr(f"PART 0  fidelity gate -- 2411851's fixed point Om={GATE_OM}")
    c, nr = eval_point(GATE_OM)
    print(f"  cnt={c}  n_raw={nr}   (job 2411851: cnt={GATE_CNT}, "
          f"n_raw={GATE_NRAW})", flush=True)
    if c != GATE_CNT or nr != GATE_NRAW:
        print("  GATE FAIL -- this run is NOT comparable to 2411851.")
        print("  Nothing below is interpretable. Stopping.", flush=True)
        sys.exit(4)
    print("  GATE PASS", flush=True)

    # ---------------- Part 1: the three tabulated labels -------------
    hdr("PART 1  cnt at the three tabulated IP-06 frequencies")
    part1 = []
    for name, Om in TABULATED:
        t1 = time.time()
        c, nr = eval_point(Om)
        part1.append((name, Om, c, nr))
        print(f"  {name:<18s} Om={Om:.6f}  cnt={c:3d}  n_raw={nr:3d}  "
              f"dt={time.time()-t1:.1f}s", flush=True)

    # ---------------- Part 2: the missing window ---------------------
    n_pts = int(round((SCAN_HI - SCAN_LO) / SCAN_STEP)) + 1
    hdr(f"PART 2  scan Om in [{SCAN_LO}, {SCAN_HI}] step {SCAN_STEP} "
        f"({n_pts} pts) -- the region 2411851 did not cover")
    scan = []
    for i in range(n_pts):
        Om = SCAN_LO + i * SCAN_STEP
        t1 = time.time()
        c, nr = eval_point(Om)
        scan.append((Om, c, nr))
        mark = "   <-- IP-06 n=36 label" if abs(Om - 1.203379) < 1e-4 else ""
        print(f"  Om={Om:.6f}  cnt={c:3d}  n_raw={nr:3d}  "
              f"dt={time.time()-t1:.1f}s{mark}", flush=True)

    # ---------------- Part 3: repeatability at the n=36 label --------
    hdr(f"PART 3  repeatability at Om={REPEAT_OM} (the n=36 label itself)")
    reps = []
    for i in range(N_REPEATS):
        t1 = time.time()
        c, nr = eval_point(REPEAT_OM)
        reps.append((c, nr))
        print(f"  rep {i+1}/{N_REPEATS}:  cnt={c:3d}  n_raw={nr:3d}  "
              f"dt={time.time()-t1:.1f}s", flush=True)

    # ---------------- reading ---------------------------------------
    hdr("READING (pre-registered)")
    scan_cnts = sorted({c for _, c, _ in scan})
    all_cnts = sorted({c for _, _, c, _ in part1} | set(scan_cnts)
                      | {c for c, _ in reps})
    rep_cnts = sorted({c for c, _ in reps})
    n36 = [c for name, Om, c, nr in part1 if "n=36" in name][0]

    print(f"  Part 1 cnts        : {[c for _,_,c,_ in part1]}")
    print(f"  Part 2 unique cnts : {scan_cnts}")
    print(f"  Part 3 rep cnts    : {[c for c,_ in reps]}")
    print(f"  cnt at the n=36 label (1.203379) = {n36}")
    print(f"  union of all cnts observed = {all_cnts}")

    if len(rep_cnts) > 1:
        verdict = ("C -- GENUINE NONDETERMINISM at Om=1.203379. Repeats "
                   "disagree. Do NOT touch the paper; investigate what in "
                   "the branch census is order/state dependent here.")
    elif all_cnts == [30]:
        verdict = ("A -- DETERMINISTIC AND CONSTANT at 30 across "
                   "[1.196, 1.205]. The historical 34 (job 2324817) does "
                   "not reproduce under s10. Appendix C.1 may be revised "
                   "and the conv-ip cell set to 30.")
    elif n36 == 34 and len(scan_cnts) == 1:
        verdict = ("D -- cnt is a stable 34 in this upper window. The "
                   "original 34 belongs to this frequency; the printed "
                   "30--34 range is correct and should be explained as "
                   "Omega-dependence, not irreproducibility.")
    else:
        verdict = ("B -- THRESHOLD EFFECT above 2411851's window. cnt "
                   "varies across the scan; report the boundary. Keep the "
                   "printed range, change the explanation from "
                   "'unresolved' to 'Omega-dependent basis fill'.")

    print(f"\n  VERDICT: {verdict}")

    if len(scan_cnts) > 1:
        print("\n  transitions in the scan:")
        for i in range(1, len(scan)):
            if scan[i][1] != scan[i - 1][1]:
                print(f"    cnt {scan[i-1][1]} -> {scan[i][1]} between "
                      f"Om={scan[i-1][0]:.6f} and {scan[i][0]:.6f}")

    print(f"\n  wall time: {time.time()-t0:.1f}s")
    print("end", flush=True)


if __name__ == "__main__":
    main()
