# -*- coding: utf-8 -*-
"""
probe_ff_oop_spurious_table.py -- DIAGNOSTIC/DATA-GENERATION ONLY. Writes
only a plain-text log; no package changes, no SOLVER_VERSION implications.

CONTEXT: PAPER1_FREEFREE_DRAFT.tex's tab:oop caption (2026-07-19 edit)
notes that its 8 physical FF-P1 OOP modes are shown but the 5 catalogued
spurious determinant zeros for this geometry are not tabulated. A 2026-07-19
external review (SuperGrok, `Project Knowledge/ReviewGrok2.txt`) asked for
those 5 zeros' raw Omega and residual-screen numbers in a supplementary
table. The only numbers on record for them
(`LESSONS_LEARNED.md` Sec 10.8, captured 2026-07-09; corroborated by
`FREE_FREE_MIDPOINT_REPORT.md` Sec 3, dated 2026-07-13, which lists the
same 5 raw Omega as "*" artifact rows retained in the reproducibility pin
alongside the CURRENT 8-physical-mode set) predate several solver-adjacent
changes made since, so this project's own methodology (validate by direct
probe, never reuse stale numbers into a paper) requires reconfirming them
fresh against the currently deployed checkpoint before they go in print --
this probe does exactly that, and nothing else.

WHAT THIS DOES: for the 13 already-known raw Omega values that make up the
FF-P1 OOP free-free spectrum at n_dofs=20 (the same "FF-P1 default" basis
size used throughout the free-free closure campaign and by
probe_paper_figures_oop.py) -- the 8 physical modes currently in
PAPER1_FREEFREE_DRAFT.tex's tab:oop, plus the 5 previously-catalogued
spurious zeros -- this probe:
  1. Reconstructs the branch set at that Omega via full_search + select_fill
     (NOT a new discovery search -- these are fixed reconstructions at
     already-validated roots, the same technique probe_paper_figures_oop.py
     already used for the mode-shape PNGs; no seeding of the residual
     screen or verdict logic itself).
  2. Runs detectors.weak_enforcement_residual_oop -- the pointwise
     edge-residual screen already validated 11/11 on this exact geometry
     (LESSONS Sec 10.7/10.8) -- to get block/rV/rM/maxW/verdict.
  3. Cross-checks solver.sigma_min at the same (Omega, sel) as an
     independent depth sanity check (cheap, reuses sel, no extra
     full_search).
  4. Compares the returned verdict against the PRE-REGISTERED expected
     label below and flags any MISMATCH loudly rather than silently
     accepting it.

PRE-REGISTERED INTERPRETATION (write this BEFORE the job runs, honor it
after):
  - If all 8 physical points return verdict=="REAL-like" and all 5
    artifact points return verdict=="ARTIFACT-like", with rM for the
    artifacts landing in the same ~75-350 range LESSONS Sec 10.8 recorded
    (vs. <0.04 for the physical points) -> the stale numbers were correct,
    were never at risk from a solver change, and this run's fresh numbers
    are what go directly into the paper's supplementary table.
  - If ANY point's verdict flips relative to the record above (e.g. a
    physical point comes back ARTIFACT-like/AMBIGUOUS, or vice versa), or
    an artifact's rM has moved into the REAL-like band -> DO NOT put the
    new numbers in the paper silently. That would mean either a solver
    behavior change since 2026-07-09/13 (unexpected, since SOLVER_VERSION
    is unchanged at 2026-07-10.s8 throughout) or a transcription error in
    this probe's OMEGA_ARTIFACT/OMEGA_PHYSICAL lists -- stop and
    re-examine both before touching the paper.
  - If any point fails to reconstruct (select_fill short of n_dofs=20, or
    "no null vector") -> report the failure plainly; do not substitute a
    neighboring Omega.

Geometry/material (FF-P1, same pinned instance used throughout, matching
probe_paper_figures_oop.py exactly):
  mat  = IsotropicMaterial(E=210e9, nu=0.30, rho=7800.0)
  geom = make_geometry(1.5, 0.5)      # r0/(2b)=1.5, 2*Theta/pi=0.5
  solver = OutOfPlaneSolver(geom, mat, M=80, n_quad=30, boundary=FreeFreeOOP())
  N_DOFS = 20, XI_MAX = 20.0 (FF-P1 default basis/search-radius elsewhere)

COST: 13 points, each ~1 full_search (~55-60s at dps=40 per the cluster-
probes skill's cost model) + 1 residual-screen build (comparable cost to a
full K build with per-item series/amplitude reconstruction, ~60-90s) +
1 cheap sigma_min reuse (~1-2s) ~= 13 x ~150s serial = ~32.5 CPU-minutes
total. Parallelized 1 point/worker across up to 13 workers, wall time
should be a few minutes plus queue overhead; budgeted generously below.
"""
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

sys.path.insert(0, os.environ.get("PKG_PATH", "."))
os.environ.setdefault("DPS", "40")

EXPECT_VER = os.environ.get("EXPECT_SOLVER_VERSION", "2026-07-10.s10")

N_DOFS = 20
XI_MAX = 20.0

# (raw Omega, Omega_lit label for readability, expected verdict) -- the 8
# physical modes match PAPER1_FREEFREE_DRAFT.tex's CURRENT tab:oop exactly;
# the 5 artifact rows match both LESSONS_LEARNED.md Sec 10.8 (2026-07-09)
# and FREE_FREE_MIDPOINT_REPORT.md Sec 3 (2026-07-13, "*" rows).
PHYSICAL = [
    (0.383695, 15.148, "REAL-like"),
    (0.616645, 24.344, "REAL-like"),
    (0.961312, 37.951, "REAL-like"),
    (1.369611, 54.070, "REAL-like"),
    (1.752188, 69.174, "REAL-like"),
    (2.320741, 91.619, "REAL-like"),
    (2.379886, 93.954, "REAL-like"),
    (2.602500, 102.743, "REAL-like"),
]
ARTIFACT = [
    (0.320670, None, "ARTIFACT-like"),
    (0.844488, None, "ARTIFACT-like"),
    (1.001419, None, "ARTIFACT-like"),
    (1.602090, None, "ARTIFACT-like"),
    (2.121685, None, "ARTIFACT-like"),
]
ALL_TARGETS = [("physical", *t) for t in PHYSICAL] + \
              [("artifact", *t) for t in ARTIFACT]


def _worker(args):
    kind, Om, Om_lit, expected = args
    t0 = time.time()
    try:
        import plate_solver as ps
        from plate_solver.detectors import (full_search, select_fill,
                                             weak_enforcement_residual_oop)

        mat = ps.IsotropicMaterial(E=210e9, nu=0.30, rho=7800.0)
        geom = ps.make_geometry(1.5, 0.5)
        solver = ps.OutOfPlaneSolver(geom, mat, M=80, n_quad=30,
                                      boundary=ps.FreeFreeOOP())
        raw = full_search(solver.fast, Om, xmax=XI_MAX)
        sel, cnt = select_fill(raw, N_DOFS)
        if cnt < N_DOFS:
            return dict(ok=False, kind=kind, Om=Om, Om_lit=Om_lit,
                        expected=expected,
                        error=f"only {cnt}/{N_DOFS} branches filled",
                        dt=time.time() - t0)

        res = weak_enforcement_residual_oop(solver, Om, sel, n_dofs=N_DOFS)
        if not res.get("ok", False):
            return dict(ok=False, kind=kind, Om=Om, Om_lit=Om_lit,
                        expected=expected,
                        error=f"residual screen failed: {res.get('reason')}",
                        dt=time.time() - t0)

        s = solver.sigma_min(Om, sel)

        verdict = res["verdict"]
        mismatch = (verdict != expected)
        return dict(ok=True, kind=kind, Om=Om, Om_lit=Om_lit,
                    expected=expected, verdict=verdict, mismatch=mismatch,
                    block=res["block"], rV=res["rV"], rM=res["rM"],
                    maxW=res["maxW"], log10_sigma_min=float(s),
                    dt=time.time() - t0)
    except Exception as exc:
        return dict(ok=False, kind=kind, Om=Om, Om_lit=Om_lit,
                    expected=expected, error=f"{exc}", dt=time.time() - t0)


def main():
    print("=" * 78)
    print("  probe_ff_oop_spurious_table -- start "
          + time.strftime("%Y-%m-%d %H:%M:%S"))
    print("=" * 78, flush=True)

    import plate_solver as ps
    if ps.SOLVER_VERSION != EXPECT_VER:
        print(f"PREFLIGHT FAIL: SOLVER_VERSION={ps.SOLVER_VERSION!r} != "
              f"{EXPECT_VER!r}")
        raise SystemExit(2)
    print(f"preflight OK: SOLVER_VERSION={ps.SOLVER_VERSION}", flush=True)

    n_workers = int(os.environ.get("N_WORKERS", "13"))
    print(f"N_WORKERS={n_workers}", flush=True)
    print(f"{len(ALL_TARGETS)} targets: {len(PHYSICAL)} physical + "
          f"{len(ARTIFACT)} artifact, N_DOFS={N_DOFS}, XI_MAX={XI_MAX}\n",
          flush=True)

    results = []
    with ProcessPoolExecutor(max_workers=min(n_workers, len(ALL_TARGETS))) as ex:
        futs = {ex.submit(_worker, t): t for t in ALL_TARGETS}
        done = 0
        for fut in as_completed(futs):
            r = fut.result()
            done += 1
            if r["ok"]:
                flag = "  <== MISMATCH" if r["mismatch"] else ""
                lit = f"{r['Om_lit']:.3f}" if r["Om_lit"] is not None else "n/a"
                print(f"  [{done}/{len(ALL_TARGETS)}] {r['kind']:8s} "
                      f"Om={r['Om']:.6f} (lit={lit}) block={r['block']:>3s} "
                      f"rV={r['rV']:.4g} rM={r['rM']:.4g} maxW={r['maxW']:.4g} "
                      f"log10smin={r['log10_sigma_min']:.4f} "
                      f"verdict={r['verdict']} (expected {r['expected']})"
                      f"{flag} [{r['dt']:.0f}s]", flush=True)
            else:
                print(f"  [{done}/{len(ALL_TARGETS)}] {r['kind']:8s} "
                      f"Om={r['Om']:.6f} FAILED: {r.get('error')} "
                      f"[{r['dt']:.0f}s]", flush=True)
            results.append(r)

    ok_results = [r for r in results if r["ok"]]
    mismatches = [r for r in ok_results if r["mismatch"]]
    failures = [r for r in results if not r["ok"]]

    print("\n" + "=" * 78)
    print(f"SUMMARY: {len(ok_results)}/{len(ALL_TARGETS)} ok, "
          f"{len(mismatches)} verdict mismatches, {len(failures)} failures")
    if mismatches:
        print("MISMATCHES (do not use in paper without re-examination):")
        for r in mismatches:
            print(f"  Om={r['Om']:.6f} kind={r['kind']} "
                  f"expected={r['expected']} got={r['verdict']} rM={r['rM']:.4g}")
    if failures:
        print("FAILURES:")
        for r in failures:
            print(f"  Om={r['Om']:.6f} kind={r['kind']}: {r.get('error')}")
    if not mismatches and not failures:
        print("All verdicts matched pre-registered expectations. Safe to "
              "transcribe rV/rM/maxW directly into the paper's "
              "supplementary spurious-zero table.")
    print("=" * 78, flush=True)

    print("\nDiagnostic/data-generation only -- no SOLVER_VERSION action, "
          "no package changes.", flush=True)
    print("End: " + time.strftime("%Y-%m-%d %H:%M:%S"), flush=True)


if __name__ == "__main__":
    t0 = time.time()
    main()
    print(f"\ntotal wall time: {(time.time() - t0)/60:.1f} min", flush=True)
