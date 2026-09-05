# -*- coding: utf-8 -*-
"""
probe_ff_oop_artifact_2p8176_v1.py -- DIAGNOSTIC/DATA-GENERATION ONLY.
Writes only a plain-text log; no package changes, no SOLVER_VERSION
implications.

EXACT COPY of probe_ff_oop_spurious_table.py (job 2325529, 13/13 clean)
with THREE changes and no others: the target list, the SOLVER_VERSION gate
(s8 -> s10), and this docstring. The worker, the residual screen call, the
verdict comparison and the summary logic are byte-identical to the run
that produced the paper's Table 5, so a verdict here is directly
comparable to that table.

WHY: the paper's Fig. 1 and Fig. 2(d) show a sixth FF-P1 OOP spurious zero
at raw Omega=2.817613 (EVEN), which sits ABOVE the validated flexural
window and is therefore not in Table 5. An earlier draft's Fig. 2 caption
told the reader to "see Sec. 5 for the residual values" of that zero --
but no r_V/r_M was ever computed for it, at any basis size, anywhere in
the project. The 2026-08-11 review pass caught the dangling pointer and
the caption was rewritten to stop promising numbers that do not exist.
This probe produces them, so the caption can instead state them.

WHAT THIS DOES: screens 2.817613 with the identical instrument used for
Table 5, alongside TWO controls carried in the same run -- one known
REAL-like mode (0.383695) and one known ARTIFACT-like zero (2.121685),
both with published rV/rM in Table 5. Controls are mandatory here (see the
plate-solver-cluster-probes skill): without them a surprising reading at
2.817613 cannot be separated from a change in the screen itself since
job 2325529 ran under s8.

PRE-REGISTERED INTERPRETATION (written BEFORE the job runs, honor it after):
  - CONTROLS FIRST. If 0.383695 does not return REAL-like with
    rM ~= 0.0025, or 2.121685 does not return ARTIFACT-like with
    rV ~= 620 / rM ~= 348 (the job 2325529 values, reproduced in the
    paper's Table 5), then the screen itself has moved between s8 and s10
    and NOTHING about 2.817613 may be concluded from this run. Report the
    control drift and stop.
  - If the controls reproduce and 2.817613 returns ARTIFACT-like with rM
    in the same ~75-350 band as Table 5's five in-window zeros -> the
    figure caption may state its rV/rM directly and describe it as a sixth
    spurious zero of the same population. This is the expected outcome.
  - If the controls reproduce but 2.817613 returns REAL-like or AMBIGUOUS
    -> do NOT quietly relabel it. That zero is currently described in the
    paper (Sec. 5, Fig. 2(d)) as spurious on the strength of its
    displacement-free reconstruction alone; a REAL-like or AMBIGUOUS
    reading would mean that description is unsupported and Sec. 5 plus
    both figure captions need revisiting. Report it plainly as a
    disconfirming result.
  - If 2.817613 fails to reconstruct (select_fill short of n_dofs=20, or
    no null vector) while both controls reconstruct fine -> report the
    failure; do not substitute a neighbouring Omega and do not fall back
    to a coarser basis.

Geometry/material: FF-P1, identical to the source probe --
  mat  = IsotropicMaterial(E=210e9, nu=0.30, rho=7800.0)
  geom = make_geometry(1.5, 0.5)      # r0/(2b)=1.5, 2*Theta/pi=0.5
  solver = OutOfPlaneSolver(geom, mat, M=80, n_quad=30, boundary=FreeFreeOOP())
  N_DOFS = 20, XI_MAX = 20.0

COST: 3 points x ~150 s serial ~= 7.5 CPU-minutes; 3 workers, a couple of
minutes wall plus queue. Budgeted at 30 min below.
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

# (raw Omega, Omega_lit label for readability, expected verdict).
# PHYSICAL/ARTIFACT here are CONTROLS with published Table-5 values; the
# unknown under test is the 2.817613 row.
PHYSICAL = [
    (0.383695, 15.148, "REAL-like"),        # control: Table 5, rV 0.0131 rM 0.00252
]
ARTIFACT = [
    (2.121685, None, "ARTIFACT-like"),      # control: Table 5, rV 620 rM 348
    (2.817613, None, "ARTIFACT-like"),      # <== THE UNKNOWN (EVEN, above window)
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
    print("  probe_ff_oop_artifact_2p8176_v1 -- start "
          + time.strftime("%Y-%m-%d %H:%M:%S"))
    print("=" * 78, flush=True)

    import plate_solver as ps
    if ps.SOLVER_VERSION != EXPECT_VER:
        print(f"PREFLIGHT FAIL: SOLVER_VERSION={ps.SOLVER_VERSION!r} != "
              f"{EXPECT_VER!r}")
        raise SystemExit(2)
    print(f"preflight OK: SOLVER_VERSION={ps.SOLVER_VERSION}", flush=True)

    n_workers = int(os.environ.get("N_WORKERS", "3"))
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
              "state 2.817613's rV/rM in the Fig. 2 caption.")
    print("=" * 78, flush=True)

    print("\nDiagnostic/data-generation only -- no SOLVER_VERSION action, "
          "no package changes.", flush=True)
    print("End: " + time.strftime("%Y-%m-%d %H:%M:%S"), flush=True)


if __name__ == "__main__":
    t0 = time.time()
    main()
    print(f"\ntotal wall time: {(time.time() - t0)/60:.1f} min", flush=True)
