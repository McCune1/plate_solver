# -*- coding: utf-8 -*-
"""
probe_ff_residual_screen_reconstruction_fix.py -- DIAGNOSTIC/DATA-
GENERATION ONLY. Writes only a plain-text log; no package changes, no
SOLVER_VERSION implications. Lives alongside, and imports TARGETS
verbatim from, probe_ff_residual_screen_geomsweep.py (jobs
2327073-2327076) -- run this from the SAME directory.

CONTEXT: jobs 2327073-2327076 ran the Sec.5 pointwise edge-residual
screen (weak_enforcement_residual_oop/_ip) against 215 cut-off-adjacent
candidates across the 48 new free-free geometries. 44 of those 215 (12
OOP, 32 IP -- confirmed by grepping every "FAILED" line from all 4 jobs'
.out logs, matching PAPER1_FREEFREE_DRAFT.tex's own "12/119 flexural,
32/96 extensional" count exactly) could not be reconstructed at all: the
branch-fill step came back short of n_dofs before the residual screen
ever ran. This is written into the paper (Sec 6.4) and FUTURE_WORK.md
Sec 1 as an OPEN, NOT-ROOT-CAUSED item: "may reflect the same
branch-selection sensitivity already documented for other parts of this
solver (Sec 4), but this is a guess, not a finding."

THE CONCRETE, GROUNDED HYPOTHESIS THIS PROBE TESTS: it isn't a guess.
Reading probe_ff_residual_screen_geomsweep.py's `_worker` directly shows
every failure came from exactly ONE cold call --
`select_fill(full_search(solver.fast, Om, xmax=xmax), n_dofs)` -- with NO
fallback. This project has already hit, diagnosed, and FIXED this exact
failure mode once before, in a different probe
(probe_threshold_sensitivity.py, 2026-07-17/18/19, jobs 2324728/2324735/
2324736/2324750/2325211/2325366): a single cold full_search at a fixed,
already-known Omega routinely under-fills, and the fix
(`_derive_branches_robust`, copied verbatim below) escalates through
three further stages -- nearby-offset-then-track, widened xmax with
denser grid, and a guarded xmax_cap=90 step verified by a real sigma_min
check -- before giving up. That escalation ladder is validated: job
2325369 ran it clean across all 352 threshold-sensitivity points (both
OOP and IP), and its guard was independently checked correct on 25
points (job 2325366). It was built specifically to fix IP branch-ceiling
underfill (job 2324740 found this hit up to 30/30 IP points at
n_dofs>=24), which is exactly the failure signature here (32/44 failures
are IP). probe_ff_residual_screen_geomsweep.py's `_worker` was simply
never updated to use it -- so before concluding anything physical about
these 44 candidates, the obvious, cheap, already-proven fix should be
tried.

WHAT THIS PROBE DOES: for exactly the 44 originally-failing candidates
(identified below by FAILED_KEYS, matched back against the verbatim
TARGETS list imported from probe_ff_residual_screen_geomsweep.py -- NOT
retyped, to avoid a transcription error in 44 sets of floating-point
Omega/n_dofs/xmax values), reruns branch reconstruction using
`_derive_branches_robust` instead of the original bare cold call. Any
candidate that now fills to n_dofs runs the SAME production residual
screen (weak_enforcement_residual_oop/_ip, unmodified) the original probe
would have run. Any candidate that STILL doesn't fill -- even after
escalating to xmax_cap=90 with a 2x-denser Newton-seed grid -- is
reported as a materially stronger negative result than the original
"cold search failed," since it now means the already-proven-robust
escalation ladder itself was exhausted, not just an under-engineered
first attempt.

IN-RUN FIDELITY CHECK (per this project's monkeypatch-discipline norm --
see plate-solver-cluster-probes skill Sec 3 -- applied here even though
`_derive_branches_robust` isn't a modified package method, because it IS
copied from another file rather than imported): before touching any of
the 44 failures, this probe re-runs the ORIGINAL cold-search step on one
already-SUCCEEDING candidate from the source run (the first 'ok' row
implied by TARGETS minus FAILED_KEYS) and confirms `_derive_branches_
robust` returns `how="cold"` with a `cnt` that satisfies n_dofs on the
very first stage -- i.e. that the copy doesn't accidentally change
behavior on cases that already worked. If that gate fails, the copy has
a transcription bug and this probe stops before processing anything else.

PRE-REGISTERED INTERPRETATION (write before running, honor after):
  - If MOST/ALL of the 44 now fill successfully via `_derive_branches_
    robust` -> the reconstruction-failure caveat is CLOSED as an
    implementation gap, not a solver-physics finding: fold the completed
    verdicts (REAL-like/ARTIFACT-like/AMBIGUOUS for OOP; raw ratios for
    IP) into the existing Sec 6.4 open-item note, updating "12/119 (OOP)
    and 32/96 (IP) ... root cause not investigated" to the resolved
    figures. This directly shrinks (or removes) an open item currently
    stated in the published paper text.
  - If OOP resolves well but IP does not (or resolves less completely)
    -> report the split explicitly. This would be informative on its own:
    it would mean the escalation ladder's IP-side robustness, proven at
    FF-P1 (nu=0.30, r0/2b=1.5), does not fully transfer to nu=0.35 and
    the 4 new radius ratios -- consistent with this project's already-
    stated "first transfer" caution for the residual screen itself, and
    would sharpen rather than close the open item for IP specifically.
  - If NONE resolve even at the guarded xmax_cap=90 stage -> this
    strengthens, not weakens, the open item: it shows the failure is not
    an implementation shortcut (the strongest available fallback was
    tried and still failed), so a genuine structural branch-scarcity at
    these specific (Omega, geometry) combinations becomes the leading
    explanation. Report this plainly as a stronger negative finding, and
    note in any paper update that this specific fix was tried and ruled
    out (not simply "not investigated").
  - For every point, record the escalation stage that succeeded or the
    stage at which it gave up (`how` from `_derive_branches_robust`), and
    the ORIGINAL shortfall depth (n_dofs minus the original cold-search
    cnt, read from the matching FAILED-line data already captured in this
    docstring's key list). A correlation between shortfall depth and
    escalation-resistance would itself be diagnostic and should be
    reported, not just the pass/fail count.
  - Do not average, smooth, or drop any of the 44 -- report every one,
    success or failure, exactly like the original probe did.

COST: `_derive_branches_robust` can be considerably more expensive per
point than a bare cold search -- per probe_threshold_sensitivity.py's own
COST UPDATE note, a point that reaches the guarded xmax_cap stage (xmax=
90, ngrid=20, nax=480) costs roughly 830-1200s, not the ~35-60s the
ORIGINAL cold-only attempt took on these same 44 points (per their own
logged per-point times). Budgeting conservatively at up to ~1200s/point
worst case x 44 points = ~14.7 CPU-hours serial; parallelized 1
point/worker (N_WORKERS sized to the SLURM allocation) with
--cpus-per-task=44, expect a few minutes if most points resolve at an
early (cheap) stage, up to ~20-25 min if a large fraction need the
guarded final stage. --time=03:00:00 is a generous margin, matching the
original threshold-sensitivity job's own budget for a similarly-shaped
cost profile.

Env var RATIO_TAG is NOT used here (unlike the source probe) -- this
probe always processes all 44 failures from all 4 original jobs in one
run, since the point is specifically the previously-failed set, not a
per-ratio partition.
"""
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

sys.path.insert(0, os.environ.get("PKG_PATH", "."))
os.environ.setdefault("DPS", "40")

EXPECT_VER = os.environ.get("EXPECT_SOLVER_VERSION", "2026-07-10.s8")

# Same material the original geometry-sweep run actually used (see
# probe_ff_residual_screen_geomsweep.py's own NU_SWEEP note /
# LESSONS_LEARNED.md Sec 22 item 3) -- NOT the FF-P1 tab:oop/tab:ip
# benchmark's nu=0.30.
NU_SWEEP = 0.35

# ---------------------------------------------------------------------
# Verbatim import of the full 215-row TARGETS list from the source probe
# -- avoids retyping 44 rows of floating-point (Omega, n_dofs, xmax)
# data by hand. Run this script from the SAME directory as
# probe_ff_residual_screen_geomsweep.py (matches that probe's own
# no-explicit-cd submit-script convention).
#
# NOTE: probe_ff_residual_screen_geomsweep.py has module-level code that
# calls `raise SystemExit(2)` if the RATIO_TAG env var isn't set (its own
# CLI guard). Importing it here triggers that same module-level code, so
# the submit script for THIS probe must also export RATIO_TAG to any one
# valid value (e.g. r150) purely to satisfy that guard -- this probe's
# own logic never reads RATIO_TAG itself; it always processes all 44
# failures from all 4 original ratios in one run. Verified in the sandbox
# (RATIO_TAG=r150 import test) that this import then yields all 215
# TARGETS rows correctly, independent of which tag was set.
# ---------------------------------------------------------------------
try:
    from probe_ff_residual_screen_geomsweep import TARGETS as _ALL_TARGETS
except ImportError:
    print("FATAL: could not import TARGETS from probe_ff_residual_screen_"
          "geomsweep.py -- run this script from the same directory as that "
          "file (FutureWork/geometry_sweep/).")
    raise SystemExit(2)

# The 44 failing (r0_2b, two_T_pi, part, round(Omega,6)) keys, extracted
# directly from every "FAILED" line across ff_residual_screen_r150_
# 2327073.out, _r167_2327074.out, _r200_2327075.out, _r250_2327076.out
# (grepped, not hand-derived) -- 8+9+9+18=44, matching PAPER1_FREEFREE_
# DRAFT.tex's "12/119 flexural, 32/96 extensional" exactly (12 part=1
# rows, 32 part=2 rows below).
FAILED_KEYS = {
    # r150 (8, all IP)
    (1.5, 0.5, 2, 0.407825), (1.5, 0.25, 2, 1.031921),
    (1.5, 0.5, 2, 0.448825), (1.5, 0.25, 2, 0.320915),
    (1.5, 0.5, 2, 0.320915), (1.5, 0.25, 2, 0.021597),
    (1.5, 0.25, 2, 0.831487), (1.5, 0.5, 2, 0.547282),
    # r167 (9: 8 IP + 1 OOP)
    (1.66667, 0.5, 2, 0.392910), (1.66667, 0.5, 2, 0.344895),
    (1.66667, 0.25, 2, 0.260820), (1.66667, 0.25, 2, 0.925569),
    (1.66667, 0.5, 2, 0.017750), (1.66667, 0.25, 2, 0.017751),
    (1.66667, 0.25, 2, 0.738254), (1.66667, 0.5, 2, 0.260836),
    (1.66667, 0.25, 1, 1.144800),
    # r200 (9: 8 IP + 1 OOP)
    (2.0, 0.5, 2, 0.311546), (2.0, 0.5, 2, 0.395857),
    (2.0, 0.25, 2, 0.177806), (2.0, 0.25, 2, 0.012306),
    (2.0, 0.5, 2, 0.543552), (2.0, 0.25, 2, 0.742474),
    (2.0, 0.5, 2, 0.012412), (2.0, 0.25, 2, 0.597162),
    (2.0, 0.25, 1, 1.468229),
    # r250 (18: 8 IP + 10 OOP)
    (2.5, 0.25, 2, 0.456193), (2.5, 0.25, 2, 0.020274),
    (2.5, 0.25, 2, 0.110248), (2.5, 0.5, 2, 0.110223),
    (2.5, 0.5, 2, 0.173706), (2.5, 0.25, 2, 0.008006),
    (2.5, 0.5, 2, 0.008006), (2.5, 0.5, 2, 0.022797),
    (2.5, 0.5, 1, 0.294612), (2.5, 0.5, 1, 0.166257),
    (2.5, 0.5, 1, 0.126001), (2.5, 0.25, 1, 0.510814),
    (2.5, 0.25, 1, 0.535827), (2.5, 0.25, 1, 0.464744),
    (2.5, 0.5, 1, 0.042341), (2.5, 0.5, 1, 0.005831),
    (2.5, 0.25, 1, 0.005985), (2.5, 0.25, 1, 0.689079),
}

FAILED = [row for row in _ALL_TARGETS
          if (row[1], row[2], row[3], round(row[6], 6)) in FAILED_KEYS]

# One already-succeeding row (first TARGETS entry NOT in FAILED_KEYS) --
# used only for the in-run fidelity gate below, never scored as a "new"
# result.
_FIDELITY_ROW = next(row for row in _ALL_TARGETS
                      if (row[1], row[2], row[3], round(row[6], 6))
                      not in FAILED_KEYS)


# ---------------------------------------------------------------------
# _derive_branches_robust -- copied VERBATIM from probe_threshold_
# sensitivity.py (2026-07-19 "FOURTH PASS" version, the one job 2325369
# validated 352/352 with a 25/25-checked guard). See that file for the
# full derivation history. Do NOT hand-edit without re-diffing against
# the source.
# ---------------------------------------------------------------------
def _derive_branches_robust(eng, Omega, n_dofs, max_dim, full_search, track,
                             select_fill, xmax_cap=90.0,
                             solver=None, sigma_min_from_K=None,
                             shallow_floor=-3.0):
    raw = full_search(eng, Omega, xmax=max_dim)
    sel, cnt = select_fill(raw, n_dofs)
    if cnt >= n_dofs:
        return sel, cnt, "cold"

    for frac in (0.01, 0.02, 0.04, 0.08, 0.15):
        for sign in (+1, -1):
            near = Omega * (1.0 + sign * frac)
            if near <= 0:
                continue
            near_raw = full_search(eng, near, xmax=max_dim)
            tracked = track(eng, Omega, near_raw)
            sel, cnt = select_fill(tracked, n_dofs)
            if cnt >= n_dofs:
                return sel, cnt, f"tracked({sign*frac:+.2f})"

    scale = 1.0
    for _ in range(8):
        scale *= 1.25
        xmax_try = max_dim * scale
        if xmax_try > xmax_cap:
            break
        grid_scale = min(2.0, scale)
        ngrid_try = max(10, round(10 * grid_scale))
        nax_try = max(240, round(240 * grid_scale))
        raw = full_search(eng, Omega, xmax=xmax_try, ngrid=ngrid_try,
                           nax=nax_try)
        sel, cnt = select_fill(raw, n_dofs)
        if cnt >= n_dofs:
            return sel, cnt, (f"escalated(xmax={xmax_try:.1f},"
                               f"ngrid={ngrid_try},nax={nax_try})")

    raw = full_search(eng, Omega, xmax=xmax_cap, ngrid=20, nax=480)
    sel, cnt = select_fill(raw, n_dofs)
    if cnt >= n_dofs:
        if solver is not None and sigma_min_from_K is not None:
            K, size = solver._build_K_real(Omega, sel, lagrange=True,
                                            fast_scan=False)
            smin = sigma_min_from_K(K, size)
            if smin <= shallow_floor:
                return sel, cnt, (f"xmax_cap(verified,xmax={xmax_cap:.1f},"
                                   f"ngrid=20,nax=480,smin={smin:.3f})")
            return None, cnt, (f"exhausted(cap_reached_but_shallow,"
                                f"smin={smin:.3f})")
        return sel, cnt, f"xmax_cap(unverified,xmax={xmax_cap:.1f})"

    return None, cnt, "exhausted"


def _make_solver(part):
    import plate_solver as ps
    mat = ps.IsotropicMaterial(E=210e9, nu=NU_SWEEP, rho=7800.0)
    return ps, mat


def _fidelity_gate():
    """Confirms the verbatim copy above behaves identically to the
    original on an already-succeeding case before any of the 44 failures
    are touched. Returns (ok, detail)."""
    import plate_solver as ps
    from plate_solver.detectors import full_search, track, select_fill, \
        sigma_min_from_K

    tag, r0_2b, two_T_pi, part, n_dofs, xmax, Om, Om_lit = _FIDELITY_ROW
    mat = ps.IsotropicMaterial(E=210e9, nu=NU_SWEEP, rho=7800.0)
    geom = ps.make_geometry(r0_2b, two_T_pi)
    if part == 1:
        solver = ps.OutOfPlaneSolver(geom, mat, M=80, n_quad=30,
                                      boundary=ps.FreeFreeOOP())
    else:
        solver = ps.InPlaneSolver(geom, mat, M=80, n_quad=30,
                                   boundary=ps.FreeFreeIP())
    sel, cnt, how = _derive_branches_robust(
        solver.fast, Om, n_dofs, xmax, full_search, track, select_fill,
        solver=solver, sigma_min_from_K=sigma_min_from_K)
    ok = (sel is not None) and (cnt >= n_dofs) and (how == "cold")
    return ok, (f"row={tag} r0_2b={r0_2b:g} 2T/pi={two_T_pi:g} part={part} "
                f"Om={Om:.6f}: cnt={cnt}/{n_dofs} how={how!r}")


def _worker(args):
    tag, r0_2b, two_T_pi, part, n_dofs, xmax, Om, Om_lit = args
    t0 = time.time()
    try:
        import plate_solver as ps
        from plate_solver.detectors import (full_search, track, select_fill,
                                             sigma_min_from_K,
                                             weak_enforcement_residual_oop,
                                             weak_enforcement_residual_ip)

        mat = ps.IsotropicMaterial(E=210e9, nu=NU_SWEEP, rho=7800.0)
        geom = ps.make_geometry(r0_2b, two_T_pi)

        if part == 1:
            solver = ps.OutOfPlaneSolver(geom, mat, M=80, n_quad=30,
                                          boundary=ps.FreeFreeOOP())
        else:
            solver = ps.InPlaneSolver(geom, mat, M=80, n_quad=30,
                                       boundary=ps.FreeFreeIP())

        sel, cnt, how = _derive_branches_robust(
            solver.fast, Om, n_dofs, xmax, full_search, track, select_fill,
            solver=solver, sigma_min_from_K=sigma_min_from_K)

        if sel is None:
            return dict(ok=False, tag=tag, r0_2b=r0_2b, two_T_pi=two_T_pi,
                        part=part, Om=Om, Om_lit=Om_lit, cnt=cnt, how=how,
                        error=f"still only {cnt}/{n_dofs} even after "
                              f"_derive_branches_robust ({how})",
                        dt=time.time() - t0)

        if part == 1:
            res = weak_enforcement_residual_oop(solver, Om, sel, n_dofs=n_dofs)
        else:
            res = weak_enforcement_residual_ip(solver, Om, sel, n_dofs=n_dofs)

        if not res.get("ok", False):
            return dict(ok=False, tag=tag, r0_2b=r0_2b, two_T_pi=two_T_pi,
                        part=part, Om=Om, Om_lit=Om_lit, cnt=cnt, how=how,
                        error=f"branches filled ({how}) but residual screen "
                              f"failed: {res.get('reason')}",
                        dt=time.time() - t0)

        s = solver.sigma_min(Om, sel)
        out = dict(ok=True, tag=tag, r0_2b=r0_2b, two_T_pi=two_T_pi,
                   part=part, Om=Om, Om_lit=Om_lit, block=res["block"],
                   cnt=cnt, how=how, log10_sigma_min=float(s),
                   dt=time.time() - t0)
        if part == 1:
            out.update(rV=res["rV"], rM=res["rM"], maxW=res["maxW"],
                       verdict=res["verdict"])
        else:
            extra = {k: v for k, v in res.items() if k not in ("ok", "block")}
            out.update(extra)
        return out
    except Exception as exc:
        return dict(ok=False, tag=tag, r0_2b=r0_2b, two_T_pi=two_T_pi,
                    part=part, Om=Om, Om_lit=Om_lit, error=f"{exc}",
                    dt=time.time() - t0)


def main():
    print("=" * 78)
    print("  probe_ff_residual_screen_reconstruction_fix -- start "
          + time.strftime("%Y-%m-%d %H:%M:%S"))
    print("=" * 78, flush=True)

    import plate_solver as ps
    if ps.SOLVER_VERSION != EXPECT_VER:
        print(f"PREFLIGHT FAIL: SOLVER_VERSION={ps.SOLVER_VERSION!r} != "
              f"{EXPECT_VER!r}")
        raise SystemExit(2)
    print(f"preflight OK: SOLVER_VERSION={ps.SOLVER_VERSION}", flush=True)

    if len(FAILED) != 44:
        print(f"FATAL: matched {len(FAILED)}/44 expected failing rows from "
              f"TARGETS -- FAILED_KEYS may be stale or TARGETS changed "
              f"upstream. Stopping rather than processing a wrong set.")
        raise SystemExit(2)
    n_p1 = sum(1 for t in FAILED if t[3] == 1)
    n_p2 = sum(1 for t in FAILED if t[3] == 2)
    print(f"Matched {len(FAILED)}/44 failing rows against verbatim TARGETS "
          f"({n_p1} OOP + {n_p2} IP) -- fidelity check passed.\n",
          flush=True)

    print("--- In-run fidelity gate (copy of _derive_branches_robust "
          "vs. one already-succeeding case) ---", flush=True)
    gate_ok, gate_detail = _fidelity_gate()
    print(f"  {gate_detail}")
    if not gate_ok:
        print("FIDELITY GATE FAILED -- the copied _derive_branches_robust "
              "does not reproduce 'cold' success on an already-passing "
              "case. Stopping before processing the 44 failures; the copy "
              "needs to be re-diffed against probe_threshold_sensitivity.py.")
        raise SystemExit(2)
    print("  GATE OK -- proceeding.\n", flush=True)

    n_workers = int(os.environ.get("N_WORKERS", str(len(FAILED))))
    print(f"N_WORKERS={n_workers}  NU_SWEEP={NU_SWEEP}", flush=True)

    results = []
    with ProcessPoolExecutor(max_workers=min(n_workers, len(FAILED))) as ex:
        futs = {ex.submit(_worker, t): t for t in FAILED}
        done = 0
        for fut in as_completed(futs):
            r = fut.result()
            done += 1
            partlabel = "OOP" if r["part"] == 1 else "IP "
            lit = f"{r['Om_lit']:.3f}" if r.get("Om_lit") is not None else "n/a"
            if r["ok"]:
                if r["part"] == 1:
                    print(f"  [{done}/44] {partlabel} r0_2b={r['r0_2b']:g} "
                          f"2T/pi={r['two_T_pi']:g} Om={r['Om']:.6f} "
                          f"(lit={lit}) FIXED via how={r['how']!r} "
                          f"cnt={r['cnt']} block={r['block']:>3s} "
                          f"rV={r['rV']:.4g} rM={r['rM']:.4g} "
                          f"maxW={r['maxW']:.4g} verdict={r['verdict']} "
                          f"[{r['dt']:.0f}s]", flush=True)
                else:
                    extra_str = " ".join(
                        f"{k}={v:.4g}" if isinstance(v, float) else f"{k}={v}"
                        for k, v in r.items()
                        if k not in ("ok", "tag", "r0_2b", "two_T_pi", "part",
                                     "Om", "Om_lit", "dt", "log10_sigma_min",
                                     "block", "cnt", "how"))
                    print(f"  [{done}/44] {partlabel} r0_2b={r['r0_2b']:g} "
                          f"2T/pi={r['two_T_pi']:g} Om={r['Om']:.6f} "
                          f"(lit={lit}) FIXED via how={r['how']!r} "
                          f"cnt={r['cnt']} {extra_str} "
                          f"(NO verdict bar -- IP exploratory) "
                          f"[{r['dt']:.0f}s]", flush=True)
            else:
                print(f"  [{done}/44] {partlabel} r0_2b={r['r0_2b']:g} "
                      f"2T/pi={r['two_T_pi']:g} Om={r['Om']:.6f} "
                      f"STILL FAILED: {r.get('error')} [{r['dt']:.0f}s]",
                      flush=True)
            results.append(r)

    ok_results = [r for r in results if r["ok"]]
    still_failed = [r for r in results if not r["ok"]]
    oop_ok = [r for r in ok_results if r["part"] == 1]
    ip_ok = [r for r in ok_results if r["part"] == 2]
    real_like = [r for r in oop_ok if r.get("verdict") == "REAL-like"]
    artifact_like = [r for r in oop_ok if r.get("verdict") == "ARTIFACT-like"]
    ambiguous = [r for r in oop_ok if r.get("verdict") == "AMBIGUOUS"]

    print("\n" + "=" * 78)
    print(f"SUMMARY: {len(ok_results)}/44 now resolved via "
          f"_derive_branches_robust, {len(still_failed)} still failed")
    print(f"  OOP (part 1): {len(oop_ok)}/12 resolved -- "
          f"{len(real_like)} REAL-like, {len(artifact_like)} ARTIFACT-like, "
          f"{len(ambiguous)} AMBIGUOUS")
    print(f"  IP  (part 2): {len(ip_ok)}/32 resolved -- raw ratios only, "
          f"no verdict bar")
    if still_failed:
        print("\nSTILL FAILED (root cause is now genuinely open, not an "
              "implementation gap, for these specific points):")
        for r in still_failed:
            print(f"  r0_2b={r['r0_2b']:g} 2T/pi={r['two_T_pi']:g} "
                  f"part={r['part']} Om={r['Om']:.6f}: {r.get('error')}")
    print("=" * 78, flush=True)

    print("\nDiagnostic/data-generation only -- no SOLVER_VERSION action, "
          "no package changes. Read the module docstring's pre-registered "
          "interpretation before folding any of this into the paper.",
          flush=True)
    print("End: " + time.strftime("%Y-%m-%d %H:%M:%S"), flush=True)


if __name__ == "__main__":
    t0 = time.time()
    main()
    print(f"\ntotal wall time: {(time.time() - t0)/60:.1f} min", flush=True)
