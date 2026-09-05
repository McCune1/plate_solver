# -*- coding: utf-8 -*-
"""
probe_ff_ip_relabeled_basis.py -- DIAGNOSTIC/DATA-GENERATION ONLY.
Direct follow-up to job 2327101 (probe_ff_residual_screen_reconstruction_
fix.py), which just ran. Read that job's summary before this docstring:
19/44 originally-failing candidates were fixed via `_derive_branches_
robust`'s escalation ladder -- ALL 12 OOP resolved (4 REAL-like, 6
ARTIFACT-like, 2 AMBIGUOUS), but only 7/32 IP resolved. The remaining 25
IP candidates split into two genuinely different failure signatures, not
one:

  GROUP A (24 points): stuck at EXACTLY n_dofs-2 branches (e.g. 18/20)
  even after escalating xmax all the way to 90.0 with a 2x-denser Newton
  grid -- job 2327101's own per-point log shows this exact "-2" shortfall
  on every single one of these 24, never a different amount. This is the
  SAME signature this project already has an established, tested
  explanation for: LESSONS_LEARNED.md Sec 21.7 found that FORCING
  select_fill to hit an exact target count reliably pulls in a worse root
  and degrades sigma_min, 32/32 times it was tried -- i.e. an
  under-fill can be protective, not a defect. Given the escalation ladder
  (already proven to recover genuine under-search failures, per job
  2327101's other 19 fixes) STILL can't find a 20th/28th/etc. branch for
  these 24 after trying up to xmax=90, the more likely explanation now is
  that 18 (or n_dofs-2 generally) IS the true achievable basis at this
  (Omega, geometry) -- exactly the already-published "used@28"/"used@36"
  basis-relabeling precedent PAPER1_FREEFREE_DRAFT.tex already applies
  elsewhere (Sec 4's convergence tables). This probe applies that SAME
  precedent here: instead of continuing to demand n_dofs branches, it
  evaluates weak_enforcement_residual_ip AT THE ACHIEVED BASIS SIZE
  (n_dofs=cnt, whatever that reproducibly turns out to be), reporting raw
  ratios honestly at the true achieved size rather than leaving these 24
  as unadjudicated failures.

  GROUP B (1 point): r0_2b=1.5, 2T/pi=0.5, Om=0.547282 -- this one
  actually DID reach the full n_dofs=20 at the guarded xmax_cap=90 stage,
  but was rejected because its verified sigma_min (-2.679) is above the
  shallow_floor=-3.0 guard, i.e. this is a SHALLOW-LOCK-ON case, not a
  capacity case: the 20-branch fill found at Om=0.547282 likely
  corresponds to an adjacent, weaker feature rather than the genuine
  determinant zero at that literature Omega. This is exactly the
  off-root-evaluation failure mode probe_ip_residual_calibration.py was
  built to avoid (LESSONS_LEARNED Sec 21.6) -- so this ONE point gets the
  seed_polish_modes treatment instead of more branch-derivation
  escalation: independently re-converge the true zero near Om=0.547282
  using a neighbour-clamped window built from this geometry's OWN sibling
  candidates in TARGETS (the same (r0_2b, two_T_pi, part) group, already
  sorted by construction), then evaluate the residual there.

PRE-REGISTERED INTERPRETATION:
  - Group A: report every achieved (n_dofs=cnt) raw ratio plainly, no
    verdict forced (IP has none). If a future session wants to fold these
    into Sec 6.4, the honest label is "evaluated at a smaller,
    reproducibly-achieved basis (n=cnt, not the originally-targeted
    n_dofs)" -- mirroring the existing used@28/used@36 convention exactly,
    not a new kind of caveat.
  - Group A ALSO reports whether cnt is reproducible run-to-run (it
    should be, since full_search/track/select_fill are deterministic
    given the same Omega/xmax/material) -- if this single point's escalation
    reproduces a DIFFERENT cnt than job 2327101 logged, that is itself a
    finding worth flagging (would suggest non-determinism somewhere in the
    pipeline, which this project has never previously observed and would
    be a genuinely new and important structural finding).
  - Group B: if seed_polish_modes converges to a materially different
    Omega than 0.547282 (i.e. drift beyond ~0.1-0.5%, the scale of drift
    already seen as "small/acceptable" in probe_ip_residual_calibration.py's
    own sanity check), that supports the off-root-lock-on diagnosis
    directly -- report the new Om_star, its achieved fill, and its
    residual. If it reconverges to essentially the SAME Omega and STILL
    comes back shallow/hard-to-fill, that is a genuinely different,
    unexplained finding for this one point -- report plainly, do not force
    a story onto it.

COST: 24 Group-A points re-run the SAME expensive escalation ladder job
2327101 already paid for (deterministic, so cost should closely match
that job's own per-point times for these exact rows -- 1150-1300s each
per that log) plus one relatively cheap residual-screen evaluation at the
achieved basis. 1 Group-B point costs one seed_polish_modes convergence
(~150-250s, per this project's standard cost model) plus escalated branch
derivation and residual evaluation. Total: ~24*1250s + ~500s =~ 8.5
CPU-hours, parallelized across 25 workers, expect roughly 20-25 min wall
time. --time=01:00:00 is generous headroom.

Env var RATIO_TAG is exported only to satisfy the imported
probe_ff_residual_screen_geomsweep module's own CLI guard (unused by this
probe's own logic) -- same as probe_ff_residual_screen_reconstruction_
fix.py. Run from FutureWork/geometry_sweep/ (same directory as that file).
"""
import os
import sys
import time
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed

sys.path.insert(0, os.environ.get("PKG_PATH", "."))
os.environ.setdefault("DPS", "40")

EXPECT_VER = os.environ.get("EXPECT_SOLVER_VERSION", "2026-07-10.s8")
NU_SWEEP = 0.35   # see probe_ff_residual_screen_geomsweep.py's own note.

try:
    from probe_ff_residual_screen_geomsweep import TARGETS as _ALL_TARGETS
except ImportError:
    print("FATAL: could not import TARGETS from probe_ff_residual_screen_"
          "geomsweep.py -- run this script from FutureWork/geometry_sweep/.")
    raise SystemExit(2)

# The 24 Group-A keys (stuck at exactly n_dofs-2 even at xmax_cap=90),
# extracted directly from job 2327101's "STILL FAILED" list -- excludes
# the 1 Group-B shallow-lock-on point (Om=0.547282), handled separately
# below via its own TARGETS-group neighbour window.
GROUP_A_KEYS = {
    (2.5, 0.5, 2, 0.173706), (2.5, 0.25, 2, 0.456193),
    (2.0, 0.5, 2, 0.395857), (2.5, 0.5, 2, 0.110223),
    (2.0, 0.5, 2, 0.311546), (2.5, 0.25, 2, 0.020274),
    (2.5, 0.25, 2, 0.110248), (1.66667, 0.5, 2, 0.344895),
    (1.66667, 0.5, 2, 0.392910), (2.0, 0.25, 2, 0.177806),
    (1.66667, 0.25, 2, 0.925569), (2.5, 0.5, 2, 0.022797),
    (1.66667, 0.25, 2, 0.260820), (1.66667, 0.5, 2, 0.260836),
    (2.5, 0.5, 2, 0.008006), (2.0, 0.25, 2, 0.012306),
    (1.66667, 0.25, 2, 0.017751), (1.66667, 0.5, 2, 0.017750),
    (2.0, 0.5, 2, 0.012412), (2.5, 0.25, 2, 0.008006),
    (2.0, 0.25, 2, 0.742474), (2.0, 0.5, 2, 0.543552),
    (1.66667, 0.25, 2, 0.738254), (2.0, 0.25, 2, 0.597162),
}
GROUP_B_KEY = (1.5, 0.5, 2, 0.547282)

GROUP_A = [row for row in _ALL_TARGETS
           if (row[1], row[2], row[3], round(row[6], 6)) in GROUP_A_KEYS]
GROUP_B = [row for row in _ALL_TARGETS
           if (row[1], row[2], row[3], round(row[6], 6)) == GROUP_B_KEY]


# ---------------------------------------------------------------------
# _derive_branches_best_effort -- a DELIBERATE VARIANT of
# _derive_branches_robust (probe_threshold_sensitivity.py, verbatim
# stages), modified ONLY in what happens on exhaustion: the original
# discards the partial fill and returns (None, cnt, "exhausted"); this
# variant tracks and returns the BEST (largest-cnt) partial fill seen
# across every stage instead, since Group A's whole point is to use the
# best achievable fill, not to require full n_dofs. Every stage's search
# logic is otherwise identical to the source -- do not hand-edit without
# re-diffing.
# ---------------------------------------------------------------------
def _derive_branches_best_effort(eng, Omega, n_dofs, max_dim, full_search,
                                  track, select_fill, xmax_cap=90.0,
                                  solver=None, sigma_min_from_K=None,
                                  shallow_floor=-3.0):
    best_sel, best_cnt, best_how = None, 0, "none"

    raw = full_search(eng, Omega, xmax=max_dim)
    sel, cnt = select_fill(raw, n_dofs)
    if cnt >= n_dofs:
        return sel, cnt, "cold"
    if cnt > best_cnt:
        best_sel, best_cnt, best_how = sel, cnt, "cold(best-effort)"

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
            if cnt > best_cnt:
                best_sel, best_cnt = sel, cnt
                best_how = f"tracked({sign*frac:+.2f})(best-effort)"

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
        if cnt > best_cnt:
            best_sel, best_cnt = sel, cnt
            best_how = (f"escalated(xmax={xmax_try:.1f},ngrid={ngrid_try},"
                        f"nax={nax_try})(best-effort)")

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
            # Full fill found but shallow -- still record it as the best
            # partial candidate (Group B's own case), caller decides how
            # to treat "full-but-shallow" vs "best-effort-partial".
            return sel, cnt, (f"xmax_cap(reached_but_shallow,smin={smin:.3f})")
        return sel, cnt, f"xmax_cap(unverified,xmax={xmax_cap:.1f})"
    if cnt > best_cnt:
        best_sel, best_cnt = sel, cnt
        best_how = f"xmax_cap(xmax={xmax_cap:.1f},ngrid=20,nax=480)(best-effort)"

    return best_sel, best_cnt, best_how


def _worker_group_a(args):
    tag, r0_2b, two_T_pi, part, n_dofs, xmax, Om, Om_lit = args
    t0 = time.time()
    try:
        import plate_solver as ps
        from plate_solver.detectors import (full_search, track, select_fill,
                                             sigma_min_from_K,
                                             weak_enforcement_residual_ip)

        mat = ps.IsotropicMaterial(E=210e9, nu=NU_SWEEP, rho=7800.0)
        geom = ps.make_geometry(r0_2b, two_T_pi)
        solver = ps.InPlaneSolver(geom, mat, M=80, n_quad=30,
                                   boundary=ps.FreeFreeIP())

        sel, cnt, how = _derive_branches_best_effort(
            solver.fast, Om, n_dofs, xmax, full_search, track, select_fill,
            solver=solver, sigma_min_from_K=sigma_min_from_K)

        if sel is None or cnt == 0:
            return dict(ok=False, tag=tag, r0_2b=r0_2b, two_T_pi=two_T_pi,
                        Om=Om, Om_lit=Om_lit, target_n_dofs=n_dofs,
                        error=f"literally no branches found ({how})",
                        dt=time.time() - t0)

        achieved_n_dofs = cnt   # relabel: evaluate at what was ACTUALLY achieved.
        res = weak_enforcement_residual_ip(solver, Om, sel, n_dofs=achieved_n_dofs)
        if not res.get("ok", False):
            return dict(ok=False, tag=tag, r0_2b=r0_2b, two_T_pi=two_T_pi,
                        Om=Om, Om_lit=Om_lit, target_n_dofs=n_dofs,
                        achieved_n_dofs=achieved_n_dofs, how=how,
                        error=f"residual screen failed at achieved basis: "
                              f"{res.get('reason')}",
                        dt=time.time() - t0)

        reproduced_2327101 = (achieved_n_dofs == n_dofs - 2)
        return dict(ok=True, tag=tag, r0_2b=r0_2b, two_T_pi=two_T_pi,
                    Om=Om, Om_lit=Om_lit, target_n_dofs=n_dofs,
                    achieved_n_dofs=achieved_n_dofs, how=how,
                    reproduced_2327101=reproduced_2327101,
                    block=res["block"], maxDisp=res["maxDisp"],
                    rTyy=res["rTyy"], rTyr=res["rTyr"],
                    dt=time.time() - t0)
    except Exception as exc:
        return dict(ok=False, tag=tag, r0_2b=r0_2b, two_T_pi=two_T_pi,
                    Om=Om, Om_lit=Om_lit, error=f"{exc}",
                    traceback=traceback.format_exc(), dt=time.time() - t0)


def _worker_group_b(args):
    """The single shallow-lock-on point: re-converge via seed_polish_modes
    using a neighbour window built from this geometry's own TARGETS
    siblings (same (r0_2b, two_T_pi, part) group, already Omega-sorted)."""
    tag, r0_2b, two_T_pi, part, n_dofs, xmax, Om, Om_lit = args
    t0 = time.time()
    try:
        import plate_solver as ps
        from plate_solver.detectors import (seed_polish_modes, full_search,
                                             track, select_fill,
                                             sigma_min_from_K,
                                             weak_enforcement_residual_ip)

        siblings = sorted(
            row[6] for row in _ALL_TARGETS
            if row[1] == r0_2b and row[2] == two_T_pi and row[3] == part
        )
        idx = siblings.index(Om)
        W_SEED = 0.12
        lo, hi = Om * (1.0 - W_SEED), Om * (1.0 + W_SEED)
        if idx > 0:
            lo = max(lo, 0.5 * (siblings[idx - 1] + Om))
        if idx < len(siblings) - 1:
            hi = min(hi, 0.5 * (Om + siblings[idx + 1]))
        lo = max(lo, 1e-4)

        mat = ps.IsotropicMaterial(E=210e9, nu=NU_SWEEP, rho=7800.0)
        geom = ps.make_geometry(r0_2b, two_T_pi)
        solver = ps.InPlaneSolver(geom, mat, M=80, n_quad=30,
                                   boundary=ps.FreeFreeIP())

        got = seed_polish_modes(solver, [Om], n_dofs=n_dofs, max_dim=xmax,
                                 part=2, scan_lo=lo, scan_hi=hi, iters=12,
                                 verbose=False)
        Om_star = got[0] if got else None
        if Om_star is None:
            return dict(ok=False, tag=tag, r0_2b=r0_2b, two_T_pi=two_T_pi,
                        Om=Om, Om_lit=Om_lit, window=(lo, hi),
                        error="seed_polish_modes returned nothing",
                        dt=time.time() - t0)

        drift_pct = 100.0 * abs(Om_star - Om) / Om

        sel, cnt, how = _derive_branches_best_effort(
            solver.fast, Om_star, n_dofs, xmax, full_search, track,
            select_fill, solver=solver, sigma_min_from_K=sigma_min_from_K)
        if sel is None or cnt == 0:
            return dict(ok=False, tag=tag, r0_2b=r0_2b, two_T_pi=two_T_pi,
                        Om=Om, Om_lit=Om_lit, Om_star=Om_star,
                        drift_pct=drift_pct, window=(lo, hi),
                        error=f"no branches at re-converged Om_star ({how})",
                        dt=time.time() - t0)

        res = weak_enforcement_residual_ip(solver, Om_star, sel, n_dofs=cnt)
        if not res.get("ok", False):
            return dict(ok=False, tag=tag, r0_2b=r0_2b, two_T_pi=two_T_pi,
                        Om=Om, Om_lit=Om_lit, Om_star=Om_star,
                        drift_pct=drift_pct, achieved_n_dofs=cnt, how=how,
                        error=f"residual screen failed: {res.get('reason')}",
                        dt=time.time() - t0)

        return dict(ok=True, tag=tag, r0_2b=r0_2b, two_T_pi=two_T_pi,
                    Om=Om, Om_lit=Om_lit, Om_star=Om_star,
                    drift_pct=drift_pct, window=(lo, hi),
                    achieved_n_dofs=cnt, target_n_dofs=n_dofs, how=how,
                    block=res["block"], maxDisp=res["maxDisp"],
                    rTyy=res["rTyy"], rTyr=res["rTyr"],
                    dt=time.time() - t0)
    except Exception as exc:
        return dict(ok=False, tag=tag, r0_2b=r0_2b, two_T_pi=two_T_pi,
                    Om=Om, Om_lit=Om_lit, error=f"{exc}",
                    traceback=traceback.format_exc(), dt=time.time() - t0)


def main():
    print("=" * 78)
    print("  probe_ff_ip_relabeled_basis -- start "
          + time.strftime("%Y-%m-%d %H:%M:%S"))
    print("=" * 78, flush=True)

    import plate_solver as ps
    if ps.SOLVER_VERSION != EXPECT_VER:
        print(f"PREFLIGHT FAIL: SOLVER_VERSION={ps.SOLVER_VERSION!r} != "
              f"{EXPECT_VER!r}")
        raise SystemExit(2)
    print(f"preflight OK: SOLVER_VERSION={ps.SOLVER_VERSION}", flush=True)

    if len(GROUP_A) != 24 or len(GROUP_B) != 1:
        print(f"FATAL: matched {len(GROUP_A)}/24 Group-A and "
              f"{len(GROUP_B)}/1 Group-B rows -- key lists may be stale.")
        raise SystemExit(2)
    print("Matched 24/24 Group-A + 1/1 Group-B rows against TARGETS.\n",
          flush=True)

    jobs = [("A", row) for row in GROUP_A] + [("B", row) for row in GROUP_B]
    n_workers = int(os.environ.get("N_WORKERS", str(len(jobs))))
    print(f"N_WORKERS={n_workers}  NU_SWEEP={NU_SWEEP}", flush=True)

    results = []
    with ProcessPoolExecutor(max_workers=min(n_workers, len(jobs))) as ex:
        futs = {}
        for grp, row in jobs:
            fn = _worker_group_a if grp == "A" else _worker_group_b
            futs[ex.submit(fn, row)] = (grp, row)
        done = 0
        for fut in as_completed(futs):
            grp, row = futs[fut]
            r = fut.result()
            r["group"] = grp
            done += 1
            if r["ok"]:
                if grp == "A":
                    print(f"  [{done}/25] GRP-A r0_2b={r['r0_2b']:g} "
                          f"2T/pi={r['two_T_pi']:g} Om={r['Om']:.6f} "
                          f"achieved n_dofs={r['achieved_n_dofs']}/"
                          f"{r['target_n_dofs']} (reproduced 2327101's "
                          f"-2 shortfall: {r['reproduced_2327101']}) "
                          f"how={r['how']!r} rTyy={r['rTyy']:.4g} "
                          f"rTyr={r['rTyr']:.4g} maxDisp={r['maxDisp']:.4g} "
                          f"[{r['dt']:.0f}s]", flush=True)
                else:
                    print(f"  [{done}/25] GRP-B r0_2b={r['r0_2b']:g} "
                          f"2T/pi={r['two_T_pi']:g} Om={r['Om']:.6f} "
                          f"-> Om_star={r['Om_star']:.6f} "
                          f"(drift={r['drift_pct']:.3f}%) "
                          f"achieved n_dofs={r['achieved_n_dofs']}/"
                          f"{r['target_n_dofs']} how={r['how']!r} "
                          f"rTyy={r['rTyy']:.4g} rTyr={r['rTyr']:.4g} "
                          f"maxDisp={r['maxDisp']:.4g} [{r['dt']:.0f}s]",
                          flush=True)
            else:
                print(f"  [{done}/25] GRP-{grp} r0_2b={r['r0_2b']:g} "
                      f"2T/pi={r['two_T_pi']:g} Om={r['Om']:.6f} FAILED: "
                      f"{r.get('error')} [{r['dt']:.0f}s]", flush=True)
            results.append(r)

    ok_a = [r for r in results if r["group"] == "A" and r["ok"]]
    ok_b = [r for r in results if r["group"] == "B" and r["ok"]]
    failed = [r for r in results if not r["ok"]]

    print("\n" + "=" * 78)
    print(f"SUMMARY: Group A {len(ok_a)}/24 evaluated at achieved basis, "
          f"Group B {len(ok_b)}/1 re-converged, {len(failed)} total failures")
    if ok_a:
        not_reproduced = [r for r in ok_a if not r["reproduced_2327101"]]
        print(f"  Group A: {len(ok_a)-len(not_reproduced)}/{len(ok_a)} "
              f"reproduced job 2327101's exact -2 shortfall "
              f"(achieved_n_dofs == target-2); "
              f"{len(not_reproduced)} did NOT (flagged below if any)")
        for r in not_reproduced:
            print(f"    NON-REPRODUCING: r0_2b={r['r0_2b']:g} "
                  f"2T/pi={r['two_T_pi']:g} Om={r['Om']:.6f} achieved="
                  f"{r['achieved_n_dofs']}/{r['target_n_dofs']}")
        ry = sorted(r["rTyy"] for r in ok_a)
        print(f"  Group A rTyy range at achieved (relabeled) basis: "
              f"[{ry[0]:.4g}, {ry[-1]:.4g}]")
    if ok_b:
        r = ok_b[0]
        print(f"  Group B: Om {r['Om']:.6f} -> Om_star {r['Om_star']:.6f} "
              f"(drift {r['drift_pct']:.3f}%), achieved "
              f"{r['achieved_n_dofs']}/{r['target_n_dofs']}, "
              f"rTyy={r['rTyy']:.4g}")
    if failed:
        print("\nFAILURES:")
        for r in failed:
            print(f"  group={r['group']} r0_2b={r['r0_2b']:g} "
                  f"2T/pi={r['two_T_pi']:g} Om={r['Om']:.6f}: "
                  f"{r.get('error')}")
    print("=" * 78, flush=True)
    print("\nDiagnostic-only -- no SOLVER_VERSION action, no package changes.")
    print("End: " + time.strftime("%Y-%m-%d %H:%M:%S"), flush=True)


if __name__ == "__main__":
    t0 = time.time()
    main()
    print(f"\ntotal wall time: {(time.time() - t0)/60:.1f} min", flush=True)
