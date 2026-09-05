# -*- coding: utf-8 -*-
"""
probe_ff_ip_relabeled_basis_s10v2.py -- DIAGNOSTIC/DATA-GENERATION ONLY.
Direct follow-up to job 2417234 (probe_ff_residual_screen_geomsweep_
s10v2_reconstruction_fix.py), which just ran under the current s10
checkpoint. Read that job's summary before this docstring: 19/44
originally-failing candidates were fixed via `_derive_branches_robust`'s
escalation ladder -- ALL 12 OOP resolved (5 REAL-like, 7 ARTIFACT-like,
0 AMBIGUOUS), but only 7/32 IP resolved. The remaining 25 IP candidates
split into two genuinely different failure signatures, not one -- EXACTLY
the same two-group structure the pre-s10 investigation found at this
identical stage (probe_ff_ip_relabeled_basis.py, job 2327405):

  GROUP A (24 points): stuck at EXACTLY n_dofs-2 branches (18/20) even
  after escalating xmax all the way to 90.0 with a 2x-denser Newton grid
  -- job 2417234's own per-point log shows this exact "-2" shortfall on
  every single one of these 24, never a different amount (confirmed by
  grepping every "STILL FAILED" line, not eyeballed). This is the SAME
  signature the pre-s10 investigation already has an established, tested
  explanation for (LESSONS_LEARNED.md Sec 21.7): forcing select_fill to
  hit an exact target count reliably pulls in a worse root and degrades
  sigma_min. Given the escalation ladder still can't find a 20th branch
  for these 24 after trying up to xmax=90 -- under BOTH checkpoints, not
  just pre-s10 -- the more likely explanation is that 18 (n_dofs-2) IS
  the true achievable basis at this (Omega, geometry), the same
  already-published "used@28"/"used@36" basis-relabeling precedent
  PAPER1_FREEFREE_DRAFT.tex already applies elsewhere (Sec 4's
  convergence tables, and Appendix C.1's IP-06 footnote). This probe
  applies that SAME precedent under s10: instead of continuing to demand
  n_dofs branches, it evaluates weak_enforcement_residual_ip AT THE
  ACHIEVED BASIS SIZE (n_dofs=cnt, whatever that reproducibly turns out
  to be), reporting raw ratios honestly at the true achieved size rather
  than leaving these 24 as unadjudicated failures.

  GROUP B (1 point): r0_2b=1.5, 2T/pi=0.5, Om=0.547221 -- job 2417234's
  own log: "still only 20/20 even after _derive_branches_robust
  (exhausted(cap_reached_but_shallow,smin=-2.681))". This reached the
  full n_dofs=20 at the guarded xmax_cap=90 stage, but was rejected
  because its verified sigma_min (-2.681) is above the shallow_floor=-3.0
  guard -- a SHALLOW-LOCK-ON case, not a capacity case, exactly like the
  pre-s10 run's own single Group-B point (r0_2b=1.5, 2T/pi=0.5,
  Om=0.547282, smin=-2.679) -- same geometry, same angle, an Omega 0.01%
  apart, an smin 0.07% apart. This is not a coincidence worth
  re-litigating; it gets the SAME seed_polish_modes treatment: independently
  re-converge the true zero near Om=0.547221 using a neighbour-clamped
  window built from this geometry's own sibling candidates in the s10v2
  manifest (same (tag, part) group, Omega-sorted), then evaluate the
  residual there.

PRE-REGISTERED INTERPRETATION (unchanged from the pre-s10 version of this
probe, since the mechanism being tested is the same mechanism, now under
a different checkpoint):
  - Group A: report every achieved (n_dofs=cnt) raw ratio plainly, no
    verdict forced (IP has none). The honest label for any that get
    folded into the paper is "evaluated at a smaller, reproducibly-
    achieved basis (n=cnt, not the originally-targeted n_dofs)" --
    mirroring the existing used@28/used@36 convention exactly.
  - Group A ALSO reports whether cnt reproduces job 2417234's exact -2
    shortfall (it should, since full_search/track/select_fill are
    deterministic given the same Omega/xmax/material) -- any point that
    does NOT reproduce would itself be a new, unexplained finding worth
    flagging plainly, not smoothed over.
  - Group B: if seed_polish_modes converges to a materially different
    Omega than 0.547221 (drift beyond ~0.1-0.5%, matching the scale
    already accepted as "small" in the pre-s10 run's own read), that
    supports the off-root-lock-on diagnosis directly. If it reconverges
    to essentially the SAME Omega and STILL comes back shallow/hard-to-
    fill, report plainly as a genuinely different, unexplained finding
    for this one point.

COST: mirrors probe_ff_ip_relabeled_basis.py's own budget almost exactly
(same 24+1 population size, same per-stage cost profile -- Group A
re-runs the same expensive escalation ladder each of these 24 points
already exhausted under job 2417234, ~1150-1300s each per that job's own
per-point times for these exact rows; Group B costs one
seed_polish_modes convergence, ~150-250s, plus one escalated branch
derivation). Total: ~24*1250s + ~500s =~ 8.5 CPU-hours, parallelized
across up to 25 workers, expect roughly 20-25 min wall time.
--time=01:00:00 is generous headroom, matching the pre-s10 job's own.

DRY-RUN VALIDATION (done before shipping, sandbox, 2026-08-13): loaded
all four figures_ff_{r150,r167,r200,r250}_s10v2/freefree_manifest.json
files, matched GROUP_A_KEYS (24) + GROUP_B_KEY (1) against a fresh
reload: 25/25 matched, 0 missing -- exactly job 2417234's own 25 "STILL
FAILED" rows, split by their own reported shortfall signature (24 at
18/20, 1 at 20/20-but-shallow). `_derive_branches_best_effort` below was
diffed byte-for-byte against probe_ff_ip_relabeled_basis.py's own copy
and confirmed identical.

Requires figures_ff_{r150,r167,r200,r250}_s10v2/freefree_manifest.json
to already exist (relative to cwd) -- the same manifests jobs
2416866-2416869 and 2417234 used.
"""
import json
import os
import sys
import time
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed

sys.path.insert(0, os.environ.get("PKG_PATH", "."))
os.environ.setdefault("DPS", "40")

EXPECT_VER = os.environ.get("EXPECT_SOLVER_VERSION", "2026-07-10.s10")

R0_2B = {"r150": 1.5, "r167": 1.66667, "r200": 2.0, "r250": 2.5}
NU_SWEEP = 0.35
TOL = 0.002

# ----------------------------------------------------------------------
# The 25 IP candidates job 2417234 still could not reconstruct, split by
# their own reported shortfall signature -- copied verbatim from that
# job's "STILL FAILED" section (24 at "only 18/20", 1 at "only 20/20 ...
# shallow"). Matched against a fresh manifest reload before shipping:
# 25/25, 0 missing (see docstring).
# ----------------------------------------------------------------------
GROUP_A_KEYS = {
    ("r250", 2, 0.233789), ("r250", 2, 0.456229), ("r250", 2, 0.300044),
    ("r200", 2, 0.940803), ("r200", 2, 0.991518), ("r200", 2, 0.311588),
    ("r167", 2, 1.130124), ("r167", 2, 1.100605), ("r200", 2, 0.395938),
    ("r250", 2, 0.173827), ("r250", 2, 0.020249), ("r167", 2, 0.344884),
    ("r250", 2, 0.395618), ("r167", 2, 0.392964), ("r167", 2, 0.925598),
    ("r167", 2, 0.488406), ("r250", 2, 0.545505), ("r167", 2, 0.738343),
    ("r200", 2, 0.543646), ("r250", 2, 0.753426), ("r200", 2, 0.742430),
    ("r200", 2, 0.576614), ("r200", 2, 0.597271), ("r167", 2, 0.639443),
}
GROUP_B_KEY = ("r150", 2, 0.547221)
N_EXPECTED_A = 24
N_EXPECTED_B = 1


def _load_targets_for_tag(tag):
    """Same data-driven extraction as the two prior s10v2 probes in this
    lineage."""
    manifest_path = f"figures_ff_{tag}_s10v2/freefree_manifest.json"
    if not os.path.isfile(manifest_path):
        print(f"FATAL: manifest not found: {manifest_path} -- unzip "
              f"figures_ff_{tag}_s10v2.zip into that directory (relative "
              f"to cwd) before submitting this job.")
        raise SystemExit(2)

    manifest = json.load(open(manifest_path))
    targets = []  # (tag, key, part, Om, n_dofs, xmax, Om_lit)

    for key, entry in manifest.items():
        if key.startswith("_"):
            continue
        diag = entry["diag"]
        n_cutoffs = len(diag["cut_offs"])
        accepted = sorted(diag["accepted"], key=lambda a: a["Omega"])
        n_take = min(n_cutoffs, len(accepted))
        if n_take == 0:
            continue

        part = diag["part"]
        n_dofs = diag["n_dofs"]
        xmax = diag["max_dim"]
        raw_list = entry.get("raw", [])
        lit_list = entry.get("omega_lit", [])
        raw_to_lit = dict(zip(raw_list, lit_list))

        for a in accepted[:n_take]:
            om = a["Omega"]
            om_lit = None
            for raw_om, lit in raw_to_lit.items():
                if abs(raw_om - om) <= TOL * max(1, om):
                    om_lit = lit
                    break
            targets.append((tag, key, part, om, n_dofs, xmax, om_lit))

    return targets


def _load_all_targets():
    all_targets = []
    for tag in ("r150", "r167", "r200", "r250"):
        all_targets.extend(_load_targets_for_tag(tag))
    return all_targets


def _select_groups(all_targets):
    group_a = [t for t in all_targets
               if (t[0], t[2], round(t[3], 6)) in GROUP_A_KEYS]
    group_b = [t for t in all_targets
               if (t[0], t[2], round(t[3], 6)) == GROUP_B_KEY]
    return group_a, group_b


# ------------------------------------------------------------------------
# _derive_branches_best_effort -- copied VERBATIM from
# probe_ff_ip_relabeled_basis.py (which is itself a deliberate variant of
# _derive_branches_robust: on exhaustion it returns the BEST partial fill
# seen across every stage instead of discarding it). Do NOT hand-edit
# without re-diffing against that source.
# ------------------------------------------------------------------------
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
    tag, key, part, Om, n_dofs, xmax, Om_lit = args
    t0 = time.time()
    try:
        import plate_solver as ps
        from plate_solver.detectors import (full_search, track, select_fill,
                                             sigma_min_from_K,
                                             weak_enforcement_residual_ip)

        mat = ps.IsotropicMaterial(E=210e9, nu=NU_SWEEP, rho=7800.0)
        two_T_pi = float(key.rsplit(" ", 1)[1])
        geom = ps.make_geometry(R0_2B[tag], two_T_pi)
        solver = ps.InPlaneSolver(geom, mat, M=80, n_quad=30,
                                   boundary=ps.FreeFreeIP())

        sel, cnt, how = _derive_branches_best_effort(
            solver.fast, Om, n_dofs, xmax, full_search, track, select_fill,
            solver=solver, sigma_min_from_K=sigma_min_from_K)

        if sel is None or cnt == 0:
            return dict(ok=False, tag=tag, key=key, Om=Om, Om_lit=Om_lit,
                        target_n_dofs=n_dofs,
                        error=f"literally no branches found ({how})",
                        dt=time.time() - t0)

        achieved_n_dofs = cnt   # relabel: evaluate at what was ACTUALLY achieved.
        res = weak_enforcement_residual_ip(solver, Om, sel, n_dofs=achieved_n_dofs)
        if not res.get("ok", False):
            return dict(ok=False, tag=tag, key=key, Om=Om, Om_lit=Om_lit,
                        target_n_dofs=n_dofs, achieved_n_dofs=achieved_n_dofs,
                        how=how,
                        error=f"residual screen failed at achieved basis: "
                              f"{res.get('reason')}",
                        dt=time.time() - t0)

        reproduced_2417234 = (achieved_n_dofs == n_dofs - 2)
        return dict(ok=True, tag=tag, key=key, Om=Om, Om_lit=Om_lit,
                    target_n_dofs=n_dofs, achieved_n_dofs=achieved_n_dofs,
                    how=how, reproduced_2417234=reproduced_2417234,
                    block=res["block"], maxDisp=res["maxDisp"],
                    rTyy=res["rTyy"], rTyr=res["rTyr"],
                    dt=time.time() - t0)
    except Exception as exc:
        return dict(ok=False, tag=tag, key=key, Om=Om, Om_lit=Om_lit,
                    error=f"{exc}", traceback=traceback.format_exc(),
                    dt=time.time() - t0)


def _worker_group_b(args):
    """The single shallow-lock-on point: re-converge via seed_polish_modes
    using a neighbour window built from this geometry's own sibling
    candidates in the s10v2 manifest (same (tag, part) group,
    Omega-sorted)."""
    tag, key, part, Om, n_dofs, xmax, Om_lit = args
    t0 = time.time()
    try:
        import plate_solver as ps
        from plate_solver.detectors import (seed_polish_modes, full_search,
                                             track, select_fill,
                                             sigma_min_from_K,
                                             weak_enforcement_residual_ip)

        siblings_all = _load_targets_for_tag(tag)
        siblings = sorted(
            t[3] for t in siblings_all if t[1] == key and t[2] == part
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
        two_T_pi = float(key.rsplit(" ", 1)[1])
        geom = ps.make_geometry(R0_2B[tag], two_T_pi)
        solver = ps.InPlaneSolver(geom, mat, M=80, n_quad=30,
                                   boundary=ps.FreeFreeIP())

        got = seed_polish_modes(solver, [Om], n_dofs=n_dofs, max_dim=xmax,
                                 part=2, scan_lo=lo, scan_hi=hi, iters=12,
                                 verbose=False)
        Om_star = got[0] if got else None
        if Om_star is None:
            return dict(ok=False, tag=tag, key=key, Om=Om, Om_lit=Om_lit,
                        window=(lo, hi),
                        error="seed_polish_modes returned nothing",
                        dt=time.time() - t0)

        drift_pct = 100.0 * abs(Om_star - Om) / Om

        sel, cnt, how = _derive_branches_best_effort(
            solver.fast, Om_star, n_dofs, xmax, full_search, track,
            select_fill, solver=solver, sigma_min_from_K=sigma_min_from_K)
        if sel is None or cnt == 0:
            return dict(ok=False, tag=tag, key=key, Om=Om, Om_lit=Om_lit,
                        Om_star=Om_star, drift_pct=drift_pct, window=(lo, hi),
                        error=f"no branches at re-converged Om_star ({how})",
                        dt=time.time() - t0)

        res = weak_enforcement_residual_ip(solver, Om_star, sel, n_dofs=cnt)
        if not res.get("ok", False):
            return dict(ok=False, tag=tag, key=key, Om=Om, Om_lit=Om_lit,
                        Om_star=Om_star, drift_pct=drift_pct,
                        achieved_n_dofs=cnt, how=how,
                        error=f"residual screen failed: {res.get('reason')}",
                        dt=time.time() - t0)

        return dict(ok=True, tag=tag, key=key, Om=Om, Om_lit=Om_lit,
                    Om_star=Om_star, drift_pct=drift_pct, window=(lo, hi),
                    achieved_n_dofs=cnt, target_n_dofs=n_dofs, how=how,
                    block=res["block"], maxDisp=res["maxDisp"],
                    rTyy=res["rTyy"], rTyr=res["rTyr"],
                    dt=time.time() - t0)
    except Exception as exc:
        return dict(ok=False, tag=tag, key=key, Om=Om, Om_lit=Om_lit,
                    error=f"{exc}", traceback=traceback.format_exc(),
                    dt=time.time() - t0)


def main():
    print("=" * 78)
    print("  probe_ff_ip_relabeled_basis_s10v2 -- start "
          + time.strftime("%Y-%m-%d %H:%M:%S"))
    print("=" * 78, flush=True)

    import plate_solver as ps
    if ps.SOLVER_VERSION != EXPECT_VER:
        print(f"PREFLIGHT FAIL: SOLVER_VERSION={ps.SOLVER_VERSION!r} != "
              f"{EXPECT_VER!r}")
        raise SystemExit(2)
    print(f"preflight OK: SOLVER_VERSION={ps.SOLVER_VERSION}", flush=True)

    all_targets = _load_all_targets()
    group_a, group_b = _select_groups(all_targets)
    if len(group_a) != N_EXPECTED_A or len(group_b) != N_EXPECTED_B:
        print(f"FATAL: matched {len(group_a)}/{N_EXPECTED_A} Group-A and "
              f"{len(group_b)}/{N_EXPECTED_B} Group-B rows against a fresh "
              f"manifest reload -- key lists may be stale, or the "
              f"manifests changed since job 2417234 ran.")
        raise SystemExit(2)
    print(f"Matched {N_EXPECTED_A}/{N_EXPECTED_A} Group-A + "
          f"{N_EXPECTED_B}/{N_EXPECTED_B} Group-B rows against a fresh "
          f"manifest reload.\n", flush=True)

    jobs = [("A", row) for row in group_a] + [("B", row) for row in group_b]
    n_total = len(jobs)
    n_workers = int(os.environ.get("N_WORKERS", str(n_total)))
    print(f"N_WORKERS={n_workers}  NU_SWEEP={NU_SWEEP}", flush=True)

    results = []
    with ProcessPoolExecutor(max_workers=min(n_workers, n_total)) as ex:
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
                    print(f"  [{done}/{n_total}] GRP-A {r['tag']} {r['key']} "
                          f"Om={r['Om']:.6f} achieved n_dofs="
                          f"{r['achieved_n_dofs']}/{r['target_n_dofs']} "
                          f"(reproduced 2417234's -2 shortfall: "
                          f"{r['reproduced_2417234']}) how={r['how']!r} "
                          f"rTyy={r['rTyy']:.4g} rTyr={r['rTyr']:.4g} "
                          f"maxDisp={r['maxDisp']:.4g} [{r['dt']:.0f}s]",
                          flush=True)
                else:
                    print(f"  [{done}/{n_total}] GRP-B {r['tag']} {r['key']} "
                          f"Om={r['Om']:.6f} -> Om_star={r['Om_star']:.6f} "
                          f"(drift={r['drift_pct']:.3f}%) achieved n_dofs="
                          f"{r['achieved_n_dofs']}/{r['target_n_dofs']} "
                          f"how={r['how']!r} rTyy={r['rTyy']:.4g} "
                          f"rTyr={r['rTyr']:.4g} maxDisp={r['maxDisp']:.4g} "
                          f"[{r['dt']:.0f}s]", flush=True)
            else:
                print(f"  [{done}/{n_total}] GRP-{grp} {r['tag']} {r['key']} "
                      f"Om={r['Om']:.6f} FAILED: {r.get('error')} "
                      f"[{r['dt']:.0f}s]", flush=True)
            results.append(r)

    ok_a = [r for r in results if r["group"] == "A" and r["ok"]]
    ok_b = [r for r in results if r["group"] == "B" and r["ok"]]
    failed = [r for r in results if not r["ok"]]

    print("\n" + "=" * 78)
    print(f"SUMMARY: Group A {len(ok_a)}/{N_EXPECTED_A} evaluated at "
          f"achieved basis, Group B {len(ok_b)}/{N_EXPECTED_B} "
          f"re-converged, {len(failed)} total failures")
    if ok_a:
        not_reproduced = [r for r in ok_a if not r["reproduced_2417234"]]
        print(f"  Group A: {len(ok_a)-len(not_reproduced)}/{len(ok_a)} "
              f"reproduced job 2417234's exact -2 shortfall "
              f"(achieved_n_dofs == target-2); "
              f"{len(not_reproduced)} did NOT (flagged below if any)")
        for r in not_reproduced:
            print(f"    NON-REPRODUCING: {r['tag']} {r['key']} "
                  f"Om={r['Om']:.6f} achieved={r['achieved_n_dofs']}/"
                  f"{r['target_n_dofs']}")
        ry = sorted(r["rTyy"] for r in ok_a)
        print(f"  Group A rTyy range at achieved (relabeled) basis: "
              f"[{ry[0]:.4g}, {ry[-1]:.4g}]")
        below = sum(1 for r in ok_a if r["rTyy"] < 0.3)
        above = sum(1 for r in ok_a if r["rTyy"] > 4)
        mid = len(ok_a) - below - above
        print(f"  Group A bimodal classification: {below} bimodal(<0.3), "
              f"{above} bimodal(>4), {mid} middle-band(0.3-4)")
    if ok_b:
        r = ok_b[0]
        print(f"  Group B: Om {r['Om']:.6f} -> Om_star {r['Om_star']:.6f} "
              f"(drift {r['drift_pct']:.3f}%), achieved "
              f"{r['achieved_n_dofs']}/{r['target_n_dofs']}, "
              f"rTyy={r['rTyy']:.4g}")
    if failed:
        print("\nFAILURES:")
        for r in failed:
            print(f"  group={r['group']} {r['tag']} {r['key']} "
                  f"Om={r['Om']:.6f}: {r.get('error')}")
    print("=" * 78, flush=True)
    print("\nDiagnostic-only -- no SOLVER_VERSION action, no package changes.")
    print("End: " + time.strftime("%Y-%m-%d %H:%M:%S"), flush=True)


if __name__ == "__main__":
    t0 = time.time()
    main()
    print(f"\ntotal wall time: {(time.time() - t0)/60:.1f} min", flush=True)
