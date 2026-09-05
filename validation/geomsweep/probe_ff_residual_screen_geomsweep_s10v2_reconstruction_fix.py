# -*- coding: utf-8 -*-
"""
probe_ff_residual_screen_geomsweep_s10v2_reconstruction_fix.py --
DIAGNOSTIC/DATA-GENERATION ONLY. Writes only a plain-text log; no package
changes, no SOLVER_VERSION implications.

CONTEXT: jobs 2416866-2416869 (probe_ff_residual_screen_geomsweep_s10v2.py)
re-screened the cut-off-adjacent population (119 flexural, 96 extensional
candidates -- same population, same counts, as the pre-s10 tally already
in PAPER1_FREEFREE_DRAFT.tex's Table 13 caption / Appendix C.3) under the
current s10 checkpoint. 107/119 flexural and 64/96 extensional
reconstructed directly; 12 flexural + 32 extensional (44 total) came back
"only n/n_dofs branches filled" -- a bare cold-search failure, because
that probe's `_worker` (like the pre-s10 probe it was modeled on) makes
exactly ONE `select_fill(full_search(solver.fast, Om, xmax=xmax), n_dofs)`
call with no fallback.

This is NOT a new problem. The pre-s10 run hit the identical 44-candidate
failure (same keys, confirmed by LESSONS_LEARNED.md Sec 18.10's own
cross-check: the s10v2 run's 107/64 direct-reconstruction split matches
the pre-s10 run's own reported "107 reconstructed directly" / "64
reconstructed directly" figures almost exactly) and closed it with
`probe_ff_residual_screen_reconstruction_fix.py` (jobs on top of
2327073-2327076), which reruns exactly these failures through
`_derive_branches_robust` -- the SAME escalation ladder validated at 352/
352 points by probe_threshold_sensitivity.py (job 2325369) and already
used successfully elsewhere in this project (probe_ff_ip_relabeled_
basis.py). `_derive_branches_robust` is copied VERBATIM below from that
pre-s10 script; it is pure branch-reconstruction algorithm, with no
SOLVER_VERSION-sensitive content (it only calls full_search/track/
select_fill/sigma_min_from_K, all unchanged names/signatures under s10 --
confirmed by grepping the deployed plate_solver/detectors.py before
writing this probe). Porting it here closes the same gap under s10 that
it closed under s8, using the identical, already-proven method rather
than inventing a new one.

WHAT THIS PROBE DOES: for exactly the 44 candidates that failed in jobs
2416866-2416869 (FAILED_KEYS below, extracted directly from those four
jobs' own printed FAILURES sections -- not retyped from memory; every
(r0_2b, two_T_pi, part, Omega) key was copied verbatim from the job logs,
then independently re-matched against a fresh data-driven re-load of each
manifest to confirm all 44 resolve to exactly 44 targets with zero
missing, before this probe was finalized -- see the dry-run validation
note at the end of this docstring), reruns branch reconstruction using
`_derive_branches_robust` instead of the bare cold call the s10v2 probe
used. Any candidate that now fills to n_dofs runs the SAME production
residual screen (weak_enforcement_residual_oop/_ip, unmodified) the
s10v2 probe would have run. Any candidate that STILL doesn't fill -- even
after escalating to xmax_cap=90 with a 2x-denser Newton-seed grid -- is a
materially stronger negative result than "cold search failed," meaning
the already-proven-robust escalation ladder itself was exhausted at this
specific (Omega, geometry) under s10, not just an under-engineered first
attempt.

Targets are loaded the SAME data-driven way probe_ff_residual_screen_
geomsweep_s10v2.py loads them (each manifest's own diag.cut_offs/
diag.accepted/diag.n_dofs/diag.max_dim, not hand-transcribed), then
filtered down to the 44 FAILED_KEYS -- this avoids retyping n_dofs/xmax
by hand for 44 rows, the same transcription-risk this project's probes
consistently avoid.

IN-RUN FIDELITY CHECK (per plate-solver-cluster-probes skill Sec 3,
applied here even though `_derive_branches_robust` isn't a modified
package method, because it IS copied from another file rather than
imported): before touching any of the 44 failures, this probe reruns the
ORIGINAL cold-search step on one already-SUCCEEDING candidate (the first
target loaded that is NOT in FAILED_KEYS) and confirms `_derive_branches_
robust` returns how="cold" with a cnt that satisfies n_dofs on the very
first stage -- i.e. the copy doesn't change behavior on cases that
already worked. If that gate fails, the copy has a transcription bug and
this probe stops before processing anything else.

HONEST SCOPE CAVEAT -- carried forward from every probe in this
lineage: weak_enforcement_residual_oop's verdict bar was calibrated ONCE,
at ONE geometry (FF-P1, r0/2b=1.5, 2Theta=0.5pi) and ONE material
(nu=0.30), then applied unchanged here to nu=0.35 and 4 new radius
ratios. weak_enforcement_residual_ip has NO calibrated verdict bar at
all -- raw ratios only, reported by angle, never pooled.

PRE-REGISTERED INTERPRETATION (write before running, honor after):
  - If MOST/ALL of the 44 now fill successfully -> fold the completed
    verdicts into Table 13's caption and Appendix C.3, producing a true
    s10-vintage 119/96-style tally (not the pre-escalation spot-check
    currently in the text) and retiring the "spot-check, escalation not
    re-run" caveat added 2026-08-13.
  - If OOP resolves well but IP resolves less completely (or vice versa)
    -> report the split explicitly, exactly as Appendix C.3 already does
    for the pre-s10 run (49/54/4 direct + 4/6/2 via escalation for OOP;
    57+7+11 bimodal / 7+3+14 middle-band for IP) -- do not force parity
    between the two motion types if the data doesn't show it.
  - If a candidate STILL fails even at the guarded xmax_cap=90 stage ->
    report it plainly as a genuine, escalation-resistant reconstruction
    failure at that specific (Omega, geometry), not folded into either
    verdict count. Compare the escalation stage (`how`) and any pattern
    against r0/2b or 2Theta/pi to the concentration already noted in
    Appendix C.3's 2026-08-13 addendum (10/30 OOP failures at r250,
    flat 8/24 IP failures at every ratio) -- does the escalation ladder
    close the same points uniformly, or does the r250 concentration
    survive escalation too?
  - Do not average, smooth, or drop any of the 44 -- report every one,
    success or failure, exactly like both prior probes in this lineage.

COST: mirrors probe_ff_residual_screen_reconstruction_fix.py's own
budget exactly (same escalation ladder, same failure population size).
Per-point cost ranges from ~35-60s (points that resolve at the bare cold
stage -- shouldn't occur here since these 44 are drawn from confirmed
cold-search failures, but the ladder always tries cold first) up to
~830-1200s for a point that reaches the guarded xmax_cap=90 stage.
Budgeting conservatively at up to ~1200s/point x 44 points = ~14.7
CPU-hours serial; parallelized 1 point/worker via ProcessPoolExecutor
across up to 44 workers, expect a few minutes if most resolve at an
early stage, up to ~20-25 min if many need the expensive final stage.
--time=03:00:00, --partition=general matches the pre-s10 reconstruction-
fix job's own margin for this identical cost profile.

DRY-RUN VALIDATION (done before shipping, sandbox, 2026-08-13): loaded
all four figures_ff_{r150,r167,r200,r250}_s10v2/freefree_manifest.json
files (from FutureWork/geometry_sweep/_s10v2_manifests.zip, the same
manifests jobs 2416866-2416869 read) via the same _load_targets_for_tag
extraction below, then matched FAILED_KEYS against the freshly-loaded
targets: 44/44 matched, 0 missing, split 8 (r150) + 9 (r167) + 9 (r200)
+ 18 (r250) -- exactly the 12 OOP / 32 IP split reported in jobs
2416866-2416869's own FAILURES sections. No FAILED_KEYS entry was
invented or approximated.

Requires figures_ff_{r150,r167,r200,r250}_s10v2/freefree_manifest.json
to already exist (relative to cwd) -- the same manifests jobs
2416866-2416869 used. If any are missing, unzip the matching
figures_ff_{tag}_s10v2.zip first (see the submit script's preflight
check).
"""
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

sys.path.insert(0, os.environ.get("PKG_PATH", "."))
os.environ.setdefault("DPS", "40")

EXPECT_VER = os.environ.get("EXPECT_SOLVER_VERSION", "2026-07-10.s10")

R0_2B = {"r150": 1.5, "r167": 1.66667, "r200": 2.0, "r250": 2.5}
NU_SWEEP = 0.35  # matches geometry_sweep's actual material -- see every
                 # prior geomsweep-residual probe's own note; NOT the
                 # FF-P1 tab:oop/tab:ip benchmark's nu=0.30.
TOL = 0.002      # relative Omega match tolerance for the raw->omega_lit
                 # lookup (same as probe_ff_residual_screen_geomsweep_
                 # s10v2.py)

# ----------------------------------------------------------------------
# The 44 candidates that failed in jobs 2416866 (r150), 2416867 (r167),
# 2416868 (r200), 2416869 (r250) -- copied verbatim from each job's own
# FAILURES section (8 + 9 + 9 + 18 = 44, matching PAPER1_FREEFREE_DRAFT.
# tex Appendix C.3's "12 flexural, 32 extensional" pre-s10 count exactly).
# part: 1 = OOP (FF-P1), 2 = IP (FF-P2). Matched against a fresh
# data-driven target reload below by (tag, part, round(Omega, 6)) --
# verified 44/44 match, 0 missing, before this probe was finalized (see
# docstring).
# ----------------------------------------------------------------------
FAILED_KEYS = {
    "r150": {
        (2, 0.407776), (2, 0.448781), (2, 1.194872), (2, 1.031950),
        (2, 1.234253), (2, 0.831576), (2, 0.547221), (2, 0.697672),
    },
    "r167": {
        (2, 1.100605), (2, 1.130124), (2, 0.392964), (2, 0.344884),
        (2, 0.925598), (2, 0.488406), (2, 0.738343), (2, 0.639443),
        (1, 1.144852),
    },
    "r200": {
        (2, 0.395938), (2, 0.940803), (2, 0.311588), (2, 0.991518),
        (2, 0.742430), (2, 0.543646), (2, 0.576614), (2, 0.597271),
        (1, 1.468258),
    },
    "r250": {
        (2, 0.173827), (2, 0.233789), (2, 0.395618), (2, 0.300044),
        (2, 0.456229), (2, 0.020249), (2, 0.545505), (2, 0.753426),
        (1, 0.294511), (1, 0.166275), (1, 0.415583), (1, 0.125977),
        (1, 0.348585), (1, 0.510870), (1, 0.535811), (1, 0.689136),
        (1, 1.361344), (1, 1.054896),
    },
}
N_EXPECTED_FAILED = sum(len(v) for v in FAILED_KEYS.values())  # 44


def _load_targets_for_tag(tag):
    """Same data-driven extraction as probe_ff_residual_screen_geomsweep_
    s10v2.py's _load_targets, generalized to take the ratio tag as an
    argument instead of reading it from RATIO_TAG (this probe processes
    all four ratios in one run, unlike that one)."""
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


def _select_failed(all_targets):
    failed = []
    for t in all_targets:
        tag, key, part, om, n_dofs, xmax, om_lit = t
        k = (part, round(om, 6))
        if k in FAILED_KEYS.get(tag, ()):
            failed.append(t)
    return failed


# ------------------------------------------------------------------------
# _derive_branches_robust -- copied VERBATIM from probe_ff_residual_
# screen_reconstruction_fix.py (which itself copied it verbatim from
# probe_threshold_sensitivity.py's 2026-07-19 "FOURTH PASS" version, job
# 2325369-validated at 352/352 points). Pure branch-reconstruction
# algorithm -- no SOLVER_VERSION-sensitive content. Do NOT hand-edit
# without re-diffing against that source.
# ------------------------------------------------------------------------
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


def _fidelity_gate(fidelity_row):
    """Confirms the verbatim copy above behaves identically to the
    original on an already-succeeding case before any of the 44 failures
    are touched. Returns (ok, detail)."""
    import plate_solver as ps
    from plate_solver.detectors import full_search, track, select_fill, \
        sigma_min_from_K

    tag, key, part, Om, n_dofs, xmax, Om_lit = fidelity_row
    mat = ps.IsotropicMaterial(E=210e9, nu=NU_SWEEP, rho=7800.0)
    two_T_pi = float(key.rsplit(" ", 1)[1])
    geom = ps.make_geometry(R0_2B[tag], two_T_pi)
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
    return ok, (f"row={tag} {key} part={part} Om={Om:.6f}: "
                f"cnt={cnt}/{n_dofs} how={how!r}")


def _worker(args):
    tag, key, part, Om, n_dofs, xmax, Om_lit = args
    t0 = time.time()
    try:
        import plate_solver as ps
        from plate_solver.detectors import (full_search, track, select_fill,
                                             sigma_min_from_K,
                                             weak_enforcement_residual_oop,
                                             weak_enforcement_residual_ip)

        mat = ps.IsotropicMaterial(E=210e9, nu=NU_SWEEP, rho=7800.0)
        two_T_pi = float(key.rsplit(" ", 1)[1])
        geom = ps.make_geometry(R0_2B[tag], two_T_pi)

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
            return dict(ok=False, tag=tag, key=key, part=part, Om=Om,
                        Om_lit=Om_lit, cnt=cnt, how=how,
                        error=f"still only {cnt}/{n_dofs} even after "
                              f"_derive_branches_robust ({how})",
                        dt=time.time() - t0)

        if part == 1:
            res = weak_enforcement_residual_oop(solver, Om, sel, n_dofs=n_dofs)
        else:
            res = weak_enforcement_residual_ip(solver, Om, sel, n_dofs=n_dofs)

        if not res.get("ok", False):
            return dict(ok=False, tag=tag, key=key, part=part, Om=Om,
                        Om_lit=Om_lit, cnt=cnt, how=how,
                        error=f"branches filled ({how}) but residual screen "
                              f"failed: {res.get('reason')}",
                        dt=time.time() - t0)

        s = solver.sigma_min(Om, sel)
        out = dict(ok=True, tag=tag, key=key, part=part, Om=Om,
                   Om_lit=Om_lit, block=res["block"], cnt=cnt, how=how,
                   log10_sigma_min=float(s), dt=time.time() - t0)
        if part == 1:
            out.update(rV=res["rV"], rM=res["rM"], maxW=res["maxW"],
                       verdict=res["verdict"])
        else:
            extra = {k: v for k, v in res.items() if k not in ("ok", "block")}
            out.update(extra)
        return out
    except Exception as exc:
        return dict(ok=False, tag=tag, key=key, part=part, Om=Om,
                    Om_lit=Om_lit, error=f"{exc}", dt=time.time() - t0)


def main():
    print("=" * 78)
    print("  probe_ff_residual_screen_geomsweep_s10v2_reconstruction_fix "
          "-- start " + time.strftime("%Y-%m-%d %H:%M:%S"))
    print("=" * 78, flush=True)

    import plate_solver as ps
    if ps.SOLVER_VERSION != EXPECT_VER:
        print(f"PREFLIGHT FAIL: SOLVER_VERSION={ps.SOLVER_VERSION!r} != "
              f"{EXPECT_VER!r}")
        raise SystemExit(2)
    print(f"preflight OK: SOLVER_VERSION={ps.SOLVER_VERSION}", flush=True)

    all_targets = _load_all_targets()
    failed = _select_failed(all_targets)
    if len(failed) != N_EXPECTED_FAILED:
        print(f"FATAL: matched {len(failed)}/{N_EXPECTED_FAILED} expected "
              f"failing rows against a fresh manifest reload -- the "
              f"manifests may have changed since jobs 2416866-2416869 ran, "
              f"or FAILED_KEYS is stale. Stopping rather than processing a "
              f"wrong set.")
        raise SystemExit(2)
    n_p1 = sum(1 for t in failed if t[2] == 1)
    n_p2 = sum(1 for t in failed if t[2] == 2)
    print(f"Matched {len(failed)}/{N_EXPECTED_FAILED} failing rows against "
          f"a fresh manifest reload ({n_p1} OOP + {n_p2} IP) -- fidelity "
          f"check passed.\n", flush=True)

    failed_key_set = {(t[0], t[2], round(t[3], 6)) for t in failed}
    fidelity_row = next(
        t for t in all_targets
        if (t[0], t[2], round(t[3], 6)) not in failed_key_set)

    print("--- In-run fidelity gate (copy of _derive_branches_robust vs. "
          "one already-succeeding case) ---", flush=True)
    gate_ok, gate_detail = _fidelity_gate(fidelity_row)
    print(f"  {gate_detail}")
    if not gate_ok:
        print("FIDELITY GATE FAILED -- the copied _derive_branches_robust "
              "does not reproduce 'cold' success on an already-passing "
              "case. Stopping before processing the 44 failures; the copy "
              "needs to be re-diffed against probe_ff_residual_screen_"
              "reconstruction_fix.py.")
        raise SystemExit(2)
    print("  GATE OK -- proceeding.\n", flush=True)

    n_workers = int(os.environ.get("N_WORKERS", str(len(failed))))
    print(f"N_WORKERS={n_workers}  NU_SWEEP={NU_SWEEP}", flush=True)

    results = []
    with ProcessPoolExecutor(max_workers=min(n_workers, len(failed))) as ex:
        futs = {ex.submit(_worker, t): t for t in failed}
        done = 0
        for fut in as_completed(futs):
            r = fut.result()
            done += 1
            partlabel = "OOP" if r["part"] == 1 else "IP "
            lit = f"{r['Om_lit']:.3f}" if r.get("Om_lit") is not None else "n/a"
            if r["ok"]:
                if r["part"] == 1:
                    print(f"  [{done}/{N_EXPECTED_FAILED}] {partlabel} "
                          f"{r['tag']} {r['key']} Om={r['Om']:.6f} "
                          f"(lit={lit}) FIXED via how={r['how']!r} "
                          f"cnt={r['cnt']} block={r['block']:>3s} "
                          f"rV={r['rV']:.4g} rM={r['rM']:.4g} "
                          f"maxW={r['maxW']:.4g} verdict={r['verdict']} "
                          f"[{r['dt']:.0f}s]", flush=True)
                else:
                    extra_str = " ".join(
                        f"{k}={v:.4g}" if isinstance(v, float) else f"{k}={v}"
                        for k, v in r.items()
                        if k not in ("ok", "tag", "key", "part", "Om",
                                     "Om_lit", "dt", "log10_sigma_min",
                                     "block", "cnt", "how"))
                    print(f"  [{done}/{N_EXPECTED_FAILED}] {partlabel} "
                          f"{r['tag']} {r['key']} Om={r['Om']:.6f} "
                          f"(lit={lit}) FIXED via how={r['how']!r} "
                          f"cnt={r['cnt']} {extra_str} "
                          f"(NO verdict bar -- IP exploratory) "
                          f"[{r['dt']:.0f}s]", flush=True)
            else:
                print(f"  [{done}/{N_EXPECTED_FAILED}] {partlabel} "
                      f"{r['tag']} {r['key']} Om={r['Om']:.6f} "
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
    print(f"SUMMARY: {len(ok_results)}/{N_EXPECTED_FAILED} now resolved via "
          f"_derive_branches_robust, {len(still_failed)} still failed")
    print(f"  OOP (part 1): {len(oop_ok)}/12 resolved -- "
          f"{len(real_like)} REAL-like, {len(artifact_like)} ARTIFACT-like, "
          f"{len(ambiguous)} AMBIGUOUS")
    print(f"  IP  (part 2): {len(ip_ok)}/32 resolved -- raw ratios only, "
          f"no verdict bar")
    if still_failed:
        print("\nSTILL FAILED (escalation-resistant -- root cause is now "
              "genuinely open, not an implementation gap, for these "
              "specific points):")
        for r in still_failed:
            print(f"  {r['tag']} {r['key']} part={r['part']} "
                  f"Om={r['Om']:.6f}: {r.get('error')}")
    if still_failed:
        by_ratio = {}
        for r in still_failed:
            by_ratio.setdefault(r["tag"], []).append(r)
        print("\nStill-failed count by ratio (compare to the spot-check's "
              "own concentration, LESSONS_LEARNED.md Sec 18.10):")
        for tag in ("r150", "r167", "r200", "r250"):
            rs = by_ratio.get(tag, [])
            n_oop = sum(1 for r in rs if r["part"] == 1)
            n_ip = sum(1 for r in rs if r["part"] == 2)
            print(f"    {tag}: {len(rs)} still failed ({n_oop} OOP, "
                  f"{n_ip} IP)")
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
