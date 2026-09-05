# -*- coding: utf-8 -*-
"""
probe_ff_residual_screen_geomsweep_s10v2.py -- DIAGNOSTIC/DATA-GENERATION
ONLY. Writes only a plain-text log; no package changes, no SOLVER_VERSION
implications.

CONTEXT: PAPER1_FREEFREE_DRAFT.tex's Table 13 caption and Appendix C.3
report a cut-off-adjacent residual-screen tally (119 flexural: 53
REAL-like, 60 ARTIFACT-like, 6 AMBIGUOUS; 96 extensional: 72/96 bimodal,
24 middle-band) that was computed under the PRE-s10 checkpoint (2026-07,
SOLVER_VERSION 2026-07-10.s8) by probe_ff_residual_screen_geomsweep.py,
using a hardcoded TARGETS list of 215 candidates extracted by hand from
figures_ff_r*/freefree_manifest.json's own diag.accepted lists.

Since then, the geomsweep was fully recaptured under the current s10
checkpoint (jobs 2414027/2414243/2414039/2411853; LESSONS_LEARNED.md
Sec 18.8) and the s10v2 candidate sets churned substantially relative to
the pre-s10 ones (new high-frequency IP modes appeared at large sector
angles; some low-Omega rigid-body-adjacent artifacts were dropped -- see
Sec 18.8's own r250 comparison). A separate follow-up (jobs 2415646/
2415645/2415648/2414038) already screened the DELTA between the old and
new candidate sets (new-only + cross-geometry-suspect candidates only,
via probe_ff_residual_screen_s10v2_growth.py) -- that is NOT the same
population as the cut-off-adjacent tally above, and does not license
updating it. PAPER1_FREEFREE_DRAFT.tex's Table 13 caption (as of
2026-08-13) explicitly flags the 119/96 tally as "from the pre-s10
candidate set and not yet independently re-screened against the
recaptured manifests" -- this script is exactly that re-screen.

WHAT THIS DOES, per (r0_2b, two_T_pi, part) key present in the s10v2
manifest figures_ff_{RATIO_TAG}_s10v2/freefree_manifest.json:
  1. Reads that key's own diag.cut_offs, diag.accepted, diag.n_dofs and
     diag.max_dim DIRECTLY from the manifest (not hardcoded, not
     re-derived) -- exactly the same extraction rule the original,
     pre-s10 script's TARGETS list was hand-built from: n_cutoffs =
     len(diag.cut_offs); take the first min(n_cutoffs, len(accepted))
     accepted candidates, sorted ascending by Omega. A dry run against
     the actual s10v2 manifests (staged from figures_ff_{tag}_s10v2.zip)
     confirms this reproduces the ORIGINAL per-ratio target counts exactly
     (53/54/54/54 = 215 total for r150/r167/r200/r250), which is expected
     -- n_cutoffs is a property of the geometry (6 angles x {5 OOP, 4 IP}
     cutoffs per angle), not of which candidates happen to sit near it.
     The actual Omega VALUES differ from the pre-s10 run; that is the
     point of re-running this under s10.
  2. Maps each target Omega back to its Omega_lit via the same key's
     'raw'/'omega_lit' parallel lists (relative tolerance 0.002), the
     same matching technique probe_ff_residual_screen_s10v2_growth.py
     already used for its own new-only/suspect population.
  3. Reconstructs the branch set at that Omega via full_search +
     select_fill using THAT KEY'S OWN diag.n_dofs / diag.max_dim (read
     from the manifest, not guessed) -- the same fixed-reconstruction
     discipline as both prior geomsweep-residual probes.
  4. Runs detectors.weak_enforcement_residual_oop (part=1, has a
     REAL-like/ARTIFACT-like/AMBIGUOUS verdict bar) or
     weak_enforcement_residual_ip (part=2, NO verdict bar -- raw
     rTyy/rTyr/maxDisp only, per that function's own docstring).
  5. Cross-checks solver.sigma_min at the same (Omega, sel) as a cheap
     independent depth sanity check (reuses sel, no extra full_search).

HONEST SCOPE CAVEAT -- read before trusting any verdict this probe
prints, carried forward UNCHANGED from both prior geomsweep-residual
probes: weak_enforcement_residual_oop's thresholds (real_rM_max=0.1,
artifact_rM_min=0.3, etc.) were calibrated ONCE, at ONE geometry (FF-P1,
r0/2b=1.5, 2Theta=0.5pi) and ONE material (nu=0.30), then applied
unchanged to r0_2b in {1.5, 1.667, 2.0, 2.5}, all 6 angles, and nu=0.35
(the geomsweep's actual material) -- a genuine extrapolation of the
thresholds' validated domain, same as before. weak_enforcement_residual_ip
has NO calibrated verdict bar at all -- report raw ratios only, and per
the growth-probe's own finding (LESSONS_LEARNED Sec 17: REAL-like at one
angle, ARTIFACT-like at another, same radius ratio), do NOT pool IP rows
into a single count; report the distribution by angle.

PRE-REGISTERED INTERPRETATION (write this BEFORE the job runs, honor it
after):
  - OOP verdict REAL-like with rM in the FF-P1 calibration's ~0.002-0.04
    real band -> supports the candidate being a genuine mode, same
    standard as tab:oop/tab:oop-spurious and the pre-s10 cut-off-adjacent
    table this replaces.
  - OOP verdict ARTIFACT-like with rM in the ~75-350 artifact band ->
    supports it being a weak-enforcement artifact, consistent with the
    rest of this project's spurious-root findings.
  - OOP AMBIGUOUS, or a verdict whose raw rM sits well outside either
    band -> do not trust the label; report the raw number.
  - A branch-reconstruction FAILURE (select_fill short of n_dofs, "no
    null vector", dimension mismatch) -> report plainly; do not drop the
    row or substitute a neighboring Omega.
  - IP (part=2): report raw rTyy/rTyr/maxDisp only, broken out by angle.
  - THE DECISION THIS JOB IS FOR: compare the resulting OOP tally
    (REAL-like/ARTIFACT-like/AMBIGUOUS counts and their sum) against the
    pre-s10 53/60/6 = 119 breakdown currently in the paper.
      * If the total is close (say within ~10%) and the REAL-like
        fraction is comparable (currently 53/119 = 44.5%), that supports
        a straightforward number-swap in Table 13's caption and
        Appendix C.3, keeping the same narrative.
      * If the total or the REAL-like fraction shifts substantially,
        that is itself the finding -- do not force the old narrative
        onto new numbers; report the shift plainly, the same way Sec 18.8
        reported the main FE-match statistic's +1.7pp shift rather than
        assuming "essentially unchanged" without checking.
      * The extensional (IP) side has no verdict bar, so there is no
        equivalent "72/96 bimodal" claim to reproduce automatically --
        that number came from a human-legible rTyy histogram inspection
        in the original write-up (Appendix C.3), and reproducing it
        requires the same by-hand inspection of this run's raw ratios,
        not an automated threshold.

COST: 215 candidates total (53 r150, 54 each r167/r200/r250 -- confirmed
by a dry-run target-count check against the actual s10v2 manifests, see
above), split across 4 jobs (one per radius ratio, matching every prior
stage of this thread's convention). Per-point cost mirrors both prior
geomsweep-residual probes' own estimate (1 full_search + 1 residual-screen
build + 1 cheap sigma_min reuse, dps=40) ~= ~150s/point serial.
Parallelized 1 point/worker via ProcessPoolExecutor sized to the SLURM
allocation -- wall time should be a few minutes per job plus queue
overhead at 64 workers; --time=01:00:00 is a generous margin (mirrors
submit_ff_residual_screen_s10v2_growth_r250.sh, which used the same
margin for a comparable-cost, larger-target-count job).

Env var RATIO_TAG (r150/r167/r200/r250) selects which ratio's manifest
this run reads -- set by the matching submit_ff_residual_screen_
geomsweep_s10v2_r*.sh script, not passed on the command line.
"""
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

sys.path.insert(0, os.environ.get("PKG_PATH", "."))
os.environ.setdefault("DPS", "40")

EXPECT_VER = os.environ.get("EXPECT_SOLVER_VERSION", "2026-07-10.s10")
RATIO_TAG = os.environ.get("RATIO_TAG")
if not RATIO_TAG:
    print("FATAL: RATIO_TAG env var not set (expected one of r150/r167/"
          "r200/r250). Set it in the submit script.")
    raise SystemExit(2)

R0_2B = {"r150": 1.5, "r167": 1.66667, "r200": 2.0, "r250": 2.5}[RATIO_TAG]
NU_SWEEP = 0.35  # matches geometry_sweep's actual material -- see both
                 # prior geomsweep-residual probes' own note; NOT the
                 # FF-P1 tab:oop/tab:ip benchmark's nu=0.30.

MANIFEST = f"figures_ff_{RATIO_TAG}_s10v2/freefree_manifest.json"
TOL = 0.002  # relative Omega match tolerance for the raw->omega_lit lookup


def _load_targets():
    """Data-driven target discovery from the s10v2 manifest itself: for
    each (r0_2b, two_T_pi, part) key, the first min(n_cutoffs,
    n_accepted) accepted candidates sorted ascending by Omega -- the same
    extraction rule the pre-s10 script's hardcoded TARGETS list encoded
    by hand. n_cutoffs = len(diag['cut_offs'])."""
    if not os.path.isfile(MANIFEST):
        print(f"FATAL: manifest not found: {MANIFEST} -- unzip "
              f"figures_ff_{RATIO_TAG}_s10v2.zip into that directory "
              f"(relative to cwd) before submitting this job.")
        raise SystemExit(2)

    manifest = json.load(open(MANIFEST))
    targets = []  # (key, part, Om, n_dofs, xmax, Om_lit)

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
            targets.append((key, part, om, n_dofs, xmax, om_lit))

    return targets


def _worker(args):
    key, part, Om, n_dofs, xmax, Om_lit = args
    t0 = time.time()
    try:
        import plate_solver as ps
        from plate_solver.detectors import (full_search, select_fill,
                                             weak_enforcement_residual_oop,
                                             weak_enforcement_residual_ip)

        mat = ps.IsotropicMaterial(E=210e9, nu=NU_SWEEP, rho=7800.0)
        two_T_pi = float(key.rsplit(" ", 1)[1])
        geom = ps.make_geometry(R0_2B, two_T_pi)

        if part == 1:
            solver = ps.OutOfPlaneSolver(geom, mat, M=80, n_quad=30,
                                          boundary=ps.FreeFreeOOP())
        else:
            solver = ps.InPlaneSolver(geom, mat, M=80, n_quad=30,
                                       boundary=ps.FreeFreeIP())

        raw = full_search(solver.fast, Om, xmax=xmax)
        sel, cnt = select_fill(raw, n_dofs)
        if cnt < n_dofs:
            return dict(ok=False, key=key, part=part, Om=Om, Om_lit=Om_lit,
                        error=f"only {cnt}/{n_dofs} branches filled",
                        dt=time.time() - t0)

        if part == 1:
            res = weak_enforcement_residual_oop(solver, Om, sel, n_dofs=n_dofs)
        else:
            res = weak_enforcement_residual_ip(solver, Om, sel, n_dofs=n_dofs)

        if not res.get("ok", False):
            return dict(ok=False, key=key, part=part, Om=Om, Om_lit=Om_lit,
                        error=f"residual screen failed: {res.get('reason')}",
                        dt=time.time() - t0)

        s = solver.sigma_min(Om, sel)

        out = dict(ok=True, key=key, part=part, Om=Om, Om_lit=Om_lit,
                   block=res["block"], log10_sigma_min=float(s),
                   dt=time.time() - t0)
        if part == 1:
            out.update(rV=res["rV"], rM=res["rM"], maxW=res["maxW"],
                       verdict=res["verdict"])
        else:
            extra = {k: v for k, v in res.items() if k not in ("ok", "block")}
            out.update(extra)
        return out
    except Exception as exc:
        return dict(ok=False, key=key, part=part, Om=Om, Om_lit=Om_lit,
                    error=f"{exc}", dt=time.time() - t0)


def main():
    print("=" * 78)
    print(f"  probe_ff_residual_screen_geomsweep_s10v2 -- RATIO_TAG={RATIO_TAG}"
          f" -- start " + time.strftime("%Y-%m-%d %H:%M:%S"))
    print("=" * 78, flush=True)

    import plate_solver as ps
    if ps.SOLVER_VERSION != EXPECT_VER:
        print(f"PREFLIGHT FAIL: SOLVER_VERSION={ps.SOLVER_VERSION!r} != "
              f"{EXPECT_VER!r}")
        raise SystemExit(2)
    print(f"preflight OK: SOLVER_VERSION={ps.SOLVER_VERSION}", flush=True)

    targets = _load_targets()
    if not targets:
        print(f"FATAL: 0 targets discovered for RATIO_TAG={RATIO_TAG!r} -- "
              f"inspect {MANIFEST} by hand; every key should contribute at "
              f"least one cut-off-adjacent candidate.")
        raise SystemExit(2)

    n_workers = int(os.environ.get("N_WORKERS", str(len(targets))))
    n_p1 = sum(1 for t in targets if t[1] == 1)
    n_p2 = sum(1 for t in targets if t[1] == 2)
    print(f"N_WORKERS={n_workers}  NU_SWEEP={NU_SWEEP}  r0_2b={R0_2B}")
    print(f"{len(targets)} targets for {RATIO_TAG}: {n_p1} OOP (part 1) + "
          f"{n_p2} IP (part 2)\n", flush=True)

    results = []
    with ProcessPoolExecutor(max_workers=min(n_workers, len(targets))) as ex:
        futs = {ex.submit(_worker, t): t for t in targets}
        done = 0
        for fut in as_completed(futs):
            r = fut.result()
            done += 1
            partlabel = "OOP" if r["part"] == 1 else "IP "
            lit = f"{r['Om_lit']:.3f}" if r.get("Om_lit") is not None else "n/a"
            if r["ok"]:
                if r["part"] == 1:
                    print(f"  [{done}/{len(targets)}] {partlabel} {r['key']} "
                          f"Om={r['Om']:.6f} (lit={lit}) block={r['block']:>3s} "
                          f"rV={r['rV']:.4g} rM={r['rM']:.4g} maxW={r['maxW']:.4g} "
                          f"log10smin={r['log10_sigma_min']:.4f} "
                          f"verdict={r['verdict']} [{r['dt']:.0f}s]", flush=True)
                else:
                    extra_str = " ".join(
                        f"{k}={v:.4g}" if isinstance(v, float) else f"{k}={v}"
                        for k, v in r.items()
                        if k not in ("ok", "key", "part", "Om", "Om_lit",
                                     "dt", "log10_sigma_min", "block"))
                    print(f"  [{done}/{len(targets)}] {partlabel} {r['key']} "
                          f"Om={r['Om']:.6f} (lit={lit}) block={r['block']:>3s} "
                          f"{extra_str} log10smin={r['log10_sigma_min']:.4f} "
                          f"(NO verdict bar -- IP exploratory) "
                          f"[{r['dt']:.0f}s]", flush=True)
            else:
                print(f"  [{done}/{len(targets)}] {partlabel} {r['key']} "
                      f"Om={r['Om']:.6f} FAILED: {r.get('error')} "
                      f"[{r['dt']:.0f}s]", flush=True)
            results.append(r)

    ok_results = [r for r in results if r["ok"]]
    failures = [r for r in results if not r["ok"]]
    oop_ok = [r for r in ok_results if r["part"] == 1]
    ip_ok = [r for r in ok_results if r["part"] == 2]
    real_like = [r for r in oop_ok if r["verdict"] == "REAL-like"]
    artifact_like = [r for r in oop_ok if r["verdict"] == "ARTIFACT-like"]
    ambiguous = [r for r in oop_ok if r["verdict"] == "AMBIGUOUS"]

    print("\n" + "=" * 78)
    print(f"SUMMARY ({RATIO_TAG}): {len(ok_results)}/{len(targets)} ok, "
          f"{len(failures)} failures")
    print(f"  OOP (part 1): {len(oop_ok)} screened -- "
          f"{len(real_like)} REAL-like, {len(artifact_like)} ARTIFACT-like, "
          f"{len(ambiguous)} AMBIGUOUS")
    print(f"  IP  (part 2): {len(ip_ok)} screened -- raw ratios only, no "
          f"verdict bar. Breakdown by angle (do not pool -- see docstring):")
    angles = sorted({r["key"].rsplit(" ", 1)[1] for r in ip_ok}, key=float)
    for ang in angles:
        rows = [r for r in ip_ok if r["key"].endswith(" " + ang)]
        below = sum(1 for r in rows if r.get("rTyy", 1e9) < 0.3)
        above = sum(1 for r in rows if r.get("rTyy", -1) > 4)
        mid = len(rows) - below - above
        print(f"    2T/pi={ang}: n={len(rows)} (rTyy<0.3: {below}, "
              f"0.3-4: {mid}, >4: {above})")
    if ambiguous:
        print("  AMBIGUOUS rows (do not trust either label without a closer look):")
        for r in ambiguous:
            print(f"    {r['key']} Om={r['Om']:.6f} rM={r['rM']:.4g}")
    if failures:
        print("FAILURES:")
        for r in failures:
            print(f"  {r['key']} Om={r['Om']:.6f}: {r.get('error')}")
    print("=" * 78, flush=True)

    print("\nDiagnostic/data-generation only -- no SOLVER_VERSION action, "
          "no package changes. Read the module docstring's HONEST SCOPE "
          "CAVEAT and PRE-REGISTERED INTERPRETATION before using any "
          "verdict in a paper claim. Compare the OOP tally above against "
          "the pre-s10 53/60/6=119 breakdown currently in "
          "PAPER1_FREEFREE_DRAFT.tex Table 13's caption / Appendix C.3.",
          flush=True)
    print("End: " + time.strftime("%Y-%m-%d %H:%M:%S"), flush=True)


if __name__ == "__main__":
    t0 = time.time()
    main()
    print(f"\ntotal wall time: {(time.time() - t0)/60:.1f} min", flush=True)
