# -*- coding: utf-8 -*-
"""
probe_ff_adjudicate_candidates.py -- DIAGNOSTIC/DATA-GENERATION ONLY.
Writes only a plain-text log + a JSON results file; no package changes,
no SOLVER_VERSION implications.

WHY THIS EXISTS: PAPER1_FREEFREE_DRAFT.tex's own reviews (2026-08-01,
GrokReviewandAnswers.txt) flagged that the free-free tables are fully
one-to-one residual-adjudicated (REAL/ARTIFACT-verdicted per candidate,
matching FF-P1's own tab:oop/tab:oop-spurious/tab:ip rigor) at exactly
ONE geometry. The 2026-07-21/22 geometry-sweep campaign (FutureWork/
geometry_sweep/, FutureWork/ansys_geomsweep_verification/) already
FE-matched 48 new (radius ratio x sector angle x part) combinations, but
only as a COARSE nearest-neighbor check, with a residual screen applied
only to each geometry's first few cut-off-adjacent candidates (see
probe_ff_residual_screen_geomsweep.py). This script closes that gap
properly for ONE chosen geometry at a time. It is GEOMETRY-AGNOSTIC (all
geometry/material parameters come from CANDIDATES_JSON) -- restored to the
package root 2026-08-08 for reuse on r200_a100_nu030 (a second radius-ratio
adjudication, R_i/R_o=0.600 at the SAME sector angle as r150_a100_nu030),
after its first use closed r150_a100_nu030 (r0/2b=1.5, 2Theta/pi=1.0, "r150/
a100", a genuinely NEW sector angle vs. FF-P1's own PHI=90deg, now archived
under archive/by_topic/r150_a100_second_geometry/).

WHAT THIS DOES: reads a candidates JSON (produced either by extracting
the already-computed accepted-Omega lists from the geometry sweep's own
freefree_manifest.json, or by a fresh probe_ff_<tag>_capture_nu030.py
production scan), and for EVERY accepted candidate (not just cut-off-
adjacent ones, unlike probe_ff_residual_screen_geomsweep.py):
  1. Reconstructs the branch set at that Omega via full_search +
     select_fill, using the SAME n_dofs/xmax the candidates file records.
  2. Runs detectors.weak_enforcement_residual_oop (part=1) or
     weak_enforcement_residual_ip (part=2) for the edge-residual verdict.
     OOP: uses the function's OWN verdict field (REAL-like/ARTIFACT-like/
     AMBIGUOUS), calibrated once at FF-P1 (nu=0.30) -- extrapolated here
     to a new angle/ratio and (for any nu=0.35 candidates file) a new
     material, exactly the honest extrapolation
     probe_ff_residual_screen_geomsweep.py already flags; report raw
     (rV, rM, maxW) alongside every verdict.
     IP: the two-sided verdict bar is mode-shape MAC on the root-owning
     parity block (ip_mac_bar.py, Rank 18, any threshold in
     (0.5801, 0.9079)). This probe does not compute MAC -- that needs an
     FE eigenvector dump and RULE_LOCALMIN_RESID reconstruction, which
     is the production IP classification path, not this residual rematch.
     What this probe still prints for IP is raw rTyy/rTyr/maxDisp plus
     the exploratory rTyy<0.3/>4.0 heuristic classify_ip(), clearly
     labeled NOT the two-sided call. Do not present that heuristic as
     a validated IP verdict in any paper table.
  3. Cross-checks solver.sigma_min at the same (Omega, sel) as a cheap
     independent depth sanity check.
  4. Matches each candidate's Omega_lit (OOP) or raw Hz (IP, via k_bar/
     w_bar) against the geometry's own independent ANSYS FE mode list
     (nearest-neighbor within a configurable tolerance), reporting the
     match error percentage alongside the residual verdict -- this is
     the actual one-to-one adjudication a coarse nearest-neighbor check
     does not provide.

PRE-REGISTERED INTERPRETATION:
  - OOP candidate: REAL-like verdict + FE match <1% -> genuinely supports
    a new tab:oop-style table row at this geometry.
  - OOP candidate: ARTIFACT-like verdict + no FE match -> supports a new
    tab:oop-spurious-style row (a weak-enforcement artifact, not a
    physical mode, exactly analogous to FF-P1's 5 spurious zeros).
  - OOP candidate: AMBIGUOUS, or a verdict whose raw rM sits well outside
    either FF-P1 calibration band (~0.002-0.04 real / ~75-350 artifact) ->
    do NOT trust the label; report raw numbers only, flag for a human
    read, exactly as probe_ff_residual_screen_geomsweep.py's own
    docstring already specifies for this situation.
  - IP candidate: report raw ratios + heuristic label + FE match
    percentage; do NOT present the heuristic label as the two-sided IP
    call. That call is ip_mac_bar.classify_mac on the root-owning block
    (LESSONS_LEARNED.md §18.22--§18.24). The residual heuristic stays
    exploratory and one-sided.
  - A reconstruction/screen FAILURE (select_fill short of n_dofs, no null
    vector) -> report plainly, do not substitute or drop silently.
  - A near-degenerate pair (two candidates both matching one FE mode, gap
    of only a few tenths of a percent in Omega_lit or Hz) -> flag it
    explicitly; this is exactly the FF-P1 ODD-pair pattern (Sec 6.3) and
    should NOT be force-assigned here either without the same multi-
    instrument treatment FF-P1's pair received.

COST: mirrors probe_ff_residual_screen_geomsweep.py's own estimate --
1 full_search (dps=40) + 1 residual-screen build + 1 cheap sigma_min
reuse, ~150s/point serial, parallelized 1 point/worker via
ProcessPoolExecutor sized to the SLURM allocation. Set N_WORKERS to (at
least) the candidate count in CANDIDATES_JSON for a few-minutes wall time.

Env vars:
  CANDIDATES_JSON  -- path to the candidates file (required)
  ANSYS_DIR        -- directory containing the FE .txt files named in the
                      candidates JSON (default: FutureWork/
                      ansys_geomsweep_verification relative to PKG_PATH's
                      parent, i.e. the package root's sibling FutureWork/)
  FE_MATCH_TOL_PCT -- nearest-neighbor FE match tolerance to report as
                      "matched" (default 3.0, matching the geomsweep
                      campaign's own coarse-match tolerance; the residual
                      verdict, not this tolerance, is what actually
                      adjudicates REAL vs ARTIFACT)
"""
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

sys.path.insert(0, os.environ.get("PKG_PATH", "."))
os.environ.setdefault("DPS", "40")
# This probe is free-free only. _worker_bc() defaults to clamped_free;
# export the token before any pool starts (Rank 18 standing action).
os.environ["R40_BC_KIND"] = "free_free"

EXPECT_VER = os.environ.get("EXPECT_SOLVER_VERSION", "2026-07-10.s10")
CANDIDATES_JSON = os.environ.get("CANDIDATES_JSON")
if not CANDIDATES_JSON:
    print("FATAL: CANDIDATES_JSON env var not set.")
    raise SystemExit(2)

with open(CANDIDATES_JSON) as fh:
    CFG = json.load(fh)

ANSYS_DIR = os.environ.get(
    "ANSYS_DIR",
    os.path.join(os.environ.get("PKG_PATH", "."), "..",
                 "FutureWork", "ansys_geomsweep_verification"))
FE_MATCH_TOL_PCT = float(os.environ.get("FE_MATCH_TOL_PCT", "3.0"))
N_WORKERS = int(os.environ.get("N_WORKERS",
                                str(os.cpu_count() or 4)))


def classify_ip(rTyy):
    """Exploratory residual heuristic ONLY. Not the Rank 18 two-sided bar.

    The two-sided IP call is ip_mac_bar.classify_mac (MAC on the
    root-owning block). This rTyy cut is the old one-sided screen and
    remains labeled unvalidated.
    """
    if rTyy < 0.3:
        return "REAL-like (unvalidated)"
    if rTyy > 4.0:
        return "ARTIFACT-like (unvalidated)"
    return "AMBIGUOUS (unvalidated)"


def _read_fe_oop(path, rigid_count):
    rows = []
    with open(path) as fh:
        lines = fh.readlines()[1:]  # skip header
    for ln in lines:
        parts = ln.split()
        if len(parts) < 3:
            continue
        rows.append(dict(mode=int(float(parts[0].rstrip("."))),
                          f_hz=float(parts[1]), omega_lit=float(parts[2])))
    return rows[rigid_count:]


def _read_fe_ip(path, rigid_count):
    rows = []
    with open(path) as fh:
        lines = fh.readlines()[1:]
    for ln in lines:
        parts = ln.split()
        if len(parts) < 2:
            continue
        rows.append(dict(mode=int(float(parts[0].rstrip("."))),
                          f_hz=float(parts[1])))
    return rows[rigid_count:]


XMAX_STEP = 8.0
XMAX_CAP = 120.0


def _worker(args):
    (part, Om, n_dofs, xmax, r0_2b, two_T_pi, nu, E, rho) = args
    t0 = time.time()
    try:
        os.environ["R40_BC_KIND"] = "free_free"
        import plate_solver as ps
        from plate_solver.detectors import (full_search, select_fill,
                                             weak_enforcement_residual_oop,
                                             weak_enforcement_residual_ip)
        from plate_solver.geometry import _omega_lit, _kbar_wbar

        mat = ps.IsotropicMaterial(E=E, nu=nu, rho=rho)
        geom = ps.make_geometry(r0_2b, two_T_pi)
        k_bar, w_bar = _kbar_wbar(geom, mat)

        if part == 1:
            solver = ps.OutOfPlaneSolver(geom, mat, M=80, n_quad=30,
                                          boundary=ps.FreeFreeOOP())
        else:
            solver = ps.InPlaneSolver(geom, mat, M=80, n_quad=30,
                                       boundary=ps.FreeFreeIP())

        # ADAPTIVE XMAX GROWTH -- the xmax recorded in the candidates JSON is
        # whatever find_modes_sigmin's scan-phase config used (often the
        # geometry's auto-derived _scan_cfg_part1/2 starting value, e.g. 14.0),
        # which was only ever proven sufficient for the n_dofs THAT capture
        # actually used. Reusing it verbatim against a DIFFERENT (often larger,
        # e.g. n_dofs=20 production-default) basis silently starves the branch
        # search -- exactly the bug this project already hit once for zone-B
        # (see LESSONS_LEARNED.md, probe_ip_zoneb_select_fill_diag_v1). Grow
        # xmax here the same way, rather than trusting the recorded value.
        cur_xmax = xmax
        raw = full_search(solver.fast, Om, xmax=cur_xmax)
        sel, cnt = select_fill(raw, n_dofs)
        while cnt < n_dofs and cur_xmax < XMAX_CAP:
            cur_xmax += XMAX_STEP
            raw = full_search(solver.fast, Om, xmax=cur_xmax)
            sel, cnt = select_fill(raw, n_dofs)
        if cnt < n_dofs:
            return dict(ok=False, part=part, Om=Om,
                        error=f"only {cnt}/{n_dofs} branches filled "
                              f"(xmax grown to cap {XMAX_CAP})",
                        dt=time.time() - t0)

        f_hz, om_lit = _omega_lit(Om, geom, mat, w_bar, k_bar, part)

        if part == 1:
            res = weak_enforcement_residual_oop(solver, Om, sel, n_dofs=n_dofs)
        else:
            res = weak_enforcement_residual_ip(solver, Om, sel, n_dofs=n_dofs)

        if not res.get("ok", False):
            return dict(ok=False, part=part, Om=Om, f_hz=float(f_hz),
                        om_lit=float(om_lit),
                        error=f"residual screen failed: {res.get('reason')}",
                        dt=time.time() - t0)

        s = solver.sigma_min(Om, sel)
        out = dict(ok=True, part=part, Om=float(Om), f_hz=float(f_hz),
                   om_lit=float(om_lit), log10_sigma_min=float(s),
                   dt=time.time() - t0)
        if part == 1:
            out.update(rV=res["rV"], rM=res["rM"], maxW=res["maxW"],
                       verdict=res["verdict"])
        else:
            extra = {k: v for k, v in res.items() if k not in ("ok", "block")}
            out.update(extra)
            out["verdict"] = classify_ip(res.get("rTyy", float("inf")))
        return out
    except Exception as exc:
        return dict(ok=False, part=part, Om=Om, error=f"{exc}",
                    dt=time.time() - t0)


def match_fe(value, fe_rows, key, tol_pct):
    best = None
    for row in fe_rows:
        if row[key] <= 0:
            continue
        err = abs(value - row[key]) / row[key] * 100.0
        if best is None or err < best[1]:
            best = (row, err)
    if best is None:
        return None
    row, err = best
    return dict(fe_mode=row["mode"], fe_value=row[key], err_pct=err,
                matched=err <= tol_pct)


def main():
    print("=" * 78)
    print(f"  probe_ff_adjudicate_candidates -- CANDIDATES_JSON={CANDIDATES_JSON}")
    print(f"  tag={CFG.get('tag')}  r0_2b={CFG['r0_2b']}  "
          f"two_T_pi={CFG['two_T_pi']}  nu={CFG['nu']}")
    print(f"  R40_BC_KIND={os.environ.get('R40_BC_KIND')!r}  "
          f"(FFFF probe; two-sided IP bar is ip_mac_bar.py, not classify_ip)")
    print("=" * 78, flush=True)

    import plate_solver as ps
    if ps.SOLVER_VERSION != EXPECT_VER:
        print(f"PREFLIGHT WARN: SOLVER_VERSION={ps.SOLVER_VERSION!r} != "
              f"{EXPECT_VER!r} (continuing; update EXPECT_SOLVER_VERSION "
              f"if intentional)")

    jobs = []
    for Om in CFG["oop_omegas"]:
        jobs.append((1, Om, CFG["oop_n_dofs"], CFG["oop_xmax"],
                     CFG["r0_2b"], CFG["two_T_pi"], CFG["nu"],
                     CFG["E"], CFG["rho"]))
    for Om in CFG["ip_omegas"]:
        jobs.append((2, Om, CFG["ip_n_dofs"], CFG["ip_xmax"],
                     CFG["r0_2b"], CFG["two_T_pi"], CFG["nu"],
                     CFG["E"], CFG["rho"]))
    print(f"  {len(jobs)} candidates total "
          f"({len(CFG['oop_omegas'])} OOP + {len(CFG['ip_omegas'])} IP), "
          f"N_WORKERS={N_WORKERS}", flush=True)

    results = []
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=N_WORKERS) as ex:
        futs = {ex.submit(_worker, j): j for j in jobs}
        for fut in as_completed(futs):
            r = fut.result()
            results.append(r)
            tag = "OOP" if r["part"] == 1 else "IP "
            if r.get("ok"):
                print(f"  [{tag}] Om={r['Om']:.6f} f_hz={r.get('f_hz', 0):.3f} "
                      f"verdict={r.get('verdict')} [{r['dt']:.0f}s]",
                      flush=True)
            else:
                print(f"  [{tag}] Om={r['Om']:.6f} FAILED: {r.get('error')} "
                      f"[{r['dt']:.0f}s]", flush=True)
    print(f"\n  all candidates done in {(time.time()-t0)/60:.1f} min",
          flush=True)

    # FE cross-match
    fe_oop_path = os.path.join(ANSYS_DIR, CFG["fe_oop_file"])
    fe_ip_path = os.path.join(ANSYS_DIR, CFG["fe_ip_file"])
    fe_oop = _read_fe_oop(fe_oop_path, CFG["fe_oop_rigid_count"])
    fe_ip = _read_fe_ip(fe_ip_path, CFG["fe_ip_rigid_count"])
    print(f"\n  FE cross-reference: {len(fe_oop)} physical OOP modes "
          f"({fe_oop_path}), {len(fe_ip)} physical IP modes "
          f"({fe_ip_path})", flush=True)

    for r in results:
        if not r.get("ok"):
            continue
        if r["part"] == 1:
            r["fe_match"] = match_fe(r["om_lit"], fe_oop, "omega_lit",
                                      FE_MATCH_TOL_PCT)
        else:
            r["fe_match"] = match_fe(r["f_hz"], fe_ip, "f_hz",
                                      FE_MATCH_TOL_PCT)

    print("\n" + "=" * 78)
    print("  ADJUDICATED TABLE")
    print("=" * 78)
    for r in sorted(results, key=lambda x: (x["part"], x.get("Om", 0))):
        tag = "OOP" if r["part"] == 1 else "IP "
        if not r.get("ok"):
            print(f"  [{tag}] Om={r['Om']:.6f}  FAILED: {r.get('error')}")
            continue
        fm = r.get("fe_match")
        fm_s = (f"FE#{fm['fe_mode']} err={fm['err_pct']:.3f}% "
                f"{'MATCH' if fm['matched'] else 'no-match'}"
                if fm else "no FE row")
        print(f"  [{tag}] Om={r['Om']:.6f}  verdict={r['verdict']:<24s} "
              f"{fm_s}")

    # near-degenerate flag: two candidates within 0.5% of the SAME FE mode
    print("\n  Near-degenerate / duplicate-FE-match check:")
    seen = {}
    flagged = False
    for r in results:
        if not r.get("ok") or not r.get("fe_match") or not r["fe_match"]["matched"]:
            continue
        key = (r["part"], r["fe_match"]["fe_mode"])
        seen.setdefault(key, []).append(r["Om"])
    for key, oms in seen.items():
        if len(oms) > 1:
            flagged = True
            tag = "OOP" if key[0] == 1 else "IP"
            print(f"    !! {tag} FE mode #{key[1]} matched by {len(oms)} "
                  f"solver candidates: {oms} -- flag for multi-instrument "
                  f"treatment like FF-P1's ODD pair (Sec 6.3), do not "
                  f"force-assign here.")
    if not flagged:
        print("    none found -- every matched FE mode has exactly one "
              "solver candidate.")

    out_name = f"ff_adjudicate_{CFG.get('tag', 'unknown')}_" \
               f"{os.environ.get('SLURM_JOB_ID', 'local')}.json"
    with open(out_name, "w") as fh:
        json.dump(dict(config=CFG, results=results), fh, indent=1)
    print(f"\n  wrote {out_name}")


if __name__ == "__main__":
    t0 = time.time()
    main()
    print(f"\ntotal wall time: {time.time() - t0:.1f}s", flush=True)
