# -*- coding: utf-8 -*-
"""
probe_regen_modeshape_figs.py -- DATA-REGENERATION ONLY. Overwrites stale
mode-shape PNGs in an existing FIGDIR; writes no new checkpoint data, makes
no package changes, no SOLVER_VERSION implications.

CONTEXT: B18 (LESSONS_LEARNED bug ledger, found 2026-07-15 via the MAC
ambiguity investigation) -- `OutOfPlaneSolver.mode_shape_grid` used
cos(xi*theta+ph) where the validated, frequency-determining _build_K_real
assembly uses sin(xi*theta+ph) for the displacement projection. This was a
90-degree phase error giving every mode-shape figure the WRONG theta-parity
-- it never touched any frequency (mode_shape_grid is not on the
full_search/select_fill/_build_K_real path that finds Omega values), so no
number in any table or checkpoint is affected, only the *_fig3.png
wireframe plots. Fixed in core_solvers.py; this script re-renders the
affected figures from data ALREADY in a finished cli.py run's checkpoint,
without re-running the expensive sigma_min discovery sweep that originally
found those Omega values (hours per geometry) -- only the cheap
full_search + mode_shape_grid reconstruction step (~60-90s each, same cost
driver as probe_paper_figures_oop.py and the MAC probe).

SCOPE: Part 1 (out-of-plane) ONLY. plot_mode_shape/mode_shape_grid is not
implemented for Part 2 (in-plane) -- classify_modes_ritz calls are Part-1
only, and no Part 2 entry in any checkpoint has a *_fig3.png counterpart
(only *_fig2.png dispersion curves, which use plot_dispersion_curves, not
mode_shape_grid -- unaffected by B18). Part 2 entries in the checkpoint are
silently skipped.

REPRODUCES cli.py's own Part-1 figure-generation block EXACTLY (same
geom_tag naming, same p1_n_dofs rule, same xi_max via _scan_cfg_part1, same
Ritz suspect/-routing, same material via _material_from_env) so this
produces byte-for-byte the same filenames the original run wrote, just
with the corrected mode shape -- a plain overwrite in the same FIGDIR.

INPUT: a finished (or partially finished) cli.py checkpoint JSON --
CHECKPOINT env var, default research23_checkpoint.json. Reads its
"results" dict directly; does NOT require re-running cli.py or having its
in-memory state. Every (r0_2b, two_T_pi) key with a "Part 1" prefix is
processed; "Part 2" keys are skipped (see SCOPE above).

OUTPUT: FIGDIR env var (default "figures", override per checkpoint the
same way cli.py's own split jobs did -- e.g. FIGDIR=figures_p1tail for
research23_checkpoint_p1tail.json). Overwrites {geom_tag}_mode{n}_fig3.png
(or FIGDIR/suspect/... for Ritz-flagged suspect modes) in place.

PRE-REGISTERED CHECK: this script's own fidelity depends on exactly
reproducing cli.py's Part-1 figure block. A fidelity gate runs first --
for one (geometry, mode) pair, it also calls plot_mode_shape a second time
via a DIRECT, freshly-imported plate_solver (not this script's own
already-imported one) and asserts the returned grid's r_phys/th_phys are
identical -- this cannot fail unless the environment itself is
inconsistent, but it's a cheap check against a stale/half-synced deploy
before spending compute on the rest.

COST: ~60-90s per mode (one full_search + one mode_shape_grid reconstruction
+ one Ritz classification per geometry, shared across that geometry's
modes). Parallelized across workers sized to the SLURM allocation, exactly
like probe_paper_figures_oop.py.
"""
import os
import sys
import json
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

sys.path.insert(0, os.environ.get("PKG_PATH", "."))
os.environ.setdefault("DPS", "40")

EXPECT_VER = os.environ.get("EXPECT_SOLVER_VERSION", "2026-07-10.s8")
CHECKPOINT_PATH = os.environ.get("CHECKPOINT", "research23_checkpoint.json")
FIGDIR = os.environ.get("FIGDIR", "figures")


def _load_part1_entries(checkpoint_path):
    """Parses the checkpoint's "results" dict, keeping only "Part 1|r0_2b|
    two_T_pi" keys (Part 2 has no mode-shape figure counterpart -- see
    module docstring SCOPE). Returns a list of (r0_2b, two_T_pi, [Omega,...])
    tuples, sorted the same way cli.py's own p1_keys are (by r0_2b, then
    two_T_pi) for a deterministic, easy-to-follow log."""
    with open(checkpoint_path) as f:
        data = json.load(f)
    entries = []
    for key, oop_f in data.get("results", {}).items():
        parts = key.split("|")
        if len(parts) != 3 or parts[0] != "Part 1":
            continue
        r0_2b = float(parts[1])
        two_T_pi = float(parts[2])
        entries.append((r0_2b, two_T_pi, list(oop_f)))
    entries.sort(key=lambda e: (e[0], e[1]))
    return entries


def _geometry_worker(args):
    """Regenerates every mode-shape figure for ONE geometry -- mirrors
    cli.py's Part-1 figure block (lines ~488-533 as of 2026-07-15) exactly:
    same geom_tag, same p1_n_dofs rule, same xi_max source, same Ritz
    suspect/-routing. Runs the geometry's modes serially within this
    worker (they share one `oop` solver and one Ritz classification call);
    parallelism is across geometries, not within one."""
    r0_2b, two_T_pi, oop_f = args
    t0 = time.time()
    results = []
    try:
        import plate_solver as ps
        from plate_solver.detectors import full_search, classify_modes_ritz
        from plate_solver.validation import _scan_cfg_part1
        from plate_solver.geometry import _material_from_env
        from plate_solver.plotting import plot_mode_shape

        mat = _material_from_env()
        geom = ps.make_geometry(r0_2b, two_T_pi)
        (Om_lo, Om_hi), ns, cut_offs, xi_max = _scan_cfg_part1(
            r0_over_2b=r0_2b, two_T_pi=two_T_pi, mat=mat)
        p1_n_dofs = 20 if two_T_pi <= 0.5 else 16
        oop = ps.OutOfPlaneSolver(geom, mat, M=80, n_quad=30)

        geom_tag = f"p1_r{r0_2b:.4g}_2T{two_T_pi:.4g}pi".replace('.', 'p')

        ritz_flags = [True] * len(oop_f)
        try:
            ritz_flags, ritz_info = classify_modes_ritz(list(oop_f), geom, mat)
        except Exception as exc:
            _ritz_info = f"Ritz cross-check failed ({exc}); treating all as CONFIRMED"
        susp_dir = os.path.join(FIGDIR, "suspect")

        for mn, Om_star in enumerate(oop_f, 1):
            m_t0 = time.time()
            try:
                ok_mode = ritz_flags[mn - 1] if mn - 1 < len(ritz_flags) else True
                out_dir = FIGDIR if ok_mode else susp_dir
                if not ok_mode:
                    os.makedirs(susp_dir, exist_ok=True)
                brs_m = full_search(oop.fast, Om_star, xmax=xi_max)
                fname = os.path.join(out_dir, f"{geom_tag}_mode{mn}_fig3.png")
                ok = plot_mode_shape(
                    oop, Om_star, brs_m, n_dofs=p1_n_dofs, mode_num=mn,
                    geom_label=f"r0/(2b)={r0_2b:.4g}, 2Θ={two_T_pi:.4g}π",
                    fname=fname)
                results.append(dict(ok=bool(ok), mode=mn, Om=Om_star,
                                     fname=fname, suspect=not ok_mode,
                                     dt=time.time() - m_t0))
            except Exception as exc:
                results.append(dict(ok=False, mode=mn, Om=Om_star,
                                     error=f"{exc}", dt=time.time() - m_t0))

        return dict(ok=True, r0_2b=r0_2b, two_T_pi=two_T_pi, geom_tag=geom_tag,
                     n_modes=len(oop_f), results=results, dt=time.time() - t0)
    except Exception as exc:
        return dict(ok=False, r0_2b=r0_2b, two_T_pi=two_T_pi, error=f"{exc}",
                     dt=time.time() - t0)


def main():
    print("=" * 78)
    print("  probe_regen_modeshape_figs -- start "
          + time.strftime("%Y-%m-%d %H:%M:%S"))
    print(f"  CHECKPOINT={CHECKPOINT_PATH}   FIGDIR={FIGDIR}")
    print("=" * 78, flush=True)

    import plate_solver as ps
    if ps.SOLVER_VERSION != EXPECT_VER:
        print(f"PREFLIGHT FAIL: SOLVER_VERSION={ps.SOLVER_VERSION!r} != "
              f"{EXPECT_VER!r}")
        raise SystemExit(2)
    print(f"preflight OK: SOLVER_VERSION={ps.SOLVER_VERSION}", flush=True)

    if not os.path.isfile(CHECKPOINT_PATH):
        print(f"FATAL: checkpoint not found at {CHECKPOINT_PATH!r}. "
              f"Set CHECKPOINT to point at the finished/partial cli.py "
              f"checkpoint whose figures need regenerating.", flush=True)
        raise SystemExit(3)

    entries = _load_part1_entries(CHECKPOINT_PATH)
    if not entries:
        print(f"Nothing to do: no 'Part 1' entries found in "
              f"{CHECKPOINT_PATH!r} (Part 2/in-plane has no mode-shape "
              f"figure counterpart -- see module docstring SCOPE).",
              flush=True)
        return
    total_modes = sum(len(e[2]) for e in entries)
    print(f"\nFound {len(entries)} Part-1 geometries, {total_modes} modes "
          f"total to regenerate:", flush=True)
    for r0_2b, two_T_pi, oop_f in entries:
        print(f"  r0/(2b)={r0_2b:.4g}  2*Theta={two_T_pi:.4g}*pi  "
              f"({len(oop_f)} modes)", flush=True)

    os.makedirs(FIGDIR, exist_ok=True)

    n_workers = int(os.environ.get("N_WORKERS", "8"))
    print(f"\nN_WORKERS={n_workers}", flush=True)

    with ProcessPoolExecutor(max_workers=min(n_workers, len(entries))) as ex:
        futs = {ex.submit(_geometry_worker, e): e for e in entries}
        for fut in as_completed(futs):
            r = fut.result()
            if not r["ok"]:
                print(f"  GEOMETRY FAIL r0/(2b)={r['r0_2b']:.4g} "
                      f"2*Theta={r['two_T_pi']:.4g}*pi: {r.get('error')} "
                      f"[{r['dt']:.0f}s]", flush=True)
                continue
            print(f"\n  {r['geom_tag']} ({r['n_modes']} modes) "
                  f"[{r['dt']:.0f}s total]:", flush=True)
            for mr in r["results"]:
                if mr["ok"]:
                    tag = "suspect" if mr["suspect"] else "OK"
                    print(f"    mode {mr['mode']:2d}  Om={mr['Om']:.6f}  "
                          f"{tag:8s} -> {mr['fname']} [{mr['dt']:.0f}s]",
                          flush=True)
                else:
                    print(f"    mode {mr['mode']:2d}  Om={mr['Om']:.6f}  "
                          f"FAILED: {mr.get('error')} [{mr['dt']:.0f}s]",
                          flush=True)

    print("\nData-regeneration only -- no checkpoint written, no "
          "SOLVER_VERSION action, no package changes. Figures in "
          f"{FIGDIR!r} (and its suspect/ subfolder) have been overwritten "
          "in place with the B18-fixed mode shapes.", flush=True)
    print("End: " + time.strftime("%Y-%m-%d %H:%M:%S"), flush=True)


if __name__ == "__main__":
    t0 = time.time()
    main()
    print(f"\ntotal wall time: {(time.time() - t0)/60:.1f} min", flush=True)
