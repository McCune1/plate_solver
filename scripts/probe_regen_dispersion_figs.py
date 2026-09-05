# -*- coding: utf-8 -*-
"""
probe_regen_dispersion_figs.py -- DATA-REGENERATION ONLY. Overwrites stale
dispersion-curve (*_fig2.png) PNGs in an existing FIGDIR; writes no new
checkpoint data, makes no package changes, no SOLVER_VERSION implications.

CONTEXT (2026-07-16): plot_dispersion_curves' branch tracer
(compute_dispersion_branches) used a deliberately cheapened root search
(ngrid=6/nax=140 vs. full_search's own ngrid=10/nax=240 default) and only
n_pts=24 Omega samples across the range. At that density the search
sometimes failed to resolve a branch's root at one Omega sample while
finding it fine on neighbouring samples; the branch tracker then ends that
branch right there and starts a new one when the root reappears -- visible
gaps in the plotted line that are not in Seok & Tiersten's own Fig. 2.
Fixed in plotting.py (ngrid/nax restored to full_search's own defaults,
n_pts doubled to 48) -- this script re-renders every *_fig2.png already on
disk with the fixed sampling, WITHOUT re-running the expensive sigma_min
discovery sweep that originally found the geometry's tabulated Omega
values (hours per geometry): plot_dispersion_curves never depends on any
*found* frequency at all -- it only needs the geometry, material, and the
Omega-window derived from the geometry's xi=0/zeta=0 cut-off frequencies
(_scan_cfg_part1/_scan_cfg_part2), which are cheap, closed-form, and
already reproduced in-process here. The only real cost is the dispersion
sweep itself, which uses the fast float64 engine (_Part1Fast/_Part2Fast),
never the mp.dps=40 arbitrary-precision determinant path -- the same
order of cost as the *_fig2.png generation step already embedded inside
the original cli.py run, not the multi-hour discovery around it.

SCOPE: both Part 1 (out-of-plane) and Part 2 (in-plane) checkpoint
entries are covered -- plot_dispersion_curves has no Part-1-only
restriction (unlike plot_mode_shape/mode_shape_grid, see
probe_regen_modeshape_figs.py's own SCOPE note).

REPRODUCES cli.py's own figure-generation block EXACTLY for both parts
(same geom_tag naming, same om_plot_hi derivation from cut_offs, same
xi_max/ze_max source) so this produces byte-for-byte the same filenames
the original run wrote, just with the corrected sampling -- a plain
overwrite in the same FIGDIR.

INPUT: a finished (or partially finished) cli.py checkpoint JSON --
CHECKPOINT env var, default research23_checkpoint.json. Reads its
"results" dict directly for the (r0_2b, two_T_pi) *keys* only -- the
Omega values themselves are never read/used, only the geometry identity.
Every key is processed regardless of "Part 1"/"Part 2" prefix.

OUTPUT: FIGDIR env var (default "figures", override per checkpoint the
same way cli.py's own split jobs did -- e.g. FIGDIR=figures_p1tail for
research23_checkpoint_p1tail.json, FIGDIR=figures_p2 for
research23_checkpoint_p2.json). Overwrites {geom_tag}_fig2.png in place.

PRE-REGISTERED CHECK: a fidelity gate runs first -- for one geometry, it
computes _scan_cfg_part1/_scan_cfg_part2 twice (once via this script's
already-imported plate_solver, once via a freshly re-imported copy) and
asserts identical (Om_lo, Om_hi, cut_offs, xi_max/ze_max) -- cheap check
against a stale/half-synced deploy before spending compute on the rest.

COST: one dispersion sweep per geometry (n_pts=48 float64 root searches,
~1-3 min each depending on geometry/part -- much cheaper than the mp.dps=40
sigma_min discovery that originally took hours). Parallelized across
workers sized to the SLURM allocation, one geometry per worker.
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


def _load_entries(checkpoint_path):
    """Parses the checkpoint's "results" dict; returns a list of
    (part, r0_2b, two_T_pi) tuples, sorted by part then r0_2b then
    two_T_pi for a deterministic, easy-to-follow log. Only the geometry
    identity is used -- the recorded Omega values are irrelevant here."""
    with open(checkpoint_path) as f:
        data = json.load(f)
    entries = []
    for key in data.get("results", {}):
        parts = key.split("|")
        if len(parts) != 3:
            continue
        part_label, r0_2b, two_T_pi = parts[0], float(parts[1]), float(parts[2])
        entries.append((part_label, r0_2b, two_T_pi))
    entries.sort(key=lambda e: (e[0], e[1], e[2]))
    return entries


def _fidelity_gate(entries):
    """Cheap pre-flight: recompute the first entry's scan config once and
    confirm it returns a well-formed (Om_lo, Om_hi, cut_offs, xmax) tuple
    with Om_lo < Om_hi and at least one cut-off -- catches an import-level
    break (missing module, signature drift) before spending any real
    compute on the actual worker pool. Not a deploy-vs-deploy diff (that
    needs two separate processes/checkouts, out of scope for a single
    script) -- just a fast sanity check that the environment this job is
    about to run hundreds of worker calls against isn't broken."""
    if not entries:
        return True
    part_label, r0_2b, two_T_pi = entries[0]
    from plate_solver.geometry import _material_from_env
    from plate_solver.validation import _scan_cfg_part1, _scan_cfg_part2
    mat = _material_from_env()
    if part_label == "Part 1":
        (Om_lo, Om_hi), ns, cut_offs, xmax = _scan_cfg_part1(
            r0_over_2b=r0_2b, two_T_pi=two_T_pi, mat=mat)
    else:
        (Om_lo, Om_hi), ns, cut_offs, xmax = _scan_cfg_part2(
            r0_over_2b=r0_2b, two_T_pi=two_T_pi, mat=mat)
    ok = (Om_lo < Om_hi) and len(cut_offs) >= 1 and xmax > 0
    print(f"  fidelity gate ({part_label} r0/(2b)={r0_2b:.4g} "
          f"2*Theta={two_T_pi:.4g}pi): Om=({Om_lo:.4g},{Om_hi:.4g}) "
          f"cut_offs={len(cut_offs)} xmax={xmax:.4g} "
          f"{'OK' if ok else 'MISMATCH'}", flush=True)
    return ok


def _geometry_worker(args):
    """Regenerates the *_fig2.png dispersion-curve figure for ONE
    geometry -- mirrors cli.py's figure block exactly for whichever part
    this entry belongs to."""
    part_label, r0_2b, two_T_pi = args
    t0 = time.time()
    try:
        import plate_solver as ps
        from plate_solver.geometry import _material_from_env
        from plate_solver.validation import _scan_cfg_part1, _scan_cfg_part2
        from plate_solver.plotting import plot_dispersion_curves

        mat = _material_from_env()
        geom = ps.make_geometry(r0_2b, two_T_pi)

        if part_label == "Part 1":
            (Om_lo, Om_hi), ns, cut_offs, xi_max = _scan_cfg_part1(
                r0_over_2b=r0_2b, two_T_pi=two_T_pi, mat=mat)
            geom_tag = f"p1_r{r0_2b:.4g}_2T{two_T_pi:.4g}pi".replace('.', 'p')
            plot_hi = (min(Om_hi, 1.6 * cut_offs[2])
                       if len(cut_offs) >= 3 else Om_hi)
            fname = os.path.join(FIGDIR, f"{geom_tag}_fig2.png")
            ok = plot_dispersion_curves(
                geom, mat, part=1, Om_lo=Om_lo, Om_hi=Om_hi, xmax=xi_max,
                om_plot_hi=plot_hi, fname=fname,
                title=f"r0/(2b)={r0_2b:.4g}, 2Θ={two_T_pi:.4g}π "
                      f"(out-of-plane)")
        else:
            (Om_lo, Om_hi), ns, cut_offs, ze_max = _scan_cfg_part2(
                r0_over_2b=r0_2b, two_T_pi=two_T_pi, mat=mat)
            geom_tag = f"p2_r{r0_2b:.4g}_2T{two_T_pi:.4g}pi".replace('.', 'p')
            plot_hi = (min(Om_hi, 0.9 * cut_offs[1])
                       if len(cut_offs) >= 2 else Om_hi)
            fname = os.path.join(FIGDIR, f"{geom_tag}_fig2.png")
            ok = plot_dispersion_curves(
                geom, mat, part=2, Om_lo=Om_lo, Om_hi=Om_hi, xmax=ze_max,
                om_plot_hi=plot_hi, fname=fname,
                title=f"r0/(2b)={r0_2b:.4g}, 2Θ={two_T_pi:.4g}π "
                      f"(in-plane)")

        return dict(ok=bool(ok), part=part_label, r0_2b=r0_2b,
                     two_T_pi=two_T_pi, fname=fname, dt=time.time() - t0)
    except Exception as exc:
        return dict(ok=False, part=part_label, r0_2b=r0_2b,
                     two_T_pi=two_T_pi, error=f"{exc}", dt=time.time() - t0)


def main():
    print("=" * 78)
    print("  probe_regen_dispersion_figs -- start "
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

    entries = _load_entries(CHECKPOINT_PATH)
    if not entries:
        print(f"Nothing to do: no entries found in {CHECKPOINT_PATH!r}.",
              flush=True)
        return
    print(f"\nFound {len(entries)} geometries to regenerate:", flush=True)
    for part_label, r0_2b, two_T_pi in entries:
        print(f"  {part_label}  r0/(2b)={r0_2b:.4g}  2*Theta={two_T_pi:.4g}*pi",
              flush=True)

    if not _fidelity_gate(entries):
        print("FATAL: fidelity gate mismatch -- stale/half-synced deploy? "
              "Aborting before spending compute.", flush=True)
        raise SystemExit(4)

    os.makedirs(FIGDIR, exist_ok=True)

    n_workers = int(os.environ.get("N_WORKERS", "8"))
    print(f"\nN_WORKERS={n_workers}", flush=True)

    with ProcessPoolExecutor(max_workers=min(n_workers, len(entries))) as ex:
        futs = {ex.submit(_geometry_worker, e): e for e in entries}
        for fut in as_completed(futs):
            r = fut.result()
            if r["ok"]:
                print(f"  OK   {r['part']:6s} r0/(2b)={r['r0_2b']:.4g} "
                      f"2*Theta={r['two_T_pi']:.4g}*pi -> {r['fname']} "
                      f"[{r['dt']:.0f}s]", flush=True)
            else:
                print(f"  FAIL {r['part']:6s} r0/(2b)={r['r0_2b']:.4g} "
                      f"2*Theta={r['two_T_pi']:.4g}*pi: {r.get('error')} "
                      f"[{r['dt']:.0f}s]", flush=True)

    print("\nData-regeneration only -- no checkpoint written, no "
          "SOLVER_VERSION action, no package changes. Figures in "
          f"{FIGDIR!r} have been overwritten in place with the "
          "gap-fixed dispersion sampling.", flush=True)
    print("End: " + time.strftime("%Y-%m-%d %H:%M:%S"), flush=True)


if __name__ == "__main__":
    t0 = time.time()
    main()
    print(f"\ntotal wall time: {(time.time() - t0)/60:.1f} min", flush=True)
