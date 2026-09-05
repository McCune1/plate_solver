# -*- coding: utf-8 -*-
"""
probe_paper_figures_oop.py -- DIAGNOSTIC/DATA-GENERATION ONLY. Writes only
PNG figures and a small text log; no package changes, no SOLVER_VERSION
implications.

CONTEXT: PAPER1_FREEFREE_DRAFT.tex Sec.5/Sec.6 carry two open \todo figure
requests (added this session, per external feedback on the draft):
  (1) a sigma_min(Omega) trace with a physical zero and a spurious zero
      annotated, plus their edge-residual profiles;
  (2) mode shapes of 3-4 representative physical modes plotted alongside
      1-2 spurious zeros at similar Omega, to show visually how the
      spurious shapes fail the pointwise edge-residual screen where the
      physical ones pass it.
Both are reconstructable from already-validated Omega values with NO new
discovery search -- plate_solver.plotting.plot_mode_shape (Fig-3 style) and
a direct sigma_min sweep already exist and only need pointing at known
Omega's. IMPORTANT SCOPE LIMIT (checked against plotting.py's own
docstring 2026-07-14): mode_shape_grid/plot_mode_shape is implemented for
OutOfPlaneSolver (Part 1) ONLY -- Part 2 (in-plane) mode-shape
reconstruction is not yet supported and would silently return False. This
probe is therefore OOP-only; an in-plane figure job is a separate,
not-yet-buildable future task.

Geometry/material (FF-P1, same pinned instance used throughout the
free-free closure campaign and the convergence-study probe):
  mat  = IsotropicMaterial(E=210e9, nu=0.30, rho=7800.0)
  geom = make_geometry(1.5, 0.5)      # r0/(2b)=1.5, 2*Theta/pi=0.5
  solver = OutOfPlaneSolver(geom, mat, M=80, n_quad=30, boundary=FreeFreeOOP())

Targets (all already-validated raw Omega, LESSONS_LEARNED Sec 10.8/20.1,
Sec 20.4 -- none seed a NEW search, each is a direct reconstruction at a
known root):
  PHYSICAL: 0.383695 (fundamental), 1.369611 (mid-spectrum), 2.602500
    (highest validated real mode) -- a representative spread across the
    8-mode validated OOP table.
  ARTIFACT: 2.817613 (EVEN parity, ARTIFACT call parity-confirmed, Sec
    20.4) -- deliberately chosen adjacent to the still-open ODD-parity
    ceiling ambiguity (raw 2.828013 / 2.838515) so the sigma_min trace
    figure shows a real mode, a confirmed artifact, AND the ambiguous
    region in one coherent window, tying the figure directly to the
    paper's own "one bounded ambiguity" discussion (Sec 6.3).

WHAT THIS DOES:
  A. For each of the 4 target Omega (3 physical + 1 artifact): full_search
     + select_fill to get the branch/basis set at n_dofs=20 (validated
     default), then plot_mode_shape -> one PNG each (4 total). This picks
     up plotting.py's current MODE_SHAPE_Z_RATIO (default 0.18 as of the
     2026-07-17 fix) automatically -- no code change needed here, this is
     an unmodified re-run to refresh figures that predate that fix.
  B. A sigma_min(Omega) trace over Omega in [2.55, 2.90] (37 points,
     n_dofs=20) spanning the last validated real mode (2.602500), the
     confirmed artifact (2.817613), and the open ambiguity pair (2.828013,
     2.838515) -- saved as a simple matplotlib line plot with all four
     Omega's annotated. (Unaffected by the mode-shape fix; regenerated
     only because it shares this script/job.)
Both (A) and (B) reuse already-validated Omega values -- no new mode
search, no threshold tuning, no SOLVER_VERSION implications. Figures are
written to ./paper_figures/ as PNGs; a plain-text log of every sigma_min
value is also written so the trace can be re-plotted/re-styled later
without re-running the (expensive) sweep.

COST: 4 mode-shape reconstructions (each ~1 full_search + 1 grid build,
~90-150s) + a 37-point sigma_min sweep (~90-150s/point, same cost driver:
one fresh full_search per point) = ~41 "full_search units" total, i.e. the
SAME per-unit cost as the free-free convergence study. Parallelized across
workers sized to the SLURM allocation (each point/target independent, own
solver built from scalars). Budget generously; not time-critical.
"""
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

sys.path.insert(0, os.environ.get("PKG_PATH", "."))
os.environ.setdefault("DPS", "40")

EXPECT_VER = os.environ.get("EXPECT_SOLVER_VERSION", "2026-07-10.s10")
OUTDIR = os.environ.get("OUTDIR", "paper_figures")

PHYSICAL_TARGETS = [
    ("phys_fundamental", 0.383695, 1),
    ("phys_midspectrum", 1.369611, 4),
    ("phys_highest", 2.602500, 8),
]
ARTIFACT_TARGETS = [
    ("artifact_even_2p8176", 2.817613, "A1"),
]
TRACE_LO, TRACE_HI, TRACE_N = 2.55, 2.90, 37
N_DOFS = 20
XI_MAX = 20.0  # matches the FF-P1 n_dofs=20 default search radius elsewhere


def _mode_shape_worker(args):
    label, Om, mode_num = args
    t0 = time.time()
    try:
        import plate_solver as ps
        from plate_solver.detectors import full_search, select_fill
        from plate_solver import plotting

        mat = ps.IsotropicMaterial(E=210e9, nu=0.30, rho=7800.0)
        geom = ps.make_geometry(1.5, 0.5)
        solver = ps.OutOfPlaneSolver(geom, mat, M=80, n_quad=30,
                                      boundary=ps.FreeFreeOOP())
        raw = full_search(solver.fast, Om, xmax=XI_MAX)
        sel, cnt = select_fill(raw, N_DOFS)
        if cnt < N_DOFS:
            return dict(ok=False, label=label, error=f"only {cnt}/{N_DOFS} "
                        f"branches filled", dt=time.time() - t0)
        fname = os.path.join(OUTDIR, f"modeshape_{label}.png")
        ok = plotting.plot_mode_shape(solver, Om, sel, N_DOFS, mode_num,
                                       "FF-P1 OOP", fname)
        return dict(ok=bool(ok), label=label, fname=fname,
                     dt=time.time() - t0)
    except Exception as exc:
        return dict(ok=False, label=label, error=f"{exc}",
                     dt=time.time() - t0)


def _sigmin_worker(Om):
    t0 = time.time()
    try:
        import plate_solver as ps
        from plate_solver.detectors import full_search, select_fill

        mat = ps.IsotropicMaterial(E=210e9, nu=0.30, rho=7800.0)
        geom = ps.make_geometry(1.5, 0.5)
        solver = ps.OutOfPlaneSolver(geom, mat, M=80, n_quad=30,
                                      boundary=ps.FreeFreeOOP())
        raw = full_search(solver.fast, Om, xmax=XI_MAX)
        sel, cnt = select_fill(raw, N_DOFS)
        if cnt < N_DOFS:
            return dict(ok=False, Om=Om, error=f"only {cnt}/{N_DOFS}",
                        dt=time.time() - t0)
        s = solver.sigma_min(Om, sel)
        return dict(ok=True, Om=Om, sigma_min=float(s), dt=time.time() - t0)
    except Exception as exc:
        return dict(ok=False, Om=Om, error=f"{exc}", dt=time.time() - t0)


def main():
    print("=" * 78)
    print("  probe_paper_figures_oop -- start " + time.strftime("%Y-%m-%d %H:%M:%S"))
    print("=" * 78, flush=True)

    import plate_solver as ps
    if ps.SOLVER_VERSION != EXPECT_VER:
        print(f"PREFLIGHT FAIL: SOLVER_VERSION={ps.SOLVER_VERSION!r} != "
              f"{EXPECT_VER!r}")
        raise SystemExit(2)
    print(f"preflight OK: SOLVER_VERSION={ps.SOLVER_VERSION}", flush=True)
    os.makedirs(OUTDIR, exist_ok=True)

    n_workers = int(os.environ.get("N_WORKERS", "8"))
    print(f"N_WORKERS={n_workers}", flush=True)

    # ---- Part A: mode-shape PNGs (physical + artifact) --------------------
    shape_jobs = [(lbl, Om, mn) for lbl, Om, mn in PHYSICAL_TARGETS] + \
                 [(lbl, Om, mn) for lbl, Om, mn in ARTIFACT_TARGETS]
    print(f"\nPart A: {len(shape_jobs)} mode-shape reconstructions "
          f"(3 physical + 1 artifact)", flush=True)
    with ProcessPoolExecutor(max_workers=min(n_workers, len(shape_jobs))) as ex:
        futs = {ex.submit(_mode_shape_worker, j): j for j in shape_jobs}
        for fut in as_completed(futs):
            r = fut.result()
            if r["ok"]:
                print(f"  OK  {r['label']:24s} -> {r['fname']} [{r['dt']:.0f}s]",
                      flush=True)
            else:
                print(f"  FAIL {r['label']:24s}: {r.get('error')} "
                      f"[{r['dt']:.0f}s]", flush=True)

    # ---- Part B: sigma_min(Omega) trace ------------------------------------
    # SKIP_TRACE=1: mode-shape PDFs only. The 37-point sweep is already
    # shipped as paper_figures/sigmin_trace_2p55_2p90.{txt,pdf}.
    if os.environ.get("SKIP_TRACE", "0") == "1":
        print("\nPart B: skipped (SKIP_TRACE=1)", flush=True)
        print("\nDiagnostic/data-generation only -- no SOLVER_VERSION action, "
              "no package changes.", flush=True)
        print("End: " + time.strftime("%Y-%m-%d %H:%M:%S"), flush=True)
        return

    Oms = [TRACE_LO + i * (TRACE_HI - TRACE_LO) / (TRACE_N - 1)
           for i in range(TRACE_N)]
    print(f"\nPart B: sigma_min trace, {TRACE_N} points in "
          f"[{TRACE_LO}, {TRACE_HI}]", flush=True)
    trace = []
    with ProcessPoolExecutor(max_workers=n_workers) as ex:
        futs = {ex.submit(_sigmin_worker, Om): Om for Om in Oms}
        done = 0
        for fut in as_completed(futs):
            r = fut.result()
            done += 1
            if r["ok"]:
                trace.append((r["Om"], r["sigma_min"]))
                print(f"  [{done}/{TRACE_N}] Om={r['Om']:.4f} "
                      f"log10(sigma_min)={r['sigma_min']:.4f} [{r['dt']:.0f}s]",
                      flush=True)
            else:
                print(f"  [{done}/{TRACE_N}] Om={r['Om']:.4f} FAILED: "
                      f"{r.get('error')} [{r['dt']:.0f}s]", flush=True)

    trace.sort(key=lambda t: t[0])
    log_path = os.path.join(OUTDIR, "sigmin_trace_2p55_2p90.txt")
    with open(log_path, "w") as f:
        f.write("# Omega  log10_sigma_min  (FF-P1 OOP, n_dofs=20)\n")
        f.write("# annotate: 2.602500=last validated REAL, "
                "2.817613=confirmed EVEN ARTIFACT, "
                "2.828013/2.838515=open ODD ambiguity pair\n")
        for Om, s in trace:
            f.write(f"{Om:.6f}  {s:.6f}\n")
    print(f"\nWrote sigma_min trace data: {log_path}", flush=True)

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(7, 4))
        xs = [t[0] for t in trace]
        ys = [t[1] for t in trace]
        ax.plot(xs, ys, "k-", linewidth=1.0)
        for Om, tag in [(2.602500, "REAL"), (2.817613, "ARTIFACT"),
                        (2.828013, "ODD cand. 1"), (2.838515, "ODD cand. 2")]:
            ax.axvline(Om, color="0.6", linestyle="--", linewidth=0.7)
            ax.annotate(tag, (Om, max(ys)), rotation=90, fontsize=7,
                        va="top", ha="right")
        ax.set_xlabel(r"$\Omega$")
        ax.set_ylabel(r"$\log_{10}\sigma_{\min}(K)$")
        ax.set_title("FF-P1 OOP: real mode, artifact, and the open ambiguity")
        fig.tight_layout()
        png_path = os.path.join(OUTDIR, "sigmin_trace_2p55_2p90.png")
        pdf_path = os.path.join(OUTDIR, "sigmin_trace_2p55_2p90.pdf")
        fig.savefig(png_path, dpi=150)
        fig.savefig(pdf_path, format="pdf", bbox_inches="tight")
        print(f"Wrote sigma_min trace figure: {png_path}", flush=True)
        print(f"Wrote sigma_min trace vector: {pdf_path}", flush=True)
    except Exception as exc:
        print(f"WARNING: could not render trace PNG ({exc}); the .txt data "
              f"file above is still usable to plot locally.", flush=True)

    print("\nDiagnostic/data-generation only -- no SOLVER_VERSION action, "
          "no package changes. Figures/data are inputs to "
          "PAPER1_FREEFREE_DRAFT.tex Sec.5/6's figure TODOs.", flush=True)
    print("End: " + time.strftime("%Y-%m-%d %H:%M:%S"), flush=True)


if __name__ == "__main__":
    t0 = time.time()
    main()
    print(f"\ntotal wall time: {(time.time() - t0)/60:.1f} min", flush=True)
