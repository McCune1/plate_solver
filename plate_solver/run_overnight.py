# -*- coding: utf-8 -*-
"""
plate_solver.run_overnight -- the full overnight batch driver (all Part 1 +
Part 2 geometries, checkpointed, run on the cluster).

Extracted verbatim (line-range provenance in LESSONS_LEARNED Sec. 22); no
numeric behaviour changed, SOLVER_VERSION not bumped.
"""
from __future__ import annotations
import os, sys, time, json, traceback
import numpy as np
from mpmath import mp, mpf, mpc, matrix

from .config import MP_DPS, N_WORKERS, PREFETCH, SOLVER_VERSION
from .geometry import IsotropicMaterial, make_geometry
from .detectors import find_modes_sigmin, paper_seeds_for, seed_polish_modes
from .validation import compare_and_collect, print_summary, PAPER_PART1, \
    PAPER_PART2, RR_CHECK_PART1, RR_CHECK_PART2, preflight, \
    _cantilever_spotcheck, _solve_freefree, compare_freefree_literature, \
    _orthotropic_scaffold_test, _print_literature_plan, _run_quality_review
from .plotting import plot_dispersion_curves, plot_mode_shape
from .workers import _StaleCheckpoint

def run_overnight(mat):
    """Top-level overnight orchestration (see header).  Best-effort throughout."""
    import time as _t, pickle, json
    t0 = _t.time()
    FAST = os.environ.get("FAST", "0") == "1"
    FIGDIR = os.environ.get("FIGDIR", "figures")
    try:
        os.makedirs(FIGDIR, exist_ok=True)
    except Exception:
        FIGDIR = "."
    ckpt_path = os.path.join(FIGDIR, "freefree_checkpoint.pkl")
    n_modes = int(os.environ.get("R40_FF_MODES", "6"))

    def _parse_geoms(env, default):
        spec = os.environ.get(env, default)
        out = []
        for tok in spec.split(","):
            tok = tok.strip()
            if not tok:
                continue
            a, b = tok.split(":")
            out.append((float(a), float(b)))
        return out

    print("\n" + "#" * 78)
    print("#  RESEARCH40 OVERNIGHT ALL-STEPS RUN")
    print(f"#  free-free modes/geom={n_modes}  FAST={FAST}  workers={N_WORKERS}")
    print("#" * 78)

    # Geometry lists parsed once so pre-flight, sweep, and review all agree.
    p1 = _parse_geoms("R40_FF_P1",
                      "1.25:0.5,1.25:1.0,1.25:1.5,1.25:0.25,1.25:0.75,1.25:1.25")
    p2 = _parse_geoms("R40_FF_P2", "1.25:0.5,1.25:1.0,1.25:0.25")

    # (0) PRE-FLIGHT — env manifest, scan plan, reference regression, parallel
    #     conditioning probe, wall-time estimate.  Catches a mis-configured or
    #     under-resolved sweep in minutes instead of after hours.  R40_PREFLIGHT=0
    #     to skip (e.g. a pure RESUME that has already passed it).
    if os.environ.get("R40_PREFLIGHT", "1") == "1":
        try:
            preflight(mat, p1, p2, n_modes, FAST)
        except Exception as exc:
            import traceback
            print(f"  WARNING: preflight failed (continuing): "
                  f"{type(exc).__name__}: {exc}")
            traceback.print_exc()

    # (1) cantilever spot-check ------------------------------------------------
    if os.environ.get("R40_SPOTCHECK", "1") == "1":
        try:
            _cantilever_spotcheck(mat)
        except Exception as exc:
            print(f"  WARNING: spot-check failed: {type(exc).__name__}: {exc}")

    # (2) free-free sweep (parallel) ------------------------------------------
    ff = {}
    if os.path.exists(ckpt_path):
        try:
            with open(ckpt_path, "rb") as f:
                ff = pickle.load(f)
            print(f"\n  (resume) loaded {len(ff)} free-free result(s) from checkpoint")
        except Exception:
            ff = {}

    if os.environ.get("R40_FREEFREE", "1") == "1":
        print("\n" + "=" * 70)
        print("  (2) FREE-FREE SPECTRA SWEEP  (both radial edges free; Step-6 parallel)")
        print("=" * 70)
        for (r0_2b, two_T) in p1:
            key = ("FF-P1", r0_2b, two_T)
            if key in ff and "error" not in ff[key]:
                print(f"  [skip] {key} already in checkpoint"); continue
            ff[key] = _solve_freefree(1, r0_2b, two_T, mat, n_modes, FAST)
            try:
                with open(ckpt_path, "wb") as f: pickle.dump(ff, f)
            except Exception as exc:
                print(f"  WARNING: checkpoint write failed: {exc}")
        for (r0_2b, two_T) in p2:
            key = ("FF-P2", r0_2b, two_T)
            if key in ff and "error" not in ff[key]:
                print(f"  [skip] {key} already in checkpoint"); continue
            ff[key] = _solve_freefree(2, r0_2b, two_T, mat, n_modes, FAST)
            try:
                with open(ckpt_path, "wb") as f: pickle.dump(ff, f)
            except Exception as exc:
                print(f"  WARNING: checkpoint write failed: {exc}")

        # consolidated free-free table (raw + Omega_lit) -----------------------
        print("\n" + "=" * 70)
        print("  FREE-FREE CONSOLIDATED TABLE")
        print("=" * 70)
        print(f"  {'part':6} {'r0/2b':>6} {'2T/pi':>6} {'mode':>4} "
              f"{'raw':>12} {'f[Hz]':>12} {'Omega_lit':>12}")
        txt_path = os.path.join(FIGDIR, "freefree_results.txt")
        lines = []
        for key in sorted(ff):
            tag, r0_2b, two_T = key
            rec = ff[key]
            if "error" in rec:
                lines.append(f"  {tag} {r0_2b} {two_T}pi  ERROR: {rec['error']}")
                continue
            for n, (rw, fz, om) in enumerate(zip(rec["raw"], rec["f_hz"],
                                                 rec["omega_lit"]), 1):
                line = (f"  {tag:6} {r0_2b:6.3g} {two_T:6.3g} {n:4d} "
                        f"{rw:12.6f} {fz:12.3f} {om:12.4f}")
                lines.append(line)
        for ln in lines:
            print(ln)
        try:
            with open(txt_path, "w") as f:
                f.write("# Research40 free-free spectra\n")
                f.write("# Omega_lit = w*R_o^2*sqrt(rho*H/D), H=2h, D=E*H^3/[12(1-nu^2)]\n")
                f.write("# part r0/2b 2T/pi mode raw f_Hz Omega_lit\n")
                f.write("\n".join(lines) + "\n")
            print(f"\n  wrote {txt_path}")
        except Exception as exc:
            print(f"  WARNING: could not write results file: {exc}")

        # post-run quality review + structured JSON manifest
        try:
            _run_quality_review(ff, FIGDIR)
        except Exception as exc:
            print(f"  WARNING: quality review failed: {type(exc).__name__}: {exc}")

    # (3) literature plan ------------------------------------------------------
    try:
        _print_literature_plan()
    except Exception as exc:
        print(f"  WARNING: literature plan print failed: {exc}")
    try:
        compare_freefree_literature(ff, mat)
    except Exception as exc:
        print(f"  WARNING: literature comparison failed: {exc}")

    # (4) orthotropic scaffold -------------------------------------------------
    if os.environ.get("R40_ORTHO_TEST", "1") == "1":
        try:
            _orthotropic_scaffold_test(mat)
        except Exception as exc:
            print(f"  WARNING: orthotropic scaffold failed: {exc}")

    print("\n" + "#" * 78)
    print(f"#  OVERNIGHT RUN COMPLETE in {(_t.time()-t0)/3600:.2f} h")
    print("#" * 78)


