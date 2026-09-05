# -*- coding: utf-8 -*-
"""
plate_solver.cli -- command-line entry point.

USABILITY ENHANCEMENT (review-doc item): adds an argparse layer over the
original environment-variable-driven __main__ block, WITHOUT changing what
the env vars do or their defaults -- every argparse flag below just sets the
corresponding os.environ value before main() runs its original (unmodified)
env-var-reading logic. This means:
  * existing submit_*.sh scripts that export env vars keep working unchanged;
  * `python -m plate_solver.cli --fast --run-p2 0` also works now;
  * `python -m plate_solver.cli` with NO flags is byte-identical to the old
    `python Research50.py` with no env vars set.

Extracted verbatim from Research50.py's `if __name__ == "__main__":` block
(now `main()`) plus rect_int.py's R40_RECT/R40_RECT_IP dispatch, inserted at
the same point rect_int.py inserted it (line-range provenance in
LESSONS_LEARNED Sec. 22). No numeric behaviour changed, SOLVER_VERSION not
bumped.
"""
from __future__ import annotations
import argparse
import os, time
import numpy as np
from mpmath import mp, mpf

from .config import SOLVER_VERSION, N_WORKERS, MAX_SCAN_PTS, PI
from .geometry import make_geometry, _material_from_env, \
    FGMPlateProperties
from .dispersion import cutoff_frequencies_part1, cutoff_frequencies_part2, \
    _Part1Fast
from .core_solvers import OutOfPlaneSolver, InPlaneSolver
from .detectors import find_modes_sigmin, paper_seeds_for, seed_polish_modes, \
    _selftest_part1, _selftest_part2, full_search, classify_modes_ritz
from .plotting import plot_dispersion_curves, plot_mode_shape
from .validation import compare_and_collect, print_summary, PAPER_PART1, \
    PAPER_PART2, run_shi_validation, run_mcgee_validation, run_rect_validation, \
    run_rect_ip_validation, _scan_cfg_part1, _scan_cfg_part2
from .run_overnight import run_overnight
from .workers import _StaleCheckpoint


# ── argparse layer: flag -> env var name, plus the arg's own env encoding ───
# Only the flags a person is actually likely to reach for at the command line
# are exposed here; every other env-var knob catalogued in LESSONS_LEARNED
# stays reachable exactly as before (`FOO=1 python -m plate_solver.cli`).
def _build_argparser():
    p = argparse.ArgumentParser(
        prog="plate_solver",
        description="Annular & Rectangular Sector Plate Vibration Solver "
                     "(Seok & Tiersten 2004 reproduction + extensions). "
                     "Every flag below is equivalent to setting the same-named "
                     "environment variable; unset flags leave existing "
                     "env-var/default behaviour untouched.")
    p.add_argument("--fast", action="store_true", default=None, help="FAST=1: quick low-resolution run (testing only)")
    p.add_argument("--run-p1", type=str, default=None, choices=["0", "1"], help="RUN_P1")
    p.add_argument("--run-p2", type=str, default=None, choices=["0", "1"], help="RUN_P2")
    p.add_argument("--geoms", type=str, default=None, help="GEOMS: comma-separated geometry subset filter")
    p.add_argument("--dps", type=int, default=None, help="DPS: mpmath decimal places (min 30)")
    p.add_argument("--n-workers", type=int, default=None, help="N_WORKERS: worker process count")
    p.add_argument("--max-scan-pts", type=int, default=None, help="MAX_SCAN_PTS")
    p.add_argument("--mat-e", type=float, default=None, help="MAT_E: Young's modulus (Pa)")
    p.add_argument("--mat-nu", type=float, default=None, help="MAT_NU: Poisson ratio")
    p.add_argument("--mat-rho", type=float, default=None, help="MAT_RHO: density (kg/m^3)")
    p.add_argument("--checkpoint", type=str, default=None, help="CHECKPOINT: path for the JSON progress checkpoint")
    p.add_argument("--resume", type=str, default=None, choices=["0", "1"], help="RESUME")
    p.add_argument("--figdir", type=str, default=None, help="FIGDIR: output directory for figures")
    p.add_argument("--make-figs", type=str, default=None, choices=["0", "1"], help="MAKE_FIGS")
    p.add_argument("--shi", action="store_true", default=None, help="R40_SHI=1: run the Shi/Lv 2016 FFFF orthotropic annular validation and exit")
    p.add_argument("--mcgee", action="store_true", default=None, help="R40_MCGEE=1: run the McGee free-free validation and exit")
    p.add_argument("--rect", action="store_true", default=None, help="R40_RECT=1: run the rectangular cantilever OOP validation and exit")
    p.add_argument("--rect-ip", action="store_true", default=None, help="R40_RECT_IP=1: run the rectangular cantilever in-plane validation and exit")
    return p


def _apply_args_to_env(args):
    mapping = {
        "run_p1": "RUN_P1", "run_p2": "RUN_P2", "geoms": "GEOMS",
        "dps": "DPS", "n_workers": "N_WORKERS", "max_scan_pts": "MAX_SCAN_PTS",
        "mat_e": "MAT_E", "mat_nu": "MAT_NU", "mat_rho": "MAT_RHO",
        "checkpoint": "CHECKPOINT", "resume": "RESUME", "figdir": "FIGDIR",
        "make_figs": "MAKE_FIGS",
    }
    for attr, env in mapping.items():
        val = getattr(args, attr)
        if val is not None:
            os.environ[env] = str(val)
    if args.fast:
        os.environ["FAST"] = "1"
    if args.shi:
        os.environ["R40_SHI"] = "1"
    if args.mcgee:
        os.environ["R40_MCGEE"] = "1"
    if args.rect:
        os.environ["R40_RECT"] = "1"
    if args.rect_ip:
        os.environ["R40_RECT_IP"] = "1"


def main(argv=None):
    _args = _build_argparser().parse_args(argv)
    _apply_args_to_env(_args)

    # ── Runtime configuration ─────────────────────────────────────────────────
    # All settings can be overridden by environment variables without editing
    # this file.  See the module docstring for the full list.
    # The float64 fast engine can overflow at low Ω / large r0_bar; track() now
    # drops the resulting garbage roots, so silence the (handled) numpy warnings
    # to keep the run log readable.
    np.seterr(over="ignore", invalid="ignore", divide="ignore")
    FAST    = os.environ.get("FAST",   "0") == "1"
    RUN_P1  = os.environ.get("RUN_P1", "1") == "1"
    RUN_P2  = os.environ.get("RUN_P2", "1") == "1"
    GEOMS_ENV = os.environ.get("GEOMS", "")
    # Mode finding is FULLY GENERAL by default: every geometry is solved from
    # scratch by locating zeros of the characteristic determinant det K(Ω)=0 —
    # the same quantity Seok & Tiersten root-find in the paper.  The published
    # PAPER_PART1/PAPER_PART2 tables are used ONLY for the post-hoc accuracy
    # report (compare_and_collect) and have ZERO influence on any computed
    # frequency.  REF_SEED=1 is an opt-in validation aid that seeds the search
    # near published values; it is OFF by default so results independently
    # reproduce the paper as if it had never been seen.
    REF_SEED = os.environ.get("REF_SEED", "0") == "1"
    # Paper-faithful detector (default ON): natural frequencies are the zeros of
    # the characteristic determinant, located as singularities of K via σ_min(K)
    # — robust, sign-free, no per-mode tuning.  PAPER_FAITHFUL=0 reverts to the
    # older log|det|-dip detector.  Ignored when REF_SEED=1.
    PAPER_FAITHFUL = os.environ.get("PAPER_FAITHFUL", "1") == "1"

    print("=" * 70)
    print("  ANNULAR SECTOR PLATE SOLVER — Seok & Tiersten 2004 (Research23)")
    if REF_SEED:
        print("  MODE: reference-seeded validation (REF_SEED=1) — search windows")
        print("        seeded from published values; set REF_SEED=0 for independent run")
    else:
        print("  MODE: fully general / independent — frequencies are zeros of the")
        print("        characteristic determinant det K(Ω)=0 found from scratch.")
        if PAPER_FAITHFUL:
            print("        Detector: σ_min(K) singularities (paper-faithful).")
        else:
            print("        Detector: log|det| dips (legacy).")
        print("        Paper tables are used only for the post-hoc accuracy report.")
    print(f"  mp.dps={mp.dps}  M=80  N_WORKERS={N_WORKERS}  "
          f"MAX_SCAN_PTS={MAX_SCAN_PTS}  FAST={FAST}")
    print("=" * 70)

    # ── Material ────────────────────────────────────────────────────────────
    # Default = the validation material (isotropic, ν=0.35) used for all 29
    # paper cases; keep it for any correctness check against the paper
    # (LESSONS §0).  MAT_E / MAT_NU / MAT_RHO let a literature-matched run
    # (e.g. ν=0.30 for the McGee/Shi FFFF comparison, LESSONS §13) be launched
    # without editing source.  Defaults reproduce ν=0.35 BIT-IDENTICALLY, so
    # existing runs are unchanged and SOLVER_VERSION is NOT bumped.
    mat = _material_from_env()

    # ── Self-tests on the canonical geometry (r0/(2b)=1.25, 2Theta=pi) ───────
    # These verify that the Frobenius recursions satisfy the governing ODEs to
    # near machine-precision; they should always pass.
    geom_ref  = make_geometry(1.25, 1.0)
    k_bar_ref = float((mpf(2)*geom_ref.b/mp.pi) *
                      mp.sqrt(3*mat.c66/(geom_ref.h**2*mat.c11_bar)))
    w_bar_ref = float((mp.pi/(2*geom_ref.b))*mp.sqrt(mat.c66/mat.rho))
    print(f"\n{geom_ref}\n{mat}")
    print(f"k_bar={k_bar_ref:.4f}  omega_bar={w_bar_ref:.4e} rad/s\n")

    # ── Demonstrate cut-off frequency computation ─────────────────────────────
    print("─" * 70)
    print("Cut-off frequency demo (ξ=0 / ζ=0 Bessel determinant roots)")
    print("─" * 70)
    r0b_ref = float(np.pi * 1.25)   # r̄₀ = π * r₀/(2b)
    _nu_d = float(mat.nu_bar)
    co_p1 = cutoff_frequencies_part1(r0b_ref, _nu_d, n_co=6)
    co_p2 = cutoff_frequencies_part2(r0b_ref, _nu_d, 2.0/(1-_nu_d), n_co=6)
    print(f"  Part 1 cut-off Ω̃: {[f'{c:.5f}' for c in co_p1]}")
    print(f"  Part 2 cut-off Ω̄: {[f'{c:.5f}' for c in co_p2]}")

    # ── Rectangular cantilever OOP validation (R40_RECT=1) ───────────────────
    #    Stage-2 integration (Paper A): exact-edge RectangularCartesianOOP +
    #    Eq.53 variational K + two-sided conjugate-pair realification.  Isolated
    #    from the annular path; dispatches before the annular ODE self-tests.
    if os.environ.get("R40_RECT", "0") == "1":
        run_rect_validation()
        raise SystemExit(0)

    # ── Rectangular cantilever IN-PLANE validation (R40_RECT_IP=1) ───────────
    #    Stage-2 (Paper A, Part 2): exact-edge RectangularCartesianIP + Eq.41
    #    unconstrained variational K; scale-invariant zeta branch-point test.
    #    Single-threaded mpmath, isolated from the annular path (no MPI/pool).
    if os.environ.get("R40_RECT_IP", "0") == "1":
        run_rect_ip_validation()
        raise SystemExit(0)

    print("─" * 70)
    print("ODE self-tests")
    print("─" * 70)
    ok1 = _selftest_part1(geom_ref, mat)
    ok2 = _selftest_part2(geom_ref, mat)
    if not (ok1 and ok2):
        print("ERROR: self-tests failed — check Frobenius recursion code.")
        raise SystemExit(1)

    # ── Overnight all-steps driver (default): cantilever spot-check +
    #    parallel free-free spectra + orthotropic scaffold.  Set
    #    R40_OVERNIGHT=0 to fall through to the legacy cantilever sweep.
    # ── Step-1 Shi/Lv 2016 FFFF orthotropic benchmark (R40_SHI=1) ────────────
    #    DIRECT annular comparison (no extrapolation); builds its own
    #    OrthotropicMaterial from SHI_* env (default = paper §3 constants).
    if os.environ.get("R40_SHI", "0") == "1":
        run_shi_validation()
        raise SystemExit(0)

    # ── Step-1 McGee solid-sector limit validation (R40_MCGEE=1) ─────────────
    if os.environ.get("R40_MCGEE", "0") == "1":
        run_mcgee_validation(mat)
        raise SystemExit(0)

    if os.environ.get("R40_OVERNIGHT", "1") == "1":
        run_overnight(mat)
        raise SystemExit(0)

    # ── Reference-seeded polish smoke test (fail fast, ~1-2 min) ─────────────
    # Exercises the seeded-convergence path on one well-conditioned mode so a
    # bug surfaces in minutes instead of after a multi-hour run.
    # Default OFF (consistent with SELFTEST_SIGMIN) for overnight/overwork
    # runs; re-enable with SELFTEST_SEED=1 when testing a fresh code change.
    # Only relevant at all when REF_SEED=1 (off by default).
    if REF_SEED and os.environ.get("SELFTEST_SEED", "0") == "1":
        print("─" * 70)
        print("Reference-seeded polish smoke test (Part 1, r0/(2b)=1.25, 2Θ=π, mode 3)")
        print("─" * 70)
        try:
            oop_ref = OutOfPlaneSolver(geom_ref, mat, M=80, n_quad=30)
            got = seed_polish_modes(oop_ref, [0.24312], n_dofs=16, max_dim=14.0,
                                    part=1, scan_lo=0.005, scan_hi=0.40,
                                    iters=8, verbose=True)
            if not got or not np.isfinite(got[0]) or not (0.18 < got[0] < 0.32):
                raise RuntimeError(f"unexpected result {got} (expected ~0.244)")
            print(f"  smoke test OK → {got[0]:.6f}\n", flush=True)
        except Exception as exc:
            print(f"ERROR: seeded-polish smoke test failed: {exc}")
            print("Aborting before the long run.  Re-run with REF_SEED=0 to use "
                  "the from-scratch detector, or SELFTEST_SEED=0 to skip this check.")
            raise SystemExit(1)

    # ── σ_min detector smoke test (fail fast, ~2-4 min) ──────────────────────
    # Full parallel σ_min path on the canonical geometry (mode 1 near Ω≈0.0332).
    # Default OFF; SELFTEST_SIGMIN=1.  Window caveat: LESSONS_LEARNED.md §Detector.
    if (not REF_SEED) and PAPER_FAITHFUL and \
            os.environ.get("SELFTEST_SIGMIN", "0") == "1":
        print("─" * 70)
        print("σ_min detector smoke test (Part 1, r0/(2b)=1.25, 2Θ=π) — expect a")
        print("mode near Ω≈0.0332 …")
        print("─" * 70)
        try:
            oop_st = OutOfPlaneSolver(geom_ref, mat, M=80, n_quad=30)
            got = find_modes_sigmin(
                oop_st, part=1, Omega_range=(0.029, 0.039), n_scan=14,
                n_dofs=16, max_dim=14.0, n_modes_wanted=1, verbose=True)
            if not got or not any(0.028 < g < 0.039 for g in got):
                raise RuntimeError(f"unexpected result {got} (expected ~0.0332)")
            print(f"  smoke test OK → {got[0]:.6f}\n", flush=True)
        except SystemExit:
            raise
        except Exception as exc:
            print(f"ERROR: σ_min detector smoke test failed: {exc}")
            print("Aborting before the long run.  Re-run with PAPER_FAITHFUL=0 to "
                  "use the legacy detector, or SELFTEST_SIGMIN=0 to skip this check.")
            raise SystemExit(1)

    # ── Branch inventory at the canonical mode-1 frequency ───────────────────
    # Prints the canonical xi-roots used to build the K matrix basis at a known
    # frequency.  Useful for debugging branch-finding if modes go missing.
    print("\n" + "─" * 70)
    print("xi-branch inventory at Omega_tilde = 0.0332 (Part 1, mode 1)")
    print("─" * 70)
    f1_ref = _Part1Fast(float(geom_ref.r0_bar), float(mat.nu_bar), M=80)
    roots  = full_search(f1_ref, 0.0332)
    print("  Canonical roots (conjugate mirrors implied):")
    print("  " + ", ".join(f"{z:.4f}" for z in roots))

    # ── Build the list of geometries to run ───────────────────────────────────
    # All (r0/(2b), 2*Theta/pi) pairs present in PAPER_PART1
    p1_keys = sorted(
        {(r0, ang) for (r0, _), d in PAPER_PART1.items() for ang in d},
        key=lambda x: (x[0], x[1])
    )
    # All pairs present in PAPER_PART2
    p2_keys = sorted(
        {(r0, ang) for (r0, _), d in PAPER_PART2.items() for ang in d},
        key=lambda x: (x[0], x[1])
    )

    # Optional geometry filter from environment variable
    if GEOMS_ENV:
        allowed = set(GEOMS_ENV.split(","))
        p1_keys = [k for k in p1_keys if f"{k[0]:.4g}_{k[1]:.4g}" in allowed]
        p2_keys = [k for k in p2_keys if f"{k[0]:.4g}_{k[1]:.4g}" in allowed]

    # summary collects results for the consolidated table printed at the end
    summary = []

    # results collects raw (Om values + geometry) for plotting (Figs. 2 & 3)
    # and for the JSON checkpoint, keyed by (part, r0_2b, two_T_pi).
    results = {}
    CHECKPOINT_PATH = os.environ.get("CHECKPOINT", "research23_checkpoint.json")
    MAKE_FIGS = os.environ.get("MAKE_FIGS", "1") == "1"
    FIGDIR = os.environ.get("FIGDIR", "figures")
    if MAKE_FIGS:
        # Check ONCE up front whether matplotlib is importable. If it isn't,
        # disable figure generation with a single clear message rather than
        # emitting a "No module named 'matplotlib'" warning for every figure
        # of every geometry (which buried the real output last run).
        try:
            import matplotlib  # noqa: F401
        except Exception:
            print("\n" + "!" * 70)
            print("  matplotlib is not installed — Figure 2/3 PNGs will be SKIPPED.")
            print("  To enable plotting, install it once and re-run:")
            print("      pip install matplotlib")
            print("  (Frequency results are unaffected; only the figures are skipped.)")
            print("!" * 70 + "\n", flush=True)
            MAKE_FIGS = False
    if MAKE_FIGS:
        try:
            os.makedirs(FIGDIR, exist_ok=True)
        except Exception as exc:
            print(f"  WARNING: could not create figure directory {FIGDIR}: {exc}")
            MAKE_FIGS = False

    def _save_checkpoint():
        """Write partial results to disk so an overnight run that crashes or
        is interrupted partway through still leaves usable output behind."""
        try:
            import json
            out = {
                "solver_version": SOLVER_VERSION,
                "summary": summary,
                "results": {
                    f"{part}|{r0_2b}|{two_T_pi}": freqs
                    for (part, r0_2b, two_T_pi), freqs in results.items()
                },
            }
            with open(CHECKPOINT_PATH, "w") as fh:
                json.dump(out, fh, indent=2)
        except Exception as exc:
            print(f"  WARNING: failed to write checkpoint: {exc}")

    # ── RESUME (default ON) ──────────────────────────────────────────────────
    # If a checkpoint from a previous (interrupted) run exists, reload its
    # completed geometries so the run PICKS UP WHERE IT LEFT OFF instead of
    # recomputing everything.  This is the practical remedy for a machine that
    # suspends/sleeps when the tab loses focus (see notes): just relaunch and
    # it continues.  Disable with RESUME=0 to force a clean from-scratch run.
    _done = set()   # set of (part_label, r0_2b, two_T_pi) already completed
    RESUME = os.environ.get("RESUME", "1") != "0"
    if RESUME and os.path.exists(CHECKPOINT_PATH):
        try:
            import json
            with open(CHECKPOINT_PATH) as fh:
                ckpt = json.load(fh)
            # VERSION GUARD: a checkpoint written by a different solver version
            # may contain results computed with now-fixed bugs.  Refuse to
            # resume from it (so the run recomputes everything with the current
            # code) instead of silently replaying stale results.  This is the
            # common "I edited the solver and reran, but it skipped everything"
            # trap — the old checkpoint is from before the edit.
            ckpt_ver = ckpt.get("solver_version")
            if ckpt_ver != SOLVER_VERSION:
                print(f"\n  RESUME: checkpoint {CHECKPOINT_PATH} was written by "
                      f"solver version '{ckpt_ver}', but this is "
                      f"'{SOLVER_VERSION}'.  IGNORING the stale checkpoint and "
                      f"recomputing from scratch (results from an older code "
                      f"version are not reused).  Delete the file or set "
                      f"RESUME=0 to silence this.\n", flush=True)
                raise _StaleCheckpoint()
            # results keys are saved as "<part_label>|<r0_2b>|<two_T_pi>" where
            # part_label is the STRING 'Part 1'/'Part 2' (matching how results
            # is keyed below) — do NOT int() it.
            for key, freqs in ckpt.get("results", {}).items():
                parts = key.split("|")
                if len(parts) != 3:
                    continue
                part_s, r_s, t_s = parts
                results[(part_s, float(r_s), float(t_s))] = freqs
            for row in ckpt.get("summary", []):
                summary.append(tuple(row))
            # mark completed geometries (by part label + geometry); a geometry
            # counts as done only if it produced at least one non-MISS row, so
            # an interrupted geometry that only logged MISS placeholders is
            # recomputed rather than skipped.
            blocks = {}
            for row in ckpt.get("summary", []):
                k = (row[0], float(row[1]), float(row[2]))
                blocks.setdefault(k, []).append(row[7])
            for k, statuses in blocks.items():
                if any(s != 'MISS' for s in statuses):
                    _done.add(k)
            if _done:
                print(f"\n  RESUME: loaded checkpoint {CHECKPOINT_PATH} — "
                      f"{len(_done)} geometry/part block(s) already complete "
                      f"will be skipped. (RESUME=0 to recompute from scratch.)\n",
                      flush=True)
        except _StaleCheckpoint:
            # version mismatch already reported above; ensure no stale state
            _done = set()
            results.clear()
            summary.clear()
        except Exception as exc:
            print(f"  WARNING: could not load checkpoint for resume: {exc}")
            _done = set()
            results.clear()
            summary.clear()

    # ── Part 1: out-of-plane frequencies ─────────────────────────────────────
    # Geometries are processed one at a time so that worker processes from one
    # geometry are fully released before the next begins.  This keeps peak RAM
    # usage bounded to a single geometry's worth of concurrent workers.
    if RUN_P1:
        print("\n" + "═" * 70)
        print("PART 1 — OUT-OF-PLANE FREQUENCIES (Omega_tilde = k_bar * Omega_bar)")
        print("  Scan ranges derived from xi=0 cut-off frequencies; mode selection")
        print("  is purely physics-based (no paper-table values used)")
        print("═" * 70)

        for geom_idx, (r0_2b, two_T_pi) in enumerate(p1_keys, 1):
            if ('Part 1', float(r0_2b), float(two_T_pi)) in _done:
                print(f"\n{'─' * 70}")
                print(f"  GEOMETRY {geom_idx}/{len(p1_keys)}  r0/(2b)={r0_2b:.4g}  "
                      f"2*Theta={two_T_pi}*pi — already in checkpoint, SKIPPING "
                      f"(RESUME).", flush=True)
                continue
            t_geom = time.time()
            try:
                geom  = make_geometry(r0_2b, two_T_pi)
                k_bar = float((mpf(2)*geom.b/mp.pi) *
                              mp.sqrt(3*mat.c66/(geom.h**2*mat.c11_bar)))
                w_bar = float((mp.pi/(2*geom.b))*mp.sqrt(mat.c66/mat.rho))

                print(f"\n{'─' * 70}")
                print(f"  GEOMETRY {geom_idx}/{len(p1_keys)}  "
                      f"r0/(2b)={r0_2b:.4g}  2*Theta={two_T_pi}*pi  |  {geom}")
                print(f"  k_bar={k_bar:.4f}  omega_bar={w_bar:.4e} rad/s")

                (Om_lo, Om_hi), ns, cut_offs, xi_max = _scan_cfg_part1(r0_over_2b=r0_2b, two_T_pi=two_T_pi, mat=mat)
                print(f"  cut-off Ω̃ (ξ=0): {[f'{c:.5f}' for c in cut_offs]}")
                if FAST:
                    Om_hi = min(Om_hi, Om_lo + (Om_hi - Om_lo) * 0.5)
                    ns    = min(ns, 20)

                oop   = OutOfPlaneSolver(geom, mat, M=80, n_quad=30)
                # Narrow sectors need more dispersion branches to converge
                # higher modes — the paper itself escalates from "C6" to
                # "C8" only for 2Θ=π/4 (its narrowest tabulated sector).
                # n_dofs=16 ≈ up to 7-8 real/imag branches or fewer complex
                # pairs; raise to 20 (≈8-9 branches) for 2Θ<=π/2, paired
                # with the wider xi_max already returned by _scan_cfg_part1.
                p1_n_dofs = 20 if two_T_pi <= 0.5 else 16
                p1_seeds = paper_seeds_for(PAPER_PART1, r0_2b, two_T_pi)
                if REF_SEED and p1_seeds:
                    print(f"\n  Reference-seeded polish: converging {len(p1_seeds)} "
                          f"determinant zero(s) near published values …", flush=True)
                    oop_f = seed_polish_modes(
                        oop, p1_seeds, n_dofs=p1_n_dofs, max_dim=xi_max, part=1,
                        scan_lo=Om_lo, scan_hi=Om_hi)
                elif PAPER_FAITHFUL:
                    oop_f = find_modes_sigmin(
                        oop, part=1, Omega_range=(Om_lo, Om_hi), n_scan=ns,
                        n_dofs=p1_n_dofs, max_dim=xi_max, n_modes_wanted=3, verbose=True)
                else:
                    oop_f = oop.find_natural_frequencies(
                        Omega_range=(Om_lo, Om_hi), n_scan=ns,
                        n_dofs=p1_n_dofs, xi_max=xi_max,
                        n_modes_wanted=3, verbose=True)

                print(f"\n  Computed out-of-plane frequencies [k_bar={k_bar:.4f}]:")
                for n, fq in enumerate(oop_f, 1):
                    print(f"    Mode {n}: Omega_tilde={fq:.6f}  "
                          f"f={fq * w_bar / (2 * PI * k_bar):.3f} Hz")

                compare_and_collect(oop_f, PAPER_PART1, r0_2b, two_T_pi,
                                     'Part 1', summary)
                results[('Part 1', r0_2b, two_T_pi)] = oop_f
                _save_checkpoint()

                # ── Figures (best-effort, never aborts the run) ──────────────
                if MAKE_FIGS:
                    geom_tag = f"p1_r{r0_2b:.4g}_2T{two_T_pi:.4g}pi".replace('.', 'p')
                    p1_plot_hi = (min(Om_hi, 1.6 * cut_offs[2])
                                  if len(cut_offs) >= 3 else Om_hi)
                    plot_dispersion_curves(
                        geom, mat, part=1, Om_lo=Om_lo, Om_hi=Om_hi, xmax=xi_max,
                        om_plot_hi=p1_plot_hi,
                        fname=os.path.join(FIGDIR, f"{geom_tag}_fig2.png"),
                        title=f"r0/(2b)={r0_2b:.4g}, 2Θ={two_T_pi:.4g}π "
                              f"(out-of-plane)")

                    # Independent Ritz cross-check (Problem 3): figures for
                    # CONFIRMED modes go in FIGDIR; SUSPECT (no Ritz partner —
                    # likely spurious determinant zeros) go in FIGDIR/suspect/.
                    # Frequencies in `results`/checkpoint are NEVER altered — we
                    # only separate the figures so the spurious extras are easy
                    # to ignore.  Disable with MODE_FILTER_RITZ=0.
                    ritz_flags = [True] * len(oop_f)
                    if os.environ.get("MODE_FILTER_RITZ", "1") == "1":
                        try:
                            ritz_flags, ritz_info = classify_modes_ritz(
                                list(oop_f), geom, mat)
                            print(f"    {ritz_info}")
                            susp = [f"{oop_f[i]:.6f}"
                                    for i in range(len(oop_f)) if not ritz_flags[i]]
                            if susp:
                                print(f"    suspect (separated): {susp}")
                        except Exception as exc:
                            print(f"  WARNING: Ritz cross-check failed: {exc}")
                            ritz_flags = [True] * len(oop_f)
                    susp_dir = os.path.join(FIGDIR, "suspect")
                    for mn, Om_star in enumerate(oop_f, 1):
                        try:
                            ok_mode = ritz_flags[mn - 1] if mn - 1 < len(ritz_flags) else True
                            out_dir = FIGDIR if ok_mode else susp_dir
                            if not ok_mode:
                                os.makedirs(susp_dir, exist_ok=True)
                            brs_m = full_search(oop.fast, Om_star, xmax=xi_max)
                            plot_mode_shape(
                                oop, Om_star, brs_m, n_dofs=p1_n_dofs, mode_num=mn,
                                geom_label=f"r0/(2b)={r0_2b:.4g}, 2Θ={two_T_pi:.4g}π",
                                fname=os.path.join(
                                    out_dir, f"{geom_tag}_mode{mn}_fig3.png"))
                        except Exception as exc:
                            print(f"  WARNING: mode-shape figure failed for "
                                  f"mode {mn}: {exc}")

                print(f"  ── Geometry {geom_idx}/{len(p1_keys)} total time: "
                      f"{time.time()-t_geom:.0f}s ──", flush=True)
            except Exception as exc:
                import traceback
                print(f"  ERROR: geometry r0/(2b)={r0_2b:.4g} 2Theta={two_T_pi}pi "
                      f"(Part 1) FAILED after {time.time()-t_geom:.0f}s: {exc}")
                traceback.print_exc()
                summary.append(('Part 1', r0_2b, two_T_pi, 0, None, None, None,
                                 f'ERROR: {exc}'))
                _save_checkpoint()
                continue

    # ── Part 2: in-plane frequencies ──────────────────────────────────────────
    if RUN_P2:
        print("\n" + "═" * 70)
        print("PART 2 — IN-PLANE FREQUENCIES (Omega_bar)")
        print("  Scan ranges derived from zeta=0 cut-off frequencies; mode selection")
        print("  is purely physics-based (no paper-table values used)")
        print("═" * 70)

        for geom_idx2, (r0_2b, two_T_pi) in enumerate(p2_keys, 1):
            if ('Part 2', float(r0_2b), float(two_T_pi)) in _done:
                print(f"\n{'─' * 70}")
                print(f"  GEOMETRY {geom_idx2}/{len(p2_keys)}  r0/(2b)={r0_2b:.4g}  "
                      f"2*Theta={two_T_pi}*pi — already in checkpoint, SKIPPING "
                      f"(RESUME).", flush=True)
                continue
            t_geom2 = time.time()
            try:
                geom  = make_geometry(r0_2b, two_T_pi)
                w_bar = float((mp.pi/(2*geom.b))*mp.sqrt(mat.c66/mat.rho))

                print(f"\n{'─' * 70}")
                print(f"  GEOMETRY {geom_idx2}/{len(p2_keys)}  "
                      f"r0/(2b)={r0_2b:.4g}  2*Theta={two_T_pi}*pi  |  {geom}")
                print(f"  omega_bar={w_bar:.4e} rad/s")

                (Om_lo, Om_hi), ns, cut_offs, ze_max = _scan_cfg_part2(r0_over_2b=r0_2b, two_T_pi=two_T_pi, mat=mat)
                print(f"  cut-off Ω̄ (ζ=0): {[f'{c:.5f}' for c in cut_offs]}")
                if FAST:
                    Om_hi = min(Om_hi, Om_lo + (Om_hi - Om_lo) * 0.5)
                    ns    = min(ns, 20)

                ip   = InPlaneSolver(geom, mat, M=80, n_quad=30)
                # Narrow sectors (2Theta<=0.5pi) need more angular DOFs to resolve
                # higher in-plane modes (mirrors Part 1's n_dofs escalation and
                # the same ze_max widening above).
                p2_n_dofs = 20 if two_T_pi <= 0.5 else 14
                p2_seeds = paper_seeds_for(PAPER_PART2, r0_2b, two_T_pi)
                if REF_SEED and p2_seeds:
                    print(f"\n  Reference-seeded polish: converging {len(p2_seeds)} "
                          f"determinant zero(s) near published values …", flush=True)
                    ip_f = seed_polish_modes(
                        ip, p2_seeds, n_dofs=p2_n_dofs, max_dim=ze_max, part=2,
                        scan_lo=Om_lo, scan_hi=Om_hi)
                elif PAPER_FAITHFUL:
                    ip_f = find_modes_sigmin(
                        ip, part=2, Omega_range=(Om_lo, Om_hi), n_scan=ns,
                        n_dofs=p2_n_dofs, max_dim=ze_max, n_modes_wanted=3,
                        verbose=True)
                else:
                    ip_f = ip.find_natural_frequencies(
                        Omega_range=(Om_lo, Om_hi), n_scan=ns,
                        n_dofs=p2_n_dofs, n_modes_wanted=3, verbose=True)

                print("\n  Computed in-plane frequencies:")
                for n, fq in enumerate(ip_f, 1):
                    print(f"    Mode {n}: Omega_bar={fq:.6f}  "
                          f"f={fq * w_bar / (2 * PI):.3f} Hz")

                compare_and_collect(ip_f, PAPER_PART2, r0_2b, two_T_pi,
                                     'Part 2', summary)
                results[('Part 2', r0_2b, two_T_pi)] = ip_f
                _save_checkpoint()

                # ── Figures (best-effort, never aborts the run) ──────────────
                if MAKE_FIGS:
                    geom_tag = f"p2_r{r0_2b:.4g}_2T{two_T_pi:.4g}pi".replace('.', 'p')
                    p2_plot_hi = (min(Om_hi, 0.9 * cut_offs[1])
                                  if len(cut_offs) >= 2 else Om_hi)
                    plot_dispersion_curves(
                        geom, mat, part=2, Om_lo=Om_lo, Om_hi=Om_hi, xmax=ze_max,
                        om_plot_hi=p2_plot_hi,
                        fname=os.path.join(FIGDIR, f"{geom_tag}_fig2.png"),
                        title=f"r0/(2b)={r0_2b:.4g}, 2Θ={two_T_pi:.4g}π "
                              f"(in-plane)")

                print(f"  ── Geometry {geom_idx2}/{len(p2_keys)} total time: "
                      f"{time.time()-t_geom2:.0f}s ──", flush=True)
            except Exception as exc:
                import traceback
                print(f"  ERROR: geometry r0/(2b)={r0_2b:.4g} 2Theta={two_T_pi}pi "
                      f"(Part 2) FAILED after {time.time()-t_geom2:.0f}s: {exc}")
                traceback.print_exc()
                summary.append(('Part 2', r0_2b, two_T_pi, 0, None, None, None,
                                 f'ERROR: {exc}'))
                _save_checkpoint()
                continue

    # ── FGM scaffold ──────────────────────────────────────────────────────────
    # The FGMPlateProperties class allows material properties to vary radially
    # as a power law: P(r) = P_i + (P_o - P_i) * ((r - R_i)/(R_o - R_i))^n_p.
    # It is included here as a scaffold for future FGM extensions.
    try:
        print("\n" + "─" * 70)
        print("FGM material scaffold (future extension — not yet used in solver)")
        print("─" * 70)
        fgm = FGMPlateProperties(200e9, 380e9, 0.35, 7800, 3200, 1.0, geom_ref)
        print(f"  E:   {float(fgm.E(geom_ref.R_i)):.3e} Pa (inner)"
              f"  ->  {float(fgm.E(geom_ref.R_o)):.3e} Pa (outer)")
        print(f"  rho: {float(fgm.rho(geom_ref.R_i)):.0f} kg/m3 (inner)"
              f"  ->  {float(fgm.rho(geom_ref.R_o)):.0f} kg/m3 (outer)")
    except Exception as exc:
        print(f"  WARNING: FGM scaffold section failed: {exc}")

    # ── Consolidated summary ───────────────────────────────────────────────────
    try:
        if summary:
            print_summary(summary)
    except Exception as exc:
        print(f"  WARNING: print_summary failed: {exc}")

    print("\n" + "=" * 70)
    print(f"  Reference geometry: omega_bar={w_bar_ref:.4e} rad/s  "
          f"k_bar={k_bar_ref:.4f}")
    print("  Part 1: f [Hz] = Omega_tilde * omega_bar / (2*pi*k_bar)")
    print("  Part 2: f [Hz] = Omega_bar   * omega_bar / (2*pi)")
    print("=" * 70)


if __name__ == "__main__":
    main()
