# -*- coding: utf-8 -*-
"""
probe_piezo_p4_lever16_freq_noise_floor_2026-09-16.py

Paper 4 lever 16 (PAPER4_LEVERS_2026-09-16.md, section A.16): non-0.01%
noise. G2 (propagated relative e31 uncertainty) is already tabulated at
0.01% and 0.1% frequency precision for F-F n=0 OC (Sec 18.175) and C-F
n=0 OC (Sec 18.192). This probe is NOT a new derivation and NOT a
geometry port (lever 11 is Liu's solid disk -- not started here). It
reuses `oc_ff_bisect` / `oc_cf_bisect` exactly as already validated and:

  (a) sweeps the SAME linear-propagation formula already used by
      probe_piezo_p4_ff_open_circuit_2026-09-16.py / cf_mixed_edge
      across a grid of assumed relative frequency errors, instead of
      the two hardcoded columns (0.01%, 0.1%), so the curve -- and its
      linearity in delta-omega/omega -- is visible directly, and
  (b) adds one lab-realistic floor sourced from the Boeringa-McCune PZT
      rod validation work (Project Knowledge/papers/BoeringaMcCune2026.pdf,
      the polished JSV-style paper -- NOT the "-2" draft, which contains
      internally inconsistent Table 4/5 numbers and explicit
      placeholder text ("[Don't include this]", "ask Stutts about
      citing exact instrumentation") and is not fit to source a number
      from).

LAB FLOOR SOURCE (read directly from the PDF via pdftotext, not
invented): BoeringaMcCune2026.pdf Table 3 reports, for the first five
axial resonances of a free-free PZT rod, both an HP4294A impedance-
analyzer frequency and an SLDV-measured frequency for the same physical
mode. The paper does not state a numeric Hz/ppm resolution spec for
either instrument (the "high-resolution sweeps ... near the expected
... resonances" sentence is qualitative), so the closest available
lab-realistic number is the two instruments' own disagreement on the
same physical resonance, which folds in real specimen/mounting/
temperature effects an electronic spec would not:

    mode   impedance (kHz)   SLDV (kHz)   |rel diff|
    1      40.516            40.563       0.116%   (paper's own printed "Imp. ref." vs SLDV rows: 0.12%)
    2      74.367            74.441       0.099%   (0.10%)
    3      112.160           112.262      0.091%   (0.09%)
    4      144.550           144.705      0.107%   (0.11%)
    5      183.750           183.917      0.091%   (0.09%)

    mean of the paper's own printed percentages (0.12, 0.10, 0.09,
    0.11, 0.09) = 0.102%.

This is used as REL_PREC_LAB = 1.02e-3 (0.102%), labelled throughout as
"Boeringa-McCune impedance-vs-SLDV floor", explicitly NOT a Boeringa-
McCune number for a piezoelectric ring or for e31 -- it is a same-
material-system (PZT), same-lab, cross-instrument frequency
disagreement used as a stand-in for what "0.01%" means physically. If
the user wants a different or tighter floor, this section says exactly
which number was used and why, so it can be swapped without re-deriving
anything else.

METHOD (linear propagation, unchanged formula from the existing
probes -- do NOT retune):
    domega/de31 computed by the same central finite difference already
    used (dp = 1% of e31_true, ITERS=40 bisection iterations, same
    brackets), at the same Duan geometry h1/2h=1/12.
    For an assumed relative frequency error rel_prec:
        d_om      = rel_prec * omega_target
        dp_prop   = d_om / abs(domega/de31)
        rel_unc   = dp_prop / abs(e31_true)
    This is linear in rel_prec by construction; the point of sweeping
    a grid (not just 0.01% and 0.1%) is to show that linearity
    explicitly on a log-log curve, not to hide it.

PRE-REGISTERED GATES:
    G0 (round-trip, both F-F and C-F): recovering e31_true from its own
       exact OC frequency (the existing bisection inverse, same
       brackets [1,10], OUTER_ITERS=26) must still give <1e-6 relative
       error at THIS geometry -- this is the existing G4 gate,
       re-verified here as a precondition, not re-derived.
    G_anchor (both F-F and C-F): the propagated uncertainty at
       rel_prec=1e-4 (0.01%) and 1e-3 (0.1%) must reproduce the
       already-published Table~tab:inv / LESSONS_LEARNED numbers
       (F-F: 1.98% / 19.8%; C-F: 2.12% / 21.2%) to within 1% relative
       -- this is a consistency check against the existing record, not
       a new finding.
    G2 (reported, not gated): the full curve, F-F and C-F, plus the
       Boeringa-McCune lab-floor row. Reported linearity check: for
       every pair of grid points, rel_unc(p2)/rel_unc(p1) must equal
       p2/p1 to a tight numerical tolerance (linear-propagation
       formula, so this is an implementation self-check, not physics).

SOLVER_VERSION is not bumped (no package change). dp (the e31 finite-
difference step) is not retuned. Brackets, ITERS, OUTER_ITERS, and
BASE_KWARGS are byte-for-byte the existing F-F/C-F probes' values.
"""
import json
import os
import sys
import time

sys.path.insert(0, os.environ.get("PKG_PATH",
                                  os.path.join(os.path.dirname(__file__), "..")))
import plate_solver as ps  # noqa: E402
from plate_solver.piezo_solver import PiezoOutOfPlaneSolver  # noqa: E402

EXPECT_SOLVER_VERSION = os.environ.get("EXPECT_SOLVER_VERSION", "2026-07-10.s10")
DPS = int(os.environ.get("PIEZO_DPS", "60"))
ITERS = int(os.environ.get("PIEZO_ITERS", "40"))
OUTER_ITERS = int(os.environ.get("PIEZO_OUTER_ITERS", "26"))
G0_TOL_REL = 1e-6

# --- unchanged from probe_piezo_p4_ff_open_circuit_2026-09-16.py /
#     probe_piezo_p4_cf_mixed_edge_2026-09-16.py ---
R_I, R_O, H = 0.1, 0.6, 0.01
E_TRUE, NU, RHO = 200e9, 0.3, 7800.0
C11E, C12E, C13E, C33E = 132e9, 71e9, 73e9, 115e9
E31_TRUE, E33 = 4.1, 14.1
X11, X33 = 7.124e-9, 5.841e-9
RHO_PZT = 7500.0
H1 = (2.0 / 12) * H

BASE_KWARGS = dict(
    r_i=R_I, r_o=R_O, h=H, E=E_TRUE, nu=NU, rho=RHO, h1=H1,
    C11E=C11E, C12E=C12E, C13E=C13E, C33E=C33E,
    e33=E33, X11=X11, X33=X33, rho_pzt=RHO_PZT, dps=DPS,
)

FF_LO, FF_HI = 720.0, 800.0
CF_LO, CF_HI = 400.0, 450.0
E31_LO, E31_HI = 1.0, 10.0

# Grid of assumed relative frequency errors (dimensionless, e.g. 1e-4 = 0.01%).
# Log-spaced from 1e-5 to 1e-2 so the linear scaling is visible over three
# decades; includes the two already-published anchor points (1e-4, 1e-3)
# exactly, plus the Boeringa-McCune lab floor.
REL_PREC_LAB = 1.02e-3  # Boeringa-McCune impedance-vs-SLDV mean, see header
REL_PREC_GRID = sorted(set([
    1e-5, 3e-5, 1e-4, 3e-4, 1e-3, REL_PREC_LAB, 3e-3, 1e-2,
]))

PUBLISHED_ANCHORS = {
    # (rel_prec) -> {"ff": pct, "cf": pct}, from PAPER4_PIEZO_DRAFT.tex
    # Table tab:inv and LESSONS_LEARNED.md Sec 18.175 / Sec 18.192.
    # 2026-09-23 Sec 18.231: re-published values under the consistent OC
    # charge arm (duan: 1e-4 ff 1.98 / cf 2.12; 1e-3 ff 19.8 / cf 21.2).
    1e-4: {"ff": 2.13, "cf": 2.28},
    1e-3: {"ff": 21.3, "cf": 22.8},
}


def _solver(e31):
    kw = dict(BASE_KWARGS)
    kw["e31"] = e31
    return PiezoOutOfPlaneSolver(**kw)


def _omega_oc_ff(e31, iters=ITERS):
    return _solver(e31).oc_ff_bisect(FF_LO, FF_HI, 0, iters=iters)


def _omega_oc_cf(e31, iters=ITERS):
    return _solver(e31).oc_cf_bisect(CF_LO, CF_HI, 0, iters=iters)


def _invert_e31(omega_fn, omega_target, lo_p, hi_p, iters_inner, iters_outer):
    g_lo = omega_fn(lo_p, iters_inner) - omega_target
    g_hi = omega_fn(hi_p, iters_inner) - omega_target
    sign_flip = (g_lo > 0) != (g_hi > 0)
    lo, hi, glo = lo_p, hi_p, g_lo
    for _ in range(iters_outer):
        mid = 0.5 * (lo + hi)
        gm = omega_fn(mid, iters_inner) - omega_target
        if (gm > 0) == (glo > 0):
            lo, glo = mid, gm
        else:
            hi = mid
    return 0.5 * (lo + hi), sign_flip


def _channel(name, omega_fn, lo_bracket, hi_bracket):
    """Run the full lever-16 pipeline for one channel (F-F or C-F)."""
    log = []
    def note(msg):
        log.append(msg)
        print("[%s] %s" % (name, msg), flush=True)

    t0 = time.time()
    omega_target = omega_fn(E31_TRUE, ITERS)
    note("omega_target(e31_true)=%.10f (%.1fs)" % (omega_target, time.time() - t0))

    # G0: round-trip inverse at the exact frequency
    t0 = time.time()
    p_hat, sign_flip = _invert_e31(
        omega_fn, omega_target, E31_LO, E31_HI, ITERS, OUTER_ITERS)
    g0_rel = abs(p_hat - E31_TRUE) / abs(E31_TRUE)
    g0_pass = sign_flip and (g0_rel < G0_TOL_REL)
    note("G0 round-trip: p_hat=%.10g rel_err=%.3e sign_flip=%s pass=%s (%.1fs)"
         % (p_hat, g0_rel, sign_flip, g0_pass, time.time() - t0))

    # forward FD sensitivity, same recipe as the existing probes
    t0 = time.time()
    dp = 0.01 * E31_TRUE
    om_plus = omega_fn(E31_TRUE + dp, ITERS)
    om_minus = omega_fn(E31_TRUE - dp, ITERS)
    domega_dp = (om_plus - om_minus) / (2 * dp)
    note("forward FD: domega/de31=%.6e rad/s per C/m^2 (%.1fs)"
         % (domega_dp, time.time() - t0))

    curve = []
    for rel_prec in REL_PREC_GRID:
        d_om = rel_prec * omega_target
        dp_prop = d_om / abs(domega_dp)
        rel_u = dp_prop / abs(E31_TRUE)
        curve.append({
            "rel_prec": rel_prec,
            "rel_prec_pct": rel_prec * 100.0,
            "propagated_de31": dp_prop,
            "rel_unc_pct": rel_u * 100.0,
            "is_lab_floor": (rel_prec == REL_PREC_LAB),
        })
        note("rel_prec=%.5f%% -> e31 unc=%.4f%%%s"
             % (rel_prec * 100, rel_u * 100,
                "  <-- Boeringa-McCune impedance/SLDV floor" if rel_prec == REL_PREC_LAB else ""))

    # G_anchor: compare against the published 0.01%/0.1% numbers
    g_anchor_pass = True
    for rp, pub in PUBLISHED_ANCHORS.items():
        row = next(r for r in curve if abs(r["rel_prec"] - rp) < 1e-12)
        pub_pct = pub[name]
        rel_diff = abs(row["rel_unc_pct"] - pub_pct) / pub_pct
        ok = rel_diff < 0.01
        g_anchor_pass = g_anchor_pass and ok
        note("G_anchor rel_prec=%.4f%%: computed=%.4f%% published=%.4f%% "
             "rel_diff=%.4e pass=%s" % (rp * 100, row["rel_unc_pct"], pub_pct,
                                         rel_diff, ok))

    # linearity self-check (reported, formula is linear by construction)
    lin_ok = True
    base = curve[0]
    for row in curve[1:]:
        expected = base["rel_unc_pct"] * (row["rel_prec"] / base["rel_prec"])
        rel_diff = abs(row["rel_unc_pct"] - expected) / expected
        if rel_diff > 1e-9:
            lin_ok = False
    note("linearity self-check pass=%s" % lin_ok)

    return {
        "omega_target": omega_target,
        "g0_p_hat": p_hat,
        "g0_rel_err": g0_rel,
        "g0_pass": g0_pass,
        "domega_de31": domega_dp,
        "curve": curve,
        "g_anchor_pass": g_anchor_pass,
        "linearity_pass": lin_ok,
        "log": log,
    }


def main():
    t_all = time.time()
    print("job start: dps=%s iters=%s outer_iters=%s" % (DPS, ITERS, OUTER_ITERS),
          flush=True)
    if ps.SOLVER_VERSION != EXPECT_SOLVER_VERSION:
        print("PREFLIGHT FAIL: SOLVER_VERSION=%r != %r"
              % (ps.SOLVER_VERSION, EXPECT_SOLVER_VERSION), flush=True)
        sys.exit(2)
    print("preflight OK: SOLVER_VERSION=%s" % ps.SOLVER_VERSION, flush=True)
    if not (hasattr(PiezoOutOfPlaneSolver, "oc_ff_bisect")
            and hasattr(PiezoOutOfPlaneSolver, "oc_cf_bisect")):
        print("PREFLIGHT FAIL: deployed PiezoOutOfPlaneSolver missing "
              "oc_ff_bisect/oc_cf_bisect.", flush=True)
        sys.exit(2)
    print("preflight OK: oc_ff_bisect + oc_cf_bisect present", flush=True)

    ff = _channel("ff", _omega_oc_ff, FF_LO, FF_HI)
    cf = _channel("cf", _omega_oc_cf, CF_LO, CF_HI)

    all_pass = (ff["g0_pass"] and ff["g_anchor_pass"] and ff["linearity_pass"]
                and cf["g0_pass"] and cf["g_anchor_pass"] and cf["linearity_pass"])

    lab_row_ff = next(r for r in ff["curve"] if r["is_lab_floor"])
    lab_row_cf = next(r for r in cf["curve"] if r["is_lab_floor"])

    results = {
        "rel_prec_grid": REL_PREC_GRID,
        "rel_prec_lab": REL_PREC_LAB,
        "lab_floor_source": (
            "BoeringaMcCune2026.pdf Table 3, mean of printed impedance-vs-"
            "SLDV frequency-discrepancy percentages across 5 axial modes "
            "of a free-free PZT rod (0.12, 0.10, 0.09, 0.11, 0.09 -> "
            "mean 0.102%). Not a Boeringa-McCune ring or e31 number; used "
            "as a same-lab, same-material-system cross-instrument "
            "frequency floor. The '-2' draft PDF was checked and rejected "
            "as a source (internally inconsistent Table 4/5, explicit "
            "placeholder text)."
        ),
        "ff": ff,
        "cf": cf,
        "lab_floor_row_ff_pct": lab_row_ff["rel_unc_pct"],
        "lab_floor_row_cf_pct": lab_row_cf["rel_unc_pct"],
        "all_pass": all_pass,
        "elapsed_s": time.time() - t_all,
    }
    out_path = os.path.join(os.path.dirname(__file__) or ".",
                            "piezo_p4_lever16_freq_noise_floor_results.json")
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)

    print("--- summary ---", flush=True)
    print("  F-F: G0=%s G_anchor=%s linearity=%s lab_floor(%.4f%%)=%.4f%%"
          % (ff["g0_pass"], ff["g_anchor_pass"], ff["linearity_pass"],
             REL_PREC_LAB * 100, lab_row_ff["rel_unc_pct"]), flush=True)
    print("  C-F: G0=%s G_anchor=%s linearity=%s lab_floor(%.4f%%)=%.4f%%"
          % (cf["g0_pass"], cf["g_anchor_pass"], cf["linearity_pass"],
             REL_PREC_LAB * 100, lab_row_cf["rel_unc_pct"]), flush=True)
    print("SENTINEL " + ("PASS_ALL" if all_pass else "FAIL_ALL")
          + " all_pass=%s" % all_pass, flush=True)


if __name__ == "__main__":
    main()
