# e31 = -4.1 COPY (LESSONS Sec 18.246) of probe_piezo_p4_inverse_recovery_2026-09-16.py.
# Written by migrate_e31m.py. Changes, all fixed before running:
#   * E31_TRUE = -4.1 (Liu/Duan Table 1 as printed, Sec 18.244);
#   * e31 inverse bracket [1, 10] -> [-10, -1] (e31bar = e31 - 8.95 is
#     monotone there; the reflection twin sits at 22.0);
#   * relative uncertainties divide by |E31_TRUE|;
#   * coupling-magnitude sanity windows (quantities ~ e31bar^2) scaled
#     by (13.05/4.85)^2 = 7.239;
#   * output files carry the suffix _e31m.
# Anchor gates that pin +4.1 manuscript numbers are reported against
# those numbers and are NOT expected to pass at -4.1.
# -*- coding: utf-8 -*-
"""
probe_piezo_p4_inverse_recovery_2026-09-16.py

Paper 4 (piezoelectric ring), PAPER4_PIEZO_ROADMAP.md Sec 6/Sec 8 step 4:
the parameter-estimation (inverse) demonstration, following Paper 3's
Southwell inverse-b/a methodology (LESSONS_LEARNED.md Sec 18.119) exactly
in kind -- self-consistency round-trip, root-quality gate (with a genuine
negative control), on-branch sensitivity, pre-registered pass/fail -- but
inverting a MATERIAL constant instead of a geometric ratio, using the now-
complete, package-integrated, worker-path-verified C-C piezo-coupled
forward model (plate_solver.piezo_solver.PiezoOutOfPlaneSolver.
cc_coupled_bisect, LESSONS_LEARNED.md Sec 18.171-18.173). C-C is used
because it is the boundary condition with literature (Duan2005 Table 4)
and cluster (Sec 18.148, job 2489140) backing -- the forward model being
inverted here is not a fresh, untrusted piece of physics.

Per the roadmap's own explicit recommendation ("recover Young's modulus
alone with the piezoelectric constants known, then recover a
piezoelectric constant alone with elastic properties known... before
attempting simultaneous multi-parameter recovery"), this probe runs TWO
independent single-parameter recovery demonstrations, not a joint fit:

  TASK 'E'   : recover the host elastic modulus E_steel, all piezo/
               dielectric constants and geometry held at their known
               (true) values.
  TASK 'e31' : recover the piezoelectric stress constant e31, all
               elastic constants and geometry held at their known (true)
               values.

Both use the SAME synthetic target mode: C-C, p=0 (axisymmetric),
h1/2h=1/12, fundamental radial mode (n_radial=0) -- exactly the point
already validated to <1e-13 relative against the standalone forward-model
probe script in Sec 18.173's TestPiezoCoupledCC gate, so the forward
model itself is not what this probe is testing.

WHY BOTH PARAMETERS ARE WORTH RUNNING TOGETHER, NOT JUST E: a sandbox
sweep this session (dps=100) already shows the two parameters are
qualitatively very different in how well this single mode constrains
them --

    domega/dE  (E in [140e9,260e9], e31 fixed at 4.1): well-conditioned,
        an 30% swing in E moves the root by ~11% (2427->3114 rad/s over
        E=140e9->260e9) -- ordinary bending-stiffness scaling, omega ~
        sqrt(E).
    domega/de31 (e31 in [1,8], E fixed at 200e9): the root moves by only
        ~0.06 rad/s (2791.82->2791.76) across the ENTIRE range -- i.e. a
        relative sensitivity around 4 orders of magnitude weaker than
        E's, consistent with LESSONS_LEARNED Sec 18.141 F's own finding
        (Duan2005 Fig. 3) that under short-circuit conditions the
        piezoelectric layer's frequency effect is dominated by ordinary
        added-layer elastic stiffness, with true electromechanical
        coupling contributing only a ~0.01-0.016% correction.

This is exactly the kind of planning input Sec 4 item 4's resolution
flagged for this demonstration, not a blocker discovered late -- so
PRE-REGISTERING it here, before running, is what "never accept a result
without knowing in advance what it would mean" requires:

PRE-REGISTERED EXPECTATION (decide before reading results):
  - E: expect GOOD conditioning -- a plausible frequency-measurement
    precision (0.01%-0.1%) should map to a comparably small (same order
    of magnitude, not orders larger) relative uncertainty in the
    recovered E.
  - e31: expect POOR conditioning -- the same 0.01%-0.1% frequency
    precision should map to a relative e31 uncertainty that is LARGE
    (order 1 or bigger, i.e. e31 is effectively unidentifiable from this
    one short-circuit mode alone). This is the honest, physically
    expected finding, not a probe-design failure -- do not tune anything
    to make it look better. If the run instead shows e31 comparably
    well-conditioned to E, THAT would be the surprising result requiring
    investigation (either a probe-design bug or a genuine new finding
    about this configuration), not e31's poor conditioning itself.

PRE-REGISTERED PASS/FAIL (both required for SENTINEL PASS_ALL; separate
from the sensitivity numbers above, which are reported either way, not
gated):
  G0 (self-consistency): recover the true parameter from its own EXACT
     forward-computed target frequency (no injected measurement noise)
     to <1e-6 relative -- matching Southwell's own G0 bar exactly
     (LESSONS_LEARNED Sec 18.119). This is a pure numerical round-trip
     test (arbitrary-precision arithmetic throughout, via mpmath dps=100
     inside PiezoOutOfPlaneSolver) and is expected to pass tightly for
     BOTH parameters regardless of physical conditioning -- poor
     conditioning affects how *measurement* error propagates, not
     whether the exact inversion round-trips.
  G1 (root quality): a genuine sign flip of the shooting residual
     G(p) = omega_root(p) - omega_target across the search bracket, PLUS
     a negative control -- G(p) evaluated at both ends of a bracket
     entirely on ONE side of the true root must show NO sign flip
     (confirms the sign-flip check itself is not trivially satisfied).

G2 (on-branch sensitivity, reported, not gated): forward finite-difference
  domega/dp at the true parameter (central difference, +-1% of p_true),
  cross-checked against an on-branch inverse finite-difference dp/domega
  (invert at omega_target +- a delta sized from the forward derivative so
  the perturbed target stays safely within the achievable omega range for
  BOTH parameters -- e31's achievable omega range from this mode is only
  ~0.06 rad/s wide, so a naive fixed absolute or relative omega
  perturbation picked for E would very likely fall outside e31's bracket
  entirely and produce a bogus/unbracketed inverse solve; sizing the
  perturbation from the forward derivative avoids that failure mode for
  either parameter). Cross-check ratio (domega/dp)*(dp/domega) should be
  approx 1.0 (implicit-function-theorem self-consistency). Then reports
  the propagated relative parameter uncertainty for two representative
  assumed frequency-measurement precisions (0.01%, 0.1% relative) --
  the actual scientific payoff of this demonstration.

Cost model (sandbox-timed this session, dps=100): one cc_coupled_bisect
call at this mode/geometry costs ~0.11-0.12s per bisection iteration
(iters=30 -> ~3.5s). This probe's default INNER_ITERS=40 and
OUTER_ITERS=26 (each far more than needed for <1e-6 relative on the
*parameter* -- bisection halves the bracket every iteration regardless of
how flat the underlying function is, so poor conditioning does not by
itself require more outer iterations, only enough INNER precision that
sign determination near the true root is not swamped by inner-solve
discretization; dps=100/INNER_ITERS=40 gives inner omega-resolution
~1e-12 absolute, five-plus orders of magnitude below the smallest signal
this probe ever needs to resolve). Per parameter: 1 forward target +
1 full G0 inversion (~26+2=28 inner bisections) + 4 G1 scan points
(forward only) + 2 G2 forward-FD evals + 2 more full inversions (G2
on-branch cross-check) = 3 full inversions + 8 forward evals per
parameter, roughly 3*28 + 8 = ~92 inner-bisect-equivalents *4.6s
(iters=40 scaling) =~ 425s =~ 7 minutes per parameter. Two parameters run
concurrently (ProcessPoolExecutor, max_workers=2) -> wall time budget
should not exceed ~15 minutes; the 30-minute SLURM request is generous
headroom, not a tight estimate, matching this project's own stated
practice for small jobs.

Standalone-vs-package note: unlike most of this project's cluster probes
(which re-derive the math standalone in pure mpmath to avoid any
PKG_PATH/SOLVER_VERSION coupling), THIS probe deliberately imports
`plate_solver.piezo_solver.PiezoOutOfPlaneSolver` directly -- the whole
point of Sec 18.171-18.173's package-integration work was to have a
single validated, reusable forward model rather than re-deriving it a
fourth time. PiezoOutOfPlaneSolver is not wired into the SOLVER_VERSION
checkpoint-compatibility system at all (it takes plain scalars, no
PlateGeometry/MaterialModel), so there is no checkpoint-resume risk from
importing it -- but the usual SOLVER_VERSION preflight print is kept
anyway, purely as a deployed-code sanity check (confirms this job is
running against the code this session actually tested, not a stale copy).
"""
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

sys.path.insert(0, os.environ.get("PKG_PATH", "."))
import plate_solver as ps  # noqa: E402
from plate_solver.piezo_solver import PiezoOutOfPlaneSolver  # noqa: E402

EXPECT_SOLVER_VERSION = os.environ.get("EXPECT_SOLVER_VERSION", "2026-07-10.s10")

DPS = int(os.environ.get("PIEZO_DPS", "100"))
INNER_ITERS = int(os.environ.get("PIEZO_INNER_ITERS", "40"))
OUTER_ITERS = int(os.environ.get("PIEZO_OUTER_ITERS", "26"))
G0_TOL_REL = 1e-6

# Duan2005 Table 1 (PZT4 piezo layer + steel host) -- same geometry/
# material this project's own cluster gate (Sec 18.148, job 2489140) and
# Sec 18.173's package tests already validate the forward model against.
R_I, R_O, H = 0.1, 0.6, 0.01
E_TRUE, NU, RHO = 200e9, 0.3, 7800.0
C11E, C12E, C13E, C33E = 132e9, 71e9, 73e9, 115e9
E31_TRUE, E33 = -4.1, 14.1
X11, X33 = 7.124e-9, 5.841e-9
RHO_PZT = 7500.0
H1 = (2.0 / 12) * H  # Duan2005 Table 4's own h1/2h=1/12 ratio

BASE_KWARGS = dict(
    r_i=R_I, r_o=R_O, h=H, nu=NU, rho=RHO, h1=H1,
    C11E=C11E, C12E=C12E, C13E=C13E, C33E=C33E,
    e33=E33, X11=X11, X33=X33, rho_pzt=RHO_PZT, dps=DPS,
)

# Known-good target-mode anchor (Sec 18.173) -- used only to seed the
# omega-bracket guesses below, not as a substitute for each trial's own
# genuine bisection.
OMEGA_ANCHOR = 2791.7866198564


def _forward_omega_E(E, iters, lo=None, hi=None):
    """omega_root(E), e31 fixed at E31_TRUE. Bracket auto-scales with
    sqrt(E) (the ordinary plate-bending scaling already confirmed by this
    session's own sandbox sweep) so this stays valid across a wide E
    sweep, not just near E_TRUE."""
    if lo is None or hi is None:
        center = OMEGA_ANCHOR * (E / E_TRUE) ** 0.5
        lo, hi = 0.80 * center, 1.20 * center
    kwargs = dict(BASE_KWARGS); kwargs["E"] = E; kwargs["e31"] = E31_TRUE
    solver = PiezoOutOfPlaneSolver(**kwargs)
    return solver.cc_coupled_bisect(lo, hi, 0, iters=iters)


def _forward_omega_e31(e31, iters, lo=2780.0, hi=2800.0):
    """omega_root(e31), E fixed at E_TRUE. Fixed absolute bracket is
    deliberate and safe: this session's sandbox sweep already confirmed
    the root only moves ~0.06 rad/s across e31 in [1,8], so a wide fixed
    bracket comfortably contains the root for any e31 this probe tries."""
    kwargs = dict(BASE_KWARGS); kwargs["E"] = E_TRUE; kwargs["e31"] = e31
    solver = PiezoOutOfPlaneSolver(**kwargs)
    return solver.cc_coupled_bisect(lo, hi, 0, iters=iters)


FORWARD_FNS = {"E": _forward_omega_E, "e31": _forward_omega_e31}
TRUE_VALUES = {"E": E_TRUE, "e31": E31_TRUE}
# Outer (parameter-space) search brackets for the inverse solve.
PARAM_BRACKETS = {"E": (140e9, 260e9), "e31": (-10.0, -1.0)}


def _invert_param(param_name, omega_target, lo_p, hi_p, iters_inner, iters_outer):
    """Bisect the shooting residual G(p) = omega_root(p) - omega_target on
    [lo_p, hi_p]. Returns (p_hat, sign_flip_seen, g_lo, g_hi)."""
    fwd = FORWARD_FNS[param_name]
    g_lo = fwd(lo_p, iters_inner) - omega_target
    g_hi = fwd(hi_p, iters_inner) - omega_target
    sign_flip = (g_lo > 0) != (g_hi > 0)
    lo, hi, glo = lo_p, hi_p, g_lo
    for _ in range(iters_outer):
        mid = 0.5 * (lo + hi)
        gm = fwd(mid, iters_inner) - omega_target
        if (gm > 0) == (glo > 0):
            lo, glo = mid, gm
        else:
            hi = mid
    return 0.5 * (lo + hi), sign_flip, g_lo, g_hi


def _worker(param_name):
    t_start = time.time()
    p_true = TRUE_VALUES[param_name]
    lo_p, hi_p = PARAM_BRACKETS[param_name]
    log = []

    def note(msg):
        log.append(msg)
        print(f"[{param_name}] {msg}", flush=True)

    # --- forward target (exact synthetic "measurement") ---
    t0 = time.time()
    omega_target = FORWARD_FNS[param_name](p_true, INNER_ITERS)
    note(f"forward target: omega({param_name}={p_true:.6g}) = "
         f"{omega_target:.10f} rad/s ({time.time()-t0:.1f}s)")

    # --- G0: self-consistency round trip ---
    t0 = time.time()
    p_hat, sign_flip_true, g_lo, g_hi = _invert_param(
        param_name, omega_target, lo_p, hi_p, INNER_ITERS, OUTER_ITERS)
    g0_rel_err = abs(p_hat - p_true) / abs(p_true)
    g0_pass = sign_flip_true and (g0_rel_err < G0_TOL_REL)
    note(f"G0 self-consistency: p_hat={p_hat:.10g} p_true={p_true:.6g} "
         f"rel_err={g0_rel_err:.3e} sign_flip={sign_flip_true} "
         f"pass={g0_pass} ({time.time()-t0:.1f}s)")

    # --- G1: root quality + negative control ---
    # Positive control already have (sign_flip_true above, full bracket).
    # Negative control: a bracket entirely on the hi-side of the true
    # root -- must show NO sign flip.
    mid_p = 0.5 * (p_hat + hi_p)
    t0 = time.time()
    g_mid = FORWARD_FNS[param_name](mid_p, INNER_ITERS) - omega_target
    g_hi_side = FORWARD_FNS[param_name](hi_p, INNER_ITERS) - omega_target
    neg_control_no_flip = (g_mid > 0) == (g_hi_side > 0)
    g1_pass = sign_flip_true and neg_control_no_flip
    note(f"G1 root quality: positive_bracket_sign_flip={sign_flip_true} "
         f"negative_control(hi-side only)_no_flip={neg_control_no_flip} "
         f"pass={g1_pass} ({time.time()-t0:.1f}s)")

    # --- G2: on-branch sensitivity (reported, not gated) ---
    t0 = time.time()
    dp_fd = 0.01 * p_true  # +-1% central difference
    om_plus = FORWARD_FNS[param_name](p_true + dp_fd, INNER_ITERS)
    om_minus = FORWARD_FNS[param_name](p_true - dp_fd, INNER_ITERS)
    domega_dp = (om_plus - om_minus) / (2 * dp_fd)
    note(f"G2 forward FD: domega/d{param_name} = {domega_dp:.6e} "
         f"rad/s per unit ({time.time()-t0:.1f}s)")

    t0 = time.time()
    # Size the on-branch omega perturbation from the forward derivative so
    # it stays safely inside the achievable range for either parameter's
    # own bracket (see module docstring).
    d_omega = abs(domega_dp) * (0.01 * p_true) * 0.5
    if d_omega == 0.0:
        d_omega = 1e-6 * abs(omega_target)
    p_plus, sf_plus, _, _ = _invert_param(
        param_name, omega_target + d_omega, lo_p, hi_p, INNER_ITERS, OUTER_ITERS)
    p_minus, sf_minus, _, _ = _invert_param(
        param_name, omega_target - d_omega, lo_p, hi_p, INNER_ITERS, OUTER_ITERS)
    dp_domega_onbranch = (p_plus - p_minus) / (2 * d_omega)
    cross_check_ratio = domega_dp * dp_domega_onbranch
    note(f"G2 on-branch inverse FD: dp/domega = {dp_domega_onbranch:.6e}, "
         f"cross_check (domega/dp)*(dp/domega) = {cross_check_ratio:.6f} "
         f"(expect ~1.0) ({time.time()-t0:.1f}s)")

    propagated = {}
    for rel_prec in (1e-4, 1e-3):  # 0.01%, 0.1% assumed frequency precision
        d_om_meas = rel_prec * omega_target
        dp_prop = d_om_meas / abs(domega_dp)
        propagated[f"{rel_prec:g}"] = {
            "assumed_omega_rel_precision": rel_prec,
            "propagated_dp": dp_prop,
            "propagated_dp_rel_to_p_true": dp_prop / abs(p_true),
        }
        note(f"G2 propagated uncertainty: {rel_prec:.4%} freq precision -> "
             f"d{param_name}={dp_prop:.6g} "
             f"({dp_prop/abs(p_true):.4%} of {param_name}_true)")

    elapsed = time.time() - t_start
    return {
        "param": param_name,
        "p_true": p_true,
        "omega_target": omega_target,
        "g0_p_hat": p_hat,
        "g0_rel_err": g0_rel_err,
        "g0_pass": g0_pass,
        "g1_positive_sign_flip": sign_flip_true,
        "g1_negative_control_no_flip": neg_control_no_flip,
        "g1_pass": g1_pass,
        "g2_domega_dp_forward": domega_dp,
        "g2_dp_domega_onbranch": dp_domega_onbranch,
        "g2_cross_check_ratio": cross_check_ratio,
        "g2_propagated": propagated,
        "elapsed_s": elapsed,
        "log": log,
    }


def main():
    print(f"job start: dps={DPS} inner_iters={INNER_ITERS} "
          f"outer_iters={OUTER_ITERS}", flush=True)
    if ps.SOLVER_VERSION != EXPECT_SOLVER_VERSION:
        print(f"PREFLIGHT FAIL: SOLVER_VERSION={ps.SOLVER_VERSION!r} != "
              f"{EXPECT_SOLVER_VERSION!r}", flush=True)
        sys.exit(2)
    print(f"preflight OK: SOLVER_VERSION={ps.SOLVER_VERSION}", flush=True)

    tasks = ["E", "e31"]
    results = {}
    with ProcessPoolExecutor(max_workers=int(os.environ.get("N_WORKERS", "2"))) as ex:
        futs = {ex.submit(_worker, t): t for t in tasks}
        for fut in as_completed(futs):
            r = fut.result()
            results[r["param"]] = r
            print(f"=== {r['param']} DONE in {r['elapsed_s']:.1f}s "
                  f"G0_pass={r['g0_pass']} G1_pass={r['g1_pass']} ===",
                  flush=True)

    all_pass = all(results[t]["g0_pass"] and results[t]["g1_pass"] for t in tasks)

    with open("piezo_p4_inverse_recovery_results_e31m.json", "w") as f:
        json.dump(results, f, indent=2)

    print("--- summary ---", flush=True)
    for t in tasks:
        r = results[t]
        print(f"  {t}: G0_rel_err={r['g0_rel_err']:.3e} "
              f"G1_pass={r['g1_pass']} "
              f"cross_check_ratio={r['g2_cross_check_ratio']:.4f} "
              f"propagated_rel_unc@0.01%={r['g2_propagated']['0.0001']['propagated_dp_rel_to_p_true']:.4%} "
              f"propagated_rel_unc@0.1%={r['g2_propagated']['0.001']['propagated_dp_rel_to_p_true']:.4%}",
              flush=True)

    print("SENTINEL " + ("PASS_ALL" if all_pass else "FAIL_ALL") +
          f" g0_g1_pass={all_pass}", flush=True)


if __name__ == "__main__":
    main()
