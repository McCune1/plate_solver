# e31 = -4.1 COPY (LESSONS Sec 18.246) of probe_piezo_p4_sequential_e_e31_2026-09-16.py.
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
probe_piezo_p4_sequential_e_e31_2026-09-16.py

Paper 4 Lever 15 (PAPER4_LEVERS_2026-09-16.md sec D.15): the inverse
chapter's own prescribed experiment, run for the first time as actually
prescribed -- SEQUENTIAL, not two independent single-parameter demos.
Roadmap Sec 6's recipe and job 2491531 (LESSONS_LEARNED.md Sec 18.174)
each recovered E and e31 INDEPENDENTLY from the same C-C short-circuit
mode, both holding the OTHER parameter at its true value. That is not
what a real experiment can do: e31 is not known in advance, so any real
recovery of e31 from a second observable must use the E already
recovered from the first observable, carrying whatever uncertainty that
estimate has -- not the untouchable true E.

THIS PROBE:
  Step 1: recover host E from C-C short-circuit flexure (job 2491531's
          own cc_coupled_bisect operator, unchanged, e31 held at its
          true value E31_TRUE -- reusing, not retuning, that job's
          result). Call the recovered value E_HAT.
  Step 2: HOLD E_HAT (not E_TRUE) and recover e31 from F-F n=0
          open-circuit (oc_ff_bisect, Sec 18.175's operator, unchanged).
          The "measurement" omega_OC used for this inversion is
          nature's TRUE open-circuit frequency, omega_OC(E_TRUE,
          e31_TRUE) -- what a real experiment on a real ring would
          report -- but the solver used to invert it is built with
          E=E_HAT, exactly as a real analyst who does not know E_TRUE
          would have to do.
  Step 3: propagate BOTH uncertainty sources into the final e31
          estimate: (a) the DIRECT sensitivity of omega_OC to e31
          itself (Sec 18.175's already-known number, ~1.98% of
          e31_TRUE at 0.01% frequency precision); and (b) the INDIRECT
          sensitivity that comes from E_HAT itself carrying step 1's
          own measurement uncertainty (implicit-function-theorem
          chain: de31/dE |omega_OC fixed =
          -(domega_OC/dE)/(domega_OC/de31), then
          delta_e31_from_E = |de31/dE| * delta_E_hat). Report both
          contributions and their quadrature-combined total -- this
          compound number, not the single-observable 1.98%, is the
          paper-honest uncertainty on a sequentially-recovered e31.
  Step 4 (negative control, a gate, not just a note): repeat the exact
          same sequential recipe but with the e31 OBSERVABLE swapped
          from F-F n=0 OC to the C-C short-circuit mode itself (i.e.
          the recipe job 2491531 already ran, this time holding E_HAT
          instead of E_TRUE) -- this MUST still fail at order-100%
          e31 uncertainty, exactly like job 2491531's 731% number.
          C-C fully-electroded OC coincides with C-C SC at leading
          order (Sec 18.175 theorem A), so the SC 6x6
          (cc_coupled_bisect) IS the C-C OC observable at this order;
          no separate C-C OC solver exists or is needed
          (PAPER4_LEVERS_2026-09-16.md standing ban: do not invert C-C
          OC for e31 -- this control confirms why, reusing the
          already-existing operator, not building a new one).

PRE-REGISTERED PASS/FAIL (SENTINEL PASS_ALL requires ALL of):
  G0 (sequential round-trip, exact data, no injected noise): the
     sequential recipe run on nature's own exact forward frequencies
     (no measurement noise anywhere) must recover e31 to <1e-6
     relative -- i.e. E_HAT's own <1e-6 self-consistency (job 2491531's
     bar) must not get amplified into a large e31 bias when chained
     into step 2. A FAIL here (>1e-6) would mean the sequential chain
     is numerically unstable on its own, independent of any
     measurement-noise question.
  G1 (root quality, both inversions): genuine sign flip of the
     shooting residual across each search bracket, PLUS a hi-side-only
     negative-control bracket with NO sign flip, for BOTH the step-1 E
     inversion and the step-2 e31 inversion.
  G_branch: the step-2 e31_hat lands on the PHYSICAL branch, i.e.
     within [1, 10] and closer to E31_TRUE=4.1 than to the reflection
     twin (C13E/C33E)*E33 ~ 8.95 (Sec 18.175 caveat) -- confirms the
     [1,10] bracket did not silently lock onto the wrong root.
  G_neg (standing prohibition, a gate not just advisory): the step-4
     negative control (C-C SC/OC as the e31 observable, E_HAT held)
     must show propagated relative e31 uncertainty at 0.01% frequency
     precision of AT LEAST 100% (order-100%, matching job 2491531's
     731%) -- i.e. the recipe visibly fails on that observable. A
     control that "passes" (small uncertainty) would mean this
     probe's C-C branch is bugged, not that C-C OC/SC secretly
     identifies e31.
  G2 (compound uncertainty, reported, not gated): the actual payoff --
     delta_e31_direct, delta_e31_from_E, and their quadrature-combined
     delta_e31_total, at both 0.01% and 0.1% assumed frequency
     precision.

Do not retune the C-C SC inverse or its bracket (PARAM_BRACKETS['E']
below is unchanged from job 2491531/probe_piezo_p4_inverse_recovery_
2026-09-16.py). Do not invert C-C OC for e31 as a real result -- the
step-4 control deliberately reuses cc_coupled_bisect AS the C-C OC
observable, per Sec 18.175 theorem A, to prove the prohibition, not to
work around it.

If lever 12 later replaces oc_ff_bisect with a same-ansatz operator,
re-score this probe on the new operator -- do not block on lever 12.

Worker-path note: this probe stays single-process. Each step's inner
bisection is already the expensive unit, and the steps below run
sequentially (step 2 needs step 1's own E_HAT, step 4 needs it too),
so the ProcessPoolExecutor + plate_solver.workers._piezo_oc_ff_root_
worker path used elsewhere in this project is not needed here -- this
probe follows probe_piezo_p4_ff_open_circuit_2026-09-16.py's serial
style, not probe_piezo_p4_inverse_recovery_2026-09-16.py's 2-worker
parallel style (that probe's two tasks were genuinely independent;
these are not).

Cost model (dps=100, sandbox-timed this session on the project
sandbox): one 40-iter cc_coupled_bisect call costs ~3.6-7s; one 50-iter
oc_ff_bisect call costs a small fraction of that (4x4 vs 6x6). Step 1
(E inversion, ~28 cc_coupled_bisect-equivalent calls) ~100-200s; step 2
(e31 inversion from OC, cheap 4x4) well under a minute; step 4
(negative control, another ~28 cc_coupled_bisect-equivalent calls)
~100-200s. Total wall time comfortably inside a 20-minute SLURM budget
(vs the sibling inverse probe's 30-minute request for two PARALLEL
7-minute tasks) -- this probe is fully serial so its total is closer to
the sum, not the max, of its stages.

SANDBOX RESULT (this session, dps=100, full production precision --
not a reduced smoke test): SENTINEL PASS_ALL. G0: E rel_err=4.470e-09,
sequential e31 rel_err=3.638e-07 (both <1e-6). G1: both stages sign-flip
+ negative control PASS. G_branch: e31_hat=4.1000015 (true 4.1, twin
8.9504). G_neg: negative control (C-C SC/OC observable) gives 730.8% at
0.01% freq (bar >=100%, matches job 2491531's 731% almost exactly --
E_HAT vs E_TRUE makes no visible difference on this already-broken
observable). G2 headline: at 0.01% frequency precision, the DIRECT
OC-frequency term alone gives 1.983% (matches Sec 18.175's 1.98%), but
the INDIRECT term from E_HAT's own uncertainty is 1.972% -- almost
exactly as large as the direct term, not a small correction -- so the
paper-honest COMPOUND sequential e31 uncertainty is 2.80% at 0.01% freq
(27.97% at 0.1% freq), not the single-observable 1.98%/19.8% that two
independent single-parameter demos would suggest.
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
DPS = int(os.environ.get("PIEZO_DPS", "100"))
SC_INNER_ITERS = int(os.environ.get("PIEZO_SC_INNER_ITERS", "40"))
OC_INNER_ITERS = int(os.environ.get("PIEZO_OC_INNER_ITERS", "50"))
OUTER_ITERS = int(os.environ.get("PIEZO_OUTER_ITERS", "26"))
G0_TOL_REL = 1e-6

# Duan2005 Table 1 (PZT4 piezo layer + steel host) -- same geometry/
# material as job 2491531 and probe_piezo_p4_ff_open_circuit_2026-09-16.py.
R_I, R_O, H = 0.1, 0.6, 0.01
E_TRUE, NU, RHO = 200e9, 0.3, 7800.0
C11E, C12E, C13E, C33E = 132e9, 71e9, 73e9, 115e9
E31_TRUE, E33 = -4.1, 14.1
X11, X33 = 7.124e-9, 5.841e-9
RHO_PZT = 7500.0
H1 = (2.0 / 12) * H  # Duan2005 Table 4's own h1/2h=1/12 ratio

REFLECTION_TWIN = (C13E / C33E) * E33  # Sec 18.175 caveat, ~8.95 C/m^2

BASE_KWARGS = dict(
    r_i=R_I, r_o=R_O, h=H, nu=NU, rho=RHO, h1=H1,
    C11E=C11E, C12E=C12E, C13E=C13E, C33E=C33E,
    e33=E33, X11=X11, X33=X33, rho_pzt=RHO_PZT, dps=DPS,
)

SC_OMEGA_ANCHOR = 2791.7866198564  # job 2491531 / Sec 18.173 anchor
OC_LO, OC_HI = 720.0, 800.0
E_LO, E_HI = 140e9, 260e9
E31_LO, E31_HI = -10.0, -1.0


# ---------- step 1: E from C-C SC (job 2491531's own operator) ----------

def _omega_sc(E, e31, iters=SC_INNER_ITERS, lo=None, hi=None):
    if lo is None or hi is None:
        center = SC_OMEGA_ANCHOR * (E / E_TRUE) ** 0.5
        lo, hi = 0.80 * center, 1.20 * center
    kw = dict(BASE_KWARGS); kw["E"] = E; kw["e31"] = e31
    return PiezoOutOfPlaneSolver(**kw).cc_coupled_bisect(lo, hi, 0, iters=iters)


def _invert_E(omega_target, e31_known, lo_p=E_LO, hi_p=E_HI,
              iters_inner=SC_INNER_ITERS, iters_outer=OUTER_ITERS):
    def fwd(E):
        return _omega_sc(E, e31_known, iters_inner)
    g_lo = fwd(lo_p) - omega_target
    g_hi = fwd(hi_p) - omega_target
    sign_flip = (g_lo > 0) != (g_hi > 0)
    lo, hi, glo = lo_p, hi_p, g_lo
    for _ in range(iters_outer):
        mid = 0.5 * (lo + hi)
        gm = fwd(mid) - omega_target
        if (gm > 0) == (glo > 0):
            lo, glo = mid, gm
        else:
            hi = mid
    return 0.5 * (lo + hi), sign_flip, g_lo, g_hi


# ---------- step 2: e31 from F-F n=0 OC, E held at an arbitrary value ----------

def _omega_oc(e31, E, iters=OC_INNER_ITERS):
    kw = dict(BASE_KWARGS); kw["E"] = E; kw["e31"] = e31
    return PiezoOutOfPlaneSolver(**kw).oc_ff_bisect(OC_LO, OC_HI, 0, iters=iters)


def _invert_e31_from_oc(omega_target, E_used, lo_p=E31_LO, hi_p=E31_HI,
                         iters_inner=OC_INNER_ITERS, iters_outer=OUTER_ITERS):
    def fwd(e31):
        return _omega_oc(e31, E_used, iters_inner)
    g_lo = fwd(lo_p) - omega_target
    g_hi = fwd(hi_p) - omega_target
    sign_flip = (g_lo > 0) != (g_hi > 0)
    lo, hi, glo = lo_p, hi_p, g_lo
    for _ in range(iters_outer):
        mid = 0.5 * (lo + hi)
        gm = fwd(mid) - omega_target
        if (gm > 0) == (glo > 0):
            lo, glo = mid, gm
        else:
            hi = mid
    return 0.5 * (lo + hi), sign_flip, g_lo, g_hi


# ---------- step 4: negative control -- e31 from C-C SC/OC, E held at E_HAT ----------

def _invert_e31_from_ccsc(omega_target, E_used, lo_p=E31_LO, hi_p=E31_HI,
                           iters_inner=SC_INNER_ITERS, iters_outer=OUTER_ITERS):
    center = SC_OMEGA_ANCHOR * (E_used / E_TRUE) ** 0.5
    lo_om, hi_om = 0.95 * center, 1.05 * center

    def fwd(e31):
        return _omega_sc(E_used, e31, iters_inner, lo=lo_om, hi=hi_om)
    g_lo = fwd(lo_p) - omega_target
    g_hi = fwd(hi_p) - omega_target
    sign_flip = (g_lo > 0) != (g_hi > 0)
    lo, hi, glo = lo_p, hi_p, g_lo
    for _ in range(iters_outer):
        mid = 0.5 * (lo + hi)
        gm = fwd(mid) - omega_target
        if (gm > 0) == (glo > 0):
            lo, glo = mid, gm
        else:
            hi = mid
    return 0.5 * (lo + hi), sign_flip, g_lo, g_hi


def main():
    t_all = time.time()
    print("job start: dps=%s sc_inner=%s oc_inner=%s outer=%s"
          % (DPS, SC_INNER_ITERS, OC_INNER_ITERS, OUTER_ITERS), flush=True)
    if ps.SOLVER_VERSION != EXPECT_SOLVER_VERSION:
        print("PREFLIGHT FAIL: SOLVER_VERSION=%r != %r"
              % (ps.SOLVER_VERSION, EXPECT_SOLVER_VERSION), flush=True)
        sys.exit(2)
    print("preflight OK: SOLVER_VERSION=%s" % ps.SOLVER_VERSION, flush=True)
    if not hasattr(PiezoOutOfPlaneSolver, "oc_ff_det"):
        print("PREFLIGHT FAIL: deployed PiezoOutOfPlaneSolver has no "
              "oc_ff_det (job 2491661's exact failure mode). Push "
              "plate_solver/piezo_solver.py to "
              "/home/ghmkfh/PythonMill/Plate_Solver_Package/plate_solver/ "
              "before re-running.", flush=True)
        sys.exit(2)
    print("preflight OK: oc_ff_det present", flush=True)

    log = []

    def note(msg):
        log.append(msg)
        print(msg, flush=True)

    # --- nature's exact frequencies (no measurement noise) ---
    t0 = time.time()
    omega_sc_true = _omega_sc(E_TRUE, E31_TRUE)
    note("nature omega_SC(E_true,e31_true) = %.10f rad/s (%.1fs)"
         % (omega_sc_true, time.time() - t0))
    t0 = time.time()
    omega_oc_true = _omega_oc(E31_TRUE, E_TRUE)
    note("nature omega_OC(e31_true; E_true) = %.10f rad/s (%.1fs)"
         % (omega_oc_true, time.time() - t0))

    # --- step 1: recover E_hat from C-C SC (e31 held known) ---
    t0 = time.time()
    E_hat, sf1, g_lo1, g_hi1 = _invert_E(omega_sc_true, E31_TRUE)
    g0_E_rel = abs(E_hat - E_TRUE) / E_TRUE
    note("STEP1 E_hat=%.10g E_true=%.6g rel_err=%.3e sign_flip=%s (%.1fs)"
         % (E_hat, E_TRUE, g0_E_rel, sf1, time.time() - t0))
    t0 = time.time()
    mid1 = 0.5 * (E_hat + E_HI)
    g_mid1 = _omega_sc(mid1, E31_TRUE) - omega_sc_true
    g_hi1b = _omega_sc(E_HI, E31_TRUE) - omega_sc_true
    neg1 = (g_mid1 > 0) == (g_hi1b > 0)
    g1_step1_pass = sf1 and neg1
    note("STEP1 G1: pos_flip=%s neg_control_no_flip=%s pass=%s (%.1fs)"
         % (sf1, neg1, g1_step1_pass, time.time() - t0))

    # --- step 2: recover e31_hat from F-F OC, holding E_hat (not E_true) ---
    t0 = time.time()
    e31_hat, sf2, g_lo2, g_hi2 = _invert_e31_from_oc(omega_oc_true, E_hat)
    g0_e31_rel = abs(e31_hat - E31_TRUE) / E31_TRUE
    note("STEP2 e31_hat=%.10g e31_true=%.6g rel_err=%.3e sign_flip=%s (%.1fs)"
         % (e31_hat, E31_TRUE, g0_e31_rel, sf2, time.time() - t0))
    t0 = time.time()
    mid2 = 0.5 * (e31_hat + E31_HI)
    g_mid2 = _omega_oc(mid2, E_hat) - omega_oc_true
    g_hi2b = _omega_oc(E31_HI, E_hat) - omega_oc_true
    neg2 = (g_mid2 > 0) == (g_hi2b > 0)
    g1_step2_pass = sf2 and neg2
    note("STEP2 G1: pos_flip=%s neg_control_no_flip=%s pass=%s (%.1fs)"
         % (sf2, neg2, g1_step2_pass, time.time() - t0))

    g0_pass = (g0_E_rel < G0_TOL_REL) and (g0_e31_rel < G0_TOL_REL)
    g1_pass = g1_step1_pass and g1_step2_pass
    g_branch_pass = ((E31_LO < e31_hat < E31_HI)
                      and (abs(e31_hat - E31_TRUE) < abs(e31_hat - REFLECTION_TWIN)))
    note("G0 pass=%s (E rel %.3e, sequential e31 rel %.3e); G1 pass=%s; "
         "G_branch pass=%s (e31_hat=%.7f vs true=4.1 vs twin=%.4f)"
         % (g0_pass, g0_E_rel, g0_e31_rel, g1_pass, g_branch_pass,
            e31_hat, REFLECTION_TWIN))

    # --- G2: compound uncertainty (reported, not gated) ---
    t0 = time.time()
    dE_fd = 0.01 * E_TRUE
    om_sc_plus = _omega_sc(E_TRUE + dE_fd, E31_TRUE)
    om_sc_minus = _omega_sc(E_TRUE - dE_fd, E31_TRUE)
    domega_sc_dE = (om_sc_plus - om_sc_minus) / (2 * dE_fd)
    note("G2 domega_SC/dE = %.6e rad/s per Pa (%.1fs)"
         % (domega_sc_dE, time.time() - t0))

    t0 = time.time()
    de31_fd = 0.01 * E31_TRUE
    om_oc_plus_e31 = _omega_oc(E31_TRUE + de31_fd, E_TRUE)
    om_oc_minus_e31 = _omega_oc(E31_TRUE - de31_fd, E_TRUE)
    domega_oc_de31 = (om_oc_plus_e31 - om_oc_minus_e31) / (2 * de31_fd)
    note("G2 domega_OC/de31 = %.6e rad/s per C/m^2 (%.1fs)"
         % (domega_oc_de31, time.time() - t0))

    t0 = time.time()
    om_oc_plus_E = _omega_oc(E31_TRUE, E_TRUE + dE_fd)
    om_oc_minus_E = _omega_oc(E31_TRUE, E_TRUE - dE_fd)
    domega_oc_dE = (om_oc_plus_E - om_oc_minus_E) / (2 * dE_fd)
    note("G2 domega_OC/dE = %.6e rad/s per Pa (%.1fs)"
         % (domega_oc_dE, time.time() - t0))

    de31_dE_implicit = -(domega_oc_dE / domega_oc_de31)
    note("G2 implicit de31/dE |omega_OC fixed = %.6e C/m^2 per Pa"
         % de31_dE_implicit)

    g2 = {}
    for rel_prec in (1e-4, 1e-3):
        dE_hat_prop = rel_prec * omega_sc_true / abs(domega_sc_dE)
        delta_e31_from_E = abs(de31_dE_implicit) * dE_hat_prop
        delta_e31_direct = rel_prec * omega_oc_true / abs(domega_oc_de31)
        delta_e31_total = (delta_e31_from_E ** 2 + delta_e31_direct ** 2) ** 0.5
        g2["%g" % rel_prec] = dict(
            assumed_omega_rel_precision=rel_prec,
            delta_E_hat=dE_hat_prop,
            delta_E_hat_rel=dE_hat_prop / E_TRUE,
            delta_e31_direct=delta_e31_direct,
            delta_e31_direct_rel=delta_e31_direct / E31_TRUE,
            delta_e31_from_E=delta_e31_from_E,
            delta_e31_from_E_rel=delta_e31_from_E / E31_TRUE,
            delta_e31_total=delta_e31_total,
            delta_e31_total_rel=delta_e31_total / E31_TRUE,
        )
        note("G2 @ %.4f%% freq: dE_hat=%.4f%% of E_true; "
             "de31_direct=%.4f%%, de31_from_E=%.4f%%, TOTAL=%.4f%% of e31_true"
             % (rel_prec * 100, dE_hat_prop / E_TRUE * 100,
                delta_e31_direct / abs(E31_TRUE) * 100,
                delta_e31_from_E / abs(E31_TRUE) * 100,
                delta_e31_total / abs(E31_TRUE) * 100))

    # --- step 4: negative control, C-C SC/OC as the e31 observable ---
    t0 = time.time()
    e31_hat_neg, sf_neg, g_lo_neg, g_hi_neg = _invert_e31_from_ccsc(
        omega_sc_true, E_hat)
    note("STEP4(neg) e31_hat=%.10g sign_flip=%s (%.1fs)"
         % (e31_hat_neg, sf_neg, time.time() - t0))
    t0 = time.time()
    de31n_fd = 0.01 * E31_TRUE
    om_sc_plus_e31 = _omega_sc(E_hat, E31_TRUE + de31n_fd)
    om_sc_minus_e31 = _omega_sc(E_hat, E31_TRUE - de31n_fd)
    domega_sc_de31 = (om_sc_plus_e31 - om_sc_minus_e31) / (2 * de31n_fd)
    note("G_neg domega_SC/de31 (at E_hat) = %.6e rad/s per C/m^2 (%.1fs)"
         % (domega_sc_de31, time.time() - t0))
    delta_e31_neg_001 = 1e-4 * omega_sc_true / abs(domega_sc_de31)
    delta_e31_neg_01 = 1e-3 * omega_sc_true / abs(domega_sc_de31)
    g_neg_pass = (delta_e31_neg_001 / abs(E31_TRUE)) >= 1.0
    note("G_neg propagated: 0.01%% freq -> %.2f%% of e31_true; 0.1%% freq -> "
         "%.2f%% (bar: >=100%% at 0.01%%) pass=%s"
         % (delta_e31_neg_001 / abs(E31_TRUE) * 100,
            delta_e31_neg_01 / abs(E31_TRUE) * 100, g_neg_pass))

    all_pass = g0_pass and g1_pass and g_branch_pass and g_neg_pass

    results = dict(
        E_true=E_TRUE, e31_true=E31_TRUE, h1_over_2h="1/12", h1=H1,
        omega_sc_true=omega_sc_true, omega_oc_true=omega_oc_true,
        E_hat=E_hat, g0_E_rel=g0_E_rel, g1_step1_pass=g1_step1_pass,
        e31_hat=e31_hat, g0_e31_rel=g0_e31_rel, g1_step2_pass=g1_step2_pass,
        g0_pass=g0_pass, g1_pass=g1_pass, g_branch_pass=g_branch_pass,
        reflection_twin=REFLECTION_TWIN,
        domega_sc_dE=domega_sc_dE, domega_oc_de31=domega_oc_de31,
        domega_oc_dE=domega_oc_dE, de31_dE_implicit=de31_dE_implicit,
        g2=g2,
        e31_hat_neg_control=e31_hat_neg,
        domega_sc_de31_at_Ehat=domega_sc_de31,
        g_neg_delta_e31_rel_0p01pct=delta_e31_neg_001 / E31_TRUE,
        g_neg_delta_e31_rel_0p1pct=delta_e31_neg_01 / E31_TRUE,
        g_neg_pass=g_neg_pass,
        elapsed_s=time.time() - t_all,
        log=log,
    )
    out_path = os.path.join(os.path.dirname(__file__) or ".",
                             "piezo_p4_sequential_e_e31_results_e31m.json")
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print("--- summary ---", flush=True)
    print("  E_hat rel_err=%.3e  sequential e31_hat rel_err=%.3e  "
          "e31 G2 direct/from_E/total @0.01%%=%.4f%%/%.4f%%/%.4f%%  "
          "G_neg@0.01%%=%.2f%%"
          % (g0_E_rel, g0_e31_rel,
             g2["0.0001"]["delta_e31_direct_rel"] * 100,
             g2["0.0001"]["delta_e31_from_E_rel"] * 100,
             g2["0.0001"]["delta_e31_total_rel"] * 100,
             delta_e31_neg_001 / abs(E31_TRUE) * 100), flush=True)
    print("SENTINEL " + ("PASS_ALL" if all_pass else "FAIL_ALL")
          + " g0=%s g1=%s g_branch=%s g_neg=%s"
          % (g0_pass, g1_pass, g_branch_pass, g_neg_pass), flush=True)


if __name__ == "__main__":
    main()
