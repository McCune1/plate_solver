# -*- coding: utf-8 -*-
"""
probe_piezo_p4_ff_open_circuit_2026-09-16.py

Paper 4, PAPER4_PIEZO_ROADMAP.md Risk 6 / Sec 6 follow-up to job 2491531
(C-C short-circuit inverse: E recoverable, e31 unidentifiable at 731%
relative uncertainty from 0.01% frequency precision).

WHAT THIS IS: the open-circuit (OC) electrical-BC counterpart of the
short-circuit (SC) sinusoidal-potential model, plus the IEEE
resonance/antiresonance observable (item 3 of the 2026-09-16 next-step
fork) computed from the same pair of frequencies.

WHY F-F n=0, NOT C-C: fully-electroded Kirchhoff flexure with a spatially
constant electrode voltage V produces a spatially constant piezo moment.
That moment does not enter the interior PDE (Q_r is a derivative of M)
and does not enter a C-C boundary (w=w'=0). Green's identity further
forces int_A nabla^2 w dA = oint w,n ds = 0 when w'=0 on both edges, so
the charge constraint Q=0 implies V=0: C-C OC coincides with C-C SC at
this order. F-F n=0 is the mode for which S = r_o w'(r_o) - r_i w'(r_i)
is nonzero, V is determined by Q=0, and V enters M_rr=0 at the free
edges. See LESSONS_LEARNED.md Sec 18.175.

MODEL (leading-order OC, linear through-thickness potential):
  phi_top    = V(t) * (z-h)/h1
  phi_bottom = -V(t) * (z+h)/h1     (same poling, parallel bimorph)
  Q_outers = 0  =>  V = - e31_bar (h+h1) h1 / (Xi33_bar A) * int nabla^2 w dA
  M_rr gets a constant extra 2 e31_bar (h + h1/2) V
  assembled as a rank-1 update on both M_rr rows of the elastic F-F 4x4,
  alpha = 4 e31_bar^2 (h+h1/2)(h+h1) h1 / (Xi33_bar (r_o^2-r_i^2))
  Mrow_OC = Mrow_elastic + alpha * S,  S_j = r_o Z_j'(r_o) - r_i Z_j'(r_i)
Sign of alpha was NOT guessed: both signs were run, and +alpha is the
one that RAISES the n=0 root (piezoelectric stiffening).

IEEE k_eff (item 3, resonance/antiresonance, no full Y(omega) sweep):
  f_r = short-circuit / constant-E frequency of this linear-potential
        model = the elastic bilayer (V=0 => E_z=0 => no electric moment)
  f_a = open-circuit frequency (this probe)
  k_eff^2 = 1 - (f_r / f_a)^2
The sine-field SC 6x6 (coupled_bisect) sits 0.001% above the elastic
bilayer; using it as f_r instead of the elastic bilayer changes k_eff
at the noise floor of this effect. Both are reported.

PRE-REGISTERED PASS/FAIL:
  G0 (elastic limit): n!=0 OC det equals elastic det at a sample omega
     (identically, not just at a root); h1=0 OC det equals elastic det.
  G_stiffen: n=0 OC root > n=0 elastic-bilayer root, relative split in
     (0.2%, 0.5%) at h1/2h=1/12 -- order-of-magnitude check against
     e31_bar^2/Xi33_bar / c11_bar * d2/(d1+d2) ~ 0.3%.
  G1 (root quality): genuine sign flip of G(e31) = omega_OC(e31) -
     omega_target across the e31 search bracket, PLUS a hi-side-only
     negative control with no flip.
  G2 (reported, not gated): forward FD d omega / d e31 at e31_true;
     propagated relative e31 uncertainty at 0.01% and 0.1% frequency
     precision. PRE-REGISTERED EXPECTATION: 0.01% freq maps to a few
     percent of e31_true (order 1-10%, NOT 700%). A large uncertainty
     here would mean the OC model did not actually unlock e31, and
     would be a probe-design or sign bug, not an expected finding.
  G3 (k_eff, reported): k_eff^2 in (0.002, 0.02) at this geometry --
     a few-tenths-of-a-percent frequency split, squared.
  G4 (inverse G0): recover e31_true from its own exact OC frequency
     to <1e-6 relative (same bar as job 2491531 / Southwell).

Cost: 4x4 Bessel, same class as elastic_bisect. Sandbox-runnable.
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

OC_LO, OC_HI = 720.0, 800.0
E31_LO, E31_HI = 1.0, 10.0


def _solver(e31):
    kw = dict(BASE_KWARGS)
    kw["e31"] = e31
    return PiezoOutOfPlaneSolver(**kw)


def _omega_oc(e31, iters=ITERS):
    return _solver(e31).oc_ff_bisect(OC_LO, OC_HI, 0, iters=iters)


def _invert_e31(omega_target, lo_p, hi_p, iters_inner, iters_outer):
    g_lo = _omega_oc(lo_p, iters_inner) - omega_target
    g_hi = _omega_oc(hi_p, iters_inner) - omega_target
    sign_flip = (g_lo > 0) != (g_hi > 0)
    lo, hi, glo = lo_p, hi_p, g_lo
    for _ in range(iters_outer):
        mid = 0.5 * (lo + hi)
        gm = _omega_oc(mid, iters_inner) - omega_target
        if (gm > 0) == (glo > 0):
            lo, glo = mid, gm
        else:
            hi = mid
    return 0.5 * (lo + hi), sign_flip, g_lo, g_hi


def main():
    t_all = time.time()
    print("job start: dps=%s iters=%s outer_iters=%s" % (DPS, ITERS, OUTER_ITERS),
          flush=True)
    if ps.SOLVER_VERSION != EXPECT_SOLVER_VERSION:
        print("PREFLIGHT FAIL: SOLVER_VERSION=%r != %r"
              % (ps.SOLVER_VERSION, EXPECT_SOLVER_VERSION), flush=True)
        sys.exit(2)
    print("preflight OK: SOLVER_VERSION=%s" % ps.SOLVER_VERSION, flush=True)
    if not hasattr(PiezoOutOfPlaneSolver, "oc_ff_det"):
        print("PREFLIGHT FAIL: deployed PiezoOutOfPlaneSolver has no "
              "oc_ff_det (job 2491661's exact failure). Push "
              "plate_solver/piezo_solver.py AND plate_solver/workers.py "
              "to /home/ghmkfh/PythonMill/Plate_Solver_Package/plate_solver/ "
              "before re-running; the probe alone is not enough.",
              flush=True)
        sys.exit(2)
    print("preflight OK: oc_ff_det present", flush=True)

    log = []
    def note(msg):
        log.append(msg)
        print(msg, flush=True)

    s = _solver(E31_TRUE)

    # --- G0: identities ---
    t0 = time.time()
    om_sample = 747.9863950952586
    det_oc_n1 = s.oc_ff_det(om_sample, 1)
    det_el_n1 = s.elastic_det(om_sample, 1)
    g0_n1 = (det_oc_n1 == det_el_n1)
    s_bare = PiezoOutOfPlaneSolver(
        r_i=R_I, r_o=R_O, h=H, E=E_TRUE, nu=NU, rho=RHO, h1=0.0, dps=DPS)
    om_bare = 726.81185755
    g0_h1 = (s_bare.oc_ff_det(om_bare, 0) == s_bare.elastic_det(om_bare, 0))
    g0_pass = g0_n1 and g0_h1
    note("G0 identities: n!=0 oc==elastic %s; h1=0 oc==elastic %s pass=%s (%.1fs)"
         % (g0_n1, g0_h1, g0_pass, time.time() - t0))

    # --- G_stiffen ---
    t0 = time.time()
    omega_el = s.elastic_bisect(OC_LO, OC_HI, 0, iters=ITERS)
    omega_oc = s.oc_ff_bisect(OC_LO, OC_HI, 0, iters=ITERS)
    drel = (omega_oc - omega_el) / omega_el
    g_stiffen = (omega_oc > omega_el) and (0.002 < drel < 0.005)
    note("G_stiffen: omega_el=%.10f omega_oc=%.10f drel=%.6e pass=%s (%.1fs)"
         % (omega_el, omega_oc, drel, g_stiffen, time.time() - t0))

    # --- G3 k_eff (item 3) ---
    t0 = time.time()
    k2_el = 1.0 - (omega_el / omega_oc) ** 2
    try:
        omega_sc_sine = s.coupled_bisect(OC_LO, OC_HI, 0, iters=ITERS)
    except Exception as e:
        omega_sc_sine = float("nan")
        note("G3 coupled_bisect failed: %s" % e)
    k2_sine = 1.0 - (omega_sc_sine / omega_oc) ** 2 if omega_sc_sine == omega_sc_sine else float("nan")
    g3_pass = 0.002 < k2_el < 0.02
    note("G3 k_eff: k2_from_elastic=%.6e k2_from_sine_SC=%.6e "
         "omega_sc_sine=%.10f pass=%s (%.1fs)"
         % (k2_el, k2_sine, omega_sc_sine, g3_pass, time.time() - t0))

    # --- G4 inverse G0 + G1 + G2 on e31 ---
    t0 = time.time()
    omega_target = omega_oc
    p_hat, sign_flip, g_lo, g_hi = _invert_e31(
        omega_target, E31_LO, E31_HI, ITERS, OUTER_ITERS)
    g4_rel = abs(p_hat - E31_TRUE) / abs(E31_TRUE)
    g4_pass = sign_flip and (g4_rel < G0_TOL_REL)
    note("G4 inverse G0: p_hat=%.10g p_true=%s rel_err=%.3e sign_flip=%s pass=%s (%.1fs)"
         % (p_hat, E31_TRUE, g4_rel, sign_flip, g4_pass, time.time() - t0))

    t0 = time.time()
    mid_p = 0.5 * (p_hat + E31_HI)
    g_mid = _omega_oc(mid_p, ITERS) - omega_target
    g_hi_side = _omega_oc(E31_HI, ITERS) - omega_target
    neg_control = (g_mid > 0) == (g_hi_side > 0)
    g1_pass = sign_flip and neg_control
    note("G1 root quality: pos_flip=%s neg_control_no_flip=%s pass=%s (%.1fs)"
         % (sign_flip, neg_control, g1_pass, time.time() - t0))

    t0 = time.time()
    dp = 0.01 * E31_TRUE
    om_plus = _omega_oc(E31_TRUE + dp, ITERS)
    om_minus = _omega_oc(E31_TRUE - dp, ITERS)
    domega_dp = (om_plus - om_minus) / (2 * dp)
    note("G2 forward FD: domega/de31 = %.6e rad/s per unit (%.1fs)"
         % (domega_dp, time.time() - t0))

    propagated = {}
    g2_ok = True
    for rel_prec in (1e-4, 1e-3):
        d_om = rel_prec * omega_target
        dp_prop = d_om / abs(domega_dp)
        rel_u = dp_prop / abs(E31_TRUE)
        propagated["%g" % rel_prec] = {
            "assumed_omega_rel_precision": rel_prec,
            "propagated_dp": dp_prop,
            "propagated_dp_rel_to_p_true": rel_u,
        }
        note("G2 propagated: %.4f%% freq -> de31=%.6g (%.4f%% of e31_true)"
             % (rel_prec * 100, dp_prop, rel_u * 100))
        if rel_prec == 1e-4 and not (0.001 < rel_u < 0.20):
            g2_ok = False

    all_pass = g0_pass and g_stiffen and g1_pass and g4_pass and g3_pass
    results = {
        "g0_pass": g0_pass,
        "g_stiffen_pass": g_stiffen,
        "g1_pass": g1_pass,
        "g3_k2_elastic": k2_el,
        "g3_k2_sine_sc": k2_sine,
        "g3_pass": g3_pass,
        "g4_p_hat": p_hat,
        "g4_rel_err": g4_rel,
        "g4_pass": g4_pass,
        "omega_elastic": omega_el,
        "omega_oc": omega_oc,
        "omega_sc_sine": omega_sc_sine,
        "drel_oc_vs_elastic": drel,
        "g2_domega_de31": domega_dp,
        "g2_propagated": propagated,
        "g2_expectation_ok": g2_ok,
        "elapsed_s": time.time() - t_all,
        "log": log,
    }
    out_path = os.path.join(os.path.dirname(__file__) or ".",
                            "piezo_p4_ff_open_circuit_results.json")
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print("--- summary ---", flush=True)
    print("  G0=%s G_stiffen=%s G1=%s G3=%s G4_rel=%s G2_0.01%%_unc=%.4f%%"
          % (g0_pass, g_stiffen, g1_pass, g3_pass, g4_rel,
             propagated["0.0001"]["propagated_dp_rel_to_p_true"] * 100),
          flush=True)
    print("SENTINEL " + ("PASS_ALL" if all_pass else "FAIL_ALL")
          + " all_pass=%s g2_expectation_ok=%s" % (all_pass, g2_ok),
          flush=True)


if __name__ == "__main__":
    main()
