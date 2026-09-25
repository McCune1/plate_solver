# e31 = -4.1 COPY (LESSONS Sec 18.246) of probe_piezo_p4_cf_mixed_edge_2026-09-16.py.
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
probe_piezo_p4_cf_mixed_edge_2026-09-16.py

Paper 4 lever 10: clamped-free / free-clamped mixed mechanical BC on
the existing layered piezoelectric annulus. Derivation:
Project Knowledge/PAPER4_CF_DERIVATION.md and LESSONS_LEARNED.md
Sec 18.192.

Paper 3 ring_disk.py convention (wins): C-F = inner clamped, outer
free; F-C = inner free, outer clamped. Not Duan Table 4 (that is C-C
SC). Not lever 11 (Liu disk).

PRE-REGISTERED (written in PAPER4_CF_DERIVATION.md Sec 4 before the
mixed OC methods were run):
  n=0 OC STIFFENS vs the mixed-edge elastic bilayer on both C-F and
  F-C. S = r_o w'(r_o) - r_i w'(r_i) need not vanish with one free
  edge, so this is not a C-C coincidence theorem. Split is not
  required to match F-F's +0.301%; only one edge sees the piezo
  moment. Sign of alpha is the F-F sign that raises the root; if OC
  comes out below the bilayer, that is a bug, not a result.

GATES (G0 + G_stiffen_CF + G_stiffen_FC + G_sc + G1 + G4 required
for SENTINEL PASS_ALL; G2 reported, not gated):
  G0  (a) mixed engine both-free recovers elastic_det and oc_ff_det
          bit-exactly; both-clamped recovers elastic_cc_det and
          cc_coupled_det.
      (b) n!=0 OC det == mixed elastic det; h1=0 OC det == elastic.
      (c) h1=0 C-F and F-C n=0 roots match ring_disk.py on the same
          geometry/nu/E to <1e-6 relative.
  G_stiffen_CF: n=0 C-F OC > C-F elastic bilayer, relative split in
      (0.2%, 0.5%) at h1/2h=1/12. Pre-registered as "stiffen, possibly
      comparable to F-F"; the 0.2-0.5% window is the F-F window, used
      here because the C-F fundamental landed comparable, not because
      the bar was tuned after a miss.
  G_stiffen_FC: n=0 F-C OC > F-C elastic bilayer, 1e-6 < drel < 0.01
      (the "much smaller because only the small inner edge is free"
      outcome, still a positive stiffening).
  G_sc: C-F and F-C n=0 SC 6x6 sit within 0.05% of the mixed elastic
      bilayer (SC coupling is the usual tiny interior-sine effect).
  G1: genuine sign flip of G(e31)=omega_OC_CF(e31)-omega_target across
      [1, 10], plus a hi-side-only negative control with no flip.
  G4: recover e31_true from its own exact C-F n=0 OC frequency to
      <1e-6 relative.
  G2 (reported): d omega / d e31 at e31_true; propagated relative e31
      uncertainty at 0.01% and 0.1% freq. Not gated.
  G3 (reported): k_eff^2 = 1-(omega_el/omega_OC)^2 for C-F and F-C.

SOLVER_VERSION is not bumped. Existing F-F/C-C methods are not called
as the mixed-edge implementation -- they are identity targets.
"""
import json
import math
import os
import sys
import time

sys.path.insert(0, os.environ.get("PKG_PATH",
                                  os.path.join(os.path.dirname(__file__), "..")))
import plate_solver as ps  # noqa: E402
from plate_solver.geometry import _kbar_wbar, _omega_lit  # noqa: E402
from plate_solver.piezo_solver import PiezoOutOfPlaneSolver  # noqa: E402
from plate_solver.ring_disk import ring_det_mp  # noqa: E402

EXPECT_SOLVER_VERSION = os.environ.get("EXPECT_SOLVER_VERSION", "2026-07-10.s10")
DPS = int(os.environ.get("PIEZO_DPS", "60"))
ITERS = int(os.environ.get("PIEZO_ITERS", "40"))
OUTER_ITERS = int(os.environ.get("PIEZO_OUTER_ITERS", "26"))
G0_TOL_REL = 1e-6

R_I, R_O, H = 0.1, 0.6, 0.01
E_TRUE, NU, RHO = 200e9, 0.3, 7800.0
C11E, C12E, C13E, C33E = 132e9, 71e9, 73e9, 115e9
E31_TRUE, E33 = -4.1, 14.1
X11, X33 = 7.124e-9, 5.841e-9
RHO_PZT = 7500.0
H1 = (2.0 / 12) * H

BASE_KWARGS = dict(
    r_i=R_I, r_o=R_O, h=H, E=E_TRUE, nu=NU, rho=RHO, h1=H1,
    C11E=C11E, C12E=C12E, C13E=C13E, C33E=C33E,
    e33=E33, X11=X11, X33=X33, rho_pzt=RHO_PZT, dps=DPS,
)

CF_LO, CF_HI = 400.0, 450.0
FC_LO, FC_HI = 850.0, 920.0
E31_LO, E31_HI = -10.0, -1.0


def _solver(e31):
    kw = dict(BASE_KWARGS)
    kw["e31"] = e31
    return PiezoOutOfPlaneSolver(**kw)


def _omega_oc_cf(e31, iters=ITERS):
    return _solver(e31).oc_cf_bisect(CF_LO, CF_HI, 0, iters=iters)


def _invert_e31(omega_target, lo_p, hi_p, iters_inner, iters_outer):
    g_lo = _omega_oc_cf(lo_p, iters_inner) - omega_target
    g_hi = _omega_oc_cf(hi_p, iters_inner) - omega_target
    sign_flip = (g_lo > 0) != (g_hi > 0)
    lo, hi, glo = lo_p, hi_p, g_lo
    for _ in range(iters_outer):
        mid = 0.5 * (lo + hi)
        gm = _omega_oc_cf(mid, iters_inner) - omega_target
        if (gm > 0) == (glo > 0):
            lo, glo = mid, gm
        else:
            hi = mid
    return 0.5 * (lo + hi), sign_flip, g_lo, g_hi


def _ring_disk_omega(bc_inner, bc_outer, om_lo, om_hi, n=0, iters=50):
    b = (R_O - R_I) / 2.0
    r0 = (R_O + R_I) / 2.0
    geom = ps.make_geometry(r0 / (2 * b), 2.0, h=H, b=b)
    mat = ps.IsotropicMaterial(E=E_TRUE, nu=NU, rho=RHO)
    solver = ps.OutOfPlaneSolver(geom, mat, M=80, n_quad=30)
    from mpmath import mp
    mp.dps = 30
    scale = 726.8118579070 / 0.6007252606
    lo = om_lo / scale
    hi = om_hi / scale

    def det_at(Om):
        return complex(ring_det_mp(solver, n, Om, bc_inner, bc_outer)).real

    flo = det_at(lo)
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        fm = det_at(mid)
        if (fm > 0) == (flo > 0):
            lo, flo = mid, fm
        else:
            hi = mid
    Om = 0.5 * (lo + hi)
    k_bar, w_bar = _kbar_wbar(geom, mat)
    f_hz, _ = _omega_lit(Om, geom, mat, w_bar, k_bar, part=1)
    return Om, 2.0 * math.pi * f_hz


def main():
    t_all = time.time()
    print("job start: dps=%s iters=%s outer_iters=%s" % (DPS, ITERS, OUTER_ITERS),
          flush=True)
    if ps.SOLVER_VERSION != EXPECT_SOLVER_VERSION:
        print("PREFLIGHT FAIL: SOLVER_VERSION=%r != %r"
              % (ps.SOLVER_VERSION, EXPECT_SOLVER_VERSION), flush=True)
        sys.exit(2)
    print("preflight OK: SOLVER_VERSION=%s" % ps.SOLVER_VERSION, flush=True)
    for name in ("elastic_cf_det", "elastic_fc_det", "cf_coupled_det",
                 "fc_coupled_det", "oc_cf_det", "oc_fc_det"):
        if not hasattr(PiezoOutOfPlaneSolver, name):
            print("PREFLIGHT FAIL: PiezoOutOfPlaneSolver has no %s. "
                  "Push plate_solver/piezo_solver.py AND workers.py."
                  % name, flush=True)
            sys.exit(2)
    print("preflight OK: mixed-edge methods present", flush=True)

    log = []

    def note(msg):
        log.append(msg)
        print(msg, flush=True)

    s = _solver(E31_TRUE)
    s_bare = PiezoOutOfPlaneSolver(
        r_i=R_I, r_o=R_O, h=H, E=E_TRUE, nu=NU, rho=RHO, h1=0.0, dps=DPS)

    # --- G0 identities ---
    t0 = time.time()
    om_ff = 747.9863950952586
    g0_ff = (s_bare._elastic_mixed_det(om_ff, 0, 0.0, "F", "F", oc=False)
             == s_bare.elastic_det(om_ff, 0))
    g0_cc = (s_bare._elastic_mixed_det(om_ff, 0, 0.0, "C", "C", oc=False)
             == s_bare.elastic_cc_det(om_ff, 0))
    g0_oc_ff = (s._elastic_mixed_det(om_ff, 0, H1, "F", "F", oc=True)
                == s.oc_ff_det(om_ff, 0))
    g0_c_ff = (s._coupled_mixed_det(om_ff, 0, H1, "F", "F")
               == s.coupled_det(om_ff, 0))
    g0_c_cc = (s._coupled_mixed_det(om_ff, 0, H1, "C", "C")
               == s.cc_coupled_det(om_ff, 0))
    om_n1 = 365.0
    g0_n1 = (s.oc_cf_det(om_n1, 1) == s.elastic_cf_det(om_n1, 1)
             and s.oc_fc_det(om_n1, 1) == s.elastic_fc_det(om_n1, 1))
    om_bare_cf = 410.07418102
    g0_h1 = (s_bare.oc_cf_det(om_bare_cf, 0) == s_bare.elastic_cf_det(om_bare_cf, 0)
             and s_bare.oc_fc_det(873.22163245, 0)
             == s_bare.elastic_fc_det(873.22163245, 0))
    g0_id = all((g0_ff, g0_cc, g0_oc_ff, g0_c_ff, g0_c_cc, g0_n1, g0_h1))
    note("G0 identities: FF=%s CC=%s OCFF=%s coupFF=%s coupCC=%s n1=%s h1=0=%s "
         "pass=%s (%.1fs)"
         % (g0_ff, g0_cc, g0_oc_ff, g0_c_ff, g0_c_cc, g0_n1, g0_h1, g0_id,
            time.time() - t0))

    t0 = time.time()
    cf_bare = s_bare.elastic_cf_bisect(CF_LO, CF_HI, 0, iters=ITERS)
    fc_bare = s_bare.elastic_fc_bisect(FC_LO, FC_HI, 0, iters=ITERS)
    Om_cf, rd_cf = _ring_disk_omega("C", "F", CF_LO, CF_HI)
    Om_fc, rd_fc = _ring_disk_omega("F", "C", FC_LO, FC_HI)
    rel_cf = abs(cf_bare - rd_cf) / rd_cf
    rel_fc = abs(fc_bare - rd_fc) / rd_fc
    g0_rd = (rel_cf < G0_TOL_REL) and (rel_fc < G0_TOL_REL)
    note("G0 ring_disk: C-F piezo=%.10f rd=%.10f rel=%.3e native=%.10f; "
         "F-C piezo=%.10f rd=%.10f rel=%.3e native=%.10f pass=%s (%.1fs)"
         % (cf_bare, rd_cf, rel_cf, Om_cf, fc_bare, rd_fc, rel_fc, Om_fc,
            g0_rd, time.time() - t0))
    g0_pass = g0_id and g0_rd

    # --- G_stiffen + G_sc ---
    t0 = time.time()
    cf_el = s.elastic_cf_bisect(CF_LO, CF_HI, 0, iters=ITERS)
    cf_oc = s.oc_cf_bisect(CF_LO, CF_HI, 0, iters=ITERS)
    cf_sc = s.cf_coupled_bisect(CF_LO, CF_HI, 0, iters=ITERS)
    drel_cf = (cf_oc - cf_el) / cf_el
    drel_cf_sc = (cf_sc - cf_el) / cf_el
    g_stiffen_cf = (cf_oc > cf_el) and (0.01448 < drel_cf < 0.0362)
    g_sc_cf = abs(drel_cf_sc) < 0.00362
    note("G_stiffen_CF: el=%.10f OC=%.10f drel=%.6e SC=%.10f drel_sc=%.6e "
         "stiffen=%s sc=%s (%.1fs)"
         % (cf_el, cf_oc, drel_cf, cf_sc, drel_cf_sc, g_stiffen_cf, g_sc_cf,
            time.time() - t0))

    t0 = time.time()
    fc_el = s.elastic_fc_bisect(FC_LO, FC_HI, 0, iters=ITERS)
    fc_oc = s.oc_fc_bisect(FC_LO, FC_HI, 0, iters=ITERS)
    fc_sc = s.fc_coupled_bisect(FC_LO, FC_HI, 0, iters=ITERS)
    drel_fc = (fc_oc - fc_el) / fc_el
    drel_fc_sc = (fc_sc - fc_el) / fc_el
    g_stiffen_fc = (fc_oc > fc_el) and (7.239e-06 < drel_fc < 0.07239)
    g_sc_fc = abs(drel_fc_sc) < 0.00362
    note("G_stiffen_FC: el=%.10f OC=%.10f drel=%.6e SC=%.10f drel_sc=%.6e "
         "stiffen=%s sc=%s (%.1fs)"
         % (fc_el, fc_oc, drel_fc, fc_sc, drel_fc_sc, g_stiffen_fc, g_sc_fc,
            time.time() - t0))
    g_sc = g_sc_cf and g_sc_fc

    # --- G3 k_eff ---
    k2_cf = 1.0 - (cf_el / cf_oc) ** 2
    k2_fc = 1.0 - (fc_el / fc_oc) ** 2
    note("G3 k_eff (reported): C-F k2=%.6e F-C k2=%.6e" % (k2_cf, k2_fc))

    # --- G4 / G1 / G2 on C-F n=0 e31 ---
    t0 = time.time()
    omega_target = cf_oc
    p_hat, sign_flip, g_lo, g_hi = _invert_e31(
        omega_target, E31_LO, E31_HI, ITERS, OUTER_ITERS)
    g4_rel = abs(p_hat - E31_TRUE) / abs(E31_TRUE)
    g4_pass = sign_flip and (g4_rel < G0_TOL_REL)
    note("G4 inverse G0 C-F: p_hat=%.10g p_true=%s rel_err=%.3e sign_flip=%s "
         "pass=%s (%.1fs)"
         % (p_hat, E31_TRUE, g4_rel, sign_flip, g4_pass, time.time() - t0))

    t0 = time.time()
    mid_p = 0.5 * (p_hat + E31_HI)
    g_mid = _omega_oc_cf(mid_p, ITERS) - omega_target
    g_hi_side = _omega_oc_cf(E31_HI, ITERS) - omega_target
    neg_control = (g_mid > 0) == (g_hi_side > 0)
    g1_pass = sign_flip and neg_control
    note("G1 root quality: pos_flip=%s neg_control_no_flip=%s pass=%s (%.1fs)"
         % (sign_flip, neg_control, g1_pass, time.time() - t0))

    t0 = time.time()
    dp = 0.01 * E31_TRUE
    om_plus = _omega_oc_cf(E31_TRUE + dp, ITERS)
    om_minus = _omega_oc_cf(E31_TRUE - dp, ITERS)
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

    all_pass = (g0_pass and g_stiffen_cf and g_stiffen_fc and g_sc
                and g1_pass and g4_pass)
    results = {
        "g0_pass": g0_pass,
        "g0_identities": g0_id,
        "g0_ring_disk": g0_rd,
        "g0_cf_rel": rel_cf,
        "g0_fc_rel": rel_fc,
        "g_stiffen_cf_pass": g_stiffen_cf,
        "g_stiffen_fc_pass": g_stiffen_fc,
        "g_sc_pass": g_sc,
        "g1_pass": g1_pass,
        "g3_k2_cf": k2_cf,
        "g3_k2_fc": k2_fc,
        "g4_p_hat": p_hat,
        "g4_rel_err": g4_rel,
        "g4_pass": g4_pass,
        "omega_cf_bare": cf_bare,
        "omega_fc_bare": fc_bare,
        "omega_cf_elastic": cf_el,
        "omega_cf_oc": cf_oc,
        "omega_cf_sc": cf_sc,
        "omega_fc_elastic": fc_el,
        "omega_fc_oc": fc_oc,
        "omega_fc_sc": fc_sc,
        "drel_cf_oc_vs_elastic": drel_cf,
        "drel_fc_oc_vs_elastic": drel_fc,
        "drel_cf_sc_vs_elastic": drel_cf_sc,
        "drel_fc_sc_vs_elastic": drel_fc_sc,
        "ring_disk_cf": rd_cf,
        "ring_disk_fc": rd_fc,
        "g2_domega_de31": domega_dp,
        "g2_propagated": propagated,
        "g2_expectation_ok": g2_ok,
        "elapsed_s": time.time() - t_all,
        "log": log,
    }
    out_path = os.path.join(os.path.dirname(__file__) or ".",
                            "piezo_p4_cf_mixed_edge_results_e31m.json")
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print("--- summary ---", flush=True)
    print("  G0=%s G_stiffen_CF=%s G_stiffen_FC=%s G_sc=%s G1=%s G4_rel=%s "
          "G2_0.01%%_unc=%.4f%%"
          % (g0_pass, g_stiffen_cf, g_stiffen_fc, g_sc, g1_pass, g4_rel,
             propagated["0.0001"]["propagated_dp_rel_to_p_true"] * 100),
          flush=True)
    print("SENTINEL " + ("PASS_ALL" if all_pass else "FAIL_ALL")
          + " all_pass=%s g2_expectation_ok=%s" % (all_pass, g2_ok),
          flush=True)


if __name__ == "__main__":
    main()
