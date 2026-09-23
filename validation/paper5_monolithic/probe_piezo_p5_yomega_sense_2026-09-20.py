# -*- coding: utf-8 -*-
"""
probe_piezo_p5_yomega_sense_2026-09-20.py

Paper 5 roadmap Sec 8 item 6, last of the four "optional-cheap extras"
(n=0 overtones and thickness sweep closed Sec 18.206/18.207; C-F/F-C
closed Sec 18.208/18.209). The voltage-driven flexural Y(omega) was
proven algebraically TRIVIAL this session (Sec 18.210: Y=j*omega*C0,
k_eff^2==0 identically, a direct corollary of the closed Phase 0
OC==SC theorem -- uniform face voltage cannot actuate flexure on this
homogeneous monolithic stacking). Per the user's explicit choice
(AskUserQuestion: "Build a force-driven sensing admittance instead"),
work redirected to a mechanically-driven, charge-out transfer
function. Derivation handed to Grok (pure math, per this project's own
grok-handoff discipline) and returned as
PAPER5_YOMEGA_SENSE_DERIVATION.md (Sec 18.211): the full-face charge
under a force drive is ALSO identically zero (the reciprocal finding,
read from the other electrical port -- same D_z(H) pointwise
cancellation as the voltage-driven case), but charge on a SEGMENT of
the still-grounded, still-fully-electroded face is not, and is the
well-posed, resonance-showing output.

New solver code this session (plate_solver/piezo_monolithic.py):
  _equilibrate_solve(A, b)     row+column-equilibrated mp.lu_solve,
                                needed because the driven 12x12's third
                                chi-cubic branch produces Bessel
                                arguments large enough (kr ~ 100 at
                                this geometry) that raw matrix entries
                                span >50 orders of magnitude and
                                mp.lu_solve reports "numerically
                                singular" without it, even though the
                                system is well-posed.
  _driven_vecs_at / driven_ff_force_sc / _driven_eval / Q_segment /
  Y_sense                      two-region (split at r_F), 12-unknown
                                linear system: 6 rim conditions
                                (M_rr=Q_r=phibar'=0 at each edge) + 6
                                interface conditions at r_F (continuity
                                of w, w', M_rr, phibar, phibar' + a
                                JUMP in Q_r of magnitude F/(2*pi*r_F),
                                derived from integrating the
                                equilibrium equation across a
                                vanishing strip, PAPER5_YOMEGA_SENSE_
                                DERIVATION.md Sec 3).

New worker this session (plate_solver/workers.py):
  _piezo_mono_driven_force_sc_worker  point evaluation of Y_sense from
                                       scalars, returns (Re, Im) as two
                                       floats (unlike every other
                                       worker in this module, this is
                                       not a root bisection).

PRE-REGISTERED PASS/FAIL (Grok's three validation gates,
PAPER5_YOMEGA_SENSE_DERIVATION.md Sec 8, re-derived/re-checked
independently before trusting):

  G1 identity(G) (hard gate, the sign check on the Q_r jump row):
     2*pi*int_{r_i}^{r_o} r*w(r) dr == -F/(A2*omega^2) exactly (up to
     mpmath quadrature/root precision). rel_err < 1e-25 at dps=60.
     This is the load-bearing proof that the jump sign in the RHS
     vector (-F/(2*pi*r_F) on the Q_r-bracket row, zero elsewhere) is
     correct, not assumed.

  G2 pole location (unique-continuation gate): the driven 12x12,
     evaluated as a homogeneous determinant (not via the RHS-solve
     path), must vanish at EXACTLY the same omega as the
     already-validated coupled_bisect F-F SC fundamental (Sec
     18.199-18.209): a homogeneous, fully-continuous 12-vector
     solution at F=0 is just a restriction of one global 6x6
     eigenmode. rel_err < 1e-10.

  G3 e31_bar -> 0 limit: Q_in -> 0 linearly in e31_bar, IF e31 and e33
     are scaled TOGETHER. This gate has a documented trap: e31_bar =
     e31 - (C13E/C33E)*e33 (PAPER5_DERIVATION.md Sec 1), so scaling
     e31 alone while leaving e33 fixed does NOT drive e31_bar to zero
     at this material point (it plateaus near -8.95 instead of 0,
     since e31_bar is dominated by the e33 cross term here). That was
     a real bug in an earlier interactive probe this session -- not in
     piezo_monolithic.py -- and is reported below (G_mu) as its own
     gate so the mistake and its resolution are both on record, not
     just the corrected result.
     Pass: log-log slope of |Q_in| vs scale over the three deepest
     decades (1e-2->1e-3->1e-4) is in [0.99, 1.01] (clean linear
     vanishing). The 1.0->1e-1 decade is NOT gated on this criterion:
     at scale=1.0, e31_bar~O(1) is not yet a small perturbation, so
     the asymptotic linear rate is not expected there (observed ratio
     ~7.5 instead of ~10 at that decade in this session -- reported,
     not gated).

  G_mu (diagnostic, explains the interrupted-session discrepancy):
     with e31 and e33 scaled together down to 1e-8, the two
     non-escaping chi-cubic branches from _branches() must match
     +-mu (the pure elastic _elastic_k(omega) value) to rel_err <
     1e-9. This is the check that would have caught the test-harness
     bug immediately: scaling e31 alone gives branches converging to
     a DIFFERENT, wrong finite value (~9.557 at this omega, not
     10.250145...) because e31_bar itself never approached zero.

  G_worker (hard gate): the new worker reproduces Y_sense from plain
     scalars BIT-IDENTICALLY (Re and Im both, mpmath value equality)
     at a generic point and near the pole (large-magnitude,
     sign-flipping region) -- required before this can run under the
     cluster's ProcessPoolExecutor.

Also reported (not gated): Y_sense(omega) swept across the fundamental
F-F SC resonance, to show the expected pole (large |Y|, sign flip).

SCOPE, per PAPER5_YOMEGA_SENSE_DERIVATION.md Sec 9/10: n=0 only (no
higher circumferential order this pass); does not restore the
voltage-driven k_eff (that path stays identically zero, Sec 18.210);
is not an OC split (Phase 0's OC==SC theorem is unaffected -- this is
a different, mechanical, drive port entirely); is not a joint
stiffness/e31 inverse; is not in-plane/membrane. No SOLVER_VERSION
bump (new unused-by-default entry points). No .tex, no FE.
"""
import json
import os
import sys
import time

sys.path.insert(0, os.environ.get("PKG_PATH",
                                  os.path.join(os.path.dirname(__file__), "..")))
import plate_solver as ps  # noqa: E402
from plate_solver.piezo_monolithic import PiezoMonolithicOutOfPlaneSolver  # noqa: E402
from plate_solver.workers import (  # noqa: E402
    _piezo_mono_driven_force_sc_worker,
)
from mpmath import mp, mpf, matrix, quad  # noqa: E402

EXPECT_SOLVER_VERSION = os.environ.get("EXPECT_SOLVER_VERSION", "2026-07-10.s10")
DPS = int(os.environ.get("PIEZO_DPS", "60"))
DPS_HI = int(os.environ.get("PIEZO_DPS_HI", "80"))
INNER_ITERS = int(os.environ.get("PIEZO_INNER_ITERS", "60"))

# Duan Table 1 / Paper 4 BASE_KWARGS (Sec 18.200), Duan radii (Sec
# 18.201), full ceramic thickness 2H = 0.02 m -- same fixed geometry/
# material as every other Paper 5 probe this session.
R_I, R_O, H = 0.1, 0.6, 0.01
C11E, C12E, C13E, C33E = 132e9, 71e9, 73e9, 115e9
RHO = 7500.0
E31, E33 = 4.1, 14.1
X11, X33 = 7.124e-9, 5.841e-9

OMEGA = 200.0          # Grok's suggested first numeric point
OMEGA_NEAR_POLE = 473.0
R_F = 0.35
R_STAR = 0.30
CP_FF_LO, CP_FF_HI = 460.0, 490.0    # bracket for the known F-F SC fundamental

BASE_KWARGS = dict(
    r_i=R_I, r_o=R_O, H=H,
    C11E=C11E, C12E=C12E, C13E=C13E, C33E=C33E, rho=RHO,
    e31=E31, e33=E33, X11=X11, X33=X33, dps=DPS,
)


def _solver(dps=None, e31=None, e33=None):
    kw = dict(BASE_KWARGS)
    if dps is not None:
        kw["dps"] = dps
    if e31 is not None:
        kw["e31"] = e31
    if e33 is not None:
        kw["e33"] = e33
    return PiezoMonolithicOutOfPlaneSolver(**kw)


def _driven12_det(s, omega, dps):
    """Homogeneous driven 12x12 determinant (F=0 path), for the
    unique-continuation pole-location gate -- built the same way as
    driven_ff_force_sc's matrix, minus the RHS."""
    with mp.workdps(dps):
        omega = mpf(omega)
        d = s._d(); A1v = s._A1(); Hm = mpf(s.H)
        e31_bar = s._e31_bar()
        K_pref = s._K_pref()  # 2026-09-23 Sec 18.231: projection-dependent
        lams = s._branches(omega)

        def vecs(r):
            return s._driven_vecs_at(0, r, lams, d, A1v, K_pref)

        w_i, wp_i, phi_i, m_i, q_i, phip_i = vecs(s.r_i)
        w_o, wp_o, phi_o, m_o, q_o, phip_o = vecs(s.r_o)
        w_F, wp_F, phi_F, m_F, q_F, phip_F = vecs(R_F)
        zero6 = [mpf(0)] * 6

        def interface_row(vec):
            return list(vec) + [-x for x in vec]

        rows = [
            m_i + zero6, q_i + zero6, phip_i + zero6,
            zero6 + m_o, zero6 + q_o, zero6 + phip_o,
            interface_row(w_F), interface_row(wp_F), interface_row(m_F),
            interface_row(q_F), interface_row(phi_F), interface_row(phip_F),
        ]
        M = matrix(rows)
        for i in range(12):
            sc = max(abs(M[i, j]) for j in range(12)) or mpf(1)
            for j in range(12):
                M[i, j] = M[i, j] / sc
        return mp.det(M)


def main():
    t_all = time.time()
    print("job start: dps=%s dps_hi=%s inner_iters=%s"
          % (DPS, DPS_HI, INNER_ITERS), flush=True)
    if ps.SOLVER_VERSION != EXPECT_SOLVER_VERSION:
        print("PREFLIGHT FAIL: SOLVER_VERSION=%r != %r"
              % (ps.SOLVER_VERSION, EXPECT_SOLVER_VERSION), flush=True)
        sys.exit(2)
    print("preflight OK: SOLVER_VERSION=%s" % ps.SOLVER_VERSION, flush=True)
    if not hasattr(PiezoMonolithicOutOfPlaneSolver, "Y_sense"):
        print("PREFLIGHT FAIL: deployed PiezoMonolithicOutOfPlaneSolver "
              "has no Y_sense. Push plate_solver/piezo_monolithic.py "
              "and plate_solver/workers.py before re-running.", flush=True)
        sys.exit(2)
    print("preflight OK: Y_sense present", flush=True)

    log = []

    def note(msg):
        log.append(msg)
        print(msg, flush=True)

    s = _solver(dps=DPS)

    # ---------- G1: identity (G), force-balance sign gate ----------
    t0 = time.time()
    res = s.driven_ff_force_sc(mpf(OMEGA), R_F, F=1.0, n=0)
    with mp.workdps(DPS):
        def w_of_r(r):
            w, _ = s._driven_eval(res, r)
            return r * w
        lhs = 2 * mp.pi * quad(w_of_r, [s.r_i, R_F, s.r_o])
        A2 = s._A2()
        rhs = -mpf(1.0) / (A2 * mpf(OMEGA) ** 2)
        rel_err_g1 = abs((lhs - rhs) / rhs)
    g1_pass = rel_err_g1 < mpf("1e-25")
    note("G1 identity(G): lhs=%s rhs=%s rel_err=%s -> %s (%.1fs)"
         % (lhs, rhs, rel_err_g1, g1_pass, time.time() - t0))

    # ---------- G2: pole location vs already-validated coupled_bisect ----------
    t0 = time.time()
    omega_coupled = s.coupled_bisect(CP_FF_LO, CP_FF_HI, 0, iters=INNER_ITERS)
    with mp.workdps(DPS):
        a, b = mpf(CP_FF_LO), mpf(CP_FF_HI)
        flo = _driven12_det(s, a, DPS)
        for _ in range(INNER_ITERS):
            m = (a + b) / 2
            fm = _driven12_det(s, m, DPS)
            if (fm.real > 0) == (flo.real > 0):
                a, flo = m, fm
            else:
                b = m
        omega_driven12 = float((a + b) / 2)
    rel_err_g2 = abs(omega_driven12 - omega_coupled) / omega_coupled
    g2_pass = rel_err_g2 < 1e-10
    note("G2 pole: coupled_bisect=%.13f driven12=%.13f rel_err=%s -> %s (%.1fs)"
         % (omega_coupled, omega_driven12, rel_err_g2, g2_pass, time.time() - t0))

    # ---------- G3: Q_in -> 0 linearly as e31_bar -> 0 (scaled correctly) ----------
    t0 = time.time()
    scales = [1.0, 1e-1, 1e-2, 1e-3, 1e-4]
    qvals = []
    for scale in scales:
        s2 = _solver(dps=DPS, e31=E31 * scale, e33=E33 * scale)
        res2 = s2.driven_ff_force_sc(mpf(OMEGA), R_F, F=1.0, n=0)
        with mp.workdps(DPS):
            q_in = s2.Q_segment(res2, s2.r_i, R_STAR)
        qval = complex(q_in).real
        qvals.append(qval)
        note("G3 scale=%.0e Q_in=%s" % (scale, qval))
    ratios = [qvals[i] / qvals[i + 1] for i in range(len(qvals) - 1)]
    note("G3 decade ratios (Q_in[k]/Q_in[k+1]): %s" % (ratios,))
    # only the three deepest decades are gated -- scale=1.0 is not yet
    # a small perturbation, see docstring.
    deep_ratios = ratios[1:]
    import math
    slope = (math.log(abs(qvals[2])) - math.log(abs(qvals[4]))) / \
            (math.log(scales[2]) - math.log(scales[4]))
    g3_pass = all(9.9 <= r <= 10.1 for r in deep_ratios)
    note("G3 log-log slope (scale=1e-2..1e-4) = %.6f -> %s (%.1fs)"
         % (slope, g3_pass, time.time() - t0))

    # ---------- G_mu: mechanical branches -> +-elastic mu (diagnostic) ----------
    t0 = time.time()
    s_hi = _solver(dps=DPS_HI)
    with mp.workdps(DPS_HI):
        mu, _ = s_hi._elastic_k(mpf(OMEGA))
    mu_f = float(mu)
    s_small = _solver(dps=DPS_HI, e31=E31 * 1e-8, e33=E33 * 1e-8)
    with mp.workdps(DPS_HI):
        lams_small = s_small._branches(mpf(OMEGA))
    lam_mech = sorted([complex(l) for _, l in lams_small],
                       key=lambda z: abs(z))[:2]
    rel_errs_mu = [abs(abs(z) - mu_f) / mu_f for z in lam_mech]
    g_mu_pass = all(e < 1e-9 for e in rel_errs_mu)
    note("G_mu elastic mu=%.15f  small-e31_bar branches=%s  rel_errs=%s "
         "-> %s (%.1fs)" % (mu_f, lam_mech, rel_errs_mu, g_mu_pass,
                             time.time() - t0))
    note("G_mu note: this is the check that catches the "
         "'scale e31 but not e33' harness bug -- that mistake gives "
         "branches converging to ~9.557 (WRONG) instead of %.6f "
         "(correct), because e31_bar = e31-(C13E/C33E)*e33 never "
         "actually reaches 0 when only e31 is scaled." % mu_f)

    # ---------- G_worker: bit-identical worker path ----------
    t0 = time.time()
    y_direct = complex(s.Y_sense(mpf(OMEGA), R_F, R_STAR, F=1.0, n=0))
    args = (R_I, R_O, H, C11E, C12E, C13E, C33E, RHO,
            E31, E33, X11, X33, DPS, OMEGA, R_F, R_STAR, 1.0, 0)
    re, im = _piezo_mono_driven_force_sc_worker(args)
    worker_pass_pt = (re == y_direct.real) and (im == y_direct.imag)
    note("G_worker at omega=%.1f: direct=%s worker=%s -> %s"
         % (OMEGA, y_direct, complex(re, im), worker_pass_pt))

    y_direct2 = complex(s.Y_sense(mpf(OMEGA_NEAR_POLE), R_F, R_STAR, F=1.0, n=0))
    args2 = (R_I, R_O, H, C11E, C12E, C13E, C33E, RHO,
             E31, E33, X11, X33, DPS, OMEGA_NEAR_POLE, R_F, R_STAR, 1.0, 0)
    re2, im2 = _piezo_mono_driven_force_sc_worker(args2)
    worker_pass_pole = (re2 == y_direct2.real) and (im2 == y_direct2.imag)
    note("G_worker near pole omega=%.1f: direct=%s worker=%s -> %s (%.1fs)"
         % (OMEGA_NEAR_POLE, y_direct2, complex(re2, im2),
            worker_pass_pole, time.time() - t0))
    worker_pass = worker_pass_pt and worker_pass_pole

    # ---------- reported only: Y_sense sweep across the F-F SC pole ----------
    t0 = time.time()
    sweep = []
    for omega_s in [460.0, 468.0, 471.0, omega_coupled - 0.5,
                    omega_coupled - 0.01, omega_coupled + 0.01,
                    omega_coupled + 0.5, 478.0, 485.0]:
        y = complex(s.Y_sense(mpf(omega_s), R_F, R_STAR, F=1.0, n=0))
        sweep.append((omega_s, y.real, y.imag))
        note("sweep omega=%.4f Y_sense=%s" % (omega_s, y))
    note("sweep done (%.1fs)" % (time.time() - t0,))

    all_pass = bool(g1_pass and g2_pass and g3_pass and g_mu_pass and worker_pass)

    results = dict(
        g1_pass=bool(g1_pass), rel_err_g1=str(rel_err_g1),
        g2_pass=bool(g2_pass), rel_err_g2=rel_err_g2,
        omega_coupled=omega_coupled, omega_driven12=omega_driven12,
        g3_pass=bool(g3_pass), qvals=qvals, ratios=ratios, slope=slope,
        g_mu_pass=bool(g_mu_pass), mu=mu_f,
        lam_mech=[[z.real, z.imag] for z in lam_mech],
        rel_errs_mu=rel_errs_mu,
        worker_pass=bool(worker_pass),
        y_sense_at_200=[y_direct.real, y_direct.imag],
        y_sense_near_pole=[y_direct2.real, y_direct2.imag],
        sweep=sweep,
        all_pass=all_pass, elapsed_s=time.time() - t_all, log=log,
    )
    out_path = os.path.join(os.path.dirname(__file__) or ".",
                             "piezo_p5_yomega_sense_results.json")
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)

    print("--- summary ---", flush=True)
    print("  G1=%s G2=%s G3=%s G_mu=%s G_worker=%s"
          % (g1_pass, g2_pass, g3_pass, g_mu_pass, worker_pass), flush=True)
    print("SENTINEL " + ("PASS_ALL" if all_pass else "FAIL_ALL")
          + " (%.1fs total)" % (time.time() - t_all), flush=True)


if __name__ == "__main__":
    main()
