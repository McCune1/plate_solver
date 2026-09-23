# -*- coding: utf-8 -*-
"""
plate_solver.piezo_solver -- Paper 4 (piezoelectric ring) F-F/C-C solver,
plus F-F open-circuit (Sec 18.175) and mixed-edge C-F/F-C (Sec 18.192).

Background: LESSONS_LEARNED.md Sec 18.141 (electroelastic derivation),
Sec 18.149 (F-F boundary rows + the corrected M_rr/M_theta_theta sign,
plus the live 5e-8-in-percent numerical cross-check against this
project's own ring_disk.py F-F elastic solver), Sec 18.148 (the C-C
cluster gate, job 2489140, PASS_ALL 27/27 Duan2005 Table 4 points), Sec
18.171 (this module's first, elastic-only F-F slice), Sec 18.172 (the
full F-F piezo-coupled 6x6 determinant), Sec 18.173 (this docstring's
newest addition: the C-C counterpart of both), PAPER4_PIEZO_ROADMAP.md
Sec 4 (the derivation) and Sec 8 (phased sequence -- step 2 is
"elastic-limit regression... before trusting any piezo-coupled number").

SCOPE: this module ports the four SC probe scripts' math into class form
-- two boundary conditions (F-F, C-C) times two model fidelities
(elastic-only, full piezo-coupled) -- plus the F-F open-circuit 4x4
(Sec 18.175, `oc_ff_det`/`oc_ff_bisect`) and the mixed-edge C-F/F-C
family (Sec 18.192 / PAPER4_CF_DERIVATION.md: elastic, SC coupled,
linear OC). C-F is inner clamped / outer free; F-C is the swap
(Paper 3 ring_disk.py convention).
`elastic_det()`/`elastic_bisect()` are the non-piezo-coupled F-F
determinant (bare host h1=0, or an elastic bilayer h1>0/e31=0 -- ordinary
two-material added stiffness, zero electromechanical coupling), matching
`Paper4_Piezo/probe_piezo_p4_elastic_ff_baseline_2026-09-15.py`'s
`det_at()`/`bisect()`. `coupled_det()`/`coupled_bisect()` are the full
3-branch piezo-coupled 6x6 F-F determinant (the chi/lambda cubic, the
electrical phi'=0 row, the piezoelectric term on M_rr/Q_r), matching
`Paper4_Piezo/probe_piezo_p4_ff_forward_model_2026-09-15.py`'s
`cubic_coeffs()`/`branches()`/`radial_quad()`/`ff_det()`/`bisect()`.
`elastic_cc_det()`/`elastic_cc_bisect()` and `cc_coupled_det()`/
`cc_coupled_bisect()` are the C-C (w=w'=phi'=0) analogues, matching
`Paper4_Piezo/probe_piezo_p4_elastic_cc_baseline_2026-09-15.py`'s
`det_at()`/`bisect()` and `Paper4_Piezo/probe_piezo_p4_cc_forward_model_
2026-09-15.py`'s `cc_det()`/`bisect()` respectively -- the C-C model is
the one with an actual literature table (Duan2005 Table 4) and its own
cluster-confirmed gate (Sec 18.148), so its port is the closer-to-
"production" of the two boundary conditions even though F-F was ported
first (F-F is this project's own primary BC target per
PAPER4_PIEZO_ROADMAP.md, C-C was always the validation vehicle). The
chi-cubic (`_cubic_coeffs`/`_branches`) and `_radial_quad` are shared,
unmodified, between F-F and C-C -- Duan2005's own derivation never makes
them depend on which mechanical edge condition is applied; only the
boundary-row construction (`coupled_det` vs `cc_coupled_det`) differs.

Architectural note (Sec 18.141 G): this lands as a NEW solver class in a
NEW module -- the same "new entry point, not a modification to the
validated OutOfPlaneSolver" pattern ring_disk.py set for Paper 3 -- so
Papers 1-3's code paths and SOLVER_VERSION stay untouched. SOLVER_VERSION
is not bumped: nothing here is on any existing default computation's path.

ELECTRICAL PROJECTION AND OC CHARGE ARM (2026-09-23, LESSONS_LEARNED.md
Sec 18.231, SC_PROJECTION_CONSISTENT_DERIVATION_2026-09-23.md Sec 3.2/4):
the default projection='consistent' is the same-weight Galerkin
projection of the electric enthalpy, with a quadratic bubble potential
in each skin (exact long-wave SC increment h1^3 e31bar^2/(12 Xi33bar)
per skin), and the variationally consistent (reciprocal) bus charge
Q = C0 V + 4 pi e31bar (h + h1/2) S, i.e. the SAME skin mid-surface arm
the moment row already uses; the interior-bubble field contributes no
bus charge (int f' dz = 0). projection='duan' is the pre-2026-09-23
model, kept bit-identical: Duan 2005 Eq. 12 charge-conservation closure
of the sine (skin SC increment 12/pi^2 high) and the pointwise
outer-face charge arm (h + h1), which overstates every OC split by
(h + h1)/(h + h1/2) -- the whole 93/90/86% FE/analytic gap of Sec 18.189.
projection='galerkin_sine' is a diagnostic (sine, consistent, 96/pi^4).
The duan option is what reproduces Duan2005's own tables as published.

Units: dimensional SI throughout (metres, Pa, kg/m^3, C/m^2, F/m, rad/s),
matching Duan2005 and both probe scripts directly -- NOT the package's
native Seok-Tiersten Omega nondimensionalization that geometry.py/
core_solvers.py use. The radial ODE each branch satisfies (Delta Z =
lambda Z, an ordinary/modified Bessel equation) is solved in closed form
(mpmath Bessel functions) rather than via the general Frobenius power
series (core_solvers.py's _series_mp) -- Sec 18.141 D/G's own point that
the piezo-coupled system decouples into Helmholtz-type radial problems,
which is why this class does not subclass or wrap OutOfPlaneSolver.

Constructor material-constant convention (Sec 18.172, superseding Sec
18.171's original choice): C11E/C12E/C13E/C33E are the piezo layer's RAW
elastic stiffness at constant E-field (Pa), i.e. Duan2005 Table 1
sourcing directly -- NOT the plane-stress-reduced c11_bar/c12_bar Sec
18.171's first cut took as C11E/C12E. The reduced constants (Sec 18.141
A: c11_bar = C11E - C13E^2/C33E, c12_bar = C12E - C13E^2/C33E, and, for
the coupled model, e31_bar/Xi33_bar/Xi11_bar) are formed internally, once,
by _c11_bar()/_c12_bar()/_e31_bar()/_Xi33_bar()/_Xi11_bar() -- exactly
mirroring how both probe scripts form them once at module scope. This
changed the elastic-only constructor's C11E/C12E meaning from Sec
18.171 (see that section for the superseded convention); the tests in
tests/test_solver.py were updated in the same session this change was
made, so no caller is left passing constants under the old convention.
"""
from __future__ import annotations

from mpmath import mp, mpf, matrix


class PiezoOutOfPlaneSolver:
    """Layered piezoelectric ring OOP flexural F-F/C-C solver for Paper 4.

    Builds the F-F/C-C x elastic/coupled families, F-F open-circuit,
    F-F driven admittance `driven_ff_qv` (Sec 18.179), and the mixed-
    edge C-F/F-C family (Sec 18.192; C-F = inner C, outer F):

    - `elastic_det()` / `elastic_bisect()`: the 4x4 PURELY ELASTIC F-F
      determinant (2 branches). h1=0 is the bare host ring (no piezo
      layer at all); h1>0 is an ordinary two-material elastic bilayer
      (added bending stiffness from the piezo layer, ZERO
      electromechanical coupling -- e31/permittivity constants are not
      needed for this path at all).
    - `coupled_det()` / `coupled_bisect()`: the full 6x6 PIEZO-COUPLED
      F-F determinant (3 branches from the chi/lambda cubic, Sec
      18.141 D). Needs h1>0 and the full piezoelectric constant set
      (e31, e33, X11, X33) in addition to the elastic ones.
    - `elastic_cc_det()` / `elastic_cc_bisect()`: the C-C (w=w'=0)
      counterpart of `elastic_det`/`elastic_bisect` -- same d1/d2/A2
      building blocks, different (w, w') boundary rows instead of
      (M_rr, Q_r). Reproduces Duan2005 Table 4's own h1=0 column and
      elastic-bilayer row (Sec 18.148's cluster gate target family).
    - `cc_coupled_det()` / `cc_coupled_bisect()`: the C-C (w=w'=phi'=0)
      counterpart of `coupled_det`/`coupled_bisect` -- reuses the same
      chi-cubic/branches/radial_quad machinery, different boundary rows.
      This is the model already CLUSTER-CONFIRMED against all 27 non-
      trivial Duan2005 Table 4 points (Sec 18.148, job 2489140) as a
      standalone script; this class just packages that same math.
    - `oc_ff_det()` / `oc_ff_bisect()`: F-F OPEN-CIRCUIT (fully electroded,
      linear through-thickness potential, Q=0 on the outer-electrode bus).
      Leading-order counterpart of short-circuit's sinusoidal-potential
      6x6: a 4x4 elastic F-F determinant plus a rank-1 update on the
      two M_rr rows, proportional to S = r_o w'(r_o) - r_i w'(r_i)
      (Green: S = (1/2pi) int_A nabla^2 w dA for n=0). n!=0 and h1=0
      reduce identically to `elastic_det` (net electrode charge cancels
      for n!=0; no layer to charge at h1=0). C-C fully-electroded OC is
      NOT a fifth mechanical-BC sibling: w'=0 on both edges forces S=0,
      so Q=0 implies V=0 and OC coincides with short-circuit / elastic
      bilayer at this order (LESSONS_LEARNED.md Sec 18.175). Do not add
      a C-C OC frequency solver expecting it to recover e31.
    - `elastic_cf_det()` / `elastic_fc_det()` and bisects: mixed-edge
      elastic 4x4. C-F = inner clamped (w=w'=0), outer free
      (M_rr=Q_r=0); F-C is the swap. Paper 3 ring_disk.py convention
      (Sec 18.192). Both-free recovers `elastic_det`; both-clamped
      recovers `elastic_cc_det`.
    - `cf_coupled_det()` / `fc_coupled_det()` and bisects: mixed-edge
      SC 6x6, same chi-cubic as `coupled_det`. Electrical row is still
      phi'=0 on both edges (Duan2005 Eq. 22).
    - `oc_cf_det()` / `oc_fc_det()` and bisects: mixed-edge linear OC.
      Rank-1 alpha*S update lives only on the FREE edge's M_rr row.
      n=0 stiffens vs the mixed elastic bilayer (not a C-C coincidence
      theorem: one edge has w' free, so S need not vanish). n!=0 and
      h1=0 reduce identically to the mixed elastic 4x4.

    Row formulas are the CORRECTED ones from LESSONS_LEARNED.md Sec
    18.149 (Duan2005 Eq. 9a, re-verified against a 400 DPI page-image
    render, not OCR -- the earlier Sec 18.141 C transcription had a
    dropped minus sign):
        M_rr-row(i) = [(d1+d2)*lambda_i + (4/pi)*h1*e31_bar*chi_i
                       + 2*A1*n^2/r^2] * Z - (2*A1/r) * dZ/dr
        Q_r-row(i)  = [(d1+d2)*lambda_i + (4/pi)*h1*e31_bar*chi_i] * dZ/dr
        phi'-row(i) = chi_i * dZ/dr
    reduced to Z, dZ/dr only via each branch's own defining ODE
    Delta(Z) = lambda_i*Z. The elastic-only path is this same M_rr/Q_r
    row with the piezoelectric term (chi_i, the phi' row entirely)
    dropped, which is exactly the e31_bar=0 limit -- except that limit
    is genuinely singular in the chi-cubic (Sec 18.149: leading and two
    subleading coefficients all vanish together), so `elastic_det` is a
    separate, from-scratch 2-branch construction rather than
    `coupled_det` called with e31_bar forced to zero, mirroring how both
    probe scripts keep their own elastic baseline in a separate script.
    """

    _PROJECTIONS = ('consistent', 'duan', 'galerkin_sine')

    def __init__(self, r_i, r_o, h, E, nu, rho,
                 h1=0.0,
                 C11E=None, C12E=None, C13E=None, C33E=None,
                 e31=None, e33=None, X11=None, X33=None,
                 rho_pzt=None, dps=100, projection='consistent'):
        """
        r_i, r_o : inner/outer radius (m).
        h        : HOST HALF-thickness (m) -- Duan2005's own symbol; full
                   host thickness is 2h.
        E, nu, rho : host material (Pa, --, kg/m^3), isotropic.
        h1       : piezo layer thickness EACH LAYER (m); 0.0 (default) is
                   the bare host, no piezo layer at all.
        C11E, C12E, C13E, C33E : the piezo layer's RAW elastic stiffness
                   at constant E-field (Pa) -- Duan2005 Table 1 sourcing
                   directly (e.g. PZT4: 132e9, 71e9, 73e9, 115e9). Needed
                   whenever h1 > 0, for both `elastic_det`/`elastic_bisect`
                   (the added-stiffness effect alone) and `coupled_det`/
                   `coupled_bisect` (the full piezo-coupled model).
        e31, e33 : piezoelectric stress constants (C/m^2). Needed ONLY
                   for `coupled_det`/`coupled_bisect` -- the elastic-only
                   path never touches them.
        X11, X33 : permittivity constants (F/m). Same scope as e31/e33.
        rho_pzt  : piezo layer density (kg/m^3); needed whenever h1 > 0.
        projection : 'consistent' (default; Galerkin quadratic bubble +
                   reciprocal OC charge arm h+h1/2), 'duan' (legacy,
                   Duan2005-faithful, bit-identical to the pre-2026-09-23
                   code), or 'galerkin_sine' (diagnostic). Sec 18.231.
        dps      : mpmath working precision for methods on this instance
                   (applied via an `mp.workdps` context -- does not leak
                   into global mpmath state, unlike the standalone probe
                   scripts' `mp.dps = ...` top-level assignment). Default
                   100 matches both probe scripts' own working precision
                   for the coupled model and the C-C production cluster
                   gate (Sec 18.148, job 2489140).
        """
        self.r_i = float(r_i)
        self.r_o = float(r_o)
        self.h = float(h)
        self.E = float(E)
        self.nu = float(nu)
        self.rho = float(rho)
        self.h1 = float(h1)
        self.C11E = None if C11E is None else float(C11E)
        self.C12E = None if C12E is None else float(C12E)
        self.C13E = None if C13E is None else float(C13E)
        self.C33E = None if C33E is None else float(C33E)
        self.e31 = None if e31 is None else float(e31)
        self.e33 = None if e33 is None else float(e33)
        self.X11 = None if X11 is None else float(X11)
        self.X33 = None if X33 is None else float(X33)
        self.rho_pzt = None if rho_pzt is None else float(rho_pzt)
        self.dps = int(dps)
        if projection not in self._PROJECTIONS:
            raise ValueError(
                "projection must be one of %r, got %r"
                % (self._PROJECTIONS, projection))
        self.projection = projection
        if self.h1 != 0.0:
            missing = [name for name, v in (
                ("C11E", self.C11E), ("C12E", self.C12E),
                ("C13E", self.C13E), ("C33E", self.C33E),
                ("rho_pzt", self.rho_pzt)) if v is None]
            if missing:
                raise ValueError(
                    "h1 > 0 (a piezo layer is present) needs %s (the "
                    "elastic-bilayer added-stiffness constants); for "
                    "the bare host case pass h1=0.0 (the default) "
                    "instead" % ", ".join(missing))

    # ---------- reduced constants (Sec 18.141 A) ----------

    def _c11_bar(self):
        C11E = mpf(self.C11E); C13E = mpf(self.C13E); C33E = mpf(self.C33E)
        return C11E - C13E ** 2 / C33E

    def _c12_bar(self):
        C12E = mpf(self.C12E); C13E = mpf(self.C13E); C33E = mpf(self.C33E)
        return C12E - C13E ** 2 / C33E

    def _e31_bar(self):
        if self.e31 is None or self.e33 is None:
            raise ValueError(
                "the piezo-coupled model needs e31 and e33 (the elastic-"
                "only path -- elastic_det/elastic_bisect -- does not)")
        e31 = mpf(self.e31); e33 = mpf(self.e33)
        C13E = mpf(self.C13E); C33E = mpf(self.C33E)
        return e31 - (C13E / C33E) * e33

    def _Xi33_bar(self):
        if self.X33 is None or self.e33 is None:
            raise ValueError(
                "the piezo-coupled model needs X33 and e33 (the elastic-"
                "only path does not)")
        X33 = mpf(self.X33); e33 = mpf(self.e33); C33E = mpf(self.C33E)
        return X33 + e33 ** 2 / C33E

    def _Xi11_bar(self):
        if self.X11 is None:
            raise ValueError(
                "the piezo-coupled model needs X11 (the elastic-only "
                "path does not)")
        return mpf(self.X11)

    # ---------- elastic building blocks (Sec 18.149 / the elastic-
    # baseline probe script's det_at(), byte-for-byte) ----------

    def _d1(self):
        """Host bending stiffness d1 = (2/3) E h^3 / (1-nu^2). Present even
        at h1=0 -- this alone is the bare-host case."""
        E = mpf(self.E); h = mpf(self.h); nu = mpf(self.nu)
        return mpf(2) / 3 * E * h ** 3 / (1 - nu ** 2)

    def _d2(self, h1):
        """Piezo layer's own added elastic bending stiffness. Zero at
        h1=0 by construction (no layer), matching the bare-host case
        exactly without needing any piezo material constants at all."""
        h1 = mpf(h1)
        if h1 == 0:
            return mpf(0)
        h = mpf(self.h)
        return mpf(2) / 3 * self._c11_bar() * ((h + h1) ** 3 - h ** 3)

    def _A1(self, d2, h1):
        """A1 = 0.5*[(1-nu)*d1 + (1 - c12_bar/c11_bar)*d2] (Duan2005's own
        A1 definition). At h1=0, d2=0 so the second term is dropped
        entirely rather than needing c11_bar/c12_bar at all."""
        d1 = self._d1()
        nu = mpf(self.nu)
        if mpf(h1) == 0:
            return mpf('0.5') * (1 - nu) * d1
        c11_bar = self._c11_bar(); c12_bar = self._c12_bar()
        return mpf('0.5') * ((1 - nu) * d1 + (1 - c12_bar / c11_bar) * d2)

    def _A2(self, h1):
        """A2 = 2*(rho_host*h + rho_pzt*h1) -- the areal inertia term."""
        h1 = mpf(h1)
        base = 2 * mpf(self.rho) * mpf(self.h)
        if h1 == 0:
            return base
        return base + 2 * mpf(self.rho_pzt) * h1

    def elastic_det(self, omega, n, h1=None):
        """4x4 elastic F-F determinant (2 branches, M_rr=Q_r=0 at r_i, r_o).

        h1=None (default) uses the instance's own h1; pass h1=0.0
        explicitly to force the bare-host case regardless of what the
        instance was constructed with. Matches
        probe_piezo_p4_elastic_ff_baseline_2026-09-15.py's det_at()
        exactly: same row formulas, same 4-column order
        [I_n, K_n, J_n, Y_n], same per-row equilibration before mp.det.
        """
        if h1 is None:
            h1 = self.h1
        with mp.workdps(self.dps):
            omega = mpf(omega)
            d1 = self._d1()
            d2 = self._d2(h1)
            dsum = d1 + d2
            A1 = self._A1(d2, h1)
            A2 = self._A2(h1)
            mu = mp.sqrt(A2 * omega ** 2 / dsum)
            k = mp.sqrt(mu)

            def rows(r):
                r = mpf(r)
                x = k * r
                I = mp.besseli(n, x); K = mp.besselk(n, x)
                J = mp.besselj(n, x); Y = mp.bessely(n, x)
                dI = k * mpf('0.5') * (mp.besseli(n - 1, x) + mp.besseli(n + 1, x))
                dK = k * mpf('-0.5') * (mp.besselk(n - 1, x) + mp.besselk(n + 1, x))
                dJ = k * mpf('0.5') * (mp.besselj(n - 1, x) - mp.besselj(n + 1, x))
                dY = k * mpf('0.5') * (mp.bessely(n - 1, x) - mp.bessely(n + 1, x))
                cols = [(I, dI, mu), (K, dK, mu), (J, dJ, -mu), (Y, dY, -mu)]
                Mrow, Qrow = [], []
                for Z, dZ, lam in cols:
                    Mrow.append((dsum * lam + 2 * A1 * n ** 2 / r ** 2) * Z
                                - (2 * A1 / r) * dZ)
                    # Kirchhoff effective shear V_r = Q_r + (1/r) dM_rtheta/dtheta
                    # (2026-09-22 fix, LESSONS Sec 18.226): twisting term added;
                    # vanishes identically at n=0, so n=0 results are unchanged.
                    Qrow.append(dsum * lam * dZ
                                - 2 * A1 * n ** 2 * (dZ / r ** 2 - Z / r ** 3))
                return Mrow, Qrow

            M_ri, Q_ri = rows(self.r_i)
            M_ro, Q_ro = rows(self.r_o)
            M = matrix([M_ri, Q_ri, M_ro, Q_ro])
            for i in range(4):
                s = max(abs(M[i, j]) for j in range(4)) or mpf(1)
                for j in range(4):
                    M[i, j] = M[i, j] / s
            return mp.det(M)

    def elastic_bisect(self, lo, hi, n, h1=None, iters=50):
        """Bisect the elastic F-F determinant's real-part sign flip on
        [lo, hi] (rad/s). Returns a plain float. Matches the sibling
        probe script's bisect() exactly (same iteration count default,
        same sign-tracking logic)."""
        if h1 is None:
            h1 = self.h1
        with mp.workdps(self.dps):
            lo = mpf(lo); hi = mpf(hi)
            flo = self.elastic_det(lo, n, h1=h1)
            for _ in range(iters):
                mid = (lo + hi) / 2
                fm = self.elastic_det(mid, n, h1=h1)
                if (fm.real > 0) == (flo.real > 0):
                    lo, flo = mid, fm
                else:
                    hi = mid
            return float((lo + hi) / 2)

    # ---------- electrical projection + OC charge arm (Sec 18.231) ----------
    # (A) (d1+d2) Lap^2 w + K_pref Lap phibar - A2 omega^2 w = 0
    # (B) a Lap phibar - phibar + c Lap w = 0   =>  lam = chi/(a chi + c)
    # 'duan' keeps its original code paths so it stays bit-identical.

    def _proj_skin_integrals(self, h1):
        """Per-skin (I0, I1, I2, I3) = int f, int z f', int f'^2, int f^2
        over one skin s in [0, h1] (derivation Sec 3.2)."""
        h1 = mpf(h1)
        pi = mp.pi
        if self.projection == 'consistent':
            return 2 * h1 / 3, -2 * h1 / 3, mpf(16) / (3 * h1), 8 * h1 / 15
        if self.projection == 'galerkin_sine':
            return 2 * h1 / pi, -2 * h1 / pi, pi ** 2 / (2 * h1), h1 / 2
        raise ValueError("'duan' is not a Galerkin projection")

    def _K_pref(self, h1):
        """Piezo moment coefficient (both skins): M_rr contains
        -K_pref*phibar."""
        h1 = mpf(h1)
        if self.projection == 'duan':
            return mpf(4) / mp.pi * h1 * self._e31_bar()
        _, I1, _, _ = self._proj_skin_integrals(h1)
        return -2 * I1 * self._e31_bar()

    def _elec_ac(self, h1):
        """(a, c) of the normalised per-skin electrical equation (B)."""
        h1 = mpf(h1)
        Xi33 = self._Xi33_bar(); Xi11 = self._Xi11_bar()
        e31_bar = self._e31_bar()
        if self.projection == 'duan':
            return (Xi11 * h1 ** 2 / (mp.pi ** 2 * Xi33),
                    e31_bar * h1 ** 2 / (2 * mp.pi * Xi33))
        _, I1, I2, I3 = self._proj_skin_integrals(h1)
        return Xi11 * I3 / (Xi33 * I2), -I1 * e31_bar / (Xi33 * I2)

    def _charge_arm(self, h1):
        """Motional bus-charge arm: Q = C0 V + 4 pi e31bar * arm * S.
        'duan': h+h1 (pointwise outer-face D_z of the linear field);
        otherwise h+h1/2, the variational (reciprocal) value = the moment
        arm of Mp = 2 e31bar (h+h1/2) V."""
        h = mpf(self.h); h1 = mpf(h1)
        if self.projection == 'duan':
            return h + h1
        return h + h1 / 2

    def _phi_charge_weight(self, h1):
        """Coefficient w_phi of the interior-field bus-charge term
        -w_phi*S_phi (lever 12). 'duan': 4 pi Xi33bar/h1 (Rayleigh-sign
        fixed, Sec 18.182); Galerkin: 0 exactly (int f' dz = 0)."""
        h1 = mpf(h1)
        if self.projection == 'duan':
            return 4 * mp.pi * self._Xi33_bar() / h1
        return mpf(0)

    def long_wave_sc_ratio(self, omega, h1=None):
        """(D_eff - d1 - d2)/[2 * h1^3 e31bar^2/(12 Xi33bar)] from the
        flexural branch (smallest |lam|), D_eff = A2 omega^2/lam^2.
        -> 1 / 96/pi^4 / 12/pi^2 as omega -> 0. Gate P2, Sec 18.231."""
        if h1 is None:
            h1 = self.h1
        with mp.workdps(self.dps):
            omega = mpf(omega); h1 = mpf(h1)
            br = self._branches(omega, h1)
            chi, lam = min(br, key=lambda cl: abs(cl[1]))
            Deff = self._A2(h1) * omega ** 2 / lam ** 2
            exact = 2 * h1 ** 3 * self._e31_bar() ** 2 / (12 * self._Xi33_bar())
            return (Deff.real - self._d1() - self._d2(h1)) / exact

    # ---------- piezo-coupled model (Sec 18.141 D/G, Sec 18.149) ----------

    def _cubic_coeffs(self, omega, h1, d2):
        """chi-cubic coefficients (Duan2005 Eq. 17a, h1^4-corrected per
        Sec 18.147). Matches the coupled-model probe script's
        cubic_coeffs() exactly -- does not depend on which mechanical
        edge condition is applied (only the boundary rows in
        `coupled_det` differ between F-F and, if ever ported, C-C).
        Non-'duan' projections: Galerkin form (Sec 18.231), same roots
        up to an overall scale when the coefficients coincide."""
        if self.projection != 'duan':
            omega = mpf(omega); h1 = mpf(h1)
            w2 = omega ** 2
            A2 = self._A2(h1)
            dsum = self._d1() + d2
            K = self._K_pref(h1)
            a, c = self._elec_ac(h1)
            return (K * a,
                    dsum + K * c - A2 * w2 * a ** 2,
                    -2 * A2 * w2 * a * c,
                    -A2 * w2 * c ** 2)
        pi = mp.pi
        omega = mpf(omega); h1 = mpf(h1)
        w2 = omega ** 2
        A2 = self._A2(h1)
        Xi11_bar = self._Xi11_bar(); Xi33_bar = self._Xi33_bar()
        e31_bar = self._e31_bar()
        d1 = self._d1()
        a = -16 * pi * Xi11_bar * Xi33_bar * e31_bar * h1 ** 3
        b = (4 * A2 * w2 * h1 ** 4 * Xi11_bar ** 2
             - 8 * pi ** 2 * Xi33_bar * e31_bar ** 2 * h1 ** 3
             - 4 * pi ** 4 * Xi33_bar ** 2 * (d1 + d2))
        c = 4 * pi * A2 * w2 * h1 ** 4 * Xi11_bar * e31_bar
        d = pi ** 2 * A2 * w2 * h1 ** 4 * e31_bar ** 2
        return a, b, c, d

    def _branches(self, omega, h1):
        """Solve the chi-cubic for the three (chi_i, lambda_i) branch
        pairs. Matches the probe script's branches() exactly."""
        pi = mp.pi
        h1 = mpf(h1)
        d2 = self._d2(h1)
        a, b, c, d = self._cubic_coeffs(omega, h1, d2)
        roots = mp.polyroots([a, b, c, d], maxsteps=300, extraprec=1500)
        if self.projection != 'duan':
            ea, ec = self._elec_ac(h1)
            return [(chi, chi / (ea * chi + ec)) for chi in roots]
        Xi33_bar = self._Xi33_bar(); Xi11_bar = self._Xi11_bar()
        e31_bar = self._e31_bar()
        out = []
        for chi in roots:
            lam = (2 * pi ** 2 * Xi33_bar * chi
                   / (h1 ** 2 * (2 * Xi11_bar * chi + e31_bar * pi)))
            out.append((chi, lam))
        return out

    @staticmethod
    def _radial_quad(n, r, lam):
        """Z1, Z2, dZ1, dZ2 for Delta w = lam*w -- the ordinary (lam<0)
        or modified (lam>=0) Bessel-family radial solution pair. Matches
        the probe script's radial_quad() exactly."""
        r = mpf(r)
        if lam.real >= 0:
            delta = mp.sqrt(lam); x = delta * r
            Z1 = mp.besseli(n, x); Z2 = mp.besselk(n, x)
            dZ1 = delta * mpf('0.5') * (mp.besseli(n - 1, x) + mp.besseli(n + 1, x))
            dZ2 = delta * mpf('-0.5') * (mp.besselk(n - 1, x) + mp.besselk(n + 1, x))
        else:
            delta = mp.sqrt(-lam); x = delta * r
            Z1 = mp.besselj(n, x); Z2 = mp.bessely(n, x)
            dZ1 = delta * mpf('0.5') * (mp.besselj(n - 1, x) - mp.besselj(n + 1, x))
            dZ2 = delta * mpf('0.5') * (mp.bessely(n - 1, x) - mp.bessely(n + 1, x))
        return Z1, Z2, dZ1, dZ2

    def coupled_det(self, omega, n, h1=None):
        """6x6 F-F piezo-coupled boundary determinant: per edge, M_rr=0,
        Q_r=0, phi'=0 (Duan2005 Eq. 25a,b; corrected M_rr/Q_r rows per
        Sec 18.149; the phi'=0 electrical row is Duan2005's own Eq. 22,
        confirmed identical for C, S, and F mechanical edge conditions).
        Matches the coupled-model probe script's ff_det() exactly.

        h1=None (default) uses the instance's own h1. h1=0 is refused --
        the bare-host/zero-coupling limit is a genuine singularity of
        the chi-cubic (Sec 18.149: leading coefficient `a` and subleading
        `c`, `d` are all proportional to e31_bar and vanish together,
        degenerating mp.polyroots's cubic solve), NOT a degenerate case
        this method can fall back to safely -- use `elastic_det(h1=0.0)`
        for that limit instead, exactly as both probe scripts keep the
        elastic baseline in a separate construction rather than trying
        to zero out e31_bar inside the coupled model.
        """
        if h1 is None:
            h1 = self.h1
        if float(h1) == 0.0:
            raise ValueError(
                "coupled_det needs h1 > 0 (a piezo layer present) -- "
                "h1=0 is a genuine singular limit of the chi-cubic "
                "(Sec 18.149), not a valid input here; use "
                "elastic_det(h1=0.0) for the bare-host limit instead")
        with mp.workdps(self.dps):
            omega = mpf(omega)
            h1 = mpf(h1)
            d1 = self._d1()
            d2 = self._d2(h1)
            lams = self._branches(omega, h1)
            A1v = self._A1(d2, h1)
            e31_bar = self._e31_bar()
            K_pref = self._K_pref(h1)

            def rows_at(r):
                r = mpf(r)
                m_row, q_row, phip_row = [], [], []
                for chi, lam in lams:
                    Z1, Z2, dZ1, dZ2 = self._radial_quad(n, r, lam)
                    Ki = (d1 + d2) * lam + K_pref * chi
                    for Z, dZ in ((Z1, dZ1), (Z2, dZ2)):
                        m_row.append((Ki + 2 * A1v * n ** 2 / r ** 2) * Z
                                    - (2 * A1v / r) * dZ)
                        # Kirchhoff effective shear V_r = Q_r + (1/r) dM_rtheta/dtheta
                        # (2026-09-22 fix, LESSONS Sec 18.226): twisting term added;
                        # vanishes identically at n=0, so n=0 results are unchanged.
                        q_row.append(Ki * dZ
                                     - 2 * A1v * n ** 2 * (dZ / r ** 2 - Z / r ** 3))
                        phip_row.append(chi * dZ)
                return m_row, q_row, phip_row

            m_ri, q_ri, phip_ri = rows_at(self.r_i)
            m_ro, q_ro, phip_ro = rows_at(self.r_o)
            M = matrix([m_ri, q_ri, phip_ri, m_ro, q_ro, phip_ro])
            for j in range(6):
                mx = max(abs(M[i, j]) for i in range(6)) or mpf(1)
                for i in range(6):
                    M[i, j] = M[i, j] / mx
            for i in range(6):
                mx = max(abs(M[i, j]) for j in range(6)) or mpf(1)
                for j in range(6):
                    M[i, j] = M[i, j] / mx
            return mp.det(M)

    def coupled_bisect(self, lo, hi, n, h1=None, iters=45):
        """Bisect the piezo-coupled F-F determinant's real-part sign flip
        on [lo, hi] (rad/s). Returns a plain float. Matches the
        coupled-model probe script's bisect() exactly (same default
        iteration count, same sign-tracking logic)."""
        if h1 is None:
            h1 = self.h1
        with mp.workdps(self.dps):
            lo = mpf(lo); hi = mpf(hi)
            flo = self.coupled_det(lo, n, h1=h1)
            for _ in range(iters):
                mid = (lo + hi) / 2
                fm = self.coupled_det(mid, n, h1=h1)
                if (fm.real > 0) == (flo.real > 0):
                    lo, flo = mid, fm
                else:
                    hi = mid
            return float((lo + hi) / 2)

    # ---------- C-C (clamped-clamped) counterparts (Sec 18.173) ----------

    def elastic_cc_det(self, omega, n, h1=None):
        """4x4 elastic C-C determinant (2 branches, w=w'=0 at r_i, r_o).
        Same d1/d2/A2 building blocks as elastic_det -- no A1 needed at
        all, since the C-C boundary never touches M_rr. Matches
        probe_piezo_p4_elastic_cc_baseline_2026-09-15.py's det_at()
        exactly (same 4-column [I_n, K_n, J_n, Y_n] order, same per-row
        equilibration before mp.det)."""
        if h1 is None:
            h1 = self.h1
        with mp.workdps(self.dps):
            omega = mpf(omega)
            d1 = self._d1()
            d2 = self._d2(h1)
            dsum = d1 + d2
            A2 = self._A2(h1)
            mu = mp.sqrt(A2 * omega ** 2 / dsum)
            k = mp.sqrt(mu)

            def row(r):
                r = mpf(r)
                x = k * r
                I = mp.besseli(n, x); K = mp.besselk(n, x)
                J = mp.besselj(n, x); Y = mp.bessely(n, x)
                dI = k * mpf('0.5') * (mp.besseli(n - 1, x) + mp.besseli(n + 1, x))
                dK = k * mpf('-0.5') * (mp.besselk(n - 1, x) + mp.besselk(n + 1, x))
                dJ = k * mpf('0.5') * (mp.besselj(n - 1, x) - mp.besselj(n + 1, x))
                dY = k * mpf('0.5') * (mp.bessely(n - 1, x) - mp.bessely(n + 1, x))
                return [I, K, J, Y], [dI, dK, dJ, dY]

            w_ri, wp_ri = row(self.r_i)
            w_ro, wp_ro = row(self.r_o)
            M = matrix([w_ri, wp_ri, w_ro, wp_ro])
            for i in range(4):
                s = max(abs(M[i, j]) for j in range(4)) or mpf(1)
                for j in range(4):
                    M[i, j] = M[i, j] / s
            return mp.det(M)

    def elastic_cc_bisect(self, lo, hi, n, h1=None, iters=50):
        """Bisect the elastic C-C determinant's real-part sign flip on
        [lo, hi] (rad/s). Returns a plain float. Matches the sibling
        probe script's bisect() exactly."""
        if h1 is None:
            h1 = self.h1
        with mp.workdps(self.dps):
            lo = mpf(lo); hi = mpf(hi)
            flo = self.elastic_cc_det(lo, n, h1=h1)
            for _ in range(iters):
                mid = (lo + hi) / 2
                fm = self.elastic_cc_det(mid, n, h1=h1)
                if (fm.real > 0) == (flo.real > 0):
                    lo, flo = mid, fm
                else:
                    hi = mid
            return float((lo + hi) / 2)

    def cc_coupled_det(self, omega, n, h1=None):
        """6x6 C-C piezo-coupled boundary determinant: per edge, w=0,
        w'=0, phi'=0 (Duan2005 Eq. 23a,b -- confirmed identical electrical
        row to F-F, Sec 18.141 F/G). Reuses `_branches`/`_radial_quad`
        unchanged (the chi-cubic does not depend on the mechanical edge
        condition) -- only the boundary-row construction differs from
        `coupled_det`. Matches
        probe_piezo_p4_cc_forward_model_2026-09-15.py's cc_det() exactly.

        h1=None (default) uses the instance's own h1. h1=0 is refused for
        the same reason as coupled_det: the bare-host limit is a genuine
        chi-cubic singularity (Sec 18.149), not a valid input here -- use
        elastic_cc_det(h1=0.0) for that limit instead.
        """
        if h1 is None:
            h1 = self.h1
        if float(h1) == 0.0:
            raise ValueError(
                "cc_coupled_det needs h1 > 0 (a piezo layer present) -- "
                "h1=0 is a genuine singular limit of the chi-cubic "
                "(Sec 18.149), not a valid input here; use "
                "elastic_cc_det(h1=0.0) for the bare-host limit instead")
        with mp.workdps(self.dps):
            omega = mpf(omega)
            h1 = mpf(h1)
            lams = self._branches(omega, h1)

            def rows_at(r):
                r = mpf(r)
                w_row, wp_row, phip_row = [], [], []
                for chi, lam in lams:
                    Z1, Z2, dZ1, dZ2 = self._radial_quad(n, r, lam)
                    w_row.append(Z1); w_row.append(Z2)
                    wp_row.append(dZ1); wp_row.append(dZ2)
                    phip_row.append(chi * dZ1); phip_row.append(chi * dZ2)
                return w_row, wp_row, phip_row

            w_ri, wp_ri, phip_ri = rows_at(self.r_i)
            w_ro, wp_ro, phip_ro = rows_at(self.r_o)
            M = matrix([w_ri, wp_ri, phip_ri, w_ro, wp_ro, phip_ro])
            for j in range(6):
                mx = max(abs(M[i, j]) for i in range(6)) or mpf(1)
                for i in range(6):
                    M[i, j] = M[i, j] / mx
            for i in range(6):
                mx = max(abs(M[i, j]) for j in range(6)) or mpf(1)
                for j in range(6):
                    M[i, j] = M[i, j] / mx
            return mp.det(M)

    def cc_coupled_bisect(self, lo, hi, n, h1=None, iters=40):
        """Bisect the piezo-coupled C-C determinant's real-part sign flip
        on [lo, hi] (rad/s). Returns a plain float. Matches the
        coupled-model C-C probe script's bisect() exactly (same default
        iteration count -- 40, not F-F's 45 -- same sign-tracking logic).
        """
        if h1 is None:
            h1 = self.h1
        with mp.workdps(self.dps):
            lo = mpf(lo); hi = mpf(hi)
            flo = self.cc_coupled_det(lo, n, h1=h1)
            for _ in range(iters):
                mid = (lo + hi) / 2
                fm = self.cc_coupled_det(mid, n, h1=h1)
                if (fm.real > 0) == (flo.real > 0):
                    lo, flo = mid, fm
                else:
                    hi = mid
            return float((lo + hi) / 2)

    # ---------- F-F open-circuit (Sec 18.175) ----------

    def oc_ff_det(self, omega, n, h1=None):
        """4x4 F-F open-circuit determinant (fully electroded, linear
        through-thickness potential, Q=0 on the outer-electrode bus).

        n!=0 or h1=0 reduce identically to elastic_det -- net charge on a
        continuous electrode vanishes for n!=0 by the theta-integral, and
        h1=0 has no piezo layer to charge. n=0 with h1>0 adds the same
        rank-1 update to BOTH M_rr rows (the piezo moment is spatially
        constant) and leaves the Q_r rows unchanged (a constant moment
        produces zero shear). Sign of the update is the one that RAISES
        the n=0 root relative to the elastic bilayer (piezoelectric
        stiffening; sandbox-confirmed, Sec 18.175).
        """
        if h1 is None:
            h1 = self.h1
        if float(h1) == 0.0 or int(n) != 0:
            return self.elastic_det(omega, n, h1=h1)
        with mp.workdps(self.dps):
            omega = mpf(omega)
            h1 = mpf(h1)
            d1 = self._d1()
            d2 = self._d2(h1)
            dsum = d1 + d2
            A1 = self._A1(d2, h1)
            A2 = self._A2(h1)
            mu = mp.sqrt(A2 * omega ** 2 / dsum)
            k = mp.sqrt(mu)
            e31b = self._e31_bar()
            Xi33b = self._Xi33_bar()
            h = mpf(self.h)
            ri = mpf(self.r_i)
            ro = mpf(self.r_o)
            alpha = (4 * e31b ** 2 * (h + h1 / 2) * self._charge_arm(h1) * h1
                     / (Xi33b * (ro ** 2 - ri ** 2)))

            def cols_at(r):
                r = mpf(r)
                x = k * r
                I = mp.besseli(n, x); K = mp.besselk(n, x)
                J = mp.besselj(n, x); Y = mp.bessely(n, x)
                dI = k * mpf('0.5') * (mp.besseli(n - 1, x) + mp.besseli(n + 1, x))
                dK = k * mpf('-0.5') * (mp.besselk(n - 1, x) + mp.besselk(n + 1, x))
                dJ = k * mpf('0.5') * (mp.besselj(n - 1, x) - mp.besselj(n + 1, x))
                dY = k * mpf('0.5') * (mp.bessely(n - 1, x) - mp.bessely(n + 1, x))
                return [(I, dI, mu), (K, dK, mu), (J, dJ, -mu), (Y, dY, -mu)]

            def m_q(r, cols):
                r = mpf(r)
                Mrow, Qrow, dZrow = [], [], []
                for Z, dZ, lam in cols:
                    Mrow.append((dsum * lam + 2 * A1 * n ** 2 / r ** 2) * Z
                                - (2 * A1 / r) * dZ)
                    # Kirchhoff effective shear V_r = Q_r + (1/r) dM_rtheta/dtheta
                    # (2026-09-22 fix, LESSONS Sec 18.226): twisting term added;
                    # vanishes identically at n=0, so n=0 results are unchanged.
                    Qrow.append(dsum * lam * dZ
                                - 2 * A1 * n ** 2 * (dZ / r ** 2 - Z / r ** 3))
                    dZrow.append(dZ)
                return Mrow, Qrow, dZrow

            cols_i = cols_at(ri)
            cols_o = cols_at(ro)
            M_ri, Q_ri, dZi = m_q(ri, cols_i)
            M_ro, Q_ro, dZo = m_q(ro, cols_o)
            S = [ro * dZo[j] - ri * dZi[j] for j in range(4)]
            for j in range(4):
                M_ri[j] = M_ri[j] + alpha * S[j]
                M_ro[j] = M_ro[j] + alpha * S[j]
            M = matrix([M_ri, Q_ri, M_ro, Q_ro])
            for i in range(4):
                s = max(abs(M[i, j]) for j in range(4)) or mpf(1)
                for j in range(4):
                    M[i, j] = M[i, j] / s
            return mp.det(M)

    def oc_ff_bisect(self, lo, hi, n, h1=None, iters=50):
        """Bisect the F-F open-circuit determinant's real-part sign flip
        on [lo, hi] (rad/s). Returns a plain float. Default iters=50
        matches elastic_bisect (this is a 4x4, not the coupled 6x6)."""
        if h1 is None:
            h1 = self.h1
        with mp.workdps(self.dps):
            lo = mpf(lo); hi = mpf(hi)
            flo = self.oc_ff_det(lo, n, h1=h1)
            for _ in range(iters):
                mid = (lo + hi) / 2
                fm = self.oc_ff_det(mid, n, h1=h1)
                if (fm.real > 0) == (flo.real > 0):
                    lo, flo = mid, fm
                else:
                    hi = mid
            return float((lo + hi) / 2)

    # ---------- F-F driven admittance (Sec 18.179) ----------

    def clamped_C0(self, h1=None):
        """Parallel-bimorph clamped capacitance, two layers:
        C0 = 2 * Xi33_bar * A / h1, A = pi (r_o^2 - r_i^2). SI farads."""
        if h1 is None:
            h1 = self.h1
        if float(h1) == 0.0:
            return 0.0
        with mp.workdps(self.dps):
            h1 = mpf(h1)
            A = mp.pi * (mpf(self.r_o) ** 2 - mpf(self.r_i) ** 2)
            return float(2 * self._Xi33_bar() * A / h1)

    def driven_ff_qv(self, omega, n, volt=1.0, h1=None):
        """Imposed-V driven F-F problem (linear-potential model).

        Returns a dict of plain floats:
          q, q_over_v, c0, s, mp
        where q is the outer-bus charge (C), s = r_o w'(r_o) - r_i w'(r_i)
        (m^2/m = m, since w is metres), mp is the constant piezo moment.

        n!=0: motional charge vanishes (theta-integral), q = C0 * volt.
        h1=0: q = 0.
        n=0, h1>0: solve the elastic F-F 4x4 with
          M_rr-row = 2*e31_bar*(h+h1/2)*volt  at BOTH edges
          Q_r-row  = 0
        then q = 4*pi*e31_bar*(h+h1)*S + C0*volt
        (LESSONS_LEARNED.md Sec 18.175 charge formula). Y = j*omega*q/volt
        is the undamped admittance. Poles of q/volt sit at the elastic
        F-F (short-circuit) roots; zeros sit at the OC roots.
        """
        if h1 is None:
            h1 = self.h1
        volt = float(volt)
        c0 = self.clamped_C0(h1=h1)
        if float(h1) == 0.0 or volt == 0.0:
            return dict(q=0.0, q_over_v=0.0, c0=c0, s=0.0, mp=0.0)
        if int(n) != 0:
            return dict(q=c0 * volt, q_over_v=c0, c0=c0, s=0.0, mp=0.0)
        with mp.workdps(self.dps):
            omega = mpf(omega)
            h1 = mpf(h1)
            V = mpf(volt)
            d1 = self._d1()
            d2 = self._d2(h1)
            dsum = d1 + d2
            A1 = self._A1(d2, h1)
            A2 = self._A2(h1)
            mu = mp.sqrt(A2 * omega ** 2 / dsum)
            k = mp.sqrt(mu)
            e31b = self._e31_bar()
            Xi33b = self._Xi33_bar()
            h = mpf(self.h)
            ri = mpf(self.r_i)
            ro = mpf(self.r_o)
            n = int(n)
            Mp = 2 * e31b * (h + h1 / 2) * V
            C0 = 2 * Xi33b * mp.pi * (ro ** 2 - ri ** 2) / h1

            def cols_at(r):
                r = mpf(r)
                x = k * r
                I = mp.besseli(n, x); K = mp.besselk(n, x)
                J = mp.besselj(n, x); Y = mp.bessely(n, x)
                dI = k * mpf('0.5') * (mp.besseli(n - 1, x) + mp.besseli(n + 1, x))
                dK = k * mpf('-0.5') * (mp.besselk(n - 1, x) + mp.besselk(n + 1, x))
                dJ = k * mpf('0.5') * (mp.besselj(n - 1, x) - mp.besselj(n + 1, x))
                dY = k * mpf('0.5') * (mp.bessely(n - 1, x) - mp.bessely(n + 1, x))
                return [(I, dI, mu), (K, dK, mu), (J, dJ, -mu), (Y, dY, -mu)]

            def m_q(r, cols):
                r = mpf(r)
                Mrow, Qrow, dZrow = [], [], []
                for Z, dZ, lam in cols:
                    Mrow.append((dsum * lam + 2 * A1 * n ** 2 / r ** 2) * Z
                                - (2 * A1 / r) * dZ)
                    # Kirchhoff effective shear V_r = Q_r + (1/r) dM_rtheta/dtheta
                    # (2026-09-22 fix, LESSONS Sec 18.226): twisting term added;
                    # vanishes identically at n=0, so n=0 results are unchanged.
                    Qrow.append(dsum * lam * dZ
                                - 2 * A1 * n ** 2 * (dZ / r ** 2 - Z / r ** 3))
                    dZrow.append(dZ)
                return Mrow, Qrow, dZrow

            cols_i = cols_at(ri)
            cols_o = cols_at(ro)
            M_ri, Q_ri, dZi = m_q(ri, cols_i)
            M_ro, Q_ro, dZo = m_q(ro, cols_o)
            A = matrix([M_ri, Q_ri, M_ro, Q_ro])
            b = matrix([Mp, mpf(0), Mp, mpf(0)])
            try:
                c = mp.lu_solve(A, b)
            except Exception:
                return dict(q=float('inf'), q_over_v=float('inf'),
                            c0=float(C0), s=float('nan'), mp=float(Mp))
            S = ro * (c[0] * dZo[0] + c[1] * dZo[1]
                      + c[2] * dZo[2] + c[3] * dZo[3]) \
                - ri * (c[0] * dZi[0] + c[1] * dZi[1]
                        + c[2] * dZi[2] + c[3] * dZi[3])
            q = 4 * mp.pi * e31b * self._charge_arm(h1) * S + C0 * V
            qv = q / V
            return dict(q=float(q), q_over_v=float(qv), c0=float(C0),
                        s=float(S), mp=float(Mp))

    # ---------- F-F "same-ansatz" open-circuit / driven (LESSONS_
    # LEARNED.md Sec 18.182, lever 12): unifies the linear bus-voltage
    # particular solution (oc_ff_det/driven_ff_qv's own ansatz) with
    # the interior sine homogeneous solution (coupled_det's own chi/
    # lambda branches) into ONE potential expansion, so SC and OC come
    # from the same operator family. See
    # PAPER4_OC_SINE_LINEAR_DERIVATION.md for the full derivation;
    # oc_ff_det/oc_ff_bisect/driven_ff_qv (the linear-only, leading-
    # order models) are kept UNCHANGED alongside these -- this is an
    # addition, not a replacement. ----------

    def oc_ff_sine_det(self, omega, n, h1=None):
        """7x7 F-F open-circuit determinant on the unified (linear +
        interior-sine) potential ansatz.

        phi_top(r,theta,z,t) = V(t)*(z-h)/h1
                                + phi_bar(r,theta,t)*sin(pi(z-h)/h1)
        (mirrored on the bottom layer, Duan2005's own convention:
        -V(t)*(z+h)/h1 + phi_bar*sin(pi(z+h)/h1) -- both outer
        electrodes wired in PARALLEL to the same bus V, matching
        clamped_C0's own "parallel-bimorph" two-layer capacitance).
        Substituting this into Duan2005's own moment/shear/Maxwell-
        integral derivation (the same one coupled_det already
        implements for the pure-sine, V=0 short-circuit case) shows:

        (1) phi_bar obeys EXACTLY the same chi-cubic/3-branch system as
            coupled_det -- V's field is spatially uniform (no r,theta
            dependence, zero z-divergence: a linear-in-z potential has
            no through-thickness Gauss's-law source), so it drops out
            of every interior PDE identically and never perturbs the
            branches, radial quadrature, or the M_rr/Q_r/phi'=0 row
            *shapes* coupled_det already has right.
        (2) V couples to the problem in exactly two places: a rigid
            additive term on BOTH M_rr boundary rows, and the OC self-
            consistency circuit condition Q=0. Written as a homogeneous
            row (M_rr_row(c) + Vcol_M*V = 0), reciprocity with the
            ALREADY-validated driven_ff_qv/oc_ff_det linear model (whose
            inhomogeneous form is M_rr_row(c) = Mp,
            Mp = 2*e31_bar*(h+h1/2)*V) fixes
                Vcol_M = -e31_bar*(2h+h1) = -Mp/V
            (the linear-only model's own Mp/V ratio, sign-flipped for
            the homogeneous-row convention -- get this backwards and
            the n=0 OC root comes out BELOW the elastic bilayer, which
            violates the piezoelectric c^D >= c^E stiffening theorem
            that must hold for ANY correctly-posed open-circuit model;
            this was caught exactly that way in-session, Sec 18.182).
        (3) The linear-only model's charge formula (oc_ff_det's implicit
            one, driven_ff_qv's explicit q = 4*pi*e31_bar*(h+h1)*S +
            C0*V) is INCOMPLETE once phi_bar != 0: it omits phi_bar's
            own contribution to the outer-electrode charge, i.e. the
            area integral of phi_bar itself (not just its radial
            derivative). For n=0, each branch's own governing ODE
            (Delta Z = lam*Z, n=0: Z'' + Z'/r = lam*Z) gives
            d/dr(r*Z') = r*lam*Z, so
                int_ri^ro r*phi_bar dr
                    = sum_i c_i*(chi_i/lam_i)*[ro*dZ_i(ro)-ri*dZ_i(ri)]
            (phi_bar_i = chi_i*w_i, Duan's own per-branch scaling,
            already used by the phi'=0 row = chi_i*dZ; this closed form
            was verified against direct mpmath numerical quadrature to
            ~50 digits on an arbitrary synthetic branch set, Sec
            18.182 -- it is a plain Bessel-calculus identity,
            independent of every electrostatics sign convention below).
            The enriched charge balance is
                Q = 4*pi*e31_bar*(h+h1)*S_w + C0*V
                    - 4*pi*Xi33_bar/h1 * S_phi
            with S_w = sum_i c_i*[ro*dZ_i(ro)-ri*dZ_i(ri)] (S,
            generalized from oc_ff_det's 4-branch elastic basis to this
            6-branch sine one) and S_phi the boundary-quadrature sum
            above. The sign of the new Xi33_bar term was fixed by the
            Rayleigh-Ritz convergence-from-above argument (Sec 18.182):
            phi_bar=0 (linear-only) is a proper subspace of this
            richer ansatz, and the n=0 OC root is a natural frequency
            of a conservative (undamped) electromechanical system, so
            enriching the electrical trial space can only move the
            fundamental-mode estimate DOWN, never up, relative to the
            restricted (linear-only) ansatz -- the minus sign gives
            750.1950... rad/s (below the 750.2394... linear-only
            anchor, h1/2h=1/12); the plus sign gives 750.2968... rad/s
            (above it) and was rejected on this basis. This sign is
            NOT independently cross-checked by a from-scratch
            electrostatics (D_z outward-normal) derivation -- two by-
            hand attempts at that gave inconsistent relative signs
            between this term and the (separately, reciprocity-
            confirmed) Vcol_M/circuit_pref pair, so the Rayleigh-Ritz
            argument plus the self-consistency check in (4) below is
            this coefficient's actual evidence, not hand electrostatics.
            Both magnitude and (this) sign only move the root by
            ~6e-5 relative (Sec 18.182) -- far below anything that
            affects a reported digit anywhere else in this project.
        (4) Internal self-consistency (necessary, not sufficient, but a
            real check on the wiring): oc_ff_sine_det's own root, fed
            into driven_ff_qv_sine, gives q_over_v/C0 ~1e-13 (dps=100,
            all three Duan2005 thickness ratios, Sec 18.182) -- i.e.
            the homogeneous 7x7's zero and the inhomogeneous-6x6-plus-
            charge-formula's Q=0 crossing land on the same frequency,
            as they must algebraically if (2) and (3) are wired
            consistently between the two methods (get either wrong in
            only one of the two methods and this check fails at
            O(0.01)-O(1) relative, not 1e-13 -- this is how the Vcol_M
            bug above was actually caught).
        (5) C-C fully-electroded OC=SC (Sec 18.175's theorem) survives
            unchanged: S_phi, evaluated on a null vector of the
            UNMODIFIED cc_coupled_det's own 6x6 matrix at one of its
            cluster-confirmed roots (2791.7866198564 rad/s, h1/2h=
            1/12, Sec 18.148), vanishes to ~150 digits at dps=150 --
            same order of exactness as S_w's algebraic zero (w'=0 is a
            literal boundary ROW at both C-C edges). This was NOT
            obvious a priori (S_phi is a differently-WEIGHTED sum
            across the 3 branches than S_w, chi_i/lam_i varies by
            branch, so w'(edge)=0 does not trivially force it to
            vanish termwise) -- it is an empirical confirmation, not
            yet a from-scratch proof; a full proof would likely follow
            from the divergence theorem on the whole piezo-layer volume
            given phi'=0 at both radial edges, but was not needed to
            settle the pre-registered question.

        Unknowns: the same 6 sine-branch coefficients coupled_det solves
        for, PLUS V (7 unknowns). Equations: the same 6 boundary rows
        (M_rr, Q_r, phi'=0 at each edge) PLUS Q=0 (7 equations).
        Nontrivial solutions require this 7x7 determinant to vanish.

        n!=0 or h1=0 reduce identically to coupled_det/elastic_det --
        the spatially uniform bus voltage V cannot couple to a non-
        axisymmetric mechanical mode (no theta-dependence to drive it),
        and there is no layer to charge at h1=0 -- mirroring
        oc_ff_det's own reduction rule exactly (verified bit-for-bit,
        Sec 18.182).
        """
        if h1 is None:
            h1 = self.h1
        if float(h1) == 0.0:
            return self.elastic_det(omega, n, h1=h1)
        if int(n) != 0:
            return self.coupled_det(omega, n, h1=h1)
        with mp.workdps(self.dps):
            omega = mpf(omega)
            h1 = mpf(h1)
            d1 = self._d1()
            d2 = self._d2(h1)
            lams = self._branches(omega, h1)
            A1v = self._A1(d2, h1)
            e31_bar = self._e31_bar()
            Xi33_bar = self._Xi33_bar()
            K_pref = self._K_pref(h1)
            h = mpf(self.h)
            ri = mpf(self.r_i)
            ro = mpf(self.r_o)
            n = int(n)

            Vcol_M = -e31_bar * (2 * h + h1)
            C0 = 2 * Xi33_bar * mp.pi * (ro ** 2 - ri ** 2) / h1
            circuit_pref = 4 * mp.pi * e31_bar * self._charge_arm(h1)
            phi_pref = -self._phi_charge_weight(h1)

            def rows_at(r):
                r = mpf(r)
                m_row, q_row, phip_row, bterm_row, col_lam = [], [], [], [], []
                for chi, lam in lams:
                    Z1, Z2, dZ1, dZ2 = self._radial_quad(n, r, lam)
                    Ki = (d1 + d2) * lam + K_pref * chi
                    for Z, dZ in ((Z1, dZ1), (Z2, dZ2)):
                        m_row.append((Ki + 2 * A1v * n ** 2 / r ** 2) * Z
                                    - (2 * A1v / r) * dZ)
                        # Kirchhoff effective shear V_r = Q_r + (1/r) dM_rtheta/dtheta
                        # (2026-09-22 fix, LESSONS Sec 18.226): twisting term added;
                        # vanishes identically at n=0, so n=0 results are unchanged.
                        q_row.append(Ki * dZ
                                     - 2 * A1v * n ** 2 * (dZ / r ** 2 - Z / r ** 3))
                        phip_row.append(chi * dZ)
                        bterm_row.append(dZ)
                        col_lam.append((chi, lam))
                return m_row, q_row, phip_row, bterm_row, col_lam

            m_ri, q_ri, phip_ri, b_ri, col_lam = rows_at(ri)
            m_ro, q_ro, phip_ro, b_ro, _ = rows_at(ro)

            bterm = [ro * b_ro[j] - ri * b_ri[j] for j in range(6)]
            circuit_row = []
            for j in range(6):
                chi_j, lam_j = col_lam[j]
                circuit_row.append(
                    bterm[j] * (circuit_pref + phi_pref * chi_j / lam_j))
            circuit_row.append(C0)

            rows = [
                list(m_ri) + [Vcol_M],
                list(q_ri) + [mpf(0)],
                list(phip_ri) + [mpf(0)],
                list(m_ro) + [Vcol_M],
                list(q_ro) + [mpf(0)],
                list(phip_ro) + [mpf(0)],
                circuit_row,
            ]
            M = matrix(rows)
            for j in range(7):
                mx = max(abs(M[i, j]) for i in range(7)) or mpf(1)
                for i in range(7):
                    M[i, j] = M[i, j] / mx
            for i in range(7):
                mx = max(abs(M[i, j]) for j in range(7)) or mpf(1)
                for j in range(7):
                    M[i, j] = M[i, j] / mx
            return mp.det(M)

    def oc_ff_sine_bisect(self, lo, hi, n, h1=None, iters=45):
        """Bisect oc_ff_sine_det's real-part sign flip on [lo, hi]
        (rad/s). Returns a plain float. Default iters=45 matches
        coupled_bisect (same 3-branch chi-cubic complexity as the
        underlying 6 sine columns this 7x7 extends)."""
        if h1 is None:
            h1 = self.h1
        with mp.workdps(self.dps):
            lo = mpf(lo); hi = mpf(hi)
            flo = self.oc_ff_sine_det(lo, n, h1=h1)
            for _ in range(iters):
                mid = (lo + hi) / 2
                fm = self.oc_ff_sine_det(mid, n, h1=h1)
                if (fm.real > 0) == (flo.real > 0):
                    lo, flo = mid, fm
                else:
                    hi = mid
            return float((lo + hi) / 2)

    def driven_ff_qv_sine(self, omega, n, volt=1.0, h1=None):
        """Imposed-V driven F-F problem on the unified (linear + sine)
        ansatz -- the same-ansatz counterpart of driven_ff_qv (the
        linear-only model), built on coupled_det's own 6x6 matrix
        instead of the elastic 4x4.

        Returns a dict of plain floats: q, q_over_v, c0, s, sphi.
        n!=0: motional charge vanishes (theta-integral), q = C0*volt,
        matching driven_ff_qv exactly. h1=0: q=0.
        n=0, h1>0: solve coupled_det's 6x6 matrix with a rigid forcing
        e31_bar*(2h+h1)*volt on BOTH M_rr rows' RHS (Q_r, phi'=0 rows
        unforced -- V never appears in them, oc_ff_sine_det's own
        docstring point (1)), then apply the same enriched charge
        formula oc_ff_sine_det's circuit row derives (point (3) there,
        including its sign-selection rationale):
            q = 4*pi*e31_bar*(h+h1)*s + C0*volt - 4*pi*Xi33_bar/h1*sphi
        Poles of q/volt sit at coupled_det's own roots (this method's
        6x6 matrix IS coupled_det's, unaffected by the RHS forcing);
        zeros sit at oc_ff_sine_det's roots -- confirmed to
        q_over_v/C0 ~1e-13 at dps=100 across all three Duan2005
        thickness ratios (Sec 18.182); this is the actual check that
        caught and then confirmed the fix to the sign errors discussed
        in oc_ff_sine_det's own docstring.
        """
        if h1 is None:
            h1 = self.h1
        volt = float(volt)
        c0 = self.clamped_C0(h1=h1)
        if float(h1) == 0.0 or volt == 0.0:
            return dict(q=0.0, q_over_v=0.0, c0=c0, s=0.0, sphi=0.0)
        if int(n) != 0:
            return dict(q=c0 * volt, q_over_v=c0, c0=c0, s=0.0, sphi=0.0)
        with mp.workdps(self.dps):
            omega = mpf(omega)
            h1 = mpf(h1)
            V = mpf(volt)
            d1 = self._d1()
            d2 = self._d2(h1)
            lams = self._branches(omega, h1)
            A1v = self._A1(d2, h1)
            e31_bar = self._e31_bar()
            Xi33_bar = self._Xi33_bar()
            K_pref = self._K_pref(h1)
            h = mpf(self.h)
            ri = mpf(self.r_i)
            ro = mpf(self.r_o)
            n = int(n)

            Mp = e31_bar * (2 * h + h1) * V
            C0 = 2 * Xi33_bar * mp.pi * (ro ** 2 - ri ** 2) / h1

            def rows_at(r):
                r = mpf(r)
                m_row, q_row, phip_row, bterm_row, col_lam = [], [], [], [], []
                for chi, lam in lams:
                    Z1, Z2, dZ1, dZ2 = self._radial_quad(n, r, lam)
                    Ki = (d1 + d2) * lam + K_pref * chi
                    for Z, dZ in ((Z1, dZ1), (Z2, dZ2)):
                        m_row.append((Ki + 2 * A1v * n ** 2 / r ** 2) * Z
                                    - (2 * A1v / r) * dZ)
                        # Kirchhoff effective shear V_r = Q_r + (1/r) dM_rtheta/dtheta
                        # (2026-09-22 fix, LESSONS Sec 18.226): twisting term added;
                        # vanishes identically at n=0, so n=0 results are unchanged.
                        q_row.append(Ki * dZ
                                     - 2 * A1v * n ** 2 * (dZ / r ** 2 - Z / r ** 3))
                        phip_row.append(chi * dZ)
                        bterm_row.append(dZ)
                        col_lam.append((chi, lam))
                return m_row, q_row, phip_row, bterm_row, col_lam

            m_ri, q_ri, phip_ri, b_ri, col_lam = rows_at(ri)
            m_ro, q_ro, phip_ro, b_ro, _ = rows_at(ro)

            A = matrix([m_ri, q_ri, phip_ri, m_ro, q_ro, phip_ro])
            b = matrix([Mp, mpf(0), mpf(0), Mp, mpf(0), mpf(0)])
            # Equilibrate columns before solving -- the 3 chi/lambda
            # branches can differ in magnitude by many orders (Sec
            # 18.182), which otherwise makes mp.lu_solve raise a
            # spurious ZeroDivisionError well short of the working
            # precision's actual resolving power.
            col_scale = []
            for j in range(6):
                mx = max(abs(A[i, j]) for i in range(6)) or mpf(1)
                col_scale.append(mx)
            Aeq = matrix(6, 6)
            for i in range(6):
                for j in range(6):
                    Aeq[i, j] = A[i, j] / col_scale[j]
            try:
                ceq = mp.lu_solve(Aeq, b)
                c = [ceq[j, 0] / col_scale[j] for j in range(6)]
            except Exception:
                return dict(q=float('inf'), q_over_v=float('inf'),
                            c0=float(C0), s=float('nan'), sphi=float('nan'))

            bterm = [ro * b_ro[j] - ri * b_ri[j] for j in range(6)]
            S = sum(c[j] * bterm[j] for j in range(6))
            Sphi = sum(c[j] * (col_lam[j][0] / col_lam[j][1]) * bterm[j]
                       for j in range(6))
            q = (4 * mp.pi * e31_bar * self._charge_arm(h1) * S + C0 * V
                 - self._phi_charge_weight(h1) * Sphi)
            qv = q / V
            return dict(q=float(q), q_over_v=float(qv), c0=float(C0),
                        s=float(S), sphi=float(Sphi))

    # ---------- mixed-edge C-F / F-C (Sec 18.192, lever 10) ----------
    # Paper 3 ring_disk.py: C-F = inner clamped, outer free; F-C = swap.
    # Existing F-F and C-C methods above are unchanged. See
    # PAPER4_CF_DERIVATION.md for the algebra and the OC pre-registration.

    @staticmethod
    def _elastic_bessel_cols(n, k, mu, r):
        """2-branch elastic columns at radius r: (Z, dZ, lam) in
        [I_n, K_n, J_n, Y_n] order. Shared by mixed-edge 4x4s only."""
        r = mpf(r)
        x = k * r
        I = mp.besseli(n, x); K = mp.besselk(n, x)
        J = mp.besselj(n, x); Y = mp.bessely(n, x)
        dI = k * mpf('0.5') * (mp.besseli(n - 1, x) + mp.besseli(n + 1, x))
        dK = k * mpf('-0.5') * (mp.besselk(n - 1, x) + mp.besselk(n + 1, x))
        dJ = k * mpf('0.5') * (mp.besselj(n - 1, x) - mp.besselj(n + 1, x))
        dY = k * mpf('0.5') * (mp.bessely(n - 1, x) - mp.bessely(n + 1, x))
        return [(I, dI, mu), (K, dK, mu), (J, dJ, -mu), (Y, dY, -mu)]

    def _elastic_mixed_det(self, omega, n, h1, inner, outer, oc=False):
        """4x4 elastic mixed-edge det. inner/outer in {'C','F'}.
        oc=True adds the Sec 18.175 rank-1 alpha*S update on FREE
        M_rr rows only (n=0 and h1>0). Both-free recovers elastic_det
        / oc_ff_det; both-clamped recovers elastic_cc_det."""
        with mp.workdps(self.dps):
            omega = mpf(omega)
            d1 = self._d1()
            d2 = self._d2(h1)
            dsum = d1 + d2
            A1 = self._A1(d2, h1)
            A2 = self._A2(h1)
            mu = mp.sqrt(A2 * omega ** 2 / dsum)
            k = mp.sqrt(mu)
            ri = mpf(self.r_i)
            ro = mpf(self.r_o)

            def w_wp_m_q_dz(r, cols):
                r = mpf(r)
                W, Wp, Mrow, Qrow, dZrow = [], [], [], [], []
                for Z, dZ, lam in cols:
                    W.append(Z)
                    Wp.append(dZ)
                    Mrow.append((dsum * lam + 2 * A1 * n ** 2 / r ** 2) * Z
                                - (2 * A1 / r) * dZ)
                    # Kirchhoff effective shear V_r = Q_r + (1/r) dM_rtheta/dtheta
                    # (2026-09-22 fix, LESSONS Sec 18.226): twisting term added;
                    # vanishes identically at n=0, so n=0 results are unchanged.
                    Qrow.append(dsum * lam * dZ
                                - 2 * A1 * n ** 2 * (dZ / r ** 2 - Z / r ** 3))
                    dZrow.append(dZ)
                return W, Wp, Mrow, Qrow, dZrow

            cols_i = self._elastic_bessel_cols(n, k, mu, ri)
            cols_o = self._elastic_bessel_cols(n, k, mu, ro)
            Wi, Wpi, Mi, Qi, dZi = w_wp_m_q_dz(ri, cols_i)
            Wo, Wpo, Mo, Qo, dZo = w_wp_m_q_dz(ro, cols_o)
            if oc and int(n) == 0 and float(h1) != 0.0:
                e31b = self._e31_bar()
                Xi33b = self._Xi33_bar()
                h = mpf(self.h)
                h1v = mpf(h1)
                alpha = (4 * e31b ** 2 * (h + h1v / 2) * self._charge_arm(h1v)
                         * h1v
                         / (Xi33b * (ro ** 2 - ri ** 2)))
                S = [ro * dZo[j] - ri * dZi[j] for j in range(4)]
                if inner == "F":
                    for j in range(4):
                        Mi[j] = Mi[j] + alpha * S[j]
                if outer == "F":
                    for j in range(4):
                        Mo[j] = Mo[j] + alpha * S[j]

            def pair(bc, W, Wp, Mrow, Qrow):
                if bc == "C":
                    return W, Wp
                if bc == "F":
                    return Mrow, Qrow
                raise ValueError(
                    "mixed-edge bc must be 'C' or 'F', got %r" % (bc,))

            r0, r1 = pair(inner, Wi, Wpi, Mi, Qi)
            r2, r3 = pair(outer, Wo, Wpo, Mo, Qo)
            M = matrix([r0, r1, r2, r3])
            for i in range(4):
                s = max(abs(M[i, j]) for j in range(4)) or mpf(1)
                for j in range(4):
                    M[i, j] = M[i, j] / s
            return mp.det(M)

    def _coupled_mixed_det(self, omega, n, h1, inner, outer):
        """6x6 SC mixed-edge det. inner/outer in {'C','F'}. Both-free
        recovers coupled_det; both-clamped recovers cc_coupled_det.
        h1=0 refused (chi-cubic singularity)."""
        if float(h1) == 0.0:
            raise ValueError(
                "coupled mixed-edge det needs h1 > 0 -- h1=0 is a "
                "genuine singular limit of the chi-cubic (Sec 18.149); "
                "use elastic_cf_det/elastic_fc_det with h1=0.0 instead")
        with mp.workdps(self.dps):
            omega = mpf(omega)
            h1 = mpf(h1)
            d1 = self._d1()
            d2 = self._d2(h1)
            lams = self._branches(omega, h1)
            A1v = self._A1(d2, h1)
            e31_bar = self._e31_bar()
            K_pref = self._K_pref(h1)

            def vecs_at(r):
                r = mpf(r)
                w_row, wp_row, m_row, q_row, phip_row = [], [], [], [], []
                for chi, lam in lams:
                    Z1, Z2, dZ1, dZ2 = self._radial_quad(n, r, lam)
                    Ki = (d1 + d2) * lam + K_pref * chi
                    for Z, dZ in ((Z1, dZ1), (Z2, dZ2)):
                        w_row.append(Z)
                        wp_row.append(dZ)
                        m_row.append((Ki + 2 * A1v * n ** 2 / r ** 2) * Z
                                     - (2 * A1v / r) * dZ)
                        # Kirchhoff effective shear V_r = Q_r + (1/r) dM_rtheta/dtheta
                        # (2026-09-22 fix, LESSONS Sec 18.226): twisting term added;
                        # vanishes identically at n=0, so n=0 results are unchanged.
                        q_row.append(Ki * dZ
                                     - 2 * A1v * n ** 2 * (dZ / r ** 2 - Z / r ** 3))
                        phip_row.append(chi * dZ)
                return w_row, wp_row, m_row, q_row, phip_row

            def triple(bc, w, wp, m_row, q_row, phip):
                if bc == "C":
                    return w, wp, phip
                if bc == "F":
                    return m_row, q_row, phip
                raise ValueError(
                    "mixed-edge bc must be 'C' or 'F', got %r" % (bc,))

            wi, wpi, mi, qi, phipi = vecs_at(self.r_i)
            wo, wpo, mo, qo, phipo = vecs_at(self.r_o)
            a0, a1, a2 = triple(inner, wi, wpi, mi, qi, phipi)
            b0, b1, b2 = triple(outer, wo, wpo, mo, qo, phipo)
            M = matrix([a0, a1, a2, b0, b1, b2])
            for j in range(6):
                mx = max(abs(M[i, j]) for i in range(6)) or mpf(1)
                for i in range(6):
                    M[i, j] = M[i, j] / mx
            for i in range(6):
                mx = max(abs(M[i, j]) for j in range(6)) or mpf(1)
                for j in range(6):
                    M[i, j] = M[i, j] / mx
            return mp.det(M)

    def _signflip_bisect(self, det_fn, lo, hi, iters):
        """Real-part sign-flip bisection used by the mixed-edge
        public bisects. Same loop as elastic_bisect."""
        with mp.workdps(self.dps):
            lo = mpf(lo); hi = mpf(hi)
            flo = det_fn(lo)
            for _ in range(iters):
                mid = (lo + hi) / 2
                fm = det_fn(mid)
                if (fm.real > 0) == (flo.real > 0):
                    lo, flo = mid, fm
                else:
                    hi = mid
            return float((lo + hi) / 2)

    def elastic_cf_det(self, omega, n, h1=None):
        """4x4 elastic C-F determinant: inner w=w'=0, outer M_rr=Q_r=0."""
        if h1 is None:
            h1 = self.h1
        return self._elastic_mixed_det(omega, n, h1, "C", "F", oc=False)

    def elastic_cf_bisect(self, lo, hi, n, h1=None, iters=50):
        if h1 is None:
            h1 = self.h1
        return self._signflip_bisect(
            lambda om, n=n, h1=h1: self.elastic_cf_det(om, n, h1=h1),
            lo, hi, iters)

    def elastic_fc_det(self, omega, n, h1=None):
        """4x4 elastic F-C determinant: inner M_rr=Q_r=0, outer w=w'=0."""
        if h1 is None:
            h1 = self.h1
        return self._elastic_mixed_det(omega, n, h1, "F", "C", oc=False)

    def elastic_fc_bisect(self, lo, hi, n, h1=None, iters=50):
        if h1 is None:
            h1 = self.h1
        return self._signflip_bisect(
            lambda om, n=n, h1=h1: self.elastic_fc_det(om, n, h1=h1),
            lo, hi, iters)

    def cf_coupled_det(self, omega, n, h1=None):
        """6x6 C-F SC coupled determinant (inner w=w'=phi'=0, outer
        M_rr=Q_r=phi'=0). Same chi-cubic as coupled_det."""
        if h1 is None:
            h1 = self.h1
        return self._coupled_mixed_det(omega, n, h1, "C", "F")

    def cf_coupled_bisect(self, lo, hi, n, h1=None, iters=45):
        if h1 is None:
            h1 = self.h1
        return self._signflip_bisect(
            lambda om, n=n, h1=h1: self.cf_coupled_det(om, n, h1=h1),
            lo, hi, iters)

    def fc_coupled_det(self, omega, n, h1=None):
        """6x6 F-C SC coupled determinant (inner M_rr=Q_r=phi'=0,
        outer w=w'=phi'=0)."""
        if h1 is None:
            h1 = self.h1
        return self._coupled_mixed_det(omega, n, h1, "F", "C")

    def fc_coupled_bisect(self, lo, hi, n, h1=None, iters=45):
        if h1 is None:
            h1 = self.h1
        return self._signflip_bisect(
            lambda om, n=n, h1=h1: self.fc_coupled_det(om, n, h1=h1),
            lo, hi, iters)

    def oc_cf_det(self, omega, n, h1=None):
        """4x4 C-F open-circuit determinant. Rank-1 alpha*S on the
        OUTER (free) M_rr row only. n!=0 and h1=0 reduce identically
        to elastic_cf_det. n=0 is pre-registered to STIFFEN vs the
        C-F elastic bilayer (PAPER4_CF_DERIVATION.md Sec 4)."""
        if h1 is None:
            h1 = self.h1
        if float(h1) == 0.0 or int(n) != 0:
            return self.elastic_cf_det(omega, n, h1=h1)
        return self._elastic_mixed_det(omega, n, h1, "C", "F", oc=True)

    def oc_cf_bisect(self, lo, hi, n, h1=None, iters=50):
        if h1 is None:
            h1 = self.h1
        return self._signflip_bisect(
            lambda om, n=n, h1=h1: self.oc_cf_det(om, n, h1=h1),
            lo, hi, iters)

    def oc_fc_det(self, omega, n, h1=None):
        """4x4 F-C open-circuit determinant. Rank-1 alpha*S on the
        INNER (free) M_rr row only. n!=0 and h1=0 reduce identically
        to elastic_fc_det."""
        if h1 is None:
            h1 = self.h1
        if float(h1) == 0.0 or int(n) != 0:
            return self.elastic_fc_det(omega, n, h1=h1)
        return self._elastic_mixed_det(omega, n, h1, "F", "C", oc=True)

    def oc_fc_bisect(self, lo, hi, n, h1=None, iters=50):
        if h1 is None:
            h1 = self.h1
        return self._signflip_bisect(
            lambda om, n=n, h1=h1: self.oc_fc_det(om, n, h1=h1),
            lo, hi, iters)
