# -*- coding: utf-8 -*-
"""
plate_solver.piezo_monolithic -- Paper 5 Option A: homogeneous
thickness-polarized piezoceramic ring, Kirchhoff, short-circuit 6x6.

Background: PAPER5_DERIVATION.md (Phase 0), PAPER5_MONOLITHIC_PIEZO_ROADMAP.md
(claim A locked on PZT-4), LESSONS_LEARNED.md Sec 18.199-18.200.

SCOPE: one ceramic spanning |z| <= H, electrodes on the two faces only.
Not a flag on PiezoOutOfPlaneSolver (Paper 4 layered sandwich). Not
Parashar d15. No open-circuit frequency solver in this module -- Phase 0
theorem: fully-electroded OC flexure coincides with SC at this ansatz.

Methods:
  elastic_det / elastic_bisect       F-F 4x4, e31 unused
  elastic_cc_det / elastic_cc_bisect C-C 4x4, e31 unused
  coupled_det / coupled_bisect       F-F 6x6 SC (chi-cubic, 3 Helmholtz branches)
  cc_coupled_det / cc_coupled_bisect C-C 6x6 SC

SOLVER_VERSION is not bumped: new unused-by-default entry point, same
pattern as piezo_solver.py / ring_disk.py.

ELECTRICAL PROJECTION (2026-09-23, LESSONS_LEARNED.md Sec 18.231,
SC_PROJECTION_CONSISTENT_DERIVATION_2026-09-23.md): the default
projection='consistent' is the same-weight Galerkin projection of the
electric enthalpy with the quadratic through-thickness potential
phi = phibar*(1 - z^2/H^2). Its long-wave SC stiffening is the exact
thin-plate value (2/3) H^3 e31bar^2/Xi33bar. projection='duan' is the
pre-2026-09-23 model (sine potential closed by Duan 2005 Eq. 12's
charge-conservation integral), kept bit-identical for provenance; it
overstates the long-wave SC increment by exactly 12/pi^2.
projection='galerkin_sine' is the consistent projection of the sine
profile (factor 96/pi^4), a diagnostic only.

Units: dimensional SI (m, Pa, kg/m^3, C/m^2, F/m, rad/s).
"""
from __future__ import annotations

from mpmath import mp, mpf, matrix


class PiezoMonolithicOutOfPlaneSolver:
    """Homogeneous thickness-polarized Kirchhoff ring (Paper 5 Option A).

    Constructor takes RAW Duan/PZT-4 stiffnesses (constant-E). Reduced
    bars are formed internally, same convention as PiezoOutOfPlaneSolver
    (Sec 18.172). H is the HALF-thickness; full thickness is 2H.

    projection: 'consistent' (default, quadratic Galerkin, exact long-wave
    SC stiffening), 'duan' (legacy, 12/pi^2 high), or 'galerkin_sine'
    (diagnostic, 96/pi^4). Elastic methods ignore it.
    """

    _PROJECTIONS = ('consistent', 'duan', 'galerkin_sine')

    def __init__(self, r_i, r_o, H,
                 C11E, C12E, C13E, C33E, rho,
                 e31=None, e33=None, X11=None, X33=None,
                 dps=100, projection='consistent'):
        self.r_i = float(r_i)
        self.r_o = float(r_o)
        self.H = float(H)
        self.C11E = float(C11E)
        self.C12E = float(C12E)
        self.C13E = float(C13E)
        self.C33E = float(C33E)
        self.rho = float(rho)
        self.e31 = None if e31 is None else float(e31)
        self.e33 = None if e33 is None else float(e33)
        self.X11 = None if X11 is None else float(X11)
        self.X33 = None if X33 is None else float(X33)
        self.dps = int(dps)
        if projection not in self._PROJECTIONS:
            raise ValueError(
                "projection must be one of %r, got %r"
                % (self._PROJECTIONS, projection))
        self.projection = projection
        if self.H <= 0.0:
            raise ValueError("H (half-thickness) must be > 0")
        if self.r_i <= 0.0 or self.r_o <= self.r_i:
            raise ValueError("need 0 < r_i < r_o")

    # ---------- reduced constants (Duan/Liu, same as Paper 4) ----------

    def _c11_bar(self):
        C11E = mpf(self.C11E); C13E = mpf(self.C13E); C33E = mpf(self.C33E)
        return C11E - C13E ** 2 / C33E

    def _c12_bar(self):
        C12E = mpf(self.C12E); C13E = mpf(self.C13E); C33E = mpf(self.C33E)
        return C12E - C13E ** 2 / C33E

    def _e31_bar(self):
        if self.e31 is None or self.e33 is None:
            raise ValueError(
                "the piezo-coupled model needs e31 and e33 (the elastic "
                "path does not)")
        e31 = mpf(self.e31); e33 = mpf(self.e33)
        C13E = mpf(self.C13E); C33E = mpf(self.C33E)
        return e31 - (C13E / C33E) * e33

    def _Xi33_bar(self):
        if self.X33 is None or self.e33 is None:
            raise ValueError(
                "the piezo-coupled model needs X33 and e33 (the elastic "
                "path does not)")
        X33 = mpf(self.X33); e33 = mpf(self.e33); C33E = mpf(self.C33E)
        return X33 + e33 ** 2 / C33E

    def _Xi11_bar(self):
        if self.X11 is None:
            raise ValueError(
                "the piezo-coupled model needs X11 (the elastic path "
                "does not)")
        return mpf(self.X11)

    def isotropic_nu(self):
        """nu_p = c12_bar / c11_bar for the elastic-limit ring_disk gate."""
        with mp.workdps(self.dps):
            return float(self._c12_bar() / self._c11_bar())

    def isotropic_E(self):
        """E = c11_bar * (1 - nu_p^2) so D_iso matches d = (2/3) c11_bar H^3."""
        with mp.workdps(self.dps):
            c11 = self._c11_bar(); c12 = self._c12_bar()
            nu = c12 / c11
            return float(c11 * (1 - nu ** 2))

    # ---------- elastic building blocks (PAPER5_DERIVATION.md Sec 1) ----------

    def _d(self):
        """Bending stiffness d = (2/3) c11_bar H^3 (full thickness 2H)."""
        H = mpf(self.H)
        return mpf(2) / 3 * self._c11_bar() * H ** 3

    def _A1(self):
        """A1 = 0.5 * (1 - nu_p) * d."""
        d = self._d()
        nu = self._c12_bar() / self._c11_bar()
        return mpf("0.5") * (1 - nu) * d

    def _A2(self):
        """A2 = 2 * rho * H."""
        return 2 * mpf(self.rho) * mpf(self.H)

    def _elastic_k(self, omega):
        d = self._d()
        A2 = self._A2()
        mu = mp.sqrt(A2 * mpf(omega) ** 2 / d)
        return mu, mp.sqrt(mu)

    @staticmethod
    def _bessel_quad(n, k, r):
        """I, K, J, Y and radial derivatives at argument k*r."""
        r = mpf(r)
        x = k * r
        I = mp.besseli(n, x); K = mp.besselk(n, x)
        J = mp.besselj(n, x); Y = mp.bessely(n, x)
        dI = k * mpf("0.5") * (mp.besseli(n - 1, x) + mp.besseli(n + 1, x))
        dK = k * mpf("-0.5") * (mp.besselk(n - 1, x) + mp.besselk(n + 1, x))
        dJ = k * mpf("0.5") * (mp.besselj(n - 1, x) - mp.besselj(n + 1, x))
        dY = k * mpf("0.5") * (mp.bessely(n - 1, x) - mp.bessely(n + 1, x))
        return I, K, J, Y, dI, dK, dJ, dY

    def elastic_det(self, omega, n):
        """4x4 elastic F-F determinant (M_rr = Q_r = 0 at r_i, r_o)."""
        with mp.workdps(self.dps):
            omega = mpf(omega)
            d = self._d()
            A1 = self._A1()
            mu, k = self._elastic_k(omega)

            def rows(r):
                I, K, J, Y, dI, dK, dJ, dY = self._bessel_quad(n, k, r)
                r = mpf(r)
                cols = [(I, dI, mu), (K, dK, mu), (J, dJ, -mu), (Y, dY, -mu)]
                Mrow, Qrow = [], []
                for Z, dZ, lam in cols:
                    Mrow.append((d * lam + 2 * A1 * n ** 2 / r ** 2) * Z
                                - (2 * A1 / r) * dZ)
                    # Kirchhoff effective shear V_r = Q_r + (1/r) dM_rtheta/dtheta
                    # (2026-09-22 fix, LESSONS Sec 18.226): twisting term added;
                    # vanishes identically at n=0, so n=0 results are unchanged.
                    Qrow.append(d * lam * dZ
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

    def elastic_bisect(self, lo, hi, n, iters=50):
        with mp.workdps(self.dps):
            lo = mpf(lo); hi = mpf(hi)
            flo = self.elastic_det(lo, n)
            for _ in range(iters):
                mid = (lo + hi) / 2
                fm = self.elastic_det(mid, n)
                if (fm.real > 0) == (flo.real > 0):
                    lo, flo = mid, fm
                else:
                    hi = mid
            return float((lo + hi) / 2)

    def elastic_cc_det(self, omega, n):
        """4x4 elastic C-C determinant (w = w' = 0 at r_i, r_o)."""
        with mp.workdps(self.dps):
            omega = mpf(omega)
            mu, k = self._elastic_k(omega)

            def row(r):
                I, K, J, Y, dI, dK, dJ, dY = self._bessel_quad(n, k, r)
                return [I, K, J, Y], [dI, dK, dJ, dY]

            w_ri, wp_ri = row(self.r_i)
            w_ro, wp_ro = row(self.r_o)
            M = matrix([w_ri, wp_ri, w_ro, wp_ro])
            for i in range(4):
                s = max(abs(M[i, j]) for j in range(4)) or mpf(1)
                for j in range(4):
                    M[i, j] = M[i, j] / s
            return mp.det(M)

    def elastic_cc_bisect(self, lo, hi, n, iters=50):
        with mp.workdps(self.dps):
            lo = mpf(lo); hi = mpf(hi)
            flo = self.elastic_cc_det(lo, n)
            for _ in range(iters):
                mid = (lo + hi) / 2
                fm = self.elastic_cc_det(mid, n)
                if (fm.real > 0) == (flo.real > 0):
                    lo, flo = mid, fm
                else:
                    hi = mid
            return float((lo + hi) / 2)

    # ---------- electrical projection (Sec 18.231) ----------
    # Plate equations for every projection (derivation file Sec 2):
    #   (A) d Lap^2 w + K_pref Lap phibar - A2 omega^2 w = 0
    #   (B) a Lap phibar - phibar + c Lap w = 0
    # so lam = chi/(a chi + c). 'duan' keeps its original code path in
    # _cubic_coeffs/_branches so it stays bit-identical.

    def _proj_integrals(self):
        """(I0, I1, I2, I3) = int f, int z f', int f'^2, int f^2 over
        |z| <= H for the Galerkin projections (derivation Sec 3.1)."""
        H = mpf(self.H)
        pi = mp.pi
        if self.projection == 'consistent':
            return 4 * H / 3, -4 * H / 3, mpf(8) / (3 * H), 16 * H / 15
        if self.projection == 'galerkin_sine':
            return 4 * H / pi, -4 * H / pi, pi ** 2 / (4 * H), H
        raise ValueError("'duan' is not a Galerkin projection")

    def _K_pref(self):
        """Piezo moment coefficient: M_rr contains -K_pref*phibar."""
        H = mpf(self.H)
        if self.projection == 'duan':
            return (4 * H / mp.pi) * self._e31_bar()
        _, I1, _, _ = self._proj_integrals()
        return -I1 * self._e31_bar()

    def _elec_ac(self):
        """(a, c) of the normalised electrical equation (B)."""
        H = mpf(self.H)
        Xi33 = self._Xi33_bar(); Xi11 = self._Xi11_bar()
        e31_bar = self._e31_bar()
        if self.projection == 'duan':
            return (4 * H ** 2 * Xi11 / (mp.pi ** 2 * Xi33),
                    2 * H ** 2 * e31_bar / (mp.pi * Xi33))
        _, I1, I2, I3 = self._proj_integrals()
        return Xi11 * I3 / (Xi33 * I2), -I1 * e31_bar / (Xi33 * I2)

    def _phi_edge_weight(self):
        """Thickness weight of the in-plane electric flux pairing on an
        edge (Paper 6 theta-edge D_theta term): Galerkin I3; 'duan'
        kept the Gauss weight int f = 4H/pi."""
        H = mpf(self.H)
        if self.projection == 'duan':
            return 4 * H / mp.pi
        return self._proj_integrals()[3]

    def _face_flux_weight(self):
        """pi*I0: Q_band = pi*I0*Xi11*[r phibar'] (derivation Sec 3.1,
        Gauss flux of the leading-order field). 'duan' = 4H exactly."""
        H = mpf(self.H)
        if self.projection == 'duan':
            return 4 * H
        return mp.pi * self._proj_integrals()[0]

    def long_wave_sc_ratio(self, omega):
        """(D_eff - d)/[(2/3) H^3 e31bar^2/Xi33bar] from the flexural
        branch of _branches (smallest |lam|), D_eff = A2 omega^2/lam^2.
        -> 1 (consistent), 96/pi^4 (galerkin_sine), 12/pi^2 (duan) as
        omega -> 0. Gate P2 of Sec 18.231."""
        self._require_coupled_consts()
        with mp.workdps(self.dps):
            omega = mpf(omega)
            br = self._branches(omega)
            chi, lam = min(br, key=lambda cl: abs(cl[1]))
            Deff = self._A2() * omega ** 2 / lam ** 2
            exact = (mpf(2) / 3 * mpf(self.H) ** 3 * self._e31_bar() ** 2
                     / self._Xi33_bar())
            return (Deff.real - self._d()) / exact

    # ---------- SC coupled 6x6 (PAPER5_DERIVATION.md Sec 2) ----------

    def _require_coupled_consts(self):
        if self.e31 is None or self.e33 is None or self.X11 is None \
                or self.X33 is None:
            raise ValueError(
                "coupled methods need e31, e33, X11, X33; use elastic_det "
                "/ elastic_cc_det for the e31=0 limit")
        if float(self.e31) == 0.0:
            raise ValueError(
                "e31=0 is a genuine singular limit of the chi-cubic "
                "(PAPER5_DERIVATION.md Sec 2); use elastic_det / "
                "elastic_cc_det instead")

    def _cubic_coeffs(self, omega):
        """a3 chi^3 + a2 chi^2 + a1 chi + a0 = 0 (PAPER5_DERIVATION.md;
        Galerkin form: SC_PROJECTION_CONSISTENT_DERIVATION Sec 2)."""
        if self.projection != 'duan':
            omega = mpf(omega)
            w2 = omega ** 2
            d = self._d()
            A2 = self._A2()
            K = self._K_pref()
            a, c = self._elec_ac()
            return (K * a,
                    d + K * c - A2 * w2 * a ** 2,
                    -2 * A2 * w2 * a * c,
                    -A2 * w2 * c ** 2)
        pi = mp.pi
        omega = mpf(omega)
        w2 = omega ** 2
        H = mpf(self.H)
        d = self._d()
        A2 = self._A2()
        e31_bar = self._e31_bar()
        Xi11 = self._Xi11_bar()
        Xi33 = self._Xi33_bar()
        beta = (4 * H / pi) * e31_bar
        alpha = 4 * H ** 2 * Xi11
        gamma = 2 * pi * H ** 2 * e31_bar
        a3 = beta * pi ** 2 * Xi33 * alpha
        a2 = (d * pi ** 4 * Xi33 ** 2
              + beta * pi ** 2 * Xi33 * gamma
              - A2 * w2 * alpha ** 2)
        a1 = -2 * A2 * w2 * alpha * gamma
        a0 = -A2 * w2 * gamma ** 2
        return a3, a2, a1, a0

    def _branches(self, omega):
        pi = mp.pi
        H = mpf(self.H)
        a3, a2, a1, a0 = self._cubic_coeffs(omega)
        roots = mp.polyroots([a3, a2, a1, a0], maxsteps=300, extraprec=1500)
        if self.projection != 'duan':
            a, c = self._elec_ac()
            return [(chi, chi / (a * chi + c)) for chi in roots]
        Xi33 = self._Xi33_bar(); Xi11 = self._Xi11_bar()
        e31_bar = self._e31_bar()
        out = []
        for chi in roots:
            lam = (pi ** 2 * Xi33 * chi
                   / (2 * H ** 2 * (2 * Xi11 * chi + pi * e31_bar)))
            out.append((chi, lam))
        return out

    @staticmethod
    def _radial_quad(n, r, lam):
        """Z1, Z2, dZ1, dZ2 for Delta w = lam w. Copied from
        PiezoOutOfPlaneSolver._radial_quad (same Bessel identities)."""
        r = mpf(r)
        if lam.real >= 0:
            delta = mp.sqrt(lam); x = delta * r
            Z1 = mp.besseli(n, x); Z2 = mp.besselk(n, x)
            dZ1 = delta * mpf("0.5") * (mp.besseli(n - 1, x) + mp.besseli(n + 1, x))
            dZ2 = delta * mpf("-0.5") * (mp.besselk(n - 1, x) + mp.besselk(n + 1, x))
        else:
            delta = mp.sqrt(-lam); x = delta * r
            Z1 = mp.besselj(n, x); Z2 = mp.bessely(n, x)
            dZ1 = delta * mpf("0.5") * (mp.besselj(n - 1, x) - mp.besselj(n + 1, x))
            dZ2 = delta * mpf("0.5") * (mp.bessely(n - 1, x) - mp.bessely(n + 1, x))
        return Z1, Z2, dZ1, dZ2

    def _equilibrate6(self, M):
        for j in range(6):
            mx = max(abs(M[i, j]) for i in range(6)) or mpf(1)
            for i in range(6):
                M[i, j] = M[i, j] / mx
        for i in range(6):
            mx = max(abs(M[i, j]) for j in range(6)) or mpf(1)
            for j in range(6):
                M[i, j] = M[i, j] / mx
        return M

    def coupled_det(self, omega, n):
        """6x6 F-F SC determinant: M_rr, Q_r, phi' = 0 at each edge."""
        self._require_coupled_consts()
        with mp.workdps(self.dps):
            omega = mpf(omega)
            d = self._d()
            A1v = self._A1()
            H = mpf(self.H)
            e31_bar = self._e31_bar()
            K_pref = self._K_pref()
            lams = self._branches(omega)

            def rows_at(r):
                r = mpf(r)
                m_row, q_row, phip_row = [], [], []
                for chi, lam in lams:
                    Z1, Z2, dZ1, dZ2 = self._radial_quad(n, r, lam)
                    Ki = d * lam + K_pref * chi
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
            return mp.det(self._equilibrate6(M))

    def coupled_bisect(self, lo, hi, n, iters=45):
        with mp.workdps(self.dps):
            lo = mpf(lo); hi = mpf(hi)
            flo = self.coupled_det(lo, n)
            for _ in range(iters):
                mid = (lo + hi) / 2
                fm = self.coupled_det(mid, n)
                if (fm.real > 0) == (flo.real > 0):
                    lo, flo = mid, fm
                else:
                    hi = mid
            return float((lo + hi) / 2)

    def cc_coupled_det(self, omega, n):
        """6x6 C-C SC determinant: w, w', phi' = 0 at each edge."""
        self._require_coupled_consts()
        with mp.workdps(self.dps):
            omega = mpf(omega)
            lams = self._branches(omega)

            def rows_at(r):
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
            return mp.det(self._equilibrate6(M))

    def cc_coupled_bisect(self, lo, hi, n, iters=40):
        with mp.workdps(self.dps):
            lo = mpf(lo); hi = mpf(hi)
            flo = self.cc_coupled_det(lo, n)
            for _ in range(iters):
                mid = (lo + hi) / 2
                fm = self.cc_coupled_det(mid, n)
                if (fm.real > 0) == (flo.real > 0):
                    lo, flo = mid, fm
                else:
                    hi = mid
            return float((lo + hi) / 2)

    # ---------- mixed-edge C-F / F-C (Paper 4 Sec 18.192 pattern) ----------
    # Same physical family as PiezoOutOfPlaneSolver._elastic_mixed_det /
    # _coupled_mixed_det (piezo_solver.py). C-F = inner clamped, outer
    # free; F-C is the swap. Both-free recovers elastic_det/coupled_det;
    # both-clamped recovers elastic_cc_det/cc_coupled_det -- this
    # reduction is the regression gate for the new row-selection logic
    # (see probe_piezo_p5_cf_boundary_2026-09-20.py's G_reduce check).
    # No open-circuit variant: Phase 0's OC=SC theorem covers this too,
    # and there is still no OC frequency solver in this module.

    def _signflip_bisect(self, det_fn, lo, hi, iters):
        """Real-part sign-flip bisection shared by the mixed-edge
        public bisects (same loop as elastic_bisect/coupled_bisect)."""
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

    def _elastic_mixed_det(self, omega, n, inner, outer):
        """4x4 elastic mixed-edge det. inner/outer in {'C','F'}.
        Both-free recovers elastic_det; both-clamped recovers
        elastic_cc_det (bit-identical row order and equilibration)."""
        with mp.workdps(self.dps):
            omega = mpf(omega)
            d = self._d()
            A1 = self._A1()
            mu, k = self._elastic_k(omega)

            def rows(r):
                I, K, J, Y, dI, dK, dJ, dY = self._bessel_quad(n, k, r)
                r = mpf(r)
                cols = [(I, dI, mu), (K, dK, mu), (J, dJ, -mu), (Y, dY, -mu)]
                W, Wp, Mrow, Qrow = [], [], [], []
                for Z, dZ, lam in cols:
                    W.append(Z)
                    Wp.append(dZ)
                    Mrow.append((d * lam + 2 * A1 * n ** 2 / r ** 2) * Z
                                - (2 * A1 / r) * dZ)
                    # Kirchhoff effective shear V_r = Q_r + (1/r) dM_rtheta/dtheta
                    # (2026-09-22 fix, LESSONS Sec 18.226): twisting term added;
                    # vanishes identically at n=0, so n=0 results are unchanged.
                    Qrow.append(d * lam * dZ
                                - 2 * A1 * n ** 2 * (dZ / r ** 2 - Z / r ** 3))
                return W, Wp, Mrow, Qrow

            Wi, Wpi, Mi, Qi = rows(self.r_i)
            Wo, Wpo, Mo, Qo = rows(self.r_o)

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

    def elastic_cf_det(self, omega, n):
        """4x4 elastic C-F determinant: inner w=w'=0, outer M_rr=Q_r=0."""
        return self._elastic_mixed_det(omega, n, "C", "F")

    def elastic_cf_bisect(self, lo, hi, n, iters=50):
        return self._signflip_bisect(
            lambda om, n=n: self.elastic_cf_det(om, n), lo, hi, iters)

    def elastic_fc_det(self, omega, n):
        """4x4 elastic F-C determinant: inner M_rr=Q_r=0, outer w=w'=0."""
        return self._elastic_mixed_det(omega, n, "F", "C")

    def elastic_fc_bisect(self, lo, hi, n, iters=50):
        return self._signflip_bisect(
            lambda om, n=n: self.elastic_fc_det(om, n), lo, hi, iters)

    def _coupled_mixed_det(self, omega, n, inner, outer):
        """6x6 SC mixed-edge det. inner/outer in {'C','F'}. Both-free
        recovers coupled_det; both-clamped recovers cc_coupled_det
        (bit-identical row order and equilibration)."""
        self._require_coupled_consts()
        with mp.workdps(self.dps):
            omega = mpf(omega)
            d = self._d()
            A1v = self._A1()
            H = mpf(self.H)
            e31_bar = self._e31_bar()
            K_pref = self._K_pref()
            lams = self._branches(omega)

            def vecs_at(r):
                r = mpf(r)
                w_row, wp_row, m_row, q_row, phip_row = [], [], [], [], []
                for chi, lam in lams:
                    Z1, Z2, dZ1, dZ2 = self._radial_quad(n, r, lam)
                    Ki = d * lam + K_pref * chi
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
            return mp.det(self._equilibrate6(M))

    def cf_coupled_det(self, omega, n):
        """6x6 C-F SC coupled determinant (inner w=w'=phi'=0, outer
        M_rr=Q_r=phi'=0). Same chi-cubic as coupled_det/cc_coupled_det."""
        return self._coupled_mixed_det(omega, n, "C", "F")

    def cf_coupled_bisect(self, lo, hi, n, iters=45):
        return self._signflip_bisect(
            lambda om, n=n: self.cf_coupled_det(om, n), lo, hi, iters)

    def fc_coupled_det(self, omega, n):
        """6x6 F-C SC coupled determinant (inner M_rr=Q_r=phi'=0,
        outer w=w'=phi'=0)."""
        return self._coupled_mixed_det(omega, n, "F", "C")

    def fc_coupled_bisect(self, lo, hi, n, iters=45):
        return self._signflip_bisect(
            lambda om, n=n: self.fc_coupled_det(om, n), lo, hi, iters)

    # ---------- force-driven sensing admittance (Sec 18.211/18.212) ----------
    # PAPER5_YOMEGA_SENSE_DERIVATION.md. F-F mechanical + SC electrical
    # only, n=0 only, this pass. A harmonic axisymmetric ring load F at
    # r_F splits the domain into two regions, each carrying the SAME
    # three chi-cubic branches as coupled_det/_branches; only the
    # Kirchhoff shear (the "q" bracket already used by coupled_det, NOT
    # w, w', M_rr, phibar, or phibar') jumps at r_F. This is NOT
    # Y(omega)=Q_full/V (Sec 18.210: voltage-driven flexural admittance
    # is algebraically Y=j*omega*C0, k_eff^2=0 identically, and the
    # reciprocal full-face force-driven charge Q_full is ALSO
    # identically zero on this stacking -- Sec 18.211 Eq (T)). The
    # well-posed, resonance-showing quantity is the charge on a
    # SEGMENT of the still-fully-electroded, still-grounded face:
    # Q_segment/Y_sense below. Do not implement Q_full/F as an
    # admittance -- it is provably, algebraically zero, not a bug.

    def _driven_vecs_at(self, n, r, lams, d, A1v, K_pref):
        """w, w', phibar, m(=-M_rr bracket), q(=-Q_r bracket),
        phip(=phibar') row vectors at radius r, on the SAME chi-cubic
        branches as coupled_det/cc_coupled_det. phi_row (=chi*Z) is
        new here -- coupled_det/cc_coupled_det never need raw phibar
        continuity, only phibar'=0 at the rims."""
        r = mpf(r)
        w_row, wp_row, phi_row, m_row, q_row, phip_row = [], [], [], [], [], []
        for chi, lam in lams:
            Z1, Z2, dZ1, dZ2 = self._radial_quad(n, r, lam)
            Ki = d * lam + K_pref * chi
            for Z, dZ in ((Z1, dZ1), (Z2, dZ2)):
                w_row.append(Z)
                wp_row.append(dZ)
                phi_row.append(chi * Z)
                m_row.append((Ki + 2 * A1v * n ** 2 / r ** 2) * Z
                              - (2 * A1v / r) * dZ)
                # Kirchhoff effective shear V_r = Q_r + (1/r) dM_rtheta/dtheta
                # (2026-09-22 fix, LESSONS Sec 18.226): twisting term added;
                # vanishes identically at n=0, so n=0 results are unchanged.
                q_row.append(Ki * dZ
                             - 2 * A1v * n ** 2 * (dZ / r ** 2 - Z / r ** 3))
                phip_row.append(chi * dZ)
        return w_row, wp_row, phi_row, m_row, q_row, phip_row

    def _equilibrate_solve(self, A, b):
        """Row+column-equilibrated mp.lu_solve. The driven 12x12's
        branch with the largest |lam| produces Bessel arguments large
        enough (kr up to ~100 at this geometry) that raw matrix entries
        span >50 orders of magnitude -- mp.lu_solve reports "numerically
        singular" on the unscaled system even off resonance. Column-
        then-row equilibration (same idea as _equilibrate6, but that
        one only needs a determinant's zero-crossing, which is scale-
        invariant; a driven SOLVE needs the unscale-back step here)."""
        n_dim = A.rows
        colscale = [mpf(1)] * n_dim
        for j in range(n_dim):
            mx = max(abs(A[i, j]) for i in range(n_dim)) or mpf(1)
            colscale[j] = mx
            for i in range(n_dim):
                A[i, j] = A[i, j] / mx
        for i in range(n_dim):
            mx = max(abs(A[i, j]) for j in range(n_dim)) or mpf(1)
            for j in range(n_dim):
                A[i, j] = A[i, j] / mx
            b[i] = b[i] / mx
        y = mp.lu_solve(A, b)
        return [y[j] / colscale[j] for j in range(n_dim)]

    def driven_ff_force_sc(self, omega, r_F, F=1.0, n=0):
        """Force-driven F-F SC response: axisymmetric ring load F at
        r_i < r_F < r_o, two-region 12x12 solve (PAPER5_YOMEGA_SENSE_
        DERIVATION.md Sec 3-4). Row/column convention: columns 0-5 are
        region I (r_i..r_F) branch constants, columns 6-11 are region
        II (r_F..r_o); every row is built as (region I trace - region
        II trace), continuity rows have RHS 0, and the Q_r-bracket jump
        row has RHS -F/(2*pi*r_F) (Sec 4's code-convention table,
        verified against the global force identity (G) in this
        session's validation -- see LESSONS_LEARNED.md Sec 18.212).

        Returns a dict of mpf/plain values sufficient to reconstruct
        w(r) and phibar'(r) anywhere via _driven_eval: c_I, c_II (each
        a length-6 list of branch constants), lams, r_F, n, H.
        n=0 only this pass; C-C/C-F/F-C sensing variants are follow-on.
        """
        self._require_coupled_consts()
        if int(n) != 0:
            raise NotImplementedError(
                "driven_ff_force_sc: n=0 only this pass "
                "(PAPER5_YOMEGA_SENSE_DERIVATION.md Sec 10)")
        if not (float(self.r_i) < float(r_F) < float(self.r_o)):
            raise ValueError("r_F must satisfy r_i < r_F < r_o")
        with mp.workdps(self.dps):
            omega = mpf(omega)
            r_F = mpf(r_F)
            F = mpf(F)
            d = self._d()
            A1v = self._A1()
            H = mpf(self.H)
            e31_bar = self._e31_bar()
            K_pref = self._K_pref()
            lams = self._branches(omega)

            def vecs(r):
                return self._driven_vecs_at(n, r, lams, d, A1v, K_pref)

            w_i, wp_i, phi_i, m_i, q_i, phip_i = vecs(self.r_i)
            w_o, wp_o, phi_o, m_o, q_o, phip_o = vecs(self.r_o)
            w_F, wp_F, phi_F, m_F, q_F, phip_F = vecs(r_F)

            zero6 = [mpf(0)] * 6

            def interface_row(vec):
                return list(vec) + [-x for x in vec]

            rows = [
                m_i + zero6,
                q_i + zero6,
                phip_i + zero6,
                zero6 + m_o,
                zero6 + q_o,
                zero6 + phip_o,
                interface_row(w_F),
                interface_row(wp_F),
                interface_row(m_F),
                interface_row(q_F),
                interface_row(phi_F),
                interface_row(phip_F),
            ]
            rhs = [mpf(0), mpf(0), mpf(0), mpf(0), mpf(0), mpf(0),
                   mpf(0), mpf(0), mpf(0),
                   -F / (2 * mp.pi * r_F),
                   mpf(0), mpf(0)]

            A_mat = matrix(rows)
            b_vec = matrix(rhs)
            c = self._equilibrate_solve(A_mat, b_vec)

            return dict(
                omega=omega, r_F=r_F, F=F, lams=lams,
                c_I=[c[j] for j in range(6)],
                c_II=[c[j] for j in range(6, 12)],
                d=d, A1v=A1v, H=H, K_pref=K_pref, n=int(n),
            )

    def _driven_eval(self, result, r):
        """w(r), phibar'(r) from a driven_ff_force_sc result, picking
        region I or II by r vs r_F (they agree at r_F by construction:
        interface rows 7/12 in driven_ff_force_sc enforce continuity)."""
        with mp.workdps(self.dps):
            r = mpf(r)
            c = result["c_I"] if r <= result["r_F"] else result["c_II"]
            n = result["n"]
            w = mpf(0)
            phip = mpf(0)
            idx = 0
            for chi, lam in result["lams"]:
                Z1, Z2, dZ1, dZ2 = self._radial_quad(n, r, lam)
                for Z, dZ in ((Z1, dZ1), (Z2, dZ2)):
                    w = w + c[idx] * Z
                    phip = phip + c[idx] * chi * dZ
                    idx += 1
            return w, phip

    def Q_segment(self, result, r_a, r_b):
        """Q_Omega = W*Xi11_bar*(r_b*phibar'(r_b) - r_a*phibar'(r_a)),
        W = pi*int f dz (_face_flux_weight: 4H for 'duan', 4 pi H/3 for
        the default quadratic projection; Sec 18.231),
        the induced SC charge on the annular band [r_a, r_b] of the
        still fully-electroded, still-grounded face (PAPER5_YOMEGA_
        SENSE_DERIVATION.md Sec 7, Eq Q_Omega). Valid for r_a, r_b in
        [r_i, r_o]; phibar/phibar' are continuous across r_F by
        construction, so this collapses to boundary values regardless
        of which side of r_F they fall on."""
        with mp.workdps(self.dps):
            H = result["H"]
            Xi11 = self._Xi11_bar()
            _, phip_a = self._driven_eval(result, r_a)
            _, phip_b = self._driven_eval(result, r_b)
            return (self._face_flux_weight() * Xi11
                    * (mpf(r_b) * phip_b - mpf(r_a) * phip_a))

    def Y_sense(self, omega, r_F, r_star, F=1.0, n=0):
        """Y_sense[r_star](omega) = Q_in(r_star; omega) / F, the
        force-driven sensing admittance of the inner partition
        [r_i, r_star] (PAPER5_YOMEGA_SENSE_DERIVATION.md Sec 7, Eq Y).
        NOT Q_full/F (Sec 18.210/18.211: that is identically zero on
        this stacking). Poles sit at the F-F SC coupled_bisect roots,
        unless w(r_F)=0 or phibar'(r_star)=0 at that mode. n=0 only."""
        res = self.driven_ff_force_sc(omega, r_F, F=F, n=n)
        with mp.workdps(self.dps):
            q_in = self.Q_segment(res, self.r_i, r_star)
            return q_in / mpf(F)
