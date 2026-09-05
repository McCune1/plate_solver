# -*- coding: utf-8 -*-
"""
plate_solver.geometry -- plate/material geometry classes, the isotropic and
polar-orthotropic material models, the lightweight worker-side scalar
reconstructions (_ScalarGeom/_ScalarMat), and geometry-table helpers.

Extracted verbatim (line-range provenance in LESSONS_LEARNED Sec. 22) from
Research50.py during the 2026-07-01 modularization; no numeric behaviour
changed, SOLVER_VERSION not bumped.
"""
from __future__ import annotations
import os
import numpy as np
from mpmath import mp, mpf, mpc, matrix

class PlateGeometry:
    """Annular-sector plate geometry in the Seok–Tiersten nondimensionalization.

    Stores physical half-thickness h, half-angle Theta, and the derived midradius
    quantities b=(R_o-R_i)/2, r_0=(R_o+R_i)/2, and r0_bar=pi*r_0/(2b) used by the
    exact-edge Frobenius engines and variational assemblers.
    """
    def __init__(self, R_i, R_o, h, Theta):
        """R_i/R_o: inner/outer radius; h: half-thickness; Theta: half sector angle."""
        self.R_i = mpf(R_i); self.R_o = mpf(R_o)
        self.h = mpf(h);     self.Theta = mpf(Theta)
        self.b = (self.R_o - self.R_i)/2
        self.r_0 = (self.R_o + self.R_i)/2
        self.r0_bar = mp.pi*self.r_0/(2*self.b)
    def __repr__(self):
        return (f"PlateGeometry(R_i={float(self.R_i):.4g}, R_o={float(self.R_o):.4g}, "
                f"2h={float(2*self.h):.4g}, 2Θ={float(2*self.Theta/mp.pi):.3g}π, "
                f"r0/(2b)={float(self.r_0/(2*self.b)):.4g})")


class MaterialModel:
    """Strategy interface (Phase-0 refactor) supplying the dimensionless ODE
    coefficient constants the exact-edge wave solvers need.

    The wave-variational method's exact-edge ODE has a small set of dimensionless
    coefficients:
      - out-of-plane (Part 1):  T, R   (enter α = 2Tξ²+R,  β = Rξ²(ξ²−2(T+R)))
      - in-plane    (Part 2):  c11, R  (the c̄'11 = c̄11/c̄66 ratio and R)
    plus the effective Poisson ratio ν̂.

    For a HOMOGENEOUS (isotropic) plate these are constants.  For a radially
    GRADED plate (Paper B) they become functions of the radial coordinate; the
    `*_is_graded` flags and the `*_series` hooks let the exact-edge solver switch
    from the constant recurrence to the graded one WITHOUT any other layer
    knowing the difference.  In Phase 0 only the constant path is implemented and
    used, so behaviour is bit-for-bit identical to before; the graded hooks are
    stubs that raise until Paper B fills them in.

    Concrete subclasses must provide: nu_bar, and the OOP constants (T, R) and IP
    constants (c11_eff_bar, R), exactly as the old IsotropicMaterial did.
    """
    #: True when material properties vary along the exact-edge coordinate (FGM).
    radially_graded = False

    def oop_constants(self):
        """Return (T, R, nu_bar) for the out-of-plane (Part 1) ODE."""
        return self.T, self.R, self.nu_bar

    def ip_constants(self):
        """Return (c11_eff_bar, R, nu_bar) for the in-plane (Part 2) ODE."""
        return self.c11_eff_bar, self.R, self.nu_bar

    # ---- graded hooks (Paper B): default to "not graded" ---------------------
    def oop_coeff_series(self, r0_bar, M):
        """Return the radial dependence of the OOP ODE coefficients as power
        series in x=(rr−r0).  Homogeneous default: constants (handled by the
        exact-edge solver's constant path), so this is unused unless
        radially_graded is True."""
        raise NotImplementedError("oop_coeff_series only needed for graded "
                                  "materials (Paper B)")

    def ip_coeff_series(self, r0_bar, M):
        """Return radial power-series coeffs of IP ODE coefficients (c11, R, nu).

        Homogeneous default: unused unless radially_graded is True (same contract
        as oop_coeff_series). Concrete graded materials override this.
        """
        raise NotImplementedError("ip_coeff_series only needed for graded "
                                  "materials (Paper B)")


class IsotropicMaterial(MaterialModel):
    """Homogeneous isotropic plate material (default validation path).

    Sets T=R=1, nu_bar=nu, and c11_eff_bar=2/(1-nu) so the exact-edge ODEs
    reduce to the isotropic Seok–Tiersten special case used by the paper tables.
    """
    radially_graded = False

    def __init__(self, E, nu, rho):
        """E [Pa], nu [-], rho [kg/m^3] — paper/validation defaults use steel-like values."""
        self.E = mpf(E); self.nu = mpf(nu); self.rho = mpf(rho)
        self.c11_bar = self.E/(1-self.nu**2)
        self.c66 = self.E/(2*(1+self.nu))
        self.nu_bar = self.nu
        self.T = mpf('1'); self.R = mpf('1')
        self.c11_eff_bar = mpf('2')/(1-self.nu)     # c̄'11 = c̄11/c̄66
    def __repr__(self):
        return (f"IsotropicMaterial(E={float(self.E):.3g}, ν={float(self.nu):.3g}, "
                f"ρ={float(self.rho):.3g})")


class FGMPlateProperties:
    """Scaffold for physical E(r)/rho(r) power-law FGM profiles.

    Distinct from RadialFGMMaterial, which grades the *dimensionless ODE ratio*
    coefficients (T, R, nu) used by the exact-edge engines. This class is a
    roadmap hook for E/rho-based formulations, not the deployed ODE grader.
    """
    def __init__(self, E_i, E_o, nu, rho_i, rho_o, n_p, geom):
        """Inner/outer E and rho, constant nu, power n_p, and host geometry."""
        self.E_i = mpf(E_i); self.E_o = mpf(E_o); self.nu = mpf(nu)
        self.rho_i = mpf(rho_i); self.rho_o = mpf(rho_o)
        self.n_p = mpf(n_p); self.geom = geom
    def _x(self, r):
        """Normalized radial coordinate (r-R_i)/(R_o-R_i) in [0,1]."""
        return (r-self.geom.R_i)/(self.geom.R_o-self.geom.R_i)
    def E(self, r):
        """Young's modulus at physical radius r (power-law blend)."""
        return self.E_i+(self.E_o-self.E_i)*self._x(r)**self.n_p
    def rho(self, r):
        """Density at physical radius r (power-law blend)."""
        return self.rho_i+(self.rho_o-self.rho_i)*self._x(r)**self.n_p


def _fgm_power_law_series(v_i, v_o, n_p):
    """Exact Taylor coefficients, in the SAME shifted/nondimensional radial
    coordinate x the exact-edge solvers already use (rr = x + r0_bar), of a
    physical power-law profile V(r) = v_i + (v_o - v_i) * ((r-R_i)/(R_o-R_i))^n_p,
    for a NONNEGATIVE INTEGER n_p (binomial expansion -- exact, no truncation
    error; non-integer n_p is rejected rather than silently approximated).

    Independent of R_i/R_o/b: the project's own r0_bar/b nondimensionalization
    already maps r in [R_i,R_o] onto rr in [r0_bar-pi/2, r0_bar+pi/2] for ANY
    annulus geometry (r0_bar = pi*r_0/(2b), b=(R_o-R_i)/2), so
        zeta(x) := (r(x)-R_i)/(R_o-R_i) = 1/2 + x/pi
    identically -- derived directly from r(x) = (x+r0_bar)*2b/pi, not assumed.
    """
    from math import comb
    n_p_int = int(round(float(n_p)))
    if abs(float(n_p) - n_p_int) > 1e-9 or n_p_int < 0:
        raise ValueError(
            "RadialFGMMaterial requires a nonnegative INTEGER grading "
            f"exponent n_p for an exact (non-truncated) series; got {n_p!r}. "
            "A fractional exponent would need a truncated Taylor expansion "
            "with a documented truncation error, not yet implemented.")
    half = mpf(1)/2
    invpi = 1/mp.pi
    zeta_coeffs = [mpf(comb(n_p_int, k)) * half**(n_p_int - k) * invpi**k
                   for k in range(n_p_int + 1)]
    dv = mpf(v_o) - mpf(v_i)
    coeffs = [dv * c for c in zeta_coeffs]
    coeffs[0] += mpf(v_i)
    return coeffs


def _poly_eval_mp(coeffs, x):
    """Horner evaluation of a coefficient list (index = power of x) at mp
    value x. Used for BOTH oop_coeff_series() itself (grading profile ->
    scalar reference value at x=0, for __repr__/compat) and by the solver's
    graded code paths (T(x)/R(x)/nu(x) evaluated pointwise at a node)."""
    v = mpf(0)
    for c in reversed(coeffs):
        v = v * x + c
    return v


class RadialFGMMaterial(MaterialModel):
    """Radially-graded polar-orthotropic OOP material (Paper B / roadmap
    "FGM" item).  Grades the DIMENSIONLESS polar-orthotropic ODE ratio
    coefficients T, R (alpha=2T xi^2+R, beta=R xi^4-2(T+R)xi^2 -- LESSONS
    Sec 8's corrected form, NOT the old "spurious extra R" form) and the
    Poisson-type constant the kernel calls "nu" (= mu_theta = D12/D11),
    each via an independent power-law profile from an inner value to an
    outer value with a shared integer exponent n_p -- the same
    v_i+(v_o-v_i)*x^n_p convention FGMPlateProperties already uses for E/rho.

    SCOPE (explicit, not silently assumed): this grades the ratio-type ODE
    coefficients, which the exact-edge recursion already isolates from the
    frequency-scale constant k_bar (k_bar is evaluated ONCE from c11_bar/
    c66 at the reference radius r0 and held fixed -- NOT graded).  This
    models a material whose polar-orthotropic RATIOS vary radially (e.g. a
    radially-varying degree of anisotropy or Poisson ratio) with a FIXED
    frequency-normalization scale -- it does NOT model the more common
    "single common profile" FGM assumption (E_r(r)=E_r0*f(r) etc. for ALL
    moduli with the SAME f(r)), which leaves T/R/mu_theta constant (they
    are ratios; f(r) cancels) but DOES vary k_bar/c11_bar with r -- that is
    a structurally similar but separate extension (same convolution
    machinery, applied to the Omega^2 term instead of the T,R terms) not
    attempted here; see FUTURE_WORK.md's FGM item for this scoping note.
    2026-08-02: this Omega^2-term grading IS now available, via the
    OPTIONAL kappa_i/kappa_o constructor args below and
    oop_kappa_series()/OutOfPlaneSolver's K_omega convolution -- it is a
    SEPARATE profile from T/R/nu (not tied to the "single common profile"
    FGM assumption either; kappa(x) is just whatever local ratio the caller
    supplies). kappa_i=kappa_o=None (the default) is an exact no-op for
    every pre-existing call site. See GrokCode/KGradingDerivation.txt for
    the derivation and FUTURE_WORK.md's FGM item for the scoping history.
    Sandbox self-tested only as of this edit -- NOT yet cluster-confirmed;
    see probe_kbar_grading_validation.py.

    Isotropic reduction: T_i=T_o=R_i=R_o=1, nu_i=nu_o=nu recovers the
    IsotropicMaterial constants exactly, and radially_graded=True still
    routes through the (more expensive) graded recursion -- for a true
    isotropic no-op, use IsotropicMaterial itself, not this class with a
    flat profile (the graded recursion is a correct but needlessly O(M^2)
    superset of the constant-coefficient O(M) path for that case).
    """
    radially_graded = True
    validated = False

    def __init__(self, T_i, T_o, R_i, R_o, nu_i, nu_o, n_p,
                 c11_bar, c66, rho, geom, c11_i=None, c11_o=None,
                 kappa_i=None, kappa_o=None):
        self.T_i = mpf(T_i); self.T_o = mpf(T_o)
        self.R_i = mpf(R_i); self.R_o = mpf(R_o)
        self.nu_i = mpf(nu_i); self.nu_o = mpf(nu_o)
        self.n_p = mpf(n_p); self.geom = geom
        self.c11_bar = mpf(c11_bar); self.c66 = mpf(c66); self.rho = mpf(rho)
        # In-plane (Part 2) c11 grading -- added 2026-07-30, LESSONS_LEARNED
        # Sec 54.4. c11_i/c11_o are OPTIONAL and default to the ungraded
        # reference value (c11_bar/c66) for BOTH endpoints when omitted, so
        # every pre-existing call site (OOP-only grading) is unaffected --
        # v_i==v_o makes _fgm_power_law_series collapse to a flat series
        # regardless of n_p, byte-identical to the old scalar c11_eff_bar.
        _c11_eff_ref = self.c11_bar / self.c66     # c̄'11 = c̄11/c̄66, ungraded ref
        self.c11_i = mpf(c11_i) if c11_i is not None else _c11_eff_ref
        self.c11_o = mpf(c11_o) if c11_o is not None else _c11_eff_ref
        # Omega^2-term ("K_omega") grading -- added 2026-08-02, see class
        # docstring above and GrokCode/KGradingDerivation.txt. OPTIONAL,
        # defaults to a flat kappa(x)=1 profile (exact no-op) when omitted.
        self.kappa_i = mpf(kappa_i) if kappa_i is not None else mpf(1)
        self.kappa_o = mpf(kappa_o) if kappa_o is not None else mpf(1)
        # Reference-radius (x=0, r=r_0) scalar values -- used by any caller
        # that only wants oop_constants()/ip_constants() (e.g. a quick
        # isotropic-reduction sanity check), and as this object's own
        # nu_bar/T/R/c11_eff_bar attributes for __repr__ and k_bar-style
        # code that expects a plain scalar.
        self._T_series = _fgm_power_law_series(self.T_i, self.T_o, self.n_p)
        self._R_series = _fgm_power_law_series(self.R_i, self.R_o, self.n_p)
        self._nu_series = _fgm_power_law_series(self.nu_i, self.nu_o, self.n_p)
        self._c11_series = _fgm_power_law_series(self.c11_i, self.c11_o, self.n_p)
        self._kappa_series = _fgm_power_law_series(self.kappa_i, self.kappa_o, self.n_p)
        self.T = _poly_eval_mp(self._T_series, mpf(0))
        self.R = _poly_eval_mp(self._R_series, mpf(0))
        self.nu_bar = _poly_eval_mp(self._nu_series, mpf(0))
        self.c11_eff_bar = _poly_eval_mp(self._c11_series, mpf(0))

    def oop_constants(self):
        """Reference-radius (r=r_0) scalar (T, R, nu_bar) -- for callers not
        yet grading-aware (e.g. a quick isotropic-reduction check). The
        actual graded solver path uses oop_coeff_series() instead, selected
        by OutOfPlaneSolver.__init__ via self.radially_graded."""
        return self.T, self.R, self.nu_bar

    def ip_constants(self):
        """Reference-radius (r=r_0) scalar (c11_eff_bar, R, nu_bar) -- for
        callers not yet grading-aware (e.g. a quick isotropic-reduction
        check). The actual graded solver path uses ip_coeff_series()
        instead, selected by InPlaneSolver.__init__ via self.radially_graded.
        Implemented 2026-07-30 (in-plane grading, LESSONS_LEARNED.md
        Sec 54.4, derivation self-tested to 1e-24 on a non-trivial profile)
        -- previously raised NotImplementedError when only the OOP path
        existed; see that history note in git blame / LESSONS_LEARNED if
        curious why this looked unfinished before."""
        return self.c11_eff_bar, self.R, self.nu_bar

    def ip_coeff_series(self, r0_bar, M):
        """Return (c11_coeffs, R_coeffs, nu_coeffs) power-series coefficient
        lists in x (rr=x+r0_bar) for the in-plane (Part 2) graded recursion
        -- mirrors oop_coeff_series() exactly (same geometry-independent
        zeta(x)=1/2+x/pi profile, both arguments accepted for interface
        symmetry only). See InPlaneSolver._series_mp's graded branch
        (core_solvers.py) for the linear-convolution recursion this feeds,
        and dispersion.py::_Part2Fast.series() for the float64 mirror used
        by the fast root search."""
        return list(self._c11_series), list(self._R_series), list(self._nu_series)

    def oop_coeff_series(self, r0_bar, M):
        """Return (T_coeffs, R_coeffs, nu_coeffs) power-series coefficient
        lists in x (rr = x+r0_bar), independent of r0_bar/M (the profile is
        already expressed in x via _fgm_power_law_series's own geometry-
        independent zeta(x)=1/2+x/pi mapping) -- both arguments accepted
        for interface symmetry with the (not yet needed) M-truncated
        non-integer-exponent case, unused here since integer-n_p series are
        exact and already finite."""
        return list(self._T_series), list(self._R_series), list(self._nu_series)

    def oop_kappa_series(self, r0_bar, M):
        """Return the Taylor-coefficient list (in x, rr=x+r0_bar) for
        kappa(x) -- the LOCAL grading of the Omega^2 (inertia/bending-
        stiffness) coefficient in the OOP radial ODE, independent of and
        in ADDITION to the T(x)/R(x)/nu(x) ratio grading above. Added
        2026-08-02. kappa_i=kappa_o=1 (the default when neither is passed)
        collapses to the constant [1] series, making
        OutOfPlaneSolver._series_mp's K_omega convolution an exact
        identity (v_n = a_n) -- a byte-for-byte no-op for every
        pre-existing RadialFGMMaterial call site that doesn't pass
        kappa_i/kappa_o. See GrokCode/KGradingDerivation.txt for the
        derivation (causal recursion, same structure as the T/R
        convolution) and core_solvers.py::OutOfPlaneSolver._series_mp's
        own K_omega/_kconv comment for how this is consumed. NOT yet
        given an in-plane (IP) counterpart -- IP's frequency term has a
        different structure (Seok & Tiersten Part 2 Eq. 48 vs Part 1's
        Eq. 43) and would need its own derivation, not this same
        recursion reused."""
        return list(self._kappa_series)

    def __repr__(self):
        return (f"RadialFGMMaterial(T:{float(self.T_i):.3g}->{float(self.T_o):.3g}, "
                f"R:{float(self.R_i):.3g}->{float(self.R_o):.3g}, "
                f"nu:{float(self.nu_i):.3g}->{float(self.nu_o):.3g}, "
                f"n_p={float(self.n_p):.3g})")


# ══════════════════════════════════════════════════════════════════════════════
#  SECTION 2 — SHARED UTILITIES
# ══════════════════════════════════════════════════════════════════════════════

class _ScalarGeom:
    """Minimal geometry object reconstructible from plain scalars.

    Field order must match the g_scalars tuples packed by solvers and
    unpacked by every ProcessPoolExecutor worker via _make_geom_mat.
    """
    def __init__(self, r0_bar, Theta, b, r_0, R_i, R_o, h):
        from mpmath import mpf, mp
        self.r0_bar = mpf(r0_bar)
        self.Theta  = mpf(Theta)
        self.b      = mpf(b)
        self.r_0    = mpf(r_0)
        self.R_i    = mpf(R_i)
        self.R_o    = mpf(R_o)
        self.h      = mpf(h)


class _ScalarMat(MaterialModel):
    """Minimal material object reconstructible from plain scalars.

    Subclasses MaterialModel so the Phase-0 refactored solver constructors
    (which pull constants via mat.oop_constants()/ip_constants() and read
    mat.radially_graded) work in the multiprocessing worker reconstruction
    path exactly as IsotropicMaterial does in-process.  The inherited
    oop_constants() returns (T, R, nu_bar) and ip_constants() returns
    (c11_eff_bar, R, nu_bar) -- all attributes this class already sets."""
    radially_graded = False
    def __init__(self, nu_bar, T, R, c11_bar, c11_eff_bar, c66, rho):
        from mpmath import mpf
        self.nu_bar      = mpf(nu_bar)
        self.T           = mpf(T)
        self.R           = mpf(R)
        self.c11_bar     = mpf(c11_bar)
        self.c11_eff_bar = mpf(c11_eff_bar)
        self.c66         = mpf(c66)
        self.rho         = mpf(rho)


def _make_geom_mat(g_scalars, m_scalars, dps):
    """Reconstruct geometry/material and set mp.dps.  Returns (geom, mat).

    Used on every process-pool worker path: sets the worker's mpmath precision
    then builds _ScalarGeom/_ScalarMat from plain Python floats (picklable).
    """
    from mpmath import mp
    mp.dps = dps
    return _ScalarGeom(*g_scalars), _ScalarMat(*m_scalars)


# Per-process guard: surface the FIRST exception each worker hits, then stay
# quiet.  Workers otherwise swallow every failure into a NaN sentinel, which
# (Phase-0 _ScalarMat regression) hid a one-line constructor-contract break
# behind 'σ_min sweep produced no finite values' across 29 geometries.  This
# keeps the NaN return identical but makes the next such bug obvious in the log.
_WORKER_EXC_SURFACED = False

def make_geometry(r0_over_2b, two_Theta_over_pi, h=0.04, b=2.0):
    """Construct PlateGeometry for a given r0/(2b) and 2*Theta/pi.

    Uses b = (R_o - R_i)/2 (defaults to 2.0 to match the paper's calibration
    geometry of R_i=3, R_o=7 when r0/(2b)=1.25).  R_i and R_o are derived:
        r_0 = r0_over_2b * 2b
        R_i = r_0 - b,  R_o = r_0 + b
    Theta = two_Theta_over_pi * pi / 2.
    h is the half-thickness (paper uses h=0.04 → full thickness 0.08).
    """
    r_0 = r0_over_2b * 2 * b
    R_i = r_0 - b
    R_o = r_0 + b
    Theta = two_Theta_over_pi * float(mp.pi) / 2.0
    return PlateGeometry(R_i=R_i, R_o=R_o, h=h, Theta=Theta)


# ── Per-geometry scan configuration ──────────────────────────────────────────
# Scan ranges come from the ξ=0/ζ=0 cut-off frequencies (no paper tables).
# Cut-off physics and the proximity-filter window: LESSONS_LEARNED.md §Scan.

class OrthotropicMaterial(MaterialModel):
    """Polar-orthotropic material SCAFFOLD (framework hook for the orthotropic
    extension; see Shi/Lv 2016 for benchmark targets).

    The exact-edge OOP ODE already carries polar-orthotropy constants (T, R);
    isotropic is the special case T=R=1.  This class supplies the isotropic-
    reduction constants exactly, so a free-free *isotropic* run can be driven
    through the orthotropic code path as a sanity check.  For TRUE orthotropy
    (E_r != E_theta) the Seok & Tiersten definitions of T and R in terms of the
    polar rigidities D_r, D_theta, D_rtheta, D_k must be filled in below before
    the numbers can be trusted -- that derivation is deliberately left as a
    flagged stub rather than guessed.
    """
    radially_graded = False
    validated = False

    def __init__(self, E_r, E_theta, nu_r, G_rtheta, rho):
        self.E_r = mpf(E_r); self.E_theta = mpf(E_theta)
        self.nu_r = mpf(nu_r); self.G_rtheta = mpf(G_rtheta); self.rho = mpf(rho)
        # reciprocal Poisson ratio from the orthotropy relation nu_theta/E_theta
        # = nu_r/E_r
        self.nu_theta = self.nu_r * self.E_theta / self.E_r
        self._isotropic = (abs(self.E_r - self.E_theta) < mpf('1e-9') * self.E_r)
        # Radial reference modulus drives the c11_bar / c66 normalization, in the
        # same way IsotropicMaterial uses E.
        self.c11_bar = self.E_r / (1 - self.nu_r * self.nu_theta)
        self.c66 = self.G_rtheta
        self.nu_bar = self.nu_r
        self.c11_eff_bar = self.c11_bar / self.c66
        # ── Polar-orthotropic OOP ODE constants (derived; LESSONS §12) ──────────
        # In terms of the polar flexural rigidities (h³/12 cancels in all ratios):
        #   D11 = E_r/(1-mu_r*mu_th), D22 = E_th/(1-mu_r*mu_th),
        #   D12 = mu_r*E_th/(1-mu_r*mu_th), D66 = G_rth.   The exact-edge OOP
        #   ODE/BCs need:
        #     R    = D22/D11 = E_th/E_r           (alpha=2T xi^2+R, beta=R xi^4-2(T+R)xi^2)
        #     T    = (D12+2D66)/D11
        #     mu_th= D12/D11 = mu_r*E_th/E_r       (the Poisson constant in every
        #            OOP moment/shear -- what the kernel calls "nu"; oop_constants
        #            returns THIS, not mu_r).
        # Verified: isotropic reduction (E_r=E_th, G=E/2(1+nu)) => R=T=1, mu_th=nu,
        # and the kernel reduces bit-for-bit to IsotropicMaterial.
        self.R = self.E_theta / self.E_r
        self.mu_theta = self.nu_theta                       # = D12/D11
        self.T = self.mu_theta + 2*self.G_rtheta*(1 - self.nu_r*self.nu_theta)/self.E_r
        # Whether this material has been validated against an external benchmark
        # (Shi/Lv 2016 FFFF).  False => run_shi_validation prints a loud caveat.
        self.validated = bool(self._isotropic)

    def oop_constants(self):
        """Return (T, R, mu_theta) for the OOP exact-edge ODE.

        The Poisson slot is mu_theta (= D12/D11), NOT mu_r. Equal for the
        isotropic reduction, so that path is unchanged.
        """
        return self.T, self.R, self.mu_theta

    def ip_constants(self):
        """Return (c11_eff_bar, R, mu_theta) for the IP exact-edge ODE.

        FIX (2026-07-26, cluster-validated against Wang/Liang/Yao/Zhang 2016
        Table 4 FFFF, 8/8 <1.4%): base MaterialModel.ip_constants returns
        nu_bar (= nu_r), wrong for orthotropy — same nu_r vs mu_theta mistake
        oop_constants already fixes on the OOP side. Returns mu_theta
        (= D12/D11 = nu_r*E_theta/E_r). Isotropic reduction is a no-op.
        """
        return self.c11_eff_bar, self.R, self.mu_theta

    def __repr__(self):
        return (f"OrthotropicMaterial(E_r={float(self.E_r):.3g}, "
                f"E_theta={float(self.E_theta):.3g}, nu_r={float(self.nu_r):.3g}, "
                f"G={float(self.G_rtheta):.3g}, iso_reduction={self._isotropic})")


def _omega_lit(fq, geom, mat, w_bar, k_bar, part):
    """Convert a dimensionless Seok-Tiersten frequency to the common literature
    parameter  Omega_lit = omega * R_o^2 * sqrt(rho*H / D),  H = 2h (full
    thickness), D = E*H^3/[12(1-nu^2)] = c11_bar*H^3/12.  part 1 = OOP (k_bar
    scaling), part 2 = IP.  Returns (f_Hz, Omega_lit)."""
    import math
    if part == 1:
        f_hz = float(fq) * w_bar / (2 * math.pi * k_bar)
    else:
        f_hz = float(fq) * w_bar / (2 * math.pi)
    omega = 2 * math.pi * f_hz
    H = 2 * float(geom.h)
    D = float(mat.c11_bar) * H ** 3 / 12.0
    Ro = float(geom.R_o)
    Om_lit = omega * Ro ** 2 * math.sqrt(float(mat.rho) * H / D)
    return f_hz, Om_lit


def _kbar_wbar(geom, mat):
    """Return (k_bar, w_bar) frequency-scale pair for literature conversions.

    k_bar scales OOP dimensionless Omega; w_bar is the reference radian
    frequency used when converting to Hz / Omega_lit (see _omega_lit).
    """
    k_bar = float((mpf(2) * geom.b / mp.pi) *
                  mp.sqrt(3 * mat.c66 / (geom.h ** 2 * mat.c11_bar)))
    w_bar = float((mp.pi / (2 * geom.b)) * mp.sqrt(mat.c66 / mat.rho))
    return k_bar, w_bar


def _r0_2b_for_ratio(eps):
    """r0/(2b) token (make_geometry uses b=2) that yields R_i/R_o = eps.

    make_geometry: R_i = 4·t − 2, R_o = 4·t + 2  ⇒  (4t−2)/(4t+2) = eps
                   ⇒  t = (1+eps)/(2(1−eps)).   eps→0 ⇒ t→0.5 (solid limit)."""
    return (1.0 + eps) / (2.0 * (1.0 - eps))


def _material_from_env():
    """Build the run material from MAT_E / MAT_NU / MAT_RHO (defaults = the
    ν=0.35 validation material).  Isolated + module-level so it is unit-testable
    and so a literature-matched run (LESSONS §13) needs no source edit."""
    E   = float(os.environ.get("MAT_E",   "210e9"))
    nu  = float(os.environ.get("MAT_NU",  "0.35"))
    rho = float(os.environ.get("MAT_RHO", "7800.0"))
    m = IsotropicMaterial(E=E, nu=nu, rho=rho)
    if abs(nu - 0.35) > 1e-12 or abs(E - 210e9) > 1.0 or abs(rho - 7800.0) > 1e-9:
        print(f"  NOTE: NON-DEFAULT material from env (E={E:.4g}, nu={nu:.4g}, "
              f"rho={rho:.4g}).  The 29-mode paper validation and the cantilever\n"
              f"        spot-check reference numbers assume the DEFAULT ν=0.35 — "
              f"treat their OK/DRIFT labels as informational for this run.")
    return m


def _radius_ratio(r0_2b):
    """R_i/R_o for a geometry token, using make_geometry's b=2 calibration."""
    g = make_geometry(r0_2b, 1.0)
    return float(g.R_i) / float(g.R_o)


