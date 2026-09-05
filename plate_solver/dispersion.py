# -*- coding: utf-8 -*-
"""
plate_solver.dispersion -- the radial/Cartesian dispersion ("exact-edge")
math: cut-off-frequency Bessel determinants (annular), the Frobenius-series
ODE engines ExactEdgeSolver/_Part1Fast/_Part2Fast (annular), and the
Cartesian counterparts RectangularCartesianOOP/IP (rectangular, ported from
rect_int.py).

BUG-DOC ITEM (mechanical cleanup, not a physics change): in
cutoff_frequencies_part1/2, T/R/nu/nu_bar/c11_eff are now cast to float ONCE
in the parent scope before the brentq/np.linspace scan, instead of being
re-cast to float on every one of the thousands of _bessel_det_* calls each
scan makes. float(float(x)) == float(x) bit-for-bit, so this changes nothing
numerically; it only removes redundant casts from a tight loop.

Extracted verbatim otherwise (line-range provenance in LESSONS_LEARNED
Sec. 22); no numeric behaviour changed, SOLVER_VERSION not bumped.
"""
from __future__ import annotations
import os, cmath as _cmath
import numpy as np
from math import comb
from mpmath import mp, mpf, mpc, matrix
from scipy.special import jv, yv, jvp, yvp
from scipy.optimize import brentq

from .config import MP_DPS, PI

def _bessel_det_p1_F(Om, r0b, nu, T=1.0, R=1.0):
    """2×2 determinant for Part 1 F-branch cut-off (ξ=0).

    At ξ=0 the Part 1 governing ODE (eq.41 in Seok-Tiersten Part 1) reduces.
    The F-branch solution (from Part 2 annular paper eq.40) satisfies:

        d²F/dx² + (1/x)dF/dx + (1 - R/x²)F = 0
        x = Ω̃(r̃+r̄₀)/√T

    This is a Bessel equation of order √R.  For isotropic T=R=1, the order is 1.

    The arc BC (paper Part 2 annular eq.43) at ξ=0:
        dF/dr̃ + ν̂·F/(r̃+r̄₀) = 0   at r̃ = ±π/2

    The 2×2 system from both arcs must have zero determinant.
    """
    from scipy.special import jv, yv, jvp, yvp
    hp = np.pi / 2.0
    # Bessel argument at inner (r̃=-π/2) and outer (r̃=+π/2) arcs
    # NOTE: T, R are expected pre-cast to float by the caller (cutoff_
    # frequencies_part1 casts once, outside the thousands-of-calls scan loop
    # brentq/np.linspace drive here); float(float(x)) == float(x), so this is
    # a pure perf hoist with zero numeric effect.
    beta = Om / np.sqrt(T)              # radial wave number scaling
    p    = np.sqrt(R)                   # Bessel order = sqrt(R) (from paper Part 2 eq.40)
    u_i  = beta * (r0b - hp)           # argument at inner arc
    u_o  = beta * (r0b + hp)           # argument at outer arc

    # BC: dF/dr̃ + ν̂·F/(r̃+r̄₀) = 0
    # dF/dr̃ = β·dF/du  (chain rule)
    # Row 0 (inner arc, r̃=-π/2, rr = r̄₀-π/2):
    rr_i = r0b - hp
    bc_Ji = beta * jvp(p, u_i, 1) + nu * jv(p, u_i) / rr_i
    bc_Yi = beta * yvp(p, u_i, 1) + nu * yv(p, u_i) / rr_i
    # Row 1 (outer arc, r̃=+π/2, rr = r̄₀+π/2):
    rr_o = r0b + hp
    bc_Jo = beta * jvp(p, u_o, 1) + nu * jv(p, u_o) / rr_o
    bc_Yo = beta * yvp(p, u_o, 1) + nu * yv(p, u_o) / rr_o

    return bc_Ji * bc_Yo - bc_Yi * bc_Jo


def _bessel_det_p1_G(Om, r0b, nu, T=1.0, R=1.0):
    """2×2 determinant for Part 1 G-branch cut-off (ξ=0).

    The G-branch solution is:
        w*(r̃) = D₁ J_1(β₂(r̃+r̄₀)) + D₂ Y_1(β₂(r̃+r̄₀))

    where β₂² = Ω̃²/R (from eq.41 of Part 1 at ξ=0).
    The arc BC (eq.44 in Part 2 annular, or Part 1 Section 3) is:
        dG/dr* - G/(r*+r̄₀) = 0   at r* = ±π/2

    This is precisely the Bessel equation BC for order 1 Bessel functions
    (J'_1(u) = J_0(u) - J_1(u)/u, so the condition dG/du·β₂ - G/rr = 0
    reduces neatly).
    """
    from scipy.special import jv, yv, jvp, yvp
    hp = np.pi / 2.0
    beta = Om / np.sqrt(R)               # R pre-cast to float by the caller
    p    = 1.0                          # always order 1 for G-branch
    u_i  = beta * (r0b - hp)
    u_o  = beta * (r0b + hp)

    rr_i = r0b - hp
    rr_o = r0b + hp
    bc_Ji = beta * jvp(p, u_i, 1) - jv(p, u_i) / rr_i
    bc_Yi = beta * yvp(p, u_i, 1) - yv(p, u_i) / rr_i
    bc_Jo = beta * jvp(p, u_o, 1) - jv(p, u_o) / rr_o
    bc_Yo = beta * yvp(p, u_o, 1) - yv(p, u_o) / rr_o

    return bc_Ji * bc_Yo - bc_Yi * bc_Jo


def cutoff_frequencies_part1(r0b, nu, T=1.0, R=1.0, n_co=8,
                              Om_lo=0.001, Om_hi=3.0, n_grid=3000):
    """Compute Part 1 (out-of-plane) cut-off frequencies by scanning ξ=0.

    At ξ=0 the 4×4 dispersion determinant factors into two independent 2×2
    Bessel determinants (F-branch and G-branch).  This function scans both
    over a fine Ω̃ grid and collects the n_co smallest sign-change roots.

    Parameters
    ----------
    r0b    : float  —  r̄₀ = π·r₀/(2b), dimensionless mid-radius
    nu     : float  —  effective Poisson ratio ν̂
    T, R   : float  —  polar orthotropy constants (isotropic: T=R=1)
    n_co   : int    —  how many cut-off values to return
    Om_lo  : float  —  lower bound of scan (avoid 0 to skip trivial root)
    Om_hi  : float  —  upper bound of scan
    n_grid : int    —  number of grid points for initial sign-change search

    Returns
    -------
    list of float — sorted cut-off Ω̃ values (smallest first)
    """
    from scipy.optimize import brentq
    # Cast once here (not inside _bessel_det_p1_F/_G, which brentq/the grid
    # scan below call thousands of times) -- perf hoist, no numeric change.
    r0b, nu, T, R = float(r0b), float(nu), float(T), float(R)
    Oms = np.linspace(Om_lo, Om_hi, n_grid)
    cutoffs = []

    for det_fn in (_bessel_det_p1_F, _bessel_det_p1_G):
        vals = np.array([det_fn(Om, r0b, nu, T, R) for Om in Oms])
        for k in range(len(vals) - 1):
            if np.isfinite(vals[k]) and np.isfinite(vals[k+1]):
                if vals[k] * vals[k+1] < 0:
                    try:
                        root = brentq(det_fn, Oms[k], Oms[k+1],
                                      args=(r0b, nu, T, R),
                                      xtol=1e-8, rtol=1e-8)
                        cutoffs.append(root)
                    except ValueError:
                        pass

    cutoffs = sorted(set(round(c, 6) for c in cutoffs))
    return cutoffs[:n_co]


def _bessel_det_p2_F(Om, r0b, nu_bar, c11_eff, R=1.0):
    """2×2 Bessel determinant for Part 2 F-branch cut-off (ζ=0).

    At ζ=0 the Part 2 coupled equations (Seok-Tiersten Part 2, eqs.40-41)
    decouple.  The F-branch satisfies:
        d²F/dx² + (1/x)dF/dx + (1 - R/x²)F = 0
    with  x = Ω̄(r*+r̄₀)/√c̄'₁₁  (eq.42).
    This is a Bessel equation of order √R.

    The arc BC at ζ=0 (eq.43 of Part 2 annular paper) is:
        dF/dr* + ν̂·F/(r*+r̄₀) = 0   at r* = ±π/2

    With x = Ω̄·rr/√c̄'₁₁,  dF/dr* = (Ω̄/√c̄'₁₁)·dF/dx, and the BC becomes:
        (Ω̄/√c̄'₁₁)·F'(x) + ν̂·F(x)/rr = 0
    """
    from scipy.special import jv, yv, jvp, yvp
    hp = np.pi / 2.0
    sc = Om / np.sqrt(c11_eff)   # scaling factor Ω̄/√c̄'₁₁
    p  = np.sqrt(R)               # Bessel order; R pre-cast to float by caller

    rr_i = r0b - hp;  u_i = sc * rr_i
    rr_o = r0b + hp;  u_o = sc * rr_o

    bc_Ji = sc * jvp(p, u_i, 1) + nu_bar * jv(p, u_i) / rr_i
    bc_Yi = sc * yvp(p, u_i, 1) + nu_bar * yv(p, u_i) / rr_i
    bc_Jo = sc * jvp(p, u_o, 1) + nu_bar * jv(p, u_o) / rr_o
    bc_Yo = sc * yvp(p, u_o, 1) + nu_bar * yv(p, u_o) / rr_o

    return bc_Ji * bc_Yo - bc_Yi * bc_Jo


def _bessel_det_p2_G(Om, r0b, nu_bar, c11_eff, R=1.0):
    """2×2 Bessel determinant for Part 2 G-branch cut-off (ζ=0).

    The G-branch satisfies:
        d²G/dZ² + (1/Z)dG/dZ + (1 - 1/Z²)G = 0
    with  Z = Ω̄(r*+r̄₀)  (order 1 Bessel equation, eq.41/42).

    The arc BC at ζ=0 (eq.44 of Part 2 annular paper) is:
        dG/dr* - G/(r*+r̄₀) = 0   at r* = ±π/2
    """
    from scipy.special import jv, yv, jvp, yvp
    hp = np.pi / 2.0
    p  = 1.0   # always order 1 for G-branch

    rr_i = r0b - hp;  u_i = Om * rr_i
    rr_o = r0b + hp;  u_o = Om * rr_o

    bc_Ji = Om * jvp(p, u_i, 1) - jv(p, u_i) / rr_i
    bc_Yi = Om * yvp(p, u_i, 1) - yv(p, u_i) / rr_i
    bc_Jo = Om * jvp(p, u_o, 1) - jv(p, u_o) / rr_o
    bc_Yo = Om * yvp(p, u_o, 1) - yv(p, u_o) / rr_o

    return bc_Ji * bc_Yo - bc_Yi * bc_Jo


def cutoff_frequencies_part2(r0b, nu_bar, c11_eff, R=1.0, n_co=8,
                              Om_lo=0.001, Om_hi=3.0, n_grid=3000):
    """Compute Part 2 (in-plane) cut-off frequencies by scanning ζ=0.

    At ζ=0 the Part 2 coupled equations decouple into two independent Bessel
    equations (Seok-Tiersten Part 2, eqs.40-44).  This function scans both the
    F-branch and G-branch 2×2 Bessel determinants over a fine Ω̄ grid.

    Parameters
    ----------
    r0b     : float  —  r̄₀ = π·r₀/(2b)
    nu_bar  : float  —  effective Poisson ratio ν̂
    c11_eff : float  —  c̄'₁₁ = c̄₁₁/c̄₆₆  (= 2/(1-ν̂) for isotropic)
    R       : float  —  polar orthotropy constant (isotropic: R=1)
    n_co    : int    —  number of cut-off values to return
    Om_lo   : float  —  lower bound of scan
    Om_hi   : float  —  upper bound of scan
    n_grid  : int    —  grid density for sign-change search

    Returns
    -------
    list of float — sorted cut-off Ω̄ values (smallest first)
    """
    from scipy.optimize import brentq
    # Cast once here (not inside _bessel_det_p2_F/_G, which brentq/the grid
    # scan below call thousands of times) -- perf hoist, no numeric change.
    r0b, nu_bar, c11_eff, R = float(r0b), float(nu_bar), float(c11_eff), float(R)
    Oms = np.linspace(Om_lo, Om_hi, n_grid)
    cutoffs = []

    for det_fn in (_bessel_det_p2_F, _bessel_det_p2_G):
        vals = np.array([det_fn(Om, r0b, nu_bar, c11_eff, R) for Om in Oms])
        for k in range(len(vals) - 1):
            if np.isfinite(vals[k]) and np.isfinite(vals[k+1]):
                if vals[k] * vals[k+1] < 0:
                    try:
                        root = brentq(det_fn, Oms[k], Oms[k+1],
                                      args=(r0b, nu_bar, c11_eff, R),
                                      xtol=1e-8, rtol=1e-8)
                        cutoffs.append(root)
                    except ValueError:
                        pass

    cutoffs = sorted(set(round(c, 6) for c in cutoffs))
    return cutoffs[:n_co]


# ══════════════════════════════════════════════════════════════════════════════
#  SECTION 3 — FAST (float64) DISPERSION ENGINES  [branch finding/tracking]
# ══════════════════════════════════════════════════════════════════════════════

class ExactEdgeSolver:
    """Strategy interface (Phase-0 refactor) for the EXACT-edge wave solution —
    the part of the method that satisfies the governing PDE and the two "exact"
    edges (the arcs r=R_i, R_o for an annulus; the two side edges for a
    rectangle) exactly, yielding the dispersion relation.

    This is the float64 "engine" consumed by `full_search`/`track`; it must
    provide:
      - series(wavenumber, Om)         : basis solutions of the exact-edge ODE
      - Lmat(wavenumber, sols)         : the 4×4 (or n×n) arc-BC matrix
      - det(wavenumber, Om)            : determinant of Lmat (dispersion function)
      - newton(z0, Om, ...)            : root polish on det

    Concrete strategies:
      - `_Part1Fast` (alias `AnnulusRadialOOP`) : annulus, out-of-plane (Part 1)
      - `_Part2Fast` (alias `AnnulusRadialIP`)  : annulus, in-plane    (Part 2)
      - (future) RectangularCartesian{OOP,IP}   : rectangle  — Paper A
      - the graded recurrence is selected INSIDE the strategy via the material's
        `radially_graded` flag — Paper B

    `geometry_kind` and `motion` are descriptive tags so the assembler/figures
    can label output; `radially_graded` records whether the ODE coefficients
    vary along the exact-edge coordinate.
    """
    geometry_kind = "abstract"
    motion = "abstract"
    radially_graded = False

    def series(self, wavenumber, Om):
        """Return exact-edge basis solutions at (wavenumber, Om)."""
        raise NotImplementedError

    def Lmat(self, wavenumber, sols):
        """Build the arc (or free-edge) BC matrix from series() output."""
        raise NotImplementedError

    def det(self, wavenumber, Om):
        """Dispersion function: det(Lmat(series(...))); full_search tracks zeros."""
        raise NotImplementedError

    def newton(self, z0, Om, itmax=40, tol=1e-13):
        """Polish a complex root of det(·, Om) starting from seed z0."""
        raise NotImplementedError


def _coef_eval_f64(coef, x):
    """float64 counterpart of core_solvers._coef_eval -- evaluate a
    scalar-or-series ODE coefficient at x. Non-graded materials pass a
    plain float (returned unchanged, zero behavior change); graded
    materials (Paper B) pass a list of float Taylor coefficients, Horner-
    evaluated pointwise. See core_solvers._coef_eval's own docstring for
    why this is exact (pointwise substitution, no r-derivative of T/R/nu
    ever appears in the arc/theta-edge algebra)."""
    if isinstance(coef, list):
        v = 0.0
        for c in reversed(coef):
            v = v * x + c
        return v
    return coef


class _Part1Fast(ExactEdgeSolver):
    """float64 Frobenius series + 4×4 arc-BC determinant for annular OOP (Part 1).

    Alias: AnnulusRadialOOP. Consumed by full_search/track for branch finding.
    Constant-coefficient path is the homogeneous annulus; when radially_graded
    is True, nu/T/R may be Taylor-coefficient lists and series() uses the graded
    recurrence (Paper B / LESSONS FGM notes).
    """
    geometry_kind = "annulus"
    motion = "out_of_plane"

    def __init__(self, r0b, nu, T=1.0, R=1.0, M=80, radially_graded=False, K_omega=1.0):
        """r0b=r0_bar; nu,T,R scalar or Taylor lists; M=Frobenius truncation.
        K_omega: OPTIONAL Omega^2-term grading (scalar or Taylor list),
        added 2026-08-02 -- see core_solvers.py's OutOfPlaneSolver.__init__
        and _series_mp's own K_omega/_kconv comment. Default 1.0 is an
        exact no-op (matches the graded-branch identity _kconv==a[n])."""
        self.r0b = float(r0b)
        self.radially_graded = bool(radially_graded)
        # Graded materials pass nu/T/R as LISTS of float Taylor coefficients
        # (Paper B) -- keep them as-is; non-graded materials pass plain
        # scalars, force-cast exactly as before (byte-for-byte unchanged).
        if self.radially_graded:
            self.nu = [float(c) for c in nu] if isinstance(nu, (list, tuple)) else float(nu)
            self.T = [float(c) for c in T] if isinstance(T, (list, tuple)) else float(T)
            self.R = [float(c) for c in R] if isinstance(R, (list, tuple)) else float(R)
            self.K_omega = ([float(c) for c in K_omega]
                             if isinstance(K_omega, (list, tuple)) else float(K_omega))
        else:
            self.nu = float(nu); self.T = float(T); self.R = float(R)
        self.M = M

    def series(self, xi, Om):
        """Corrected ordinary-point power series: 4 solutions with unit IC."""
        xi = complex(xi); Ot2 = complex(Om)**2
        T, R, r0 = self.T, self.R, self.r0b
        M = self.M
        if not isinstance(T, list):
            # ---- EXACT existing constant-coefficient path, byte-for-byte
            # unchanged ----
            al = 2*T*xi**2 + R
            # β: correct polar-ORTHOTROPIC form  R·ξ⁴ − 2(T+R)·ξ²  (LESSONS §12);
            # byte-identical to the old R·ξ²(ξ²−2(T+R)) at R=1 (isotropic).  The old
            # form carried a spurious extra R on the ξ² term when R≠1.
            be = xi**2*(R*xi**2 - 2*(T+R))
            sols = []
            for q in range(4):
                a = np.zeros(M+5, complex); a[q] = 1.0
                for s in range(4, M+5):
                    rhs = 0.0+0.0j
                    for j in range(1, 5):
                        m = s-j
                        if m >= 0:
                            ff = m*(m-1)*(m-2)*(m-3)
                            if ff: rhs += comb(4, j)*r0**(4-j)*ff*a[m]
                    for j in range(0, 4):
                        m = s-1-j
                        if m >= 0:
                            ff = m*(m-1)*(m-2)
                            if ff: rhs += 2*comb(3, j)*r0**(3-j)*ff*a[m]
                    for j in range(0, 3):
                        m = s-2-j
                        if m >= 0:
                            ff = m*(m-1)
                            if ff: rhs += -al*comb(2, j)*r0**(2-j)*ff*a[m]
                    for j in range(0, 2):
                        m = s-3-j
                        if m >= 1:
                            rhs += al*comb(1, j)*r0**(1-j)*m*a[m]
                    if s-4 >= 0: rhs += be*a[s-4]
                    for j in range(0, 5):
                        m = s-4-j
                        if m >= 0:
                            rhs += -Ot2*comb(4, j)*r0**(4-j)*a[m]
                    a[s] = -rhs/(r0**4 * s*(s-1)*(s-2)*(s-3))
                sols.append(a)
            return sols
        # ---- Radially-graded path (Paper B): mirrors core_solvers.py's mp
        # _series_mp graded branch exactly (same corrected, LINEAR-in-T(x)/
        # R(x) beta(x) formula, same term-by-term widening) -- see that
        # method's own comment for the full derivation/cross-check. ----
        def _poly_mul(p, q):
            """Multiply two power-series coefficient lists (convolution)."""
            n = len(p) + len(q) - 1
            out = [0.0+0.0j] * n
            for ip, pv in enumerate(p):
                if pv == 0:
                    continue
                for iq, qv in enumerate(q):
                    out[ip+iq] += pv*qv
            return out
        def _binom_pow(n):
            """Binomial expansion of (r0+x)^n as a coefficient list in x."""
            return [comb(n, k)*r0**(n-k) + 0.0j for k in range(n+1)]
        Tl = [complex(c) for c in T]; Rl = [complex(c) for c in R]
        maxlen = max(len(Tl), len(Rl))
        Tl += [0.0+0.0j]*(maxlen-len(Tl)); Rl += [0.0+0.0j]*(maxlen-len(Rl))
        xi2 = xi**2
        Phi = [2*xi2*Tk + Rk for Tk, Rk in zip(Tl, Rl)]
        Psi = [xi2**2*Rk - 2*xi2*Tk - 2*xi2*Rk for Tk, Rk in zip(Tl, Rl)]
        rx2 = _binom_pow(2); rx1 = _binom_pow(1)
        rx2_Phi = _poly_mul(rx2, Phi)
        rx1_Phi = _poly_mul(rx1, Phi)
        # Omega^2-term ("K_omega") grading -- mirrors core_solvers.py's mp
        # _series_mp graded branch exactly (same _kconv convolution); see
        # that method's comment for the derivation pointer.
        Kw_raw = getattr(self, 'K_omega', None)
        if isinstance(Kw_raw, list):
            Kw = [complex(c) for c in Kw_raw]
        elif Kw_raw is not None:
            Kw = [complex(Kw_raw)]
        else:
            Kw = [1.0+0.0j]
        def _kconv(a, n):
            if n < 0:
                return 0.0+0.0j
            kmax = n if n < len(Kw) - 1 else len(Kw) - 1
            v = 0.0+0.0j
            for k in range(kmax + 1):
                v += Kw[k] * a[n - k]
            return v
        sols = []
        for q in range(4):
            a = np.zeros(M+5, complex); a[q] = 1.0
            for s in range(4, M+5):
                rhs = 0.0+0.0j
                for j in range(1, 5):
                    m = s-j
                    if m >= 0:
                        ff = m*(m-1)*(m-2)*(m-3)
                        if ff: rhs += comb(4, j)*r0**(4-j)*ff*a[m]
                for j in range(0, 4):
                    m = s-1-j
                    if m >= 0:
                        ff = m*(m-1)*(m-2)
                        if ff: rhs += 2*comb(3, j)*r0**(3-j)*ff*a[m]
                for k in range(len(rx2_Phi)):
                    m = s-2-k
                    if m >= 0:
                        ff = m*(m-1)
                        if ff: rhs += -rx2_Phi[k]*ff*a[m]
                for k in range(len(rx1_Phi)):
                    m = s-3-k
                    if m >= 1:
                        rhs += rx1_Phi[k]*m*a[m]
                for k in range(len(Psi)):
                    m = s-4-k
                    if m >= 0:
                        rhs += Psi[k]*a[m]
                for j in range(0, 5):
                    m = s-4-j
                    if m >= 0:
                        rhs += -Ot2*comb(4, j)*r0**(4-j)*_kconv(a, m)
                a[s] = -rhs/(r0**4 * s*(s-1)*(s-2)*(s-3))
            sols.append(a)
        return sols

    @staticmethod
    def _ev(a, x, d=0):
        """Evaluate d-th derivative of power series Σ a[m] x^m at x."""
        v = 0.0+0.0j
        for m in range(len(a)-1, d-1, -1):
            ff = 1.0
            for t in range(d): ff *= (m-t)
            if ff: v += a[m]*ff*x**(m-d)
        return v

    def Lmat(self, xi, sols):
        """4×4 free-arc BC matrix (moment/shear-type rows) at xi from series sols."""
        nu, T, R, r0 = self.nu, self.T, self.R, self.r0b
        hp = PI/2; xi2 = complex(xi)**2
        L = np.zeros((4, 4), complex)
        for q in range(4):
            a = sols[q]
            for row, x in enumerate((-hp, hp)):
                rr = x + r0
                nu_x = _coef_eval_f64(nu, x); T_x = _coef_eval_f64(T, x); R_x = _coef_eval_f64(R, x)
                G = self._ev(a, x, 0); Gp = self._ev(a, x, 1)
                Gpp = self._ev(a, x, 2); Gppp = self._ev(a, x, 3)
                L[row, q] = rr**2*Gpp + nu_x*(rr*Gp - xi2*G)
                L[row+2, q] = (rr*(rr**2*Gppp + rr*Gpp - R_x*Gp)
                               + xi2*((2*T_x+R_x-nu_x)*G - (2*T_x-nu_x)*rr*Gp))
        return L

    def det(self, xi, Om):
        """Row-scaled det(Lmat); max-abs scaling avoids float64 overflow on large Im(xi)."""
        L = self.Lmat(xi, self.series(xi, Om))
        # Scale each row by its maximum absolute entry before computing the
        # determinant.  np.linalg.norm squares the entries, which overflows
        # float64 when the Frobenius series has large imaginary-part roots
        # (e.g. xi ~ 1 + 9j).  Max-abs scaling avoids any squaring.
        for i in range(L.shape[0]):
            scale = np.max(np.abs(L[i]))
            if scale > 0:
                L[i] /= scale
        return np.linalg.det(L)

    def newton(self, z0, Om, itmax=40, tol=1e-13):
        """Complex Newton polish of det(z, Om)=0 starting at seed z0."""
        z = complex(z0)
        for _ in range(itmax):
            f = self.det(z, Om)
            h = 1e-7*max(1.0, abs(z))
            fp = (self.det(z+h, Om) - self.det(z-h, Om))/(2*h)
            if fp == 0: return None
            dz = f/fp
            z -= dz
            if abs(dz) < tol*max(1.0, abs(z)):
                return z
        return z if abs(self.det(z, Om)) < 1e-9 else None


class _Part2Fast(ExactEdgeSolver):
    """float64 coupled F/G power series + 4×4 arc-BC determinant for annular IP.

    Alias: AnnulusRadialIP. Part 2 (in-plane) exact-edge engine for full_search/
    track. Supports constant-coefficient scalars and (when radially_graded) Taylor
    lists for nu/c11/R — float64 mirror of InPlaneSolver._series_mp graded path
    (LESSONS §54.5). Keep the graded branch in sync with core_solvers.py.
    """
    geometry_kind = "annulus"
    motion = "in_plane"

    def __init__(self, r0b, nu, c11p, R=1.0, M=80, radially_graded=False):
        """r0b=r0_bar; nu,c11p,R scalar or Taylor lists; M=series length."""
        self.r0b = float(r0b)
        self.radially_graded = bool(radially_graded)
        # Graded materials pass nu/c11p/R as LISTS of float Taylor
        # coefficients (LESSONS_LEARNED Sec 54.4) -- keep them as-is;
        # non-graded materials pass plain scalars, force-cast exactly as
        # before (byte-for-byte unchanged). Mirrors _Part1Fast.__init__.
        if self.radially_graded:
            self.nu = [float(c) for c in nu] if isinstance(nu, (list, tuple)) else float(nu)
            self.c11 = [float(c) for c in c11p] if isinstance(c11p, (list, tuple)) else float(c11p)
            self.R = [float(c) for c in R] if isinstance(R, (list, tuple)) else float(R)
        else:
            self.nu = float(nu); self.c11 = float(c11p); self.R = float(R)
        self.M = M

    def series(self, ze, Om):
        """Corrected coupled recursion.  4 unit-IC solutions over
        (F(0), F'(0), G(0), G'(0))."""
        ze = complex(ze); O2 = complex(Om)**2
        r0 = self.r0b
        M = self.M
        if not isinstance(self.c11, list):
            # ---- EXACT existing constant-coefficient path, byte-for-byte
            # unchanged ----
            c11, R = self.c11, self.R
            cF = c11*R + ze**2
            cG = 1 + c11*R*ze**2
            coup = (c11*self.nu + 1)*ze
            cst = (c11*R + 1)*ze
            sols = []
            for q in range(4):
                a = np.zeros(M+3, complex); b = np.zeros(M+3, complex)
                if q == 0: a[0] = 1
                elif q == 1: a[1] = 1
                elif q == 2: b[0] = 1
                else: b[1] = 1
                for p in range(0, M+1):
                    am2 = a[p-2] if p >= 2 else 0.0
                    am1 = a[p-1] if p >= 1 else 0.0
                    bm2 = b[p-2] if p >= 2 else 0.0
                    bm1 = b[p-1] if p >= 1 else 0.0
                    rhsF = (c11*(2*r0*(p+1)*p*a[p+1] + p*(p-1)*a[p]
                                 + r0*(p+1)*a[p+1] + p*a[p])
                            + O2*(am2 + 2*r0*am1 + r0**2*a[p])
                            - cF*a[p]
                            - coup*(r0*(p+1)*b[p+1] + p*b[p])
                            + cst*b[p])
                    a[p+2] = -rhsF/(c11*r0**2*(p+2)*(p+1))
                    rhsG = ((2*r0*(p+1)*p*b[p+1] + p*(p-1)*b[p]
                             + r0*(p+1)*b[p+1] + p*b[p])
                            + O2*(bm2 + 2*r0*bm1 + r0**2*b[p])
                            - cG*b[p]
                            + coup*(r0*(p+1)*a[p+1] + p*a[p])
                            + cst*a[p])
                    b[p+2] = -rhsG/(r0**2*(p+2)*(p+1))
                sols.append((a, b))
            return sols
        # ---- Radially-graded path (Paper B / LESSONS_LEARNED Sec 54.4):
        # float64 mirror of core_solvers.py InPlaneSolver._series_mp's graded
        # branch, same linear-convolution formula -- keep both in sync. ----
        c11s = [complex(c) for c in self.c11]
        Rs   = [complex(c) for c in self.R]
        nus  = [complex(c) for c in self.nu] if isinstance(self.nu, list) else [complex(self.nu)]
        N = M + 2

        def _conv_full(p, q, n):
            """Full linear convolution of two Taylor lists up to degree n."""
            out = [0.0+0.0j] * (n + 1)
            for j in range(n + 1):
                s = 0.0+0.0j
                for k in range(0, j + 1):
                    pk = p[k] if k < len(p) else 0.0+0.0j
                    qk = q[j - k] if (j - k) < len(q) else 0.0+0.0j
                    s += pk * qk
                out[j] = s
            return out

        c11R = _conv_full(c11s, Rs, N)
        c11nu = _conv_full(c11s, nus, N)
        cF   = [c11R[j] + (ze**2 if j == 0 else 0.0+0.0j) for j in range(N + 1)]
        cG   = [(1.0+0.0j if j == 0 else 0.0+0.0j) + c11R[j]*ze**2 for j in range(N + 1)]
        coup = [(c11nu[j] + (1.0+0.0j if j == 0 else 0.0+0.0j))*ze for j in range(N + 1)]
        cst  = [(c11R[j]  + (1.0+0.0j if j == 0 else 0.0+0.0j))*ze for j in range(N + 1)]

        def _dot(cs, xs, pp):
            """Dot material Taylor list cs with kinematic series xs at order pp."""
            # cs may be a raw (short) material Taylor list, not zero-padded
            # to N+1 -- guard both operands defensively (xs is always
            # p+1 long by construction, but the bounds check costs nothing).
            s = 0.0+0.0j
            for k in range(pp + 1):
                ck = cs[k] if k < len(cs) else 0.0+0.0j
                xk = xs[pp - k] if (pp - k) < len(xs) else 0.0+0.0j
                s += ck*xk
            return s

        sols = []
        for q in range(4):
            a = np.zeros(M+3, complex); b = np.zeros(M+3, complex)
            if q == 0: a[0] = 1
            elif q == 1: a[1] = 1
            elif q == 2: b[0] = 1
            else: b[1] = 1
            for p in range(0, M+1):
                am2 = a[p-2] if p >= 2 else 0.0
                am1 = a[p-1] if p >= 1 else 0.0
                bm2 = b[p-2] if p >= 2 else 0.0
                bm1 = b[p-1] if p >= 1 else 0.0
                r2A_p = am2 + 2*r0*am1 + r0**2*a[p]
                r2B_p = bm2 + 2*r0*bm1 + r0**2*b[p]
                opA = [r0**2*(j+2)*(j+1)*a[j+2] + (2*r0*(j+1)*j + r0*(j+1))*a[j+1]
                       + (j*(j-1)+j)*a[j] for j in range(p+1)]
                opB = [r0**2*(j+2)*(j+1)*b[j+2] + (2*r0*(j+1)*j + r0*(j+1))*b[j+1]
                       + (j*(j-1)+j)*b[j] for j in range(p+1)]
                rBp = [r0*(j+1)*b[j+1] + j*b[j] for j in range(p+1)]
                rAp = [r0*(j+1)*a[j+1] + j*a[j] for j in range(p+1)]
                known_resF = (_dot(c11s, opA, p) + O2*r2A_p - _dot(cF, a[:p+1], p)
                              - _dot(coup, rBp, p) + _dot(cst, b[:p+1], p))
                a[p+2] = -known_resF/(c11s[0]*r0**2*(p+2)*(p+1))
                known_resG = (opB[p] + O2*r2B_p - _dot(cG, b[:p+1], p)
                              + _dot(coup, rAp, p) + _dot(cst, a[:p+1], p))
                b[p+2] = -known_resG/(r0**2*(p+2)*(p+1))
            sols.append((a, b))
        return sols

    @staticmethod
    def _ev(c, x, d=0):
        """Evaluate d-th derivative of power series Σ c[m] x^m at x (float64)."""
        v = 0.0+0.0j
        for m in range(len(c)-1, d-1, -1):
            ff = 1.0
            for t in range(d): ff *= (m-t)
            if ff: v += c[m]*ff*x**(m-d)
        return v

    def Lmat(self, ze, sols):
        """4×4 free-arc traction BC matrix (t_rr, t_rθ); nu may be graded via _coef_eval_f64."""
        nu, r0 = self.nu, self.r0b
        hp = PI/2
        L = np.zeros((4, 4), complex)
        for q, (fa, fb) in enumerate(sols):
            for row, x in enumerate((-hp, hp)):
                rr = x + r0
                F = self._ev(fa, x, 0); Fp = self._ev(fa, x, 1)
                G = self._ev(fb, x, 0); Gp = self._ev(fb, x, 1)
                nu_x = _coef_eval_f64(nu, x)
                L[row, q] = Fp + nu_x*(F - ze*G)/rr        # t_rr = 0
                L[row+2, q] = Gp + (ze*F - G)/rr          # t_rθ = 0
        return L

    def det(self, ze, Om):
        """Row-scaled det(Lmat); same max-abs scaling as Part 1."""
        L = self.Lmat(ze, self.series(ze, Om))
        # Same max-abs row scaling as Part 1 — avoids norm overflow.
        for i in range(L.shape[0]):
            scale = np.max(np.abs(L[i]))
            if scale > 0:
                L[i] /= scale
        return np.linalg.det(L)

    def newton(self, z0, Om, itmax=40, tol=1e-13):
        """Complex Newton polish of det(z, Om)=0 for the IP annular engine."""
        z = complex(z0)
        for _ in range(itmax):
            f = self.det(z, Om)
            h = 1e-7*max(1.0, abs(z))
            fp = (self.det(z+h, Om) - self.det(z-h, Om))/(2*h)
            if fp == 0: return None
            dz = f/fp
            z -= dz
            if abs(dz) < tol*max(1.0, abs(z)):
                return z
        return z if abs(self.det(z, Om)) < 1e-9 else None


# ── Branch-degeneracy guard (rejects the geometry-invariant 0.4542 artifact) ──
# A spurious singularity appears when two canonical branch roots ξ/ζ nearly
# coincide (near-parallel basis columns → K near-singular for EVERY 2Θ).  Reject
# a dip whose basis min pairwise root spacing < SIGMIN_DEGEN_TOL (real ≈ 0.8,
# artifact ≈ 0.05).  Cause-based, not a hard-coded drop.  See §Detector.
SIGMIN_DEGEN_SENTINEL = 99.0

AnnulusRadialOOP = _Part1Fast
AnnulusRadialIP = _Part2Fast


# ══════════════════════════════════════════════════════════════════════════════
#  SECTION 4 — PART 1: OUT-OF-PLANE SOLVER (mp precision K / det)
# ══════════════════════════════════════════════════════════════════════════════

# ==============================================================================
#  BOUNDARY CONDITIONS  (Phase-0 Step 3)
#  The annular sector has two straight (radial) theta-edges at theta=+/-Theta.
#  Their conditions define the problem; the arc (Frobenius/exact-edge) machinery
#  is BC-independent.  This axis is ORTHOGONAL to the material's radially_graded
#  flag (Paper B), which lives on the radial ODE coefficients: a graded free-free
#  plate is simply boundary=FreeFree* together with a radially_graded material,
#  set independently.
# ==============================================================================

def _rect_dom(d):
    """Signed dominant component of a complex value (phase-agnostic real face of
    the dispersion determinant): the antisymmetric class is purely imaginary on
    the axes, so full_search's det(...).real scan needs this to see the roots."""
    return d.real if abs(d.real) >= abs(d.imag) else d.imag


class RectangularCartesianOOP(ExactEdgeSolver):
    """Out-of-plane exact-edge engine for the rectangular cantilever plate
    (Seok, Tiersten & Scarton 2004, Part 1), fitted to the ExactEdgeSolver
    contract.

    NOTE ON THE FREQUENCY ARGUMENT: this strategy is parameterised by
    Lambda = k_bar * Omega_bar (the Table-3 convention).  The `Om` argument of
    det/newton/series IS Lambda here (the rectangular exact-edge problem depends
    on material only through (T,R,nu_bar,mu_bar), not on b/h -- working in Lambda
    is what makes the dispersion "a function of Poisson's ratio only").

    NOTE ON SYMMETRY: the symmetric (sym=True) and antisymmetric (sym=False)
    waves are separate eigenproblems; `sym` is engine state so det keeps the
    annular (wavenumber, Om) signature.  Use two instances for the two classes.

    `det(xi, Om)` returns the DOMINANT-COMPONENT real-valued dispersion function
    (so the shared full_search .real sign-change scan finds both classes), with
    eta-branch-point candidates (min|eta|~0 or eta1^2~eta2^2) pushed away from
    zero so they are rejected by full_search's |det|>1e-6 filter rather than
    mistaken for plate modes.
    """
    geometry_kind = "rectangle"
    motion = "out_of_plane"

    def __init__(self, nu_b=0.3, R=1.0, T=1.0, sym=True, eta_tol=0.05):
        """nu_b, R, T: material ratios; sym selects even/odd class; eta_tol for filters."""
        self.nu = float(nu_b)
        self.R = float(R)
        self.T = float(T)
        self.mu = 2.0 * self.T - self.nu          # Eq. 28
        self.sym = bool(sym)
        self.eta_tol = float(eta_tol)
        self.radially_graded = False

    # exact-edge transverse wavenumbers (Eq. 34, dimensionless) ----------------
    def _eta_sq(self, xi, Lam):
        """Return (eta1^2, eta2^2) for longitudinal wavenumber xi and Lambda."""
        T, R = self.T, self.R
        disc = _cmath.sqrt((T * T - R) * xi**4 + R * R * Lam * Lam)
        return (-T * xi * xi - disc) / R, (-T * xi * xi + disc) / R

    def _etas(self, xi, Lam):
        """Return (eta1, eta2) square-roots of _eta_sq."""
        e1, e2 = self._eta_sq(xi, Lam)
        return _cmath.sqrt(e1), _cmath.sqrt(e2)

    def _amp_ratio(self, xi, Lam):
        """Amplitude ratios (Eq. 40), for this engine's symmetry class.  Used by
        the Stage-2 variational assembler (run_rect_validation)."""
        R, nu = self.R, self.nu
        e1, e2 = self._etas(complex(xi), complex(Lam))
        x2 = complex(xi) * complex(xi)
        ph = (0.0 if not self.sym else _cmath.pi / 2.0)
        C1 = -(R * e2 * e2 + nu * x2) * _cmath.sin(_cmath.pi * e2 / 2.0 + ph)
        C2 = (R * e1 * e1 + nu * x2) * _cmath.sin(_cmath.pi * e1 / 2.0 + ph)
        return C1, C2

    # ExactEdgeSolver contract -------------------------------------------------
    def series(self, wavenumber, Om):
        """Closed-form rectangular exact-edge solution: there is no Frobenius
        series (constant-coefficient PDE).  Returns the transverse wavenumbers as
        the 'basis descriptor' consumed by Lmat."""
        return self._etas(complex(wavenumber), complex(Om))

    def Lmat(self, wavenumber, sols):
        """2x2 free-edge matrix (Eq. 36); `sols` are (eta1, eta2) from series."""
        R, nu, mu = self.R, self.nu, self.mu
        xi = complex(wavenumber)
        e1, e2 = sols
        x2 = xi * xi
        phi = 0.0 if not self.sym else _cmath.pi / 2.0   # n=1 antisym, n=2 sym
        a1 = _cmath.pi * e1 / 2.0 + phi
        a2 = _cmath.pi * e2 / 2.0 + phi
        P1 = R * e1 * e1 + nu * x2
        P2 = R * e2 * e2 + nu * x2
        Q1 = e1 * (R * e1 * e1 + mu * x2)
        Q2 = e2 * (R * e2 * e2 + mu * x2)
        return [[P1 * _cmath.sin(a1), P2 * _cmath.sin(a2)],
                [Q1 * _cmath.cos(a1), Q2 * _cmath.cos(a2)]]

    def _dispersion(self, xi, Lam):
        """Complex free-edge 2×2 determinant (Eq. 36) before dominant-component map."""
        e1, e2 = self._etas(complex(xi), complex(Lam))
        a1 = _cmath.pi * e1 / 2.0 + (0.0 if not self.sym else _cmath.pi / 2.0)
        a2 = _cmath.pi * e2 / 2.0 + (0.0 if not self.sym else _cmath.pi / 2.0)
        if abs(a1.imag) > 300 or abs(a2.imag) > 300:
            return complex(1e300, 0.0)
        R, nu, mu = self.R, self.nu, self.mu
        x2 = complex(xi) * complex(xi)
        P1 = R * e1 * e1 + nu * x2
        P2 = R * e2 * e2 + nu * x2
        Q1 = e1 * (R * e1 * e1 + mu * x2)
        Q2 = e2 * (R * e2 * e2 + mu * x2)
        try:
            return (P1 * _cmath.sin(a1) * Q2 * _cmath.cos(a2)
                    - P2 * _cmath.sin(a2) * Q1 * _cmath.cos(a1))
        except (OverflowError, ValueError):
            return complex(1e300, 0.0)

    def det(self, wavenumber, Om):
        """Phase-agnostic real-valued dispersion function (dominant component),
        with eta-branch-point candidates lifted away from zero so full_search
        rejects them.  `Om` is Lambda (see class docstring)."""
        xi = complex(wavenumber); Lam = complex(Om)
        val = _rect_dom(self._dispersion(xi, Lam))
        # eta branch points (xi=sqrt(Lam) etc.) are representation singularities,
        # not plate modes: lift |det| so full_search's |det|>1e-6 filter drops it.
        e1, e2 = self._etas(xi, Lam)
        if min(abs(e1), abs(e2)) < self.eta_tol or abs(e1*e1 - e2*e2) < self.eta_tol:
            return complex(1.0, 0.0)
        return complex(val, 0.0)

    def newton(self, z0, Om, itmax=40, tol=1e-13):
        """Newton polish of the rectangular OOP dispersion det (Lambda as Om)."""
        z = complex(z0)
        for _ in range(itmax):
            f = self.det(z, Om).real
            h = 1e-7 * max(1.0, abs(z))
            fp = (self.det(z + h, Om).real - self.det(z - h, Om).real) / (2 * h)
            if fp == 0:
                return None
            dz = f / fp
            z -= dz
            if abs(dz) < tol * max(1.0, abs(z)):
                return z
        return z if abs(self.det(z, Om).real) < 1e-9 else None


# -- Rectangular Cartesian IN-PLANE exact-edge strategy (Paper A, Part 2) -------
# Seok, Tiersten & Scarton, JSV 271 (2004) 147-158 (in-plane motion).  Gated by
# R40_RECT_IP; the annular path and the rect OOP path are unaffected when off.
# Dimensionless per Eqs. 22-23: gamma_bar=(2b/pi)*gamma (the x1 wavenumber, the
# `wavenumber`/`xi` argument here), zeta_bar=(2b/pi)*zeta (transverse), and the
# frequency argument `Om` IS Omega_bar = omega/omega_bar DIRECTLY (no k_bar
# factor, unlike the OOP Lambda convention).  Two symmetry classes (sym=engine
# state, like RectangularCartesianOOP); phase s' = (s-1)pi/2, s=2 sym -> pi/2,
# s=1 antisym -> 0.  The free-edge determinant Eq.27 is purely imaginary on the
# real-gamma axis for sym and purely real for antisym, so `det` returns the
# dominant component (_rect_dom) for the shared full_search .real scan -- exactly
# as in the OOP engine.


class RectangularCartesianIP(ExactEdgeSolver):
    """In-plane exact-edge engine for the rectangular cantilever plate
    (Seok, Tiersten & Scarton 2004, Part 2), fitted to the ExactEdgeSolver
    contract.

    Two-component displacement (u1,u2): for a given x1-wavenumber gamma_bar and
    frequency Omega_bar the transverse wavenumber zeta_bar^2 solves a biquadratic
    (Eq.24 determinant); each of the two roots zeta_(1),(2) carries amplitude
    ratios A1bar,A2bar (Eq.25-26).  The free-edge (x2=+/-b) 2x2 determinant
    (Eq.27) is the dispersion relation.  Material enters only through
    c11s = c*11/c66, nu_hat = c12/c11, R = c22/c11 (isotropic nu: c11s=2/(1-nu),
    nu_hat=nu, R=1).
    """
    geometry_kind = "rectangle"
    motion = "in_plane"

    def __init__(self, nu_b=0.3, R=1.0, c11s=None, sym=True, zeta_tol=0.05):
        """nu_b=nu_hat; R=c22/c11; c11s defaults to isotropic 2/(1-nu); sym class."""
        self.nu = float(nu_b)                       # nu_hat = c12/c11
        self.R = float(R)                           # c22/c11
        self.c11s = (2.0 / (1.0 - self.nu)) if c11s is None else float(c11s)
        self.sym = bool(sym)
        self.zeta_tol = float(zeta_tol)
        self.radially_graded = False

    @property
    def _phase(self):
        """Symmetry phase s'=(s-1)pi/2: pi/2 symmetric, 0 antisymmetric."""
        return (_cmath.pi / 2.0) if self.sym else 0.0   # s'=(s-1)pi/2

    # transverse wavenumbers: biquadratic in zeta_bar^2 (Eq.24 det=0) -----------
    def _zeta_sq(self, g, Om):
        """Return (zeta1^2, zeta2^2) for x1-wavenumber g and Omega_bar."""
        g = complex(g); Om = complex(Om)
        c11s, nu, R = self.c11s, self.nu, self.R
        g2 = g * g; Om2 = Om * Om
        a = c11s * R
        b = (c11s * g2 - Om2) * c11s * R + (g2 - Om2) - (c11s * nu + 1.0) ** 2 * g2
        c = (c11s * g2 - Om2) * (g2 - Om2)
        disc = _cmath.sqrt(b * b - 4 * a * c)
        return (-b - disc) / (2 * a), (-b + disc) / (2 * a)

    def _zetas(self, g, Om):
        """Return (zeta1, zeta2) square-roots of _zeta_sq."""
        z1, z2 = self._zeta_sq(g, Om)
        return _cmath.sqrt(z1), _cmath.sqrt(z2)

    def _amp(self, g, Om):
        """Per-branch amplitude ratios A1bar^(n), A2bar^(n) (Eq.25-26)."""
        g = complex(g); Om = complex(Om)
        c11s, nu = self.c11s, self.nu
        z1, z2 = self._zetas(g, Om)
        z1s, z2s = self._zeta_sq(g, Om)
        A1 = ((-(c11s * nu + 1.0) * g * z1), (-(c11s * nu + 1.0) * g * z2))
        A2 = ((c11s * g * g + z1s - Om * Om), (c11s * g * g + z2s - Om * Om))
        return (z1, z2), A1, A2

    def _Cratio(self, g, Om):
        """Inter-branch amplitude ratios C^(1),C^(2) (Eq.28-29) from the
        free-edge consistency: C^(1) = -(nu*g*A1^(2)+R*z2*A2^(2)) sin(pi z2/2+s'),
        C^(2) = (nu*g*A1^(1)+R*z1*A2^(1)) sin(pi z1/2+s')."""
        nu, R = self.nu, self.R
        (z1, z2), A1, A2 = self._amp(g, Om)
        ph = self._phase
        s1 = _cmath.sin(_cmath.pi * z1 / 2 + ph)
        s2 = _cmath.sin(_cmath.pi * z2 / 2 + ph)
        C1 = -(nu * complex(g) * A1[1] + R * z2 * A2[1]) * s2
        C2 = (nu * complex(g) * A1[0] + R * z1 * A2[0]) * s1
        return C1, C2

    def H(self, g, Om):
        """Transverse-profile coefficients H1[q], H2[q] = C^(q) * An^(q)
        (Eq.39), the per-branch packed amplitudes for the Stage-2 assembler."""
        (z1, z2), A1, A2 = self._amp(g, Om)
        C1, C2 = self._Cratio(g, Om)
        H1 = (C1 * A1[0], C2 * A1[1])
        H2 = (C1 * A2[0], C2 * A2[1])
        return (z1, z2), H1, H2

    # ExactEdgeSolver contract -------------------------------------------------
    def series(self, wavenumber, Om):
        """Exact-edge basis descriptor: (zeta1, zeta2) for this (gamma, Omega)."""
        return self._zetas(complex(wavenumber), complex(Om))

    def _dispersion(self, g, Om):
        """Free-edge 2x2 determinant Eq.27 (full complex)."""
        g = complex(g); Om = complex(Om)
        nu, R = self.nu, self.R
        (z1, z2), A1, A2 = self._amp(g, Om)
        ph = self._phase
        a1 = _cmath.pi * z1 / 2 + ph; a2 = _cmath.pi * z2 / 2 + ph
        if abs(a1.imag) > 300 or abs(a2.imag) > 300:
            return complex(1e300, 0.0)
        try:
            M11 = (nu * g * A1[0] + R * z1 * A2[0]) * _cmath.sin(a1)
            M12 = (nu * g * A1[1] + R * z2 * A2[1]) * _cmath.sin(a2)
            M21 = (z1 * A1[0] + g * A2[0]) * _cmath.cos(a1)
            M22 = (z2 * A1[1] + g * A2[1]) * _cmath.cos(a2)
            return M11 * M22 - M12 * M21
        except (OverflowError, ValueError):
            return complex(1e300, 0.0)

    def is_branch_point(self, g, Om, rel_tol=5e-3, abs_min=1e-3):
        """Scale-invariant zeta branch-point test.  The two in-plane transverse
        wavenumbers are physically near-degenerate for isotropic material, with
        |z1^2 - z2^2| ~ const while |z^2| grows with Omega; an ABSOLUTE
        |z1^2-z2^2| threshold therefore wrongly rejects genuine evanescent
        branches at low Omega / high l/b (verified: imag-gamma branch at
        Om=0.0749, l/b=2.5 was lost).  Use the RELATIVE separation instead, plus
        a tiny absolute floor on min|zeta| to catch the true confluent point."""
        z1, z2 = self._zetas(complex(g), complex(Om))
        z1s = z1 * z1; z2s = z2 * z2
        denom = abs(z1s) + abs(z2s)
        rel = abs(z1s - z2s) / denom if denom > 0 else 0.0
        return (min(abs(z1), abs(z2)) < abs_min) or (rel < rel_tol)

    def det(self, wavenumber, Om):
        """Phase-agnostic real dispersion (dominant component) with scale-
        invariant zeta branch-point rejection.  `Om` is Omega_bar."""
        g = complex(wavenumber); Om = complex(Om)
        val = _rect_dom(self._dispersion(g, Om))
        if self.is_branch_point(g, Om):
            return complex(1.0, 0.0)
        return complex(val, 0.0)

    def newton(self, z0, Om, itmax=40, tol=1e-13):
        """Newton polish on the real dominant-component dispersion (Eq. 27 face)."""
        z = complex(z0)
        for _ in range(itmax):
            f = self.det(z, Om).real
            h = 1e-7 * max(1.0, abs(z))
            fp = (self.det(z + h, Om).real - self.det(z - h, Om).real) / (2 * h)
            if fp == 0:
                return None
            dz = f / fp
            z -= dz
            if abs(dz) < tol * max(1.0, abs(z)):
                return z
        return z if abs(self.det(z, Om).real) < 1e-9 else None

    def cutoffs(self, n=4):
        """Closed-form cut-offs at gamma=0 (Eqs.34-35)."""
        import math as _math
        s = _math.sqrt(self.c11s * self.R)
        if not self.sym:                       # antisymmetric
            a = [2 * k - 1 for k in range(1, n + 1)]
            b = [s * (2 * k) for k in range(1, n + 1)]
        else:                                  # symmetric
            a = [2 * k for k in range(1, n + 1)]
            b = [s * (2 * k - 1) for k in range(1, n + 1)]
        return sorted(a + b)



# R40_RECT gate: descriptive registry entry only (Stage 1 = engine strategy).
# The rectangular variational K-assembly (Stage 2) is intentionally NOT wired in
# yet, so the validated annular assembler stays frozen.
_RECT_ENABLED = os.environ.get("R40_RECT", "0") == "1"


# ══════════════════════════════════════════════════════════════════════════════
#  SECTION 4 — PART 1: OUT-OF-PLANE SOLVER (mp precision K / det)
# ══════════════════════════════════════════════════════════════════════════════

# ==============================================================================
#  BOUNDARY CONDITIONS  (Phase-0 Step 3)
#  The annular sector has two straight (radial) theta-edges at theta=+/-Theta.
#  Their conditions define the problem; the arc (Frobenius/exact-edge) machinery
#  is BC-independent.  This axis is ORTHOGONAL to the material's radially_graded
#  flag (Paper B), which lives on the radial ODE coefficients: a graded free-free
#  plate is simply boundary=FreeFree* together with a radially_graded material,
#  set independently.
# ==============================================================================

def _rect_G(g, d):
    """Green-type 1D kernel integral used by rectangular variational assembly."""
    if abs(g) < mpf(10) ** -18:
        return mp.pi * mp.cos(d)
    return 2 * mp.cos(d) * mp.sin(g * mp.pi / 2) / g


def _rect_J(a, pa, b, pb):
    """Combination of _rect_G terms for two-wavenumber rectangular products."""
    return (_rect_G(a - b, pa - pb) - _rect_G(a + b, pa + pb)) / 2


