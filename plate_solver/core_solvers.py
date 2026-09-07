# -*- coding: utf-8 -*-
"""
plate_solver.core_solvers -- the big variational-assembly solver classes:
OutOfPlaneSolver / InPlaneSolver (annular, Parts 1 & 2) and
RectOOPAssembler / RectIPAssembler (rectangular, ported from rect_int.py).

BUG FIX (review-doc item 2): InPlaneSolver.find_natural_frequencies used to
initialize `freqs = []` twice in a row (harmless but confusing dead code);
the duplicate line has been removed.

MECHANICAL CLEANUP (review-doc item, mode_shape_grid): the SVD fallback used
when the arbitrary-precision null-vector extraction is disabled/fails now
has its OWN narrow try/except around just the np.linalg.svd call (instead
of relying on the method's single outer blanket except), so a
non-convergent fallback degrades to `c = None` (and is reported once)
without being silently conflated with unrelated earlier failures in the
same method.  Both changes are plotting/reporting-path only; no computed
frequency is affected, SOLVER_VERSION not bumped.

Extracted verbatim otherwise (line-range provenance in LESSONS_LEARNED
Sec. 22).
"""
from __future__ import annotations
import os, sys, time, cmath as _cmath
import numpy as np
from math import comb
from mpmath import mp, mpf, mpc, matrix
import concurrent.futures

from .config import MP_DPS, N_WORKERS, PREFETCH, PI, MAX_SCAN_PTS, \
    _P1_LD_GENUINE, _P1_PROM, _P2_LD_GENUINE, SOLVER_VERSION
from .geometry import PlateGeometry, MaterialModel, IsotropicMaterial, \
    OrthotropicMaterial, FGMPlateProperties, RadialFGMMaterial, \
    _ScalarGeom, _ScalarMat, _make_geom_mat
from .boundary import EdgeKind, EdgeSpec, BoundaryCondition, ClampedFreeOOP, \
    ClampedFreeIP, FreeFreeOOP, FreeFreeIP, make_bc, _worker_bc, \
    _OOP_CORNER_COEFF, _BC_REGISTRY, _BC_REGISTRY_TOKENS
from .dispersion import ExactEdgeSolver, _Part1Fast, _Part2Fast, \
    RectangularCartesianOOP, RectangularCartesianIP, _rect_dom, _rect_G, _rect_J
from .detectors import canon, real_basis_items, dof_count, select_fill, \
    set_signature, mp_gauss_legendre, mp_proj, equilibrated_logdet, \
    sigma_min_from_K, equilibrated_nullvec_mp, _nullvec_mp_rect, \
    _theta_collocation_coeffs, _quick_shape_stats, ritz_spectrum_oop, \
    classify_modes_ritz, _min_root_spacing, _degen_guard, full_search, track, \
    find_modes_sigmin, SIGMIN_DEGEN_SENTINEL, _rect_absdet_axis_roots, \
    _rect_cnewton, rect_resolve_branches, rect_select_branches, \
    _rect_ip_cnewton, _rect_ip_axis_roots, rect_ip_resolve_branches, \
    rect_ip_select_branches
from .workers import _p1_logdet_worker, _p1_bisect_worker, _p1_golden_worker, \
    _p1_ld_at_worker, _p1_sigmin_worker, _p1_sigmin_golden_worker, \
    _p2_logdet_worker, _p2_bisect_worker, _p2_golden_worker, _p2_ld_at_worker, \
    _p2_sigmin_worker, _p2_sigmin_golden_worker, _surface_worker_exc, \
    _throttled_map

AnnulusRadialOOP = _Part1Fast
AnnulusRadialIP = _Part2Fast


def _coef_eval(coef, x):
    """Evaluate a scalar-or-series ODE coefficient at x (Paper B / FGM).

    Non-graded materials pass a plain mpf/mpc scalar (self.T/self.R/self.nu
    are never lists for IsotropicMaterial/OrthotropicMaterial) -- this
    returns it UNCHANGED, so every existing call site is a byte-for-byte
    no-op. Graded materials (RadialFGMMaterial) pass a list of Taylor
    coefficients (index = power of x, rr=x+r0); this Horner-evaluates it
    pointwise. Used ONLY in the arc/theta-edge algebra (_Lmat_mp,
    _build_K_real's Tyy/Qy/Mc), which Seok & Tiersten Eq.(1a)/(6) shows are
    evaluated locally at a fixed r (the theta-edge terms differentiate only
    w.r.t. theta, and the corner/arc terms are point evaluations at r=Ri/Ro
    -- an r-derivative of T/R/nu never appears in either) -- so pointwise
    substitution is exact, no extra product-rule terms needed. The radial
    ODE itself (_series_mp) DOES differentiate w.r.t. r and needs the full
    widened convolution recursion instead; see that method's own graded
    branch. Full derivation: FUTURE_WORK.md's FGM item."""
    if isinstance(coef, list):
        v = mpf(0)
        for c in reversed(coef):
            v = v * x + c
        return v
    return coef


class OutOfPlaneSolver:
    """Flexural vibration, Seok & Tiersten Part 1.  Frequencies are Ω̃=k̄Ω̄."""
    def __init__(self, geom, mat, M=80, n_quad=30, boundary=None):
        self.geom = geom; self.mat = mat; self.M = M
        self.r0b = geom.r0_bar
        # theta-edge boundary condition (Phase-0 Step 3).  Defaults to the
        # clamped-free cantilever so the worker reconstruction path (which
        # builds the solver with no `boundary` arg) is unchanged -- same
        # safe-default pattern as radially_graded=False.
        self.bc = boundary if boundary is not None else ClampedFreeOOP()
        # Pull the OOP ODE constants through the MaterialModel interface (Phase-0
        # refactor).  For IsotropicMaterial these are exactly mat.T, mat.R,
        # mat.nu_bar as before, so frequencies are unchanged; a graded material
        # would supply graded coefficients here without this class changing.
        if mat.radially_graded:
            self.T, self.R, self.nu = mat.oop_coeff_series(self.r0b, M)
            # K_omega: OPTIONAL Omega^2-term grading, added 2026-08-02 --
            # see geometry.py::RadialFGMMaterial.oop_kappa_series's
            # docstring and GrokCode/KGradingDerivation.txt. getattr
            # fallback keeps this a no-op for any radially_graded
            # material that predates this addition (or doesn't supply
            # the method), since [mpf(1)] makes the _series_mp
            # convolution below an exact identity.
            _kappa_fn = getattr(mat, 'oop_kappa_series', None)
            self.K_omega = (list(_kappa_fn(self.r0b, M)) if _kappa_fn is not None
                            else [mpf(1)])
        else:
            self.T, self.R, self.nu = mat.oop_constants()
            self.K_omega = [mpf(1)]
        self.Theta = geom.Theta
        self.k_bar = float((mpf(2)*geom.b/mp.pi) *
                           mp.sqrt(3*mat.c66/(geom.h**2*mat.c11_bar)))
        # float() cast only for the non-graded (scalar) case -- forcing a
        # graded material's T/R/nu (Taylor-coefficient LISTS) through
        # float() would raise; pass the list through unchanged instead, and
        # AnnulusRadialOOP's own radially_graded branch handles it.
        _fast_nu = [float(c) for c in self.nu] if mat.radially_graded else float(self.nu)
        _fast_T = [float(c) for c in self.T] if mat.radially_graded else float(self.T)
        _fast_R = [float(c) for c in self.R] if mat.radially_graded else float(self.R)
        _fast_Kw = [float(c) for c in self.K_omega] if mat.radially_graded else 1.0
        self.fast = AnnulusRadialOOP(float(self.r0b), _fast_nu,
                                     _fast_T, _fast_R, M,
                                     radially_graded=mat.radially_graded,
                                     K_omega=_fast_Kw)
        nd, wt = mp_gauss_legendre(n_quad)
        self.nodes = [x*mp.pi/2 for x in nd]
        self.wts = [w*mp.pi/2 for w in wt]

    # ---------- mp series / arc-BC machinery ----------
    def _series_mp(self, xi, Om):
        xi = mpc(xi); Ot2 = mpc(Om)**2
        T, R, r0 = self.T, self.R, self.r0b
        if not isinstance(T, list):
            # ---- EXACT existing constant-coefficient path, byte-for-byte
            # unchanged (isotropic/orthotropic materials never take the
            # graded branch below -- self.T/self.R are plain mpf scalars). ----
            al = 2*T*xi**2 + R
            # β: polar-orthotropic form (see float64 series); R=1 ⇒ byte-identical.
            be = xi**2*(R*xi**2 - 2*(T+R))
            M = self.M
            sols = []
            for q in range(4):
                a = [mpc(0)]*(M+5); a[q] = mpc(1)
                for s in range(4, M+5):
                    rhs = mpc(0)
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
        # ---- Radially-graded path (Paper B): T, R are Taylor-coefficient
        # lists in x (rr=x+r0). alpha(x)=2*T(x)*xi^2+R(x), beta(x)=
        # R(x)*xi^4-2*(T(x)+R(x))*xi^2 -- coefficient-wise generalization of
        # LESSONS Sec 8's corrected polar-orthotropic form (linear in
        # T(x)/R(x), NOT quadratic -- probe_fgm_graded_recursion_selftest.py's
        # own beta(x)=R(x)*xi^2*(xi^2-2T(x)-2R(x)) formula was checked and
        # found to only agree with this at R=1; it reintroduces exactly the
        # "spurious extra R" bug LESSONS Sec 8 already fixed once, so it is
        # NOT used here -- see FUTURE_WORK.md's FGM item for the full
        # symbolic cross-check). Every OTHER block below (the r0-binomial
        # (x+r0)^n expansions with no T/R at all) is untouched structurally
        # from the constant-coefficient branch, just widened from a fixed
        # 2-4 term window to a sum over however many terms alpha(x)'s own
        # series has -- verified to collapse EXACTLY to the branch above
        # when T,R are length-1 lists (single-term series).
        M = self.M
        def _poly_mul(p, q):
            n = len(p) + len(q) - 1
            out = [mpc(0)] * n
            for ip, pv in enumerate(p):
                if pv == 0:
                    continue
                for iq, qv in enumerate(q):
                    out[ip+iq] += pv*qv
            return out
        def _binom_pow(n):
            return [mpc(comb(n, k)) * r0**(n-k) for k in range(n+1)]
        Tl = [mpc(c) for c in T]; Rl = [mpc(c) for c in R]
        maxlen = max(len(Tl), len(Rl))
        Tl += [mpc(0)]*(maxlen-len(Tl)); Rl += [mpc(0)]*(maxlen-len(Rl))
        xi2 = xi**2
        Phi = [2*xi2*Tk + Rk for Tk, Rk in zip(Tl, Rl)]
        Psi = [xi2**2*Rk - 2*xi2*Tk - 2*xi2*Rk for Tk, Rk in zip(Tl, Rl)]
        rx2 = _binom_pow(2); rx1 = _binom_pow(1)
        rx2_Phi = _poly_mul(rx2, Phi)
        rx1_Phi = _poly_mul(rx1, Phi)
        # Omega^2-term ("K_omega") grading -- added 2026-08-02, see
        # geometry.py::RadialFGMMaterial.oop_kappa_series's docstring and
        # GrokCode/KGradingDerivation.txt for the derivation. Kw defaults
        # to [1] (getattr fallback) for any graded material that doesn't
        # supply K_omega, making _kconv(a, n) == a[n] exactly -- a byte-
        # for-byte no-op for every T/R-only-graded call site that
        # predates this addition.
        Kw_raw = getattr(self, 'K_omega', None)
        Kw = [mpc(c) for c in Kw_raw] if Kw_raw is not None else [mpc(1)]
        def _kconv(a, n):
            if n < 0:
                return mpc(0)
            kmax = n if n < len(Kw) - 1 else len(Kw) - 1
            v = mpc(0)
            for k in range(kmax + 1):
                v += Kw[k] * a[n - k]
            return v
        sols = []
        for q in range(4):
            a = [mpc(0)]*(M+5); a[q] = mpc(1)
            for s in range(4, M+5):
                rhs = mpc(0)
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
    def _ev_mp(a, x, d=0):
        """Evaluate d-th derivative of Σ a[m]·x^m using Horner's rule.
        Horner gives ~12× speedup over the x**(m-d) loop for d=0,1,2.
        """
        N = len(a) - 1
        if d == 0:
            v = a[N]
            for m in range(N-1, -1, -1):
                v = v*x + a[m]
            return v
        elif d == 1:
            v = mpc(0)
            for m in range(N, 0, -1):
                v = v*x + m*a[m]
            return v
        elif d == 2:
            v = mpc(0)
            for m in range(N, 1, -1):
                v = v*x + m*(m-1)*a[m]
            return v
        else:  # d==3 (only used in Lmat)
            v = mpc(0)
            for m in range(N, 2, -1):
                v = v*x + m*(m-1)*(m-2)*a[m]
            return v

    def _Lmat_mp(self, xi, sols):
        nu, T, R, r0 = self.nu, self.T, self.R, self.r0b
        hp = mp.pi/2; xi2 = mpc(xi)**2
        L = matrix(4, 4)
        for q in range(4):
            a = sols[q]
            for row, x in enumerate((-hp, hp)):
                rr = x + r0
                # _coef_eval is a no-op (returns nu/T/R unchanged) unless the
                # material is radially_graded (Paper B), in which case it
                # Horner-evaluates the Taylor series at this arc position x.
                nu_x = _coef_eval(nu, x); T_x = _coef_eval(T, x); R_x = _coef_eval(R, x)
                G = self._ev_mp(a, x, 0); Gp = self._ev_mp(a, x, 1)
                Gpp = self._ev_mp(a, x, 2); Gppp = self._ev_mp(a, x, 3)
                L[row, q] = rr**2*Gpp + nu_x*(rr*Gp - xi2*G)
                L[row+2, q] = (rr*(rr**2*Gppp + rr*Gpp - R_x*Gp)
                               + xi2*((2*T_x+R_x-nu_x)*G - (2*T_x-nu_x)*rr*Gp))
        return L

    def _det_mp(self, xi, Om):
        L = self._Lmat_mp(xi, self._series_mp(xi, Om))
        for i in range(4):
            s = max(abs(L[i, j]) for j in range(4))
            if s > 0:
                for j in range(4): L[i, j] /= s
        return mp.det(L)

    def _refine_root_mp(self, z, Om, iters=4):
        """A few mp Newton steps to reach mp-level null residual."""
        z = mpc(z)
        for _ in range(iters):
            f = self._det_mp(z, Om)
            h = mpf(10)**(-(mp.dps//2))*max(mpf(1), abs(z))
            fp = (self._det_mp(z+h, Om) - self._det_mp(z-h, Om))/(2*h)
            if fp == 0: break
            dz = f/fp
            z -= dz
            if abs(dz) < mpf(10)**(-(mp.dps-8)): break
        return z

    def _amp_mp(self, xi, Om, sols):
        """Null vector of L (pivot column + 3 best rows), normalised W(0)=1."""
        L = self._Lmat_mp(xi, sols)
        cn = [sum(abs(L[r, q])**2 for r in range(4)) for q in range(4)]
        pc = max(range(4), key=lambda q: cn[q])
        oc = [q for q in range(4) if q != pc]
        rn = [sum(abs(L[r, q])**2 for q in range(4)) for r in range(4)]
        ur = sorted(range(4), key=lambda r: rn[r], reverse=True)[:3]
        sub = matrix(3, 3); rhs = matrix(3, 1)
        for i, row in enumerate(ur):
            rhs[i, 0] = -L[row, pc]
            for j, col in enumerate(oc):
                sub[i, j] = L[row, col]
        try:
            t = mp.lu_solve(sub, rhs)
            A = [mpc(0)]*4; A[pc] = mpc(1)
            for j, col in enumerate(oc): A[col] = t[j, 0]
        except Exception:
            A = [mpc(1), mpc(0), mpc(0), mpc(0)]
        W0 = sum(A[q]*self._ev_mp(sols[q], mpf(0), 0) for q in range(4))
        if abs(W0) > mpf('1e-30'):
            A = [A[q]/W0 for q in range(4)]
        return A

    def _W_mp(self, x, sols, A, d=0):
        return sum(A[q]*self._ev_mp(sols[q], x, d) for q in range(4))

    # ---------- mp K matrix in the REAL basis ----------
    def _build_K_real(self, Om_f, brs, fast_scan=False):
        """Real N_dof × N_dof variational matrix from Eq.(43).

        fast_scan=True skips the mp Newton root-refinement (_refine_root_mp)
        and uses the float64 branch roots directly.  This gives ~8× speedup
        for scan-phase evaluations where sign detection — not high accuracy —
        is the goal.  Always use fast_scan=False for final polish evaluations.
        """
        Om = mpf(Om_f)
        nu, T, R, r0 = self.nu, self.T, self.R, self.r0b
        Th = self.Theta
        items = real_basis_items(brs)
        n = len(items)

        # cache series/amp per distinct root
        cache = {}
        for z in {it[0] for it in items}:
            if fast_scan:
                zr = mpc(z)   # skip mp Newton refinement
            else:
                zr = self._refine_root_mp(mpc(z), Om)
            sols = self._series_mp(zr, Om)
            A = self._amp_mp(zr, Om, sols)
            cache[z] = (zr, sols, A)

        # ---- BC-driven theta-edge assembly (Phase-0 Step 3) -------------------
        # Each theta-edge contributes, by kind:
        #   FREE    -> integral(trial-stress . test-disp) + corner(coeff -2)
        #   CLAMPED -> integral(trial-disp . test-stress) + corner(coeff +1)
        # Cantilever bc = [FREE@+Theta, CLAMPED@-Theta] reproduces the original
        # K = s + jm - 2*jp exactly (verified by verify_step3.py).
        edges = self.bc.edges
        ne = len(edges)
        corner_coeff = [_OOP_CORNER_COEFF[e.kind] for e in edges]
        edge_free    = [e.kind == EdgeKind.FREE for e in edges]
        # Edge-orientation factor (2026-07-08 free-free fix, extended
        # 2026-07-28 -- "Geo-2", E.2/LESSONS_LEARNED Sec 38/46): the boundary
        # functional arises from [.]_{-Theta}^{+Theta}, so ANY edge's coded
        # form (originally derived at -Theta) must flip sign when
        # transplanted to +Theta -- this is a geometric fact about which
        # wall the edge sits at, not a property of FREE vs CLAMPED. The
        # previous CLAMPED-always-+1 rule (B8-era) was correct only because
        # every validated geometry to date clamped at -Theta (orient=1
        # trivially); it was never actually derived for CLAMPED@+Theta, and
        # per the E.2 "reverse cantilever" investigation does not hold
        # there. Fixed to orient = e.sign uniformly, matching the FREE rule.
        # CLAMPED's pairing term below also needed a matching sign split
        # (sA=+1, sB=-1) -- selected empirically from a 4-way scan (job
        # 2330901) after this orient rule alone did not suffice (job
        # 2334118); the combination passes the real 41-point mirror-
        # symmetry acceptance test exactly (job 2334147) and the full
        # 20-row published OOP cantilever regression (job 2335764, 19/20
        # OK + 1 root-caused as a probe-window artifact, true error 1.3%;
        # LESSONS_LEARNED Sec 46). NOTE: InPlaneSolver._build_K_real below
        # has its OWN, separate _orient/pairing construction (its CLAMPED
        # term also carries a Lagrange constraint) -- this fix is OOP-only
        # and does not touch it; IP's own CLAMPED-orientation question is
        # unexamined and NOT assumed to carry over.
        _orient = [ e.sign for e in edges ]

        hp = mp.pi/2

        # W, Wp, Wpp at quadrature nodes and (W,Wp) at arcs per distinct root
        # (unchanged; shared across edges).
        node_vals = {}
        arc_vals  = {}
        for z in {it[0] for it in items}:
            xi, sols, A = cache[z]
            node_vals[z] = []
            for x in self.nodes:
                W   = self._W_mp(x, sols, A, 0)
                Wp  = self._W_mp(x, sols, A, 1)
                Wpp = self._W_mp(x, sols, A, 2)
                node_vals[z].append((W, Wp, Wpp))
            arc_vals[z] = []
            for x in (-hp, hp):
                Wc  = self._W_mp(x, sols, A, 0)
                Wpc = self._W_mp(x, sols, A, 1)
                arc_vals[z].append((Wc, Wpc))

        # Per item, per edge: projected stress (a=Qy, b=Tyy) and displacement
        # (A=W, B=W/rr) node vectors, plus complex Kirchhoff corner factors
        # (cM=Mry.ang, cW=W.ang) at [inner, outer] arc, using THIS edge's sign.
        edata = []
        for (z, q, P) in items:
            xi, sols, A = cache[z]
            xi2 = xi**2; ph = q*mp.pi/2
            per_edge = []
            for e in edges:
                sg = e.sign
                ssg = mp.sin(sg*xi*Th + ph); csg = mp.cos(sg*xi*Th + ph)
                a_=[]; b_=[]; A_=[]; B_=[]
                for k, x in enumerate(self.nodes):
                    rr = x + r0
                    # _coef_eval is a no-op unless the material is
                    # radially_graded (Paper B) -- see that helper's own
                    # docstring for why pointwise substitution here is exact.
                    nu_x = _coef_eval(nu, x); T_x = _coef_eval(T, x); R_x = _coef_eval(R, x)
                    W, Wp, Wpp = node_vals[z][k]
                    Tyy = nu_x*Wpp + R_x*Wp/rr - R_x*xi2*W/rr**2
                    Qy  = ((2*T_x-nu_x)*Wpp/rr + (R_x-2*T_x+2*nu_x)*Wp/rr**2
                           + (2*T_x-2*nu_x-R_x*xi2)*W/rr**3)
                    a_.append(mp_proj(xi*csg*Qy,   P))   # stress Qy
                    b_.append(mp_proj(ssg*Tyy,     P))   # stress Tyy
                    A_.append(mp_proj(ssg*W,       P))   # disp   W
                    B_.append(mp_proj(xi*csg*W/rr, P))   # disp   W/rr
                cM=[]; cW=[]
                for arc_idx, x in enumerate((-hp, hp)):
                    rr = x + r0
                    T_x = _coef_eval(T, x); nu_x = _coef_eval(nu, x)
                    Wc, Wpc = arc_vals[z][arc_idx]
                    # SOLID-DISK GUARD (2026-07-27, geometry-generalization
                    # future work, LESSONS_LEARNED Sec 37.4/38): this corner
                    # curvature combination (1/r)dW/dr - W/r^2 has a finite
                    # limit for any REGULAR Frobenius solution as r->0
                    # (W ~ a0 + a2*r^2 + ... near the origin), but the bare
                    # division is undefined exactly at the inner arc's r=0
                    # (rr==0), which occurs only for the not-yet-supported
                    # solid-disk limit (Ri=0 -- confirmed by direct
                    # reproduction, job 2334126: r0_bar collapses to exactly
                    # hp=pi/2 there, so rr=x+r0_bar hits bit-exact 0 at
                    # x=-hp). Guarded rather than left to raise
                    # ZeroDivisionError. DERIVED EXACT LIMIT (2026-08-01,
                    # LESSONS_LEARNED Sec 68.6, `Small-Bounded.txt`): the
                    # Laurent expansion of (1/r)dW/dr - W/r^2 is
                    # -a0/r^2 + 0/r + a2 + O(r) -- the 1/r term vanishes
                    # identically for every branch, so the finite part as
                    # r->0 is exactly a2 (the r^2 Taylor coefficient), NOT
                    # the earlier placeholder's guessed 2*a2. sols[q][2] is
                    # that coefficient, already computed by _series_mp and
                    # already in scope here. Odd-q branches (q=1,3) have
                    # a2=0 identically, so this is a correct no-op for them;
                    # even-q branches (q=0,2) generally do not. Uses exact
                    # equality, not a magnitude threshold: every currently-
                    # validated annular geometry has Ri>0, so rr is never
                    # bit-exactly 0 there and this is a provable no-op for
                    # all of them (SOLVER_VERSION not bumped).
                    if rr == 0:
                        Mc = (T_x-nu_x)*sols[q][2]
                    else:
                        Mc = (T_x-nu_x)*(Wpc/rr - Wc/rr**2)
                    cM.append(xi*csg*Mc)
                    cW.append(ssg*Wc)
                per_edge.append((a_, b_, A_, B_, cM, cW))
            edata.append(per_edge)

        K = matrix(n, n)
        for i in range(n):
            Pi = items[i][2]
            for j in range(n):
                Pj = items[j][2]
                # Integral part: per node, sum the per-edge weak pairing
                # (FREE: stress_j.disp_i ; CLAMPED: disp_j.stress_i).  Node-outer
                # / edge-inner ordering matches the original interleaved sum.
                s = mpf(0)
                for k, wt in enumerate(self.wts):
                    acc = mpf(0)
                    for ei in range(ne):
                        aj,bj,Aj,Bj,_,_ = edata[j][ei]
                        ai,bi,Ai,Bi,_,_ = edata[i][ei]
                        if edge_free[ei]:
                            acc += _orient[ei]*(aj[k]*Ai[k] - bj[k]*Bi[k])
                        else:
                            # Geo-2 CLAMPED pairing (sA=+1, sB=-1 relative to
                            # the FREE-edge sign convention, folded directly
                            # into the "+" here) -- see the _orient comment
                            # above for provenance. Was: "Aj[k]*ai[k] -
                            # Bj[k]*bi[k]" with no _orient multiplier
                            # (equivalent to orient forced to +1, sB=+1).
                            acc += _orient[ei]*(Aj[k]*ai[k] + Bj[k]*bi[k])
                    s += wt*acc
                # Kirchhoff corner part: per edge, coeff * jump(proj_j(Mry.ang) *
                # proj_i(W.ang)) over [inner -> outer].
                corner = mpf(0)
                for ei in range(ne):
                    _,_,_,_,cMj,_ = edata[j][ei]
                    _,_,_,_,_,cWi = edata[i][ei]
                    fj_in  = mp_proj(cMj[0], Pj); fj_out = mp_proj(cMj[1], Pj)
                    gi_in  = mp_proj(cWi[0], Pi); gi_out = mp_proj(cWi[1], Pi)
                    corner += _orient[ei]*corner_coeff[ei]*(fj_out*gi_out - fj_in*gi_in)
                K[i, j] = s + corner

        return K, n


    def logdet(self, Om_f, brs, fast_scan=False):
        K, n = self._build_K_real(Om_f, brs, fast_scan=fast_scan)
        return equilibrated_logdet(K, n)

    def sigma_min(self, Om_f, brs, fast_scan=False):
        """log10 of the smallest singular value of K (mode indicator).
        σ_min → 0 (log10 → very negative) exactly at a natural frequency."""
        K, n = self._build_K_real(Om_f, brs, fast_scan=fast_scan)
        return sigma_min_from_K(K, n)

    def _golden_sigmin(self, a, b, n_dofs, xi_max, iters,
                       n_quad_iter=None, seed_brs=None):
        """Golden-section minimise log10 σ_min(K) over [a, b] using a FRESH
        branch set at every evaluation (σ_min depth is basis-sensitive, so a
        tracked basis cannot be trusted here).  Returns (Ω*, log10 σ_min)."""
        g = (np.sqrt(5) - 1) / 2

        def f(Om):
            raw = full_search(self.fast, Om, xmax=xi_max)
            sel, cnt = select_fill(raw, n_dofs)
            return float(self.sigma_min(Om, sel, fast_scan=False))

        x1 = b - g * (b - a);  x2 = a + g * (b - a)
        f1, f2 = f(x1), f(x2)
        for _ in range(iters):
            if f1 < f2:
                b, x2, f2 = x2, x1, f1
                x1 = b - g * (b - a);  f1 = f(x1)
            else:
                a, x1, f1 = x1, x2, f2
                x2 = a + g * (b - a);  f2 = f(x2)
        Om_star = (a + b) / 2
        # Final fresh evaluation (one full_search; same cost as f(Om_star)) and
        # branch-degeneracy guard against the geometry-invariant artifact.
        raw_s = full_search(self.fast, Om_star, xmax=xi_max)
        sel_s, _ = select_fill(raw_s, n_dofs)
        s_star = float(self.sigma_min(Om_star, sel_s, fast_scan=False))
        s_star = _degen_guard(sel_s, s_star)
        return (Om_star, s_star)

    # ---------- helper: pack solver parameters for worker dispatch -----------
    def _p1_worker_params(self):
        """Return (g_scalars, m_scalars) tuples for worker reconstruction."""
        if self.mat.radially_graded:
            raise NotImplementedError(
                "worker-pool reconstruction (_p1_worker_params / _ScalarMat) "
                "does not yet support radially-graded materials -- only "
                "scalar T/R/nu are carried through the worker tuple. Graded "
                "runs must use the single-process solver path directly "
                "(matching this project's established pattern for new, "
                "not-yet-worker-integrated capabilities -- see "
                "FUTURE_WORK.md's FGM item).")
        g = (float(self.geom.r0_bar), float(self.geom.Theta),
             float(self.geom.b),      float(self.geom.r_0),
             float(self.geom.R_i),    float(self.geom.R_o),
             float(self.geom.h))
        # Pack the OOP Poisson constant self.nu (= oop_constants()[2] = μ_θ =
        # D₁₂/D₁₁), NOT mat.nu_bar — equal for isotropic (byte-identical) but
        # different for orthotropy, where nu_bar=μ_r would mis-solve workers
        # (LESSONS §1/§3).
        m = (float(self.nu),              float(self.mat.T),
             float(self.mat.R),           float(self.mat.c11_bar),
             float(self.mat.c11_eff_bar), float(self.mat.c66),
             float(self.mat.rho))
        return g, m

    # ---------- frequency search (parallel scan + parallel polish) -----------
    def find_natural_frequencies(self, Omega_range=(0.02, 0.35), n_scan=67,
                                  n_dofs=16, xi_max=14.0, dip_orders=2.0,
                                  polish_iters=8, n_quad_iter=None,
                                  n_modes_wanted=3,
                                  verbose=True, n_workers=None):
        """Find Part 1 natural frequencies with parallel scan and polish.

        FULLY GENERAL (Research22+): no paper-table reference values are
        used anywhere in this method.  Mode discovery and acceptance rely
        only on:
          (a) a dense coarse + fine scan covering the ENTIRE Omega_range
              (the fine scan is no longer anchored to a known-mode window —
              it now sweeps the full range with sequential, small-step
              branch tracking, Fix H from the lessons-learned doc, applied
              everywhere rather than only near references),
          (b) a physics-based genuine/spurious classifier (Phase 4): a
              candidate is genuine if log|det K| has a true local minimum
              at Om_star (3-point check) AND is deep enough to represent a
              near-singular K (not a generic dip),
          (c) selection of the lowest `n_modes_wanted` genuine modes, and
          (d) an iterative gap-refinement pass (Phase 5) that re-scans more
              densely between accepted modes (and below the lowest one)
              whenever fewer than n_modes_wanted genuine modes were found —
              this replaces the old reference-anchored "Phase 3c" targeted
              rescan with a purely geometric criterion (mode count vs.
              n_modes_wanted), so it generalizes to any geometry.

        n_workers : override N_WORKERS for this call (None → use global).
        """
        # Step 6: thread this solver's boundary condition to every worker pool
        # via the inherited R40_BC_KIND env var (workers rebuild it with
        # _worker_bc()).  Covers ALL pools, including the shared scan/refine
        # helpers; defaults to clamped_free so cantilever runs are unchanged.
        if self.bc.token not in _BC_REGISTRY_TOKENS:
            raise NotImplementedError(
                "parallel find_natural_frequencies has no registered worker BC "
                "for boundary=%r (token=%r)" % (self.bc, self.bc.token))
        os.environ['R40_BC_KIND'] = self.bc.token
        nw = n_workers if n_workers is not None else N_WORKERS
        t0 = time.time()
        Om_lo, Om_hi = Omega_range
        print(f"\n{'='*68}")
        print(f"[Part 1 – Out-of-Plane]  Ω̃ ∈ [{Om_lo:.4f},{Om_hi:.4f}]")
        print(f"  n_dofs={n_dofs}, M={self.M}, dps={mp.dps}, k̄={self.k_bar:.4f}")
        print(f"  N_WORKERS={nw}  n_modes_wanted={n_modes_wanted}")
        print(f"{'='*68}")

        Oms = np.linspace(Om_lo, Om_hi, n_scan)

        # Fine scan: dense rescan with sequential, small-step branch tracking
        # (Fix H) — now covers the FULL [Om_lo, Om_hi] range instead of a
        # reference-anchored sub-window.  This is the key generalization:
        # narrow sign-change/dip features (e.g. a 0.001-wide window) can
        # occur anywhere in the range for an unknown geometry, not just near
        # a known mode.  Resolution is capped by MAX_FINE; for wide ranges
        # the step widens proportionally (same trade-off as before, but now
        # applied uniformly rather than only to a sub-range).
        coarse_step = (Om_hi - Om_lo) / max(n_scan - 1, 1)
        MAX_FINE    = max(40, MAX_SCAN_PTS // 2)   # at least 40 fine points

        fine_lo   = Om_lo
        fine_hi   = Om_hi
        fine_step = min(coarse_step * 0.4, 0.002)

        # Enforce point cap: widen fine_step if needed
        raw_n_fine = int((fine_hi - fine_lo) / fine_step) + 1
        if raw_n_fine > MAX_FINE:
            fine_step = (fine_hi - fine_lo) / (MAX_FINE - 1)

        # ── Phase 1: sequential float64 branch continuation ──────────────────
        print("  Phase 1: branch tracking (sequential, fast) …", flush=True)
        prev = full_search(self.fast, float(Oms[0]), xmax=xi_max)
        branch_sets = []
        signatures  = []
        for Om in Oms:
            prev = track(self.fast, float(Om), prev)
            sel, cnt = select_fill(prev, n_dofs)
            if cnt < n_dofs:
                prev = full_search(self.fast, float(Om), xmax=xi_max)
                sel, cnt = select_fill(prev, n_dofs)
            branch_sets.append(list(sel))
            signatures.append(set_signature(sel))

        # ── Phase 2: parallel mp logdet evaluations ───────────────────────────
        # Coarse scan args: use n_quad=15 for scan phase (fast but reliable signs)
        # Polish args will use the full n_quad=30 via p_args below.
        g_sc, m_sc = self._p1_worker_params()
        dps = mp.dps
        n_quad_scan = max(15, len(self.nodes) // 2)  # 15 nodes for scan, full for polish
        coarse_args = [
            (float(Om), list(brs), g_sc, m_sc, self.M, n_quad_scan, dps)
            for Om, brs in zip(Oms, branch_sets)
        ]

        # Fine scan: track branch sets sequentially along the FINE grid itself (step ≤
        # 0.002), not seeded from the nearest coarse point — one Newton step in track()
        # can otherwise skip a fast-moving narrow root (missed Part1 0.75π m1, 1.25π
        # m1&3).  See LESSONS_LEARNED.md §Detector.
        fine_Oms = np.arange(fine_lo, fine_hi + fine_step * 0.5, fine_step)
        coarse_Oms = list(Oms)
        fine_args = []
        if len(fine_Oms) > 0:
            # Start the fine walk from the nearest coarse point's seed (just
            # for the initial condition), then track sequentially along the
            # fine grid with small steps.
            idx0 = int(np.argmin([abs(float(fine_Oms[0]) - c) for c in Oms]))
            fine_prev = list(branch_sets[idx0])
            for fine_Om in fine_Oms:
                fine_prev = track(self.fast, float(fine_Om), fine_prev)
                sel_fine, cnt_fine = select_fill(fine_prev, n_dofs)
                if cnt_fine < n_dofs:
                    fine_prev = full_search(self.fast, float(fine_Om), xmax=xi_max)
                    sel_fine, _ = select_fill(fine_prev, n_dofs)
                    fine_prev = list(sel_fine)
                fine_args.append(
                    (float(fine_Om), list(sel_fine), g_sc, m_sc, self.M, n_quad_scan, dps)
                )

        all_scan_args = coarse_args + fine_args
        n_coarse = len(coarse_args)
        n_fine   = len(fine_args)

        print(f"  Phase 2: mp logdet sweep ({n_coarse} coarse + {n_fine} fine = "
              f"{len(all_scan_args)} points, {nw} workers) …", flush=True)

        # All phases (scan + polish) share a single pool so that process
        # spawning/teardown happens only once per call to find_natural_frequencies.
        # _throttled_map keeps at most PREFETCH jobs ahead of completion,
        # preventing RAM exhaustion from pre-serialising the full argument list.
        all_scan_results = []
        t_scan = time.time()
        with concurrent.futures.ProcessPoolExecutor(max_workers=nw) as pool:
            n_total_scan = len(all_scan_args)
            n_done = 0
            for result in _throttled_map(pool, _p1_logdet_worker, all_scan_args):
                all_scan_results.append(result)
                n_done += 1
                if not verbose and n_done % max(1, n_total_scan // 10) == 0:
                    pct = 100 * n_done // n_total_scan
                    elapsed = time.time() - t_scan
                    rate = n_done / elapsed if elapsed > 0 else 0
                    eta = (n_total_scan - n_done) / rate if rate > 0 else 0
                    print(f"  scan {n_done}/{n_total_scan} ({pct}%)  "
                          f"elapsed {elapsed:.0f}s  ETA {eta:.0f}s", flush=True)

            # Split coarse and fine results; sort each by Om.
            coarse_results = sorted(all_scan_results[:n_coarse], key=lambda r: r[0])
            fine_results   = sorted(all_scan_results[n_coarse:], key=lambda r: r[0])

            rows = []
            for idx, ((Om_f, sgn, ld), sig) in enumerate(zip(coarse_results, signatures)):
                rows.append((Om_f, sgn, ld, sig))
                if verbose:
                    br = ",".join(f"{z:.2f}" for z in branch_sets[idx])
                    print(f"  {Om_f:.5f}  sign={sgn:+d}  "
                          f"log10|det|={ld:+8.3f}  [{br}]", flush=True)

            fine_rows = []
            for (Om_f, sgn, ld) in fine_results:
                fine_rows.append((Om_f, sgn, ld))
                if verbose:
                    print(f"  [fine] {Om_f:.5f}  sign={sgn:+d}  "
                          f"log10|det|={ld:+8.3f}", flush=True)

            t_scan_done = time.time() - t_scan
            print(f"  Phase 2 done: {n_coarse} coarse + {n_fine} fine pts in {t_scan_done:.1f}s", flush=True)

            # ── Phase 2b: candidate detection ─────────────────────────────────
            cands = []
            seen  = set()

            # Sign changes in the coarse scan.
            for k in range(1, len(rows)):
                Om_a, Om_b = rows[k-1][0], rows[k][0]
                key = (round(Om_a, 4), round(Om_b, 4))
                if rows[k-1][1] * rows[k][1] < 0:
                    if abs(rows[k-1][2] - rows[k][2]) < 4.0 and key not in seen:
                        cands.append((Om_a, Om_b, 'sign'))
                        seen.add(key)

            # Local minima in the coarse scan (catches modes outside fine range,
            # and modes that touch zero without producing a sign change).
            # Use the looser threshold always: the Phase 4 physics-based
            # genuine-mode filter (deep + locally-minimal log|det|) removes
            # false positives downstream, so we can afford to over-detect
            # candidates here.
            lds = [r[2] for r in rows]
            med = float(np.median(lds))
            coarse_dip_thresh = min(1.5, dip_orders)
            for k in range(1, len(rows) - 1):
                Om_a, Om_b = rows[k-1][0], rows[k+1][0]
                key = (round(Om_a, 4), round(Om_b, 4))
                if (lds[k] < lds[k-1] and lds[k] < lds[k+1]
                        and lds[k] < med - coarse_dip_thresh
                        and key not in seen):
                    cands.append((Om_a, Om_b, 'dip'))
                    seen.add(key)

            # Sign changes in the fine scan.
            for k in range(1, len(fine_rows)):
                Om_a, Om_b = fine_rows[k-1][0], fine_rows[k][0]
                key = (round(Om_a, 4), round(Om_b, 4))
                if fine_rows[k-1][1] * fine_rows[k][1] < 0:
                    if abs(fine_rows[k-1][2] - fine_rows[k][2]) < 4.0 and key not in seen:
                        cands.append((Om_a, Om_b, 'fine_sign'))
                        seen.add(key)

            # Local minima in the fine scan.  Threshold 0.6 (lowered from 0.8)
            # to catch narrow double-zero modes that barely dip below surroundings.
            fine_lds = [r[2] for r in fine_rows]
            fine_med = float(np.median(fine_lds)) if fine_lds else 0.0
            for k in range(1, len(fine_rows) - 1):
                Om_a, Om_b = fine_rows[k-1][0], fine_rows[k+1][0]
                key = (round(Om_a, 4), round(Om_b, 4))
                if (fine_lds[k] < fine_lds[k-1] and fine_lds[k] < fine_lds[k+1]
                        and fine_lds[k] < fine_med - 0.6
                        and key not in seen):
                    cands.append((Om_a, Om_b, 'fine_dip'))
                    seen.add(key)

            print(f"\n  candidates: "
                  f"{[(f'{a:.4f}', f'{b:.4f}', t) for a,b,t in cands]}")

            # ── Phase 2c: cheap pre-filter by candidate depth ───────────────
            # Cap polished candidates to the MAX_POLISH deepest (by the log|det| already
            # measured at Phase-2 brackets) plus any past DEPTH_AUTOKEEP.  No extra evals.
            ld_lookup = {}
            for (Om_f, sgn, ld, sig) in rows:
                ld_lookup[round(Om_f, 5)] = ld
            for (Om_f, sgn, ld) in fine_rows:
                ld_lookup[round(Om_f, 5)] = ld

            def _cand_depth(a, b):
                la = ld_lookup.get(round(a, 5))
                lb = ld_lookup.get(round(b, 5))
                vals = [v for v in (la, lb) if v is not None and np.isfinite(v)]
                return min(vals) if vals else 0.0

            DEPTH_AUTOKEEP = -4.0  # always keep candidates at least this deep
            MAX_POLISH = max(2 * n_modes_wanted + 4, 12)
            cand_depths = [_cand_depth(a, b) for a, b, _ in cands]
            order = sorted(range(len(cands)), key=lambda i: cand_depths[i])
            keep_idx = set()
            for rank, i in enumerate(order):
                if cand_depths[i] < DEPTH_AUTOKEEP or rank < MAX_POLISH:
                    keep_idx.add(i)
            if len(keep_idx) < len(cands):
                print(f"  Phase 2c: pre-filter kept {len(keep_idx)}/{len(cands)} "
                      f"candidates (deepest log|det|, cap={MAX_POLISH})", flush=True)
            cands = [cands[i] for i in sorted(keep_idx)]

            # ── Phase 3: candidate polishing (same pool) ───────────────────────
            # sign/fine_sign candidates have a confirmed bracket → bisection.
            # dip/fine_dip candidates have no sign change → golden-section.
            worker_fn_map = {
                'sign':      _p1_bisect_worker,
                'dip':       _p1_golden_worker,
                'fine_sign': _p1_bisect_worker,
                'fine_dip':  _p1_golden_worker,
            }
            nq_fast = (n_quad_iter if n_quad_iter is not None
                       else max(8, len(self.nodes) // 2))
            cand_types = [t for _, _, t in cands]
            cand_lo    = [a for a, _, _ in cands]
            cand_hi    = [b for _, b, _ in cands]

            # For each candidate bracket [a, b], find the nearest coarse branch set
            # to use as a Newton seed — avoids 34 s full_search in each polish worker.
            def _nearest_brs(Om_mid):
                idx = int(np.argmin([abs(Om_mid - c) for c in coarse_Oms]))
                return branch_sets[idx]

            p_args = [
                (a, b, g_sc, m_sc, self.M, len(self.nodes),
                 xi_max, n_dofs, polish_iters, nq_fast, dps,
                 _nearest_brs((a + b) / 2))
                for a, b in zip(cand_lo, cand_hi)
            ]

            print(f"  Phase 3: polishing {len(cands)} candidate(s) …", flush=True)
            polish_results = []
            t_polish = time.time()
            if p_args:
                fns = [worker_fn_map[t] for t in cand_types]
                futures = [pool.submit(fn, arg) for fn, arg in zip(fns, p_args)]
                for i, fut in enumerate(futures):
                    try:
                        Om_star, ld_star = fut.result()
                        polish_results.append(
                            (cand_types[i], cand_lo[i], cand_hi[i], Om_star, ld_star))
                        print(f"  polish {i+1}/{len(cands)}  "
                              f"Ω̃={Om_star:.6f}  logdet={ld_star:+.2f}  "
                              f"({time.time()-t_polish:.0f}s)", flush=True)
                    except Exception as exc:
                        print(f"  WARNING: polish [{cand_lo[i]:.4f},"
                              f"{cand_hi[i]:.4f}] failed: {exc}")

        # Pool is closed here; all worker processes have exited.

        # ── Phase 3b: post-polish verification for all candidates ──────────
        # Re-run branch search at each converged frequency to confirm the
        # branch count is complete; refresh ld_star with the verified value.
        verified = []  # list of dicts: typ, a, b, Om_star, ld_star
        for typ, a, b, Om_star, ld_star in polish_results:
            chk_brs = full_search(self.fast, Om_star, xmax=xi_max)
            chk_sel, chk_cnt = select_fill(chk_brs, n_dofs)
            if chk_cnt < n_dofs:
                print(f"  WARNING: only {chk_cnt}/{n_dofs} dofs at "
                      f"Omega_tilde={Om_star:.6f} — result may be unreliable")
            else:
                _, ld_star_verified = self.logdet(Om_star, chk_sel)
                ld_star = float(ld_star_verified)
            verified.append({'typ': typ, 'a': a, 'b': b,
                              'Om_star': Om_star, 'ld_star': ld_star})

        freqs = []

        # ── Phase 4: physics-based genuine-mode filter ──────────────────────
        # GENUINE = log|det K| has a prominent local min at Om_star (≥ P1_PROMINENCE
        # below both neighbours) AND depth ≤ P1_LD_GENUINE.  History: §Classifier.
        LD_GENUINE   = _P1_LD_GENUINE
        P1_PROM      = _P1_PROM
        d_probe = max(1e-3, 0.5 * fine_step)
        g_sc, m_sc = self._p1_worker_params()

        def _nearest_brs_full(Om_mid):
            idx = int(np.argmin([abs(Om_mid - c) for c in coarse_Oms]))
            return branch_sets[idx]

        filter_points = [
            Om_star + delta
            for v in verified
            for Om_star, delta in ((v['Om_star'], -d_probe),
                                    (v['Om_star'], 0.0),
                                    (v['Om_star'], +d_probe))
        ]
        filter_args = [
            (Om_pt, g_sc, m_sc, self.M, n_quad_scan, xi_max, n_dofs, dps,
             _nearest_brs_full(Om_pt))
            for Om_pt in filter_points
        ]
        filter_lds = []
        if filter_args:
            with concurrent.futures.ProcessPoolExecutor(max_workers=nw) as pool4:
                for ld_val in _throttled_map(pool4, _p1_ld_at_worker, filter_args):
                    filter_lds.append(ld_val)

        for idx, v in enumerate(verified):
            if filter_lds:
                lo_v, mid_v, hi_v = filter_lds[idx*3], filter_lds[idx*3+1], filter_lds[idx*3+2]
            else:
                lo_v = mid_v = hi_v = v['ld_star']
            # Prominent local minimum: central dip at least P1_PROM orders below
            # both shoulders.  Falls back to a plain local-min if the shoulder
            # evaluations were unavailable.
            prominent = (mid_v <= lo_v - P1_PROM and mid_v <= hi_v - P1_PROM)
            local_min = (mid_v <= lo_v and mid_v <= hi_v)
            shape_ok = prominent or (local_min and not filter_lds)
            v['genuine'] = bool(shape_ok and v['ld_star'] < LD_GENUINE)
            if verbose:
                status = 'OK' if v['genuine'] else 'SPURIOUS'
                print(f"  polished [{v['a']:.4f},{v['b']:.4f}] ({v['typ']}) → "
                      f"Omega_tilde={v['Om_star']:.6f}  "
                      f"log10|det|={v['ld_star']:+.2f}  {status}")

        # ── Mode selection ──────────────────────────────────────────────────
        # Collect genuine candidates, dedupe near-duplicates, sort by Om, and
        # keep the lowest n_modes_wanted — these are the lowest natural
        # frequencies, which is what every geometry's "first n modes" means.
        def _merge_genuine(verified_list, accepted):
            for v in verified_list:
                if not v['genuine']:
                    continue
                Om_star = v['Om_star']
                if all(abs(Om_star - f['Om_star']) > 1e-3 for f in accepted):
                    accepted.append(v)
            accepted.sort(key=lambda v: v['Om_star'])
            return accepted

        accepted = _merge_genuine(verified, [])

        # ── Phase 5: gap-refinement (replaces reference-anchored Phase 3c) ──
        # If fewer than n_modes_wanted genuine modes were found, the global
        # scan likely lost a fast-moving branch across one tracking step
        # (see lessons-learned Fix H/I).  Run a denser, fresh-seeded local
        # rescan in the "gaps" — below the lowest accepted mode, between
        # consecutive accepted modes, and up to Om_hi — purely driven by the
        # mode COUNT (a geometric criterion), never by a known answer.
        MAX_REFINE_PASSES = 3
        refine_pass = 0
        while len(accepted) < n_modes_wanted and refine_pass < MAX_REFINE_PASSES:
            refine_pass += 1
            gaps = []
            edges = [Om_lo] + [v['Om_star'] for v in accepted] + [Om_hi]
            for k in range(len(edges) - 1):
                lo_g, hi_g = edges[k], edges[k+1]
                if hi_g - lo_g > 1e-6:
                    gaps.append((lo_g, hi_g))
            if not gaps:
                break
            print(f"  Phase 5: gap-refinement pass {refine_pass} — "
                  f"{len(accepted)}/{n_modes_wanted} genuine modes found, "
                  f"rescanning {len(gaps)} gap(s) "
                  f"{[(f'{a:.4f}', f'{b:.4f}') for a,b in gaps]} …", flush=True)
            new_verified = []
            for (lo_g, hi_g) in gaps:
                new_verified.extend(
                    self._local_rescan(lo_g, hi_g, n_dofs, xi_max, polish_iters,
                                        nq_fast, n_quad_scan, dps, g_sc, m_sc,
                                        nw, verbose, label=f"gap[{lo_g:.4f},{hi_g:.4f}]"))
            if not new_verified:
                break
            n_before = len(accepted)
            accepted = _merge_genuine(new_verified, accepted)
            if len(accepted) == n_before:
                # No new genuine modes found this pass — stop to avoid
                # looping forever on a geometry that simply has fewer than
                # n_modes_wanted modes in this range.
                break

        for v in accepted:
            freqs.append(v['Om_star'])

        freqs.sort()
        freqs = freqs[:n_modes_wanted]
        print(f"  elapsed {time.time()-t0:.0f}s")
        return freqs

    # ------------------------------------------------------------------
    # Local rescan of a sub-window [lo, hi]: fresh-seeded sequential
    # tracking with small steps (Fix G/H), candidate detection, polish, and
    # the same genuine-mode 3-point check as the main scan.  Used by Phase 5
    # gap-refinement.  Returns a list of verified dicts (with 'genuine' set).
    # ------------------------------------------------------------------
    def _local_rescan(self, lo, hi, n_dofs, xi_max, polish_iters, nq_fast,
                       n_quad_scan, dps, g_sc, m_sc, nw, verbose, label=""):
        n_pts = max(15, min(30, MAX_SCAN_PTS // 2))
        t_Oms = np.linspace(lo, hi, n_pts)
        prev_seed = full_search(self.fast, float(t_Oms[0]), xmax=xi_max)
        t_rows, t_brs = [], []
        for Om_t in t_Oms:
            prev_seed = track(self.fast, float(Om_t), prev_seed)
            sel_t, cnt_t = select_fill(prev_seed, n_dofs)
            if cnt_t < n_dofs:
                prev_seed = full_search(self.fast, float(Om_t), xmax=xi_max)
                sel_t, _ = select_fill(prev_seed, n_dofs)
            sgn_t, ld_t = self._logdet_nq(float(Om_t), sel_t, n_quad_scan, fast_scan=False)
            t_rows.append((float(Om_t), float(sgn_t), float(ld_t)))
            t_brs.append(list(sel_t))
            if verbose:
                print(f"    [{label}] {Om_t:.5f}  sign={int(sgn_t):+d}  "
                      f"log10|det|={ld_t:+8.3f}", flush=True)

        cands = []
        for k in range(1, len(t_rows)):
            if t_rows[k-1][1] * t_rows[k][1] < 0:
                cands.append((t_rows[k-1][0], t_rows[k][0], 'sign', k))
        lds_t = [r[2] for r in t_rows]
        med_t = float(np.median(lds_t))
        for k in range(1, len(t_rows) - 1):
            if (lds_t[k] < lds_t[k-1] and lds_t[k] < lds_t[k+1]
                    and lds_t[k] < med_t - 1.5):
                cands.append((t_rows[k-1][0], t_rows[k+1][0], 'dip', k))

        verified = []
        for (a_t, b_t, typ, k) in cands:
            seed_mid = t_brs[k]
            if typ == 'sign':
                Om_star_t, ld_star_t = self._bisect(
                    a_t, b_t, n_dofs, xi_max, polish_iters,
                    n_quad_iter=nq_fast, seed_brs=seed_mid)
            else:
                Om_star_t, ld_star_t = self._golden(
                    a_t, b_t, n_dofs, xi_max, polish_iters,
                    n_quad_iter=nq_fast, seed_brs=seed_mid)

            chk_brs = full_search(self.fast, Om_star_t, xmax=xi_max)
            chk_sel, chk_cnt = select_fill(chk_brs, n_dofs)
            if chk_cnt < n_dofs:
                print(f"  WARNING: only {chk_cnt}/{n_dofs} dofs at "
                      f"Omega_tilde={Om_star_t:.6f} ({label}) — "
                      f"result may be unreliable")
            else:
                _, ld_v = self.logdet(Om_star_t, chk_sel)
                ld_star_t = float(ld_v)

            # 3-point genuine check, using a fresh local seed.
            d_probe = max(1e-3, 0.5 * (hi - lo) / n_pts)
            seeds = {}
            for delta in (-d_probe, 0.0, d_probe):
                Om_p = Om_star_t + delta
                p_seed = track(self.fast, Om_p, [complex(z) for z in seed_mid])
                sel_p, cnt_p = select_fill(p_seed, n_dofs)
                if cnt_p < n_dofs:
                    p_seed = full_search(self.fast, Om_p, xmax=xi_max)
                    sel_p, _ = select_fill(p_seed, n_dofs)
                seeds[delta] = float(self.logdet(Om_p, sel_p, fast_scan=False)[1])
            local_min = (seeds[0.0] <= seeds[-d_probe] and seeds[0.0] <= seeds[d_probe])
            prominent = (seeds[0.0] <= seeds[-d_probe] - _P1_PROM and
                         seeds[0.0] <= seeds[d_probe] - _P1_PROM)
            genuine = bool((prominent or local_min) and ld_star_t < _P1_LD_GENUINE)

            verified.append({'typ': f'local_{typ}', 'a': a_t, 'b': b_t,
                              'Om_star': Om_star_t, 'ld_star': ld_star_t,
                              'genuine': genuine})
            if verbose:
                status = 'OK' if genuine else 'SPURIOUS'
                print(f"  local [{label}] [{a_t:.4f},{b_t:.4f}] ({typ}) → "
                      f"Omega_tilde={Om_star_t:.6f}  "
                      f"log10|det|={ld_star_t:+.2f}  {status}", flush=True)
        return verified

    # ------------------------------------------------------------------
    # Evaluate log|det| at an arbitrary Omega (for genuine-mode 3-point
    # checks).  Uses track() from seed_brs with a full_search fallback.
    # ------------------------------------------------------------------
    def _ld_at(self, Om, n_dofs, xi_max, seed_brs=None):
        if seed_brs:
            raw = track(self.fast, Om, [complex(z) for z in seed_brs])
            if not raw:
                raw = full_search(self.fast, Om, xmax=xi_max)
        else:
            raw = full_search(self.fast, Om, xmax=xi_max)
        sel, cnt = select_fill(raw, n_dofs)
        if cnt < n_dofs:
            raw = full_search(self.fast, Om, xmax=xi_max)
            sel, _ = select_fill(raw, n_dofs)
        return float(self.logdet(Om, sel, fast_scan=False)[1])

    # ------------------------------------------------------------------
    # Fast-quadrature logdet: uses a coarser node set during polishing
    # iterations so each evaluation is ~(n_quad_full/n_quad_fast)^? faster.
    # The full-resolution logdet is only called once at the final point.
    # ------------------------------------------------------------------
    def _logdet_nq(self, Om_f, brs, n_quad_override=None, fast_scan=False):
        """logdet with optional quadrature-node override (for cheap iterations).
        fast_scan=True skips mp root refinement — suitable for bisect/golden iterations.
        """
        if n_quad_override is None or n_quad_override == len(self.nodes):
            return self.logdet(Om_f, brs, fast_scan=fast_scan)
        nd, wt = mp_gauss_legendre(n_quad_override)
        saved_nodes, saved_wts = self.nodes, self.wts
        self.nodes = [x * mp.pi / 2 for x in nd]
        self.wts   = [w * mp.pi / 2 for w in wt]
        try:
            result = self.logdet(Om_f, brs, fast_scan=fast_scan)
        finally:
            self.nodes, self.wts = saved_nodes, saved_wts
        return result

    def _bisect(self, a, b, n_dofs, xi_max, iters, n_quad_iter=None, seed_brs=None):
        """Bisect sign-change bracket.

        seed_brs : optional float64 branch list from the coarse scan to use as
                   Newton seed instead of calling full_search (saves 34 s per call).
        Iterations use fast_scan=True; only the final evaluation uses full precision.
        """
        n_quad_fast = n_quad_iter if n_quad_iter is not None else max(8, len(self.nodes) // 2)

        # Seed: use provided coarse branches (track from them) or fall back to full_search
        if seed_brs:
            prev_a = track(self.fast, a, [complex(z) for z in seed_brs])
            if not prev_a:
                prev_a = full_search(self.fast, a, xmax=xi_max)
            prev_b = track(self.fast, b, [complex(z) for z in seed_brs])
            if not prev_b:
                prev_b = full_search(self.fast, b, xmax=xi_max)
        else:
            prev_a = full_search(self.fast, a, xmax=xi_max)
            prev_b = full_search(self.fast, b, xmax=xi_max)

        sel_a, cnt = select_fill(prev_a, n_dofs)
        if cnt < n_dofs:
            prev_a = full_search(self.fast, a, xmax=xi_max)
            sel_a, _ = select_fill(prev_a, n_dofs)
        sel_b, cnt = select_fill(prev_b, n_dofs)
        if cnt < n_dofs:
            prev_b = full_search(self.fast, b, xmax=xi_max)
            sel_b, _ = select_fill(prev_b, n_dofs)

        sa, _ = self._logdet_nq(a, sel_a, n_quad_fast, fast_scan=False)
        prev_lo, prev_hi = prev_a, prev_b

        for _ in range(iters):
            m = (a + b) / 2
            if abs(m - a) <= abs(m - b):
                prev_m = track(self.fast, m, prev_lo)
            else:
                prev_m = track(self.fast, m, prev_hi)
            sel_m, cnt = select_fill(prev_m, n_dofs)
            if cnt < n_dofs:
                prev_m = full_search(self.fast, m, xmax=xi_max)
                sel_m, _ = select_fill(prev_m, n_dofs)
            sm_, _ = self._logdet_nq(m, sel_m, n_quad_fast, fast_scan=False)
            if sm_ == sa:
                a = m; prev_lo = prev_m
            else:
                b = m; prev_hi = prev_m

        m = (a + b) / 2
        prev_m = track(self.fast, m, prev_lo)
        sel_m, cnt = select_fill(prev_m, n_dofs)
        if cnt < n_dofs:
            prev_m = full_search(self.fast, m, xmax=xi_max)
            sel_m, _ = select_fill(prev_m, n_dofs)
        _, ld_final = self.logdet(m, sel_m, fast_scan=False)   # full precision
        return (m, float(ld_final))

    def _golden(self, a, b, n_dofs, xi_max, iters, n_quad_iter=None, seed_brs=None):
        """Golden-section minimise log|det K|.
        seed_brs : coarse branch list to track from instead of full_search.
        Iterations use fast_scan=True; final evaluation uses full precision.
        """
        n_quad_fast = n_quad_iter if n_quad_iter is not None else max(8, len(self.nodes) // 2)
        g = (np.sqrt(5) - 1) / 2

        if seed_brs:
            prev_x1 = track(self.fast, a, [complex(z) for z in seed_brs]) or \
                       full_search(self.fast, a, xmax=xi_max)
            prev_x2 = track(self.fast, b, [complex(z) for z in seed_brs]) or \
                       full_search(self.fast, b, xmax=xi_max)
        else:
            prev_x1 = full_search(self.fast, a, xmax=xi_max)
            prev_x2 = full_search(self.fast, b, xmax=xi_max)

        def f_tracked(Om, prev_seed):
            prev = track(self.fast, Om, prev_seed)
            sel, cnt = select_fill(prev, n_dofs)
            if cnt < n_dofs:
                prev = full_search(self.fast, Om, xmax=xi_max)
                sel, _ = select_fill(prev, n_dofs)
            ld = float(self._logdet_nq(Om, sel, n_quad_fast, fast_scan=False)[1])
            return ld, prev

        x1 = b - g * (b - a);  x2 = a + g * (b - a)
        f1, prev_x1 = f_tracked(x1, prev_x1)
        f2, prev_x2 = f_tracked(x2, prev_x2)

        for _ in range(iters):
            if f1 < f2:
                b, x2, f2, prev_x2 = x2, x1, f1, prev_x1
                x1 = b - g * (b - a)
                f1, prev_x1 = f_tracked(x1, prev_x2)
            else:
                a, x1, f1, prev_x1 = x1, x2, f2, prev_x2
                x2 = a + g * (b - a)
                f2, prev_x2 = f_tracked(x2, prev_x1)

        Om_star = (a + b) / 2
        prev_f = track(self.fast, Om_star, prev_x1)
        sel_f, cnt = select_fill(prev_f, n_dofs)
        if cnt < n_dofs:
            prev_f = full_search(self.fast, Om_star, xmax=xi_max)
            sel_f, _ = select_fill(prev_f, n_dofs)
        _, ld_final = self.logdet(Om_star, sel_f, fast_scan=False)   # full precision
        return (Om_star, float(ld_final))

    # ------------------------------------------------------------------
    # Mode shape reconstruction (Figure 3 style).  Builds K at the
    # converged (Om_star, brs), finds its null vector c via SVD, then
    # evaluates W(r_tilde, theta) = sum_j c_j * phi_j(r_tilde, theta) on a
    # grid for plotting.  Best-effort: returns None on any failure so a
    # plotting error never aborts the main solve.
    # ------------------------------------------------------------------
    def mode_shape_grid(self, Om_star, brs, n_dofs, n_r=25, n_th=25):
        try:
            sel, cnt = select_fill(brs, n_dofs)
            if cnt < n_dofs:
                return None
            K, n = self._build_K_real(Om_star, sel, fast_scan=False)

            items = real_basis_items(sel)
            # Pre-build per-item (xi, sols, A) using fast (non-refined) roots
            # — sufficient for a mode-shape PLOT (visual accuracy only).
            cache = {}
            for z in {it[0] for it in items}:
                zr = mpc(z)
                sols = self._series_mp(zr, mpf(Om_star))
                A = self._amp_mp(zr, mpf(Om_star), sols)
                cache[z] = (zr, sols, A)

            hp = float(mp.pi/2)
            Th = float(self.Theta)

            # ── Mode coefficients c (null vector of K) ───────────────────────
            # Smallest right singular vector of K at mp precision with column
            # equilibration (float64 SVD is wrecked by the evanescent dynamic range).
            # MODE_SHAPE_MP_NULLVEC=0 forces float64.  Weak-clamp note: §Mode-shapes.
            c = None
            if os.environ.get("MODE_SHAPE_MP_NULLVEC", "1") == "1":
                c = equilibrated_nullvec_mp(K, n)
            if c is None:
                try:
                    Kn = np.array([[complex(K[i, j]) for j in range(n)]
                                    for i in range(n)], dtype=np.complex128).real
                    _, _, Vt = np.linalg.svd(Kn)
                    c = np.asarray(Vt[-1, :], dtype=float)
                except np.linalg.LinAlgError:
                    # float64 SVD didn't converge on this (typically
                    # ill-conditioned) K. Surface it once instead of letting
                    # it fall through indistinguishably from unrelated
                    # earlier failures in this method's outer except.
                    print(f"  [mode_shape_grid] float64 SVD fallback did not "
                          f"converge at Om_star={Om_star}; skipping this mode shape.")
                    return None

            # ── Boundary-collocation refinement (experimental; Problem 2) ────
            # Tries coefficients that satisfy the θ-edge BCs pointwise; accepts only if it
            # flattens the clamp without the spike, else keeps the weak-form shape (never
            # worse).  Self-rejects with the current small basis.  MODE_SHAPE_COLLOCATION=1.
            if os.environ.get("MODE_SHAPE_COLLOCATION", "0") == "1":
                try:
                    cl0, _in0 = _quick_shape_stats(self, items, cache, c)
                    best = None; best_cl = cl0
                    for vmode in ('sin', 'xicos'):
                        cc = _theta_collocation_coeffs(self, Om_star, items,
                                                       cache, n, v_theta=vmode)
                        if cc is None:
                            continue
                        cl, inn = _quick_shape_stats(self, items, cache, cc)
                        # accept only a clear clamp improvement with no spike
                        if cl < 0.25 and inn > 0.30 and cl < best_cl - 0.05:
                            best = cc; best_cl = cl
                    if best is not None:
                        c = best
                except Exception:
                    pass

            r_vals = np.linspace(-hp, hp, n_r)
            th_vals = np.linspace(-Th, Th, n_th)
            Wgrid = np.zeros((n_r, n_th))
            for ir, x in enumerate(r_vals):
                # Evaluate W_j(x) for each distinct root once per radial point.
                Wcache = {}
                for z in {it[0] for it in items}:
                    xi, sols, A = cache[z]
                    Wcache[z] = complex(self._W_mp(mpf(x), sols, A, 0))
                for it_idx, (z, q, P) in enumerate(items):
                    Wj = Wcache[z]
                    for jt, th in enumerate(th_vals):
                        xi, sols, A = cache[z]
                        ph = q*float(mp.pi)/2
                        # angular factor: sin(xi*theta+ph) -- matches the
                        # DISPLACEMENT projection actually used in
                        # _build_K_real's validated edge assembly (ssg =
                        # sin(sg*xi*Th+ph) multiplying W for both the
                        # A_ node vector and Tyy; see B18). This was
                        # cos(xi*theta+ph) before 2026-07-15, giving every
                        # mode_shape_grid reconstruction the WRONG parity
                        # relative to items' actual q assignment (figures
                        # only -- this never touched _build_K_real/
                        # full_search/select_fill, so no frequency this
                        # project has published was ever affected).
                        ang_val = np.sin(complex(xi)*th + ph)
                        phi = Wj * ang_val
                        val = phi.real if P in ('re', 'asis') else phi.imag
                        Wgrid[ir, jt] += c[it_idx] * val

            # Normalize for plotting (paper figures show unit-scale mode
            # shapes; the eigenvector is only defined up to a scale factor).
            mx = np.max(np.abs(Wgrid))
            if mx > 0:
                Wgrid = Wgrid / mx
            # Physical planform.  r_tilde in [-pi/2,pi/2] maps to the physical
            # radius r in [R_i, R_o] via r = r_0 + (2b/pi)*r_tilde, and theta is
            # kept in TRUE radians [-Theta,+Theta].  The old code returned a
            # normalised radius in [0,1] and a normalised angle in [-1,1], which
            # made plot_mode_shape draw a unit PIE wedge from the origin at a
            # fixed ~114° spread for every geometry — not the actual annular
            # (ring) sector of width 2Theta the paper shows in Fig. 3.
            r0_phys = float(self.geom.r_0)
            b_phys  = float(self.geom.b)
            r_phys  = r0_phys + (2.0 * b_phys / float(mp.pi)) * r_vals  # [R_i,R_o]
            th_phys = th_vals                                           # radians
            return r_phys, th_phys, Wgrid
        except Exception:
            return None


# ══════════════════════════════════════════════════════════════════════════════
#  SECTION 5 — PART 2: IN-PLANE SOLVER (mp precision K̂ / det)
# ══════════════════════════════════════════════════════════════════════════════

class InPlaneSolver:
    """In-plane vibration, Seok & Tiersten Part 2.  Frequencies are Ω̄."""
    def __init__(self, geom, mat, M=80, n_quad=30, boundary=None):
        self.geom = geom; self.mat = mat; self.M = M
        self.r0b = geom.r0_bar
        # theta-edge boundary condition (Phase-0 Step 3).  Defaults to the
        # clamped-free cantilever so the worker reconstruction path (which
        # builds the solver with no `boundary` arg) is unchanged -- same
        # safe-default pattern as radially_graded=False.
        self.bc = boundary if boundary is not None else ClampedFreeIP()
        # Pull the IP ODE constants through the MaterialModel interface (Phase-0
        # refactor); identical values to mat.c11_eff_bar, mat.R, mat.nu_bar for
        # IsotropicMaterial, so frequencies are unchanged. A radially_graded
        # material (RadialFGMMaterial, LESSONS_LEARNED Sec 54.4) supplies
        # Taylor-coefficient LISTS via ip_coeff_series() instead of scalars --
        # mirrors OutOfPlaneSolver's own oop_coeff_series()/oop_constants()
        # branch above.
        if mat.radially_graded:
            self.c11, self.R, self.nu = mat.ip_coeff_series(self.r0b, M)
        else:
            self.c11, self.R, self.nu = mat.ip_constants()
        self.Theta = geom.Theta
        # float() cast only for the non-graded (scalar) case -- forcing a
        # graded material's c11/R/nu (Taylor-coefficient LISTS) through
        # float() would raise; pass the list through unchanged instead, and
        # _Part2Fast's own radially_graded branch handles it (same pattern
        # as OutOfPlaneSolver's _fast_T/_fast_R/_fast_nu above).
        _fast_nu  = [float(c) for c in self.nu]  if mat.radially_graded else float(self.nu)
        _fast_c11 = [float(c) for c in self.c11] if mat.radially_graded else float(self.c11)
        _fast_R   = [float(c) for c in self.R]   if mat.radially_graded else float(self.R)
        self.fast = AnnulusRadialIP(float(self.r0b), _fast_nu,
                                    _fast_c11, _fast_R, M,
                                    radially_graded=mat.radially_graded)
        nd, wt = mp_gauss_legendre(n_quad)
        self.nodes = [x*mp.pi/2 for x in nd]
        self.wts = [w*mp.pi/2 for w in wt]

    def _series_mp(self, ze, Om):
        ze = mpc(ze); O2 = mpc(Om)**2
        r0 = self.r0b
        M = self.M
        if not isinstance(self.c11, list):
            # ---- EXACT existing constant-coefficient path, byte-for-byte
            # unchanged ----
            c11, R, nu = self.c11, self.R, self.nu
            cF = c11*R + ze**2
            cG = 1 + c11*R*ze**2
            coup = (c11*nu + 1)*ze
            cst = (c11*R + 1)*ze
            sols = []
            for q in range(4):
                a = [mpc(0)]*(M+3); b = [mpc(0)]*(M+3)
                if q == 0: a[0] = mpc(1)
                elif q == 1: a[1] = mpc(1)
                elif q == 2: b[0] = mpc(1)
                else: b[1] = mpc(1)
                for p in range(0, M+1):
                    am2 = a[p-2] if p >= 2 else mpc(0)
                    am1 = a[p-1] if p >= 1 else mpc(0)
                    bm2 = b[p-2] if p >= 2 else mpc(0)
                    bm1 = b[p-1] if p >= 1 else mpc(0)
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
        # ---- Radially-graded path (LESSONS_LEARNED Sec 54.4; derived and
        # self-tested standalone by Grok, FGMInPlane.txt, then wired in here).
        # c11/R/nu are Taylor-coefficient LISTS in x (rr=x+r0); every scalar
        # product in the constant-coefficient recursion above becomes a
        # discrete LINEAR convolution of the material series with the
        # kinematic-operator series -- explicitly NOT quadratic in any
        # material series (the exact bug already fixed once for the OOP
        # grading work, see RadialFGMMaterial's docstring). The highest-order
        # unknown (a[p+2]/b[p+2]) only ever appears multiplied by the material
        # series' OWN CONSTANT term (c11s[0] for the F-equation; the G-
        # equation's leading operator has no material multiplier at all --
        # "OpB multiplier == 1"), so the recursion stays causal with the same
        # denominator structure as the scalar path. Mirrored in float64 by
        # dispersion.py::_Part2Fast.series(); keep both in sync. ----
        c11s = [mpc(c) for c in self.c11]
        Rs   = [mpc(c) for c in self.R]
        nus  = [mpc(c) for c in self.nu] if isinstance(self.nu, list) else [mpc(self.nu)]
        N = M + 2

        def _conv_full(p, q, n):
            out = [mpc(0)] * (n + 1)
            for j in range(n + 1):
                s = mpc(0)
                for k in range(0, j + 1):
                    pk = p[k] if k < len(p) else mpc(0)
                    qk = q[j - k] if (j - k) < len(q) else mpc(0)
                    s += pk * qk
                out[j] = s
            return out

        c11R = _conv_full(c11s, Rs, N)
        c11nu = _conv_full(c11s, nus, N)
        cF   = [c11R[j] + (ze**2 if j == 0 else mpc(0)) for j in range(N + 1)]
        cG   = [(mpc(1) if j == 0 else mpc(0)) + c11R[j]*ze**2 for j in range(N + 1)]
        coup = [(c11nu[j] + (mpc(1) if j == 0 else mpc(0)))*ze for j in range(N + 1)]
        cst  = [(c11R[j]  + (mpc(1) if j == 0 else mpc(0)))*ze for j in range(N + 1)]

        def _dot(cs, xs, pp):
            # cs may be a raw (short) material Taylor list, not zero-padded
            # to N+1 -- guard both operands defensively (xs is always
            # p+1 long by construction, but the bounds check costs nothing).
            s = mpc(0)
            for k in range(pp + 1):
                ck = cs[k] if k < len(cs) else mpc(0)
                xk = xs[pp - k] if (pp - k) < len(xs) else mpc(0)
                s += ck*xk
            return s

        sols = []
        for q in range(4):
            a = [mpc(0)]*(M+3); b = [mpc(0)]*(M+3)
            if q == 0: a[0] = mpc(1)
            elif q == 1: a[1] = mpc(1)
            elif q == 2: b[0] = mpc(1)
            else: b[1] = mpc(1)
            for p in range(0, M+1):
                am2 = a[p-2] if p >= 2 else mpc(0)
                am1 = a[p-1] if p >= 1 else mpc(0)
                bm2 = b[p-2] if p >= 2 else mpc(0)
                bm1 = b[p-1] if p >= 1 else mpc(0)
                r2A_p = am2 + 2*r0*am1 + r0**2*a[p]
                r2B_p = bm2 + 2*r0*bm1 + r0**2*b[p]
                # Op series j=0..p (a[p+2]/b[p+2] read as the still-zero
                # placeholder from the pre-allocated array -- exactly the
                # "still-unknown=0" convention the scalar path uses too).
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
    def _ev_mp(c, x, d=0):
        """Horner evaluation of Σ c[m]·x^m and its derivatives (d=0 or 1)."""
        N = len(c) - 1
        if d == 0:
            v = c[N]
            for m in range(N-1, -1, -1):
                v = v*x + c[m]
            return v
        else:  # d==1
            v = mpc(0)
            for m in range(N, 0, -1):
                v = v*x + m*c[m]
            return v

    def _Lmat_mp(self, ze, sols):
        nu, r0 = self.nu, self.r0b
        hp = mp.pi/2
        L = matrix(4, 4)
        for q, (fa, fb) in enumerate(sols):
            for row, x in enumerate((-hp, hp)):
                rr = x + r0
                F = self._ev_mp(fa, x, 0); Fp = self._ev_mp(fa, x, 1)
                G = self._ev_mp(fb, x, 0); Gp = self._ev_mp(fb, x, 1)
                # _coef_eval is a no-op (returns nu unchanged) unless the
                # material is radially_graded (Paper B) -- same convention as
                # OutOfPlaneSolver._Lmat_mp's nu_x/T_x/R_x above.
                nu_x = _coef_eval(nu, x)
                L[row, q] = Fp + nu_x*(F - ze*G)/rr
                L[row+2, q] = Gp + (ze*F - G)/rr
        return L

    def _det_mp(self, ze, Om):
        L = self._Lmat_mp(ze, self._series_mp(ze, Om))
        for i in range(4):
            s = max(abs(L[i, j]) for j in range(4))
            if s > 0:
                for j in range(4): L[i, j] /= s
        return mp.det(L)

    def _refine_root_mp(self, z, Om, iters=4):
        z = mpc(z)
        for _ in range(iters):
            f = self._det_mp(z, Om)
            h = mpf(10)**(-(mp.dps//2))*max(mpf(1), abs(z))
            fp = (self._det_mp(z+h, Om) - self._det_mp(z-h, Om))/(2*h)
            if fp == 0: break
            dz = f/fp
            z -= dz
            if abs(dz) < mpf(10)**(-(mp.dps-8)): break
        return z

    def _amp_mp(self, ze, Om, sols):
        L = self._Lmat_mp(ze, sols)
        cn = [sum(abs(L[r, q])**2 for r in range(4)) for q in range(4)]
        pc = max(range(4), key=lambda q: cn[q])
        oc = [q for q in range(4) if q != pc]
        rn = [sum(abs(L[r, q])**2 for q in range(4)) for r in range(4)]
        ur = sorted(range(4), key=lambda r: rn[r], reverse=True)[:3]
        sub = matrix(3, 3); rhs = matrix(3, 1)
        for i, row in enumerate(ur):
            rhs[i, 0] = -L[row, pc]
            for j, col in enumerate(oc):
                sub[i, j] = L[row, col]
        try:
            t = mp.lu_solve(sub, rhs)
            A = [mpc(0)]*4; A[pc] = mpc(1)
            for j, col in enumerate(oc): A[col] = t[j, 0]
        except Exception:
            A = [mpc(1), mpc(0), mpc(0), mpc(0)]
        s = mp.sqrt(sum(abs(x)**2 for x in A))
        return [x/s for x in A]

    def _F_mp(self, x, sols, A, d=0):
        return sum(A[q]*self._ev_mp(sols[q][0], x, d) for q in range(4))
    def _G_mp(self, x, sols, A, d=0):
        return sum(A[q]*self._ev_mp(sols[q][1], x, d) for q in range(4))

    def _build_K_real(self, Om_f, brs, lagrange=True, fast_scan=False):
        """Real augmented K̂ (Eq.53 + Lagrange row Eq.52), corrected pairing (B7).
        u_r = F sin(ζy+φ), u_θ = G cos(ζy+φ).
        fast_scan=True skips mp root refinement for 8× scan speedup.
        """
        Om = mpf(Om_f)
        c11, nu, R, r0 = self.c11, self.nu, self.R, self.r0b
        Th = self.Theta

        # Build the per-root cache directly from `brs` (not via
        # real_basis_items) so we can determine, for each pure-imaginary
        # root, which projection pair gives nonzero basis functions.
        cache = {}
        for z in set(brs):
            if fast_scan:
                zr = mpc(z)
            else:
                zr = self._refine_root_mp(mpc(z), Om)
            sols = self._series_mp(zr, Om)
            A = self._amp_mp(zr, Om, sols)
            cache[z] = (zr, sols, A)

        # BUGFIX (Research18): for a pure-imaginary zeta=ia root, F and G are
        # each purely real OR purely imaginary (never mixed) — but WHICH one
        # is real flips with Omega_bar (a branch-crossing phenomenon). The
        # projection pair must therefore be chosen per-root, per-frequency:
        # whichever of F,G is real at that root determines the ('re','im')
        # vs ('im','re') pairing that yields nonzero basis contributions.
        # A static choice (Research17/early Research18) gives two
        # structurally-zero rows/columns whenever the "wrong" parity occurs
        # — the root cause of the Part2 0.25pi mode3 '-inf plateau'.
        imag_proj = {}
        for z, (zr, sols, A) in cache.items():
            if abs(z.real) < 1e-9 and abs(z.imag) > 1e-9:
                x0 = self.nodes[0]
                F0 = self._F_mp(x0, sols, A, 0)
                # If F is (numerically) real at this root, q=0 with P='re'
                # and q=1 with P='im' give the nonzero pair; if F is
                # imaginary (so G is real), the swap applies.
                f_is_real = abs(mp.im(F0)) < abs(mp.re(F0)) * 1e-6 + mpf('1e-30')
                imag_proj[z] = ('im', 're') if f_is_real else ('re', 'im')

        items = []
        for z in brs:
            if abs(z.imag) < 1e-9:
                items += [(z, 0, 'asis'), (z, 1, 'asis')]
            elif abs(z.real) < 1e-9:
                p0, p1 = imag_proj[z]
                items += [(z, 0, p0), (z, 1, p1)]
            else:
                items += [(z, 0, 're'), (z, 0, 'im'), (z, 1, 're'), (z, 1, 'im')]
        n = len(items)

        # Pre-compute F, Fp, G, Gp at all quadrature nodes for each distinct root.
        node_FG = {}  # z -> list of (F, Fp, G, Gp) per node
        for z in {it[0] for it in items}:
            ze, sols, A = cache[z]
            node_FG[z] = []
            for x in self.nodes:
                F  = self._F_mp(x, sols, A, 0); Fp = self._F_mp(x, sols, A, 1)
                G  = self._G_mp(x, sols, A, 0); Gp = self._G_mp(x, sols, A, 1)
                node_FG[z].append((F, Fp, G, Gp))

        # ---- BC-driven theta-edge assembly (Phase-0 Step 3) -------------------
        # FREE edge: integral(trial-stress . test-disp); CLAMPED edge:
        # integral(trial-disp . test-stress) plus a weak Lagrange constraint on
        # the edge-normal displacement (u_r = F = 0).  Cantilever bc =
        # [FREE@+Theta, CLAMPED@-Theta] reproduces the original assembly exactly.
        edges = self.bc.edges
        ne = len(edges)
        edge_free = [e.kind == EdgeKind.FREE for e in edges]
        clamp_eidx = [ei for ei in range(ne) if not edge_free[ei]]
        # Edge-orientation factor (2026-07-08 free-free fix): the boundary
        # functional arises from [.]_{-Theta}^{+Theta}, so a FREE edge's
        # transplanted +Theta functional must flip sign at -Theta. CLAMPED
        # edges keep +1 (their coded form was derived at -Theta; see B8).
        # Cantilever (free edge at +Theta) is byte-identical under this.
        _orient = [ (e.sign if edge_free[_k] else 1)
                    for _k, e in enumerate(edges) ]

        # Per item, per edge: projected stress (p_tyy, p_tyr) and displacement
        # (p_G, p_F) node vectors using THIS edge's angular sign.
        edata = []
        for (z, q, P) in items:
            ze, sols, A = cache[z]
            ph = q*mp.pi/2
            per_edge = []
            for e in edges:
                sg = e.sign
                ssg = mp.sin(sg*ze*Th + ph); csg = mp.cos(sg*ze*Th + ph)
                ptyy=[]; ptyr=[]; pG=[]; pF=[]
                for k, x in enumerate(self.nodes):
                    rr = x + r0
                    F, Fp, G, Gp = node_FG[z][k]
                    # _coef_eval is a no-op unless the material is
                    # radially_graded (Paper B); pointwise substitution is
                    # exact here (theta-edge terms differentiate only w.r.t.
                    # theta, never r) -- same convention as _Lmat_mp above and
                    # OutOfPlaneSolver's own Tyy/Qy/Mc evaluation.
                    c11_x = _coef_eval(c11, x); nu_x = _coef_eval(nu, x)
                    R_x = _coef_eval(R, x)
                    tyy = c11_x*(nu_x*Fp + R_x*(F - ze*G)/rr)
                    tyr = Gp - G/rr + ze*F/rr
                    ptyy.append(mp_proj(tyy*ssg, P))   # stress Tyy
                    ptyr.append(mp_proj(tyr*csg, P))   # stress Tyr
                    pG.append(mp_proj(G*csg, P))       # disp   G
                    pF.append(mp_proj(F*ssg, P))       # disp   F
                per_edge.append((ptyy, ptyr, pG, pF))
            edata.append(per_edge)

        n_clamp = len(clamp_eidx) if lagrange else 0
        size = n + n_clamp
        K = matrix(size, size)
        for i in range(n):
            for j in range(n):
                s = mpf(0)
                for k, wt in enumerate(self.wts):
                    acc = mpf(0)
                    for ei in range(ne):
                        ptyy_j, ptyr_j, pG_j, pF_j = edata[j][ei]
                        ptyy_i, ptyr_i, pG_i, pF_i = edata[i][ei]
                        if edge_free[ei]:
                            acc += _orient[ei]*(ptyy_j[k]*pG_i[k] + ptyr_j[k]*pF_i[k])
                        else:
                            acc += -(pF_j[k]*ptyr_i[k] + pG_j[k]*ptyy_i[k])
                    s += wt*acc
                K[i, j] = s
        # One Lagrange row/col per clamped edge: weakly impose u_r=F=0 there.
        if lagrange:
            for row, ei in enumerate(clamp_eidx):
                rr_idx = n + row
                for j in range(n):
                    _, _, _, pF_j = edata[j][ei]
                    s = mpf(0)
                    for k, wt in enumerate(self.wts):
                        s += wt*pF_j[k]
                    K[rr_idx, j] = s; K[j, rr_idx] = s
                for row2, ei2 in enumerate(clamp_eidx):
                    K[rr_idx, n + row2] = mpf(0)

        return K, size

    def logdet(self, Om_f, brs, lagrange=True, fast_scan=False):
        K, size = self._build_K_real(Om_f, brs, lagrange, fast_scan=fast_scan)
        return equilibrated_logdet(K, size)

    def sigma_min(self, Om_f, brs, lagrange=True, fast_scan=False):
        """log10 of the smallest singular value of the (constrained) in-plane K.
        σ_min → 0 exactly at an in-plane natural frequency."""
        K, size = self._build_K_real(Om_f, brs, lagrange, fast_scan=fast_scan)
        return sigma_min_from_K(K, size)

    def _golden_sigmin(self, a, b, n_dofs, ze_max, iters,
                       n_quad_iter=None, seed_brs=None):
        """Golden-section minimise log10 σ_min of the UNCONSTRAINED in-plane K
        using a FRESH branch set at every evaluation.

        Returns (Ω*, log10 σ_min_unconstrained, log10 σ_min_constrained).

        PAPER PHYSICS (corrected this pass): the in-plane natural frequencies
        are the DOUBLE ROOTS of the unconstrained variational equation Eq.
        (48) — that is the real mode condition, so the UNCONSTRAINED matrix is
        the primary signal here.  The Lagrange-constrained matrix K# (Eq.
        52-53) is only an auxiliary the paper uses to split those double roots
        for amplitude determination; about half its roots are spurious.
        Verified on (1.25, 0.25π): the unconstrained σ_min dips deeply at all
        three paper modes (-3.9/-4.0/-4.4) and is shallow at the spurious
        constrained root Ω=0.6445 (-1.8) — the exact opposite of the
        constrained σ_min, which is why a constrained-primary detector both
        missed real modes and accepted spurious ones.  The constrained value
        is still returned so find_modes_sigmin can, if desired, note where the
        two coincide (true mode) vs. not.
        """
        g = (np.sqrt(5) - 1) / 2

        def f(Om):
            raw = full_search(self.fast, Om, xmax=ze_max)
            sel, cnt = select_fill(raw, n_dofs)
            return float(self.sigma_min(Om, sel, lagrange=False, fast_scan=False))

        x1 = b - g * (b - a);  x2 = a + g * (b - a)
        f1, f2 = f(x1), f(x2)
        for _ in range(iters):
            if f1 < f2:
                b, x2, f2 = x2, x1, f1
                x1 = b - g * (b - a);  f1 = f(x1)
            else:
                a, x1, f1 = x1, x2, f2
                x2 = a + g * (b - a);  f2 = f(x2)
        Om_star = (a + b) / 2

        # Single fresh basis at Ω*, used for BOTH σ_min values (the constrained
        # cross-check is then free — no extra full_search).  SIGMIN_P2_CROSSCHECK
        # (default 1) computes the constrained σ_min for the informational
        # coincidence note; set 0 to skip it entirely.
        raw_s = full_search(self.fast, Om_star, xmax=ze_max)
        sel_s, _ = select_fill(raw_s, n_dofs)
        s_star = float(self.sigma_min(Om_star, sel_s, lagrange=False,
                                      fast_scan=False))   # unconstrained (primary)
        s_star = _degen_guard(sel_s, s_star)              # branch-collision guard
        if os.environ.get("SIGMIN_P2_CROSSCHECK", "1") == "1":
            s_constr = float(self.sigma_min(Om_star, sel_s, lagrange=True,
                                            fast_scan=False))
        else:
            s_constr = float('nan')

        return (Om_star, s_star, s_constr)

    # ---------- helper: pack solver parameters for worker dispatch -----------
    def _p2_worker_params(self):
        """Return (g_scalars, m_scalars) tuples for worker reconstruction."""
        g = (float(self.geom.r0_bar), float(self.geom.Theta),
             float(self.geom.b),      float(self.geom.r_0),
             float(self.geom.R_i),    float(self.geom.R_o),
             float(self.geom.h))
        # FIX (2026-07-26): pack the IP coupling constant self.nu (=
        # ip_constants()[2] = mu_theta = D12/D11 for orthotropy), NOT
        # mat.nu_bar -- the exact same "second definition of what a
        # material is" mistake the OOP path's _p1_worker_params already
        # avoids (see its comment, three lines above this class in the
        # file). Equal for isotropic (byte-identical, since nu_bar==mu_theta
        # there), so this is a no-op for every already-validated isotropic
        # table. For orthotropy this was a live bug: _ScalarMat inherits the
        # base MaterialModel.ip_constants(), which returns nu_bar verbatim,
        # so every worker-pool IP search (find_modes_sigmin/
        # find_natural_frequencies, i.e. every real cluster run) silently
        # used mu_r instead of mu_theta regardless of what
        # OrthotropicMaterial.ip_constants() itself returns or how it's
        # monkeypatched in-process -- workers never see the original mat
        # object, only this scalar tuple. Discovered when landing the
        # OrthotropicMaterial.ip_constants() mu_theta fix (LESSONS_LEARNED
        # Sec 24.2/35): the prior Wang Table 4 cluster validation (job
        # 2332092) patched mat.ip_constants in the driver process only, so
        # it actually re-tested the OLD nu_bar behavior end-to-end without
        # anyone noticing -- that "8/8 PASS" result does not yet validate
        # mu_theta and needs to be re-run now that this is fixed.
        m = (float(self.nu),              float(self.mat.T),
             float(self.mat.R),           float(self.mat.c11_bar),
             float(self.mat.c11_eff_bar), float(self.mat.c66),
             float(self.mat.rho))
        return g, m

    def find_natural_frequencies(self, Omega_range=(0.02, 0.30), n_scan=57,
                                  n_dofs=14, ze_max=14.0, dip_orders=1.0,
                                  polish_iters=8, n_quad_iter=None,
                                  n_modes_wanted=3,
                                  verbose=True, n_workers=None):
        """Find Part 2 natural frequencies with parallel scan and polish.

        FULLY GENERAL (Research22+): no paper-table reference values are
        used anywhere in this method.  Candidate acceptance relies only on
        the physics-based genuine-mode filter (Phase 4: unconstrained-det
        local dip or deep constrained/unconstrained zero), and mode
        selection takes the lowest n_modes_wanted genuine candidates.  A
        Phase 5 gap-refinement pass (mirroring Part 1) rescans gaps between
        accepted modes if fewer than n_modes_wanted genuine modes are found.

        FIX (E1): Dual-matrix detection (constrained + unconstrained) retained.

        PARALLELISM (Research8):
          Phase 1 — sequential float64 branch tracking (fast, ~1 s).
          Phase 2 — all mp logdet evaluations (both constrained and
                    unconstrained) dispatched to the worker pool in one batch.
          Phase 3 — all polish candidates submitted to the pool simultaneously.
          Phase 4 — genuine-mode filter evaluations (3 points per candidate)
                    also run in parallel across all candidates.

        n_workers : override N_WORKERS for this call (None → use global).
        """
        # Step 6: thread this solver's boundary condition to every worker pool
        # via the inherited R40_BC_KIND env var (workers rebuild it with
        # _worker_bc()).  Covers ALL pools, including the shared scan/refine
        # helpers; defaults to clamped_free so cantilever runs are unchanged.
        if self.bc.token not in _BC_REGISTRY_TOKENS:
            raise NotImplementedError(
                "parallel find_natural_frequencies has no registered worker BC "
                "for boundary=%r (token=%r)" % (self.bc, self.bc.token))
        os.environ['R40_BC_KIND'] = self.bc.token
        nw = n_workers if n_workers is not None else N_WORKERS
        t0 = time.time()
        print(f"\n{'='*68}")
        print(f"[Part 2 – In-Plane]  Ω̄ ∈ [{Omega_range[0]:.4f},{Omega_range[1]:.4f}]")
        print(f"  n_dofs={n_dofs}, M={self.M}, dps={mp.dps}")
        print(f"  N_WORKERS={nw}  n_modes_wanted={n_modes_wanted}")
        print(f"{'='*68}")

        Om_lo, Om_hi = Omega_range
        Oms = np.linspace(Om_lo, Om_hi, n_scan)

        # ── Phase 1: sequential float64 branch tracking ───────────────────────
        print("  Phase 1: branch tracking (sequential, fast) …", flush=True)
        prev = full_search(self.fast, float(Oms[0]), xmax=ze_max)
        branch_sets = []
        signatures  = []
        for Om in Oms:
            prev = track(self.fast, float(Om), prev)
            sel_t, cnt = select_fill(prev, n_dofs)
            if cnt < n_dofs:
                prev = full_search(self.fast, float(Om), xmax=ze_max)
                sel_t, cnt = select_fill(prev, n_dofs)
            branch_sets.append(list(sel_t))
            signatures.append(set_signature(sel_t))

        # Phase 2: parallel mp logdet evaluations
        g_sc, m_sc = self._p2_worker_params()
        dps = mp.dps
        n_quad_scan = max(15, len(self.nodes) // 2)  # 15 nodes for scan
        scan_args = [
            (float(Om), list(brs), g_sc, m_sc, self.M, n_quad_scan,
             ze_max, n_dofs, dps)
            for Om, brs in zip(Oms, branch_sets)
        ]

        print(f"  Phase 2: mp logdet sweep ({n_scan} points, "
              f"{nw} workers) …", flush=True)

        # All three phases share one pool to avoid the overhead and memory
        # spike of repeated process spawning within a single geometry run.
        rows = []
        scan_results = []
        t_scan2 = time.time()
        with concurrent.futures.ProcessPoolExecutor(max_workers=nw) as pool:
            n_done2 = 0
            for result in _throttled_map(pool, _p2_logdet_worker, scan_args):
                scan_results.append(result)
                n_done2 += 1
                if not verbose and n_done2 % max(1, n_scan // 10) == 0:
                    pct = 100 * n_done2 // n_scan
                    elapsed2 = time.time() - t_scan2
                    rate2 = n_done2 / elapsed2 if elapsed2 > 0 else 0
                    eta2 = (n_scan - n_done2) / rate2 if rate2 > 0 else 0
                    print(f"  scan {n_done2}/{n_scan} ({pct}%)  "
                          f"elapsed {elapsed2:.0f}s  ETA {eta2:.0f}s", flush=True)
            scan_results.sort(key=lambda r: r[0])
            print(f"  Phase 2 done: {n_scan} pts in {time.time()-t_scan2:.1f}s", flush=True)

            for idx, ((Om_f, sgn_c, ld_c, ld_u), sig) in enumerate(
                    zip(scan_results, signatures)):
                rows.append((Om_f, sgn_c, ld_c, ld_u, sig))
                if verbose:
                    br = ",".join(f"{z:.2f}" for z in branch_sets[idx])
                    print(f"  {Om_f:.5f}  sign={sgn_c:+d}  "
                          f"logC={ld_c:+8.3f}  logU={ld_u:+8.3f}  [{br}]",
                          flush=True)

            # ── Phase 2b: collect candidates ───────────────────────────────────
            cands = []
            seen  = set()

            # Sign changes in the constrained matrix.
            for k in range(1, len(rows)):
                Om_a, Om_b = rows[k-1][0], rows[k][0]
                key = (round(Om_a, 4), round(Om_b, 4))
                if rows[k-1][1] * rows[k][1] < 0:
                    if abs(rows[k-1][2] - rows[k][2]) < 4.0 and key not in seen:
                        cands.append((Om_a, Om_b, 'sign_C'))
                        seen.add(key)

            # Local minima in the unconstrained matrix.  This is the primary
            # detector for Part 2: genuine modes show a wide, reliable dip in
            # the unconstrained det that is often missed as a sign change in the
            # constrained det (the sign-change window can be as narrow as 0.002).
            # Use the loose threshold always — the Phase 4 physics-based
            # genuine-mode filter removes false positives downstream.
            lds_u = [r[3] for r in rows]
            med_u = float(np.median(lds_u))
            dip_thresh_u = min(0.7, dip_orders)
            for k in range(1, len(rows) - 1):
                Om_a, Om_b = rows[k-1][0], rows[k+1][0]
                key = (round(Om_a, 4), round(Om_b, 4))
                if (lds_u[k] < lds_u[k-1] and lds_u[k] < lds_u[k+1]
                        and lds_u[k] < med_u - dip_thresh_u
                        and key not in seen):
                    cands.append((Om_a, Om_b, 'dip_U'))
                    seen.add(key)

            # Local minima in the constrained matrix (backup detector).
            # Requires a stable branch signature across the three points to
            # avoid flagging topology-change artefacts as candidates.
            lds_c = [r[2] for r in rows]
            med_c = float(np.median(lds_c))
            dip_thresh_c = min(0.7, dip_orders)
            for k in range(1, len(rows) - 1):
                Om_a, Om_b = rows[k-1][0], rows[k+1][0]
                key = (round(Om_a, 4), round(Om_b, 4))
                if (lds_c[k] < lds_c[k-1] and lds_c[k] < lds_c[k+1]
                        and rows[k-1][4] == rows[k][4] == rows[k+1][4]
                        and lds_c[k] < med_c - dip_thresh_c
                        and key not in seen):
                    cands.append((Om_a, Om_b, 'dip_C'))
                    seen.add(key)

            print(f"\n  candidates: "
                  f"{[(f'{a:.4f}', f'{b:.4f}', t) for a,b,t in cands]}")

            # ── Phase 2c: cheap pre-filter by candidate depth ───────────────
            # Cap polished candidates to the MAX_POLISH deepest (by the log|det| already
            # measured at Phase-2 brackets) plus any past DEPTH_AUTOKEEP.  No extra evals.
            ldu_lookup = {round(r[0], 5): r[3] for r in rows}

            def _cand_depth_p2(a, b):
                la = ldu_lookup.get(round(a, 5))
                lb = ldu_lookup.get(round(b, 5))
                vals = [v for v in (la, lb) if v is not None and np.isfinite(v)]
                return min(vals) if vals else 0.0

            DEPTH_AUTOKEEP = -4.0
            MAX_POLISH = max(2 * n_modes_wanted + 4, 12)
            cand_depths = [_cand_depth_p2(a, b) for a, b, _ in cands]
            order = sorted(range(len(cands)), key=lambda i: cand_depths[i])
            keep_idx = set()
            for rank, i in enumerate(order):
                if cand_depths[i] < DEPTH_AUTOKEEP or rank < MAX_POLISH:
                    keep_idx.add(i)
            if len(keep_idx) < len(cands):
                print(f"  Phase 2c: pre-filter kept {len(keep_idx)}/{len(cands)} "
                      f"candidates (deepest logU, cap={MAX_POLISH})", flush=True)
            cands = [cands[i] for i in sorted(keep_idx)]

            # ── Phase 3: polishing (same pool) ────────────────────────────────
            nq_fast = (n_quad_iter if n_quad_iter is not None
                       else max(8, len(self.nodes) // 2))
            worker_fn_map = {
                'sign_C': _p2_bisect_worker,
                'dip_U':  _p2_golden_worker,
                'dip_C':  _p2_golden_worker,
            }
            cand_types = [t for _, _, t in cands]
            cand_lo    = [a for a, _, _ in cands]
            cand_hi    = [b for _, b, _ in cands]

            coarse_Oms_p2 = list(Oms)
            def _nearest_brs_p2(Om_mid):
                idx = int(np.argmin([abs(Om_mid - c) for c in coarse_Oms_p2]))
                return branch_sets[idx]

            p_args = [
                (a, b, g_sc, m_sc, self.M, len(self.nodes),
                 ze_max, n_dofs, polish_iters, nq_fast, dps,
                 _nearest_brs_p2((a + b) / 2))
                for a, b in zip(cand_lo, cand_hi)
            ]

            print(f"  Phase 3: polishing {len(cands)} candidate(s) …", flush=True)
            polish_results = []
            t_p2_polish = time.time()
            if p_args:
                fns = [worker_fn_map[t] for t in cand_types]
                futures = [pool.submit(fn, arg) for fn, arg in zip(fns, p_args)]
                for i, fut in enumerate(futures):
                    try:
                        Om_star, ld_star = fut.result()
                        polish_results.append(
                            (cand_types[i], cand_lo[i], cand_hi[i], Om_star, ld_star))
                        print(f"  polish {i+1}/{len(cands)}  "
                              f"Ω̄={Om_star:.6f}  logC={ld_star:+.2f}  "
                              f"({time.time()-t_p2_polish:.0f}s)", flush=True)
                    except Exception as exc:
                        print(f"  WARNING: polish [{cand_lo[i]:.4f},"
                              f"{cand_hi[i]:.4f}] failed: {exc}")

            # ── Phase 4: genuine-mode filter (same pool) ───────────────────────
            # Evaluate the unconstrained log|det| at three points around each
            # polished candidate to confirm it is a genuine local minimum, not
            # a topology artefact.  Only 3 * len(polish_results) evaluations.
            print("  Phase 4: genuine-mode filter …", flush=True)
            d = 5e-3
            filter_points = [
                Om_star + delta
                for _, _, _, Om_star, _ in polish_results
                for delta in (-d, 0.0, d)
            ]
            filter_args = [
                (Om_pt, g_sc, m_sc, self.M, n_quad_scan,
                 ze_max, n_dofs, False, dps,
                 _nearest_brs_p2(Om_pt))
                for Om_pt in filter_points
            ]
            filter_lds = []
            if filter_args:
                for ld_val in _throttled_map(pool, _p2_ld_at_worker, filter_args):
                    filter_lds.append(ld_val)

        # Pool closed; all worker processes have exited.

        # ── Verify all candidates with the genuine/proximity physics filter ──
        verified = []  # dicts: typ, a, b, Om_star, ld_star(logC), logU, genuine_physics
        for idx, (typ, a, b, Om_star, ld_star) in enumerate(polish_results):
            if filter_lds:
                lo_u  = filter_lds[idx * 3]
                mid_u = filter_lds[idx * 3 + 1]
                hi_u  = filter_lds[idx * 3 + 2]
            else:
                lo_u = mid_u = hi_u = 0.0

            # A candidate is genuine if the unconstrained det has a local dip
            # at Om_star, OR either det is extremely deep (near-zero), OR it is
            # a CONFIRMED sign change (a real bracketed zero) whose logC is at
            # least as deep as the floor.  The last clause is essential: a
            # genuine mode can manifest as a C sign-change with only moderate
            # depth and no U-dip (e.g. real mode 0.347 at 2Θ=0.25π, logC=-4.06),
            # which the dip/deep tests alone wrongly discarded.
            unc_dip = (mid_u < lo_u and mid_u < hi_u)
            deep_c  = (ld_star < -8.0)
            deep_u  = (mid_u < -5.5)
            sign_genuine = (str(typ).startswith('sign') and ld_star < _P2_LD_GENUINE)
            genuine_physics = unc_dip or deep_c or deep_u or sign_genuine

            verified.append({'typ': typ, 'a': a, 'b': b, 'Om_star': Om_star,
                              'ld_star': ld_star, 'logU': mid_u,
                              'genuine_physics': genuine_physics})

        freqs = []

        for v in verified:
            if verbose:
                status = 'GENUINE' if v['genuine_physics'] else 'SPURIOUS'
                logU_str = f"{v['logU']:+.2f}" if v['logU'] is not None else "  n/a"
                print(f"  polished [{v['a']:.4f},{v['b']:.4f}] ({v['typ']}) → "
                      f"Omega_bar={v['Om_star']:.6f}  logC={v['ld_star']:+.2f}  "
                      f"logU={logU_str}  {status}")

        # ── Mode selection ──────────────────────────────────────────────────
        # Collect genuine candidates, dedupe near-duplicates, sort by Ω̄, and
        # keep the lowest n_modes_wanted — purely a frequency-ordering
        # criterion, valid for any geometry.
        def _merge_genuine(verified_list, accepted):
            for v in verified_list:
                if not v['genuine_physics']:
                    continue
                Om_star = v['Om_star']
                if all(abs(Om_star - f['Om_star']) > 1e-3 for f in accepted):
                    accepted.append(v)
            accepted.sort(key=lambda v: v['Om_star'])
            return accepted

        accepted = _merge_genuine(verified, [])

        # ── Phase 5: gap-refinement (replaces reference-anchored Phase 4b) ──
        # If fewer than n_modes_wanted genuine modes were found, rescan the
        # gaps — below the lowest accepted mode, between consecutive
        # accepted modes, and up to Omega_range[1] — with a fresh-seeded,
        # small-step local scan.  Driven purely by mode COUNT, never by a
        # known answer.
        MAX_REFINE_PASSES = 3
        refine_pass = 0
        while len(accepted) < n_modes_wanted and refine_pass < MAX_REFINE_PASSES:
            refine_pass += 1
            gaps = []
            edges = [Om_lo] + [v['Om_star'] for v in accepted] + [Om_hi]
            for k in range(len(edges) - 1):
                lo_g, hi_g = edges[k], edges[k+1]
                if hi_g - lo_g > 1e-6:
                    gaps.append((lo_g, hi_g))
            if not gaps:
                break
            print(f"  Phase 5: gap-refinement pass {refine_pass} — "
                  f"{len(accepted)}/{n_modes_wanted} genuine modes found, "
                  f"rescanning {len(gaps)} gap(s) "
                  f"{[(f'{a:.4f}', f'{b:.4f}') for a,b in gaps]} …", flush=True)
            new_verified = []
            for (lo_g, hi_g) in gaps:
                new_verified.extend(
                    self._local_rescan(lo_g, hi_g, n_dofs, ze_max, polish_iters,
                                        nq_fast, n_quad_scan, dps, g_sc, m_sc,
                                        nw, verbose, label=f"gap[{lo_g:.4f},{hi_g:.4f}]"))
            if not new_verified:
                break
            n_before = len(accepted)
            accepted = _merge_genuine(new_verified, accepted)
            if len(accepted) == n_before:
                break

        for v in accepted:
            freqs.append(v['Om_star'])

        freqs.sort()
        freqs = freqs[:n_modes_wanted]
        print(f"  elapsed {time.time()-t0:.0f}s")
        return freqs

    # ------------------------------------------------------------------
    # Local rescan of a sub-window [lo, hi]: fresh-seeded sequential
    # tracking with small steps (Fix G/H), candidate detection (sign change
    # in logC, or deepest local minimum of logU), polish, and the same
    # genuine-mode physics check (unconstrained dip / deep zero) as the
    # main scan.  Used by Phase 5 gap-refinement.  Returns a list of
    # verified dicts (with 'genuine_physics' set).
    # ------------------------------------------------------------------
    def _local_rescan(self, lo, hi, n_dofs, ze_max, polish_iters, nq_fast,
                       n_quad_scan, dps, g_sc, m_sc, nw, verbose, label=""):
        n_pts = max(15, min(30, MAX_SCAN_PTS // 2))
        t_Oms = np.linspace(lo, hi, n_pts)
        prev_seed = full_search(self.fast, float(t_Oms[0]), xmax=ze_max)
        t_rows, t_brs = [], []  # rows: (Om, sign_c, logC, logU)
        for Om_t in t_Oms:
            prev_seed = track(self.fast, float(Om_t), prev_seed)
            sel_t, cnt_t = select_fill(prev_seed, n_dofs)
            if cnt_t < n_dofs:
                prev_seed = full_search(self.fast, float(Om_t), xmax=ze_max)
                sel_t, _ = select_fill(prev_seed, n_dofs)
            sgn_c, ld_c = self._logdet_nq(float(Om_t), sel_t, lagrange=True,
                                           n_quad_override=n_quad_scan, fast_scan=False)
            _, ld_u = self._logdet_nq(float(Om_t), sel_t, lagrange=False,
                                       n_quad_override=n_quad_scan, fast_scan=False)
            t_rows.append((float(Om_t), float(sgn_c), float(ld_c), float(ld_u)))
            t_brs.append(list(sel_t))
            if verbose:
                print(f"    [{label}] {Om_t:.5f}  sign={int(sgn_c):+d}  "
                      f"logC={ld_c:+8.3f}  logU={ld_u:+8.3f}", flush=True)

        cands = []  # (a, b, typ, k)
        for k in range(1, len(t_rows)):
            if t_rows[k-1][1] * t_rows[k][1] < 0:
                cands.append((t_rows[k-1][0], t_rows[k][0], 'sign_C', k))

        lds_u_t = [r[3] for r in t_rows]
        med_u_t = float(np.median(lds_u_t))
        for k in range(1, len(t_rows) - 1):
            if (lds_u_t[k] < lds_u_t[k-1] and lds_u_t[k] < lds_u_t[k+1]
                    and lds_u_t[k] < med_u_t - 0.7):
                cands.append((t_rows[k-1][0], t_rows[k+1][0], 'dip_U', k))

        verified = []
        for (a_t, b_t, typ, k) in cands:
            seed_mid = t_brs[k]
            if typ == 'sign_C':
                Om_star_t, ld_star_t = self._bisect(
                    a_t, b_t, n_dofs, ze_max, polish_iters,
                    n_quad_iter=nq_fast, seed_brs=seed_mid)
            else:
                Om_star_t, ld_star_t = self._golden(
                    a_t, b_t, n_dofs, ze_max, polish_iters,
                    n_quad_iter=nq_fast, seed_brs=seed_mid)

            # 3-point genuine_physics check using the unconstrained det.
            d_probe = max(1e-3, 0.5 * (hi - lo) / n_pts)
            logU_vals = {}
            for delta in (-d_probe, 0.0, d_probe):
                Om_p = Om_star_t + delta
                p_seed = track(self.fast, Om_p, [complex(z) for z in seed_mid])
                sel_p, cnt_p = select_fill(p_seed, n_dofs)
                if cnt_p < n_dofs:
                    p_seed = full_search(self.fast, Om_p, xmax=ze_max)
                    sel_p, _ = select_fill(p_seed, n_dofs)
                _, ld_u_p = self.logdet(Om_p, sel_p, lagrange=False, fast_scan=False)
                logU_vals[delta] = float(ld_u_p)
            mid_u = logU_vals[0.0]
            unc_dip = (mid_u < logU_vals[-d_probe] and mid_u < logU_vals[d_probe])
            deep_c  = (ld_star_t < -8.0)
            deep_u  = (mid_u < -5.5)
            sign_genuine = (str(typ).startswith('sign') and ld_star_t < _P2_LD_GENUINE)
            genuine = bool(unc_dip or deep_c or deep_u or sign_genuine)

            verified.append({'typ': f'local_{typ}', 'a': a_t, 'b': b_t,
                              'Om_star': Om_star_t, 'ld_star': ld_star_t,
                              'logU': mid_u, 'genuine_physics': genuine})
            if verbose:
                status = 'GENUINE' if genuine else 'SPURIOUS'
                print(f"  local [{label}] [{a_t:.4f},{b_t:.4f}] ({typ}) → "
                      f"Omega_bar={Om_star_t:.6f}  logC={ld_star_t:+.2f}  "
                      f"logU={mid_u:+.2f}  {status}", flush=True)
        return verified

    def _ld_at(self, Om, n_dofs, ze_max, lagrange=True, seed_brs=None):
        """Evaluate log|det| at Om. Uses track from seed_brs or full_search fallback."""
        if seed_brs:
            raw = track(self.fast, Om, [complex(z) for z in seed_brs])
            if not raw:
                raw = full_search(self.fast, Om, xmax=ze_max)
        else:
            raw = full_search(self.fast, Om, xmax=ze_max)
        sel, _ = select_fill(raw, n_dofs)
        return float(self.logdet(Om, sel, lagrange, fast_scan=False)[1])

    # ------------------------------------------------------------------
    # Fast-quadrature logdet helper for cheap polishing iterations
    # ------------------------------------------------------------------
    def _logdet_nq(self, Om_f, brs, lagrange=True, n_quad_override=None, fast_scan=False):
        """logdet with optional quadrature-node override."""
        if n_quad_override is None or n_quad_override == len(self.nodes):
            return self.logdet(Om_f, brs, lagrange, fast_scan=fast_scan)
        nd, wt = mp_gauss_legendre(n_quad_override)
        saved_nodes, saved_wts = self.nodes, self.wts
        self.nodes = [x * mp.pi / 2 for x in nd]
        self.wts   = [w * mp.pi / 2 for w in wt]
        try:
            result = self.logdet(Om_f, brs, lagrange, fast_scan=fast_scan)
        finally:
            self.nodes, self.wts = saved_nodes, saved_wts
        return result

    def _bisect(self, a, b, n_dofs, ze_max, iters, n_quad_iter=None, seed_brs=None):
        """Bisect sign change in constrained K̂ determinant.
        seed_brs : coarse branch list for Newton seed (avoids full_search).
        Iterations use fast_scan=True; final eval uses full precision.
        """
        n_quad_fast = n_quad_iter if n_quad_iter is not None else max(8, len(self.nodes) // 2)

        if seed_brs:
            prev_a = track(self.fast, a, [complex(z) for z in seed_brs]) or \
                     full_search(self.fast, a, xmax=ze_max)
            prev_b = track(self.fast, b, [complex(z) for z in seed_brs]) or \
                     full_search(self.fast, b, xmax=ze_max)
        else:
            prev_a = full_search(self.fast, a, xmax=ze_max)
            prev_b = full_search(self.fast, b, xmax=ze_max)

        sel_a, cnt = select_fill(prev_a, n_dofs)
        if cnt < n_dofs:
            prev_a = full_search(self.fast, a, xmax=ze_max)
            sel_a, _ = select_fill(prev_a, n_dofs)

        sa, _ = self._logdet_nq(a, sel_a, lagrange=True,
                                n_quad_override=n_quad_fast, fast_scan=False)
        prev_lo, prev_hi = prev_a, prev_b

        for _ in range(iters):
            m = (a + b) / 2
            if abs(m - a) <= abs(m - b):
                prev_m = track(self.fast, m, prev_lo)
            else:
                prev_m = track(self.fast, m, prev_hi)
            sel_m, cnt = select_fill(prev_m, n_dofs)
            if cnt < n_dofs:
                prev_m = full_search(self.fast, m, xmax=ze_max)
                sel_m, _ = select_fill(prev_m, n_dofs)
            sm_, _ = self._logdet_nq(m, sel_m, lagrange=True,
                                     n_quad_override=n_quad_fast, fast_scan=False)
            if sm_ == sa:
                a = m; prev_lo = prev_m
            else:
                b = m; prev_hi = prev_m

        m = (a + b) / 2
        prev_m = track(self.fast, m, prev_lo)
        sel_m, cnt = select_fill(prev_m, n_dofs)
        if cnt < n_dofs:
            prev_m = full_search(self.fast, m, xmax=ze_max)
            sel_m, _ = select_fill(prev_m, n_dofs)
        _, ld_final = self.logdet(m, sel_m, lagrange=True, fast_scan=False)
        return (m, float(ld_final))

    def _golden(self, a, b, n_dofs, ze_max, iters, n_quad_iter=None, seed_brs=None):
        """Golden-section minimise UNCONSTRAINED log|det|.
        seed_brs : coarse branch list for Newton seed.
        Iterations use fast_scan=True; final constrained eval uses full precision.
        """
        n_quad_fast = n_quad_iter if n_quad_iter is not None else max(8, len(self.nodes) // 2)
        g = (np.sqrt(5) - 1) / 2

        if seed_brs:
            prev_x1 = track(self.fast, a, [complex(z) for z in seed_brs]) or \
                       full_search(self.fast, a, xmax=ze_max)
            prev_x2 = track(self.fast, b, [complex(z) for z in seed_brs]) or \
                       full_search(self.fast, b, xmax=ze_max)
        else:
            prev_x1 = full_search(self.fast, a, xmax=ze_max)
            prev_x2 = full_search(self.fast, b, xmax=ze_max)

        def f_tracked(Om, prev_seed):
            prev = track(self.fast, Om, prev_seed)
            sel, cnt = select_fill(prev, n_dofs)
            if cnt < n_dofs:
                prev = full_search(self.fast, Om, xmax=ze_max)
                sel, _ = select_fill(prev, n_dofs)
            ld = float(self._logdet_nq(Om, sel, lagrange=False,
                                       n_quad_override=n_quad_fast,
                                       fast_scan=False)[1])
            return ld, prev

        x1 = b - g * (b - a);  x2 = a + g * (b - a)
        f1, prev_x1 = f_tracked(x1, prev_x1)
        f2, prev_x2 = f_tracked(x2, prev_x2)

        for _ in range(iters):
            if f1 < f2:
                b, x2, f2, prev_x2 = x2, x1, f1, prev_x1
                x1 = b - g * (b - a)
                f1, prev_x1 = f_tracked(x1, prev_x2)
            else:
                a, x1, f1, prev_x1 = x1, x2, f2, prev_x2
                x2 = a + g * (b - a)
                f2, prev_x2 = f_tracked(x2, prev_x1)

        Om_star = (a + b) / 2
        prev_f = track(self.fast, Om_star, prev_x1)
        sel_f, cnt = select_fill(prev_f, n_dofs)
        if cnt < n_dofs:
            prev_f = full_search(self.fast, Om_star, xmax=ze_max)
            sel_f, _ = select_fill(prev_f, n_dofs)
        _, ld_c = self.logdet(Om_star, sel_f, lagrange=True, fast_scan=False)
        return (Om_star, float(ld_c))


# ══════════════════════════════════════════════════════════════════════════════
# ══════════════════════════════════════════════════════════════════════════════
#  SECTION 5b — MODULE-LEVEL PARALLEL WORKER FUNCTIONS
#  Defined at module level so multiprocessing can pickle them (spawn on
#  Win/macOS).  Each worker re-inits mp.dps (per-process state) and rebuilds a
#  lightweight solver from picklable scalars.  See LESSONS_LEARNED.md §Workers.
# ══════════════════════════════════════════════════════════════════════════════

class RectOOPAssembler:
    """mpmath rectangular Eq.53 variational assembler + equilibrated sigma_min
    detector with two-sided conjugate-pair realification.  Uses two
    RectangularCartesianOOP engines (one per symmetry class) for branch/amplitude
    evaluation; the symmetry class is selected per assemble() call."""

    def __init__(self, nu_b=0.3, R=1.0, T=1.0, dps=30, bc='clamped_free'):
        mp.dps = dps
        self.dps = dps
        self.eng_sym = RectangularCartesianOOP(nu_b, R, T, sym=True)
        self.eng_anti = RectangularCartesianOOP(nu_b, R, T, sym=False)
        self.nu = mpf(str(nu_b)); self.R = mpf(str(R)); self.T = mpf(str(T))
        self.mu = 2 * self.T - self.nu
        # BC at the x1=-L end -- 'clamped_free' (default, BYTE-IDENTICAL to
        # every pre-2026-08-02 call site) or 'free_free'. The free_free
        # corner jump is the 2026-09-03 closed-contour checkerboard
        # (assemble() else-branch); sandbox FE-checked against 4 ANSYS
        # anchors, not yet a cluster production table. x2=+/-b stay free
        # in both cases.
        if bc not in ('clamped_free', 'free_free'):
            raise ValueError(f"RectOOPAssembler: unknown bc={bc!r}")
        self.bc = bc

    def eng(self, sym):
        return self.eng_sym if sym else self.eng_anti

    def _pterms(self, e1, e2, H1, H2, phin, m):
        ph = phin + m * mp.pi / 2
        return [(H1 * e1 ** m, e1, ph), (H2 * e2 ** m, e2, ph)]

    def _pint(self, A, B):
        s = mpc(0)
        for ca, ea, pa in A:
            for cb, eb, pb in B:
                s += ca * cb * _rect_J(ea, pa, eb, pb)
        return s

    def assemble(self, full, Lam, sym, lob):
        nu, _R, T, mu = self.nu, self.R, self.T, self.mu
        phin = mp.pi / 2 if sym else mpf(0)
        L = mp.pi * mpf(str(lob)) / 2
        Lf = float(Lam)
        eng = self.eng(sym)
        bd = []
        for xi in full:
            e1, e2 = eng._etas(complex(xi), Lf)
            H1, H2 = eng._amp_ratio(complex(xi), Lf)
            e1 = mpc(e1); e2 = mpc(e2); H1 = mpc(H1); H2 = mpc(H2)
            psi = {m: self._pterms(e1, e2, H1, H2, phin, m) for m in (0, 1, 2, 3)}
            bd.append((mpc(xi), psi))

        def U(xi, r, m, x1):
            return xi ** m * mp.sin(xi * x1 + (r - 1) * mp.pi / 2 + m * mp.pi / 2)

        def peval(terms, x2):
            return sum(c * mp.sin(e * x2 + p) for c, e, p in terms)

        P = len(full)
        K = mp.zeros(2 * P, 2 * P)
        for pi_, (xip, psip) in enumerate(bd):
            for ppi, (xipp, psipp) in enumerate(bd):
                I = {(m, k): self._pint(psip[m], psipp[k])
                     for m in (0, 1, 2, 3) for k in (0, 1, 2, 3)}
                p1t = peval(psip[1], mp.pi / 2); p1b = peval(psip[1], -mp.pi / 2)
                q0t = peval(psipp[0], mp.pi / 2); q0b = peval(psipp[0], -mp.pi / 2)
                for r in (1, 2):
                    for rp in (1, 2):
                        Ut3 = U(xip, r, 3, L); Ut1 = U(xip, r, 1, L)
                        Ut2 = U(xip, r, 2, L); Ut0 = U(xip, r, 0, L)
                        Vt0 = U(xipp, rp, 0, L); Vt1 = U(xipp, rp, 1, L)
                        tip = (Ut3 * Vt0 * I[(0, 0)] + mu * Ut1 * Vt0 * I[(2, 0)]
                               - Ut2 * Vt1 * I[(0, 0)] - nu * Ut0 * Vt1 * I[(2, 0)])
                        Uw0 = U(xip, r, 0, -L); Uw1 = U(xip, r, 1, -L)
                        Vw0 = U(xipp, rp, 0, -L); Vw1 = U(xipp, rp, 1, -L)
                        if self.bc == 'clamped_free':
                            # ---- EXACT pre-2026-08-02 path, byte-for-byte unchanged ----
                            Vw2 = U(xipp, rp, 2, -L)
                            wall = (Uw0 * Vw2 * I[(0, 1)] + T * Uw0 * Vw0 * I[(0, 3)]
                                    - Uw1 * Vw2 * I[(0, 0)] - nu * Uw1 * Vw0 * I[(0, 2)]
                                    - (T - nu) * Uw0 * Vw1 * I[(1, 1)])
                            corner = (T - nu) * (p1t * Uw1 * q0t * Vw0
                                                 + p1b * Uw1 * q0b * Vw0
                                                 - 2 * p1t * Ut1 * q0t * Vt0)
                            K[2 * ppi + (rp - 1), 2 * pi_ + (r - 1)] = tip + wall + corner
                        else:
                            # ---- 'free_free' Kirchhoff corner jump ----
                            # left = orientation-reversed primal at x1=-L (same
                            # formula as tip, L -> -L, minus from reversed n).
                            # corner_ff is the closed-contour Eq.(16) jump at all
                            # four free corners (LESSONS Sec 18.75/18.81-18.83,
                            # 2026-09-03 Grok trace). Naive same-sign 4-term sum
                            # vanishes identically (p1t*q0t + p1b*q0b == 0).
                            # Checkerboard: TB and WT flip because the inward-axis
                            # Jacobian is -1 there (F_in - F_out = -2 t_12).
                            # Under the parity identity this is exactly 2*FIX.
                            # Sandbox FE check (probe_rect_ff_oop_derived_vs_fix_
                            # 2026-09-03.py): 4/4 ANSYS anchors closer than FIX
                            # and CURRENT. clamped_free `corner` above is
                            # untouched. SOLVER_VERSION not bumped: default
                            # bc remains clamped_free; prior free_free scans
                            # (zero corner term) are stale.
                            Uw3 = U(xip, r, 3, -L); Uw2 = U(xip, r, 2, -L)
                            left = -(Uw3 * Vw0 * I[(0, 0)] + mu * Uw1 * Vw0 * I[(2, 0)]
                                     - Uw2 * Vw1 * I[(0, 0)] - nu * Uw0 * Vw1 * I[(2, 0)])
                            corner_ff = (T - nu) * (-2) * (
                                p1t * Ut1 * q0t * Vt0 - p1b * Ut1 * q0b * Vt0
                                - p1t * Uw1 * q0t * Vw0 + p1b * Uw1 * q0b * Vw0)
                            K[2 * ppi + (rp - 1), 2 * pi_ + (r - 1)] = tip + left + corner_ff
        return K

    @staticmethod
    def _conjugate_pair_cols(full, tol=1e-9):
        pairs = []; used = set()
        for p, z in enumerate(full):
            if p in used or abs(z.imag) < tol:
                continue
            for q in range(p + 1, len(full)):
                if q in used:
                    continue
                w = full[q]
                if abs(w.real - z.real) < tol and abs(w.imag + z.imag) < tol:
                    used.add(p); used.add(q)
                    for r in (0, 1):
                        pairs.append((2 * p + r, 2 * q + r))
                    break
        return pairs

    def _realify_conjugate_pairs(self, K, full):
        n = K.cols; I = mpc(0, 1); s2 = mp.sqrt(2)
        for (a, b) in self._conjugate_pair_cols(full):
            for k in range(n):
                xa, xb = K[k, a], K[k, b]
                K[k, a] = (xa + xb) / s2
                K[k, b] = (xa - xb) / (I * s2)
            for k in range(n):
                xa, xb = K[a, k], K[b, k]
                K[a, k] = (xa + xb) / s2
                K[b, k] = (xa - xb) / (I * s2)

    def equil_sigma(self, K, full=None):
        """Column-equilibrated smallest/largest singular-value ratio, with
        two-sided realification of conjugate-pair blocks when `full` is given.
        REQUIRED for complex branches (the Eq.53 form is non-symmetric, so a
        conjugate pair is a [[A,B],[B*,A*]] block made real only by transforming
        BOTH rows and columns)."""
        K = K.copy()
        if full is not None:
            self._realify_conjugate_pairs(K, full)
        n = K.cols
        for j in range(n):
            cn = mp.sqrt(sum(abs(K[i, j]) ** 2 for i in range(n)))
            if cn > 0:
                for i in range(n):
                    K[i, j] /= cn
        S = mp.svd(K, compute_uv=False)
        vals = [S[i] for i in range(S.rows)]
        return float(min(vals) / max(vals))


class RectIPAssembler:
    """mpmath rectangular in-plane Eq.41 UNCONSTRAINED variational assembler +
    equilibrated sigma_min detector with two-sided conjugate-pair realification.
    Frequencies are double roots of det K (Eq.36) found from scratch.

    bc='clamped_free' (default) is the Part 2 cantilever: wall at x1=-L plus
    free tip at x1=+L. bc='free_free' makes both short edges natural (same
    traction bilinear form at +/-L). No Kirchhoff corner jump exists for
    in-plane motion (Part 2 Eq.(1)/(36); 2026-09-05 confirmation). Default
    path is unchanged; SOLVER_VERSION is not bumped."""

    def __init__(self, nu_b=0.3, R=1.0, c11s=None, dps=30, bc='clamped_free'):
        mp.dps = dps
        self.dps = dps
        self.eng_sym = RectangularCartesianIP(nu_b, R, c11s, sym=True)
        self.eng_anti = RectangularCartesianIP(nu_b, R, c11s, sym=False)
        self.nu = mpf(str(nu_b)); self.R = mpf(str(R))
        self.c11s = mpf(str(self.eng_sym.c11s))
        if bc not in ('clamped_free', 'free_free'):
            raise ValueError(f"RectIPAssembler: unknown bc={bc!r}")
        self.bc = bc

    def eng(self, sym):
        return self.eng_sym if sym else self.eng_anti

    def assemble(self, full, Om, sym, lob):
        """2P x 2P unconstrained K, indexed K[test, trial].
        clamped_free: Eq.41 wall+tip. free_free: free-edge form at both x1=+/-L,
        no corner term."""
        nu, _R, c11s = self.nu, self.R, self.c11s
        sph = mp.pi / 2 if sym else mpf(0)
        L = mp.pi * mpf(str(lob)) / 2
        Of = float(Om)
        eng = self.eng(sym)

        # Per branch: x1-wavenumber g (mpc), zeta pair, packed transverse
        # coefficients H1[q],H2[q] (Eq.39), all promoted to mp precision.
        bd = []
        for g in full:
            (z1, z2), H1, H2 = eng.H(complex(g), Of)
            bd.append((mpc(g), (mpc(z1), mpc(z2)),
                       (mpc(H1[0]), mpc(H1[1])), (mpc(H2[0]), mpc(H2[1]))))

        HALF = mp.pi / 2

        # x2-profile terms as lists of (coeff, freq, phase) so products reduce
        # via _rect_J.  cos(.) represented as sin(.+pi/2).
        def prof_u1(zs, H):                 # Psi1 = sum H1[q] sin(z_q x2 + s')
            return [(H[0], zs[0], sph), (H[1], zs[1], sph)]

        def prof_u2(zs, H):                 # Psi2 = sum H2[q] cos(z_q x2 + s')
            return [(H[0], zs[0], sph + HALF), (H[1], zs[1], sph + HALF)]

        def prof_du1(zs, H):                # dPsi1/dx2 = sum H1[q] z_q cos(...)
            return [(H[0] * zs[0], zs[0], sph + HALF), (H[1] * zs[1], zs[1], sph + HALF)]

        def prof_du2(zs, H):                # dPsi2/dx2 = sum -H2[q] z_q sin(...)
            return [(-H[0] * zs[0], zs[0], sph), (-H[1] * zs[1], zs[1], sph)]

        def pint(A, B):
            s = mpc(0)
            for ca, fa, pa in A:
                for cb, fb, pb in B:
                    s += ca * cb * _rect_J(fa, pa, fb, pb)
            return s

        def xfac(g, r, X, kind):
            """x1 scalar at x1=X.  kind: 'u1'=cos(gX+r'), 'u1_1'=-g sin,
            'u2'=sin(gX+r'), 'u2_1'=+g cos."""
            ang = g * X + (r - 1) * mp.pi / 2
            if kind == "u1":
                return mp.cos(ang)
            if kind == "u1_1":
                return -g * mp.sin(ang)
            if kind == "u2":
                return mp.sin(ang)
            if kind == "u2_1":
                return g * mp.cos(ang)
            raise ValueError(kind)

        P = len(full)
        K = mp.zeros(2 * P, 2 * P)
        for pj, (gj, zj, H1j, H2j) in enumerate(bd):          # trial = column
            U1j = prof_u1(zj, H1j); U2j = prof_u2(zj, H2j)
            dU1j = prof_du1(zj, H1j); dU2j = prof_du2(zj, H2j)
            for pi_, (gi, zi, H1i, H2i) in enumerate(bd):     # test = row
                U1i = prof_u1(zi, H1i); U2i = prof_u2(zi, H2i)
                dU1i = prof_du1(zi, H1i); dU2i = prof_du2(zi, H2i)
                # x2-integrals (trial profile . test profile) needed below:
                I_u1j_u1i = pint(U1j, U1i)   # wall t1: u1_trial * u1_test,1 (x1 part)
                I_u1j_du2i = pint(U1j, dU2i) # wall t1: u1_trial * u2_test,2
                I_u2j_du1i = pint(U2j, dU1i) # wall t2: u2_trial * u1_test,2
                I_u2j_u2i = pint(U2j, U2i)   # wall t2: u2_trial * u2_test,1
                I_du2j_u1i = pint(dU2j, U1i) # tip  t3: u2_trial,2 * u1_test
                I_u1j_u1i_t = I_u1j_u1i      # tip  t3: u1_trial,1 * u1_test (same x2 prod)
                I_du1j_u2i = pint(dU1j, U2i) # tip  t4: u1_trial,2 * u2_test
                I_u2j_u2i_t = I_u2j_u2i      # tip  t4: u2_trial,1 * u2_test
                for r in (1, 2):                              # trial r' index
                    for rr in (1, 2):                         # test  r' index
                        # TIP at x1=+L (free): -(trial traction * test value)
                        # Unchanged for both bc values (x2=+/-b and x1=+L stay free).
                        tip = -(c11s * (xfac(gj, r, L, "u1_1") * xfac(gi, rr, L, "u1") * I_u1j_u1i_t
                                        + nu * xfac(gj, r, L, "u2") * xfac(gi, rr, L, "u1") * I_du2j_u1i)
                                + (xfac(gj, r, L, "u1") * xfac(gi, rr, L, "u2") * I_du1j_u2i
                                   + xfac(gj, r, L, "u2_1") * xfac(gi, rr, L, "u2") * I_u2j_u2i_t))
                        if self.bc == 'clamped_free':
                            # WALL at x1=-L (clamped): trial value * test derivative
                            # EXACT pre-bc path, formula unchanged.
                            wall = (c11s * (xfac(gj, r, -L, "u1") * xfac(gi, rr, -L, "u1_1") * I_u1j_u1i
                                            + nu * xfac(gj, r, -L, "u1") * xfac(gi, rr, -L, "u2") * I_u1j_du2i)
                                    + (xfac(gj, r, -L, "u2") * xfac(gi, rr, -L, "u1") * I_u2j_du1i
                                       + xfac(gj, r, -L, "u2") * xfac(gi, rr, -L, "u2_1") * I_u2j_u2i))
                            K[2 * pi_ + (rr - 1), 2 * pj + (r - 1)] = wall + tip
                        else:
                            # free_free: both short edges natural. Same traction
                            # bilinear form as tip, at x1=-L. Contour sign: both
                            # free ends contribute -int t_1b delta u_b dx2
                            # (Part 2 Eq.(1); n=+/-e1 is in n_a t_ab). No corner.
                            left = -(c11s * (xfac(gj, r, -L, "u1_1") * xfac(gi, rr, -L, "u1") * I_u1j_u1i_t
                                            + nu * xfac(gj, r, -L, "u2") * xfac(gi, rr, -L, "u1") * I_du2j_u1i)
                                    + (xfac(gj, r, -L, "u1") * xfac(gi, rr, -L, "u2") * I_du1j_u2i
                                       + xfac(gj, r, -L, "u2_1") * xfac(gi, rr, -L, "u2") * I_u2j_u2i_t))
                            K[2 * pi_ + (rr - 1), 2 * pj + (r - 1)] = left + tip
        return K

    @staticmethod
    def _conjugate_pair_cols(full, tol=1e-9):
        pairs = []; used = set()
        for p, z in enumerate(full):
            if p in used or abs(z.imag) < tol:
                continue
            for q in range(p + 1, len(full)):
                if q in used:
                    continue
                w = full[q]
                if abs(w.real - z.real) < tol and abs(w.imag + z.imag) < tol:
                    used.add(p); used.add(q)
                    for r in (0, 1):
                        pairs.append((2 * p + r, 2 * q + r))
                    break
        return pairs

    def _realify_conjugate_pairs(self, K, full):
        n = K.cols; I = mpc(0, 1); s2 = mp.sqrt(2)
        for (a, b) in self._conjugate_pair_cols(full):
            for k in range(n):
                xa, xb = K[k, a], K[k, b]
                K[k, a] = (xa + xb) / s2
                K[k, b] = (xa - xb) / (I * s2)
            for k in range(n):
                xa, xb = K[a, k], K[b, k]
                K[a, k] = (xa + xb) / s2
                K[b, k] = (xa - xb) / (I * s2)

    def equil_sigma(self, K, full=None):
        K = K.copy()
        if full is not None:
            self._realify_conjugate_pairs(K, full)
        n = K.cols
        for j in range(n):
            cn = mp.sqrt(sum(abs(K[i, j]) ** 2 for i in range(n)))
            if cn > 0:
                for i in range(n):
                    K[i, j] /= cn
        S = mp.svd(K, compute_uv=False)
        vals = [S[i] for i in range(S.rows)]
        return float(min(vals) / max(vals))

    # -- Lagrange-constrained detector (Eq.40/41) ------------------------------
    # The rectangular center-origin unconstrained K (Eq.36) fully DECOUPLES into
    # two identical sub-systems (paper p.154, footnote 3): its singular values
    # come in frozen pairs and the rank drops by TWO at a true mode, so a plain
    # sigma_min/sigma_max cannot localize modes for higher l/b (verified: flat at
    # l/b=2.5).  The paper removes this double root with an auxiliary constraint
    # -- the MEAN displacement of component d at the wall x1=-L vanishes, d=1
    # (u1) for symmetric, d=2 (u2) for antisymmetric (Eq.40).  The augmented
    # (2P+1)x(2P+1) matrix (Eq.42-43) has a simple root at each true mode; about
    # half of its roots are spurious (constraint-induced) and are discarded by
    # cross-checking against the unconstrained dip / a Ritz proximity cue.

    def _wall_mean_vector(self, full, Om, sym, lob):
        """Mean over x2 in [-pi/2,pi/2] of the constrained displacement component
        (d=1 u1 for sym, d=2 u2 for antisym) at the wall x1=-L, per basis column
        (branch p, r).  Returns a length-2P list of mpc."""
        eng = self.eng(sym)
        sph = mp.pi / 2 if sym else mpf(0)
        L = mp.pi * mpf(str(lob)) / 2
        HALF = mp.pi / 2
        Of = float(Om)

        def Msin(z, p):                     # int_{-pi/2}^{pi/2} sin(z x2 + p) dx2
            if abs(z) < mpf(10) ** -18:
                return mp.sin(p) * mp.pi
            return (-mp.cos(z * HALF + p) + mp.cos(-z * HALF + p)) / z

        vec = []
        for g in full:
            (z1, z2), H1, H2 = eng.H(complex(g), Of)
            z1 = mpc(z1); z2 = mpc(z2)
            if sym:                          # u1 = Psi1 cos(g x1 + r'); Psi1 = sum H1 sin(z x2+sph)
                Hq = (mpc(H1[0]), mpc(H1[1])); base_ph = sph
                x1kind = "cos"
            else:                            # u2 = Psi2 sin(g x1 + r'); Psi2 = sum H2 cos(z x2+sph)
                Hq = (mpc(H2[0]), mpc(H2[1])); base_ph = sph + HALF
                x1kind = "sin"
            mean_x2 = Hq[0] * Msin(z1, base_ph) + Hq[1] * Msin(z2, base_ph)
            for r in (1, 2):
                ang = mpc(g) * (-L) + (r - 1) * mp.pi / 2
                x1fac = mp.cos(ang) if x1kind == "cos" else mp.sin(ang)
                vec.append(x1fac * mean_x2 / mp.pi)
        return vec

    def assemble_constrained(self, full, Om, sym, lob):
        """(2P+1)x(2P+1) augmented matrix: unconstrained Eq.41 block + one
        symmetric Lagrange border (mean-displacement constraint, Eq.40).
        Cantilever-only: the constraint is a wall-mean at x1=-L."""
        if self.bc != 'clamped_free':
            raise ValueError(
                "RectIPAssembler.assemble_constrained is cantilever-only "
                f"(bc='clamped_free'), not bc={self.bc!r}")
        K2 = self.assemble(full, Om, sym, lob)
        P = len(full)
        n = 2 * P
        Kc = mp.zeros(n + 1, n + 1)
        for i in range(n):
            for j in range(n):
                Kc[i, j] = K2[i, j]
        vec = self._wall_mean_vector(full, Om, sym, lob)
        for j in range(n):
            Kc[n, j] = vec[j]
            Kc[j, n] = vec[j]
        Kc[n, n] = mpf(0)
        return Kc

    def equil_sigma_constrained(self, Kc, full):
        """sigma_min/sigma_max of the augmented matrix with two-sided
        realification of the conjugate-pair columns in the 2P block ONLY (the
        Lagrange index is a single real dof and is left untouched)."""
        Kc = Kc.copy()
        self._realify_conjugate_pairs(Kc, full)   # acts on 2P-indexed pairs only
        n = Kc.cols
        for j in range(n):
            cn = mp.sqrt(sum(abs(Kc[i, j]) ** 2 for i in range(n)))
            if cn > 0:
                for i in range(n):
                    Kc[i, j] /= cn
        S = mp.svd(Kc, compute_uv=False)
        vals = [S[i] for i in range(S.rows)]
        return float(min(vals) / max(vals))


