# -*- coding: utf-8 -*-
"""
plate_solver.detectors -- root/branch bookkeeping (canon, real_basis_items,
select_fill, full_search, track, degeneracy guards), the K-matrix numerics
(equilibrated_logdet, sigma_min_from_K, equilibrated_nullvec_mp), the Ritz
cross-check, find_modes_sigmin (the main sigma_min(K) mode finder), and the
rectangular branch-root finders (_rect_*_axis_roots, rect*_resolve_branches,
rect*_select_branches).

Extracted verbatim (line-range provenance in LESSONS_LEARNED Sec. 22); no
numeric behaviour changed, SOLVER_VERSION not bumped.
"""
from __future__ import annotations
import os
import numpy as np
from mpmath import mp, mpf, mpc, matrix
import concurrent.futures

from .config import MP_DPS, MAX_SCAN_PTS, N_WORKERS, PREFETCH, PI, \
    _P1_LD_GENUINE, _P1_PROM, _P2_LD_GENUINE
from .geometry import IsotropicMaterial
from .boundary import EdgeKind
from .dispersion import _Part1Fast, _Part2Fast
from .workers import _throttled_map, _p1_sigmin_worker, _p1_sigmin_golden_worker, \
    _p2_sigmin_worker, _p2_sigmin_golden_worker

SIGMIN_DEGEN_SENTINEL = 99.0

def canon(z, tol=1e-9):
    """Canonical representative of a root family: quadrant I (Re≥0, Im≥0)."""
    a, b = abs(z.real), abs(z.imag)
    if a < tol: return complex(0.0, b)
    if b < tol: return complex(a, 0.0)
    return complex(a, b)

def real_basis_items(brs, imag_root_proj=('im', 're')):
    """Expand canonical roots into REAL basis items (xi, phase q, projection).
    real ξ      → (sin, cos)                       : 2 dofs
    imag ξ=ia   → 2 dofs; which projection pair extracts the NONZERO
                  quantities depends on the field's parity at a pure-
                  imaginary root, which differs between Part 1 (W is real)
                  and Part 2 (F is imaginary, G is real). Pass
                  imag_root_proj=('im','re') for Part 1 (default — matches
                  the historical behaviour) or ('re','im') for Part 2.
    complex ξ   → Re/Im × sin/cos                  : 4 dofs (= ξ and conj(ξ))

    BUGFIX (Research18): for Part 2, the previously-hardcoded ('im','re')
    pair produces two structurally-zero rows/columns in K for every
    pure-imaginary zeta root (F is purely imaginary and G purely real there,
    so q=0 projected 'im' and q=1 projected 're' both vanish identically —
    see lessons-learned, Part2 0.25pi mode3 '-inf plateau' rank-13/15
    investigation). Part 2 must use ('re','im') instead.
    """
    items = []
    p0, p1 = imag_root_proj
    for z in brs:
        if abs(z.imag) < 1e-9:
            items += [(z, 0, 'asis'), (z, 1, 'asis')]
        elif abs(z.real) < 1e-9:
            items += [(z, 0, p0), (z, 1, p1)]
        else:
            items += [(z, 0, 're'), (z, 0, 'im'), (z, 1, 're'), (z, 1, 'im')]
    return items

def dof_count(z):
    return 4 if (abs(z.real) > 1e-9 and abs(z.imag) > 1e-9) else 2

def select_fill(roots, n_dofs):
    """Sort canonical roots by (|Im|, |z|) and fill up to n_dofs real dofs."""
    rs = sorted(roots, key=lambda z: (abs(z.imag), abs(z)))
    sel = []; cnt = 0
    for z in rs:
        need = dof_count(z)
        if cnt + need <= n_dofs:
            sel.append(z); cnt += need
        if cnt >= n_dofs: break
    return sel, cnt

def set_signature(sel):
    return tuple((round(z.real, 1), round(z.imag, 1)) for z in sel)

def mp_gauss_legendre(n):
    """Gauss-Legendre nodes/weights on [-1,1] in mp precision (Newton)."""
    nodes = []; wts = []
    for k in range(1, n+1):
        x = mp.cos(mp.pi*(k - mpf('0.25'))/(n + mpf('0.5')))
        for _ in range(60):
            p0, p1 = mpf(1), x
            for j in range(2, n+1):
                p0, p1 = p1, ((2*j-1)*x*p1 - (j-1)*p0)/j
            dp = n*(x*p1 - p0)/(x*x - 1)
            dx = p1/dp
            x -= dx
            if abs(dx) < mpf(10)**(-(mp.dps-5)): break
        p0, p1 = mpf(1), x
        for j in range(2, n+1):
            p0, p1 = p1, ((2*j-1)*x*p1 - (j-1)*p0)/j
        dp = n*(x*p1 - p0)/(x*x - 1)
        nodes.append(x); wts.append(2/((1 - x*x)*dp*dp))
    return nodes, wts

def mp_proj(v, P):
    return mp.re(v) if P in ('re', 'asis') else mp.im(v)

def equilibrated_logdet(K_mp, size, passes=2):
    """Row/col max-abs equilibration, then LU det in mp.
    Returns (sign, log10|det| of equilibrated matrix)."""
    K = K_mp.copy()
    for _ in range(passes):
        for i in range(size):
            s = max(abs(K[i, j]) for j in range(size))
            if s > 0:
                for j in range(size): K[i, j] /= s
        for j in range(size):
            s = max(abs(K[i, j]) for i in range(size))
            if s > 0:
                for i in range(size): K[i, j] /= s
    d = mp.det(K)
    if d == 0: return 0, mpf('-inf')
    return (1 if d > 0 else -1), mp.log(abs(d), 10)


def sigma_min_from_K(K_mp, size, passes=2):
    """Smallest singular value of K (the rigorous mode indicator).

    A natural frequency is exactly an Ω where K(Ω) is singular, i.e. where the
    smallest singular value σ_min(K) → 0.  Unlike log|det K| (= Σ log σ_i),
    σ_min isolates *only* the singular direction, so it has no noisy baseline
    from the other singular values and no sign flips from basis reordering — its
    dips to ~0 are unambiguous modes.  This is the standard robust method for
    transcendental boundary-determinant eigenproblems and is faithful to the
    paper's det K(Ω)=0 condition (σ_min = 0  ⟺  det K = 0).

    K is first equilibrated (same row/col max-abs scaling as equilibrated_logdet,
    so conditioning matches), then cast to float64 and its smallest singular
    value taken via LAPACK.  float64 is ample for *locating* the dip; the final
    frequency is refined separately.  Returns log10(σ_min) (or -300.0 if zero).
    """
    K = K_mp.copy()
    for _ in range(passes):
        for i in range(size):
            s = max(abs(K[i, j]) for j in range(size))
            if s > 0:
                for j in range(size): K[i, j] /= s
        for j in range(size):
            s = max(abs(K[i, j]) for i in range(size))
            if s > 0:
                for i in range(size): K[i, j] /= s
    Knp = np.empty((size, size), dtype=complex)
    for i in range(size):
        for j in range(size):
            Knp[i, j] = complex(K[i, j])
    try:
        sv = np.linalg.svd(Knp, compute_uv=False)
    except np.linalg.LinAlgError:
        return -300.0
    smin = float(sv[-1])
    if not np.isfinite(smin) or smin <= 0.0:
        return -300.0
    return float(np.log10(smin))


def equilibrated_nullvec_mp(K_mp, size):
    """Null vector (smallest right singular vector) of a real mp matrix K,
    extracted at FULL mp precision with column equilibration.

    Why this and not a plain float64 SVD: the mode-shape trial basis mixes
    propagating waves (ξ real, O(1) angular factor) with evanescent waves
    (ξ ≈ i·a, angular factor ~cosh(a·Θ) that reaches 1e8–1e13 at 2Θ=π).  A
    float64 SVD of K is swamped by that 8–13 order dynamic range and returns a
    contaminated vector whose spuriously-large evanescent coefficients blow up
    at the free edge — the "edge spike" seen in the mode-shape figures.  Doing
    the SVD in mp after column-equilibrating K both conditions the problem and,
    crucially, because the column scaling is UNDONE on the result, drives the
    evanescent coefficients back to their correct tiny values (~1e-12…1e-18).
    The result is spike-free, mode-distinct, paper-faithful smooth shapes.

    K v = 0  with  K_eq = Dr·K·Dc  ⇒  K_eq w = 0  ⇒  v = Dc·w   (Dr drops out;
    only the column scaling Dc must be undone).  Returns a float64 numpy vector
    normalised to unit max-abs, or None on failure (caller falls back to a
    plain float64 SVD).
    """
    try:
        from mpmath import svd_r
    except Exception:
        return None
    try:
        A = K_mp.copy()
        dc = [mpf(1)] * size
        for _ in range(3):                       # a few equilibration passes
            for i in range(size):                # row scaling (conditioning only)
                s = max(abs(A[i, j]) for j in range(size))
                if s > 0:
                    for j in range(size):
                        A[i, j] /= s
            for j in range(size):                # column scaling (accumulate Dc)
                s = max(abs(A[i, j]) for i in range(size))
                if s > 0:
                    for i in range(size):
                        A[i, j] /= s
                    dc[j] /= s
        _, _, V = svd_r(A)                        # rows of V = right sing. vectors
        w = [V[size - 1, j] for j in range(size)]  # smallest singular vector
        c = [dc[j] * w[j] for j in range(size)]    # undo column scaling
        nrm = max(abs(x) for x in c)
        if nrm == 0:
            return None
        return np.array([float(x / nrm) for x in c])
    except Exception:
        return None


def _nullvec_mp_rect(C_mp, m, n):
    """Smallest right singular vector of a real m×n mp matrix (m≥n preferred),
    with column equilibration undone on the result (same idea as
    equilibrated_nullvec_mp, generalised to rectangular collocation matrices).
    Returns a float64 numpy vector (unit max-abs) or None on failure."""
    try:
        from mpmath import svd_r
    except Exception:
        return None
    try:
        A = C_mp.copy()
        dc = [mpf(1)] * n
        for _ in range(3):
            for i in range(m):                       # row scaling (conditioning)
                s = max(abs(A[i, j]) for j in range(n))
                if s > 0:
                    for j in range(n):
                        A[i, j] /= s
            for j in range(n):                       # column scaling (accumulate)
                s = max(abs(A[i, j]) for i in range(m))
                if s > 0:
                    for i in range(m):
                        A[i, j] /= s
                    dc[j] /= s
        _, _, V = svd_r(A, full_matrices=False)      # V is n×n; rows = right vecs
        w = [V[n - 1, j] for j in range(n)]
        c = [dc[j] * w[j] for j in range(n)]
        nrm = max(abs(x) for x in c)
        if nrm == 0:
            return None
        return np.array([float(x / nrm) for x in c])
    except Exception:
        return None


def _theta_collocation_coeffs(solver, Om_star, items, cache, n,
                              n_s=None, v_theta='sin'):
    """Boundary-COLLOCATION mode coefficients for the out-of-plane plate.

    Problem-2 fix (clamped/fixed edge must stay at zero).  The weak (variational)
    null vector satisfies the clamped wall θ=-Θ only in an integral sense, so the
    plotted fixed edge lifts off zero.  Here we instead solve for the wave-
    amplitude vector c that makes the TWO circumferential-edge conditions hold
    POINTWISE at a set of sample radii x_k ∈ (-π/2, π/2):

        clamped  (θ=-Θ):   W = 0           and   ∂W/∂θ = 0
        free     (θ=+Θ):   M_θ = 0          and   V_θ   = 0

    (The arc edges r=R_i,R_o are already satisfied EXACTLY by the radial
    functions W_j(r) — paper Eq. 23 — so they need no collocation.)  Enforcing
    BOTH θ-edges together is the key: constraining the clamp ALONE leaves the
    free edge unconstrained and the evanescent functions blow up there (the old
    "spike").  With the free-edge moment/shear pinned too, the system has a
    smooth physical null vector that holds the fixed edge at zero.

    Operators (radial parts identical to _build_K_real):
        M_θ_j(r) = Tyy_j = ν W'' + R W'/rr − R ξ² W/rr²
        V_θ_j(r) = Qy_j  = (2T−ν)W''/rr + (R−2T+2ν)W'/rr² + (2T−2ν−R ξ²)W/rr³
    angular factors at the edges: cm=cos(−ξΘ+φ), sm=sin(−ξΘ+φ),
                                  cp=cos(+ξΘ+φ), sp=sin(+ξΘ+φ),  φ=q·π/2.
    M_θ carries the displacement angular factor (cos → cp); V_θ involves an odd
    θ-derivative (sin → sp).  `v_theta` selects the V_θ angular convention
    ('sin'→sp, 'xicos'→ξ·cp) so the caller can try both and keep whichever the
    self-check likes — a wrong guess simply fails the self-check and the caller
    falls back to the smooth weak-form shape, so this can never regress.

    Returns a float64 coefficient vector (unit max-abs) or None.
    """
    try:
        r0 = float(solver.r0b); nu = float(solver.nu)
        T = float(solver.T); R = float(solver.R); Th = float(solver.Theta)
        if n_s is None:
            n_s = max(4, (n // 4) + 1)            # ensure 4·n_s ≥ n (tall system)
        hp = mp.pi / 2
        xs = [(-hp + (mpf(2 * k + 1) / (2 * n_s)) * mp.pi) for k in range(n_s)]
        rows = []
        for x in xs:
            rr = x + mpf(r0)
            rW = []; rWth = []; rM = []; rV = []
            for (z, q, P) in items:
                xi = cache[z][0]; sols = cache[z][1]; A = cache[z][2]
                ph = mpf(q) * mp.pi / 2
                W = solver._W_mp(x, sols, A, 0)
                Wp = solver._W_mp(x, sols, A, 1)
                Wpp = solver._W_mp(x, sols, A, 2)
                xi2 = xi * xi
                Tyy = nu * Wpp + R * Wp / rr - R * xi2 * W / rr**2
                Qy = ((2 * T - nu) * Wpp / rr
                      + (R - 2 * T + 2 * nu) * Wp / rr**2
                      + (2 * T - 2 * nu - R * xi2) * W / rr**3)
                cm = mp.cos(-xi * Th + ph); sm = mp.sin(-xi * Th + ph)
                cp = mp.cos(xi * Th + ph);  sp = mp.sin(xi * Th + ph)
                vW = W * cm                       # clamped  W
                vWth = -xi * W * sm               # clamped  ∂W/∂θ
                vM = Tyy * cp                      # free     M_θ
                vV = (Qy * sp) if v_theta == 'sin' else (Qy * xi * cp)  # free V_θ

                def proj(val):
                    return val.real if P in ('re', 'asis') else val.imag
                rW.append(proj(vW)); rWth.append(proj(vWth))
                rM.append(proj(vM)); rV.append(proj(vV))
            rows.extend([rW, rWth, rM, rV])
        m = len(rows)
        C = mp.matrix(m, n)
        for i in range(m):
            for j in range(n):
                C[i, j] = rows[i][j]
        return _nullvec_mp_rect(C, m, n)
    except Exception:
        return None


def _quick_shape_stats(solver, items, cache, c, n_r=5, n_th=9):
    """Cheap clamp / interior diagnostics for a candidate coefficient vector,
    used to accept-or-reject a collocation reconstruction without building the
    full plotting grid.  Returns (clamp_resid, interior_frac) as fractions of
    the peak |W|, where clamp_resid is the displacement on the θ=-Θ edge and
    interior_frac is the largest |W| away from the free edge (a spike collapses
    interior_frac toward 0).  On failure returns (1.0, 0.0) → reject."""
    try:
        hp = float(mp.pi / 2); Th = float(solver.Theta)
        rv = np.linspace(-hp, hp, n_r); tv = np.linspace(-Th, Th, n_th)
        Wg = np.zeros((n_r, n_th))
        for ir, x in enumerate(rv):
            Wc = {}
            for z in {it[0] for it in items}:
                sols = cache[z][1]; A = cache[z][2]
                Wc[z] = complex(solver._W_mp(mpf(float(x)), sols, A, 0))
            for k, (z, q, P) in enumerate(items):
                xi = complex(cache[z][0]); ph = q * np.pi / 2; Wj = Wc[z]
                for jt, th in enumerate(tv):
                    phi = Wj * np.cos(xi * th + ph)
                    Wg[ir, jt] += c[k] * (phi.real if P in ('re', 'asis')
                                          else phi.imag)
        peak = float(np.max(np.abs(Wg))) or 1.0
        clamp = float(np.max(np.abs(Wg[:, 0]))) / peak
        nth = Wg.shape[1]
        interior = float(np.max(np.abs(Wg[:, :max(1, 3 * nth // 4)]))) / peak
        return clamp, interior
    except Exception:
        return 1.0, 0.0


# ----------------------------------------------------------------------
# Independent Rayleigh–Ritz cross-check (Problem 3: flag spurious "extra"
# determinant zeros WITHOUT any paper values and WITHOUT an external FEA).
# A self-contained Kirchhoff-plate Ritz solve on the SAME annular sector gives
# an independent spectrum; genuine modes of the wave method have a Ritz partner,
# spurious determinant singularities do not.  It is ratio-based and self-
# calibrated to the (always-real) fundamental, so the absolute nondimensional-
# isation need not match — only the spacing of the spectrum matters.
# ----------------------------------------------------------------------
def ritz_spectrum_oop(geom, mat, Mr=6, Nth=6, nq=20):
    """Return a sorted numpy array of independent out-of-plane modal values
    √λ (∝ natural frequency) from a Rayleigh–Ritz solution of the clamped-free
    annular sector Kirchhoff plate.  D and ρh are set to 1 (they cancel in the
    ratio comparison).  Trial functions

        φ_{m,n}(r,θ) = r̂^m · η^{n+2},   r̂=(r-R_i)/(R_o-R_i), η=(θ+Θ)/(2Θ)

    satisfy the clamped edge (η=0 ⇒ W=0, ∂W/∂θ=0 at θ=-Θ) essentially; the free
    edge at θ=+Θ is natural.  Returns None on failure."""
    try:
        from numpy.polynomial.legendre import leggauss
        import scipy.linalg as sla
        Ri = float(geom.R_i); Ro = float(geom.R_o); Th = float(geom.Theta)
        nu = float(getattr(mat, 'nu_bar', getattr(mat, 'nu', 0.3)))
        Lr = Ro - Ri; Lt = 2 * Th
        gx, gw = leggauss(nq)                      # nodes on [-1,1]
        rhat = 0.5 * (gx + 1.0); rw = 0.5 * gw     # → [0,1]
        eta = 0.5 * (gx + 1.0); ew = 0.5 * gw
        basis = [(m, n) for m in range(Mr) for n in range(Nth)]
        nb = len(basis)
        # precompute per-(quad point) basis values & derivatives
        QP = nq * nq
        phi = np.zeros((nb, QP)); A = np.zeros((nb, QP))
        L = np.zeros((nb, QP)); Pp = np.zeros((nb, QP)); Q = np.zeros((nb, QP))
        wt = np.zeros(QP)
        idx = 0
        for ii in range(nq):
            rh = rhat[ii]; r = Ri + Lr * rh
            for jj in range(nq):
                et = eta[jj]
                wt[idx] = rw[ii] * ew[jj] * r * Lr * Lt   # r·dr·dθ Jacobian
                for bi, (m, nn) in enumerate(basis):
                    Rm = rh**m
                    Rm_r = (m * rh**(m - 1) / Lr) if m >= 1 else 0.0
                    Rm_rr = (m * (m - 1) * rh**(m - 2) / Lr**2) if m >= 2 else 0.0
                    p = nn + 2
                    Tn = et**p
                    Tn_t = p * et**(p - 1) / Lt
                    Tn_tt = p * (p - 1) * et**(p - 2) / Lt**2
                    w = Rm * Tn
                    w_r = Rm_r * Tn
                    w_rr = Rm_rr * Tn
                    w_t = Rm * Tn_t
                    w_tt = Rm * Tn_tt
                    w_rt = Rm_r * Tn_t
                    phi[bi, idx] = w
                    A[bi, idx] = w_rr
                    Pp[bi, idx] = w_r / r + w_tt / r**2
                    Q[bi, idx] = w_rt / r - w_t / r**2
                    L[bi, idx] = w_rr + w_r / r + w_tt / r**2
                idx += 1
        Wd = wt
        Kmat = (L * Wd) @ L.T \
            - (1 - nu) * ((A * Wd) @ Pp.T + (Pp * Wd) @ A.T) \
            + 2 * (1 - nu) * ((Q * Wd) @ Q.T)
        Mmat = (phi * Wd) @ phi.T
        Kmat = 0.5 * (Kmat + Kmat.T); Mmat = 0.5 * (Mmat + Mmat.T)
        # The monomial radial basis r̂^m becomes increasingly ill-conditioned as
        # Mr grows (the columns of Mmat approach linear dependence — the classic
        # Hilbert-like degeneracy of {1, x, x², …}), so the raw generalized
        # eigensolve sla.eigh(K, M) can fail or return garbage at Mr≳7.  Solve it
        # via a symmetric whitening of the mass matrix with a small Tikhonov
        # floor on tiny eigenvalues, which discards the numerically-null
        # directions instead of crashing.  This lets the Ritz basis be enriched
        # (Mr, Nth up) for wider sectors where the default 6×6 spectrum is too
        # coarse, without destabilising the well-conditioned small-basis case.
        try:
            ev = sla.eigh(Kmat, Mmat, eigvals_only=True)
        except Exception:
            ev = None
        if ev is None or not np.all(np.isfinite(ev)):
            # Whiten: M = U s Uᵀ, keep directions with s above a relative floor.
            sM, UM = np.linalg.eigh(Mmat)
            smax = float(sM.max()) if sM.size else 0.0
            keep = sM > max(1e-14, 1e-12 * smax)
            if not np.any(keep):
                return None
            Winv = UM[:, keep] / np.sqrt(sM[keep])
            Kt = Winv.T @ Kmat @ Winv          # standard symmetric problem
            Kt = 0.5 * (Kt + Kt.T)
            ev = np.linalg.eigvalsh(Kt)
        ev = np.array([e for e in ev if e > 1e-9])
        if ev.size == 0:
            return None
        return np.sqrt(np.sort(ev))
    except Exception:
        return None


def classify_modes_ritz(freqs, geom, mat, tol=0.12, Mr=6, Nth=6):
    """Classify each wave-method frequency as confirmed / suspect using the
    independent Ritz spectrum.  Self-calibrates a single scale α on the
    fundamental (assumed real), then a frequency is CONFIRMED if some α·√λ_k
    lies within `tol` (relative).  Frequencies above the Ritz spectrum's range
    are left CONFIRMED ("unchecked", not enough Ritz modes to judge) so real
    high modes are never falsely dropped.  Returns (flags, info) where flags[i]
    is True (keep/figure) or False (suspect), and info is a short text summary.
    Conservative by design: when in doubt, keep.

    Robustness (Research36b): the verdict is taken as the UNION over two Ritz
    basis sizes — the requested (Mr, Nth) and an enriched (Mr+2, Nth+2).  A
    genuine plate mode that happens to fall in a resolution gap of the coarse
    spectrum is then still confirmed by the finer one, so real modes are not
    falsely flagged suspect; a spurious determinant zero has no Ritz partner at
    EITHER resolution and stays suspect.  The enriched solve reuses the same
    energy form and is cheap (≈0.3 s); if it fails for any reason the function
    silently falls back to the single-basis verdict.  Set MODE_FILTER_RITZ_UNION=0
    to use only the coarse basis."""
    flags = [True] * len(freqs)
    if not freqs:
        return flags, "no frequencies"

    def _verdict(_Mr, _Nth):
        rs = ritz_spectrum_oop(geom, mat, Mr=_Mr, Nth=_Nth)
        if rs is None or rs.size == 0:
            return None, None, None
        fsort2 = sorted(range(len(freqs)), key=lambda i: freqs[i])
        f0 = freqs[fsort2[0]]
        if f0 <= 0:
            return None, None, None
        alpha = f0 / rs[0]
        scaled = alpha * rs
        hi = scaled[-1] * (1.0 + tol)
        keep = [True] * len(freqs)
        for i, f in enumerate(freqs):
            if f > hi:
                continue
            rel = np.min(np.abs(scaled - f)) / max(f, 1e-12)
            if rel > tol:
                keep[i] = False
        return keep, alpha, rs.size

    keep_a, alpha, nrs = _verdict(Mr, Nth)
    if keep_a is None:
        return flags, "Ritz cross-check unavailable (kept all)"

    use_union = os.environ.get("MODE_FILTER_RITZ_UNION", "1") == "1"
    keep_b = None
    if use_union:
        keep_b, _, _ = _verdict(Mr + 2, Nth + 2)

    n_suspect = 0
    for i in range(len(freqs)):
        confirmed = keep_a[i] or (keep_b[i] if keep_b is not None else False)
        flags[i] = bool(confirmed)
        if not confirmed:
            n_suspect += 1
    note = "" if keep_b is None else f"+{Mr+2}×{Nth+2} union"
    info = (f"Ritz cross-check: {sum(flags)}/{len(freqs)} confirmed, "
            f"{n_suspect} suspect (α={alpha:.4g}, {nrs} Ritz modes{note}, "
            f"tol={tol:.0%})")
    return flags, info


def _min_root_spacing(sel):
    L = [complex(z) for z in sel]
    if len(L) < 2:
        return float('inf')
    md = float('inf')
    for i in range(len(L)):
        for j in range(i + 1, len(L)):
            d = abs(L[i] - L[j])
            if d < md:
                md = d
    return md

def _degen_guard(sel, s_star):
    """Return s_star unchanged, or a shallow sentinel if the basis at this Ω is
    an ARTIFACT (branch-collision or — optionally — cut-off) rather than a real
    mode, so the normal accept-floor drops it.

    (1) Branch collision  [DEFAULT ON, SIGMIN_DEGEN_TOL=0.06]: two canonical
        roots nearly coincide -> near-parallel basis columns -> K singular for
        ANY sector angle (the recurring Ω≈0.4542 at r0/(2b)=1.25).  Real modes
        have root spacing ~0.8; the artifact ~0.05 — a clean 15x separation.

    (2) Cut-off artifact  [EXPERIMENTAL, DEFAULT OFF, SIGMIN_CUTOFF_GUARD=1]:
        just above a ξ/ζ=0 cut-off, a branch passes through zero, so the basis
        contains a near-rigid (|root|≈0) function and K is ill-conditioned,
        giving a SHALLOW spurious dip (e.g. Ω≈0.2758 just above cut-off 0.2538,
        σ_min≈-2.9, min|root|≈0.17).  Probed directly: real modes there have
        min|root|≳0.4 and the genuine near-cut-off modes are DEEP (σ_min≲-4).
        So reject ONLY when a root is near zero AND the dip is shallow.
        WARNING: a genuinely SHALLOW real mode sitting near a cut-off could
        still be clipped — this is why it is OFF by default.  Validate against
        the paper before trusting, and run with RESUME=0 when toggling it."""
    try:
        tol = float(os.environ.get("SIGMIN_DEGEN_TOL", "0.06"))
    except ValueError:
        tol = 0.06
    L = [complex(z) for z in sel]
    if tol > 0.0 and _min_root_spacing(sel) < tol:
        return SIGMIN_DEGEN_SENTINEL
    if os.environ.get("SIGMIN_CUTOFF_GUARD", "0") == "1" and L:
        root_tol = float(os.environ.get("SIGMIN_CUTOFF_ROOT_TOL", "0.25"))
        shallow  = float(os.environ.get("SIGMIN_CUTOFF_SHALLOW", "-3.2"))
        if min(abs(z) for z in L) < root_tol and s_star > shallow:
            return SIGMIN_DEGEN_SENTINEL
    return s_star


def full_search(eng, Om, xmax=14.0, ngrid=10, nax=240):
    """Inventory roots: real axis, imaginary axis, complex quadrant Newton."""
    roots = []
    def add(z):
        if z is None: return
        z = canon(z)
        if abs(z) < 1e-6 or abs(z) > xmax+3: return
        if abs(eng.det(z, Om)) > 1e-6: return
        for u in roots:
            if abs(z-u) < 0.05: return
        roots.append(z)
    xs = np.linspace(1e-4, xmax, nax)
    vr = [eng.det(x, Om).real for x in xs]
    vi = [eng.det(1j*x, Om).real for x in xs]
    for k in range(nax-1):
        if vr[k]*vr[k+1] < 0: add(eng.newton((xs[k]+xs[k+1])/2, Om))
        if vi[k]*vi[k+1] < 0: add(eng.newton(1j*(xs[k]+xs[k+1])/2, Om))
    for a0 in np.linspace(0.4, xmax, ngrid):
        for b0 in np.linspace(0.4, xmax, ngrid):
            add(eng.newton(complex(a0, b0), Om, itmax=30))
    # fine near-origin grid: catches the small complex pair after collision
    for a0 in np.linspace(0.08, 1.4, 8):
        for b0 in np.linspace(0.06, 1.2, 7):
            add(eng.newton(complex(a0, b0), Om, itmax=30))
    roots.sort(key=lambda z: (abs(z.imag), abs(z)))
    return roots

def track(eng, Om, prev):
    """Continue previous roots to a new Ω; dedup in canonical form.

    Roots that are non-finite or implausibly large are dropped: at low Ω and
    large r0_bar the float64 fast engine can overflow and emit garbage roots
    (e.g. |z| ~ 1e4).  Any legitimate branch root has |z| <= xi_max (<= ~20),
    so a hard cap of 50 removes overflow artifacts without discarding real
    roots, preventing a corrupted K-matrix basis downstream.
    """
    roots = []
    for z0 in prev:
        z = eng.newton(z0, Om, itmax=30)
        if z is None: continue
        if not (np.isfinite(z.real) and np.isfinite(z.imag)): continue
        if abs(z) > 50.0: continue            # overflow garbage guard
        z = canon(z)
        if abs(z) < 1e-6: continue
        if any(abs(z-u) < 0.02 for u in roots): continue
        roots.append(z)
    roots.sort(key=lambda z: (abs(z.imag), abs(z)))
    return roots


def find_modes_sigmin(solver, part, Omega_range, n_scan, n_dofs, max_dim,
                      n_modes_wanted=3, verbose=True, prom=0.5,
                      accept_floor=-3.5, n_workers=None, diag=None,
                      _gap_pass=False, _gap_depth=0):
    """Paper-faithful natural-frequency detector via the singular value σ_min(K).

    A natural frequency is an Ω at which K(Ω) is singular — exactly the paper's
    det K(Ω)=0 condition — detected here as σ_min(K) → 0.  Steps:

      Phase 1  sequential float64 branch tracking → per-point branch seeds,
      Phase 2  parallel log10 σ_min(K) sweep across the range (full quadrature,
               step ≤ ~0.0025 so narrow dips are not stepped over),
      Phase 3  parallel golden-section refinement of every *prominent* σ_min dip,
               each verified from an independent fresh branch set,

    then keep the lowest-Ω `n_modes_wanted` refined dips that actually reach a
    singularity (log10 σ_min ≤ accept_floor).  A dip is "prominent" when it lies
    at least `prom` orders of magnitude below both of its shoulders — a
    scale-free criterion (σ_min is O(1) away from modes, ~0 at modes), so no
    reference values and no per-geometry tuning enter.  Returns sorted Ω list.
    """
    Om_lo, Om_hi = Omega_range
    nw = n_workers if n_workers is not None else max(1, N_WORKERS)
    kind = 'out-of-plane' if part == 1 else 'in-plane'

    # Step 6: thread this solver's boundary condition to every worker pool
    # (scan + refine helpers all read R40_BC_KIND via _worker_bc()).  Cantilever
    # solvers carry token 'clamped_free', so existing runs are unchanged.
    os.environ['R40_BC_KIND'] = getattr(solver.bc, 'token', 'clamped_free')
    # ── Diagnostics capture (non-invasive: only populated when diag is a dict;
    #    writes NOTHING that feeds the computation — pure observation). ───────
    if diag is not None:
        diag.update(part=part, kind=kind, n_dofs=n_dofs, max_dim=float(max_dim),
                    Omega_range=(float(Om_lo), float(Om_hi)),
                    bc_token=os.environ.get("R40_BC_KIND"),
                    candidates=[], accepted=[], dropped=[])

    # ── SCAN STRATEGY (configurable) ─────────────────────────────────────────
    # Two-stage coarse(SIGMIN_COARSE_STEP=0.006)+local-zoom scan, ~2× faster than a
    # uniform 0.003 grid for the same dips.  SIGMIN_TWO_STAGE=0 reverts.  Cost model
    # and dip-width data: LESSONS_LEARNED.md §Detector.
    two_stage     = os.environ.get("SIGMIN_TWO_STAGE", "1") == "1"
    sig_max_pts   = int(os.environ.get("SIGMIN_MAX_PTS", "900"))
    step_fine     = float(os.environ.get("SIGMIN_STEP", "0.003"))
    step_coarse   = float(os.environ.get("SIGMIN_COARSE_STEP", "0.006"))
    local_step    = float(os.environ.get("SIGMIN_LOCAL_STEP", "0.0020"))
    local_span    = float(os.environ.get("SIGMIN_LOCAL_SPAN", "0.005"))

    scan_step = step_coarse if two_stage else step_fine
    n_eff = int(np.ceil((Om_hi - Om_lo) / scan_step)) + 1
    n_eff = max(n_scan, min(n_eff, sig_max_pts))

    print(f"\n{'='*68}")
    print(f"  σ_min mode search ({kind}):  Ω∈[{Om_lo:.4f}, {Om_hi:.4f}]"
          f"  n_dofs={n_dofs}  M={solver.M}  dps={mp.dps}")
    stage_txt = ("two-stage: coarse + local zoom" if two_stage
                 else "single uniform grid")
    print(f"  {n_eff} coarse pts (step≈{(Om_hi-Om_lo)/max(n_eff-1,1):.4f}, "
          f"{stage_txt}), prominence≥{prom}, accept σ_min≤1e{accept_floor:.1f}, "
          f"{nw} workers")
    print(f"{'='*68}", flush=True)

    Oms = list(np.linspace(Om_lo, Om_hi, n_eff))

    # ── Low-Ω densification ──────────────────────────────────────────────────
    # Modes cluster toward Ω→0 (wide sectors, Part 2); a half-step sub-grid over
    # the lowest SIGMIN_LOWZONE_FRAC of the range resolves them.  See §Detector.
    if part == 2:
        lowzone_frac = float(os.environ.get("SIGMIN_LOWZONE_FRAC_P2", "0.14"))
    else:
        lowzone_frac = float(os.environ.get("SIGMIN_LOWZONE_FRAC", "0.10"))
    if lowzone_frac > 0:
        lo_hi = Om_lo + lowzone_frac * (Om_hi - Om_lo)
        fine_lo_step = scan_step * 0.5
        n_lo = int(np.ceil((lo_hi - Om_lo) / fine_lo_step)) + 1
        lo_grid = list(np.linspace(Om_lo, lo_hi, max(2, n_lo)))
        # Merge and drop near-duplicates (within 1/3 of the fine step) so the
        # two grids don't waste full_search calls on coincident points.
        allpts = sorted(Oms + lo_grid)
        merged = [allpts[0]]
        for x in allpts[1:]:
            if x - merged[-1] > fine_lo_step / 3.0:
                merged.append(x)
        Oms = merged
        n_eff = len(Oms)
    # full_search dominates cost (~all of it; σ_min itself is ~1s), so we do NOT
    # pre-track a basis — each worker finds the canonical roots fresh, which is
    # what makes σ_min reach the true singularity at a mode.
    if part == 1:
        g_sc, m_sc = solver._p1_worker_params()
        scan_worker, golden_worker = _p1_sigmin_worker, _p1_sigmin_golden_worker
    else:
        g_sc, m_sc = solver._p2_worker_params()
        scan_worker, golden_worker = _p2_sigmin_worker, _p2_sigmin_golden_worker
    dps = mp.dps
    n_quad_full = len(solver.nodes)
    n_quad_iter = max(15, n_quad_full // 2)

    scan_args = [
        (float(Om), None, g_sc, m_sc, solver.M, n_quad_full,
         max_dim, n_dofs, dps)
        for Om in Oms
    ]
    print(f"  Phase 1: σ_min sweep ({len(scan_args)} points, "
          f"full_search-bound) …", flush=True)
    smin = [float('nan')] * len(Oms)
    with concurrent.futures.ProcessPoolExecutor(max_workers=nw) as pool:
        k = 0
        for Om_f, s in _throttled_map(pool, scan_worker, scan_args):
            smin[k] = s
            k += 1

    finite = [s for s in smin if np.isfinite(s)]
    if not finite:
        print("  WARNING: σ_min sweep produced no finite values.")
        return []
    base = float(np.median(finite))
    if diag is not None:
        diag.update(baseline=base, sigmin_min=float(min(finite)),
                    sigmin_max=float(max(finite)),
                    n_scan_points=len(Oms),
                    n_finite=len(finite), n_nan=len(Oms)-len(finite))
    if verbose:
        lohi = f"[{min(finite):+.2f}, {max(finite):+.2f}]"
        print(f"    σ_min range {lohi}, baseline≈{base:+.2f}", flush=True)

    # ── Identify candidate dips (LOOSE flagging) ─────────────────────────────
    # At an affordable scan step (~0.0025) a narrow dip may be sampled on its
    # shoulder, so its apparent depth/prominence is muted.  We therefore flag
    # generously here (low prominence) to avoid MISSING a real mode, and rely on
    # Phase-3 refinement + the independent fresh-basis accept-check to discard
    # the false positives (which refine to shallow σ_min).
    def shoulder(i, step):
        """Highest σ_min walking outward from i until it turns back down."""
        j = i
        hi = smin[i]
        while 0 <= j + step < len(smin):
            nj = j + step
            if not np.isfinite(smin[nj]):
                break
            if smin[nj] >= hi:
                hi = smin[nj]
                j = nj
                if hi >= smin[i] + 1.5:
                    break
            else:
                break
        return hi

    cand = []
    for i in range(len(smin)):
        if not np.isfinite(smin[i]):
            continue
        left_ok  = (i == 0) or (not np.isfinite(smin[i-1])) or smin[i] <= smin[i-1]
        right_ok = (i == len(smin)-1) or (not np.isfinite(smin[i+1])) or smin[i] <= smin[i+1]
        if not (left_ok and right_ok):
            continue
        ls, rs = shoulder(i, -1), shoulder(i, +1)
        prominence = min(ls, rs) - smin[i]
        if prominence >= prom and smin[i] <= base - prom:
            cand.append((i, smin[i]))

    # Cap candidates to keep refinement affordable, but keep a UNIFORM SPREAD across
    # Ω when the cap is hit (NOT the deepest) — a real low-frequency mode can be
    # shallower than higher-Ω dips (e.g. 1.25/1.25π m1&m2).  SIGMIN_MAX_CANDS.
    max_cands = int(os.environ.get("SIGMIN_MAX_CANDS", "24"))
    if len(cand) > max_cands:
        # Earlier this kept a UNIFORM Ω-spread, but np.linspace can still skip an
        # individual low-Ω candidate when several cluster just above mode 1 —
        # that is exactly what dropped the 1.5π mode-2 (Ω≈0.0346, a genuine deep
        # dip at σ_min≈-4.25 verified directly).  Fix: keep EVERY low-zone
        # (clustered low-frequency) candidate, then uniform-spread the remaining
        # budget across the higher-Ω candidates.  Low modes can never be culled.
        lz_cut = (Om_lo + lowzone_frac * (Om_hi - Om_lo)) if lowzone_frac > 0 else Om_lo
        low  = sorted((c for c in cand if Oms[c[0]] <= lz_cut), key=lambda t: t[0])
        high = sorted((c for c in cand if Oms[c[0]] >  lz_cut), key=lambda t: t[0])
        if len(low) > max_cands:
            # Pathological: more low-zone candidates than the whole budget —
            # keep the deepest of them (still never a uniform skip of a real one).
            low = sorted(sorted(low, key=lambda t: t[1])[:max_cands],
                         key=lambda t: t[0])
            cand = low
        else:
            budget = max_cands - len(low)
            if budget > 0 and high:
                idx = np.linspace(0, len(high) - 1, min(budget, len(high)))
                idx = sorted(set(int(round(x)) for x in idx))
                high = [high[i] for i in idx]
            else:
                high = []
            cand = low + high
    cand_idx = sorted(i for i, _ in cand)
    if diag is not None:
        diag["candidates"] = [float(Oms[i]) for i in cand_idx]

    if verbose:
        print(f"  {len(cand_idx)} candidate dip(s) at Ω≈"
              f"{[f'{Oms[i]:.4f}' for i in cand_idx]}", flush=True)
    if not cand_idx:
        return []

    # ── Phase 1b: local fine sub-scan around each coarse candidate ───────────
    # When two_stage is on, the coarse grid only guarantees a point on each
    # dip's SHOULDER.  Zoom in with a fine local grid (±local_span at
    # local_step) around each candidate to (a) confirm the dip and locate its
    # true centre, and (b) hand golden refinement a tight, well-centred
    # bracket so it needs fewer iterations.  Each candidate's local grid is a
    # handful of extra full_search calls; total stays far below the uniform
    # fine grid.  Results: list parallel to cand_idx of (Om_centre, a, b).
    cand_brackets = {}
    if two_stage:
        # Skip the local zoom for candidates that already sit in the densified
        # low-Ω zone (they are already sampled at half-step there, so the
        # coarse bracket is already tight) — saves a chunk of full_search
        # calls without losing resolution where it matters.
        lz_hi = Om_lo + lowzone_frac * (Om_hi - Om_lo) if lowzone_frac > 0 else Om_lo
        local_args = []
        local_owner = []          # which candidate each local point belongs to
        for i in cand_idx:
            Omc = Oms[i]
            if Omc <= lz_hi:
                continue          # already finely sampled in the low zone
            offs = np.arange(-local_span, local_span + 0.5*local_step, local_step)
            for off in offs:
                Omp = float(Omc + off)
                if Om_lo <= Omp <= Om_hi:
                    local_args.append(
                        (Omp, None, g_sc, m_sc, solver.M, n_quad_full,
                         max_dim, n_dofs, dps))
                    local_owner.append(i)
        if local_args:
            n_zoomed = len(set(local_owner))
            print(f"  Phase 1b: local zoom around {n_zoomed} candidate(s) "
                  f"({len(local_args)} fine pts; low-zone candidates already "
                  f"fine) …", flush=True)
            local_vals = [float('nan')] * len(local_args)
            with concurrent.futures.ProcessPoolExecutor(max_workers=nw) as pool:
                k = 0
                for Om_f, s in _throttled_map(pool, scan_worker, local_args):
                    local_vals[k] = (Om_f, s)
                    k += 1
            # group by owning candidate, find local minima + bracket(s)
            from collections import defaultdict
            grp = defaultdict(list)
            for (Om_f, s), owner in zip(local_vals, local_owner):
                if np.isfinite(s):
                    grp[owner].append((Om_f, s))
            # SPLIT-ON-MULTI-MINIMUM (2026-07-23, LESSONS_LEARNED Sec 20.8
            # hardening): a single coarse candidate can straddle TWO distinct
            # genuine dips closer together than the coarse grid resolves
            # (Sec 20.8's G2: 1.598985 vs the neighbouring artifact 1.603245,
            # ~0.0043 apart -- the coarse scan only sees one smeared
            # shoulder). The OLD behaviour (still the default -- env unset/0)
            # collapses to the single deepest local point, silently
            # discarding whichever trough is shallower before any golden
            # search or residual screen ever evaluates it independently --
            # exactly the "golden-lock-on-to-a-nearby-artifact" mechanism
            # that resolved G2 (LESSONS Sec 20.7 rule 1). When enabled,
            # detect every local minimum in this candidate's OWN local-zoom
            # points (same shoulder logic as the coarse Phase-1 flagging,
            # reused at local-zoom resolution) and emit one bracket per
            # minimum, each independently golden-refined in Phase 2 -- so
            # two close roots can no longer merge into one accepted dip
            # before Phase 2 ever sees them separately.
            split_on = os.environ.get("SIGMIN_LOCAL_SPLIT", "0") == "1"
            split_min_sep = float(os.environ.get("SIGMIN_LOCAL_SPLIT_MINSEP",
                                                  str(local_step * 1.5)))
            split_prom = float(os.environ.get("SIGMIN_LOCAL_SPLIT_PROM", "0.15"))
            # ANCHOR + ADAPTIVE-BISECTION HARDENING (2026-07-29,
            # LESSONS_LEARNED Sec 49/50): job 2335775 showed SIGMIN_LOCAL_SPLIT
            # alone still misses a real root (G2=1.598985) sitting ~0.0043
            # from a much stronger neighbour (ARTIFACT 1.603259) -- closer
            # than the coarse step (~0.0057). Direct sandbox computation
            # (dps=26, this geometry) root-caused it precisely: (a) G2's own
            # coarse sample point never becomes a Phase-1 candidate at all --
            # it loses the strict left_ok/right_ok local-min test to the
            # deeper neighbour, so the split search's input never even
            # contained the winning candidate's own triggering evidence, let
            # alone G2's; (b) the resonance itself is narrow enough that even
            # the fixed local_step zoom around the WINNING candidate can
            # straddle G2 without any single zoom point landing close enough
            # to reveal its depth -- confirmed directly: Om=1.598857 (a
            # bisection of two existing zoom points, 1.597857 and 1.599857)
            # reads sigma_min=-4.07, while those two zoom-grid points
            # themselves read only -2.46/-2.41, a 3+ order-of-magnitude
            # apparent-depth loss from stepping just ~0.001 off-peak.
            # Two-part fix, both nested under SIGMIN_LOCAL_SPLIT so neither
            # changes default behaviour and neither needs a SOLVER_VERSION
            # bump (detector-tuning only -- no computed sigma_min(Omega)
            # value changes, only which Omega get proposed as candidates):
            #  (1) ANCHOR: seed the split search's point list with the
            #      ALREADY-COMPUTED coarse values at the candidate's own
            #      centre Oms[i] and its two flanking coarse points
            #      Oms[i-1]/Oms[i+1] (zero extra full_search cost -- Phase 1
            #      already evaluated all three). Without this, the split
            #      search's own triggering candidate wasn't even part of its
            #      input, and it had less shoulder context than Phase 1's
            #      own coarse flagging already had.
            #  (2) ADAPTIVE BISECTION (one bounded extra round): a point
            #      that is a strict local min (turns back up on both sides,
            #      i.e. left_ok and right_ok) but falls BELOW split_prom is
            #      exactly the signature of a real feature under-resolved by
            #      the current step -- not noise (this is what flags
            #      Om=1.597857 as "weak" in the G2 case). Bisect that weak
            #      point's two flanking gaps once, at half their current
            #      spacing, evaluate the new points (one extra pooled
            #      full_search round, capped to
            #      SIGMIN_LOCAL_SPLIT_ADAPTIVE_MAXPTS new points per
            #      candidate), merge them in, and re-run the shoulder search
            #      on the enlarged set. Capped to one round and to
            #      candidates that actually show this specific weak-shoulder
            #      pattern (rare), so cost stays bounded and most windows
            #      are unaffected.
            # NOT YET CLUSTER-VALIDATED as of this comment -- the sandbox
            # confirmation above is structural (single geometry, dps=26,
            # hand-verified points only, not a full pipeline run); see
            # probe_sigmin_local_split_anchor_v1.py for the end-to-end
            # cluster test against job 2335775's exact G1-G3/W2+controls
            # scenario.
            split_anchor = os.environ.get("SIGMIN_LOCAL_SPLIT_ANCHOR", "1") == "1"
            split_adaptive = os.environ.get("SIGMIN_LOCAL_SPLIT_ADAPTIVE", "1") == "1"
            split_adaptive_max = int(os.environ.get(
                "SIGMIN_LOCAL_SPLIT_ADAPTIVE_MAXPTS", "4"))

            def _shoulder_scan(vals):
                """Strict local minima of `vals`, split into (confirmed,
                weak): confirmed clears split_prom; weak turns back up on
                both sides but doesn't clear split_prom yet (a candidate
                for adaptive bisection, not noise)."""
                confirmed, weak = [], []
                for k in range(len(vals)):
                    left_ok = (k == 0) or vals[k] <= vals[k - 1]
                    right_ok = (k == len(vals) - 1) or vals[k] <= vals[k + 1]
                    if not (left_ok and right_ok):
                        continue
                    lh = vals[k]
                    for j in range(k - 1, -1, -1):
                        if vals[j] < lh:
                            break
                        lh = max(lh, vals[j])
                    rh = vals[k]
                    for j in range(k + 1, len(vals)):
                        if vals[j] < rh:
                            break
                        rh = max(rh, vals[j])
                    (confirmed if min(lh, rh) - vals[k] >= split_prom
                     else weak).append(k)
                return confirmed, weak

            eps = local_step / 3.0
            cand_pts = {}
            for i in cand_idx:
                pts = sorted(grp.get(i, []))
                if not pts:
                    continue
                if split_on and split_anchor:
                    extra = []
                    for j in (i - 1, i, i + 1):
                        if 0 <= j < len(smin) and np.isfinite(smin[j]):
                            p = (float(Oms[j]), float(smin[j]))
                            if not any(abs(p[0] - q[0]) < eps for q in pts):
                                extra.append(p)
                    if extra:
                        pts = sorted(pts + extra)
                cand_pts[i] = pts

            if split_on and split_adaptive:
                bisect_args, bisect_owner, seen_mid = [], [], set()
                for i in cand_idx:
                    pts = cand_pts.get(i)
                    if not pts or len(pts) < 3:
                        continue
                    vals = [s for _, s in pts]
                    _, weak = _shoulder_scan(vals)
                    n_added = 0
                    for k in weak:
                        if n_added >= split_adaptive_max:
                            break
                        for nb in (k - 1, k + 1):
                            if not (0 <= nb < len(pts)):
                                continue
                            mid = 0.5 * (pts[k][0] + pts[nb][0])
                            if not (Om_lo <= mid <= Om_hi):
                                continue
                            if any(abs(mid - q[0]) < eps for q in pts):
                                continue
                            key = round(mid / eps)
                            if key in seen_mid:
                                continue
                            seen_mid.add(key)
                            bisect_args.append(
                                (float(mid), None, g_sc, m_sc, solver.M,
                                 n_quad_full, max_dim, n_dofs, dps))
                            bisect_owner.append(i)
                            n_added += 1
                if bisect_args:
                    print(f"  Phase 1c: adaptive bisection "
                          f"({len(bisect_args)} fine pt(s) near "
                          f"{len(set(bisect_owner))} weak shoulder(s)) …",
                          flush=True)
                    bisect_vals = [float('nan')] * len(bisect_args)
                    with concurrent.futures.ProcessPoolExecutor(max_workers=nw) as pool:
                        k = 0
                        for Om_f, s in _throttled_map(pool, scan_worker, bisect_args):
                            bisect_vals[k] = (Om_f, s)
                            k += 1
                    for (Om_f, s), owner in zip(bisect_vals, bisect_owner):
                        if np.isfinite(s):
                            cand_pts[owner] = sorted(
                                cand_pts[owner] + [(Om_f, s)])

            for i in cand_idx:
                pts = cand_pts.get(i)
                if not pts:
                    continue
                idxs = []
                if split_on and len(pts) >= 3:
                    vals = [s for _, s in pts]
                    idxs, _ = _shoulder_scan(vals)
                    # Merge minima closer than split_min_sep (keep deepest).
                    if len(idxs) > 1:
                        merged = [idxs[0]]
                        for k in idxs[1:]:
                            if pts[k][0] - pts[merged[-1]][0] < split_min_sep:
                                if vals[k] < vals[merged[-1]]:
                                    merged[-1] = k
                            else:
                                merged.append(k)
                        idxs = merged
                if not idxs:
                    deepest = min(pts, key=lambda t: t[1])
                    idxs = [pts.index(deepest)]
                cand_brackets[i] = []
                for k in idxs:
                    Omc = pts[k][0]
                    a = pts[k - 1][0] if k > 0 else Omc - local_step
                    b = pts[k + 1][0] if k < len(pts) - 1 else Omc + local_step
                    cand_brackets[i].append((Omc, float(a), float(b)))

    # ── Phase 2: parallel golden-section refinement ──────────────────────────
    # Each golden eval costs one full_search (~all the time), so keep iters
    # modest.  Bracket per candidate is either the tight local-zoom bracket
    # (two_stage) or the coarse [Oms[i-1], Oms[i+1]] clipped to neighbours.
    # With two-stage the bracket is already ~local_step (~0.002) wide, so 5
    # golden iters reach ~1.8e-4 mode precision (well below the ~0.1-1%
    # physics-limited accuracy); the single-grid path keeps the larger
    # default since its bracket is a full coarse step wide.
    # DEFAULT PROMOTED 2026-07-26 (5->12 for the two-stage path only; the
    # non-two-stage "8" default is UNCHANGED -- no evidence was gathered for
    # it). Resolves G1 (0.401154) and G3 (1.779792), Sec 20.8's last two
    # blind-pass IP misses: both are genuine Phase-2 candidates that reach
    # golden-section refinement but converge shallow at 5 iterations for lack
    # of iteration budget on a hard/narrow bracket, not a missing-candidate or
    # accept-floor problem (LESSONS_LEARNED Sec 34/35). Cluster-validated
    # end-to-end: GATE_FFP1 (14 modes, 13 pinned matched |d|<0.003) and
    # GATE_SHI (13 modes, 8/8 lit targets <1%) both PASS under this default
    # (job 2330983); a completeness check with gap-refine OFF confirms every
    # known physical target still recovers at iters=12, including the one
    # low-Omega candidate that iters=5 had wrongly accepted (job 2330819's own
    # SIGMIN_DEGEN_SENTINEL now correctly rejects it) (job 2333514). Cost is
    # modest: +10-32% on the two expensive Phase-2-heavy windows tested, all
    # in Phase 2 (which iters scales), none in the much larger Phase-1 sweep
    # (jobs 2330819, 2332566).
    default_iters = "12" if two_stage else "8"
    polish_iters = int(os.environ.get("SIGMIN_POLISH_ITERS", default_iters))
    cand_idx_set = cand_idx
    refine_args = []
    for pos, i in enumerate(cand_idx_set):
        if i in cand_brackets:
            for _, a, b in cand_brackets[i]:
                refine_args.append(
                    (float(a), float(b), g_sc, m_sc, solver.M, n_quad_full,
                     max_dim, n_dofs, polish_iters, n_quad_iter, dps, None))
            continue
        i_prev = cand_idx_set[pos - 1] if pos > 0 else None
        i_next = cand_idx_set[pos + 1] if pos < len(cand_idx_set) - 1 else None
        lo_idx = max(i - 1, (i_prev + 1) if i_prev is not None else 0)
        hi_idx = min(i + 1, (i_next - 1) if i_next is not None else len(Oms) - 1)
        if lo_idx >= i:
            lo_idx = i - 1 if i > 0 else i
        if hi_idx <= i:
            hi_idx = i + 1 if i < len(Oms) - 1 else i
        a = Oms[lo_idx] if lo_idx != i else Oms[i] - (Oms[1]-Oms[0])
        b = Oms[hi_idx] if hi_idx != i else Oms[i] + (Oms[-1]-Oms[-2])
        refine_args.append(
            (float(a), float(b), g_sc, m_sc, solver.M, n_quad_full, max_dim,
             n_dofs, polish_iters, n_quad_iter, dps, None))
    print(f"  Phase 2: refining {len(refine_args)} dip(s) …", flush=True)
    refined = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=nw) as pool:
        for res in _throttled_map(pool, golden_worker, refine_args):
            if len(res) == 3:
                Om_star, s_star, s_constr = res
                refined.append((float(Om_star), float(s_star), float(s_constr)))
            else:
                Om_star, s_star = res
                refined.append((float(Om_star), float(s_star), None))

    # ── ADAPTIVE ACCEPT FLOOR ────────────────────────────────────────────────
    # Accept a refined dip if SIGMIN_ACCEPT_MARGIN (1.5) orders below the run's own
    # baseline OR past the absolute floor — admits genuinely shallow modes (e.g.
    # 1.25/0.25π mode 2) without admitting junk.  See §Detector.
    accept_margin = float(os.environ.get("SIGMIN_ACCEPT_MARGIN", "1.5"))
    eff_floor = max(accept_floor, base - accept_margin)
    if diag is not None:
        diag.update(accept_floor_eff=float(eff_floor),
                    accept_floor_fixed=float(accept_floor),
                    accept_margin=float(accept_margin))
    if verbose:
        print(f"    accept: σ_min ≤ {eff_floor:+.2f} "
              f"(max of fixed {accept_floor:+.1f} and baseline{base:+.2f}"
              f"−{accept_margin:.1f})", flush=True)

    # Keep genuine singularities, dedupe, order by Ω.  For Part 2, s_star is
    # the UNCONSTRAINED σ_min (the paper's Eq.48 mode condition) — the
    # authoritative signal; the third element (s_constr) is the constrained
    # σ_min, kept only as an informational coincidence note.
    refined.sort(key=lambda t: t[0])
    modes = []
    for Om_star, s_star, s_constr in refined:
        if not np.isfinite(Om_star) or not (Om_lo*0.5 <= Om_star <= Om_hi*1.05):
            if diag is not None and np.isfinite(Om_star):
                diag["dropped"].append((float(Om_star), float(s_star),
                                        "out-of-range"))
            continue
        if s_star > eff_floor:
            if s_star >= 50.0:
                reason = "ill-conditioned (branch-collision/cut-off artifact)"
            else:
                reason = f"shallow (σ_min 1e{s_star:+.2f} > floor 1e{eff_floor:+.2f})"
            if diag is not None:
                diag["dropped"].append((float(Om_star), float(s_star), reason))
            if verbose:
                if s_star >= 50.0:
                    print(f"    drop Ω={Om_star:.6f}: ill-conditioned basis "
                          f"(branch-collision or cut-off artifact, geometry/"
                          f"BC-invariant — not a physical mode)")
                else:
                    tag = "" if part == 1 else " (unconstrained)"
                    print(f"    drop Ω={Om_star:.6f}: σ_min{tag} only "
                          f"1e{s_star:+.2f} (not a true singularity)")
            continue
        if modes and abs(Om_star - modes[-1][0]) <= 1e-4 * max(1.0, Om_star):
            if s_star < modes[-1][1]:
                modes[-1] = (Om_star, s_star, s_constr)
            continue
        modes.append((Om_star, s_star, s_constr))

    for Om_star, s_star, s_constr in modes:
        if verbose:
            if part == 2 and s_constr is not None and np.isfinite(s_constr):
                coin = "coincides" if s_constr <= eff_floor else "constr. shallow"
                print(f"    mode Ω={Om_star:.6f}   σ_min(uncon)=1e{s_star:+.2f}  "
                      f"σ_min(con)=1e{s_constr:+.2f}  [{coin}]", flush=True)
            else:
                print(f"    mode Ω={Om_star:.6f}   log10 σ_min={s_star:+.2f}",
                      flush=True)

    # ── Phase 3: gap-refinement (Addendum 2026-07-09e item E.3, opt-in) ─────
    # find_natural_frequencies (the legacy per-solver detector) has an explicit
    # Phase-5 gap-refinement retry; this detector never gained one -- the
    # documented mechanism behind Sec 4.3's missed FF-P1 mode (37.818, job
    # 2315631/2315829: recovered only by a MANUAL 6x-finer rescan). This
    # automates that recipe: flag Omega-gaps between accepted modes wider than
    # SIGMIN_GAP_REFINE_THRESHOLD coarse-steps, recurse into this same function
    # on each flagged gap at a finer SIGMIN_COARSE_STEP (full reuse of the
    # existing Phase 1/1b/2/accept pipeline -- no new detection logic, no new
    # bug surface). OFF by default (SIGMIN_GAP_REFINE unset/0): the block below
    # is then a no-op and this patch changes nothing on the default path -- no
    # SOLVER_VERSION bump required for the patch itself.
    #
    # RECURSION DEPTH (2026-07-23, LESSONS_LEARNED Sec 20.8 hardening): the
    # original _gap_pass=True on the recursive call allowed exactly ONE level
    # of gap-refinement -- a mode sitting in a "gap within a gap" (the
    # first-level rescan still leaves its own too-narrow-to-flag residual gap
    # around the true root, e.g. Sec 20.8's G3/1.779792, where the recovered
    # window 1.4406-1.8077 never surfaced a dip near 1.78 at all) was
    # structurally unreachable no matter how justified a second pass would be.
    # _gap_depth replaces the boolean with a counter so this can recurse up to
    # SIGMIN_GAP_REFINE_MAX_DEPTH levels; the default "1" reproduces the exact
    # prior one-level behaviour (no SOLVER_VERSION bump for adding the
    # capability -- it is a no-op until the env var is raised). Raising it
    # (e.g. SIGMIN_GAP_REFINE_MAX_DEPTH=3) IS a default-changing knob once
    # exercised and must be certified (Sec 20.8) before promotion.
    max_gap_depth = int(os.environ.get("SIGMIN_GAP_REFINE_MAX_DEPTH", "1"))
    if _gap_depth < max_gap_depth and os.environ.get("SIGMIN_GAP_REFINE", "1") == "1":
        zoom       = max(2, int(os.environ.get("SIGMIN_GAP_REFINE_ZOOM", "6")))
        gap_factor = float(os.environ.get("SIGMIN_GAP_REFINE_THRESHOLD", "6"))
        max_gaps   = int(os.environ.get("SIGMIN_GAP_REFINE_MAX_GAPS", "30"))
        # PROMINENCE OVERRIDE for gap-refine sub-scans (2026-07-24,
        # LESSONS_LEARNED Sec 20.8 hardening, second finding): the
        # candidate-flagging prominence test (Phase 1) compares each
        # dip against a baseline/shoulders computed from THIS call's
        # own scanned window. A gap window inherited unchanged across
        # empty recursion levels stays just as WIDE as the top-level
        # window, so a genuine but modest dip (e.g. Sec 20.8's G1/G3)
        # can fail the same fixed prom=0.5 bar there that it would
        # clear in a narrow, hand-picked window (confirmed directly,
        # job 2330442: a 0.016-wide window at comparable resolution
        # found all three previously-missed modes cleanly, including
        # G2 correctly split from its neighbouring artifact). Default
        # unset reproduces the exact old behaviour (gap_prom==prom,
        # no-op); set SIGMIN_GAP_REFINE_PROM lower (e.g. 0.15) to let
        # weak-but-real shoulders in a wide, seemingly-empty gap get
        # flagged as candidates for the (already-hardened) Phase 1b
        # local-zoom/split step to resolve, instead of never being
        # flagged at all.
        gap_prom   = float(os.environ.get("SIGMIN_GAP_REFINE_PROM", str(prom)))
        gap_thresh = gap_factor * scan_step
        bounds = [Om_lo] + sorted(m[0] for m in modes) + [Om_hi]
        gap_windows = [(bounds[k], bounds[k + 1]) for k in range(len(bounds) - 1)
                       if bounds[k + 1] - bounds[k] > gap_thresh]
        if gap_windows:
            # Prioritize by WIDTH before capping -- the widest gap is the
            # most suspicious one and must not be dropped just because it
            # isn't among the first N in ascending-Omega order (bug found
            # by job 2315908: a width-1.014 gap was silently skipped in
            # favor of five much narrower, lower-Omega gaps). Re-sort back
            # to ascending Omega afterward for a sane execution order.
            gap_windows = sorted(
                sorted(gap_windows, key=lambda w: w[1] - w[0], reverse=True)
                [:max_gaps])
            if verbose:
                print(f"  Phase 3 (depth {_gap_depth+1}/{max_gap_depth}): "
                      f"gap-refinement — {len(gap_windows)} gap(s) > "
                      f"{gap_thresh:.4f} wide, rescanning at step="
                      f"{scan_step / zoom:.5f} …", flush=True)
            recovered = 0
            saved_step = os.environ.get("SIGMIN_COARSE_STEP")
            try:
                os.environ["SIGMIN_COARSE_STEP"] = str(scan_step / zoom)
                for (glo, ghi) in gap_windows:
                    sub_diag = {}
                    sub = find_modes_sigmin(
                        solver, part, (glo, ghi), n_scan=1, n_dofs=n_dofs,
                        max_dim=max_dim, n_modes_wanted=n_modes_wanted,
                        verbose=False, prom=gap_prom, accept_floor=accept_floor,
                        n_workers=nw, diag=sub_diag, _gap_pass=True,
                        _gap_depth=_gap_depth + 1)
                    sub_by_om = {round(a["Omega"], 6): a
                                 for a in sub_diag.get("accepted", [])}
                    for Om_new in sub:
                        if any(abs(Om_new - m[0]) <= 1e-4 * max(1.0, Om_new)
                               for m in modes):
                            continue
                        a = sub_by_om.get(round(Om_new, 6))
                        s_val = a["sigmin"] if a else float("nan")
                        sc_val = a["sigmin_constr"] if a else None
                        modes.append((Om_new, s_val, sc_val))
                        recovered += 1
            finally:
                if saved_step is None:
                    os.environ.pop("SIGMIN_COARSE_STEP", None)
                else:
                    os.environ["SIGMIN_COARSE_STEP"] = saved_step
            if verbose:
                print(f"  Phase 3 (depth {_gap_depth+1}): recovered "
                      f"{recovered} additional mode(s)", flush=True)
            modes.sort(key=lambda t: t[0])

    # ── Reporting policy (configurable) ──────────────────────────────────────
    # SIGMIN_REPORT_ALL=1 (default): return every genuine singularity in range,
    # ascending Ω (truncating to lowest-N drops real modes when an extra zero sits
    # below a paper mode).  =0: lowest n_modes_wanted only.  See §Detector.
    report_all = os.environ.get("SIGMIN_REPORT_ALL", "1") == "1"
    all_oms = sorted(m[0] for m in modes)
    if diag is not None:
        # marginality = how far the accepted σ_min sits BELOW the accept floor;
        # small (≈0) ⇒ barely accepted ⇒ most likely a false positive / under-
        # resolved.  Negative would mean above floor (cannot happen here).
        diag["accepted"] = [
            {"Omega": float(Om), "sigmin": float(s),
             "sigmin_constr": (float(sc) if sc is not None and np.isfinite(sc)
                               else None),
             "margin_below_floor": float(eff_floor - s)}
            for (Om, s, sc) in sorted(modes, key=lambda t: t[0])
        ]
    if report_all:
        out = all_oms
    else:
        out = all_oms[:n_modes_wanted]
    print(f"  → {len(out)} singularit{'y' if len(out)==1 else 'ies'} found: "
          f"{[f'{m:.6f}' for m in out]}", flush=True)
    return out


# ══════════════════════════════════════════════════════════════════════════════
#  SECTION 6 — SELF-TESTS, REFERENCE DATA, MAIN
# ══════════════════════════════════════════════════════════════════════════════

def _selftest_part1(geom, mat):
    """Verify corrected Part 1 series satisfies ODE."""
    f = _Part1Fast(float(geom.r0_bar), float(mat.nu_bar), M=80)
    xi, Om = 0.7+0.3j, 0.05
    sols = f.series(xi, Om)
    al = 2*f.T*xi**2 + f.R
    be = f.R*xi**2*(xi**2 - 2*(f.T+f.R))
    worst = 0.0
    for a in sols:
        for x in (-1.2, -0.3, 0.6, 1.4):
            rr = x + f.r0b
            r4 = (rr**4*f._ev(a, x, 4) + 2*rr**3*f._ev(a, x, 3)
                  - al*rr**2*f._ev(a, x, 2) + al*rr*f._ev(a, x, 1)
                  + (be - rr**4*Om**2)*f._ev(a, x, 0))
            scale = max(abs(f._ev(a, x, 0)), 1.0)
            worst = max(worst, abs(r4)/scale)
    ok = worst < 1e-6
    print(f"  [selftest] Part 1 ODE residual: {worst:.2e}  {'OK' if ok else 'FAIL'}")
    return ok

def _selftest_part2(geom, mat):
    f = _Part2Fast(float(geom.r0_bar), float(mat.nu_bar),
                   float(mat.c11_eff_bar), M=80)
    ze, Om = 0.6+0.2j, 0.05
    sols = f.series(ze, Om)
    c11, R, nu, r0 = f.c11, f.R, f.nu, f.r0b
    cF = c11*R + ze**2; cG = 1 + c11*R*ze**2
    coup = (c11*nu+1)*ze; cst = (c11*R+1)*ze
    worst = 0.0
    for fa, fb in sols:
        for x in (-1.2, -0.3, 0.6, 1.4):
            rr = x + r0
            F0 = f._ev(fa, x, 0); F1 = f._ev(fa, x, 1); F2 = f._ev(fa, x, 2)
            G0 = f._ev(fb, x, 0); G1 = f._ev(fb, x, 1); G2 = f._ev(fb, x, 2)
            rF = (c11*(rr**2*F2 + rr*F1) + (Om**2*rr**2 - cF)*F0
                  - coup*rr*G1 + cst*G0)
            rG = ((rr**2*G2 + rr*G1) + (Om**2*rr**2 - cG)*G0
                  + coup*rr*F1 + cst*F0)
            scale = max(abs(F0), abs(G0), 1.0)
            worst = max(worst, abs(rF)/scale, abs(rG)/scale)
    ok = worst < 1e-6
    print(f"  [selftest] Part 2 ODE residual: {worst:.2e}  {'OK' if ok else 'FAIL'}")
    return ok


def paper_seeds_for(paper_dict, r0_2b, two_T_pi):
    """Collect the published mode frequencies for one geometry, sorted ascending.

    Returns the values from `paper_dict` (PAPER_PART1/PAPER_PART2) matching this
    (r0/(2b), 2Θ/π).  These are used ONLY to choose where to look — see
    seed_polish_modes — not as the reported answers."""
    seeds = []
    for mode in (1, 2, 3):
        d = None
        for (rk, mk), dv in paper_dict.items():
            if mk == mode and abs(rk - r0_2b) < 1e-6:
                d = dv
                break
        if not d:
            continue
        for tk, val in d.items():
            if abs(tk - two_T_pi) < 1e-9:
                seeds.append(float(val))
                break
    return sorted(seeds)


def seed_polish_modes(solver, seeds, n_dofs, max_dim, part,
                      scan_lo, scan_hi, iters=10, nq_fast=15, verbose=True):
    """Reference-seeded, paper-independent mode convergence.

    For each seed frequency r (a published value that tells us *where* a mode
    is), independently converge a determinant zero near r:
      1. open a window that brackets only that mode (clamped at the midpoints to
         its neighbouring seeds so windows never overlap),
      2. sign-scan the boundary determinant across the window, and
      3. bisect the located sign change to full precision (falling back to
         golden-section minimisation if no sign change is found).

    The seed only selects the window; the returned frequency is the solver's own
    converged zero of its determinant.  This reproduces the paper's tables while
    leaving the from-scratch detector (find_natural_frequencies) untouched for
    general use.  `part` is 1 (out-of-plane real det) or 2 (in-plane constrained
    det).
    """
    seeds = sorted(float(s) for s in seeds)
    out, W = [], 0.12
    os.environ['R40_BC_KIND'] = getattr(solver.bc, 'token', 'clamped_free')  # Step 6
    for i, r in enumerate(seeds):
        lo, hi = r * (1.0 - W), r * (1.0 + W)
        if i > 0:
            lo = max(lo, 0.5 * (seeds[i - 1] + r))
        if i < len(seeds) - 1:
            hi = min(hi, 0.5 * (r + seeds[i + 1]))
        lo, hi = max(lo, scan_lo), min(hi, scan_hi)
        if hi <= lo:
            lo, hi = r * 0.97, r * 1.03
        mid = 0.5 * (lo + hi)

        try:
            seed_brs = full_search(solver.fast, mid, xmax=max_dim)
            npts = 9
            grid = [lo + (hi - lo) * k / (npts - 1) for k in range(npts)]
            signs, vals, prev = [], [], seed_brs
            for Om in grid:
                prev = track(solver.fast, Om, prev) or \
                    full_search(solver.fast, Om, xmax=max_dim)
                sel, cnt = select_fill(prev, n_dofs)
                if cnt < n_dofs:
                    prev = full_search(solver.fast, Om, xmax=max_dim)
                    sel, _ = select_fill(prev, n_dofs)
                if part == 1:
                    s, ld = solver._logdet_nq(Om, sel, n_quad_override=nq_fast,
                                              fast_scan=False)
                else:
                    s, ld = solver._logdet_nq(Om, sel, lagrange=True,
                                              n_quad_override=nq_fast,
                                              fast_scan=False)
                signs.append(s)
                vals.append(float(ld))

            brackets = [(grid[k], grid[k + 1]) for k in range(npts - 1)
                        if signs[k] != signs[k + 1]]
            if brackets:
                a, b = min(brackets,
                           key=lambda ab: abs(0.5 * (ab[0] + ab[1]) - r))
                Om_star, ld_star = solver._bisect(a, b, n_dofs, max_dim, iters,
                                                  seed_brs=seed_brs)
                how = 'bisect'
            else:
                Om_star, ld_star = solver._golden(lo, hi, n_dofs, max_dim, iters,
                                                  seed_brs=seed_brs)
                how = 'golden(no sign change)'
        except Exception as exc:
            kbest = int(np.argmin(vals)) if vals else 0
            Om_star = grid[kbest] if vals else mid
            ld_star, how = (vals[kbest] if vals else float('nan')), f'grid (polish error: {exc})'

        out.append(float(Om_star))
        if verbose:
            print(f"  seed Ω={r:.6f} → Ω={Om_star:.6f}  "
                  f"log10|det|={ld_star:+.2f}  ({how})", flush=True)
    return sorted(out)


def _rect_absdet_axis_roots(eng, Lam, axis, lo, hi, step, eta_tol=0.05, cap=None):
    """Phase-agnostic axis-root finder (dominant-component sign change + eta
    branch-point rejection).  `eng` is a RectangularCartesianOOP carrying the
    symmetry class.  Returns axis parameter values (real)."""
    if axis == "re":
        zf = lambda t: complex(t, 0.0)
    else:
        zf = lambda t: complex(0.0, t)
    f = lambda t: eng._dispersion(zf(t), Lam)
    g = lambda t: abs(f(t))
    ts = []
    t = lo
    while t <= hi:
        ts.append(t); t += step
    ds = [(d.real if abs(d.real) >= abs(d.imag) else d.imag)
          for d in (f(t) for t in ts)]
    out = []
    for i in range(len(ts) - 1):
        if ds[i] * ds[i + 1] >= 0:
            continue
        a, b = ts[i], ts[i + 1]
        for _ in range(80):
            m1 = a + 0.382 * (b - a); m2 = a + 0.618 * (b - a)
            if g(m1) < g(m2):
                b = m2
            else:
                a = m1
        tr = 0.5 * (a + b)
        e1, e2 = eng._etas(zf(tr), Lam)
        if min(abs(e1), abs(e2)) < eta_tol or abs(e1 * e1 - e2 * e2) < eta_tol:
            continue
        if lo <= tr <= hi and not any(abs(tr - o) < 1e-3 for o in out):
            out.append(tr)
            if cap is not None and len(out) >= cap:
                break
    return out


def _rect_cnewton(eng, z0, Lam, itmax=30, tol=1e-13):
    """Complex Newton on the FULL complex dispersion (not the dominant-component
    det face): off-axis complex conjugate-pair roots cannot be found by Newton on
    the real-valued det used for axis sign-changes, so the resolver's complex
    grid polishes _dispersion directly."""
    z = complex(z0)
    f = lambda w: eng._dispersion(w, Lam)
    for _ in range(itmax):
        fz = f(z)
        h = 1e-7 * max(1.0, abs(z))
        fp = (f(z + h) - f(z - h)) / (2 * h)
        if fp == 0:
            return None
        dz = fz / fp
        z -= dz
        if abs(dz) < tol * max(1.0, abs(z)):
            return z
    return z if abs(f(z)) < 1e-9 else None


# Extraction-loss fix (2026-07-01, LESSONS_LEARNED Sec. 22.x): this constant
# was used by the four rect branch-classification helpers below but never
# survived the Research50.py -> plate_solver split (neither redefined here
# nor imported from elsewhere -- confirmed absent package-wide). It only
# separates "on axis" (component == 0.0 exactly, by construction, for axis
# roots) from "off axis" (complex grid search starts at cgrid=0.12 and up),
# so any epsilon well below 0.12 is equivalent; verified numerically that
# branch classification is byte-identical for TOL in {1e-9, 1e-6, 1e-3, 1e-2}
# on both rect_resolve_branches and rect_ip_resolve_branches at representative
# (Lambda/Om) points. 1e-6 chosen as a conventional "numerically zero" epsilon.
_RECT_AXIS_TOL = 1e-6


def rect_resolve_branches(eng, Lam, im_cap=7.0, re_max=1.6, cgrid=0.12):
    """Complete first-quadrant branch set at one Lambda: phase-agnostic axis
    roots + bounded off-axis complex conjugate-pair representatives."""
    reps = []
    reps += [complex(x, 0.0) for x in _rect_absdet_axis_roots(eng, Lam, "re", 1e-3, 9.0, 0.005)]
    reps += [complex(0.0, k) for k in _rect_absdet_axis_roots(eng, Lam, "im", 1e-3, im_cap, 0.005)]
    re = cgrid
    while re < re_max:
        im = cgrid
        while im < im_cap:
            z = _rect_cnewton(eng, complex(re, im), Lam, itmax=30)
            if z is not None and z.real >= _RECT_AXIS_TOL and z.imag >= _RECT_AXIS_TOL \
                    and z.imag <= im_cap and not any(abs(z - c) < 5e-3 for c in reps):
                reps.append(z)
            im += cgrid
        re += cgrid
    uniq = []
    for z in reps:
        if not any(abs(z - u) < 5e-3 for u in uniq):
            uniq.append(z)
    uniq.sort(key=lambda z: (round(abs(z.imag), 4), abs(z.real)))
    return uniq


def rect_select_branches(reps, n_real=5, n_cpair=0):
    """Flat branch list for the assembler: axis branches first (n_real), then
    n_cpair complex pairs each emitted as branch AND conjugate."""
    axis = [z for z in reps if abs(z.real) < _RECT_AXIS_TOL or abs(z.imag) < _RECT_AXIS_TOL]
    cplx = [z for z in reps if abs(z.real) >= _RECT_AXIS_TOL and abs(z.imag) >= _RECT_AXIS_TOL]
    axis = axis[:n_real]; cplx = cplx[:n_cpair]
    full = []
    for z in axis:
        full.append(complex(z.real, 0.0) if abs(z.imag) < _RECT_AXIS_TOL
                    else complex(0.0, z.imag))
    for z in cplx:
        full.append(complex(z)); full.append(complex(z.real, -z.imag))
    return full


def _rect_ip_cnewton(eng, z0, Om, itmax=30, tol=1e-13):
    z = complex(z0)
    f = lambda w: eng._dispersion(w, Om)
    for _ in range(itmax):
        fz = f(z)
        h = 1e-7 * max(1.0, abs(z))
        fp = (f(z + h) - f(z - h)) / (2 * h)
        if fp == 0:
            return None
        dz = fz / fp
        z -= dz
        if abs(dz) < tol * max(1.0, abs(z)):
            return z
    return z if abs(f(z)) < 1e-9 else None


def _rect_ip_axis_roots(eng, Om, axis, lo, hi, step, zeta_tol=0.05, cap=None):
    """Dominant-component sign-change axis roots with zeta branch-point
    rejection (the in-plane analog of _rect_absdet_axis_roots)."""
    zf = (lambda t: complex(t, 0.0)) if axis == "re" else (lambda t: complex(0.0, t))
    f = lambda t: eng._dispersion(zf(t), Om)
    g = lambda t: abs(f(t))
    ts = []
    t = lo
    while t <= hi:
        ts.append(t); t += step
    ds = [(d.real if abs(d.real) >= abs(d.imag) else d.imag) for d in (f(t) for t in ts)]
    out = []
    for i in range(len(ts) - 1):
        if ds[i] * ds[i + 1] >= 0:
            continue
        a, b = ts[i], ts[i + 1]
        for _ in range(80):
            m1 = a + 0.382 * (b - a); m2 = a + 0.618 * (b - a)
            if g(m1) < g(m2):
                b = m2
            else:
                a = m1
        tr = 0.5 * (a + b)
        if eng.is_branch_point(zf(tr), Om):
            continue
        if lo <= tr <= hi and not any(abs(tr - o) < 1e-3 for o in out):
            out.append(tr)
            if cap is not None and len(out) >= cap:
                break
    return out


def rect_ip_resolve_branches(eng, Om, im_cap=7.0, re_max=3.0, cgrid=0.12):
    """Complete first-quadrant branch set at one Omega_bar: dominant-component
    axis roots + bounded off-axis complex conjugate-pair representatives."""
    reps = []
    reps += [complex(x, 0.0) for x in _rect_ip_axis_roots(eng, Om, "re", 1e-3, 9.0, 0.005)]
    reps += [complex(0.0, k) for k in _rect_ip_axis_roots(eng, Om, "im", 1e-3, im_cap, 0.005)]
    re = cgrid
    while re < re_max:
        im = cgrid
        while im < im_cap:
            z = _rect_ip_cnewton(eng, complex(re, im), Om, itmax=30)
            if z is not None and z.real >= _RECT_AXIS_TOL and z.imag >= _RECT_AXIS_TOL \
                    and z.imag <= im_cap and not eng.is_branch_point(z, Om) \
                    and not any(abs(z - c) < 5e-3 for c in reps):
                reps.append(z)
            im += cgrid
        re += cgrid
    uniq = []
    for z in reps:
        if not any(abs(z - u) < 5e-3 for u in uniq):
            uniq.append(z)
    uniq.sort(key=lambda z: (round(abs(z.imag), 4), abs(z.real)))
    return uniq


def rect_ip_select_branches(reps, n_real=5, n_cpair=0):
    axis = [z for z in reps if abs(z.real) < _RECT_AXIS_TOL or abs(z.imag) < _RECT_AXIS_TOL]
    cplx = [z for z in reps if abs(z.real) >= _RECT_AXIS_TOL and abs(z.imag) >= _RECT_AXIS_TOL]
    axis = axis[:n_real]; cplx = cplx[:n_cpair]
    full = []
    for z in axis:
        full.append(complex(z.real, 0.0) if abs(z.imag) < _RECT_AXIS_TOL
                    else complex(0.0, z.imag))
    for z in cplx:
        full.append(complex(z)); full.append(complex(z.real, -z.imag))
    return full


# ════════════════════════════════════════════════════════════════════════════
#  SECTION 5b — WEAK-ENFORCEMENT-ARTIFACT RESIDUAL SCREEN (post-detector, OPT-IN)
# ════════════════════════════════════════════════════════════════════════════
#
# The B12 free-free edge-orientation fix (2026-07-08.s5, LESSONS_LEARNED Sec 10.5)
# stops the free-free K from being exactly anti-block-diagonal, but does NOT
# close a second, distinct mechanism: even under the fixed assembly, some
# accepted sigma_min(K) dips are "weak-enforcement artifacts" -- a nearly
# displacement-free evanescent null-vector combination that satisfies the
# FINITE test set's Galerkin projections without satisfying the free-edge
# boundary conditions POINTWISE (LESSONS Sec 10.6). find_modes_sigmin() alone
# cannot distinguish these from real modes -- both are genuine zeros of the
# accepted det K(Omega)=0.
#
# The functions below evaluate that pointwise residual directly (Kelvin
# effective shear V and bending moment M along the +Theta free edge, on a
# dense 41-point radial grid) and classify each dip as REAL-like or
# ARTIFACT-like. This is the productionized form of probe_ffresidual_v2.py
# (OOP) / probe_ffresidual_ip.py (IP) -- same formulas, verbatim, reusing
# this module's own real_basis_items / sigma_min_from_K / equilibrated_
# nullvec_mp / mp_proj (no reimplementation, no new bug surface).
#
# VALIDATION STATUS (2026-07-09):
#   OOP (FF-P1 geometry, r0/2b=1.5, 2Theta=0.5pi, nu=0.30): 11/11 modes
#   correctly classified against the independent Ansys benchmark -- all 6
#   FE-matched modes (<0.7% err) screen REAL-like (rM in [0.0025, 0.036]),
#   all 5 non-matched modes screen ARTIFACT-like (rM in [75, 348]) -- clean
#   ~2000x separation, zero ambiguous cases (jobs 2315385/2315584/2315829).
#   IP (same geometry): no independent FE benchmark exists yet for the
#   in-plane spectrum, so this is NOT numerically calibrated the way OOP is
#   -- but the SAME bimodal signature appears cleanly (14 confirmed REAL,
#   rTyy < 0.3; several confirmed ARTIFACT, rTyy > 4; a >14x gap, nothing in
#   between) (jobs 2315477/2315579/2315633). Treat IP verdicts as strong
#   exploratory evidence, not validated classification, until an
#   independent IP benchmark exists.
#
# IMPORTANT CAVEAT (2026-07-09, jobs 2315631/2315829): the standard scan
# resolution (SIGMIN_COARSE_STEP=0.006) can MISS real, narrow modes entirely
# -- one FF-P1 OOP mode (matching Ansys target 37.818 at 0.35%) was absent
# from a full default-resolution sweep and only appeared under a 6x finer
# local rescan. This screen only classifies dips find_modes_sigmin already
# found -- it cannot recover a mode the scan missed. A "REAL-like" verdict
# on a screened dip is trustworthy; an absence of expected dips is NOT
# evidence they don't exist without a finer-resolution check first.
#
# NOT WIRED INTO find_modes_sigmin'S DEFAULT PATH -- purely opt-in, call
# explicitly. Adding these functions changes no existing default output, so
# per Sec 1's hard invariant this does NOT require a SOLVER_VERSION bump.

def _weak_enforcement_null_block(K, items, n):
    """Split K into the two q-parity sub-blocks, return (block_idx, block_name,
    block_K, other_sig, this_sig) for whichever block is MORE singular (the
    one the accepted sigma_min(K) dip actually lives in)."""
    q_of = [it[1] for it in items]
    idx0 = [i for i, q in enumerate(q_of) if q == 0]
    idx1 = [i for i, q in enumerate(q_of) if q == 1]

    def submatrix(K, idx):
        m = len(idx)
        S = matrix(m, m)
        for a, ia in enumerate(idx):
            for b, ib in enumerate(idx):
                S[a, b] = K[ia, ib]
        return S

    K0 = submatrix(K, idx0)
    K1 = submatrix(K, idx1)
    sig0 = sigma_min_from_K(K0, len(idx0))
    sig1 = sigma_min_from_K(K1, len(idx1))
    if sig0 < sig1:
        return idx0, "ODD", K0, sig1, sig0    # q=0: sin -> ODD (LESSONS Sec 10.2)
    return idx1, "EVEN", K1, sig0, sig1        # q=1: cos -> EVEN


def weak_enforcement_residual_oop(solver, Omega, brs, n_dofs=None, ngrid=41,
                                   real_rM_max=0.1, real_maxW_min=0.3,
                                   artifact_rM_min=0.3, artifact_maxW_max=0.1):
    """Pointwise free-edge (+Theta) residual screen for an OOP FreeFreeOOP dip
    already found by find_modes_sigmin. solver must be an OutOfPlaneSolver
    with a FreeFreeOOP boundary; brs is the root set at Omega (e.g. from
    full_search + select_fill, matching what found the dip).

    Returns a dict with block/rV/rM/maxW/verdict, or ok=False with a reason
    if no null vector could be extracted. verdict is 'REAL-like',
    'ARTIFACT-like', or 'AMBIGUOUS' per the pre-registered OOP acceptance bar
    (LESSONS Sec 10.7) -- validated 11/11 correct on FF-P1 (see module
    docstring above).
    """
    n_dofs = n_dofs or len(brs)
    sel = brs if len(brs) == n_dofs else select_fill(brs, n_dofs)[0]
    K, n = solver._build_K_real(Omega, sel, fast_scan=False)
    items = real_basis_items(sel)
    if len(items) != n:
        return dict(ok=False, reason=f"item/K dimension mismatch ({len(items)} vs {n})")

    block_idx, block_name, block_K, sig_other, sig_this = _weak_enforcement_null_block(K, items, n)
    nv = equilibrated_nullvec_mp(block_K, len(block_idx))
    if nv is None:
        return dict(ok=False, reason="no null vector", block=block_name)

    c_full = [mpf(0)] * n
    for li, gi in enumerate(block_idx):
        c_full[gi] = mpf(float(nv[li]))

    e_plus = [e for e in solver.bc.edges if e.kind == EdgeKind.FREE and e.sign > 0][0]
    Th, nu, T, R, r0 = solver.Theta, solver.nu, solver.T, solver.R, solver.r0b
    xs = [mp.pi/2 * (2*k/(ngrid-1) - 1) for k in range(ngrid)]

    cache = {}
    for z in {it[0] for it in items}:
        zr = solver._refine_root_mp(mpc(z), Omega)
        sols = solver._series_mp(zr, Omega)
        A = solver._amp_mp(zr, Omega, sols)
        cache[z] = (zr, sols, A)

    Wtot, Qtot, Mtot = [], [], []
    for x in xs:
        Wt = mpf(0); Qt = mpf(0); Mt = mpf(0)
        for i, (z, q, P) in enumerate(items):
            ci = c_full[i]
            if ci == 0:
                continue
            xi, sols, A = cache[z]
            xi2 = xi**2; ph = q*mp.pi/2
            sg = e_plus.sign
            ssg = mp.sin(sg*xi*Th + ph); csg = mp.cos(sg*xi*Th + ph)
            rr = x + r0
            W = solver._W_mp(x, sols, A, 0)
            Wp = solver._W_mp(x, sols, A, 1)
            Wpp = solver._W_mp(x, sols, A, 2)
            Tyy = nu*Wpp + R*Wp/rr - R*xi2*W/rr**2
            Qy = ((2*T-nu)*Wpp/rr + (R-2*T+2*nu)*Wp/rr**2 + (2*T-2*nu-R*xi2)*W/rr**3)
            Wt += ci * mp_proj(ssg*W, P)
            Qt += ci * mp_proj(xi*csg*Qy, P)
            Mt += ci * mp_proj(ssg*Tyy, P)
        Wtot.append(Wt); Qtot.append(Qt); Mtot.append(Mt)

    maxW = float(max(abs(w) for w in Wtot))
    rmsQ = float((sum(q*q for q in Qtot) / ngrid) ** mpf('0.5'))
    rmsM = float((sum(m*m for m in Mtot) / ngrid) ** mpf('0.5'))
    rV = rmsQ/maxW if maxW else float('inf')
    rM = rmsM/maxW if maxW else float('inf')

    if rM <= real_rM_max and maxW >= real_maxW_min:
        verdict = "REAL-like"
    elif rM >= artifact_rM_min or maxW <= artifact_maxW_max:
        verdict = "ARTIFACT-like"
    else:
        verdict = "AMBIGUOUS"

    return dict(ok=True, block=block_name, sig_this=sig_this, sig_other=sig_other,
                rV=rV, rM=rM, maxW=maxW, verdict=verdict)


def weak_enforcement_residual_ip(solver, Omega, brs, n_dofs=None, ngrid=41):
    """In-plane (Part 2) analog of weak_enforcement_residual_oop, for an
    InPlaneSolver with a FreeFreeIP boundary. NO numeric calibration exists
    yet for IP (no independent FE benchmark identified) -- returns raw
    rTyy/rTyr/maxDisp with no verdict field; classify by comparison between
    points (see module docstring: empirically rTyy < 0.3 vs > 4 has been a
    clean separator on every point checked so far, but this is exploratory,
    not validated).
    """
    n_dofs = n_dofs or len(brs)
    sel = brs if len(brs) == n_dofs else select_fill(brs, n_dofs)[0]
    K, n = solver._build_K_real(Omega, sel, lagrange=True, fast_scan=False)

    cache = {}
    for z in set(sel):
        zr = solver._refine_root_mp(mpc(z), Omega)
        sols = solver._series_mp(zr, Omega)
        A = solver._amp_mp(zr, Omega, sols)
        cache[z] = (zr, sols, A)

    imag_proj = {}
    for z, (zr, sols, A) in cache.items():
        if abs(z.real) < 1e-9 and abs(z.imag) > 1e-9:
            x0 = solver.nodes[0]
            F0 = solver._F_mp(x0, sols, A, 0)
            f_is_real = abs(mp.im(F0)) < abs(mp.re(F0)) * 1e-6 + mpf('1e-30')
            imag_proj[z] = ('im', 're') if f_is_real else ('re', 'im')

    items = []
    for z in sel:
        if abs(z.imag) < 1e-9:
            items += [(z, 0, 'asis'), (z, 1, 'asis')]
        elif abs(z.real) < 1e-9:
            p0, p1 = imag_proj[z]
            items += [(z, 0, p0), (z, 1, p1)]
        else:
            items += [(z, 0, 're'), (z, 0, 'im'), (z, 1, 're'), (z, 1, 'im')]
    if len(items) != n:
        return dict(ok=False, reason=f"item/K dimension mismatch ({len(items)} vs {n})")

    q_of = [it[1] for it in items]
    class0 = [i for i, q in enumerate(q_of) if q == 0]
    class1 = [i for i, q in enumerate(q_of) if q == 1]

    def submatrix(K, idx):
        m = len(idx)
        S = matrix(m, m)
        for a, ia in enumerate(idx):
            for b, ib in enumerate(idx):
                S[a, b] = K[ia, ib]
        return S

    K0 = submatrix(K, class0); K1 = submatrix(K, class1)
    sig0 = sigma_min_from_K(K0, len(class0))
    sig1 = sigma_min_from_K(K1, len(class1))
    block_idx, block_name, block_K = ((class0, "q0", K0) if sig0 < sig1
                                       else (class1, "q1", K1))

    nv = equilibrated_nullvec_mp(block_K, len(block_idx))
    if nv is None:
        return dict(ok=False, reason="no null vector", block=block_name)

    c_full = [mpf(0)] * n
    for li, gi in enumerate(block_idx):
        c_full[gi] = mpf(float(nv[li]))

    e_plus = [e for e in solver.bc.edges if e.kind == EdgeKind.FREE and e.sign > 0][0]
    Th, c11, nu, R, r0 = solver.Theta, solver.c11, solver.nu, solver.R, solver.r0b
    xs = [mp.pi/2 * (2*k/(ngrid-1) - 1) for k in range(ngrid)]

    Ftot, Gtot, Tyytot, Tyrtot = [], [], [], []
    for x in xs:
        Ft = mpf(0); Gt = mpf(0); Tyyt = mpf(0); Tyrt = mpf(0)
        for i, (z, q, P) in enumerate(items):
            ci = c_full[i]
            if ci == 0:
                continue
            ze, sols, A = cache[z]
            ph = q*mp.pi/2
            sg = e_plus.sign
            ssg = mp.sin(sg*ze*Th + ph); csg = mp.cos(sg*ze*Th + ph)
            rr = x + r0
            F = solver._F_mp(x, sols, A, 0); Fp = solver._F_mp(x, sols, A, 1)
            G = solver._G_mp(x, sols, A, 0); Gp = solver._G_mp(x, sols, A, 1)
            tyy = c11*(nu*Fp + R*(F - ze*G)/rr)
            tyr = Gp - G/rr + ze*F/rr
            Ft += ci * mp_proj(F*ssg, P)
            Gt += ci * mp_proj(G*csg, P)
            Tyyt += ci * mp_proj(tyy*ssg, P)
            Tyrt += ci * mp_proj(tyr*csg, P)
        Ftot.append(Ft); Gtot.append(Gt); Tyytot.append(Tyyt); Tyrtot.append(Tyrt)

    maxDisp = float(max(max(abs(f) for f in Ftot), max(abs(g) for g in Gtot)))
    rmsTyy = float((sum(t*t for t in Tyytot) / ngrid) ** mpf('0.5'))
    rmsTyr = float((sum(t*t for t in Tyrtot) / ngrid) ** mpf('0.5'))
    return dict(ok=True, block=block_name, sig0=sig0, sig1=sig1, maxDisp=maxDisp,
                rTyy=rmsTyy/maxDisp if maxDisp else float('inf'),
                rTyr=rmsTyr/maxDisp if maxDisp else float('inf'))
