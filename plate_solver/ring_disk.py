# -*- coding: utf-8 -*-
"""
plate_solver.ring_disk -- closed-ring (2 Theta = 2 pi) and solid-disk
(R_i = 0) entry points for Paper 3.

This module does not change default sector assembly. Frequencies still
come from the existing `_series_mp` / `_Lmat_mp` / `_det_mp` operators
in core_solvers.py. SOLVER_VERSION is not bumped.

Ring
----
2 pi-periodicity makes circumferential order an integer n. Each
frequency is a zero of the 4x4 free-edge (or mixed-edge) radial
determinant already coded as `_Lmat_mp`. OutOfPlaneSolver._Lmat_mp
rows 0,1 = M_r inner/outer; 2,3 = V_r inner/outer. InPlaneSolver._Lmat_mp
is N_r, N_r theta at both arcs. Searching det L(n, Omega) = 0 at fixed
integer n IS the ring.

The sector depth/gap bar does not travel here (LESSONS_LEARNED.md
Sec. 18.38). Adjudication is the real-part sign-flip ladder with
positive and negative controls in the same call.

Disk
----
Boundedness at r = 0 kills Y_n and K_n. Two remaining coefficients, two
outer-edge conditions: a 2x2. Job 2339453: evaluating the *4x4* at
R_i = 0 gives rank 2 at n = 0 and rank 3 (not 2) at n >= 1. Frequencies
from that matrix are not Leissa's disk.

    DISK_4X4_AT_RI0_IS_NOT_THE_DISK

Do not "fix" the rank-3 result. Production disk_search / disk_L_mp is a
2x2 on the regular pair. A tiny-hole 4x4 at b/a = 1e-3 may exist later
as a *control*, not as the production disk.

Phase D implements OOP F and C on each ring edge (F-F, F-C, C-F,
C-C). SS still raises UnsupportedRingBC (Phase F). F-F still calls
the validated `_det_mp`. Mixed-edge replaces the clamped Mr,Vr rows
by W, W'. Disk remains free-outer only. In-plane mixed-edge and IP
disk 2x2 are not Phase D.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Optional

import numpy as np
from mpmath import mp, mpf, mpc, matrix

from .core_solvers import OutOfPlaneSolver, InPlaneSolver
from .geometry import (
    IsotropicMaterial, make_geometry, _kbar_wbar, _omega_lit, _r0_2b_for_ratio,
)


# Job 2339453. The 4x4 at R_i=0 is NOT the disk. Rank 2 at n=0, rank 3
# (not 2) at n>=1. Do not "fix" the rank-3 result. Production disk is 2x2.
DISK_4X4_AT_RI0_IS_NOT_THE_DISK = True

# Sec. 18.38 confirmed pair (job 2444601). If an anchor ever moves, do
# not edit it -- that is a changed computed frequency.
SECTION_18_38_RADIUS_RATIO = 0.5
SECTION_18_38_NU = 0.30
SECTION_18_38_N2_OM = 0.1081885366
SECTION_18_38_N0_OM = 0.2359132157
SECTION_18_38_N2_LOGABSDET = -9.6838868232130650
SECTION_18_38_N0_LOGABSDET = -11.4312477453016985
SECTION_18_38_NEG_N0_OM = 0.2689680157
SECTION_18_38_NEG_N1_OM = 0.2689680157

# Paper 1 FF-P1 confirmed-real mode 7 (archive / Sec. 17 units writeup).
FFP1_MODE7_NATIVE = 1.369611
FFP1_MODE7_ANSYS_LIT = 54.0755

SIGNFLIP_DELTAS = (1.0e-2, 1.0e-3, 1.0e-4, 1.0e-5, 1.0e-6)
SIGNFLIP_KS = (-2, -1, 1, 2)

_RANK_TOL = 1e-12  # job 2339453 matrix_rank_mp


class UnsupportedRingBC(NotImplementedError):
    """SS (and IP mixed-edge) row-swaps are not Phase D."""


class DiskPathError(RuntimeError):
    """Raised if a disk routine would use the illegal 4x4 at R_i=0."""


# ---------------------------------------------------------------------------
# Conversion. Factor 39.478... is geometry-specific (Ri/Ro=0.5 => 4 pi^2),
# never an identity and never a universal constant (job 2406948).
# ---------------------------------------------------------------------------

def flexural_lambda2(Omega_native, geom, mat):
    """Map Seok Omega_native to lambda^2 = omega a^2 sqrt(rho / D).

    a = R_o, H = 2h (full thickness), D = c11_bar H^3 / 12. Delegates to
    `_kbar_wbar` / `_omega_lit` (part=1). Do not bake 4*pi**2 in.
    """
    k_bar, w_bar = _kbar_wbar(geom, mat)
    _f_hz, om_lit = _omega_lit(Omega_native, geom, mat, w_bar, k_bar, part=1)
    return float(om_lit)


def native_from_flexural_lambda2(lambda2, geom, mat):
    """Inverse of flexural_lambda2. Factor is recomputed from geom/mat."""
    factor = flexural_lambda2(1.0, geom, mat)
    if factor == 0.0:
        raise ZeroDivisionError("flexural conversion factor is zero")
    return float(lambda2) / factor


def inplane_lambda_irie(Omega_native, geom, mat):
    """Irie 1984 lambda = omega a sqrt(rho (1-nu^2) / E).

    This is not flexural lambda^2. `_omega_lit` part=2 still returns a
    flexural-style Omega_lit as its second value; that number is discarded.
    Physical frequency comes from the part=2 f_Hz conversion only.
    """
    k_bar, w_bar = _kbar_wbar(geom, mat)
    f_hz, _flexural_misuse = _omega_lit(
        Omega_native, geom, mat, w_bar, k_bar, part=2)
    omega = 2.0 * math.pi * float(f_hz)
    a = float(geom.R_o)
    E = float(mat.E)
    nu = float(mat.nu)
    rho = float(mat.rho)
    return omega * a * math.sqrt(rho * (1.0 - nu * nu) / E)


# ---------------------------------------------------------------------------
# Geometry helpers
# ---------------------------------------------------------------------------

def make_annulus_solver(radius_ratio, nu=0.30, motion="oop", M=80, n_quad=30,
                        two_theta=0.5, E=210e9, rho=7800.0):
    """Isotropic annular solver at given R_i/R_o. Angle is unused by L."""
    t = _r0_2b_for_ratio(float(radius_ratio))
    geom = make_geometry(t, two_theta)
    mat = IsotropicMaterial(E=E, nu=nu, rho=rho)
    motion = _norm_motion(motion)
    cls = OutOfPlaneSolver if motion == "oop" else InPlaneSolver
    solver = cls(geom, mat, M=M, n_quad=n_quad)
    return solver, geom, mat


def make_disk_solver(nu=0.30, motion="oop", M=80, n_quad=30,
                     two_theta=0.5, E=210e9, rho=7800.0):
    """Solid disk: R_i = 0 via _r0_2b_for_ratio(0) -> r0/(2b) = 0.5."""
    return make_annulus_solver(0.0, nu=nu, motion=motion, M=M, n_quad=n_quad,
                               two_theta=two_theta, E=E, rho=rho)


def _norm_motion(motion):
    m = str(motion).lower()
    if m in ("oop", "out-of-plane", "flexural", "part1", "1"):
        return "oop"
    if m in ("ip", "in-plane", "extensional", "part2", "2"):
        return "ip"
    raise ValueError("motion must be 'oop' or 'ip', got %r" % (motion,))


def _norm_bc(bc):
    s = str(bc).upper().strip()
    if s in ("F", "FREE"):
        return "F"
    if s in ("C", "CL", "CLAMPED"):
        return "C"
    if s in ("S", "SS", "SIMPLY", "SIMPLY-SUPPORTED", "SIMPLY_SUPPORTED"):
        return "S"
    raise ValueError("bc must be F, C, or S; got %r" % (bc,))


def _require_ring_bc(bc_inner, bc_outer, motion="oop"):
    """OOP F/C on each edge (Phase D). F-F always allowed. SS is Phase F."""
    inner = _norm_bc(bc_inner)
    outer = _norm_bc(bc_outer)
    if inner == "S" or outer == "S":
        raise UnsupportedRingBC(
            "SS row-swaps are Paper 3 Phase F, not Phase D; "
            "got inner=%r outer=%r" % (bc_inner, bc_outer))
    if (inner, outer) == ("F", "F"):
        return inner, outer
    motion = _norm_motion(motion)
    if motion != "oop":
        raise UnsupportedRingBC(
            "IP mixed-edge is not Phase D; got motion=%r inner=%r outer=%r"
            % (motion, bc_inner, bc_outer))
    if inner in ("F", "C") and outer in ("F", "C"):
        return inner, outer
    raise UnsupportedRingBC(
        "Phase D OOP allows F/C on each edge; got inner=%r outer=%r"
        % (bc_inner, bc_outer))


def _solver_motion(solver):
    if isinstance(solver, OutOfPlaneSolver):
        return "oop"
    if isinstance(solver, InPlaneSolver):
        return "ip"
    raise TypeError("expected OutOfPlaneSolver or InPlaneSolver, got %r"
                    % type(solver))


def _require_ri0(solver):
    ri = float(solver.geom.R_i)
    if abs(ri) > 1e-14:
        raise ValueError("disk path requires R_i=0, got R_i=%s" % ri)


# ---------------------------------------------------------------------------
# Ring 4x4 (validated path: solver._det_mp)
# ---------------------------------------------------------------------------

def ring_L_mp(solver, n, Om, bc_inner="F", bc_outer="F"):
    """4x4 characteristic matrix at integer n.

    F-F uses `_Lmat_mp` as-is (Mr, Vr at both arcs). A clamped edge
    replaces that arc's Mr, Vr rows by W, W' (G, G' at x=±π/2).
    """
    motion = _solver_motion(solver)
    inner, outer = _require_ring_bc(bc_inner, bc_outer, motion)
    xi = mpc(n)
    sols = solver._series_mp(xi, mpf(Om))
    L = solver._Lmat_mp(xi, sols)
    if (inner, outer) == ("F", "F"):
        return L
    hp = mp.pi / 2
    for edge, x, row_w, row_wp in (
            (inner, -hp, 0, 2),
            (outer, hp, 1, 3)):
        if edge != "C":
            continue
        for q in range(4):
            a = sols[q]
            L[row_w, q] = solver._ev_mp(a, x, 0)
            L[row_wp, q] = solver._ev_mp(a, x, 1)
    return L


def ring_det_mp(solver, n, Om, bc_inner="F", bc_outer="F"):
    """Equilibrated det of the 4x4. F-F calls the validated `_det_mp`."""
    motion = _solver_motion(solver)
    inner, outer = _require_ring_bc(bc_inner, bc_outer, motion)
    if (inner, outer) == ("F", "F"):
        return solver._det_mp(mpc(n), mpf(Om))
    L = ring_L_mp(solver, n, Om, inner, outer)
    for i in range(4):
        s = max(abs(L[i, j]) for j in range(4))
        if s > 0:
            for j in range(4):
                L[i, j] /= s
    return mp.det(L)


def ring_logabsdet(solver, n, Om, bc_inner="F", bc_outer="F"):
    """log10 |det L|. -300 if the determinant is literally zero."""
    d = ring_det_mp(solver, n, Om, bc_inner, bc_outer)
    ad = abs(complex(d))
    if ad <= 0.0:
        return -300.0
    return float(math.log10(ad))


# ---------------------------------------------------------------------------
# Disk 2x2. Never uses the 4x4 determinant at R_i=0.
# ---------------------------------------------------------------------------

def _regular_column_weights(solver, sols, n):
    """4x2 weights spanning the regular (bounded-at-r=0) subspace.

    The four `_series_mp` columns are mid-radius Taylor ICs, not the
    Bessel {J, Y, I, K} pair. A truncated series is finite at r=0, so
    "smallest |W|" does not isolate {J_n, I_n} (job 2461504: that SVD
    mix had no 2x2 zero). Regular solutions behave as r^{|n|}, so the
    inner Taylor jet supplies two homogeneous conditions whose 2-D
    kernel is the regular pair (sandbox 2026-09-08, Narita 2(b) natives):

      n = 0 : W' = W''' = 0   (even)
      n = 1 : W  = W''  = 0   (odd)
      n >= 2: W  = W'   = 0   (r^n, n>=2)

    Do not replace this by the inner M_r, V_r rows of the 4x4 (that is
    the illegal disk, job 2339453 rank 3).
    """
    x_inner = -mp.pi / 2
    n_int = int(round(abs(float(np.real(n)))))

    def col(d):
        return [complex(solver._ev_mp(sols[q], x_inner, d)) for q in range(4)]

    if n_int == 0:
        rows = (col(1), col(3))
    elif n_int == 1:
        rows = (col(0), col(2))
    else:
        rows = (col(0), col(1))
    A = np.array(rows, dtype=np.complex128)
    _u, _s, vh = np.linalg.svd(A, full_matrices=True)
    null = vh[A.shape[0]:]
    if null.shape[0] < 2:
        null = vh[-2:]
    C = null[:2]
    return [[mpc(complex(C[j, q])) for j in range(2)] for q in range(4)]


def disk_L_mp(solver, n, Om, bc_outer="F"):
    """2x2 outer-edge matrix on the regular pair. Not the 4x4 at R_i=0."""
    if not DISK_4X4_AT_RI0_IS_NOT_THE_DISK:
        raise DiskPathError("DISK_4X4_AT_RI0_IS_NOT_THE_DISK must stay True")
    _require_ri0(solver)
    if str(bc_outer).upper() != "F":
        raise UnsupportedRingBC(
            "Phase 0 disk implements free outer only; got bc_outer=%r"
            % (bc_outer,))
    xi = mpc(n)
    sols = solver._series_mp(xi, mpf(Om))
    # Outer rows of the 4x4: row 1 = M_r (OOP) / N_r (IP) at outer;
    # row 3 = V_r / N_r theta at outer. Inner rows are regularity, not BCs.
    L4 = solver._Lmat_mp(xi, sols)
    C = _regular_column_weights(solver, sols, n)
    L2 = matrix(2, 2)
    for j in range(2):
        for i, row4 in enumerate((1, 3)):
            acc = mpc(0)
            for q in range(4):
                acc += L4[row4, q] * C[q][j]
            L2[i, j] = acc
    if L2.rows != 2 or L2.cols != 2:
        raise DiskPathError(
            "FAIL_DISK_IS_4x4: disk_L_mp produced %sx%s, not 2x2"
            % (L2.rows, L2.cols))
    return L2


def disk_det_mp(solver, n, Om, bc_outer="F"):
    """Equilibrated det of the 2x2. Does not call `_det_mp` (that is 4x4)."""
    L = disk_L_mp(solver, n, Om, bc_outer=bc_outer)
    if L.rows != 2 or L.cols != 2:
        raise DiskPathError("FAIL_DISK_IS_4x4")
    for i in range(2):
        s = max(abs(L[i, j]) for j in range(2))
        if s > 0:
            for j in range(2):
                L[i, j] /= s
    return mp.det(L)


def disk_logabsdet(solver, n, Om, bc_outer="F"):
    d = disk_det_mp(solver, n, Om, bc_outer=bc_outer)
    ad = abs(complex(d))
    if ad <= 0.0:
        return -300.0
    return float(math.log10(ad))


def illegal_ri0_4x4_rank(solver, n, Om, tol=_RANK_TOL):
    """Numerical rank of the 4x4 at R_i=0. This is NOT the disk.

    Job 2339453 criterion (verbatim): one row-then-column max-abs
    equilibration, rank = #{s > tol * s[0]} with tol=1e-12. Expected:
    rank 2 at n=0, rank 3 (not 2) at n>=1. Do not 'fix' the rank-3 result.
    """
    _require_ri0(solver)
    L_mp = solver._Lmat_mp(mpc(n), solver._series_mp(mpc(n), mpf(Om)))
    nsz = L_mp.rows
    K = L_mp.copy()
    for i in range(nsz):
        s = max(abs(K[i, j]) for j in range(nsz))
        if s > 0:
            inv = mpf(1) / s
            for j in range(nsz):
                K[i, j] *= inv
    for j in range(nsz):
        s = max(abs(K[i, j]) for i in range(nsz))
        if s > 0:
            inv = mpf(1) / s
            for i in range(nsz):
                K[i, j] *= inv
    A = np.array([[complex(K[i, j]) for j in range(nsz)]
                  for i in range(nsz)], dtype=np.complex128)
    s = np.linalg.svd(A, compute_uv=False)
    s = np.sort(np.abs(s))[::-1]
    if s[0] <= 0.0:
        return 0, [-300.0] * nsz
    rank = int(np.sum(s > tol * s[0]))
    log_s = [-300.0 if v <= 0 else float(np.log10(v)) for v in s]
    return rank, log_s


# ---------------------------------------------------------------------------
# Golden section (verbatim contract of the Sec. 18.38 ring probes)
# ---------------------------------------------------------------------------

def golden_min(f, a, b, tol=2e-5, maxit=50):
    """Minimise scalar f on [a, b]. Returns (x_best, f_best, n_evals)."""
    phi = (1.0 + 5.0 ** 0.5) / 2.0
    invphi = 1.0 / phi
    invphi2 = invphi * invphi
    a, b = float(a), float(b)
    if b < a:
        a, b = b, a
    h = b - a
    if h <= tol:
        fa = f(a)
        return a, fa, 1
    n_plan = int(round(math.log(tol / h) / math.log(invphi))) + 1
    n_plan = min(max(n_plan, 2), maxit)
    c = a + invphi2 * h
    d = a + invphi * h
    fc = f(c)
    fd = f(d)
    evals = 2
    for _it in range(1, n_plan + 1):
        if fc < fd:
            b = d
            d = c
            fd = fc
            h = invphi * h
            c = a + invphi2 * h
            fc = f(c)
            evals += 1
        else:
            a = c
            c = d
            fc = fd
            h = invphi * h
            d = a + invphi * h
            fd = f(d)
            evals += 1
        if h < tol:
            break
    if fc < fd:
        x_best, f_best = c, fc
    else:
        x_best, f_best = d, fd
    for x_end in (a, b):
        f_end = f(x_end)
        evals += 1
        if f_end < f_best:
            x_best, f_best = x_end, f_end
    return x_best, f_best, evals


# ---------------------------------------------------------------------------
# Sign-flip ladder (Sec. 18.38 / job 2444601)
# ---------------------------------------------------------------------------

def _sgn(x):
    if x == 0.0:
        return 0
    return 1 if x > 0.0 else -1


@dataclass
class SignFlipRung:
    delta: float
    flip_in: bool
    flip_out: bool
    mag_inner: float
    re_minus: float
    re_plus: float


@dataclass
class SignFlipResult:
    n: int
    Omega: float
    rungs: list
    n_in: int
    n_out: int
    shrinking: bool
    clean: bool

    @property
    def n_deltas(self):
        return len(self.rungs)


@dataclass
class ControlSpec:
    kind: str  # 'POS' or 'NEG'
    n: int
    Omega: float
    solver: object = None  # default: the search solver
    disk: bool = False  # True only for a 2x2 disk control, never the 4x4 at R_i=0


@dataclass
class ControlResult:
    kind: str
    n: int
    Omega: float
    ladder: SignFlipResult
    ok: bool


@dataclass
class RingSearchResult:
    n: int
    Omega: float
    log_abs_det: float
    lambda2: Optional[float]
    ladder: SignFlipResult
    controls: list
    controls_ok: bool
    n_evals: int
    motion: str
    bc_inner: str
    bc_outer: str
    disk: bool = False
    window: tuple = field(default_factory=tuple)
    # Window-edge guard, set by _search_impl. "EXACT" if the returned
    # Omega IS a window end, "NEAR" if within one coarse step of one,
    # else "NONE". A boundary value is NOT a located root -- log|det|
    # merely fell monotonically to the edge, and the same edge comes
    # back for every n, which is the tell. Its lambda2 must not be
    # quoted as a frequency. See LESSONS_LEARNED.md Sec. 18.120.
    edge_flag: str = "NONE"
    edge_dist: float = float("inf")


def signflip_ladder(det_fn, n, Om0, deltas=SIGNFLIP_DELTAS):
    """Real-part sign-flip ladder, delta = 1e-2 .. 1e-6.

    A genuine simple root: Re(det) flips across Om0 at every delta,
    |Re(det)| shrinks toward Om0. A coincidental dip does not flip and
    converges to a nonzero constant (Sec. 18.38 negatives).

    det_fn(Om) -> complex. Clean = inner 5/5 (k = -1, +1).
    """
    Om0 = float(Om0)
    rungs = []
    for dl in deltas:
        vals = {}
        for k in SIGNFLIP_KS:
            z = complex(det_fn(Om0 + k * dl))
            vals[k] = (z.real, z.imag)
        s_in = (_sgn(vals[-1][0]), _sgn(vals[1][0]))
        s_out = (_sgn(vals[-2][0]), _sgn(vals[2][0]))
        flip_in = (s_in[0] != 0 and s_in[1] != 0 and s_in[0] != s_in[1])
        flip_out = (s_out[0] != 0 and s_out[1] != 0 and s_out[0] != s_out[1])
        mag = min(abs(vals[-1][0]), abs(vals[1][0]))
        rungs.append(SignFlipRung(
            delta=float(dl), flip_in=flip_in, flip_out=flip_out,
            mag_inner=float(mag), re_minus=float(vals[-1][0]),
            re_plus=float(vals[1][0])))
    n_in = sum(1 for r in rungs if r.flip_in)
    n_out = sum(1 for r in rungs if r.flip_out)
    mags = [r.mag_inner for r in rungs]
    shrinking = all(mags[i + 1] <= mags[i] * 1.5
                    for i in range(len(mags) - 1))
    clean = (n_in == len(rungs))
    return SignFlipResult(
        n=int(n), Omega=Om0, rungs=rungs, n_in=n_in, n_out=n_out,
        shrinking=shrinking, clean=clean)


def section_18_38_controls(solver_ba05):
    """Two POS and two NEG points from job 2444601, on the b/a=0.5 solver."""
    return [
        ControlSpec("POS", 0, SECTION_18_38_N0_OM, solver=solver_ba05),
        ControlSpec("POS", 2, SECTION_18_38_N2_OM, solver=solver_ba05),
        ControlSpec("NEG", 0, SECTION_18_38_NEG_N0_OM, solver=solver_ba05),
        ControlSpec("NEG", 1, SECTION_18_38_NEG_N1_OM, solver=solver_ba05),
    ]


def _det_fn_for(solver, n, bc_inner, bc_outer, disk):
    """Complex det(Om) on the 4x4 ring or the 2x2 disk; never mixes them."""
    if disk:
        def fn(Om, _s=solver, _n=n, _bo=bc_outer):
            return complex(disk_det_mp(_s, _n, Om, bc_outer=_bo))
        return fn

    def fn(Om, _s=solver, _n=n, _bi=bc_inner, _bo=bc_outer):
        return complex(ring_det_mp(_s, _n, Om, _bi, _bo))
    return fn


def _eval_controls(controls, search_solver, bc_inner, bc_outer):
    """Run the sign-flip ladder on each control. POS must be clean; NEG not.

    Each ControlSpec uses its own solver (default: the search solver) and
    its own disk flag. Section 18.38 POS/NEG are ring 4x4 even when the
    search itself is a disk 2x2.
    """
    out = []
    all_ok = True
    for spec in controls:
        kind = str(spec.kind).upper()
        slv = spec.solver if spec.solver is not None else search_solver
        det_fn = _det_fn_for(slv, spec.n, bc_inner, bc_outer, spec.disk)
        ladder = signflip_ladder(det_fn, spec.n, spec.Omega)
        if kind == "POS":
            ok = bool(ladder.clean)
        elif kind == "NEG":
            ok = (not ladder.clean)
        else:
            raise ValueError("control kind must be POS or NEG, got %r" % kind)
        out.append(ControlResult(kind=kind, n=int(spec.n),
                                 Omega=float(spec.Omega),
                                 ladder=ladder, ok=ok))
        all_ok = all_ok and ok
    return out, all_ok


def _coarse_deepest(f, lo, hi, step):
    """Grid f on [lo, hi]; return (x_deep, f_deep, n_evals, bracket)."""
    lo, hi = float(lo), float(hi)
    if hi < lo:
        lo, hi = hi, lo
    step = float(step)
    xs = []
    x = lo
    # inclusive end
    nmax = int(math.floor((hi - lo) / step + 0.5)) + 1
    for i in range(max(nmax, 2)):
        xi = lo + i * step
        if xi > hi + 0.25 * step:
            break
        xs.append(min(xi, hi))
    if xs[-1] < hi:
        xs.append(hi)
    vals = [(xi, f(xi)) for xi in xs]
    i_best = min(range(len(vals)), key=lambda i: vals[i][1])
    x_best, f_best = vals[i_best]
    # bracket for golden: neighbors of the deepest grid point
    i0 = max(0, i_best - 1)
    i1 = min(len(vals) - 1, i_best + 1)
    bracket = (vals[i0][0], vals[i1][0])
    if bracket[0] == bracket[1]:
        # single-point window: expand by one step if possible
        bracket = (max(lo, x_best - step), min(hi, x_best + step))
    return x_best, f_best, len(vals), bracket


def _search_impl(n, Om_window, solver, geom, mat, bc_inner, bc_outer,
                 motion, disk, coarse_step, polish, controls):
    n = int(n)
    lo, hi = Om_window
    det_fn = _det_fn_for(solver, n, bc_inner, bc_outer, disk)

    def logabs(Om):
        z = complex(det_fn(Om))
        ad = abs(z)
        return -300.0 if ad <= 0.0 else float(math.log10(ad))

    n_evals = 0
    x0, f0, n_grid, bracket = _coarse_deepest(logabs, lo, hi, coarse_step)
    n_evals += n_grid
    if polish:
        x0, f0, n_g = golden_min(logabs, bracket[0], bracket[1])
        n_evals += n_g
    # Window-edge guard (Sec. 18.120). Labels x0; never moves it, so no
    # computed frequency changes and SOLVER_VERSION does not bump.
    _d_edge = min(abs(float(x0) - float(lo)), abs(float(hi) - float(x0)))
    if _d_edge <= 1e-9:
        _edge = "EXACT"
    elif _d_edge <= float(coarse_step):
        _edge = "NEAR"
    else:
        _edge = "NONE"
    ladder = signflip_ladder(det_fn, n, x0)
    n_evals += 4 * len(ladder.rungs)
    ctrl_results = []
    controls_ok = True
    if controls:
        ctrl_results, controls_ok = _eval_controls(
            controls, solver, bc_inner, bc_outer)
        n_evals += sum(4 * len(c.ladder.rungs) for c in ctrl_results)
    lam2 = None
    if (not disk) and motion == "oop" and geom is not None and mat is not None:
        lam2 = flexural_lambda2(x0, geom, mat)
    elif disk and motion == "oop" and geom is not None and mat is not None:
        lam2 = flexural_lambda2(x0, geom, mat)
    return RingSearchResult(
        n=n, Omega=float(x0), log_abs_det=float(f0), lambda2=lam2,
        ladder=ladder, controls=ctrl_results, controls_ok=controls_ok,
        n_evals=n_evals, motion=motion, bc_inner=bc_inner, bc_outer=bc_outer,
        disk=bool(disk), window=(float(lo), float(hi)),
        edge_flag=_edge, edge_dist=float(_d_edge))


def ring_search(n, Om_window, bc_inner="F", bc_outer="F", motion="oop",
                solver=None, geom=None, mat=None, radius_ratio=None,
                nu=0.30, M=80, n_quad=30, coarse_step=0.0005, polish=True,
                controls=None):
    """Fix integer n, search Omega on det L = 0 of the 4x4.

    Om_window : (lo, hi) in native Seok Omega.
    controls  : optional sequence of ControlSpec; POS and NEG are
                evaluated with the same sign-flip ladder in this call.

    OOP F/C on each edge (Phase D). F-F still uses `_det_mp`. Does not
    use the sector residual/depth bar.
    """
    motion = _norm_motion(motion)
    bc_inner, bc_outer = _require_ring_bc(bc_inner, bc_outer, motion)
    if solver is None:
        if radius_ratio is None:
            if geom is None:
                raise ValueError(
                    "ring_search needs solver, geom, or radius_ratio")
            mat = mat or IsotropicMaterial(E=210e9, nu=nu, rho=7800.0)
            cls = OutOfPlaneSolver if motion == "oop" else InPlaneSolver
            solver = cls(geom, mat, M=M, n_quad=n_quad)
        else:
            solver, geom, mat = make_annulus_solver(
                radius_ratio, nu=nu, motion=motion, M=M, n_quad=n_quad)
    else:
        geom = geom if geom is not None else getattr(solver, "geom", None)
        mat = mat if mat is not None else getattr(solver, "mat", None)
    return _search_impl(n, Om_window, solver, geom, mat, bc_inner, bc_outer,
                        motion, False, coarse_step, polish, controls)


def disk_search(n, Om_window, bc_outer="F", motion="oop",
                solver=None, geom=None, mat=None, nu=0.30,
                M=80, n_quad=30, coarse_step=0.0005, polish=True,
                controls=None):
    """Fix integer n, search Omega on det L = 0 of the 2x2 regular disk.

    Asserts the 4x4 at R_i=0 is not used. OOP F (outer) only in Phase 0.
    """
    motion = _norm_motion(motion)
    if motion != "oop":
        raise NotImplementedError(
            "IP disk 2x2 is Paper 3 Phase E; Phase 0 is OOP only")
    if str(bc_outer).upper() != "F":
        raise UnsupportedRingBC(
            "Phase 0 disk implements free outer only; got bc_outer=%r"
            % (bc_outer,))
    if not DISK_4X4_AT_RI0_IS_NOT_THE_DISK:
        raise DiskPathError("FAIL_DISK_IS_4x4")
    if solver is None:
        solver, geom, mat = make_disk_solver(
            nu=nu, motion=motion, M=M, n_quad=n_quad)
    else:
        _require_ri0(solver)
        geom = geom if geom is not None else getattr(solver, "geom", None)
        mat = mat if mat is not None else getattr(solver, "mat", None)
    return _search_impl(n, Om_window, solver, geom, mat, "F", "F",
                        motion, True, coarse_step, polish, controls)
