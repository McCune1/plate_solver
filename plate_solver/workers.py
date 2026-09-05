# -*- coding: utf-8 -*-
"""
plate_solver.workers -- multiprocessing worker-side functions.

These are the functions ProcessPoolExecutor actually pickles and ships to
subprocesses. They MUST stay top-level (not nested) and MUST NOT rely on
anything the worker process can't cheaply reconstruct -- see the
worker-reconstruction-path warning in LESSONS_LEARNED and the
plate-solver-sandbox-probe skill's regression gate. This is exactly why they
were split into their own module rather than left scattered: it makes the
worker/in-process boundary explicit instead of implicit.

Extracted verbatim (line-range provenance in LESSONS_LEARNED Sec. 22); no
numeric behaviour changed, SOLVER_VERSION not bumped.
"""
from __future__ import annotations

from .config import PREFETCH
from .geometry import _make_geom_mat
from .boundary import _worker_bc


# Lazy accessor to avoid a module-load-time circular import: core_solvers.py
# imports several *_worker functions from THIS module, so this module must
# not import core_solvers at module scope. Workers that reconstruct a solver
# call this once, at call time (inside the already-spawned worker process),
# by which point the whole plate_solver package is fully initialized.
def _get_core_solvers():
    from . import core_solvers as _cs
    return _cs.OutOfPlaneSolver, _cs.InPlaneSolver


# Second lazy accessor, same rationale: several workers call full_search/
# select_fill/track/_min_root_spacing, which live in detectors.py -- and
# detectors.py imports worker functions FROM this module at module scope,
# so this module must not import detectors at module scope either.
def _get_detectors():
    from . import detectors as _dt
    return _dt.full_search, _dt.select_fill, _dt.track, _dt._min_root_spacing

_WORKER_EXC_SURFACED = False

class _StaleCheckpoint(Exception):
    """Raised internally when a checkpoint's solver_version does not match the
    current SOLVER_VERSION, so RESUME ignores it instead of replaying stale
    results from an older code version."""
    pass

# Worker pool size.  Each worker holds an mp K matrix (~100-200 MB at dps=40,
# 16-dof Part 1); keep N_WORKERS × 200 MB well inside job RAM.  os.cpu_count()
# over-counts a shared node, so size from the SLURM allocation / CPU-affinity
# mask (one single-threaded worker per allocated core).  See §Cluster.
def _throttled_map(pool, fn, args_list, prefetch=None):
    """Submit jobs to a pool one at a time, keeping at most prefetch ahead.

    Unlike pool.map(), which serialises the entire argument list into the IPC
    pipe before any work starts, this function submits up to (N_WORKERS +
    prefetch) jobs and then yields results as each completes.  This bounds
    peak memory to a small multiple of the per-job cost regardless of how
    many total jobs there are.

    Results are returned in submission order.  If a job raises an exception
    the exception is re-raised here so the caller can handle it.

    Parameters
    ----------
    pool     : ProcessPoolExecutor already entered as a context manager
    fn       : module-level callable accepted by the pool
    args_list: iterable of single argument objects passed to fn
    prefetch : maximum jobs to have queued ahead of completion (default: N_WORKERS)

    Yields
    ------
    The return value of fn(arg) for each arg, in submission order.
    """
    if prefetch is None:
        prefetch = PREFETCH

    # Use a deque as a sliding window of in-flight futures.
    from collections import deque
    window = deque()
    it = iter(args_list)

    # Pre-fill the window up to prefetch slots.
    for arg in it:
        window.append(pool.submit(fn, arg))
        if len(window) >= prefetch:
            break

    # As each future finishes, yield its result and submit the next job.
    for arg in it:
        fut = window.popleft()
        yield fut.result()
        window.append(pool.submit(fn, arg))

    # Drain remaining futures.
    for fut in window:
        yield fut.result()


# ══════════════════════════════════════════════════════════════════════════════
# ══════════════════════════════════════════════════════════════════════════════
#  SECTION 2b — CUT-OFF FREQUENCY COMPUTATION  (ξ=0 / ζ=0 limits)
#  At ξ=0/ζ=0 the radial ODEs decouple into Bessel equations with closed-form
#  2×2 determinants; their roots are exact, geometry/ν-only cut-off frequencies
#  that anchor the dispersion branches.  See LESSONS_LEARNED.md §Cut-offs.
# ══════════════════════════════════════════════════════════════════════════════

def _surface_worker_exc(exc, where=""):
    """Print this worker process's first exception (repr + traceback) once."""
    global _WORKER_EXC_SURFACED
    if _WORKER_EXC_SURFACED:
        return
    _WORKER_EXC_SURFACED = True
    import sys, os, traceback
    tag = f' [{where}]' if where else ''
    sys.stderr.write(f'  [worker-exc pid={os.getpid()}]{tag} first failure: '
                     f'{type(exc).__name__}: {exc}\n')
    traceback.print_exc(file=sys.stderr)
    sys.stderr.flush()


# ── Part 1 scan worker ────────────────────────────────────────────────────────

def _p1_logdet_worker(args):
    """Evaluate Part 1 log|det K| at one (Om, branch-set) pair.

    args = (Om_f, brs_list, g_scalars, m_scalars, M, n_quad, dps)
    brs_list : list of complex numbers (canonical float64 branch roots) used as
               Newton seed — track() is called from these to Om_f.
               If None, a fresh full_search is performed (legacy, slow path).
    Returns  : (Om_f, sign:int, log10|det|:float)
    """
    Om_f, brs_list, g_scalars, m_scalars, M, n_quad, dps = args
    geom, mat = _make_geom_mat(g_scalars, m_scalars, dps)
    try:
        OutOfPlaneSolver, InPlaneSolver = _get_core_solvers()
        full_search, select_fill, track, _min_root_spacing = _get_detectors()
        solver = OutOfPlaneSolver(geom, mat, M=M, n_quad=n_quad, boundary=_worker_bc('out_of_plane'))
        if brs_list is None:
            # Legacy slow path: fresh full_search (kept for compatibility)
            xi_max = 14.0
            raw = full_search(solver.fast, Om_f, xmax=xi_max)
            sel, _ = select_fill(raw, 16)
            brs_list = list(sel)
        else:
            # Fast path: track from provided seed (a single Newton step per root)
            tracked = track(solver.fast, Om_f, [complex(z) for z in brs_list])
            sel, cnt = select_fill(tracked, 16)
            if cnt < 16:
                # Fall back to provided seed directly if tracking lost roots
                sel, _ = select_fill([complex(z) for z in brs_list], 16)
            brs_list = list(sel)
        sgn, ld = solver.logdet(Om_f, brs_list, fast_scan=False)
        return (Om_f, int(sgn), float(ld))
    except Exception as _e:
        _surface_worker_exc(_e, 'logdet_worker')
        return (Om_f, 0, float('nan'))


# ── Part 1 polish workers ─────────────────────────────────────────────────────

def _p1_bisect_worker(args):
    """Bisect a Part 1 sign-change bracket in a worker process.

    args = (a, b, g_scalars, m_scalars, M, n_quad, xi_max, n_dofs, iters, nq_fast, dps, seed_brs)
    Returns : (Om_star:float, ld_star:float)
    """
    a, b, g_scalars, m_scalars, M, n_quad, xi_max, n_dofs, iters, nq_fast, dps, seed_brs = args
    geom, mat = _make_geom_mat(g_scalars, m_scalars, dps)
    OutOfPlaneSolver, InPlaneSolver = _get_core_solvers()
    solver = OutOfPlaneSolver(geom, mat, M=M, n_quad=n_quad, boundary=_worker_bc('out_of_plane'))
    return solver._bisect(a, b, n_dofs, xi_max, iters,
                          n_quad_iter=nq_fast, seed_brs=seed_brs)


def _p1_golden_worker(args):
    """Golden-section minimise Part 1 log|det K| in a worker process.

    args = (a, b, g_scalars, m_scalars, M, n_quad, xi_max, n_dofs, iters, nq_fast, dps, seed_brs)
    Returns : (Om_star:float, ld_star:float)
    """
    a, b, g_scalars, m_scalars, M, n_quad, xi_max, n_dofs, iters, nq_fast, dps, seed_brs = args
    geom, mat = _make_geom_mat(g_scalars, m_scalars, dps)
    OutOfPlaneSolver, InPlaneSolver = _get_core_solvers()
    solver = OutOfPlaneSolver(geom, mat, M=M, n_quad=n_quad, boundary=_worker_bc('out_of_plane'))
    return solver._golden(a, b, n_dofs, xi_max, iters,
                          n_quad_iter=nq_fast, seed_brs=seed_brs)


def _p1_ld_at_worker(args):
    """Evaluate Part 1 log|det| at a single point (genuine-mode 3-point filter).

    args = (Om, g_scalars, m_scalars, M, n_quad, xi_max, n_dofs, dps, seed_brs)
    Returns : float(ld)
    """
    Om, g_scalars, m_scalars, M, n_quad, xi_max, n_dofs, dps, seed_brs = args
    geom, mat = _make_geom_mat(g_scalars, m_scalars, dps)
    OutOfPlaneSolver, InPlaneSolver = _get_core_solvers()
    solver = OutOfPlaneSolver(geom, mat, M=M, n_quad=n_quad, boundary=_worker_bc('out_of_plane'))
    return solver._ld_at(Om, n_dofs, xi_max, seed_brs=seed_brs)


def _p1_sigmin_worker(args):
    """Evaluate Part 1 log10 σ_min(K) at one (Om, branch-set) pair (scan phase).

    args = (Om_f, brs_list, g_scalars, m_scalars, M, n_quad, xi_max, n_dofs, dps)
    Returns : (Om_f, log10σ_min:float)  — large negative ⇒ near a mode.
    """
    Om_f, brs_list, g_scalars, m_scalars, M, n_quad, xi_max, n_dofs, dps = args
    geom, mat = _make_geom_mat(g_scalars, m_scalars, dps)
    try:
        OutOfPlaneSolver, InPlaneSolver = _get_core_solvers()
        full_search, select_fill, track, _min_root_spacing = _get_detectors()
        solver = OutOfPlaneSolver(geom, mat, M=M, n_quad=n_quad, boundary=_worker_bc('out_of_plane'))
        if brs_list is None:
            raw = full_search(solver.fast, Om_f, xmax=xi_max)
            sel, _ = select_fill(raw, n_dofs)
        else:
            tracked = track(solver.fast, Om_f, [complex(z) for z in brs_list])
            sel, cnt = select_fill(tracked, n_dofs)
            if cnt < n_dofs:
                sel, _ = select_fill([complex(z) for z in brs_list], n_dofs)
        s = solver.sigma_min(Om_f, list(sel), fast_scan=False)
        return (Om_f, float(s))
    except Exception as _e:
        _surface_worker_exc(_e, 'p1_sigmin_worker')
        return (Om_f, float('nan'))


def _p1_sigmin_golden_worker(args):
    """Refine a Part 1 σ_min minimum (golden-section) in a worker process.

    args = (a, b, g_scalars, m_scalars, M, n_quad, xi_max, n_dofs, iters, nq_fast, dps, seed_brs)
    Returns : (Om_star:float, log10σ_min:float)
    """
    a, b, g_scalars, m_scalars, M, n_quad, xi_max, n_dofs, iters, nq_fast, dps, seed_brs = args
    geom, mat = _make_geom_mat(g_scalars, m_scalars, dps)
    OutOfPlaneSolver, InPlaneSolver = _get_core_solvers()
    solver = OutOfPlaneSolver(geom, mat, M=M, n_quad=n_quad, boundary=_worker_bc('out_of_plane'))
    return solver._golden_sigmin(a, b, n_dofs, xi_max, iters,
                                 n_quad_iter=nq_fast, seed_brs=seed_brs)


# ── Part 2 scan worker ────────────────────────────────────────────────────────

def _p2_logdet_worker(args):
    """Evaluate Part 2 log|det K̂| (constrained + unconstrained) at one scan point.

    args = (Om_f, brs_list, g_scalars, m_scalars, M, n_quad, ze_max, n_dofs, dps)
    brs_list : pre-computed tracked branch set (for constrained eval)
    Returns  : (Om_f, sgn_c:int, ld_c:float, ld_u:float)
    """
    Om_f, brs_list, g_scalars, m_scalars, M, n_quad, ze_max, n_dofs, dps = args
    geom, mat = _make_geom_mat(g_scalars, m_scalars, dps)
    try:
        OutOfPlaneSolver, InPlaneSolver = _get_core_solvers()
        solver = InPlaneSolver(geom, mat, M=M, n_quad=n_quad, boundary=_worker_bc('in_plane'))
        # Constrained matrix with full mp root refinement (fast_scan=False).
        brs = [complex(z) for z in brs_list]
        sgn_c, ld_c = solver.logdet(Om_f, brs, lagrange=True, fast_scan=False)
        # Unconstrained: track from the same seed instead of full_search
        full_search, select_fill, track, _min_root_spacing = _get_detectors()
        tracked = track(solver.fast, Om_f, brs)
        sel_f, _ = select_fill(tracked if tracked else brs, n_dofs)
        _, ld_u = solver.logdet(Om_f, list(sel_f), lagrange=False, fast_scan=False)
        return (Om_f, int(sgn_c), float(ld_c), float(ld_u))
    except Exception as _e:
        _surface_worker_exc(_e, 'logdet_conuncon_worker')
        return (Om_f, 0, float('nan'), float('nan'))


# ── Part 2 polish workers ─────────────────────────────────────────────────────

def _p2_bisect_worker(args):
    """Bisect a Part 2 constrained sign-change bracket in a worker process.

    args = (a, b, g_scalars, m_scalars, M, n_quad, ze_max, n_dofs, iters, nq_fast, dps, seed_brs)
    Returns : (Om_star:float, ld_star:float)
    """
    a, b, g_scalars, m_scalars, M, n_quad, ze_max, n_dofs, iters, nq_fast, dps, seed_brs = args
    geom, mat = _make_geom_mat(g_scalars, m_scalars, dps)
    OutOfPlaneSolver, InPlaneSolver = _get_core_solvers()
    solver = InPlaneSolver(geom, mat, M=M, n_quad=n_quad, boundary=_worker_bc('in_plane'))
    return solver._bisect(a, b, n_dofs, ze_max, iters,
                          n_quad_iter=nq_fast, seed_brs=seed_brs)


def _p2_golden_worker(args):
    """Golden-section minimise Part 2 UNCONSTRAINED log|det K| in a worker process.

    args = (a, b, g_scalars, m_scalars, M, n_quad, ze_max, n_dofs, iters, nq_fast, dps, seed_brs)
    Returns : (Om_star:float, ld_c:float)
    """
    a, b, g_scalars, m_scalars, M, n_quad, ze_max, n_dofs, iters, nq_fast, dps, seed_brs = args
    geom, mat = _make_geom_mat(g_scalars, m_scalars, dps)
    OutOfPlaneSolver, InPlaneSolver = _get_core_solvers()
    solver = InPlaneSolver(geom, mat, M=M, n_quad=n_quad, boundary=_worker_bc('in_plane'))
    return solver._golden(a, b, n_dofs, ze_max, iters,
                          n_quad_iter=nq_fast, seed_brs=seed_brs)


def _p2_ld_at_worker(args):
    """Evaluate Part 2 log|det| at a single point (used in genuine-mode filter).

    args = (Om, g_scalars, m_scalars, M, n_quad, ze_max, n_dofs, lagrange, dps, seed_brs)
    Returns : float(ld)
    """
    Om, g_scalars, m_scalars, M, n_quad, ze_max, n_dofs, lagrange, dps, seed_brs = args
    geom, mat = _make_geom_mat(g_scalars, m_scalars, dps)
    OutOfPlaneSolver, InPlaneSolver = _get_core_solvers()
    solver = InPlaneSolver(geom, mat, M=M, n_quad=n_quad, boundary=_worker_bc('in_plane'))
    return solver._ld_at(Om, n_dofs, ze_max, lagrange, seed_brs=seed_brs)


def _p2_sigmin_worker(args):
    """Evaluate Part 2 log10 σ_min of the UNCONSTRAINED in-plane K at one point.

    PAPER PHYSICS (corrected this pass): the paper's natural frequencies are
    the DOUBLE ROOTS of the UNCONSTRAINED variational equation Eq. (48) — that
    is the actual mode condition.  The Lagrange-CONSTRAINED matrix K# (Eq.
    52-53) is an auxiliary device the paper uses only to split those double
    roots so the eigenvector amplitudes can be determined; about half of K#'s
    roots are spurious and discarded.  Verified directly on (1.25, 0.25π): at
    every paper mode the UNCONSTRAINED σ_min dips deeply (m1 -3.9, m2 -4.0,
    m3 -4.4) while the constrained σ_min is SHALLOW there (m1 -3.3, m2 -2.8,
    m3 -3.7) — and at a spurious point (Ω=0.6445) it's the reverse
    (constrained -5.5 deep, unconstrained -1.8 shallow).  So the previous
    constrained-as-primary detector both MISSED real modes (shallow
    constrained σ_min fell below the accept floor) and ACCEPTED spurious ones.
    The unconstrained matrix is now the primary signal.

    args = (Om_f, brs_list, g_scalars, m_scalars, M, n_quad, ze_max, n_dofs, dps)
    Returns : (Om_f, log10σ_min_unconstrained:float)
    """
    Om_f, brs_list, g_scalars, m_scalars, M, n_quad, ze_max, n_dofs, dps = args
    geom, mat = _make_geom_mat(g_scalars, m_scalars, dps)
    try:
        OutOfPlaneSolver, InPlaneSolver = _get_core_solvers()
        full_search, select_fill, track, _min_root_spacing = _get_detectors()
        solver = InPlaneSolver(geom, mat, M=M, n_quad=n_quad, boundary=_worker_bc('in_plane'))
        if brs_list is None:
            raw = full_search(solver.fast, Om_f, xmax=ze_max)
            sel, _ = select_fill(raw, n_dofs)
        else:
            tracked = track(solver.fast, Om_f, [complex(z) for z in brs_list])
            sel, cnt = select_fill(tracked, n_dofs)
            if cnt < n_dofs:
                sel, _ = select_fill([complex(z) for z in brs_list], n_dofs)
        s = solver.sigma_min(Om_f, list(sel), lagrange=False, fast_scan=False)
        return (Om_f, float(s))
    except Exception as _e:
        _surface_worker_exc(_e, 'p2_sigmin_worker')
        return (Om_f, float('nan'))


def _p2_sigmin_golden_worker(args):
    """Refine a Part 2 σ_min minimum (golden-section) in a worker process.

    args = (a, b, g_scalars, m_scalars, M, n_quad, ze_max, n_dofs, iters, nq_fast, dps, seed_brs)
    Returns : (Om_star:float, log10σ_min_constrained:float, log10|det_unconstrained|:float)
    """
    a, b, g_scalars, m_scalars, M, n_quad, ze_max, n_dofs, iters, nq_fast, dps, seed_brs = args
    geom, mat = _make_geom_mat(g_scalars, m_scalars, dps)
    OutOfPlaneSolver, InPlaneSolver = _get_core_solvers()
    solver = InPlaneSolver(geom, mat, M=M, n_quad=n_quad, boundary=_worker_bc('in_plane'))
    return solver._golden_sigmin(a, b, n_dofs, ze_max, iters,
                                 n_quad_iter=nq_fast, seed_brs=seed_brs)


def _preflight_probe_worker(args):
    """float64 conditioning probe at one Ω (BC-independent: the branch roots
    come from the exact-edge arcs, not the θ-edge BC).  Returns root-count and
    conditioning signals so under-resolution / degeneracy is caught up front."""
    part, key, Om, g_sc, m_sc, M, xmax, n_dofs, dps = args
    geom, mat = _make_geom_mat(g_sc, m_sc, dps)
    try:
        # Both names must be fetched unconditionally: InPlaneSolver was only
        # bound inside the part==1 branch, so the part==2 (else) branch hit
        # UnboundLocalError (Python treats a name assigned anywhere in a
        # function as local to the whole function). Fixed 2026-07-01.
        OutOfPlaneSolver, InPlaneSolver = _get_core_solvers()
        if part == 1:
            solver = OutOfPlaneSolver(geom, mat, M=M, n_quad=20,
                                      boundary=_worker_bc("out_of_plane"))
        else:
            solver = InPlaneSolver(geom, mat, M=M, n_quad=20,
                                   boundary=_worker_bc("in_plane"))
        full_search, select_fill, track, _min_root_spacing = _get_detectors()
        raw = full_search(solver.fast, float(Om), xmax=float(xmax))
        nr = len(raw)
        sel, filled = select_fill(raw, n_dofs)
        sp = _min_root_spacing(raw) if nr >= 2 else float("inf")
        ab = min(abs(complex(z)) for z in raw) if raw else float("nan")
        return (part, key, float(Om), int(nr), int(filled), float(sp), float(ab))
    except Exception as _e:
        _surface_worker_exc(_e, "preflight_probe")
        return (part, key, float(Om), -1, -1, float("nan"), float("nan"))


