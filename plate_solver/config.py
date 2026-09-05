# -*- coding: utf-8 -*-
"""
plate_solver.config — centralized runtime configuration.

BUG FIX (review-doc item 1, 2026-07-01): the thread-limiting env vars used to
be set with os.environ.setdefault(...). On a SLURM allocation the scheduler
injects OMP_NUM_THREADS (and friends) into the environment *before* this
script even starts, so setdefault silently left the scheduler's value in
place and BLAS still spawned node-core threads inside every worker -- the
exact nested-multithreading slowdown this code exists to prevent. Changed to
a direct assignment so this process's single-thread policy always wins,
regardless of what the scheduler pre-set. This MUST run before numpy/scipy
are imported anywhere in the package, which is why it lives at the very top
of this module and __init__ imports config first.

This module intentionally keeps the ORIGINAL flat module-level names
(MP_DPS, N_WORKERS, PREFETCH, PI, _P1_LD_GENUINE, ...) as aliases onto the
SolverConfig dataclass fields, so every extracted module below can still say
`from .config import MP_DPS` exactly as it could when this was one file --
zero behavioural change, no SOLVER_VERSION bump.

Scope note: this dataclass centralizes the ~10 tunables that are read at
import time and referenced by bare name from multiple modules. The other
~60 SIGMIN_*/FIG2_*/MODE_SHAPE_*/R40_*/SHI_* knobs documented in
LESSONS_LEARNED are each read locally, once, inside the single function that
uses them -- moving all 70 into one dataclass and rewiring every call site
was judged too high-risk to do blind in this pass (see LESSONS_LEARNED
Sec. 22); they remain plain `os.environ.get(...)` reads in place, still
fully overridable, just not yet centralized.
"""
from __future__ import annotations
import os
from dataclasses import dataclass

# -- Single-thread the BLAS backend BEFORE numpy/scipy are imported --------
# Parallelism is the ProcessPool (one mp worker/core); the tiny SVDs gain
# nothing from threaded BLAS, which instead spawns node-core threads inside
# EVERY worker and starves the mp work (the cluster-slowdown cause). Direct
# assignment (not setdefault) so a scheduler-injected value never wins.
for _thr_var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                 "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS",
                 "BLIS_NUM_THREADS"):
    os.environ[_thr_var] = "1"

import numpy as np
from mpmath import mp


def _allocated_cpus():
    """CPUs actually available to THIS process, scheduler/cgroup aware."""
    for _var in ("SLURM_CPUS_PER_TASK", "SLURM_CPUS_ON_NODE"):
        _v = os.environ.get(_var)
        if _v and _v.isdigit() and int(_v) > 0:
            return int(_v)
    try:
        n = len(os.sched_getaffinity(0))   # respects cpuset/cgroup limits
        if n > 0:
            return n
    except (AttributeError, OSError):
        pass
    return os.cpu_count() or 2


@dataclass(frozen=True)
class SolverConfig:
    """Globally-shared, cross-module tunables. See module docstring for the
    scope note on why this doesn't (yet) cover every env-var knob."""
    dps: int = 40                 # DPS   - mpmath decimal places (min 30)
    max_scan_pts: int = 80        # MAX_SCAN_PTS
    n_workers: int = 1            # N_WORKERS (computed from allocation below)
    prefetch: int = 1             # PREFETCH
    p1_ld_genuine: float = -2.0   # P1_LD_GENUINE
    p1_prominence: float = 0.3    # P1_PROMINENCE
    p2_ld_genuine: float = -3.0   # P2_LD_GENUINE
    solver_version: str = "2026-07-10.s10"

    @classmethod
    def from_env(cls) -> "SolverConfig":
        cpu = _allocated_cpus()
        n_workers = int(os.environ.get("N_WORKERS", str(max(1, min(64, cpu)))))
        return cls(
            dps=int(os.environ.get("DPS", "40")),
            max_scan_pts=int(os.environ.get("MAX_SCAN_PTS", "80")),
            n_workers=n_workers,
            prefetch=int(os.environ.get("PREFETCH", str(max(1, n_workers)))),
            p1_ld_genuine=float(os.environ.get("P1_LD_GENUINE", "-2.0")),
            p1_prominence=float(os.environ.get("P1_PROMINENCE", "0.3")),
            p2_ld_genuine=float(os.environ.get("P2_LD_GENUINE", "-3.0")),
            # BUMPED 2026-07-26 s8->s9: three default-affecting changes land
            # together -- (1) detectors.py's two-stage SIGMIN_POLISH_ITERS
            # default 5->12 (recovers G1=0.401154/G3=1.779792, previously
            # missed at iters=5; cluster-gated via jobs 2330819/2330983/
            # 2332566/2333514, GATE_FFP1+GATE_SHI both PASS, no known target
            # regresses -- this part IS cluster-validated); (2) geometry.py's
            # OrthotropicMaterial.ip_constants() now returns mu_theta instead
            # of the inherited nu_bar, the IP-side analogue of the fix
            # oop_constants() already carries, derived from scratch and
            # matched term-for-term to the governing equations + Wang et al.
            # 2016's own reciprocity relation (LESSONS_LEARNED Sec 24.2);
            # (3) core_solvers.py's InPlaneSolver._p2_worker_params() now
            # packs self.nu (mu_theta) instead of mat.nu_bar, fixing a
            # second, independent bug in the multiprocessing worker-scalar
            # path that (1) alone does not reach -- every worker-pool IP
            # search was silently using mu_r regardless of what
            # ip_constants() returns, because _ScalarMat never sees the
            # original material object. Both (2) and (3) are bit-identical
            # no-ops for isotropic materials (nu_bar==mu_theta there), so no
            # isotropic table is affected. NOTE: (2)+(3) together are NOT
            # yet cluster-validated against Wang Table 4 -- job 2332092's
            # "8/8 PASS" predates fix (3) and, because it went through the
            # worker pool, actually re-exercised the OLD nu_bar behavior
            # despite an in-process ip_constants() monkeypatch that never
            # reached the workers; that job must be re-run post-(3) before
            # IP orthotropy can be called cluster-validated. All three
            # changes can shift a computed frequency for orthotropic/
            # hard-bracket cases relative to s8, so per the invariant above
            # this is a real bump, not a refactor.
            # BUMPED 2026-07-28 s9->s10: core_solvers.py's
            # OutOfPlaneSolver._build_K_real CLAMPED-edge orientation rule
            # ("Geo-2", E.2/LESSONS_LEARNED Sec 38/46) changed from a
            # hardcoded orient=+1 for every CLAMPED edge to orient=e.sign
            # uniformly (matching the FREE-edge rule), with a matching
            # CLAMPED pairing-term sign split (sA=+1, sB=-1). Legacy folded
            # orient=+1 directly into the CLAMPED code path regardless of
            # which wall the edge sat at, so this DOES change every CLAMPED
            # computation's matrix values (not a no-op) -- verified this
            # stays within the published table's own tolerance for all 20
            # OOP Part-1 cantilever rows (job 2335764: 19/20 within
            # 0.0-0.5%, 1/20 within 1.3% after correcting a too-narrow
            # probe scan window, LESSONS_LEARNED Sec 46) but the actual
            # sigma_min values at a given Omega DO shift, so per the
            # invariant above this is a real bump. InPlaneSolver's own,
            # separate CLAMPED-orientation code is UNCHANGED (not examined
            # for this fix) and IP results are unaffected by this bump.
            solver_version="2026-07-10.s10",
        )


CONFIG = SolverConfig.from_env()
mp.dps = CONFIG.dps

# -- Backward-compatible flat aliases (unchanged names/values/semantics) ----
MP_DPS         = CONFIG.dps
MAX_SCAN_PTS   = CONFIG.max_scan_pts
N_WORKERS      = CONFIG.n_workers
PREFETCH       = CONFIG.prefetch
PI             = np.pi
_P1_LD_GENUINE = CONFIG.p1_ld_genuine
_P1_PROM       = CONFIG.p1_prominence
_P2_LD_GENUINE = CONFIG.p2_ld_genuine

# Solver-logic version stamp. Written into every checkpoint and checked on
# RESUME. BUMP ONLY when a code change alters a computed (isotropic)
# frequency by default -- never for comments, figures, or refactors (this
# module split is exactly such a refactor: NOT bumped).
SOLVER_VERSION = CONFIG.solver_version
