# -*- coding: utf-8 -*-
"""
plate_solver.boundary -- boundary-condition registry (EdgeKind/EdgeSpec,
BoundaryCondition and the four concrete Clamped/Free x OOP/IP classes) and
the worker-side BC token lookup.

Extracted verbatim from Research50.py; no numeric behaviour changed.
"""
from __future__ import annotations
import os

class EdgeKind:
    CLAMPED = "clamped"   # essential: W=dW/dtheta=0 (OOP) / u=0 (IP), weakly imposed
    FREE    = "free"      # natural: moment/shear (OOP) or traction (IP) vanish

class EdgeSpec:
    """One straight theta-edge.  sign=+1 -> theta=+Theta, sign=-1 -> theta=-Theta."""
    __slots__ = ("sign", "kind")
    def __init__(self, sign, kind):
        self.sign = int(sign); self.kind = kind
    def __repr__(self):
        return f"EdgeSpec(theta={'+' if self.sign>0 else '-'}Theta, {self.kind})"

# OOP Kirchhoff theta-corner coefficient by edge kind.  The cantilever variational
# form (Seok & Tiersten Eq.43) gives the clamped wall +1 and the free edge -2;
# this dict reproduces that and lets a free-free plate use -2 on BOTH free edges.
_OOP_CORNER_COEFF = {EdgeKind.CLAMPED: 1, EdgeKind.FREE: -2}

class BoundaryCondition:
    """Describes the two theta-edges.  `edges` order also fixes the assembly
    order (kept as [+Theta, -Theta] so the cantilever reproduces the original
    summation bit-stably)."""
    name = "abstract"
    validated = True          # False for configs without a reference table
    def __init__(self, edges):
        self.edges = list(edges)
    def clamped_edges(self):
        return [e for e in self.edges if e.kind == EdgeKind.CLAMPED]
    @property
    def has_clamped(self):
        return any(e.kind == EdgeKind.CLAMPED for e in self.edges)
    def __repr__(self):
        return f"{type(self).__name__}({self.name})"

class ClampedFreeOOP(BoundaryCondition):
    token = "clamped_free"
    name = "clamped(-Theta)-free(+Theta)"
    def __init__(self):
        super().__init__([EdgeSpec(+1, EdgeKind.FREE), EdgeSpec(-1, EdgeKind.CLAMPED)])

class ClampedFreeIP(BoundaryCondition):
    token = "clamped_free"
    name = "clamped(-Theta)-free(+Theta)"
    def __init__(self):
        super().__init__([EdgeSpec(+1, EdgeKind.FREE), EdgeSpec(-1, EdgeKind.CLAMPED)])

class FreeFreeOOP(BoundaryCondition):
    token = "free_free"
    # NOTE: no Seok & Tiersten reference table exists for the free-free annulus
    # (the source papers never treat this configuration) -- so there is no
    # literature table to reproduce. The symmetric corner coefficient (-2 on
    # both free edges) is FE-validated instead: 8/8 isotropic targets <=0.64%
    # against a mesh-converged ANSYS SHELL281 model and 8/8 orthotropic targets
    # <=1.08% against an independent ANSYS benchmark (Shi, Liang & Wang material
    # properties). Closed 2026-07-12/13 -- see LESSONS_LEARNED.md Sec 20 and the
    # paper's Sec 6.2-6.3.
    name = "free(-Theta)-free(+Theta)"
    validated = True
    def __init__(self):
        super().__init__([EdgeSpec(+1, EdgeKind.FREE), EdgeSpec(-1, EdgeKind.FREE)])

class FreeFreeIP(BoundaryCondition):
    token = "free_free"
    # NOTE: no Seok & Tiersten reference table exists for the free-free annulus
    # -- FE-validated instead: all 27 independent ANSYS PLANE183 modes below
    # 1120 Hz matched to <=0.052% (26 distinct solver roots). Closed 2026-07-12/13
    # -- see LESSONS_LEARNED.md Sec 20 and the paper's Sec 6.3.
    name = "free(-Theta)-free(+Theta)"
    validated = True
    def __init__(self):
        super().__init__([EdgeSpec(+1, EdgeKind.FREE), EdgeSpec(-1, EdgeKind.FREE)])

_BC_REGISTRY = {
    ("clamped_free", "out_of_plane"): ClampedFreeOOP,
    ("clamped_free", "in_plane"):     ClampedFreeIP,
    ("free_free",    "out_of_plane"): FreeFreeOOP,
    ("free_free",    "in_plane"):     FreeFreeIP,
}

def make_bc(kind, motion):
    """Factory: kind in {'clamped_free','free_free'}, motion in
    {'out_of_plane','in_plane'}.  Lets a BC token be threaded through worker
    scalars (Step 6) so parallel free-free reconstruction picks the right BC."""
    try:
        return _BC_REGISTRY[(kind, motion)]()
    except KeyError:
        raise ValueError(f"no boundary condition for kind={kind!r} motion={motion!r}")

_BC_REGISTRY_TOKENS = {t for (t, _m) in _BC_REGISTRY}

def _worker_bc(motion):
    """Reconstruct this solve's boundary condition inside a worker process
    from the inherited R40_BC_KIND environment variable.  Defaults to
    clamped_free, so any pool launched without a BC set (every existing
    cantilever run, and all figure-layer pools) behaves exactly as before.
    The env var is inherited by workers at fork/spawn on all platforms, so
    this reaches EVERY pool -- including the shared scan/refine helpers --
    without threading a token through any worker argument tuple."""
    return make_bc(os.environ.get('R40_BC_KIND', 'clamped_free'), motion)


