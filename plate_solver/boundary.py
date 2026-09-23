# -*- coding: utf-8 -*-
"""
plate_solver.boundary -- boundary-condition registry (EdgeKind/EdgeSpec,
BoundaryCondition and the Clamped/Free x OOP/IP classes) and the
worker-side BC token lookup.

Extracted verbatim from Research50.py; no numeric behaviour changed.
P2-easy (2026-09-11): free_clamped (reverse cantilever) and
clamped_clamped (both theta-edges clamped, arcs still free) added.
Corner coefficients are the already-validated CLAMPED=+1 / FREE=-2
pair; SOLVER_VERSION is not bumped. IP clamped-at-+Theta orientation
is still the pre-E.2 rule -- see probe_p2easy_reverse_identity.
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

class FreeClampedOOP(BoundaryCondition):
    """Reverse cantilever: FREE at -Theta, CLAMPED at +Theta.

    OOP Geo-2 (E.2, SOLVER_VERSION s10) made this the mirror of
    ClampedFreeOOP. validated=False until the identity probe closes it
    through the worker path; it is not a new frequency table.
    """
    token = "free_clamped"
    name = "free(-Theta)-clamped(+Theta)"
    validated = False
    def __init__(self):
        super().__init__([EdgeSpec(+1, EdgeKind.CLAMPED), EdgeSpec(-1, EdgeKind.FREE)])

class FreeClampedIP(BoundaryCondition):
    """Reverse cantilever, in-plane.

    IP E.2 (2026-09-14, SOLVER_VERSION s10 -- LESSONS_LEARNED Sec 18.135)
    ported Geo-2 to InPlaneSolver._build_K_real: CLAMPED edges now use
    orient=e.sign uniformly (was hardcoded +1, and unused for CLAMPED
    edges at all), with NO internal sign split of the CLAMPED pairing
    term needed (unlike OOP's sA=+1/sB=-1 -- confirmed empirically, not
    assumed to carry over). Sandbox-verified (n_dofs=8, dps=26,
    r0/2b=1.25, 2T/pi in {0.25,0.5,1.0}): the published cantilever
    (ClampedFreeIP, CLAMPED@-Theta) is bit-identical to pre-fix, and
    std/rev sigma_min now match EXACTLY (delta=0.0000) at every point
    tried. validated=False until the cluster reverse-identity probe
    (re-run of probe_p2easy_reverse_identity_2026-09-11.py against this
    fix) confirms G2 PASS on real hardware -- see that probe's own
    pre-registered G2 reading.
    """
    token = "free_clamped"
    name = "free(-Theta)-clamped(+Theta)"
    validated = False
    def __init__(self):
        super().__init__([EdgeSpec(+1, EdgeKind.CLAMPED), EdgeSpec(-1, EdgeKind.FREE)])

class ClampedClampedOOP(BoundaryCondition):
    """Both radial walls clamped; inner/outer arcs stay free.

    Corner coefficient +1 on both walls (same pairing as the cantilever
    mixed corner). Not a published table in this project.
    """
    token = "clamped_clamped"
    name = "clamped(-Theta)-clamped(+Theta)"
    validated = False
    def __init__(self):
        super().__init__([EdgeSpec(+1, EdgeKind.CLAMPED), EdgeSpec(-1, EdgeKind.CLAMPED)])

class ClampedClampedIP(BoundaryCondition):
    """Both radial walls clamped, in-plane.

    IP E.2 (2026-09-14) changes THIS BC's own spectrum, not just the
    reverse cantilever's: one of the two clamped walls sits at +Theta, and
    that wall's contribution sign-flips under the fix (the -Theta wall is
    bit-identical to pre-fix; see FreeClampedIP). Sandbox spot-check
    (n_dofs=8, ffp1 geometry) shows sigma_min at the OLD 15-matched-mode
    Omega_native values shifts by a non-trivial amount at that basis size
    -- the roots likely move and the P2_EASY_STATUS.md Sec 5.2 "15 matched
    C-C IP Hz" table is NOT assumed unchanged. A cluster re-search of the
    C-C IP window against this fix, re-scored against cc_ip_mesh64/96.txt
    FE data the same way as Sec 5.2, is required before quoting a revised
    (or confirmed-unchanged) C-C IP table. Scout only until that job runs
    AND the reverse-cantilever identity (FreeClampedIP) is cluster-
    confirmed.
    """
    token = "clamped_clamped"
    name = "clamped(-Theta)-clamped(+Theta)"
    validated = False
    def __init__(self):
        super().__init__([EdgeSpec(+1, EdgeKind.CLAMPED), EdgeSpec(-1, EdgeKind.CLAMPED)])

_BC_REGISTRY = {
    ("clamped_free", "out_of_plane"): ClampedFreeOOP,
    ("clamped_free", "in_plane"):     ClampedFreeIP,
    ("free_free",    "out_of_plane"): FreeFreeOOP,
    ("free_free",    "in_plane"):     FreeFreeIP,
    ("free_clamped", "out_of_plane"): FreeClampedOOP,
    ("free_clamped", "in_plane"):     FreeClampedIP,
    ("clamped_clamped", "out_of_plane"): ClampedClampedOOP,
    ("clamped_clamped", "in_plane"):     ClampedClampedIP,
}

def make_bc(kind, motion):
    """Factory: kind in {'clamped_free','free_free','free_clamped',
    'clamped_clamped'}, motion in {'out_of_plane','in_plane'}.  Lets a BC
    token be threaded through worker scalars (Step 6) so parallel
    reconstruction picks the right BC."""
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


