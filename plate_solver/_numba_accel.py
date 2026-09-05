# -*- coding: utf-8 -*-
"""
plate_solver._numba_accel -- optional, opt-in Numba acceleration scaffold.

STATUS: infrastructure only, NOT wired into the float64 branch trackers.

The review doc's suggestion ("explore JIT compilation via Numba for the fast
float64 branch trackers") is directionally reasonable -- ExactEdgeSolver.series
/ .newton in dispersion.py are exactly the kind of tight complex128 loop that
benefits from Numba. But actually wiring it in requires:

  1. Refactoring those methods from bound methods reading `self.M`,
     `self.r0_bar`, etc. into free functions with explicit, Numba-typeable
     arguments (nopython mode cannot JIT arbitrary Python objects/self
     attribute lookups) -- a real code change to the numerical hot path,
     not a decorator drop-in.
  2. Verifying, numerically, that the JIT'd version reproduces the existing
     float64 branch seeds bit-for-bit (or close enough not to perturb which
     branches get tracked/selected) -- exactly the kind of claim this
     project's own rules (LESSONS_LEARNED, the sandbox-probe skill) say must
     be checked by direct numerical probe, never by reasoning alone.

This sandbox has no network access and Numba is not installed here (checked
2026-07-01: `pip install numba` fails with "no matching distribution", and
`import numba` fails with ModuleNotFoundError) -- so step 2 above cannot be
done in this environment at all. Shipping a refactored hot path that was
never actually run would be exactly the "hopeful attempt that fails later"
this project's own guidance warns against, so it was deliberately NOT done.

What IS here: a safe, opt-in decorator so that WHEN a Numba-enabled
environment (i.e. the cluster, if numba is added there) wants to try this,
the plumbing exists and degrades to a plain no-op everywhere else --
importing this module never fails and never changes behaviour unless
R40_NUMBA=1 AND numba is actually importable.

To finish this properly on the cluster: pip install numba there, refactor
ExactEdgeSolver.series/newton (or a hand-picked hot subroutine of them) into
free functions decorated with `@maybe_jit(nopython=True)`, then run the
sandbox-probe skill's known-good branch set at Ω̃≈0.033 (see
plate-solver-sandbox-probe skill, "Known-good probe targets") through BOTH
the JIT'd and un-JIT'd path and diff the branch roots before trusting it.
"""
from __future__ import annotations
import os

NUMBA_REQUESTED = os.environ.get("R40_NUMBA", "0") == "1"

try:
    import numba as _numba
    NUMBA_AVAILABLE = True
except ImportError:
    _numba = None
    NUMBA_AVAILABLE = False

NUMBA_ACTIVE = NUMBA_REQUESTED and NUMBA_AVAILABLE

if NUMBA_REQUESTED and not NUMBA_AVAILABLE:
    print("  [plate_solver] R40_NUMBA=1 was set but `numba` is not "
          "importable in this environment; continuing with plain Python/"
          "numpy (no acceleration, no behaviour change).")


def maybe_jit(*jit_args, **jit_kwargs):
    """Decorator: apply numba.njit(*jit_args, **jit_kwargs) IFF R40_NUMBA=1
    and numba is importable; otherwise return the function unchanged.
    Nothing in this package currently uses this on a hot path -- see the
    module docstring for why, and what's needed to change that safely."""
    def deco(fn):
        if NUMBA_ACTIVE:
            return _numba.njit(*jit_args, **jit_kwargs)(fn)
        return fn
    return deco
