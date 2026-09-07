# -*- coding: utf-8 -*-
"""
Unit tests for plate_solver.

Uses stdlib unittest (NOT pytest): this sandbox has no network access and
pytest is not installed here (`pip install pytest` fails -- no matching
distribution). unittest is stdlib and needs nothing extra, so `python -m
unittest discover` works in any environment including a bare cluster node.

Test tiers, cheapest first (see plate-solver-sandbox-probe skill for why the
sandbox can't run more than this):

  * TestSelfTests    -- ODE residual self-tests (~1s, dps=30)
  * TestCutoffs      -- Bessel cut-off frequency computation (~1s)
  * TestRegressionGates -- the four known-good sigma_min(K) anchors from the
    sandbox-probe skill, at dps=26 (~1-2 min total; single-point probes only,
    exactly the pattern the skill mandates -- NOT a full find_modes_sigmin
    sweep, which the sandbox cannot finish, see TestWorkerPath below).
  * TestRectangularFreeFree -- rectangular FFFF: two single-point
    sigma_min anchors (out-of-plane and in-plane, ~2s total) plus a
    source-level guard that the four-corner Kirchhoff jump is the
    checkerboard and not the naive same-sign sum that vanishes.
  * TestWorkerPath   -- multiprocessing worker-side scalar reconstruction
    (fast, no ProcessPoolExecutor) is always run; the ProcessPoolExecutor
    end-to-end variant is skipped by default (SLOW_TESTS=1 to include it)
    because even a narrow window can exceed sandbox time, per the skill.

Run with:  python -m unittest tests.test_solver -v
Slow tier: SLOW_TESTS=1 python -m unittest tests.test_solver -v

The summary line's test count is interpreter-dependent and is NOT a
completeness signal. Python 3.11 counts the one @skipUnless test in
"Ran N"; Python 3.12 does not, so the same 14 collected tests report
"Ran 14" and "Ran 13" respectively, both with skipped=1 (measured
2026-09-07 on 3.11.15 and 3.12.1). To ask whether every gate is actually
present on a given machine, run diag_test_collection.sh at the package
root -- it enumerates what the loader collects -- rather than reading the
count off this summary.

ANCHOR PROVENANCE (all three annular gates re-derived 2026-09-07). Every
value below reproduces to ~15 significant digits on two independent
environments -- scipy 1.17.1 / numpy 2.4.4 / mpmath 1.3.0 and 1.4.1, and the
Mill cluster venv at scipy 1.18.0 / numpy 2.5.1 / mpmath 1.4.1 -- so the
generous `places` tolerances have several orders of magnitude of headroom and
these are not library-sensitive.

Until 2026-09-07 all three asserted stale values, in two different ways, and
the suite had not been run since. Both are recorded in LESSONS_LEARNED
Sec. 1 / Sec. 16 / Sec. 18.34 / Sec. 18.37:

  * The two free-free gates asserted -2.553504 and -1.1217, which Sec. 1
    labels verbatim as "PRE-fix historical values (any un-patched tree)".
    They were corrected once already, on 2026-07-09, and had drifted back.
    The post-fix values -2.2307545049 and -0.8946416553 are what the solver
    has produced since 2026-07-08.s5.
  * The cantilever gate asserted -1.418728, captured at 2026-07-10.s8 and
    annotated "must NEVER move". It moved deliberately: the s9->s10 bump was
    the E.2/Geo-2 deployment of 2026-07-28, which changed
    OutOfPlaneSolver._build_K_real's CLAMPED orientation -- and CLAMPED is
    what a cantilever has. Paper 1's Tables 3 and 4 were regenerated at s10
    for exactly this reason (job 2443456, all 29 rows moved). The two
    free-free gates, which have no clamped edge, did NOT move, which is the
    signature that deployment predicts.

Do NOT edit a stored anchor merely to make a gate pass. Change one only when
an independent record says what the new value should be and why, as above.
"""
from __future__ import annotations
import os
import sys
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from mpmath import mp, mpf
import plate_solver as ps
from plate_solver.detectors import (
    _selftest_part1, _selftest_part2, full_search, select_fill,
    sigma_min_from_K, find_modes_sigmin,
    rect_resolve_branches, rect_select_branches,
    rect_ip_resolve_branches, rect_ip_select_branches,
)
from plate_solver.boundary import FreeFreeOOP, FreeFreeIP
from plate_solver.dispersion import cutoff_frequencies_part1, cutoff_frequencies_part2
from plate_solver.geometry import _make_geom_mat, MaterialModel

SLOW = os.environ.get("SLOW_TESTS", "0") == "1"


def _canonical_geom_mat():
    geom = ps.make_geometry(1.25, 1.0)
    mat = ps.IsotropicMaterial(E=210e9, nu=0.35, rho=7800.0)
    return geom, mat


class TestSelfTests(unittest.TestCase):
    """ODE-residual self-tests. Reference values from LESSONS_LEARNED /
    plate-solver-sandbox-probe skill: Part 1 residual ~2.42e-13, Part 2
    ~1.95e-15 (both must read 'OK', i.e. < 1e-6)."""

    def setUp(self):
        mp.dps = 30

    def test_part1_ode_residual(self):
        geom, mat = _canonical_geom_mat()
        self.assertTrue(_selftest_part1(geom, mat))

    def test_part2_ode_residual(self):
        geom, mat = _canonical_geom_mat()
        self.assertTrue(_selftest_part2(geom, mat))


class TestCutoffs(unittest.TestCase):
    """Cut-off frequency computation (ξ=0 / ζ=0 Bessel determinants). Only
    checks basic sanity (positive, sorted, right count) -- these values are
    geometry/material-derived, not paper-tabulated, so there's no external
    reference to gate against; the real check is that cutoff_frequencies_*
    still runs after the float-cast hoist (see dispersion.py)."""

    def test_part1_cutoffs_isotropic(self):
        r0b = float(ps.make_geometry(1.25, 1.0).r0_bar)
        cutoffs = cutoff_frequencies_part1(r0b, 0.35, T=1.0, R=1.0, n_co=4)
        self.assertGreaterEqual(len(cutoffs), 1)
        self.assertTrue(all(c > 0 for c in cutoffs))
        self.assertEqual(cutoffs, sorted(cutoffs))

    def test_part2_cutoffs_isotropic(self):
        r0b = float(ps.make_geometry(1.25, 1.0).r0_bar)
        c11_eff = 2.0 / (1 - 0.35)
        cutoffs = cutoff_frequencies_part2(r0b, 0.35, c11_eff, R=1.0, n_co=4)
        self.assertGreaterEqual(len(cutoffs), 1)
        self.assertTrue(all(c > 0 for c in cutoffs))

    def test_part1_cutoffs_accept_mpf_inputs(self):
        """The float-cast hoist (dispersion.py) must not break callers that
        still pass mpf/mpmath scalars instead of plain floats."""
        from mpmath import mpf
        r0b = ps.make_geometry(1.25, 1.0).r0_bar   # mpf, not float
        cutoffs = cutoff_frequencies_part1(r0b, mpf("0.35"), T=mpf(1), R=mpf(1), n_co=3)
        self.assertGreaterEqual(len(cutoffs), 1)


class TestRegressionGates(unittest.TestCase):
    """The known-good sigma_min(K) anchors of LESSONS_LEARNED Sec. 1 / Sec. 16,
    at SOLVER_VERSION 2026-07-10.s10 (see the module docstring for provenance).
    Three of the four listed there are exercised here; the fourth, the FF-P1
    geometry point -2.0383, needs a different geometry and is not included. Each is a SINGLE-POINT probe (full_search + _build_K_real +
    sigma_min_from_K at one Omega), matching the skill's mandated pattern --
    not a full mode search. dps=26 per the skill's sandbox timing guidance."""

    def setUp(self):
        mp.dps = 26

    def _sigma_min_at(self, solver, Om, n_dofs=16):
        brs = full_search(solver.fast, Om, xmax=14.0)
        sel, n = select_fill(brs, n_dofs)
        K, sz = solver._build_K_real(Om, sel, fast_scan=False)
        return sigma_min_from_K(K, sz)

    def test_cantilever_oop_gate(self):
        geom, mat = _canonical_geom_mat()
        oop = ps.OutOfPlaneSolver(geom, mat, M=80, n_quad=30)
        val = self._sigma_min_at(oop, 0.12)
        self.assertAlmostEqual(val, -1.6497881851, places=5)

    def test_freefree_oop_gate(self):
        geom, mat = _canonical_geom_mat()
        oop = ps.OutOfPlaneSolver(geom, mat, M=80, n_quad=30, boundary=FreeFreeOOP())
        val = self._sigma_min_at(oop, 0.12)
        self.assertAlmostEqual(val, -2.2307545049, places=5)

    def test_ip_freefree_gate(self):
        geom, mat = _canonical_geom_mat()
        ip = ps.InPlaneSolver(geom, mat, M=80, n_quad=30, boundary=FreeFreeIP())
        val = self._sigma_min_at(ip, 0.12)
        self.assertAlmostEqual(val, -0.8946416553, places=5)


class TestWorkerPath(unittest.TestCase):
    """Worker-side scalar reconstruction (_make_geom_mat / _ScalarMat) must
    subclass MaterialModel and reproduce the in-process material's oop/ip
    constants -- see LESSONS_LEARNED's repeated worker-reconstruction-path
    warning and the sandbox-probe skill Sec. 4. This is the fast, always-run
    check; the full ProcessPoolExecutor variant is opt-in (SLOW_TESTS=1)."""

    def test_scalar_material_matches_in_process(self):
        geom, mat = _canonical_geom_mat()
        g_scalars = (float(geom.r0_bar), float(geom.Theta), float(geom.b),
                     float(geom.r_0), float(geom.R_i), float(geom.R_o), float(geom.h))
        m_scalars = (float(mat.nu_bar), float(mat.T), float(mat.R),
                     float(mat.c11_bar), float(mat.c11_eff_bar), float(mat.c66),
                     float(mat.rho))
        g, m = _make_geom_mat(g_scalars, m_scalars, 26)
        self.assertIsInstance(m, MaterialModel)
        self.assertFalse(m.radially_graded)
        oop_a, oop_b = mat.oop_constants(), m.oop_constants()
        for a, b in zip(oop_a, oop_b):
            self.assertAlmostEqual(float(a), float(b), places=9)

    @unittest.skipUnless(SLOW, "set SLOW_TESTS=1 to run the ProcessPoolExecutor gate "
                                "(can exceed a constrained sandbox's time budget; "
                                "always safe/fast on the cluster)")
    def test_find_modes_sigmin_worker_pool(self):
        geom, mat = _canonical_geom_mat()
        mp.dps = 26
        oop = ps.OutOfPlaneSolver(geom, mat, M=80, n_quad=30)
        freqs = find_modes_sigmin(oop, 1, (0.028, 0.038), n_scan=6, n_dofs=16,
                                   max_dim=16, n_modes_wanted=1, verbose=False,
                                   n_workers=2)
        self.assertTrue(len(freqs) >= 1)
        self.assertAlmostEqual(freqs[0], 0.033238, places=2)


class TestRectangularFreeFree(unittest.TestCase):
    """Rectangular completely-free (FFFF) gates for the rectangular companion
    paper: two single-point sigma_min anchors plus a source-level guard on the
    four-corner Kirchhoff jump.

    The source guard is not redundant with the numeric anchors -- it names the
    specific failure it is guarding against. Before 2026-09-03 the free-free
    corner term was assembled as the naive same-sign four-term sum
    TT + TB + WT + WB. That sum vanishes identically under the parity identity
    (Supplementary S.1.3), so it assembles without raising and silently yields
    the pre-corner spectrum: at the anchor below it moves sigma_min from
    2.3836e-04 to 3.9507e-03, a factor of 16.6. A released tree carrying the
    old formula therefore looks healthy and reproduces nothing.

    Both anchors are stable to every printed digit at dps = 26, 30 and 40, and
    each costs well under a second, so this class adds ~2 s to the suite.
    """

    CHECKERBOARD = "p1t * Ut1 * q0t * Vt0 - p1b * Ut1 * q0b * Vt0"
    NAIVE_SUM = "p1t * Ut1 * q0t * Vt0 + p1b * Ut1 * q0b * Vt0"

    def setUp(self):
        mp.dps = 26

    @staticmethod
    def _sigma(asm, sym, lob, x, n_real, n_cpair, im_cap, in_plane=False):
        """Single-point sigma_min through the production path: resolve
        branches, select the basis, assemble, take the equilibrated smallest
        singular value. Mirrors p5_rect_ff_lib.sigma_at / sigma_at_ip."""
        if in_plane:
            reps = rect_ip_resolve_branches(asm.eng(sym), x, im_cap=im_cap)
            full = rect_ip_select_branches(reps, n_real=n_real,
                                           n_cpair=n_cpair)
        else:
            reps = rect_resolve_branches(asm.eng(sym), x, im_cap=im_cap)
            full = rect_select_branches(reps, n_real=n_real, n_cpair=n_cpair)
        K = asm.assemble(full, mpf(str(round(x, 6))), sym, lob)
        return float(asm.equil_sigma(K, full)), len(full)

    def test_free_free_corner_is_the_checkerboard(self):
        import plate_solver.core_solvers as cs
        with open(cs.__file__, encoding="utf-8") as fh:
            src = fh.read()
        self.assertTrue(
            self.CHECKERBOARD in src,
            "free_free corner_ff is not the closed-contour checkerboard "
            "(TT - TB - WT + WB); this tree cannot reproduce the "
            "rectangular FFFF flexural tables.")
        self.assertFalse(
            self.NAIVE_SUM in src,
            "free_free corner_ff is the naive same-sign four-term sum, which "
            "vanishes identically under the parity identity; this is the "
            "pre-2026-09-03 assembler.")

    def test_rect_ff_oop_gate(self):
        """l/b = 1.5, SYM, Lambda* = 0.964 (first symmetric flexural mode,
        production basis n_real=6 / n_cpair=0)."""
        asm = ps.RectOOPAssembler(nu_b=0.3, R=1.0, T=1.0, dps=mp.dps,
                                  bc="free_free")
        sig, n = self._sigma(asm, True, 1.5, 0.964, 6, 0, 30.0)
        self.assertEqual(n, 3)
        self.assertAlmostEqual(sig / 2.3835918992545717e-04, 1.0, places=9)

    def test_rect_ff_ip_gate(self):
        """l/b = 2.0, ANTI, Omega-bar* = 0.5200 (the Bardell-matched first
        antisymmetric extensional mode, basis n_real=5 / n_cpair=3)."""
        asm = ps.RectIPAssembler(nu_b=0.3, R=1.0, dps=mp.dps, bc="free_free")
        sig, n = self._sigma(asm, False, 2.0, 0.5200, 5, 3, 7.0,
                             in_plane=True)
        self.assertEqual(n, 4)
        self.assertAlmostEqual(sig / 3.900364925307852e-08, 1.0, places=9)

    def test_rect_ip_bc_switch_is_live(self):
        """The in-plane free-free option must exist and must actually change
        the assembly. A tree whose RectIPAssembler has no bc argument raises
        TypeError here rather than silently returning cantilever numbers."""
        ff = ps.RectIPAssembler(nu_b=0.3, R=1.0, dps=mp.dps, bc="free_free")
        cf = ps.RectIPAssembler(nu_b=0.3, R=1.0, dps=mp.dps,
                                bc="clamped_free")
        s_ff, _ = self._sigma(ff, False, 2.0, 0.5200, 5, 3, 7.0,
                              in_plane=True)
        s_cf, _ = self._sigma(cf, False, 2.0, 0.5200, 5, 3, 7.0,
                              in_plane=True)
        self.assertGreater(abs(s_cf - s_ff) / s_cf, 0.5)
        with self.assertRaises(ValueError):
            ps.RectIPAssembler(nu_b=0.3, R=1.0, dps=mp.dps, bc="nonsense")

if __name__ == "__main__":
    unittest.main()
