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
  * TestWorkerPath   -- multiprocessing worker-side scalar reconstruction
    (fast, no ProcessPoolExecutor) is always run; the ProcessPoolExecutor
    end-to-end variant is skipped by default (SLOW_TESTS=1 to include it)
    because even a narrow window can exceed sandbox time, per the skill.

Run with:  python -m unittest tests.test_solver -v
Slow tier: SLOW_TESTS=1 python -m unittest tests.test_solver -v
"""
from __future__ import annotations
import os
import sys
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from mpmath import mp
import plate_solver as ps
from plate_solver.detectors import (
    _selftest_part1, _selftest_part2, full_search, select_fill,
    sigma_min_from_K, find_modes_sigmin,
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
    """The four known-good sigma_min(K) anchors from the sandbox-probe
    skill. Each is a SINGLE-POINT probe (full_search + _build_K_real +
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
        self.assertAlmostEqual(val, -1.418728, places=5)

    def test_freefree_oop_gate(self):
        geom, mat = _canonical_geom_mat()
        oop = ps.OutOfPlaneSolver(geom, mat, M=80, n_quad=30, boundary=FreeFreeOOP())
        val = self._sigma_min_at(oop, 0.12)
        self.assertAlmostEqual(val, -2.553504, places=5)

    def test_ip_freefree_gate(self):
        geom, mat = _canonical_geom_mat()
        ip = ps.InPlaneSolver(geom, mat, M=80, n_quad=30, boundary=FreeFreeIP())
        val = self._sigma_min_at(ip, 0.12)
        self.assertAlmostEqual(val, -1.1217, places=3)


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


if __name__ == "__main__":
    unittest.main()
