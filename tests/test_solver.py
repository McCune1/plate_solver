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
  * TestRingDisk -- Paper 3 Phase 0/D: disk path is 2x2; Sec. 18.38
    F-F anchors; Phase D F-C/C-F vs Table 2.16; SS still raises.
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
import math
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
from plate_solver.ring_disk import (
    DISK_4X4_AT_RI0_IS_NOT_THE_DISK,
    FFP1_MODE7_ANSYS_LIT, FFP1_MODE7_NATIVE,
    SECTION_18_38_N0_LOGABSDET, SECTION_18_38_N0_OM,
    SECTION_18_38_N2_LOGABSDET, SECTION_18_38_N2_OM,
    UnsupportedRingBC,
    disk_L_mp, flexural_lambda2, illegal_ri0_4x4_rank, inplane_lambda_irie,
    make_annulus_solver, make_disk_solver, native_from_flexural_lambda2,
    ring_L_mp, ring_logabsdet, ring_search,
)

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

class TestRingDisk(unittest.TestCase):
    """Closed-ring / solid-disk gates for Paper 3 Phase 0 and Phase D.

    Pre-registered verdicts (do not retune M, dps, or the screen to make
    a row pass):

      PASS              all checks below hold
      FAIL_DISK_IS_4x4  disk path is the 4x4 at R_i=0, or disk_L_mp is not 2x2
      FAIL_ANCHOR       a Sec. 18.38 |det L| anchor moved -- do not edit it
      FAIL_IP_LMAT      in-plane `_Lmat_mp` at integer n is not finite
      FAIL_PHASE_D      mixed F-C/C-F miss Table 2.16, or SS no longer raises

    Anchors are identical at dps 26/30/40 (captured 2026-09-08). If one
    ever moves, that is a changed computed frequency: do not edit the
    stored value.
    """

    def test_disk_path_is_2x2_not_4x4_at_ri0(self):
        import plate_solver.ring_disk as rd
        with open(rd.__file__, encoding="utf-8") as fh:
            src = fh.read()
        self.assertTrue(
            DISK_4X4_AT_RI0_IS_NOT_THE_DISK,
            "FAIL_DISK_IS_4x4: the rank-3 job-2339453 flag was flipped")
        self.assertIn("DISK_4X4_AT_RI0_IS_NOT_THE_DISK", src)
        self.assertIn("matrix(2, 2)", src)
        self.assertIn("is NOT the disk", src)
        self.assertNotRegex(
            src,
            r"def disk_det_mp[\s\S]{0,400}solver\._det_mp",
            "FAIL_DISK_IS_4x4: disk_det_mp calls the 4x4 _det_mp")
        mp.dps = 26
        solver, _g, _m = make_disk_solver(nu=0.30, motion="oop")
        self.assertEqual(float(solver.geom.R_i), 0.0)
        L2 = disk_L_mp(solver, 2, 0.5)
        self.assertEqual(L2.rows, 2, "FAIL_DISK_IS_4x4")
        self.assertEqual(L2.cols, 2, "FAIL_DISK_IS_4x4")
        rank0, _ = illegal_ri0_4x4_rank(solver, 0, 0.5)
        rank1, _ = illegal_ri0_4x4_rank(solver, 1, 0.5)
        self.assertEqual(rank0, 2, "job 2339453: 4x4 at R_i=0 is rank 2 at n=0")
        self.assertEqual(
            rank1, 3,
            "FAIL_DISK_IS_4x4 / job 2339453: 4x4 at R_i=0 is rank 3 "
            "(not 2) at n>=1; do not 'fix' this")

    def test_ring_oop_ff_n2_anchor(self):
        """OOP F-F n=2, Omega* = 0.1081885366 (Sec. 18.38 confirmed pair)."""
        for dps in (26, 30, 40):
            mp.dps = dps
            solver, _g, _m = make_annulus_solver(0.5, nu=0.30, motion="oop")
            logd = ring_logabsdet(solver, 2, SECTION_18_38_N2_OM)
            self.assertAlmostEqual(
                logd / SECTION_18_38_N2_LOGABSDET, 1.0, places=9,
                msg="FAIL_ANCHOR n=2 at dps=%s: log10|det|=%s stored=%s"
                    % (dps, logd, SECTION_18_38_N2_LOGABSDET))

    def test_ring_oop_ff_n0_anchor(self):
        """OOP F-F n=0, Omega* = 0.2359132157 (Sec. 18.38 confirmed pair)."""
        for dps in (26, 30, 40):
            mp.dps = dps
            solver, _g, _m = make_annulus_solver(0.5, nu=0.30, motion="oop")
            logd = ring_logabsdet(solver, 0, SECTION_18_38_N0_OM)
            self.assertAlmostEqual(
                logd / SECTION_18_38_N0_LOGABSDET, 1.0, places=9,
                msg="FAIL_ANCHOR n=0 at dps=%s: log10|det|=%s stored=%s"
                    % (dps, logd, SECTION_18_38_N0_LOGABSDET))

    def test_ip_lmat_live_at_integer_n(self):
        """In-plane `_Lmat_mp` at integer n is finite. Not a full IP table."""
        mp.dps = 26
        solver, _g, _m = make_annulus_solver(0.5, nu=0.30, motion="ip")
        L = ring_L_mp(solver, 2, 0.5)
        self.assertEqual(L.rows, 4)
        self.assertEqual(L.cols, 4)
        for i in range(4):
            for j in range(4):
                z = complex(L[i, j])
                self.assertTrue(
                    math.isfinite(z.real) and math.isfinite(z.imag),
                    "FAIL_IP_LMAT: L[%s,%s]=%s" % (i, j, z))

    def test_flexural_conversion_ffp1_mode7(self):
        """FF-P1 mode 7 native 1.369611 -> Ansys Omega_lit 54.0755 at ~0.01%.

        The factor 39.478... is geometry-specific (job 2406948), not a
        universal constant: b/a=0.3 must not reuse it. IP Irie lambda is
        a different symbol and must not reuse the flexural factor.
        """
        import math as _math
        _solver, geom, mat = make_annulus_solver(0.5, nu=0.30, motion="oop")
        om_lit = flexural_lambda2(FFP1_MODE7_NATIVE, geom, mat)
        rel = abs(om_lit - FFP1_MODE7_ANSYS_LIT) / FFP1_MODE7_ANSYS_LIT
        self.assertLess(rel, 2.0e-4,  # 0.02%; measured 0.010%
                        "FF-P1 mode 7 conversion %s vs Ansys %s rel=%s"
                        % (om_lit, FFP1_MODE7_ANSYS_LIT, rel))
        factor_ba05 = flexural_lambda2(1.0, geom, mat)
        _s3, geom3, mat3 = make_annulus_solver(0.3, nu=0.30, motion="oop")
        factor_ba03 = flexural_lambda2(1.0, geom3, mat3)
        self.assertGreater(abs(factor_ba05 - factor_ba03) / factor_ba05, 0.3)
        self.assertAlmostEqual(factor_ba05 / (4.0 * _math.pi ** 2), 1.0,
                               places=9)
        lam_ip = inplane_lambda_irie(1.0, geom, mat)
        self.assertGreater(lam_ip, 0.0)
        self.assertGreater(abs(lam_ip - factor_ba05) / factor_ba05, 0.5)

    def test_ss_and_ip_mixed_still_unsupported(self):
        """SS and IP mixed-edge stay Phase F / not Phase D."""
        mp.dps = 26
        slv, _g, _m = make_annulus_solver(0.5, nu=0.30, motion="oop")
        with self.assertRaises(UnsupportedRingBC):
            ring_L_mp(slv, 2, 0.5, bc_inner="S", bc_outer="F")
        with self.assertRaises(UnsupportedRingBC):
            ring_L_mp(slv, 2, 0.5, bc_inner="F", bc_outer="SS")
        slv_ip, _gi, _mi = make_annulus_solver(0.5, nu=0.30, motion="ip")
        with self.assertRaises(UnsupportedRingBC):
            ring_L_mp(slv_ip, 2, 0.5, bc_inner="C", bc_outer="F")

    def test_fc_cf_lmat_is_4x4_finite(self):
        """Phase D mixed L is 4x4 and finite; F-F L is unchanged in shape."""
        mp.dps = 26
        slv, _g, _m = make_annulus_solver(0.5, nu=0.30, motion="oop")
        Lff = ring_L_mp(slv, 2, 0.5, "F", "F")
        self.assertEqual((Lff.rows, Lff.cols), (4, 4))
        for bi, bo in (("F", "C"), ("C", "F"), ("C", "C")):
            L = ring_L_mp(slv, 2, 0.5, bi, bo)
            self.assertEqual((L.rows, L.cols), (4, 4), (bi, bo))
            for i in range(4):
                for j in range(4):
                    z = complex(L[i, j])
                    self.assertTrue(
                        math.isfinite(z.real) and math.isfinite(z.imag),
                        "L[%s,%s]=%s for %s-%s" % (i, j, z, bi, bo))

    def test_phase_d_table216_axisym_fc_cf(self):
        """Table 2.16 exact n=0 b/a=0.5 nu=1/3: F-C 17.51 and C-F 13.05.

        F-C = inner F, outer C (Leissa 'Clamped, Free'). C-F = inner C,
        outer F (Leissa 'Free, Clamped'). Rel < 1% vs the exact column.
        """
        mp.dps = 26
        slv, geom, mat = make_annulus_solver(0.5, nu=1.0 / 3.0, motion="oop")
        for inner, outer, printed in (("F", "C", 17.51), ("C", "F", 13.05)):
            native = native_from_flexural_lambda2(printed, geom, mat)
            res = ring_search(
                0, (max(1e-4, native - 0.03), native + 0.03),
                bc_inner=inner, bc_outer=outer, motion="oop",
                solver=slv, geom=geom, mat=mat,
                coarse_step=0.0005, polish=True)
            rel = abs(res.lambda2 - printed) / printed
            self.assertLess(
                rel, 0.01,
                "Table 2.16 %s-%s printed=%s lam2=%s rel=%s log|det|=%s"
                % (inner, outer, printed, res.lambda2, rel, res.log_abs_det))

if __name__ == "__main__":
    unittest.main()
