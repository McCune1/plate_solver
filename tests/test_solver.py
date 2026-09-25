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
  * TestRingDisk -- Paper 3 Phase 0/D/F: disk path is 2x2 (4x4 at R_i=0
    is not the disk); two full-ring OOP F-F |det L| anchors from Sec.
    18.38; in-plane `_Lmat_mp` liveness at integer n; flexural conversion
    of FF-P1 mode 7; Phase D F-C/C-F and Phase F C-C/S-S vs Table 2.16.
    Cheap 4x4/2x2 probes, not a K-build.
  * TestWorkerPath   -- multiprocessing worker-side scalar reconstruction
    (fast, no ProcessPoolExecutor) is always run; the ProcessPoolExecutor
    end-to-end variant is skipped by default (SLOW_TESTS=1 to include it)
    because even a narrow window can exceed sandbox time, per the skill.

Run with:  python -m unittest tests.test_solver -v
Slow tier: SLOW_TESTS=1 python -m unittest tests.test_solver -v

The summary line's test count is interpreter-dependent and is NOT a
completeness signal. Python 3.11 counts the one @skipUnless test in
"Ran N"; Python 3.12 does not. After TestRingDisk (2026-09-08) the
loader collects 22 tests (21 run + 1 skip on 3.14). To ask whether every gate is actually present on a
given machine, run diag_test_collection.sh at the package root -- it
enumerates what the loader collects -- rather than reading the count
off this summary.

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
from plate_solver.boundary import (
    FreeFreeOOP, FreeFreeIP, FreeClampedOOP, FreeClampedIP,
    ClampedClampedOOP, ClampedClampedIP, make_bc, EdgeKind,
)
from plate_solver.dispersion import cutoff_frequencies_part1, cutoff_frequencies_part2
from plate_solver.geometry import _make_geom_mat, MaterialModel, \
    _kbar_wbar, _omega_lit
from plate_solver.piezo_solver import PiezoOutOfPlaneSolver
from plate_solver.piezo_disk import (
    PiezoDiskSC, DiskPathError as PiezoDiskPathError,
    make_liu_disk, LIU_DISK_KWARGS,
)
from plate_solver.piezo_monolithic import PiezoMonolithicOutOfPlaneSolver
from plate_solver.workers import (
    _piezo_elastic_ff_root_worker,
    _piezo_coupled_ff_root_worker,
    _piezo_elastic_cc_root_worker,
    _piezo_coupled_cc_root_worker,
    _piezo_oc_ff_root_worker,
    _piezo_oc_sine_ff_root_worker,
    _piezo_elastic_cf_root_worker,
    _piezo_oc_cf_root_worker,
    _piezo_mono_elastic_ff_root_worker,
    _piezo_mono_coupled_ff_root_worker,
    _piezo_mono_elastic_cc_root_worker,
    _piezo_mono_coupled_cc_root_worker,
    _piezo_mono_elastic_cf_root_worker,
    _piezo_mono_elastic_fc_root_worker,
    _piezo_mono_coupled_cf_root_worker,
    _piezo_mono_coupled_fc_root_worker,
    _piezo_mono_driven_force_sc_worker,
)
from plate_solver.ring_disk import (
    DISK_4X4_AT_RI0_IS_NOT_THE_DISK,
    FFP1_MODE7_LAM2, FFP1_MODE7_NATIVE,
    SECTION_18_38_N0_LOGABSDET, SECTION_18_38_N0_OM,
    SECTION_18_38_N2_LOGABSDET, SECTION_18_38_N2_OM,
    UnsupportedRingBC,
    disk_L_mp, flexural_lambda2, illegal_ri0_4x4_rank, inplane_lambda_irie,
    make_annulus_solver, make_disk_solver, native_from_flexural_lambda2,
    ring_L_mp, ring_det_mp, ring_logabsdet, ring_search,
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


class TestP2EasyBoundary(unittest.TestCase):
    """Registry and edge-list liveness for reverse-cantilever and
    clamped-clamped. No K-build: the identity/C-C frequencies are
    cluster probes. Does not bump SOLVER_VERSION and does not touch
    the three regression-gate numbers above."""

    def test_make_bc_new_tokens(self):
        for kind in ("free_clamped", "clamped_clamped"):
            for motion in ("out_of_plane", "in_plane"):
                bc = make_bc(kind, motion)
                self.assertEqual(bc.token, kind)
                self.assertFalse(bc.validated)
                self.assertEqual(len(bc.edges), 2)

    def test_reverse_swaps_the_cantilever_walls(self):
        rev = FreeClampedOOP()
        kinds = [(e.sign, e.kind) for e in rev.edges]
        self.assertEqual(kinds, [(+1, EdgeKind.CLAMPED), (-1, EdgeKind.FREE)])
        self.assertEqual(make_bc("free_clamped", "in_plane").token, "free_clamped")

    def test_cc_both_walls_clamped(self):
        cc = ClampedClampedOOP()
        self.assertTrue(all(e.kind == EdgeKind.CLAMPED for e in cc.edges))
        self.assertEqual(len(cc.clamped_edges()), 2)
        ip = ClampedClampedIP()
        self.assertEqual(len(ip.clamped_edges()), 2)

    def test_legacy_tokens_unchanged(self):
        self.assertEqual(make_bc("clamped_free", "out_of_plane").token, "clamped_free")
        self.assertEqual(make_bc("free_free", "in_plane").token, "free_free")
        with self.assertRaises(ValueError):
            make_bc("no_such_bc", "out_of_plane")


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
    """Closed-ring / solid-disk gates for Paper 3 Phase 0, D and F.

    Pre-registered verdicts (do not retune M, dps, or the screen to make
    a row pass):

      PASS              all checks below hold
      FAIL_DISK_IS_4x4  disk path is the 4x4 at R_i=0, or disk_L_mp is not 2x2
      FAIL_ANCHOR       a Sec. 18.38 |det L| anchor moved -- do not edit it
      FAIL_IP_LMAT      in-plane `_Lmat_mp` at integer n is not finite
      FAIL_PHASE_D      mixed F-C/C-F miss Table 2.16, or IP C/F L is not 4x4
      FAIL_PHASE_F      C-C or S-S miss Table 2.16 exact n=0 b/a=0.5 nu=1/3

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
        for bo in ("C", "S"):
            Lcs = disk_L_mp(solver, 2, 0.5, bc_outer=bo)
            self.assertEqual((Lcs.rows, Lcs.cols), (2, 2), bo)
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
        """FF-P1 mode 7 native 1.369611 -> Kirchhoff lambda^2 54.070.

        Paper 1 Table 6 prints FEM 53.828 for this mode (0.45%). That
        finite-element residual is not this check. The factor 39.478...
        is geometry-specific (job 2406948), not a universal constant:
        b/a=0.3 must not reuse it. IP Irie lambda is a different symbol
        and must not reuse the flexural factor.
        """
        import math as _math
        _solver, geom, mat = make_annulus_solver(0.5, nu=0.30, motion="oop")
        om_lit = flexural_lambda2(FFP1_MODE7_NATIVE, geom, mat)
        rel = abs(om_lit - FFP1_MODE7_LAM2) / FFP1_MODE7_LAM2
        self.assertLess(rel, 2.0e-5,
                        "FF-P1 mode 7 conversion %s vs Kirchhoff %s rel=%s"
                        % (om_lit, FFP1_MODE7_LAM2, rel))
        factor_ba05 = flexural_lambda2(1.0, geom, mat)
        _s3, geom3, mat3 = make_annulus_solver(0.3, nu=0.30, motion="oop")
        factor_ba03 = flexural_lambda2(1.0, geom3, mat3)
        self.assertGreater(abs(factor_ba05 - factor_ba03) / factor_ba05, 0.3)
        self.assertAlmostEqual(factor_ba05 / (4.0 * _math.pi ** 2), 1.0,
                               places=9)
        lam_ip = inplane_lambda_irie(1.0, geom, mat)
        self.assertGreater(lam_ip, 0.0)
        self.assertGreater(abs(lam_ip - factor_ba05) / factor_ba05, 0.5)

    def test_ip_mixed_cf_is_4x4_ip_ss_still_raises(self):
        """IP C/F is a live row-swap; IP S/G still raise. OOP S and G live."""
        mp.dps = 26
        slv, _g, _m = make_annulus_solver(0.5, nu=0.30, motion="oop")
        Lsf = ring_L_mp(slv, 2, 0.5, bc_inner="S", bc_outer="F")
        self.assertEqual((Lsf.rows, Lsf.cols), (4, 4))
        Lfs = ring_L_mp(slv, 2, 0.5, bc_inner="F", bc_outer="SS")
        self.assertEqual((Lfs.rows, Lfs.cols), (4, 4))
        Lgg = ring_L_mp(slv, 2, 0.5, bc_inner="G", bc_outer="G")
        self.assertEqual((Lgg.rows, Lgg.cols), (4, 4))
        slv_ip, _gi, _mi = make_annulus_solver(0.5, nu=0.30, motion="ip")
        for bi, bo in (("C", "F"), ("F", "C"), ("C", "C")):
            L = ring_L_mp(slv_ip, 2, 0.5, bc_inner=bi, bc_outer=bo)
            self.assertEqual((L.rows, L.cols), (4, 4), (bi, bo))
            for i in range(4):
                for j in range(4):
                    z = complex(L[i, j])
                    self.assertTrue(
                        math.isfinite(z.real) and math.isfinite(z.imag),
                        "IP L[%s,%s]=%s for %s-%s" % (i, j, z, bi, bo))
        with self.assertRaises(UnsupportedRingBC):
            ring_L_mp(slv_ip, 2, 0.5, bc_inner="S", bc_outer="S")
        with self.assertRaises(UnsupportedRingBC):
            ring_L_mp(slv_ip, 2, 0.5, bc_inner="G", bc_outer="F")

    def test_fc_cf_lmat_is_4x4_finite(self):
        """Phase D mixed L is 4x4 and finite; F-F L is unchanged in shape."""
        mp.dps = 26
        slv, _g, _m = make_annulus_solver(0.5, nu=0.30, motion="oop")
        Lff = ring_L_mp(slv, 2, 0.5, "F", "F")
        self.assertEqual((Lff.rows, Lff.cols), (4, 4))
        for bi, bo in (("F", "C"), ("C", "F"), ("C", "C"),
                       ("S", "S"), ("S", "F"), ("F", "S"),
                       ("C", "S"), ("S", "C"),
                       ("G", "G"), ("G", "F"), ("F", "G")):
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

    def test_phase_f_table216_axisym_cc_ss(self):
        """Table 2.16 exact n=0 b/a=0.5 nu=1/3: C-C 89.30 and S-S 40.01.

        Leissa r=a is outer, r=b is inner. Rel < 1% vs the exact column.
        """
        mp.dps = 26
        slv, geom, mat = make_annulus_solver(0.5, nu=1.0 / 3.0, motion="oop")
        for inner, outer, printed in (("C", "C", 89.30), ("S", "S", 40.01)):
            native = native_from_flexural_lambda2(printed, geom, mat)
            res = ring_search(
                0, (max(1e-4, native - 0.03), native + 0.03),
                bc_inner=inner, bc_outer=outer, motion="oop",
                solver=slv, geom=geom, mat=mat,
                coarse_step=0.0005, polish=True)
            rel = abs(res.lambda2 - printed) / printed
            self.assertLess(
                rel, 0.01,
                "FAIL_PHASE_F Table 2.16 %s-%s printed=%s lam2=%s rel=%s "
                "log|det|=%s" % (inner, outer, printed, res.lambda2, rel,
                                 res.log_abs_det))

class TestPiezoElasticLimit(unittest.TestCase):
    """Paper 4 (piezoelectric ring) elastic-limit regression --
    PAPER4_PIEZO_ROADMAP.md Sec 8 step 2 ("elastic-limit regression...
    before trusting any piezo-coupled number"), LESSONS_LEARNED.md Sec
    18.149's own live ~5e-10-relative cross-check.

    With the piezo layer thickness/coupling set to zero (h1=0.0, the bare-
    host case), plate_solver.piezo_solver.PiezoOutOfPlaneSolver must
    reproduce this project's OWN validated F-F elastic ring solver
    (ring_disk.py's make_annulus_solver path / OutOfPlaneSolver._det_mp,
    Paper 3, literature-matched against Vogel/Leissa/Narita/Irie) at the
    exact geometry/material this Paper 4 thread targets: r_i=0.1, r_o=0.6,
    h=0.01, steel E=200 GPa/nu=0.3/rho=7800, n=0. This is a LIVE
    cross-check, not a hardcoded magic number alone: both solvers' own
    root-finds run inside this test.

    Pre-registered verdict: PASS if the two independent constructions
    (closed-form Bessel here, Frobenius series in core_solvers.py) agree
    to <1e-6 relative -- three orders of magnitude of headroom above the
    ~5e-10 relative agreement already demonstrated
    (Paper4_Piezo/probe_piezo_p4_elastic_ff_baseline_2026-09-15.py CHECK 1,
    LESSONS_LEARNED Sec 18.149 cross-check 2).
    """

    R_I, R_O, H = 0.1, 0.6, 0.01
    E_STEEL, NU_STEEL, RHO_STEEL = 200e9, 0.3, 7800.0

    def _ring_disk_bare_host_root(self):
        """Live root-find through THIS project's own validated F-F path
        (ring_disk.py's own geometry/material construction ->
        OutOfPlaneSolver._det_mp), converted to physical rad/s via
        geometry._omega_lit -- independent of piezo_solver.py entirely,
        the same construction LESSONS_LEARNED Sec 18.149 used."""
        b = (self.R_O - self.R_I) / 2.0
        r0 = (self.R_O + self.R_I) / 2.0
        geom = ps.make_geometry(r0 / (2 * b), 2.0, h=self.H, b=b)
        mat = ps.IsotropicMaterial(E=self.E_STEEL, nu=self.NU_STEEL,
                                   rho=self.RHO_STEEL)
        solver = ps.OutOfPlaneSolver(geom, mat, M=80, n_quad=30)
        mp.dps = 30
        from mpmath import mpc

        def det_at(Om):
            return complex(solver._det_mp(mpc(0), mpf(Om)))

        lo, hi = 0.55, 0.65
        flo = det_at(lo).real
        for _ in range(50):
            mid = (lo + hi) / 2
            fm = det_at(mid).real
            if (fm > 0) == (flo > 0):
                lo, flo = mid, fm
            else:
                hi = mid
        om_native = (lo + hi) / 2
        k_bar, w_bar = _kbar_wbar(geom, mat)
        f_hz, _om_lit = _omega_lit(om_native, geom, mat, w_bar, k_bar, part=1)
        return om_native, 2.0 * math.pi * f_hz

    def test_bare_host_matches_ring_disk_live(self):
        om_native, omega_ring_disk = self._ring_disk_bare_host_root()
        # Sanity: matches the recorded Sec 18.149 anchor (do not retune).
        self.assertAlmostEqual(om_native, 0.6007252606, places=8)

        solver = PiezoOutOfPlaneSolver(
            r_i=self.R_I, r_o=self.R_O, h=self.H,
            E=self.E_STEEL, nu=self.NU_STEEL, rho=self.RHO_STEEL,
            h1=0.0, dps=60)
        omega_piezo = solver.elastic_bisect(700.0, 800.0, 0)

        rel = abs(omega_piezo - omega_ring_disk) / omega_ring_disk
        self.assertLess(
            rel, 1e-6,
            "FAIL_PIEZO_ELASTIC_LIMIT: piezo_solver bare-host root=%s vs "
            "live ring_disk.py root=%s, rel=%s"
            % (omega_piezo, omega_ring_disk, rel))

    def test_elastic_bilayer_matches_probe_script_anchors(self):
        """CHECK 2 of probe_piezo_p4_elastic_ff_baseline_2026-09-15.py --
        ordinary two-material elastic stiffening (e31_bar=0 by
        construction), all three Duan2005 Table 4 thickness ratios. Not
        the zero-coupling limit itself, but the other elastic building
        block this module ports; anchors captured 2026-09-15/09-16, do
        not retune.

        Constants passed here are the RAW Duan2005 Table 1 values (Sec
        18.172 convention, superseding Sec 18.171's original pre-reduced
        c11_bar/c12_bar convention) -- PiezoOutOfPlaneSolver now forms
        c11_bar/c12_bar internally via _c11_bar()/_c12_bar()."""
        C11E, C12E, C13E, C33E = 132e9, 71e9, 73e9, 115e9
        rho_pzt = 7500.0
        targets = {
            2.0 / 12: 747.9863950952503,
            2.0 / 8: 764.1133267197268,
            2.0 / 5: 800.302692202868,
        }
        for frac, target in targets.items():
            h1 = frac * self.H
            solver = PiezoOutOfPlaneSolver(
                r_i=self.R_I, r_o=self.R_O, h=self.H,
                E=self.E_STEEL, nu=self.NU_STEEL, rho=self.RHO_STEEL,
                h1=h1, C11E=C11E, C12E=C12E, C13E=C13E, C33E=C33E,
                rho_pzt=rho_pzt, dps=60)
            root = solver.elastic_bisect(700.0, 850.0, 0)
            rel = abs(root - target) / target
            self.assertLess(rel, 1e-6, "h1/2h=%s root=%s target=%s rel=%s"
                            % (frac / 2, root, target, rel))

    def test_worker_scalar_reconstruction_matches_in_process(self):
        """Fast, always-run worker-path check (mirrors TestWorkerPath's own
        convention): call _piezo_elastic_ff_root_worker directly (no
        ProcessPoolExecutor) and confirm the scalar-reconstructed solver
        matches the in-process one exactly."""
        solver = PiezoOutOfPlaneSolver(
            r_i=self.R_I, r_o=self.R_O, h=self.H,
            E=self.E_STEEL, nu=self.NU_STEEL, rho=self.RHO_STEEL,
            h1=0.0, dps=60)
        in_process = solver.elastic_bisect(700.0, 800.0, 0)

        args = (self.R_I, self.R_O, self.H, self.E_STEEL, self.NU_STEEL,
                self.RHO_STEEL, 0.0, None, None, None, 60, 0, 700.0, 800.0, 50)
        via_worker = _piezo_elastic_ff_root_worker(args)
        self.assertEqual(in_process, via_worker)

    @unittest.skipUnless(SLOW, "set SLOW_TESTS=1 to run the ProcessPoolExecutor "
                                "gate (always safe/fast on the cluster)")
    def test_worker_pool_end_to_end(self):
        """Full multiprocessing round-trip: pickle the worker function and
        its plain-scalar args, run in a real subprocess, unpickle the
        result -- the actual path production jobs would use, not just the
        fast direct-call check above."""
        import concurrent.futures as _cf
        args = (self.R_I, self.R_O, self.H, self.E_STEEL, self.NU_STEEL,
                self.RHO_STEEL, 0.0, None, None, None, 60, 0, 700.0, 800.0, 50)
        with _cf.ProcessPoolExecutor(max_workers=1) as ex:
            root = ex.submit(_piezo_elastic_ff_root_worker, args).result()
class TestPiezoCoupledModel(unittest.TestCase):
    """Paper 4 (piezoelectric ring) full 3-branch piezo-coupled F-F
    determinant -- LESSONS_LEARNED.md Sec 18.172, direct follow-up to
    TestPiezoElasticLimit above (the elastic-limit regression that must
    pass, and does, before trusting any piezo-coupled number per
    PAPER4_PIEZO_ROADMAP.md Sec 8 step 2/step 3).

    Anchors below were captured this session by a live, fresh re-run of
    Paper4_Piezo/probe_piezo_p4_ff_forward_model_2026-09-15.py's own
    bisect() at full float precision (not copied from that script's
    module-docstring RESULTS TABLE, which rounds to 4 decimal places) --
    do not retune. coupled_det()/coupled_bisect() reproduce every one of
    these to <1e-12 relative in the sandbox (machine precision at
    dps=100), three-plus orders of magnitude inside the 1e-6 gate used
    here and throughout this test module.
    """

    R_I, R_O, H = 0.1, 0.6, 0.01
    E_STEEL, NU_STEEL, RHO_STEEL = 200e9, 0.3, 7800.0
    # Duan2005 Table 1, PZT4 piezo layer.
    C11E, C12E, C13E, C33E = 132e9, 71e9, 73e9, 115e9
    E31, E33 = 4.1, 14.1
    X11, X33 = 7.124e-9, 5.841e-9
    RHO_PZT = 7500.0

    def _make_solver(self, h1, dps=100):
        return PiezoOutOfPlaneSolver(
            r_i=self.R_I, r_o=self.R_O, h=self.H,
            E=self.E_STEEL, nu=self.NU_STEEL, rho=self.RHO_STEEL,
            h1=h1, C11E=self.C11E, C12E=self.C12E,
            C13E=self.C13E, C33E=self.C33E,
            e31=self.E31, e33=self.E33, X11=self.X11, X33=self.X33,
            rho_pzt=self.RHO_PZT, dps=dps)

    def test_coupled_matches_forward_model_anchors_p0(self):
        """All three Duan2005 Table 4 thickness ratios, p=0 (axisymmetric)
        -- the coupled model's headline RESULTS TABLE row."""
        targets = {
            # 2026-09-23 Sec 18.231 (consistent projection): old duan targets 747.9929429355723,
            # 764.1335816242521, 800.3737601605251.
            2.0 / 12: (747.9, 748.1, 747.9917804877718),
            2.0 / 8: (760.0, 770.0, 764.1299858014145),
            2.0 / 5: (795.0, 810.0, 800.3611443131172),
        }
        for frac, (lo, hi, target) in targets.items():
            h1 = frac * self.H
            solver = self._make_solver(h1)
            root = solver.coupled_bisect(lo, hi, 0, iters=40)
            rel = abs(root - target) / target
            self.assertLess(
                rel, 1e-6,
                "FAIL_PIEZO_COUPLED p0 h1/2h=%s root=%s target=%s rel=%s"
                % (frac / 2, root, target, rel))

    def test_coupled_matches_forward_model_anchor_p1(self):
        """p=1 (n=1 circumferential) robustness check at h1/2h=1/12 --
        confirms the 6x6 determinant behaves sensibly across
        circumferential order, not just at the one axisymmetric point."""
        h1 = 2.0 * self.H / 12
        solver = self._make_solver(h1)
        # 2026-09-22 (LESSONS Sec 18.226): target and bracket changed after the
        # free-edge row was corrected from Q_r to the Kirchhoff effective shear
        # V_r. The old value 1573.492573301136 came from the Q_r-only row and
        # was not a root of the physical problem. New value: fixed code at
        # dps=100; its elastic limit is cross-checked against ring_disk in
        # TestPiezoFreeEdgeEffectiveShear.
        root = solver.coupled_bisect(1740.0, 1770.0, 1, iters=40)
        target = 1756.0121287010907  # 2026-09-23 Sec 18.231 (consistent projection): duan 1756.0147046511106
        rel = abs(root - target) / target
        self.assertLess(
            rel, 1e-6,
            "FAIL_PIEZO_COUPLED p1 root=%s target=%s rel=%s"
            % (root, target, rel))

    def test_coupled_det_refuses_h1_zero(self):
        """coupled_det's own documented invariant (Sec 18.172): h1=0 is a
        genuine chi-cubic singularity, not a valid input -- must raise
        ValueError rather than silently returning a degenerate value.
        elastic_det(h1=0.0) is the documented alternative for that
        limit (exercised by TestPiezoElasticLimit above)."""
        solver = self._make_solver(2.0 * self.H / 12)
        with self.assertRaises(ValueError):
            solver.coupled_det(750.0, 0, h1=0.0)

    def test_worker_scalar_reconstruction_matches_in_process(self):
        """Fast, always-run worker-path check (mirrors
        TestPiezoElasticLimit's own convention): call
        _piezo_coupled_ff_root_worker directly (no ProcessPoolExecutor)
        and confirm the scalar-reconstructed solver matches the
        in-process one exactly."""
        h1 = 2.0 * self.H / 12
        solver = self._make_solver(h1)
        in_process = solver.coupled_bisect(747.9, 748.1, 0, iters=40)

        args = (self.R_I, self.R_O, self.H, self.E_STEEL, self.NU_STEEL,
                self.RHO_STEEL, h1,
                self.C11E, self.C12E, self.C13E, self.C33E,
                self.E31, self.E33, self.X11, self.X33, self.RHO_PZT,
                100, 0, 747.9, 748.1, 40)
        via_worker = _piezo_coupled_ff_root_worker(args)
        self.assertEqual(in_process, via_worker)

    @unittest.skipUnless(SLOW, "set SLOW_TESTS=1 to run the ProcessPoolExecutor "
                                "gate (always safe/fast on the cluster)")
    def test_worker_pool_end_to_end(self):
        """Full multiprocessing round-trip: pickle the worker function and
        its plain-scalar args, run in a real subprocess, unpickle the
        result -- the actual path production jobs would use, not just the
        fast direct-call check above."""
        import concurrent.futures as _cf
        h1 = 2.0 * self.H / 12
        args = (self.R_I, self.R_O, self.H, self.E_STEEL, self.NU_STEEL,
                self.RHO_STEEL, h1,
                self.C11E, self.C12E, self.C13E, self.C33E,
                self.E31, self.E33, self.X11, self.X33, self.RHO_PZT,
                100, 0, 747.9, 748.1, 40)
        with _cf.ProcessPoolExecutor(max_workers=1) as ex:
            root = ex.submit(_piezo_coupled_ff_root_worker, args).result()
        target = 747.9917804877718  # 2026-09-23 Sec 18.231 (consistent projection): duan 747.9929429355723
        rel = abs(root - target) / target
        self.assertLess(rel, 1e-6, "pooled worker root=%s target=%s rel=%s"
                        % (root, target, rel))


class TestPiezoElasticCC(unittest.TestCase):
    """Paper 4 (piezoelectric ring) C-C elastic building blocks --
    LESSONS_LEARNED.md Sec 18.173, C-C counterpart of TestPiezoElasticLimit.
    C-C is the boundary condition with an actual literature table
    (Duan2005 Table 4) and its own cluster-confirmed gate (Sec 18.148,
    job 2489140) -- these anchors are this project's own live re-
    derivation of that same table via
    Paper4_Piezo/probe_piezo_p4_elastic_cc_baseline_2026-09-15.py, not a
    hardcoded copy of the docstring's rounded values; do not retune.
    """

    R_I, R_O, H = 0.1, 0.6, 0.01
    E_STEEL, NU_STEEL, RHO_STEEL = 200e9, 0.3, 7800.0
    C11E, C12E, C13E, C33E = 132e9, 71e9, 73e9, 115e9
    RHO_PZT = 7500.0

    def test_bare_host_matches_table4_h1_zero_column(self):
        """Table 4's own h1=0 column, p=0, n_radial=0,1,2 -- Duan2005's
        printed 2718/7520/14783 rad/s, matched here (as the elastic
        baseline script itself documents) to ~0.2%, not exactly -- this
        test's own anchors are the SCRIPT's live roots, not Duan2005's
        printed numbers, so the assertion gate is tight (<1e-6) even
        though the physical/literature agreement is only ~0.2%."""
        targets = {
            0: (2712.3830136590, 2700.0, 2730.0),
            1: (7503.0602047435, 7480.0, 7560.0),
            2: (14749.3542504426, 14700.0, 14860.0),
        }
        solver = PiezoOutOfPlaneSolver(
            r_i=self.R_I, r_o=self.R_O, h=self.H,
            E=self.E_STEEL, nu=self.NU_STEEL, rho=self.RHO_STEEL,
            h1=0.0, dps=30)
        for n_radial, (target, lo, hi) in targets.items():
            root = solver.elastic_cc_bisect(lo, hi, 0)
            rel = abs(root - target) / target
            self.assertLess(rel, 1e-6, "n_radial=%s root=%s target=%s rel=%s"
                            % (n_radial, root, target, rel))

    def test_elastic_bilayer_matches_probe_script_anchors(self):
        """CHECK 2 of probe_piezo_p4_elastic_cc_baseline_2026-09-15.py --
        ordinary two-material elastic stiffening (e31_bar=0 by
        construction), all three Duan2005 Table 4 thickness ratios."""
        targets = {
            2.0 / 12: 2791.7640219016,
            2.0 / 8: 2852.1259090756,
            2.0 / 5: 2987.4935180437,
        }
        for frac, target in targets.items():
            h1 = frac * self.H
            solver = PiezoOutOfPlaneSolver(
                r_i=self.R_I, r_o=self.R_O, h=self.H,
                E=self.E_STEEL, nu=self.NU_STEEL, rho=self.RHO_STEEL,
                h1=h1, C11E=self.C11E, C12E=self.C12E,
                C13E=self.C13E, C33E=self.C33E,
                rho_pzt=self.RHO_PZT, dps=30)
            root = solver.elastic_cc_bisect(2700.0, 3150.0, 0)
            rel = abs(root - target) / target
            self.assertLess(rel, 1e-6, "h1/2h=%s root=%s target=%s rel=%s"
                            % (frac / 2, root, target, rel))

    def test_worker_scalar_reconstruction_matches_in_process(self):
        """Fast, always-run worker-path check, mirroring
        TestPiezoElasticLimit's own convention."""
        solver = PiezoOutOfPlaneSolver(
            r_i=self.R_I, r_o=self.R_O, h=self.H,
            E=self.E_STEEL, nu=self.NU_STEEL, rho=self.RHO_STEEL,
            h1=0.0, dps=30)
        in_process = solver.elastic_cc_bisect(2700.0, 2730.0, 0)

        args = (self.R_I, self.R_O, self.H, self.E_STEEL, self.NU_STEEL,
                self.RHO_STEEL, 0.0, None, None, None, 30, 0, 2700.0, 2730.0, 50)
        via_worker = _piezo_elastic_cc_root_worker(args)
        self.assertEqual(in_process, via_worker)

    @unittest.skipUnless(SLOW, "set SLOW_TESTS=1 to run the ProcessPoolExecutor "
                                "gate (always safe/fast on the cluster)")
    def test_worker_pool_end_to_end(self):
        """Full multiprocessing round-trip, mirroring TestPiezoElasticLimit's
        own convention."""
        import concurrent.futures as _cf
        args = (self.R_I, self.R_O, self.H, self.E_STEEL, self.NU_STEEL,
                self.RHO_STEEL, 0.0, None, None, None, 30, 0, 2700.0, 2730.0, 50)
        with _cf.ProcessPoolExecutor(max_workers=1) as ex:
            root = ex.submit(_piezo_elastic_cc_root_worker, args).result()
        target = 2712.3830136590
        rel = abs(root - target) / target
        self.assertLess(rel, 1e-6, "pooled worker root=%s target=%s rel=%s"
                        % (root, target, rel))


class TestPiezoCoupledCC(unittest.TestCase):
    """Paper 4 (piezoelectric ring) full 3-branch piezo-coupled C-C
    determinant -- LESSONS_LEARNED.md Sec 18.173, C-C counterpart of
    TestPiezoCoupledModel. C-C is the boundary condition already
    CLUSTER-CONFIRMED against all 27 non-trivial Duan2005 Table 4 points
    as a standalone script (Sec 18.148, job 2489140) -- this class
    verifies the package port reproduces that same standalone script's
    own live roots (captured fresh this session at dps=100, not copied
    from the module docstring's dps=200 rounded VALIDATION TABLE) to
    machine precision.
    """

    R_I, R_O, H = 0.1, 0.6, 0.01
    E_STEEL, NU_STEEL, RHO_STEEL = 200e9, 0.3, 7800.0
    C11E, C12E, C13E, C33E = 132e9, 71e9, 73e9, 115e9
    E31, E33 = 4.1, 14.1
    X11, X33 = 7.124e-9, 5.841e-9
    RHO_PZT = 7500.0

    def _make_solver(self, h1, dps=100):
        return PiezoOutOfPlaneSolver(
            r_i=self.R_I, r_o=self.R_O, h=self.H,
            E=self.E_STEEL, nu=self.NU_STEEL, rho=self.RHO_STEEL,
            h1=h1, C11E=self.C11E, C12E=self.C12E,
            C13E=self.C13E, C33E=self.C33E,
            e31=self.E31, e33=self.E33, X11=self.X11, X33=self.X33,
            rho_pzt=self.RHO_PZT, dps=dps)

    def test_coupled_matches_forward_model_anchors(self):
        """Spans three radial modes (n_radial=0,1,2), two circumferential
        orders (p=0,1), and three piezo-thickness ratios -- the same
        six-point VALIDATION TABLE the standalone C-C forward-model probe
        script itself uses, per its own docstring."""
        targets = [
            # 2026-09-23 Sec 18.231 (consistent projection): old duan targets 2791.7866198564,
            # 7722.7083975136, 15181.1339298614, 2852.1957347552,
            # 2987.7380652758, 2928.2222884367 (Duan Table 4 still 27/27
            # within 0.15% under the consistent default, worst -0.049%).
            ('p0_n0_1_12', 0, 2.0 / 12, 2780.0, 2800.0, 2791.7826080022496),
            ('p0_n1_1_12', 0, 2.0 / 12, 7700.0, 7750.0, 7722.697300567961),
            ('p0_n2_1_12', 0, 2.0 / 12, 15150.0, 15220.0, 15181.112117626035),
            ('p0_n0_1_8', 0, 2.0 / 8, 2830.0, 2880.0, 2852.1833386189246),
            ('p0_n0_1_5', 0, 2.0 / 5, 2960.0, 3020.0, 2987.6946524184314),
            ('p1_n0_1_12', 1, 2.0 / 12, 2900.0, 2960.0, 2928.2180805536336),
        ]
        for label, p, frac, lo, hi, target in targets:
            h1 = frac * self.H
            solver = self._make_solver(h1)
            root = solver.cc_coupled_bisect(lo, hi, p, iters=35)
            rel = abs(root - target) / target
            self.assertLess(
                rel, 1e-6,
                "FAIL_PIEZO_COUPLED_CC %s root=%s target=%s rel=%s"
                % (label, root, target, rel))

    def test_cc_coupled_det_refuses_h1_zero(self):
        """cc_coupled_det's own documented invariant (Sec 18.173, mirrors
        coupled_det's Sec 18.172 one): h1=0 is a genuine chi-cubic
        singularity, must raise ValueError."""
        solver = self._make_solver(2.0 * self.H / 12)
        with self.assertRaises(ValueError):
            solver.cc_coupled_det(2790.0, 0, h1=0.0)

    def test_worker_scalar_reconstruction_matches_in_process(self):
        """Fast, always-run worker-path check, mirroring
        TestPiezoCoupledModel's own convention."""
        h1 = 2.0 * self.H / 12
        solver = self._make_solver(h1)
        in_process = solver.cc_coupled_bisect(2780.0, 2800.0, 0, iters=35)

        args = (self.R_I, self.R_O, self.H, self.E_STEEL, self.NU_STEEL,
                self.RHO_STEEL, h1,
                self.C11E, self.C12E, self.C13E, self.C33E,
                self.E31, self.E33, self.X11, self.X33, self.RHO_PZT,
                100, 0, 2780.0, 2800.0, 35)
        via_worker = _piezo_coupled_cc_root_worker(args)
        self.assertEqual(in_process, via_worker)

    @unittest.skipUnless(SLOW, "set SLOW_TESTS=1 to run the ProcessPoolExecutor "
                                "gate (always safe/fast on the cluster)")
    def test_worker_pool_end_to_end(self):
        """Full multiprocessing round-trip, mirroring TestPiezoCoupledModel's
        own convention."""
        import concurrent.futures as _cf
        h1 = 2.0 * self.H / 12
        args = (self.R_I, self.R_O, self.H, self.E_STEEL, self.NU_STEEL,
                self.RHO_STEEL, h1,
                self.C11E, self.C12E, self.C13E, self.C33E,
                self.E31, self.E33, self.X11, self.X33, self.RHO_PZT,
                100, 0, 2780.0, 2800.0, 35)
        with _cf.ProcessPoolExecutor(max_workers=1) as ex:
            root = ex.submit(_piezo_coupled_cc_root_worker, args).result()
        target = 2791.7826080022496  # 2026-09-23 Sec 18.231 (consistent projection): duan 2791.7866198564
        rel = abs(root - target) / target
        self.assertLess(rel, 1e-6, "pooled worker root=%s target=%s rel=%s"
                        % (root, target, rel))


class TestPiezoOpenCircuitFF(unittest.TestCase):
    """Paper 4 F-F open-circuit (fully electroded, linear potential,
    Q=0). LESSONS_LEARNED.md Sec 18.175. Elastic-limit and stiffening
    gates must pass before any inverse-e31 number is trusted.

    Anchors captured this session from a live sandbox run of
    oc_ff_bisect at dps=40 (sign check) then confirmed at dps=60:
    elastic bilayer 747.9863950952586 rad/s, OC n=0 750.2394210586135
    rad/s (+0.301%). Do not retune.
    """

    R_I, R_O, H = 0.1, 0.6, 0.01
    E_STEEL, NU_STEEL, RHO_STEEL = 200e9, 0.3, 7800.0
    C11E, C12E, C13E, C33E = 132e9, 71e9, 73e9, 115e9
    E31, E33 = 4.1, 14.1
    X11, X33 = 7.124e-9, 5.841e-9
    RHO_PZT = 7500.0
    ELASTIC_N0 = 747.9863950952586
    # 2026-09-23 Sec 18.231 (consistent projection): OC charge arm (h+h1) -> (h+h1/2);
    # duan value 750.2394210586135 (+0.301%), now +0.280% = FE 0.280%.
    OC_N0 = 750.0794987774267

    def _make(self, dps=60):
        h1 = 2.0 * self.H / 12
        return PiezoOutOfPlaneSolver(
            r_i=self.R_I, r_o=self.R_O, h=self.H,
            E=self.E_STEEL, nu=self.NU_STEEL, rho=self.RHO_STEEL,
            h1=h1, C11E=self.C11E, C12E=self.C12E, C13E=self.C13E,
            C33E=self.C33E, e31=self.E31, e33=self.E33,
            X11=self.X11, X33=self.X33, rho_pzt=self.RHO_PZT, dps=dps)

    def test_h1_zero_matches_elastic(self):
        """h1=0 OC is identically elastic_det (no layer to charge)."""
        s = PiezoOutOfPlaneSolver(
            r_i=self.R_I, r_o=self.R_O, h=self.H,
            E=self.E_STEEL, nu=self.NU_STEEL, rho=self.RHO_STEEL,
            h1=0.0, dps=40)
        omega = 726.81185755
        self.assertEqual(s.oc_ff_det(omega, 0), s.elastic_det(omega, 0))

    def test_n_nonzero_matches_elastic(self):
        """n!=0 OC is identically elastic_det (net electrode charge
        cancels). Checked at an n=1 frequency, not just a root."""
        s = self._make(dps=40)
        omega = 1575.0
        self.assertEqual(s.oc_ff_det(omega, 1), s.elastic_det(omega, 1))

    def test_n0_stiffens_vs_elastic_bilayer(self):
        """n=0 OC root is ABOVE the elastic bilayer (stiffening) and
        matches the sandbox-captured anchor to 1e-6 relative."""
        s = self._make(dps=60)
        elastic = s.elastic_bisect(720.0, 780.0, 0)
        oc = s.oc_ff_bisect(720.0, 800.0, 0)
        self.assertGreater(oc, elastic)
        rel_el = abs(elastic - self.ELASTIC_N0) / self.ELASTIC_N0
        rel_oc = abs(oc - self.OC_N0) / self.OC_N0
        self.assertLess(rel_el, 1e-6, "elastic=%s anchor=%s rel=%s"
                        % (elastic, self.ELASTIC_N0, rel_el))
        self.assertLess(rel_oc, 1e-6, "oc=%s anchor=%s rel=%s"
                        % (oc, self.OC_N0, rel_oc))
        drel = (oc - elastic) / elastic
        self.assertGreater(drel, 0.002)
        self.assertLess(drel, 0.005)

    def test_worker_scalar_reconstruction_matches_in_process(self):
        s = self._make(dps=60)
        in_process = s.oc_ff_bisect(720.0, 800.0, 0)
        h1 = 2.0 * self.H / 12
        args = (self.R_I, self.R_O, self.H, self.E_STEEL, self.NU_STEEL,
                self.RHO_STEEL, h1,
                self.C11E, self.C12E, self.C13E, self.C33E,
                self.E31, self.E33, self.X11, self.X33, self.RHO_PZT,
                60, 0, 720.0, 800.0, 50)
        via_worker = _piezo_oc_ff_root_worker(args)
        self.assertEqual(in_process, via_worker)

    @unittest.skipUnless(SLOW, "set SLOW_TESTS=1 to run the ProcessPoolExecutor "
                                "gate (always safe/fast on the cluster)")
    def test_worker_pool_end_to_end(self):
        import concurrent.futures as _cf
        h1 = 2.0 * self.H / 12
        args = (self.R_I, self.R_O, self.H, self.E_STEEL, self.NU_STEEL,
                self.RHO_STEEL, h1,
                self.C11E, self.C12E, self.C13E, self.C33E,
                self.E31, self.E33, self.X11, self.X33, self.RHO_PZT,
                60, 0, 720.0, 800.0, 50)
        with _cf.ProcessPoolExecutor(max_workers=1) as ex:
            root = ex.submit(_piezo_oc_ff_root_worker, args).result()
        rel = abs(root - self.OC_N0) / self.OC_N0
        self.assertLess(rel, 1e-6, "pooled worker root=%s target=%s rel=%s"
                        % (root, self.OC_N0, rel))


class TestPiezoDrivenAdmittanceFF(unittest.TestCase):
    """Paper 4 F-F driven Y(omega) = j omega Q/V, linear-potential
    model. LESSONS_LEARNED.md Sec 18.179. Poles at elastic F-F (SC),
    zeros at OC. Sign of driven_ff_qv was sandbox-confirmed against
    the Sec 18.175 OC root (q/v ~ 4e-17 at omega_OC).
    """

    R_I, R_O, H = 0.1, 0.6, 0.01
    E_STEEL, NU_STEEL, RHO_STEEL = 200e9, 0.3, 7800.0
    C11E, C12E, C13E, C33E = 132e9, 71e9, 73e9, 115e9
    E31, E33 = 4.1, 14.1
    X11, X33 = 7.124e-9, 5.841e-9
    RHO_PZT = 7500.0
    ELASTIC_N0 = 747.9863950952586
    # 2026-09-23 Sec 18.231 (consistent projection): OC charge arm (h+h1) -> (h+h1/2);
    # duan value 750.2394210586135 (+0.301%), now +0.280% = FE 0.280%.
    OC_N0 = 750.0794987774267

    def _make(self, dps=40):
        h1 = 2.0 * self.H / 12
        return PiezoOutOfPlaneSolver(
            r_i=self.R_I, r_o=self.R_O, h=self.H,
            E=self.E_STEEL, nu=self.NU_STEEL, rho=self.RHO_STEEL,
            h1=h1, C11E=self.C11E, C12E=self.C12E, C13E=self.C13E,
            C33E=self.C33E, e31=self.E31, e33=self.E33,
            X11=self.X11, X33=self.X33, rho_pzt=self.RHO_PZT, dps=dps)

    def test_n_nonzero_is_clamped_C0(self):
        s = self._make()
        c0 = s.clamped_C0()
        r = s.driven_ff_qv(750.0, 1)
        self.assertEqual(r["q_over_v"], c0)
        self.assertEqual(r["s"], 0.0)

    def test_zero_at_oc_root(self):
        s = self._make()
        c0 = s.clamped_C0()
        r = s.driven_ff_qv(self.OC_N0, 0)
        self.assertLess(abs(r["q_over_v"]) / c0, 1e-8)

    def test_pole_sign_flip_across_sc_root(self):
        s = self._make()
        c0 = s.clamped_C0()
        lo = s.driven_ff_qv(self.ELASTIC_N0 - 0.5, 0)["q_over_v"]
        hi = s.driven_ff_qv(self.ELASTIC_N0 + 0.5, 0)["q_over_v"]
        self.assertGreater(abs(lo) / c0, 2.0)
        self.assertGreater(abs(hi) / c0, 2.0)
        self.assertLess(lo * hi, 0.0)

    def test_above_antiresonance_returns_toward_C0(self):
        s = self._make()
        c0 = s.clamped_C0()
        r = s.driven_ff_qv(800.0, 0)
        ratio = r["q_over_v"] / c0
        self.assertGreater(ratio, 0.5)
        self.assertLess(ratio, 1.0)

class TestPiezoOCSineFF(unittest.TestCase):
    """Paper 4 F-F open-circuit, unified linear+interior-sine ansatz
    (LESSONS_LEARNED.md Sec 18.182, lever 12). Same-ansatz counterpart
    of TestPiezoOpenCircuitFF above -- oc_ff_det/oc_ff_bisect (the
    linear-only, leading-order model tested there) are UNCHANGED; this
    class exercises the new oc_ff_sine_det/oc_ff_sine_bisect added
    alongside them.

    Anchors below were captured this session at dps=100 (production
    precision, matching TestPiezoCoupledModel's own convention), all
    three Duan2005 Table 4 thickness ratios. Do not retune -- if a
    future change to _branches/_radial_quad/coupled_det moves these,
    that is a regression to investigate, not a reason to update the
    anchor.
    """

    R_I, R_O, H = 0.1, 0.6, 0.01
    E_STEEL, NU_STEEL, RHO_STEEL = 200e9, 0.3, 7800.0
    C11E, C12E, C13E, C33E = 132e9, 71e9, 73e9, 115e9
    E31, E33 = 4.1, 14.1
    X11, X33 = 7.124e-9, 5.841e-9
    RHO_PZT = 7500.0

    # (h1_fraction_of_2h, elastic_bilayer, sc_coupled_det, oc_sine)
    RATIO_ANCHORS = [
        # 2026-09-23 Sec 18.231 (consistent projection): duan sc/oc_sine were 747.9929429355705 /
        # 750.1950280038534, 764.1335816242539 / 767.4613367586278,
        # 800.3737601605269 / 805.7204479759378. Consistent OC-sine =
        # SC stiffening + linear-OC stiffening (additive), no S_phi term.
        (2.0 / 12, 747.9863950952503, 747.99178048777, 750.0848429101234),
        (2.0 / 8, 764.1133267197267, 764.1299858014199, 767.2255947436147),
        (2.0 / 5, 800.302692202868, 800.3611443131117, 805.1698609187677),
    ]

    def _make(self, h1, dps=100):
        return PiezoOutOfPlaneSolver(
            r_i=self.R_I, r_o=self.R_O, h=self.H,
            E=self.E_STEEL, nu=self.NU_STEEL, rho=self.RHO_STEEL,
            h1=h1, C11E=self.C11E, C12E=self.C12E, C13E=self.C13E,
            C33E=self.C33E, e31=self.E31, e33=self.E33,
            X11=self.X11, X33=self.X33, rho_pzt=self.RHO_PZT, dps=dps)

    def test_h1_zero_matches_elastic(self):
        """h1=0 same-ansatz OC is identically elastic_det (no layer to
        charge) -- bit-for-bit, mirroring oc_ff_det's own invariant."""
        s = PiezoOutOfPlaneSolver(
            r_i=self.R_I, r_o=self.R_O, h=self.H,
            E=self.E_STEEL, nu=self.NU_STEEL, rho=self.RHO_STEEL,
            h1=0.0, dps=40)
        omega = 726.81185755
        self.assertEqual(s.oc_ff_sine_det(omega, 0), s.elastic_det(omega, 0))

    def test_n_nonzero_matches_coupled(self):
        """n!=0 same-ansatz OC is identically coupled_det (the uniform
        bus voltage cannot couple to a non-axisymmetric mode) -- bit-
        for-bit, checked at n=1 and n=2, not just a root."""
        h1 = 2.0 * self.H / 12
        s = self._make(h1, dps=40)
        for n in (1, 2):
            omega = 1600.0
            self.assertEqual(s.oc_ff_sine_det(omega, n),
                             s.coupled_det(omega, n))

    def test_n0_stiffens_and_matches_anchors(self):
        """n=0 same-ansatz OC root is ABOVE both the elastic bilayer
        AND the SC (coupled_det) root at all three thickness ratios --
        the piezoelectric c^D >= c^E stiffening theorem must hold for
        any correctly-posed open-circuit model (Sec 18.182) -- and
        matches the dps=100 anchors above to 1e-6 relative."""
        for frac, elastic_anchor, sc_anchor, oc_anchor in self.RATIO_ANCHORS:
            h1 = frac * self.H
            s = self._make(h1)
            elastic = s.elastic_bisect(700.0, 850.0, 0)
            sc = s.coupled_bisect(sc_anchor - 2.0, sc_anchor + 2.0, 0, iters=40)
            oc = s.oc_ff_sine_bisect(oc_anchor - 30.0, oc_anchor + 30.0, 0, iters=45)
            self.assertGreater(oc, elastic,
                "h1/2h=%s: OC=%s not above elastic=%s" % (frac / 2, oc, elastic))
            self.assertGreater(oc, sc,
                "h1/2h=%s: OC=%s not above SC=%s" % (frac / 2, oc, sc))
            for val, anchor, label in (
                    (elastic, elastic_anchor, "elastic"),
                    (sc, sc_anchor, "sc"), (oc, oc_anchor, "oc")):
                rel = abs(val - anchor) / anchor
                self.assertLess(rel, 1e-6,
                    "h1/2h=%s %s=%s anchor=%s rel=%s"
                    % (frac / 2, label, val, anchor, rel))

    def test_self_consistent_with_driven_admittance(self):
        """oc_ff_sine_det's own root, fed into driven_ff_qv_sine, gives
        q_over_v ~0 -- the two methods' independent formula copies must
        agree algebraically if wired consistently (Sec 18.182 point 4;
        this is the check that actually caught the sign bugs fixed in
        that section)."""
        for frac, _, _, oc_anchor in self.RATIO_ANCHORS:
            h1 = frac * self.H
            s = self._make(h1)
            oc = s.oc_ff_sine_bisect(oc_anchor - 30.0, oc_anchor + 30.0, 0, iters=45)
            r = s.driven_ff_qv_sine(oc, 0)
            c0 = s.clamped_C0(h1=h1)
            self.assertLess(abs(r["q_over_v"]) / c0, 1e-8,
                "h1/2h=%s q_over_v/c0=%s" % (frac / 2, r["q_over_v"] / c0))

    def test_worker_scalar_reconstruction_matches_in_process(self):
        h1 = 2.0 * self.H / 12
        s = self._make(h1)
        in_process = s.oc_ff_sine_bisect(720.0, 800.0, 0, iters=45)
        args = (self.R_I, self.R_O, self.H, self.E_STEEL, self.NU_STEEL,
                self.RHO_STEEL, h1,
                self.C11E, self.C12E, self.C13E, self.C33E,
                self.E31, self.E33, self.X11, self.X33, self.RHO_PZT,
                100, 0, 720.0, 800.0, 45)
        via_worker = _piezo_oc_sine_ff_root_worker(args)
        self.assertEqual(in_process, via_worker)

    @unittest.skipUnless(SLOW, "set SLOW_TESTS=1 to run the ProcessPoolExecutor "
                                "gate (always safe/fast on the cluster)")
    def test_worker_pool_end_to_end(self):
        import concurrent.futures as _cf
        h1 = 2.0 * self.H / 12
        args = (self.R_I, self.R_O, self.H, self.E_STEEL, self.NU_STEEL,
                self.RHO_STEEL, h1,
                self.C11E, self.C12E, self.C13E, self.C33E,
                self.E31, self.E33, self.X11, self.X33, self.RHO_PZT,
                100, 0, 720.0, 800.0, 45)
        with _cf.ProcessPoolExecutor(max_workers=1) as ex:
            root = ex.submit(_piezo_oc_sine_ff_root_worker, args).result()
        target = self.RATIO_ANCHORS[0][3]
        rel = abs(root - target) / target
        self.assertLess(rel, 1e-6, "pooled worker root=%s target=%s rel=%s"
                        % (root, target, rel))


class TestPiezoDrivenAdmittanceSineFF(unittest.TestCase):
    """Paper 4 F-F driven Y(omega) = j omega Q/V on the same-ansatz
    model (LESSONS_LEARNED.md Sec 18.182). Same-ansatz counterpart of
    TestPiezoDrivenAdmittanceFF above -- driven_ff_qv (linear-only) is
    UNCHANGED; this exercises the new driven_ff_qv_sine.

    Poles at coupled_det's own roots (NOT the elastic bilayer -- this
    model's 6x6 matrix IS coupled_det's), zeros at oc_ff_sine_det's
    roots.
    """

    R_I, R_O, H = 0.1, 0.6, 0.01
    E_STEEL, NU_STEEL, RHO_STEEL = 200e9, 0.3, 7800.0
    C11E, C12E, C13E, C33E = 132e9, 71e9, 73e9, 115e9
    E31, E33 = 4.1, 14.1
    X11, X33 = 7.124e-9, 5.841e-9
    RHO_PZT = 7500.0
    SC_N0 = 747.99178048777  # 2026-09-23 Sec 18.231 (consistent projection): duan 747.9929429355705
    OC_N0 = 750.0848429101234  # duan 750.1950280038534

    def _make(self, dps=60):
        h1 = 2.0 * self.H / 12
        return PiezoOutOfPlaneSolver(
            r_i=self.R_I, r_o=self.R_O, h=self.H,
            E=self.E_STEEL, nu=self.NU_STEEL, rho=self.RHO_STEEL,
            h1=h1, C11E=self.C11E, C12E=self.C12E, C13E=self.C13E,
            C33E=self.C33E, e31=self.E31, e33=self.E33,
            X11=self.X11, X33=self.X33, rho_pzt=self.RHO_PZT, dps=dps)

    def test_n_nonzero_is_clamped_C0(self):
        s = self._make()
        c0 = s.clamped_C0()
        r = s.driven_ff_qv_sine(750.0, 1)
        self.assertEqual(r["q_over_v"], c0)
        self.assertEqual(r["s"], 0.0)

    def test_zero_at_oc_sine_root(self):
        s = self._make()
        c0 = s.clamped_C0()
        r = s.driven_ff_qv_sine(self.OC_N0, 0)
        self.assertLess(abs(r["q_over_v"]) / c0, 1e-8)

    def test_pole_sign_flip_across_sc_root(self):
        """Pole sits at coupled_det's SC root -- NOT the elastic
        bilayer, since this model's 6x6 IS coupled_det's own matrix."""
        s = self._make()
        c0 = s.clamped_C0()
        lo = s.driven_ff_qv_sine(self.SC_N0 - 0.5, 0)["q_over_v"]
        hi = s.driven_ff_qv_sine(self.SC_N0 + 0.5, 0)["q_over_v"]
        self.assertGreater(abs(lo) / c0, 2.0)
        self.assertGreater(abs(hi) / c0, 2.0)
        self.assertLess(lo * hi, 0.0)


class TestPiezoElasticCF(unittest.TestCase):
    """Paper 4 lever 10 mixed-edge elastic 4x4. C-F = inner C, outer F
    (Paper 3 ring_disk.py). h1=0 must match ring_disk live; mixed
    both-free / both-clamped recover the existing F-F / C-C dets.
    Anchors from the 2026-09-16 sandbox scan, dps=50 then confirmed
    in the SENTINEL probe; do not retune.
    """

    R_I, R_O, H = 0.1, 0.6, 0.01
    E_STEEL, NU_STEEL, RHO_STEEL = 200e9, 0.3, 7800.0
    C11E, C12E, C13E, C33E = 132e9, 71e9, 73e9, 115e9
    RHO_PZT = 7500.0
    CF_BARE = 410.07418102201177
    FC_BARE = 873.2216324549811
    CF_BILAYER = 421.89468129782034
    FC_BILAYER = 899.0613175542029

    def _bare(self, dps=50):
        return PiezoOutOfPlaneSolver(
            r_i=self.R_I, r_o=self.R_O, h=self.H,
            E=self.E_STEEL, nu=self.NU_STEEL, rho=self.RHO_STEEL,
            h1=0.0, dps=dps)

    def _bilayer(self, dps=50):
        h1 = 2.0 * self.H / 12
        return PiezoOutOfPlaneSolver(
            r_i=self.R_I, r_o=self.R_O, h=self.H,
            E=self.E_STEEL, nu=self.NU_STEEL, rho=self.RHO_STEEL,
            h1=h1, C11E=self.C11E, C12E=self.C12E, C13E=self.C13E,
            C33E=self.C33E, rho_pzt=self.RHO_PZT, dps=dps)

    def _ring_disk_omega(self, bc_inner, bc_outer, om_lo, om_hi):
        b = (self.R_O - self.R_I) / 2.0
        r0 = (self.R_O + self.R_I) / 2.0
        geom = ps.make_geometry(r0 / (2 * b), 2.0, h=self.H, b=b)
        mat = ps.IsotropicMaterial(E=self.E_STEEL, nu=self.NU_STEEL,
                                   rho=self.RHO_STEEL)
        solver = ps.OutOfPlaneSolver(geom, mat, M=80, n_quad=30)
        mp.dps = 30
        scale = 726.8118579070 / 0.6007252606
        lo = om_lo / scale
        hi = om_hi / scale

        def det_at(Om):
            return complex(ring_det_mp(solver, 0, Om, bc_inner, bc_outer)).real

        flo = det_at(lo)
        for _ in range(50):
            mid = 0.5 * (lo + hi)
            fm = det_at(mid)
            if (fm > 0) == (flo > 0):
                lo, flo = mid, fm
            else:
                hi = mid
        Om = 0.5 * (lo + hi)
        k_bar, w_bar = _kbar_wbar(geom, mat)
        f_hz, _ = _omega_lit(Om, geom, mat, w_bar, k_bar, part=1)
        return 2.0 * math.pi * f_hz

    def test_mixed_engine_recovers_ff_and_cc(self):
        s = self._bare(dps=40)
        om = 726.81185755
        self.assertEqual(
            s._elastic_mixed_det(om, 0, 0.0, "F", "F", oc=False),
            s.elastic_det(om, 0))
        self.assertEqual(
            s._elastic_mixed_det(om, 0, 0.0, "C", "C", oc=False),
            s.elastic_cc_det(om, 0))

    def test_bare_host_matches_ring_disk_live(self):
        s = self._bare(dps=50)
        cf = s.elastic_cf_bisect(400.0, 450.0, 0)
        fc = s.elastic_fc_bisect(850.0, 920.0, 0)
        rd_cf = self._ring_disk_omega("C", "F", 400.0, 450.0)
        rd_fc = self._ring_disk_omega("F", "C", 850.0, 920.0)
        self.assertLess(abs(cf - rd_cf) / rd_cf, 1e-6, "C-F %s vs %s" % (cf, rd_cf))
        self.assertLess(abs(fc - rd_fc) / rd_fc, 1e-6, "F-C %s vs %s" % (fc, rd_fc))
        self.assertLess(abs(cf - self.CF_BARE) / self.CF_BARE, 1e-6)
        self.assertLess(abs(fc - self.FC_BARE) / self.FC_BARE, 1e-6)

    def test_elastic_bilayer_anchors(self):
        s = self._bilayer(dps=50)
        cf = s.elastic_cf_bisect(400.0, 450.0, 0)
        fc = s.elastic_fc_bisect(850.0, 920.0, 0)
        self.assertLess(abs(cf - self.CF_BILAYER) / self.CF_BILAYER, 1e-6)
        self.assertLess(abs(fc - self.FC_BILAYER) / self.FC_BILAYER, 1e-6)

    def test_worker_scalar_reconstruction_matches_in_process(self):
        s = self._bare(dps=50)
        in_process = s.elastic_cf_bisect(400.0, 450.0, 0)
        args = (self.R_I, self.R_O, self.H, self.E_STEEL, self.NU_STEEL,
                self.RHO_STEEL, 0.0, None, None, None, None, None, 50,
                0, 400.0, 450.0, 50)
        via_worker = _piezo_elastic_cf_root_worker(args)
        self.assertEqual(in_process, via_worker)


class TestPiezoCoupledCF(unittest.TestCase):
    """Mixed-edge SC 6x6. Elastic limit of the coupled root is the
    mixed bilayer; h1=0 is refused. F-F/C-C identities of the mixed
    engine recover coupled_det / cc_coupled_det.
    """

    R_I, R_O, H = 0.1, 0.6, 0.01
    E_STEEL, NU_STEEL, RHO_STEEL = 200e9, 0.3, 7800.0
    C11E, C12E, C13E, C33E = 132e9, 71e9, 73e9, 115e9
    E31, E33 = 4.1, 14.1
    X11, X33 = 7.124e-9, 5.841e-9
    RHO_PZT = 7500.0
    CF_SC = 421.89823972456895  # 2026-09-23 Sec 18.231 (consistent projection): duan 421.89900782092377
    # 2026-09-22 (Sec 18.226): was 364.894915861246 under the Q_r-only free-edge
    # row. Corrected V_r row gives the value below (elastic C-F n=1 bilayer
    # 380.26952903233564).
    CF_N1_SC = 380.2721396843914  # 2026-09-23 Sec 18.231 (consistent projection): duan 380.272703201043

    def _make(self, dps=50):
        h1 = 2.0 * self.H / 12
        return PiezoOutOfPlaneSolver(
            r_i=self.R_I, r_o=self.R_O, h=self.H,
            E=self.E_STEEL, nu=self.NU_STEEL, rho=self.RHO_STEEL,
            h1=h1, C11E=self.C11E, C12E=self.C12E, C13E=self.C13E,
            C33E=self.C33E, e31=self.E31, e33=self.E33,
            X11=self.X11, X33=self.X33, rho_pzt=self.RHO_PZT, dps=dps)

    def test_mixed_engine_recovers_ff_and_cc(self):
        s = self._make(dps=40)
        om = 747.9863950952586
        self.assertEqual(
            s._coupled_mixed_det(om, 0, s.h1, "F", "F"),
            s.coupled_det(om, 0))
        self.assertEqual(
            s._coupled_mixed_det(om, 0, s.h1, "C", "C"),
            s.cc_coupled_det(om, 0))

    def test_cf_coupled_det_refuses_h1_zero(self):
        s = PiezoOutOfPlaneSolver(
            r_i=self.R_I, r_o=self.R_O, h=self.H,
            E=self.E_STEEL, nu=self.NU_STEEL, rho=self.RHO_STEEL,
            h1=0.0, dps=40)
        with self.assertRaises(ValueError):
            s.cf_coupled_det(410.0, 0, h1=0.0)
        with self.assertRaises(ValueError):
            s.fc_coupled_det(873.0, 0, h1=0.0)

    def test_n0_and_n1_near_elastic_bilayer(self):
        s = self._make(dps=50)
        el = s.elastic_cf_bisect(400.0, 450.0, 0)
        sc = s.cf_coupled_bisect(400.0, 450.0, 0)
        self.assertLess(abs(sc - self.CF_SC) / self.CF_SC, 1e-6)
        self.assertLess(abs(sc - el) / el, 5e-4)
        sc_n1 = s.cf_coupled_bisect(350.0, 400.0, 1)
        self.assertLess(abs(sc_n1 - self.CF_N1_SC) / self.CF_N1_SC, 1e-6)


class TestPiezoOpenCircuitCF(unittest.TestCase):
    """Mixed-edge linear OC. n=0 C-F stiffens vs the C-F bilayer
    (comparable to F-F +0.301%); n=0 F-C stiffens much less (only the
    inner edge is free). n!=0 and h1=0 reduce identically to elastic.
    Pre-registered in PAPER4_CF_DERIVATION.md Sec 4 before the OC
    methods were run: stiffen, not coincide.
    """

    R_I, R_O, H = 0.1, 0.6, 0.01
    E_STEEL, NU_STEEL, RHO_STEEL = 200e9, 0.3, 7800.0
    C11E, C12E, C13E, C33E = 132e9, 71e9, 73e9, 115e9
    E31, E33 = 4.1, 14.1
    X11, X33 = 7.124e-9, 5.841e-9
    RHO_PZT = 7500.0
    CF_EL = 421.89468129782034
    CF_OC = 423.0008413290095  # 2026-09-23 Sec 18.231 (consistent projection): duan 423.0853515180034
    FC_EL = 899.0613175542029
    FC_OC = 899.1096715294408  # 2026-09-23 Sec 18.231 (consistent projection): duan 899.1133893217434

    def _make(self, dps=50):
        h1 = 2.0 * self.H / 12
        return PiezoOutOfPlaneSolver(
            r_i=self.R_I, r_o=self.R_O, h=self.H,
            E=self.E_STEEL, nu=self.NU_STEEL, rho=self.RHO_STEEL,
            h1=h1, C11E=self.C11E, C12E=self.C12E, C13E=self.C13E,
            C33E=self.C33E, e31=self.E31, e33=self.E33,
            X11=self.X11, X33=self.X33, rho_pzt=self.RHO_PZT, dps=dps)

    def test_h1_zero_matches_elastic(self):
        s = PiezoOutOfPlaneSolver(
            r_i=self.R_I, r_o=self.R_O, h=self.H,
            E=self.E_STEEL, nu=self.NU_STEEL, rho=self.RHO_STEEL,
            h1=0.0, dps=40)
        self.assertEqual(s.oc_cf_det(410.07418102, 0),
                         s.elastic_cf_det(410.07418102, 0))
        self.assertEqual(s.oc_fc_det(873.22163245, 0),
                         s.elastic_fc_det(873.22163245, 0))

    def test_n_nonzero_matches_elastic(self):
        s = self._make(dps=40)
        om = 365.0
        self.assertEqual(s.oc_cf_det(om, 1), s.elastic_cf_det(om, 1))
        self.assertEqual(s.oc_fc_det(om, 1), s.elastic_fc_det(om, 1))

    def test_n0_stiffens_vs_elastic_bilayer(self):
        s = self._make(dps=50)
        cf_el = s.elastic_cf_bisect(400.0, 450.0, 0)
        cf_oc = s.oc_cf_bisect(400.0, 450.0, 0)
        self.assertGreater(cf_oc, cf_el)
        self.assertLess(abs(cf_el - self.CF_EL) / self.CF_EL, 1e-6)
        self.assertLess(abs(cf_oc - self.CF_OC) / self.CF_OC, 1e-6)
        drel_cf = (cf_oc - cf_el) / cf_el
        self.assertGreater(drel_cf, 0.002)
        self.assertLess(drel_cf, 0.005)
        fc_el = s.elastic_fc_bisect(850.0, 920.0, 0)
        fc_oc = s.oc_fc_bisect(850.0, 920.0, 0)
        self.assertGreater(fc_oc, fc_el)
        self.assertLess(abs(fc_el - self.FC_EL) / self.FC_EL, 1e-6)
        self.assertLess(abs(fc_oc - self.FC_OC) / self.FC_OC, 1e-6)
        drel_fc = (fc_oc - fc_el) / fc_el
        self.assertGreater(drel_fc, 1e-6)
        self.assertLess(drel_fc, 0.01)
        self.assertLess(drel_fc, drel_cf)

    def test_worker_scalar_reconstruction_matches_in_process(self):
        s = self._make(dps=50)
        in_process = s.oc_cf_bisect(400.0, 450.0, 0)
        h1 = 2.0 * self.H / 12
        args = (self.R_I, self.R_O, self.H, self.E_STEEL, self.NU_STEEL,
                self.RHO_STEEL, h1,
                self.C11E, self.C12E, self.C13E, self.C33E,
                self.E31, self.E33, self.X11, self.X33, self.RHO_PZT,
                50, 0, 400.0, 450.0, 50)
        via_worker = _piezo_oc_cf_root_worker(args)
        self.assertEqual(in_process, via_worker)


class TestPiezoMonoElastic(unittest.TestCase):
    """Paper 5 Option A elastic 4x4. PZT-4 Duan Table 1 / Paper 4
    BASE_KWARGS, H=0.01 (full thickness 0.02 m), r_i=0.1, r_o=0.6.
    Elastic F-F must match live ring_disk.py under the isotropic
    reduction of PAPER5_DERIVATION.md Sec 5. Anchors captured
    2026-09-19 at dps=60; do not retune.
    """

    R_I, R_O, H = 0.1, 0.6, 0.01
    C11E, C12E, C13E, C33E = 132e9, 71e9, 73e9, 115e9
    RHO = 7500.0
    EL_FF = 462.3976497355817
    EL_CC = 1726.889795840675

    def _make(self, dps=60):
        return PiezoMonolithicOutOfPlaneSolver(
            r_i=self.R_I, r_o=self.R_O, H=self.H,
            C11E=self.C11E, C12E=self.C12E, C13E=self.C13E,
            C33E=self.C33E, rho=self.RHO, dps=dps)

    def test_ff_matches_ring_disk_live(self):
        s = self._make(dps=60)
        omega_el = s.elastic_bisect(430.0, 520.0, 0)
        self.assertLess(abs(omega_el - self.EL_FF) / self.EL_FF, 1e-8)

        E = s.isotropic_E()
        nu = s.isotropic_nu()
        b = (self.R_O - self.R_I) / 2.0
        r0 = (self.R_O + self.R_I) / 2.0
        geom = ps.make_geometry(r0 / (2 * b), 2.0, h=self.H, b=b)
        mat = ps.IsotropicMaterial(E=E, nu=nu, rho=self.RHO)
        solver = ps.OutOfPlaneSolver(geom, mat, M=80, n_quad=30)
        mp.dps = 30
        from mpmath import mpc

        def det_at(Om):
            return complex(solver._det_mp(mpc(0), mpf(Om)))

        lo, hi = 0.55, 0.65
        flo = det_at(lo).real
        for _ in range(50):
            mid = (lo + hi) / 2
            fm = det_at(mid).real
            if (fm > 0) == (flo > 0):
                lo, flo = mid, fm
            else:
                hi = mid
        om_native = (lo + hi) / 2
        k_bar, w_bar = _kbar_wbar(geom, mat)
        f_hz, _ = _omega_lit(om_native, geom, mat, w_bar, k_bar, part=1)
        omega_rd = 2.0 * math.pi * f_hz
        rel = abs(omega_el - omega_rd) / omega_rd
        self.assertLess(
            rel, 1e-6,
            "FAIL_P5_ELASTIC_LIMIT: mono elastic FF=%s vs live "
            "ring_disk=%s rel=%s" % (omega_el, omega_rd, rel))

    def test_cc_anchor(self):
        s = self._make(dps=60)
        omega_cc = s.elastic_cc_bisect(1600.0, 1850.0, 0)
        self.assertLess(abs(omega_cc - self.EL_CC) / self.EL_CC, 1e-8)

    def test_worker_scalar_reconstruction_matches_in_process(self):
        s = self._make(dps=60)
        in_process = s.elastic_bisect(430.0, 520.0, 0)
        args = (self.R_I, self.R_O, self.H, self.C11E, self.C12E,
                self.C13E, self.C33E, self.RHO, 60, 0, 430.0, 520.0, 50)
        via_worker = _piezo_mono_elastic_ff_root_worker(args)
        self.assertEqual(in_process, via_worker)


class TestPiezoMonoCoupled(unittest.TestCase):
    """Paper 5 Option A short-circuit 6x6. Pre-register: SC above
    elastic by a paper-visible amount (~2% at this PZT-4 point);
    e31=0 refused; three chi-cubic branches. Anchors 2026-09-19
    dps=60. No OC frequency solver.
    """

    R_I, R_O, H = 0.1, 0.6, 0.01
    C11E, C12E, C13E, C33E = 132e9, 71e9, 73e9, 115e9
    E31, E33 = 4.1, 14.1
    X11, X33 = 7.124e-9, 5.841e-9
    RHO = 7500.0
    EL_FF = 462.3976497355817
    # 2026-09-23 Sec 18.231 (consistent projection): duan CP_FF 473.2496088764243
    # (+2.347%), CP_CC 1764.2691696505472 (+2.165%).
    CP_FF = 471.36460219209835
    EL_CC = 1726.889795840675
    CP_CC = 1757.694304873462

    def _make(self, dps=60):
        return PiezoMonolithicOutOfPlaneSolver(
            r_i=self.R_I, r_o=self.R_O, H=self.H,
            C11E=self.C11E, C12E=self.C12E, C13E=self.C13E,
            C33E=self.C33E, rho=self.RHO,
            e31=self.E31, e33=self.E33, X11=self.X11, X33=self.X33,
            dps=dps)

    def test_sc_stiffens_vs_elastic_ff(self):
        s = self._make(dps=60)
        el = s.elastic_bisect(430.0, 520.0, 0)
        cp = s.coupled_bisect(460.0, 490.0, 0)
        self.assertLess(abs(el - self.EL_FF) / self.EL_FF, 1e-8)
        self.assertLess(abs(cp - self.CP_FF) / self.CP_FF, 1e-8)
        drel = (cp - el) / el
        self.assertGreater(drel, 0.015)
        self.assertLess(drel, 0.04)

    def test_sc_stiffens_vs_elastic_cc(self):
        s = self._make(dps=60)
        el = s.elastic_cc_bisect(1600.0, 1850.0, 0)
        cp = s.cc_coupled_bisect(1740.0, 1810.0, 0)
        self.assertLess(abs(el - self.EL_CC) / self.EL_CC, 1e-8)
        self.assertLess(abs(cp - self.CP_CC) / self.CP_CC, 1e-8)
        drel = (cp - el) / el
        self.assertGreater(cp, el)
        self.assertGreater(drel, 0.015)
        self.assertLess(drel, 0.04)

    def test_three_branches(self):
        s = self._make(dps=40)
        br = s._branches(self.CP_FF)
        self.assertEqual(len(br), 3)

    def test_coupled_refuses_missing_e31(self):
        s = PiezoMonolithicOutOfPlaneSolver(
            r_i=self.R_I, r_o=self.R_O, H=self.H,
            C11E=self.C11E, C12E=self.C12E, C13E=self.C13E,
            C33E=self.C33E, rho=self.RHO, dps=40)
        with self.assertRaises(ValueError):
            s.coupled_det(470.0, 0)

    def test_worker_scalar_reconstruction_matches_in_process(self):
        s = self._make(dps=60)
        in_process = s.coupled_bisect(460.0, 490.0, 0)
        args = (self.R_I, self.R_O, self.H, self.C11E, self.C12E,
                self.C13E, self.C33E, self.RHO,
                self.E31, self.E33, self.X11, self.X33,
                60, 0, 460.0, 490.0, 45)
        via_worker = _piezo_mono_coupled_ff_root_worker(args)
        self.assertEqual(in_process, via_worker)


class TestPiezoMonoMixedEdge(unittest.TestCase):
    """Paper 5 Option A mixed-edge C-F/F-C (Sec 18.208, same pattern as
    PiezoOutOfPlaneSolver's Sec 18.192 mixed-edge family). C-F = inner
    clamped, outer free; F-C = inner free, outer clamped. No open-
    circuit variant -- Phase 0's OC=SC theorem still covers this.

    Two kinds of gate here. (1) Reduction: _elastic_mixed_det/
    _coupled_mixed_det with both edges the same label must reproduce
    the already-validated elastic_det/elastic_cc_det/coupled_det/
    cc_coupled_det BIT-IDENTICALLY (same row order, same
    equilibration) -- this is the real correctness proof for the new
    row-selection logic, not just an assertion. (2) Physical anchors
    for the genuine C-F/F-C roots, captured 2026-09-20 at dps=60; do
    not retune. Same PZT-4 Duan Table 1 geometry as
    TestPiezoMonoElastic/TestPiezoMonoCoupled (r_i=0.1, r_o=0.6,
    H=0.01).
    """

    R_I, R_O, H = 0.1, 0.6, 0.01
    C11E, C12E, C13E, C33E = 132e9, 71e9, 73e9, 115e9
    E31, E33 = 4.1, 14.1
    X11, X33 = 7.124e-9, 5.841e-9
    RHO = 7500.0
    EL_CF = 1605.4427637454191
    CP_CF = 1637.9221255770667  # 2026-09-23 Sec 18.231 (consistent projection): duan 1644.8419769729373
    EL_FC = 556.8844158380273
    CP_FC = 564.9003213461377  # 2026-09-23 Sec 18.231 (consistent projection): duan 566.6020939564801

    def _make(self, dps=60):
        return PiezoMonolithicOutOfPlaneSolver(
            r_i=self.R_I, r_o=self.R_O, H=self.H,
            C11E=self.C11E, C12E=self.C12E, C13E=self.C13E,
            C33E=self.C33E, rho=self.RHO,
            e31=self.E31, e33=self.E33, X11=self.X11, X33=self.X33,
            dps=dps)

    def test_elastic_mixed_reduces_to_ff_and_cc(self):
        s = self._make(dps=40)
        for n in (0, 2):
            omega = 500.0 if n == 0 else 1500.0
            self.assertEqual(
                s.elastic_det(omega, n),
                s._elastic_mixed_det(omega, n, "F", "F"))
            self.assertEqual(
                s.elastic_cc_det(omega, n),
                s._elastic_mixed_det(omega, n, "C", "C"))

    def test_coupled_mixed_reduces_to_ff_and_cc(self):
        s = self._make(dps=40)
        for n in (0, 2):
            omega = 500.0 if n == 0 else 1500.0
            self.assertEqual(
                s.coupled_det(omega, n),
                s._coupled_mixed_det(omega, n, "F", "F"))
            self.assertEqual(
                s.cc_coupled_det(omega, n),
                s._coupled_mixed_det(omega, n, "C", "C"))

    def test_cf_anchor_stiffens(self):
        s = self._make(dps=60)
        el = s.elastic_cf_bisect(1590.0, 1620.0, 0)
        cp = s.cf_coupled_bisect(1636.0, 1640.0, 0)
        self.assertLess(abs(el - self.EL_CF) / self.EL_CF, 1e-8)
        self.assertLess(abs(cp - self.CP_CF) / self.CP_CF, 1e-8)
        drel = (cp - el) / el
        self.assertGreater(drel, 0.0)
        self.assertLess(drel, 0.05)

    def test_fc_anchor_stiffens(self):
        s = self._make(dps=60)
        el = s.elastic_fc_bisect(545.0, 570.0, 0)
        cp = s.fc_coupled_bisect(564.0, 565.5, 0)
        self.assertLess(abs(el - self.EL_FC) / self.EL_FC, 1e-8)
        self.assertLess(abs(cp - self.CP_FC) / self.CP_FC, 1e-8)
        drel = (cp - el) / el
        self.assertGreater(drel, 0.0)
        self.assertLess(drel, 0.05)

    def test_coupled_cf_refuses_missing_e31(self):
        s = PiezoMonolithicOutOfPlaneSolver(
            r_i=self.R_I, r_o=self.R_O, H=self.H,
            C11E=self.C11E, C12E=self.C12E, C13E=self.C13E,
            C33E=self.C33E, rho=self.RHO, dps=30)
        with self.assertRaises(ValueError):
            s.cf_coupled_det(1644.0, 0)

    def test_worker_scalar_reconstruction_matches_in_process(self):
        s = self._make(dps=60)

        in_el_cf = s.elastic_cf_bisect(1590.0, 1620.0, 0, iters=50)
        args_el_cf = (self.R_I, self.R_O, self.H, self.C11E, self.C12E,
                      self.C13E, self.C33E, self.RHO,
                      60, 0, 1590.0, 1620.0, 50)
        self.assertEqual(in_el_cf,
                          _piezo_mono_elastic_cf_root_worker(args_el_cf))

        in_el_fc = s.elastic_fc_bisect(545.0, 570.0, 0, iters=50)
        args_el_fc = (self.R_I, self.R_O, self.H, self.C11E, self.C12E,
                      self.C13E, self.C33E, self.RHO,
                      60, 0, 545.0, 570.0, 50)
        self.assertEqual(in_el_fc,
                          _piezo_mono_elastic_fc_root_worker(args_el_fc))

        in_cp_cf = s.cf_coupled_bisect(1636.0, 1640.0, 0, iters=45)
        args_cp_cf = (self.R_I, self.R_O, self.H, self.C11E, self.C12E,
                      self.C13E, self.C33E, self.RHO,
                      self.E31, self.E33, self.X11, self.X33,
                      60, 0, 1636.0, 1640.0, 45)
        self.assertEqual(in_cp_cf,
                          _piezo_mono_coupled_cf_root_worker(args_cp_cf))

        in_cp_fc = s.fc_coupled_bisect(564.0, 565.5, 0, iters=45)
        args_cp_fc = (self.R_I, self.R_O, self.H, self.C11E, self.C12E,
                      self.C13E, self.C33E, self.RHO,
                      self.E31, self.E33, self.X11, self.X33,
                      60, 0, 564.0, 565.5, 45)
        self.assertEqual(in_cp_fc,
                          _piezo_mono_coupled_fc_root_worker(args_cp_fc))




class TestPiezoMonoDrivenForceSensing(unittest.TestCase):
    """Paper 5 Option A force-driven sensing admittance (Sec 18.21x,
    PAPER5_YOMEGA_SENSE_DERIVATION.md). Redirected here from a
    voltage-driven Y(omega) after Sec 18.210 proved that path
    identically trivial (Y=j*omega*C0, k_eff^2==0) on this homogeneous
    monolithic stacking -- a direct corollary of the closed Phase 0
    OC==SC theorem. This is the reciprocal finding read from the other
    electrical port: a harmonic radial ring load F at r_i < r_F < r_o
    on the still-grounded, still-fully-electroded F-F SC ring drives a
    nonzero SEGMENT charge Q_Omega = 4*H*Xi11_bar*(r_b*phibar'(r_b) -
    r_a*phibar'(r_a)) even though the FULL-FACE charge is identically
    zero (same D_z(H) pointwise cancellation as the voltage-driven
    case, read backwards).

    Three pre-registered gates, all captured 2026-09-20 at dps=60,
    same PZT-4 Duan Table 1 geometry as TestPiezoMonoCoupled
    (r_i=0.1, r_o=0.6, H=0.01), r_F=0.35, r_star=0.30:

    (1) Global force-balance identity (G):
        2*pi*int_{r_i}^{r_o} r*w(r) dr = -F/(A2*omega^2)
    is the load-bearing SIGN gate on the Q_r jump row (RHS =
    -F/(2*pi*r_F) on the Q_r-bracket row only). rel_err ~ 1e-61 at
    dps=60 -- essentially exact, not a fit.

    (2) The driven 12x12, evaluated as a homogeneous determinant (not
    via the RHS-solve path), vanishes at EXACTLY the same omega as the
    already-validated coupled_bisect F-F SC fundamental (473.2496...,
    Sec 18.199-18.209): a unique-continuation argument, since a
    homogeneous fully-continuous 12-vector solution at F=0 is just a
    restriction of one global 6x6 eigenmode. rel_err = 0.0 exactly.

    (3) As e31_bar -> 0, Q_in -> 0 linearly. IMPORTANT: e31_bar =
    e31 - (C13E/C33E)*e33 (PAPER5_DERIVATION.md Sec 1), so scaling e31
    alone while leaving e33 fixed does NOT drive e31_bar to zero at
    this material point -- it plateaus near -8.95 instead of 0,
    because e31_bar is dominated by the e33 cross term. That was a bug
    in an earlier interactive probe (not in this module); e31 and e33
    must be scaled TOGETHER to test this limit. Once done correctly,
    the two "mechanical" chi-cubic branches converge to exactly
    +-elastic mu (10.250145149520318... at omega=200, matching
    _elastic_k to 17 digits) and Q_in vanishes with log-log slope ->
    1.0 in the deep-perturbative regime (scale <= 1e-2 here; at
    scale=1.0 e31_bar is O(1), not a small perturbation, so the first
    decade alone does not yet show the asymptotic linear rate).
    """

    R_I, R_O, H = 0.1, 0.6, 0.01
    C11E, C12E, C13E, C33E = 132e9, 71e9, 73e9, 115e9
    E31, E33 = 4.1, 14.1
    X11, X33 = 7.124e-9, 5.841e-9
    RHO = 7500.0
    OMEGA = 200.0
    R_F = 0.35
    R_STAR = 0.30
    # 2026-09-23 Sec 18.231 (consistent projection): duan Y_SENSE_AT_200 6.050054141033202e-12,
    # CP_FF 473.2496088764243; Q_segment weight 4H -> 4 pi H/3.
    Y_SENSE_AT_200 = 5.00405914569838e-12
    CP_FF = 471.36460219209835

    def _make(self, e31=None, e33=None, dps=60):
        return PiezoMonolithicOutOfPlaneSolver(
            r_i=self.R_I, r_o=self.R_O, H=self.H,
            C11E=self.C11E, C12E=self.C12E, C13E=self.C13E,
            C33E=self.C33E, rho=self.RHO,
            e31=self.E31 if e31 is None else e31,
            e33=self.E33 if e33 is None else e33,
            X11=self.X11, X33=self.X33, dps=dps)

    def test_identity_G_force_balance(self):
        """2*pi*int r*w dr = -F/(A2*omega^2), the sign gate on the
        Q_r jump row. Exact to dps precision, not a fit."""
        from mpmath import mp, mpf, quad
        s = self._make(dps=60)
        res = s.driven_ff_force_sc(mpf(self.OMEGA), self.R_F, F=1.0, n=0)
        with mp.workdps(60):
            def w_of_r(r):
                w, _ = s._driven_eval(res, r)
                return r * w
            lhs = 2 * mp.pi * quad(w_of_r, [s.r_i, self.R_F, s.r_o])
            A2 = s._A2()
            rhs = -mpf(1.0) / (A2 * mpf(self.OMEGA) ** 2)
            rel_err = abs((lhs - rhs) / rhs)
        self.assertLess(float(rel_err), 1e-30)

    def test_pole_matches_coupled_bisect(self):
        """Driven 12x12 homogeneous det vanishes at the already-
        validated F-F SC fundamental (unique-continuation gate)."""
        from mpmath import mp, mpf, matrix
        s = self._make(dps=60)
        cp = s.coupled_bisect(460.0, 490.0, 0)
        self.assertLess(abs(cp - self.CP_FF) / self.CP_FF, 1e-8)

        def driven12_det(omega):
            with mp.workdps(60):
                omega = mpf(omega)
                d = s._d(); A1v = s._A1(); Hm = mpf(s.H)
                e31_bar = s._e31_bar()
                K_pref = s._K_pref()  # Sec 18.231
                lams = s._branches(omega)

                def vecs(r):
                    return s._driven_vecs_at(0, r, lams, d, A1v, K_pref)

                w_i, wp_i, phi_i, m_i, q_i, phip_i = vecs(s.r_i)
                w_o, wp_o, phi_o, m_o, q_o, phip_o = vecs(s.r_o)
                w_F, wp_F, phi_F, m_F, q_F, phip_F = vecs(self.R_F)
                zero6 = [mpf(0)] * 6

                def interface_row(vec):
                    return list(vec) + [-x for x in vec]

                rows = [
                    m_i + zero6, q_i + zero6, phip_i + zero6,
                    zero6 + m_o, zero6 + q_o, zero6 + phip_o,
                    interface_row(w_F), interface_row(wp_F),
                    interface_row(m_F), interface_row(q_F),
                    interface_row(phi_F), interface_row(phip_F),
                ]
                M = matrix(rows)
                for i in range(12):
                    sc = max(abs(M[i, j]) for j in range(12)) or mpf(1)
                    for j in range(12):
                        M[i, j] = M[i, j] / sc
                return mp.det(M)

        with mp.workdps(60):
            a, b = mpf(460.0), mpf(490.0)
            flo = driven12_det(a)
            for _ in range(60):
                m = (a + b) / 2
                fm = driven12_det(m)
                if (fm.real > 0) == (flo.real > 0):
                    a, flo = m, fm
                else:
                    b = m
            omega_driven12 = float((a + b) / 2)
        self.assertLess(abs(omega_driven12 - cp) / cp, 1e-10)

    def test_qin_vanishes_linearly_as_e31_bar_to_zero(self):
        """Scaling e31 AND e33 together drives e31_bar -> 0 and Q_in
        -> 0 with log-log slope -> 1 (linear), once past the O(1)
        non-perturbative point at scale=1.0."""
        from mpmath import mp, mpf
        qvals = []
        for scale in (1e-1, 1e-2, 1e-3, 1e-4):
            s2 = self._make(e31=self.E31 * scale, e33=self.E33 * scale,
                             dps=60)
            res2 = s2.driven_ff_force_sc(mpf(self.OMEGA), self.R_F,
                                          F=1.0, n=0)
            with mp.workdps(60):
                q_in = s2.Q_segment(res2, s2.r_i, self.R_STAR)
            qvals.append(complex(q_in).real)
        ratios = [qvals[i] / qvals[i + 1] for i in range(len(qvals) - 1)]
        for r in ratios:
            self.assertGreater(r, 9.9)
            self.assertLess(r, 10.1)

    def test_mechanical_branches_match_elastic_mu_as_e31_bar_to_zero(self):
        """The two non-escaping chi-cubic branches converge to exactly
        +-mu (the pure elastic _elastic_k value) as e31_bar -> 0, once
        e31 and e33 are scaled together correctly."""
        from mpmath import mp, mpf
        s0 = self._make(dps=80)
        with mp.workdps(80):
            mu, _ = s0._elastic_k(mpf(self.OMEGA))
        s2 = self._make(e31=self.E31 * 1e-8, e33=self.E33 * 1e-8, dps=80)
        with mp.workdps(80):
            lams = s2._branches(mpf(self.OMEGA))
        lam_small = sorted([complex(l) for _, l in lams],
                            key=lambda z: abs(z))[:2]
        got = sorted(abs(z) for z in lam_small)
        self.assertLess(abs(got[0] - float(mu)) / float(mu), 1e-9)
        self.assertLess(abs(got[1] - float(mu)) / float(mu), 1e-9)

    def test_worker_matches_in_process(self):
        from mpmath import mpf
        s = self._make(dps=60)
        y_direct = complex(s.Y_sense(mpf(self.OMEGA), self.R_F,
                                      self.R_STAR, F=1.0, n=0))
        args = (self.R_I, self.R_O, self.H, self.C11E, self.C12E,
                self.C13E, self.C33E, self.RHO,
                self.E31, self.E33, self.X11, self.X33,
                60, self.OMEGA, self.R_F, self.R_STAR, 1.0, 0)
        re, im = _piezo_mono_driven_force_sc_worker(args)
        self.assertEqual(re, y_direct.real)
        self.assertEqual(im, y_direct.imag)
        self.assertAlmostEqual(re, self.Y_SENSE_AT_200, delta=1e-20)

    def test_driven_force_sc_ff_reduces_to_driven_ff(self):
        """G_reduce: driven_force_sc('F','F') matches driven_ff_force_sc
        bit for bit. Same rim-row order and the same equilibrated solve,
        so the branch constants and the segment charge agree exactly."""
        from mpmath import mpf
        s = self._make(dps=40)
        ff = s.driven_ff_force_sc(mpf(self.OMEGA), self.R_F, F=1.0, n=0)
        gen = s.driven_force_sc(
            mpf(self.OMEGA), self.R_F, "F", "F", F=1.0, n=0)
        self.assertEqual(list(ff["c_I"]), list(gen["c_I"]))
        self.assertEqual(list(ff["c_II"]), list(gen["c_II"]))
        self.assertEqual(
            s.Q_segment(ff, s.r_i, self.R_STAR),
            s.Q_segment(gen, s.r_i, self.R_STAR))

    def test_driven_force_sc_clamped_changes_the_solution(self):
        """A clamped rim has to change the branch constants. This fails
        if inner/outer are ignored and the builder always emits F-F rows."""
        from mpmath import mpf
        s = self._make(dps=25)
        ff = s.driven_force_sc(
            mpf(self.OMEGA), self.R_F, "F", "F", F=1.0, n=0)
        cc = s.driven_force_sc(
            mpf(self.OMEGA), self.R_F, "C", "C", F=1.0, n=0)
        self.assertNotEqual(list(ff["c_I"]), list(cc["c_I"]))

    def test_driven_force_sc_rejects_bad_bc_and_n(self):
        s = self._make(dps=15)
        with self.assertRaises(ValueError):
            s.driven_force_sc(self.OMEGA, self.R_F, "X", "F")
        with self.assertRaises(NotImplementedError):
            s.driven_force_sc(self.OMEGA, self.R_F, "F", "F", n=1)


class TestPiezoFreeEdgeEffectiveShear(unittest.TestCase):
    """2026-09-22 (LESSONS Sec 18.226): the free-edge shear row must be the
    Kirchhoff effective shear V_r = Q_r + (1/r) dM_rtheta/dtheta, not Q_r.
    Q_r alone is exact only at n=0, which is why no n=0 gate caught it. This
    checks both piezo solvers' ELASTIC paths against ring_disk's independent
    Frobenius _Lmat_mp (FE-confirmed, Sec 18.38) at n=2, for F-F / C-F / F-C:
      * F-F n=0 and n=2 roots, in native Omega, equal to the Sec 18.38
        anchors to 1e-6 (the Q_r-only row gave 0.1001086 at n=2, -7.5%);
      * C-F and F-C n=2 roots are clean ring_disk sign-flip roots.
    Isotropic R_i/R_o=0.5, nu=0.3, b=2, h=0.04 (HALF-thickness), E=rho=1.
    """

    @classmethod
    def setUpClass(cls):
        solver, geom, mat = make_annulus_solver(0.5, nu=0.30)
        cls.ring = solver
        nu = mpf("0.3")
        c11 = 1 / (1 - nu ** 2)
        c66 = 1 / (2 * (1 + nu))
        L = 2 * geom.b / mp.pi
        cls.conv = ((mp.sqrt(c66) / L)
                    / (L * mp.sqrt(3 * c66 / (geom.h ** 2 * c11))))
        cls.mono = PiezoMonolithicOutOfPlaneSolver(
            geom.R_i, geom.R_o, geom.h, C11E=c11, C12E=nu * c11, C13E=0,
            C33E=1, rho=1, dps=30)
        cls.lay = PiezoOutOfPlaneSolver(
            float(geom.R_i), float(geom.R_o), float(geom.h), 1.0, 0.3, 1.0,
            h1=0.0, dps=30)

    def _root(self, f, n, lo, hi):
        conv = self.conv
        a, b = mpf(lo), mpf(hi)
        fa = f(a * conv, n).real
        self.assertNotEqual(fa > 0, f(b * conv, n).real > 0,
                            "no sign change in [%s, %s]" % (lo, hi))
        for _ in range(50):
            m = (a + b) / 2
            fm = f(m * conv, n).real
            if (fm > 0) == (fa > 0):
                a, fa = m, fm
            else:
                b = m
        return (a + b) / 2

    def test_ff_matches_section_18_38_anchors(self):
        from plate_solver.ring_disk import (SECTION_18_38_N0_OM,
                                            SECTION_18_38_N2_OM)
        for f in (self.mono.elastic_det,
                  lambda w, n: self.lay.elastic_det(w, n, 0.0)):
            for n, ref in ((0, SECTION_18_38_N0_OM), (2, SECTION_18_38_N2_OM)):
                om = self._root(f, n, ref * 0.97, ref * 1.03)
                self.assertLess(abs(float(om) / ref - 1), 1e-6,
                                "F-F n=%d root %s vs ring_disk %s" % (n, om, ref))

    def test_mixed_edges_are_ring_disk_roots_at_n2(self):
        from plate_solver.ring_disk import signflip_ladder, _det_fn_for
        cases = (("C", "F", self.mono.elastic_cf_det, 0.35, 0.40),
                 ("F", "C", self.mono.elastic_fc_det, 0.78, 0.84),
                 ("C", "F", lambda w, n: self.lay.elastic_cf_det(w, n, 0.0), 0.35, 0.40),
                 ("F", "C", lambda w, n: self.lay.elastic_fc_det(w, n, 0.0), 0.78, 0.84))
        for bi, bo, f, lo, hi in cases:
            om = self._root(f, 2, lo, hi)
            lad = signflip_ladder(_det_fn_for(self.ring, 2, bi, bo, False), 2,
                                  float(om), deltas=(1e-3, 1e-4, 1e-5))
            self.assertTrue(lad.clean, "%s-%s n=2 root %s is not a ring_disk root"
                            % (bi, bo, om))

class TestPiezoConsistentProjection(unittest.TestCase):
    """2026-09-23, LESSONS Sec 18.231 /
    SC_PROJECTION_CONSISTENT_DERIVATION_2026-09-23.md. Pre-registered gates:
    P2 long-wave SC factor (consistent 1, galerkin_sine 96/pi^4, duan
    12/pi^2) for BOTH piezo modules; P0 legacy 'duan' still reproduces the
    pre-fix anchors; OC charge arm (h+h1) -> (h+h1/2) scales the linear OC
    split by (h+h1/2)/(h+h1) to first order; Q_segment weight."""

    K5 = dict(r_i=0.1, r_o=0.6, H=0.01, C11E=132e9, C12E=71e9, C13E=73e9,
              C33E=115e9, rho=7500.0, e31=4.1, e33=14.1, X11=7.124e-9,
              X33=5.841e-9)
    K4 = dict(r_i=0.1, r_o=0.6, h=0.01, E=200e9, nu=0.3, rho=7800.0,
              h1=0.01 / 6, C11E=132e9, C12E=71e9, C13E=73e9, C33E=115e9,
              e31=4.1, e33=14.1, X11=7.124e-9, X33=5.841e-9, rho_pzt=7500.0)

    def _targets(self):
        from mpmath import mp
        return {'consistent': 1.0, 'galerkin_sine': float(96 / mp.pi ** 4),
                'duan': float(12 / mp.pi ** 2)}

    def test_long_wave_factor_monolithic(self):
        for pr, tgt in self._targets().items():
            s = PiezoMonolithicOutOfPlaneSolver(projection=pr, dps=50, **self.K5)
            got = float(s.long_wave_sc_ratio(1e-5))
            self.assertLess(abs(got - tgt) / tgt, 1e-8, "%s %s" % (pr, got))

    def test_long_wave_factor_layered(self):
        for pr, tgt in self._targets().items():
            s = PiezoOutOfPlaneSolver(projection=pr, dps=50, **self.K4)
            got = float(s.long_wave_sc_ratio(1e-3))
            self.assertLess(abs(got - tgt) / tgt, 1e-8, "%s %s" % (pr, got))

    def test_duan_legacy_reproduces_pre_fix_anchors(self):
        m = PiezoMonolithicOutOfPlaneSolver(projection='duan', dps=60, **self.K5)
        cp = m.coupled_bisect(460.0, 490.0, 0)
        self.assertLess(abs(cp - 473.2496088764243) / 473.2496088764243, 1e-8)
        p = PiezoOutOfPlaneSolver(projection='duan', dps=60, **self.K4)
        oc = p.oc_ff_bisect(720.0, 800.0, 0)
        self.assertLess(abs(oc - 750.2394210586135) / 750.2394210586135, 1e-6)
        ocs = p.oc_ff_sine_bisect(720.0, 780.0, 0, iters=45)
        self.assertLess(abs(ocs - 750.1950280038534) / 750.1950280038534, 1e-6)

    def test_oc_charge_arm_ratio(self):
        el = 747.9863950952586
        dn = PiezoOutOfPlaneSolver(projection='duan', dps=40, **self.K4)
        cs = PiezoOutOfPlaneSolver(dps=40, **self.K4)
        r = ((cs.oc_ff_bisect(720.0, 800.0, 0) - el)
             / (dn.oc_ff_bisect(720.0, 800.0, 0) - el))
        h, h1 = 0.01, 0.01 / 6
        arm = (h + h1 / 2) / (h + h1)
        self.assertLess(abs(r - arm) / arm, 2e-3, "ratio %s vs %s" % (r, arm))

    def test_q_segment_weight(self):
        import math
        c = PiezoMonolithicOutOfPlaneSolver(dps=30, **self.K5)
        d = PiezoMonolithicOutOfPlaneSolver(projection='duan', dps=30, **self.K5)
        self.assertAlmostEqual(float(c._face_flux_weight()),
                               4 * math.pi * 0.01 / 3, delta=1e-15)
        self.assertAlmostEqual(float(d._face_flux_weight()), 0.04, delta=1e-15)

    def test_bad_projection_refused(self):
        with self.assertRaises(ValueError):
            PiezoMonolithicOutOfPlaneSolver(projection='sine', **self.K5)
        with self.assertRaises(ValueError):
            PiezoOutOfPlaneSolver(projection='sine', **self.K4)



class TestPiezoDiskSC(unittest.TestCase):
    """Paper 4 lever 11 -- Liu 2002 solid-disk SC CPT (Tables 2 & 6).

    Full Gate 2 (12 modes, from-scratch scan, dps=100):
    validation/paper4_piezo/probe_piezo_p4_liu_disk_gate2_2026-09-24.py.
    These package tests use dps=40 for CI smoke; tolerance is the Gate 2
    bar (0.1% vs Liu x_kir; 0.1% vs Paper 3 lambda^2 for elastic), plus a
    tight regression on the default ('consistent') roots and on the
    legacy 'duan' roots.
    Disk API is main-process only (not wired into workers).
    """

    # Paper 3 / Gate 2a lambda^2 refs (Leissa flexural)
    LAMBDA2_P3 = {
        ('C', 0, 1): 10.215815,
        ('C', 1, 1): 21.260405,
        ('S', 0, 1): 4.935127,
        ('S', 1, 1): 13.898176,
    }
    # Gate 2 SC roots at dps=100 with LIU_DISK_KWARGS, re-anchored 2026-09-24
    # for e31 = -4.1 (Liu Table 1 as printed, LESSONS Sec 18.244; the +4.1
    # anchors were 902.4165163374934 / 902.4186619938962 and
    # 435.59981690040496 / 435.60128907622254).
    # (table, bc, n, m): (omega_consistent, omega_duan, omega_kir, lo, hi)
    LIU_CPT = {
        (2, 'C', 0, 1): (902.4785335853126, 902.494065168716, 902.5, 880.0, 930.0),
        (6, 'S', 0, 1): (435.6423672713382, 435.6530233280846, 435.6, 420.0, 450.0),
    }
    R0 = 0.6
    H_FULL = 0.02  # 2h
    E, NU, RHO = 200e9, 0.3, 7800.0

    @staticmethod
    def _lambda2_from_omega(omega, a=0.6, H=0.02, E=200e9, nu=0.3, rho=7800.0):
        D = E * H ** 3 / (12.0 * (1.0 - nu * nu))
        return float(omega) * (a * a) * ((rho * H) / D) ** 0.5

    def test_flag_disk_6x6_at_ri0_is_not_the_disk(self):
        self.assertTrue(PiezoDiskSC.DISK_6X6_AT_RI0_IS_NOT_THE_DISK)

    def test_refuse_ri_gt_0(self):
        disk = make_liu_disk(dps=30, h1=0.0)
        disk.r_i = 0.1
        with self.assertRaises(PiezoDiskPathError):
            disk.elastic_disk_det(100.0, 0, 'C')
        disk2 = make_liu_disk(dps=30)
        disk2.r_i = 0.05
        with self.assertRaises(PiezoDiskPathError):
            disk2.sc_disk_det(900.0, 0, 'C')

    def test_annular_solver_refuses_ri_zero(self):
        from plate_solver.piezo_solver import PiezoOutOfPlaneSolver
        for ri in (0.0, -0.1):
            with self.assertRaises(ValueError):
                PiezoOutOfPlaneSolver(r_i=ri, r_o=0.6, **LIU_DISK_KWARGS)

    def test_sc_disk_refuses_h1_zero(self):
        disk = make_liu_disk(dps=30, h1=0.0)
        with self.assertRaises(ValueError):
            disk.sc_disk_det(900.0, 0, 'C')

    def test_elastic_disk_lambda2_vs_paper3(self):
        """Bare-host elastic 2x2 recovers Paper 3 lambda^2 (couple C/S)."""
        disk = make_liu_disk(dps=40, h1=0.0)
        cases = [
            ('C', 0, 1, 800.0, 950.0),
            ('C', 1, 1, 1700.0, 1950.0),
            ('S', 0, 1, 380.0, 460.0),
            ('S', 1, 1, 1100.0, 1280.0),
        ]
        for bc, n, m, lo, hi in cases:
            w = disk.elastic_disk_bisect(lo, hi, n, bc_outer=bc, h1=0.0, iters=40)
            lam = self._lambda2_from_omega(w)
            ref = self.LAMBDA2_P3[(bc, n, m)]
            rel = abs(lam - ref) / ref
            self.assertLess(
                rel, 1e-3,
                "FAIL_ELASTIC_DISK %s n=%s m=%s lam=%s ref=%s rel=%s"
                % (bc, n, m, lam, ref, rel))

    def test_liu_table2_and_table6_n0_m1_within_0p1pct(self):
        """Primary Gate 2 CPT points: Table 2 and Table 6, n=0 m=1."""
        disk = make_liu_disk(dps=40)  # default projection='consistent'
        self.assertEqual(disk.projection, 'consistent')
        duan = make_liu_disk(dps=40, projection='duan')
        for (table, bc, n, m), (omega_c, omega_d, omega_kir, lo, hi) in self.LIU_CPT.items():
            w = disk.sc_disk_bisect(lo, hi, n, bc_outer=bc, iters=40)
            rel_kir = abs(w - omega_kir) / omega_kir
            self.assertLess(
                rel_kir, 1e-3,
                "FAIL_LIU_CPT T%s %s n=%s m=%s omega=%s kir=%s rel=%s"
                % (table, bc, n, m, w, omega_kir, rel_kir))
            # Regression on the dps=100 Gate 2 roots (both projections)
            rel_c = abs(w - omega_c) / omega_c
            self.assertLess(
                rel_c, 1e-7,
                "FAIL_GATE2_REGRESSION T%s omega=%s ref=%s rel=%s"
                % (table, w, omega_c, rel_c))
            wd = duan.sc_disk_bisect(lo, hi, n, bc_outer=bc, iters=40)
            rel_d = abs(wd - omega_d) / omega_d
            self.assertLess(
                rel_d, 1e-7,
                "FAIL_GATE2_DUAN_REGRESSION T%s omega=%s ref=%s rel=%s"
                % (table, wd, omega_d, rel_d))


class TestE31PrintedSign(unittest.TestCase):
    """2026-09-24, LESSONS Sec 18.244/18.245: the series material is Liu 2002 /
    Duan 2005 Table 1 PZT-4 exactly as printed, e31 = -4.1 C/m^2.

    The older piezo test classes above keep e31 = +4.1 on purpose: they pin
    code paths to the anchors captured when those paths were written, and the
    code does not depend on the sign. THIS class pins the physics of the sign:
      (1) the default bundle carries -4.1;
      (2) Duan 2005 Table 4 (C-C, p=0, first radial mode, projection 'duan',
          Duan's own closure) is reproduced within 0.015% at all three
          h1/2h ratios at -4.1, where +4.1 misses 1/5 by 0.042%;
      (3) the -4.1 headline numbers FE-confirmed by job 2532086: Paper 5 F-F
          SC root 519.8388116516741 rad/s (split +12.42%, FE +12.43%) and the
          Paper 4 F-F OC/SC split at 1/12, 1.9443% (FE 1.9439%)."""

    MAT = dict(C11E=132e9, C12E=71e9, C13E=73e9, C33E=115e9, e33=14.1,
               X11=7.124e-9, X33=5.841e-9)

    def _p4(self, den, e31, projection='consistent', dps=30):
        return PiezoOutOfPlaneSolver(
            r_i=0.1, r_o=0.6, h=0.01, E=200e9, nu=0.3, rho=7800.0,
            h1=0.02 / den, e31=e31, rho_pzt=7500.0, dps=dps,
            projection=projection, **self.MAT)

    def test_default_bundle_is_printed_sign(self):
        self.assertEqual(LIU_DISK_KWARGS['e31'], -4.1)

    def test_duan_table4_cc_fundamental_prefers_printed_sign(self):
        targets = {12: 2792.0, 8: 2853.0, 5: 2989.0}
        for den, tg in targets.items():
            w = self._p4(den, -4.1, projection='duan').cc_coupled_bisect(
                tg * 0.97, tg * 1.03, 0, iters=40)
            self.assertLess(abs(w / tg - 1.0), 1.5e-4,
                            "FAIL_DUAN_T4 1/%d omega=%s target=%s" % (den, w, tg))
        wp = self._p4(5, 4.1, projection='duan').cc_coupled_bisect(
            2989.0 * 0.97, 2989.0 * 1.03, 0, iters=40)
        self.assertGreater(abs(wp / 2989.0 - 1.0), 3.5e-4,
                           "+4.1 unexpectedly matches Duan 1/5: %s" % wp)

    def test_paper5_ff_sc_root_at_printed_sign(self):
        s = PiezoMonolithicOutOfPlaneSolver(
            0.1, 0.6, 0.01, 132e9, 71e9, 73e9, 115e9, 7500.0, e31=-4.1,
            e33=14.1, X11=7.124e-9, X33=5.841e-9, dps=30)
        w = s.coupled_bisect(515.0, 525.0, 0)
        self.assertLess(abs(w - 519.8388116516741) / 519.8388116516741, 1e-9)

    def test_paper4_ff_oc_sc_split_at_printed_sign(self):
        s = self._p4(12, -4.1)
        el = s.elastic_bisect(600.0, 900.0, 0, iters=50)
        sc = s.coupled_bisect(el * 0.99, el * 1.01, 0, iters=45)
        oc = s.oc_ff_bisect(el * 0.999, el * 1.10, 0, iters=50)
        self.assertLess(abs((oc / sc - 1.0) - 0.019443190102619523), 1e-8)


if __name__ == "__main__":
    unittest.main()


