# -*- coding: utf-8 -*-
"""
validation/test_validated_tables.py -- automated regression gate for the VALIDATED
result tables (LESSONS_LEARNED Sec 17 item 18 / Sec 14 item 4).
Shipped path is validation/, not tests/ (github_repo/tests/ has only test_solver.py).

WHY: silent-revert drift has now bitten THREE times (job 2315355's un-patched
tree, tests/test_solver.py's stale gates, job 2315890's missing residual
screen -- Sec 13.4). tests/test_solver.py pins single-point sigma_min anchors;
NOTHING pins the validated TABLES, so a refactor can revert them without any
test noticing. This file adds:

  ALWAYS-RUN (seconds; run after EVERY deploy, alongside test_solver):
    * SOLVER_VERSION pin (EXPECT_SOLVER_VERSION env to bump intentionally)
    * detectors presence: both residual-screen functions AND the
      `from .boundary import EdgeKind` import line (the A.21 narrow drift)
    * gap-refinement defaults in source: default-on ("1"), cap 30 (A.24-26)
    * Shi orthotropic constants T/R/mu_theta (Sec 16)
    * free-free BC registry intact
  GATE_FFP1=1 (~1.5-3 h on >=16 cores; the gate that CAUGHT the s7
  regression, job 2316094):
    * fully-default find_modes_sigmin on FF-P1 OOP must reproduce ALL 13
      pinned raw Omegas (incl. the gap-refine-recovered 0.961312 AND the
      2026-07-11 addendum's two window-extension modes, 2.320741/2.602500)
      to |d| < 0.003, nearest-neighbour.
  GATE_SHI=1 (~3-5 h on >=16 cores):
    * fully-default Shi orthotropic FFFF discovery must contain modes
      nearest-matching all 8 closed REAL targets' computed values (Sec 16,
      job 2316501, UPDATED 2026-07-11 addendum Sec B: 66.76 resolved REAL,
      job 2316744) to <1.0% in Omega_lit. Table is now 8/8 REAL; no
      AMBIGUOUS candidate remains to track separately.

  GATE_FFP1_IP=1 (~5-7 h on >=16 cores; ADDED 2026-07-12 from job 2317336's
  verbatim capture):
    * fully-default find_modes_sigmin on FF-P1 IP (window extended to 2.85
      with density-preserving n_scan, exactly the capture config) must
      reproduce ALL 43 pinned scan locations to |d| < 0.003.

  TODO (deliberately NOT implemented rather than faked):
    * pinning the rect-IP 6-fundamental table needs value capture from
      run_rect_ip_validation, which currently prints rather than returns its
      comparison -- add a return path first, then pin here.

Run:            python -m unittest tests.test_validated_tables -v
Full gates:     GATE_FFP1=1 GATE_SHI=1 GATE_FFP1_IP=1 python -m unittest tests.test_validated_tables -v
(submit_validated_tables_gate.sh wraps this for the cluster. All three
gates passed together on 2026-07-12, job 2318033 + capture job 2317336.)

Lesson baked in from the s6->s7->s8 chain (Sec 17 item 12): a completeness/
threshold change that passes ONE geometry's gate is not validated -- BOTH
gates here must pass before any detector-completeness change ships.
"""
from __future__ import annotations
import os
import sys
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from mpmath import mp
import plate_solver as ps
from plate_solver import detectors as det

EXPECT_VER = os.environ.get("EXPECT_SOLVER_VERSION", "2026-07-10.s10")
GATE_FFP1 = os.environ.get("GATE_FFP1", "0") == "1"
GATE_SHI = os.environ.get("GATE_SHI", "0") == "1"
GATE_FFP1_IP = os.environ.get("GATE_FFP1_IP", "0") == "1"

# ---- pinned tables (post-hoc references ONLY; never seed a computation) ----
# Sec 10.8, jobs 2315385/2315829, s5-values re-confirmed under s8 (2316500).
# UPDATED 2026-07-11 addendum Sec C.1 (jobs 2316746/2316935/2316939): the
# window-2.80 extension found two further confirmed REAL modes (raw
# 2.320741 -> Omega_lit 91.6192, Ansys 91.452 0.18%; raw 2.602500 -> Omega_lit
# 102.7426, Ansys 102.136 0.594%) -- all 8 Ansys iso-FFFF targets now matched.
# NOTE (final status 2026-07-12, jobs 2317684/2318031): the ceiling item is
# CLOSED as a documented TWO-CANDIDATE ASSIGNMENT AMBIGUITY above the pinned
# window. The unmatched FE mode Omega_lit=111.4627 is ODD (half-sector Ansys,
# calibration 4/4); the solver has TWO parity-consistent ODD dips there --
# raw 2.828013 (0.16% away, rM=0.119 marginal) and raw 2.838515 (0.52% away,
# rM=0.027 clean) -- while 2.817613 is EVEN (its ARTIFACT call is thereby
# parity-CONFIRMED). Five independent techniques were exhausted (basis bump,
# sensitivity map, null-vector gap, parity readout, ODD-block topology scan);
# per the pre-registered O2 outcome, NEITHER candidate is forced into the
# table. The FE mode is NOT missing from the solver spectrum; only the
# assignment between the two ODD dips is open. No further cluster time
# without a genuinely new lever.
FFP1_OOP_RAW = [0.320670, 0.383695, 0.616645, 0.844488, 0.961312, 1.001419,
                1.369611, 1.602090, 1.752188, 2.121685, 2.320741, 2.379886,
                2.602500]
# Sec 16 final Shi table (job 2316501, s8): computed Omega_lit of the 7 REAL
# matches. UPDATED 2026-07-11 addendum Sec B (job 2316744): the 8th target
# (66.76), formerly stably AMBIGUOUS across 3 pre-session runs, resolved REAL
# via a larger-basis (n_dofs 20->28) + null-vector-anatomy probe -- combined
# dps=60/n28 run: Omega_lit=67.2235 (0.69% err). Shi orthotropic FFFF table is
# now 8/8; `SHI_AMBIG_LIT` is retired (folded into SHI_REAL_LIT below).
SHI_REAL_LIT = [12.5335, 15.6298, 36.3107, 43.8826, 67.2235, 88.0668, 90.3233,
                96.0312]

# ---- IP (Part 2) free-free: PINNED 2026-07-12 ------------------------------
# Source: job 2317336's verbatim capture (default s8 path, n_dofs=20, window
# (0.0071, 2.85) with density-preserving n_scan) -- the exact 43 accepted
# scan locations, NEVER reconstructed from memory. This is a REPRODUCIBILITY
# pin: what a fully-default rerun must find, INCLUDING dips later adjudicated
# as artifacts (they are stable features of the n20 landscape). Verdict
# provenance: n20 screens (job 2317336) + n28/n36 adjudications (jobs
# 2317685/2318032 + zone truths from jobs 2316938/2317313).
FFP1_IP_RAW_SCAN = [
    0.311045,  # ARTIFACT (Task E, permanently closed; Sec D.5)
    0.401154,  # REAL (FE 161.4 Hz, 0.03%)
    0.441229,  # ARTIFACT (n28 rTyy=17)
    0.541978,  # ARTIFACT
    0.687446,  # REAL -- NEW 2026-07-12 (FE 276.5 Hz, 0.007%)
    0.744954,  # ARTIFACT
    0.764697,  # REAL (FE 307.5 Hz, 0.03%)
    1.021209,  # ARTIFACT (n28 rTyy=11)
    1.142287,  # REAL (FE 459.5 Hz, 0.006% -- near-degenerate pair member)
    1.145158,  # REAL (FE 460.6 Hz, 0.006% -- the solver SPLITS this pair)
    1.197854,  # zone A: true mode at 1.197600 @n28 (FE 481.7 Hz, 0.005%)
    1.310594,  # ARTIFACT
    1.416689,  # REAL (FE 569.7 Hz, 0.03%)
    1.440568,  # REAL (FE 579.3 Hz, 0.03%)
    1.454636,  # ARTIFACT
    1.523230,  # REAL -- NEW 2026-07-12 (n28 1.522330; FE 612.3 Hz, 0.007%)
    1.598985,  # REAL -- NEW 2026-07-12 (FE 643.1 Hz, 0.012%)
    1.613555,  # zone B: true mode at 1.613123 @n28 (FE 648.8 Hz, 0.010%)
    1.616291,  # ARTIFACT (zone B partner; worsens with basis -- Sec D.2)
    1.664592,  # ARTIFACT (n28 rTyy=7.4)
    1.695387,  # ARTIFACT (true mode 1.690020 found by batch-5 n28 scan, FE 679.8 Hz 0.001%)
    1.779792,  # REAL -- NEW 2026-07-12 (FE 715.9 Hz, 0.001%)
    1.806063,  # ARTIFACT
    1.832694,  # REAL -- NEW 2026-07-12 (FE 736.8 Hz, 0.05%)
    1.843853,  # ARTIFACT
    1.975027,  # REAL (n28 1.972966; FE 793.6 Hz, 0.001%)
    1.987365,  # n20 shoulder of TRUE mode 1.987180 (batch-5; FE 799.3 Hz 0.003%; the earlier n28 golden mislocated to 1.988009)
    1.997627,  # REAL (n36 1.996235; FE 803.0 Hz, 0.004%)
    2.077957,  # ARTIFACT
    2.169751,  # REAL (n36 2.166815; FE 871.6 Hz, 0.002%)
    2.242363,  # ARTIFACT
    2.294399,  # REAL (n28 2.293214; FE 922.3 Hz, 0.013%)
    2.302746,  # REAL (n28 2.301512; FE 925.7 Hz, 0.007%)
    2.389113,  # REAL (FE 960.7 Hz, 0.03%)
    2.456143,  # n20 artifact (batch-5 n28 scan: only root in window is 2.450795, FE 985.7 Hz 0.011%)
    2.477799,  # ARTIFACT
    2.524851,  # zone D: true mode at 2.521489 @n28 (FE 1014.1 Hz, 0.014%)
    2.534250,  # ARTIFACT (n28/n36 rTyy 18/22 despite deep double-root sigma)
    2.570686,  # ARTIFACT (true mode 2.566362 found by batch-5 n28 scan, FE 1032.0 Hz 0.029%)
    2.626042,  # ARTIFACT
    2.679942,  # n20 artifact (batch-5 n28 scan: only root in window is 2.676026, FE 1075.9 Hz 0.043%)
    2.710103,  # ARTIFACT
    2.779320,  # ARTIFACT (old "Task C row" -- n20 mislocation; true mode below)
]
# Physical REAL table (adjudicated locations; documentation, NOT the gate):
# 0.401154, 0.687446, 0.764697, 1.142287, 1.145158, 1.197600, 1.416689,
# 1.440568, 1.522330, 1.598985, 1.613123, 1.690020, 1.779792, 1.832694,
# 1.972966, 1.987180, 1.996235, 2.166815, 2.293214, 2.301512, 2.389113,
# 2.450795, 2.521489, 2.566362, 2.676026, 2.771955.
# 26 REAL modes; ALL 27 FE PLANE183 modes (job 2316743) matched at <=0.052%,
# zero spurious REAL rows. The five 2026-07-12 batch-5 promotions (job
# 2318784: 1.690020/1.987180/2.450795/2.566362/2.676026, FE err
# 0.001-0.043%) are n28 fine-scan double roots that are INVISIBLE to the
# default n20 first pass -- the n20 landscape shows only nearby
# weak-enforcement artifacts (1.695387, 2.456143, 2.570686, 2.679942) or a
# mislocated shoulder (1.987365). This is the documented n20 completeness
# limitation driving the pending "blind first-pass at n28" default-config
# decision (probe_ip_blind_firstpass_n28.py; LESSONS Sec 20.8). The FE
# near-degenerate pair 1114.9/1115.3 Hz remains reproduced as ONE root
# (2.771955, Task-C relocation); the 459.5/460.6 Hz pair as TWO.


def _detectors_source():
    return open(os.path.join(os.path.dirname(det.__file__), "detectors.py"),
                "rb").read().decode("utf-8", "replace")


# Paper Tables 3--4 (annular cantilever, SOLVER_VERSION s10 regen).
# Always-run pin of the PRINTED numbers so a silent table rewrite is
# caught without a cluster discovery gate. These are post-hoc references
# only; they never seed a computation. Rectangular seok3/seok4 is out of
# Paper 1's printed scope and is not pinned here.
PAPER_CANTILEVER_OOP = [
    # r0/2b, 2Theta/pi, mode, computed, paper, err%
    (1.25, 0.25, 1, 0.345324, 0.341120, 1.232),
    (1.25, 0.25, 2, 1.018530, 1.013440, 0.502),
    (1.25, 0.25, 3, 1.722434, 1.718720, 0.216),
    (1.25, 0.50, 1, 0.099419, 0.099060, 0.363),
    (1.25, 0.50, 2, 0.343742, 0.340670, 0.902),
    (1.25, 0.50, 3, 0.605904, 0.604220, 0.279),
    (1.25, 0.75, 1, 0.050717, 0.050830, 0.221),
    (1.25, 0.75, 2, 0.151458, 0.149870, 1.060),
    (1.25, 0.75, 3, 0.376503, 0.376280, 0.059),
    (1.25, 1.00, 1, 0.033153, 0.033220, 0.201),
    (1.25, 1.00, 2, 0.080358, 0.079770, 0.738),
    (1.25, 1.00, 3, 0.243521, 0.243120, 0.165),
    (1.25, 1.25, 1, 0.024652, 0.024640, 0.048),
    (1.25, 1.25, 2, 0.049473, 0.049340, 0.269),
    (1.25, 1.25, 3, 0.151555, 0.151160, 0.261),
    (1.25, 1.50, 1, 0.019790, 0.019730, 0.305),
    (1.25, 1.50, 2, 0.034523, 0.034550, 0.078),
    (1.25, 1.50, 3, 0.097470, 0.097260, 0.216),
    (5 / 3, 1.00, 1, 0.018023, 0.018050, 0.148),
    (5 / 2, 1.00, 1, 0.007762, 0.007770, 0.102),
]
PAPER_CANTILEVER_IP = [
    (1.25, 0.25, 1, 0.341243, 0.345890, 1.343),
    (1.25, 0.25, 2, 0.808929, 0.819210, 1.255),
    (1.25, 0.25, 3, 0.940032, 0.940120, 0.009),
    (1.25, 0.50, 1, 0.120677, 0.120880, 0.168),
    (1.25, 0.50, 2, 0.347033, 0.351910, 1.386),
    (1.25, 0.50, 3, 0.529239, 0.529300, 0.012),
    (1.25, 1.00, 1, 0.039466, 0.039480, 0.035),
    (1.25, 1.00, 2, 0.106529, 0.105560, 0.918),
    (1.25, 1.00, 3, 0.263049, 0.261290, 0.673),
]


class TestPaperCantileverTables(unittest.TestCase):
    """Pin Paper 1 Tables 3--4 as printed. Cheap; no detector run."""

    def test_twenty_nine_rows_and_err_column(self):
        self.assertEqual(len(PAPER_CANTILEVER_OOP), 20)
        self.assertEqual(len(PAPER_CANTILEVER_IP), 9)
        comps = [r[3] for r in PAPER_CANTILEVER_OOP + PAPER_CANTILEVER_IP]
        self.assertEqual(len(set(comps)), 29, "duplicate computed frequency")
        for r0, tpi, mode, comp, pap, err in (
                PAPER_CANTILEVER_OOP + PAPER_CANTILEVER_IP):
            rec = 100.0 * abs(comp - pap) / abs(pap)
            self.assertLessEqual(
                abs(rec - err), 0.002,
                msg=f"r0={r0} 2T/pi={tpi} mode={mode}: "
                    f"printed err={err} recomputed={rec:.4f}")
        oop_max = max(r[5] for r in PAPER_CANTILEVER_OOP)
        ip_max = max(r[5] for r in PAPER_CANTILEVER_IP)
        self.assertAlmostEqual(oop_max, 1.232, places=3)
        self.assertAlmostEqual(ip_max, 1.386, places=3)


class TestDeployPresence(unittest.TestCase):
    """Cheap presence/version pins -- the Sec 13.4 drift class. Verify
    PRESENCE, not just version (job 2315894's lesson: version can be right
    while one import line is missing)."""

    def test_solver_version(self):
        self.assertEqual(ps.SOLVER_VERSION, EXPECT_VER,
                         f"SOLVER_VERSION drifted (got {ps.SOLVER_VERSION!r}, "
                         f"expected {EXPECT_VER!r}; set EXPECT_SOLVER_VERSION "
                         f"when bumping intentionally)")

    def test_residual_screen_functions_present(self):
        for fn in ("weak_enforcement_residual_oop",
                   "weak_enforcement_residual_ip",
                   "_weak_enforcement_null_block"):
            self.assertTrue(hasattr(det, fn),
                            f"detectors.{fn} missing -- residual screen not "
                            f"deployed (job 2315890 drift class)")

    def test_detectors_edgekind_import_line(self):
        # the A.21 NARROW drift: whole block present but this line missing
        self.assertIn("from .boundary import EdgeKind", _detectors_source(),
                      "detectors.py lost its EdgeKind import (job 2315894)")

    def test_gap_refine_defaults(self):
        src = _detectors_source().replace("'", '"')
        self.assertIn('"SIGMIN_GAP_REFINE", "1"', src,
                      "gap-refinement no longer default-on (s6 promotion "
                      "reverted?)")
        self.assertIn('"SIGMIN_GAP_REFINE_MAX_GAPS", "30"', src,
                      "gap cap default != 30 -- A.24-26: a cap of 5 silently "
                      "dropped the single most important gap in BOTH known "
                      "geometries; do not lower without re-running both gates")

    def test_shi_orthotropic_constants(self):
        m = ps.OrthotropicMaterial(E_r=40e9, E_theta=70e9, nu_r=0.3,
                                   G_rtheta=3.51e9, rho=7850)
        T, R, mu = (float(x) for x in m.oop_constants())
        self.assertAlmostEqual(T, 0.6729, delta=5e-3)
        self.assertAlmostEqual(R, 1.75, places=6)
        self.assertAlmostEqual(mu, 0.525, places=6)

    def test_free_free_bc_registry(self):
        for motion in ("out_of_plane", "in_plane"):
            bc = ps.make_bc("free_free", motion)
            self.assertEqual(sorted(e.sign for e in bc.edges), [-1, 1])
            self.assertTrue(all(e.kind == ps.EdgeKind.FREE for e in bc.edges))


def _nearest_err(x, targets):
    best = min(targets, key=lambda t: abs(t - x))
    return best, abs(best - x)


@unittest.skipUnless(GATE_FFP1, "set GATE_FFP1=1 (cluster, >=16 cores, "
                                "~1.5-3 h) -- the gate that caught the s7 "
                                "regression, job 2316094")
class TestFFP1TableGate(unittest.TestCase):
    """Fully-default detector on FF-P1 OOP must reproduce the 13-mode
    spectrum (the original Sec 10.8 11 modes, INCLUDING the gap-refine-
    recovered 0.961312, PLUS the 2026-07-11 addendum's two window-extension
    modes 2.320741/2.602500). If this fails, suspect
    SIGMIN_GAP_REFINE_MAX_GAPS truncation FIRST (Sec 16), not the detector
    logic."""

    def test_ffp1_oop_11_modes(self):
        mp.dps = int(os.environ.get("DPS", "40"))
        mat = ps.IsotropicMaterial(E=210e9, nu=0.30, rho=7800.0)
        geom = ps.make_geometry(1.5, 0.5)
        from plate_solver.validation import _scan_cfg_part1
        (lo, hi), ns, co, xmax = _scan_cfg_part1(1.5, 0.5, mat, n_modes_wanted=12)
        # UPDATED 2026-07-11 (addendum Sec C.1): window widened past the old
        # 2.45 ceiling so the two newly-confirmed modes (raw 2.320741,
        # 2.602500) are in scan range; spectrum now spans raw 0.32-2.60.
        hi = max(hi, 2.65)
        solver = ps.OutOfPlaneSolver(geom, mat, M=80, n_quad=30,
                                     boundary=ps.FreeFreeOOP())
        freqs = ps.find_modes_sigmin(
            solver, part=1, Omega_range=(lo, hi), n_scan=ns, n_dofs=20,
            max_dim=xmax, n_modes_wanted=16, verbose=True)
        self.assertGreaterEqual(len(freqs), 13,
                                f"only {len(freqs)} modes accepted: {freqs}")
        misses = []
        for ref in FFP1_OOP_RAW:
            _, d = _nearest_err(ref, freqs)
            if d >= 0.003:
                misses.append((ref, d))
        self.assertFalse(misses,
                         f"pinned FF-P1 modes not reproduced (ref, |d|): "
                         f"{misses}; computed={freqs}")


@unittest.skipUnless(GATE_FFP1_IP, "set GATE_FFP1_IP=1 (cluster, >=16 "
                                   "cores, ~5-7 h)")
class TestFFP1IPTableGate(unittest.TestCase):
    """Fully-default detector on FF-P1 IP (capture config: window extended
    to 2.85 with density-preserving n_scan) must reproduce all 43 pinned
    scan locations from job 2317336 -- artifacts included, since this pins
    the n20 landscape's stability, not physics. Verdicts are NOT asserted.
    If this fails, suspect SIGMIN_GAP_REFINE_MAX_GAPS truncation FIRST
    (Sec 16), then window/n_scan config drift (the density-preserving
    extension below must match probe_ip_baseline_capture.py exactly)."""

    def test_ffp1_ip_scan_table(self):
        mp.dps = int(os.environ.get("DPS", "40"))
        mat = ps.IsotropicMaterial(E=210e9, nu=0.30, rho=7800.0)
        geom = ps.make_geometry(1.5, 0.5)
        from plate_solver.validation import _scan_cfg_part2
        (lo, hi), ns, co, zmax = _scan_cfg_part2(1.5, 0.5, mat,
                                                 n_modes_wanted=22)
        hi_ext = max(hi, 2.85)
        ns_ext = max(ns, int(round(ns * (hi_ext - lo) / (hi - lo))))
        solver = ps.InPlaneSolver(geom, mat, M=80, n_quad=30,
                                  boundary=ps.FreeFreeIP())
        freqs = ps.find_modes_sigmin(
            solver, part=2, Omega_range=(lo, hi_ext), n_scan=ns_ext,
            n_dofs=20, max_dim=zmax, n_modes_wanted=40, verbose=True)
        self.assertGreaterEqual(len(freqs), 43,
                                f"only {len(freqs)} dips accepted: {freqs}")
        misses = []
        for ref in FFP1_IP_RAW_SCAN:
            _, d = _nearest_err(ref, freqs)
            if d >= 0.003:
                misses.append((ref, d))
        self.assertFalse(misses,
                         f"pinned FF-P1 IP scan locations not reproduced "
                         f"(ref, |d|): {misses}; computed={freqs}")


@unittest.skipUnless(GATE_SHI, "set GATE_SHI=1 (cluster, >=16 cores, "
                               "~3-5 h)")
class TestShiTableGate(unittest.TestCase):
    """Fully-default Shi orthotropic FFFF discovery must contain the final
    8/8 table's computed values (UPDATED 2026-07-11: all 8 now REAL, the
    former stably-AMBIGUOUS 66.76 candidate resolved REAL via job 2316744 --
    addendum Sec B) to <1.0% nearest-neighbour in Omega_lit. Verdicts are NOT
    asserted here -- this pins mode EXISTENCE/location against silent
    detector regressions."""

    def test_shi_ortho_table(self):
        mp.dps = int(os.environ.get("DPS", "40"))
        omat = ps.OrthotropicMaterial(E_r=40e9, E_theta=70e9, nu_r=0.3,
                                      G_rtheta=3.51e9, rho=7850)
        geom = ps.make_geometry(1.5, 0.5)
        from plate_solver.validation import _scan_cfg_part1
        from plate_solver.geometry import _omega_lit, _kbar_wbar
        (lo, hi), ns, co, xmax = _scan_cfg_part1(1.5, 0.5, omat, n_modes_wanted=14)
        hi = max(hi, 2.45)   # highest pinned lit 96.03 -> raw ~2.43
        solver = ps.OutOfPlaneSolver(geom, omat, M=80, n_quad=30,
                                     boundary=ps.FreeFreeOOP())
        freqs = ps.find_modes_sigmin(
            solver, part=1, Omega_range=(lo, hi), n_scan=ns, n_dofs=20,
            max_dim=xmax, n_modes_wanted=18, verbose=True)
        k_bar, w_bar = _kbar_wbar(geom, omat)
        lits = [_omega_lit(f, geom, omat, w_bar, k_bar, 1)[1] for f in freqs]
        misses = []
        for ref in SHI_REAL_LIT:
            best, d = _nearest_err(ref, lits)
            if 100.0 * d / ref >= 1.0:
                misses.append((ref, best))
        self.assertFalse(misses,
                         f"pinned Shi table values not reproduced <1% "
                         f"(ref, nearest): {misses}; computed lits={lits}")


if __name__ == "__main__":
    unittest.main()
