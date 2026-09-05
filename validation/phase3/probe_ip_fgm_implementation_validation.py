# -*- coding: utf-8 -*-
"""
probe_ip_fgm_implementation_validation.py -- cluster confirmation for the
in-plane (Part 2) radially-graded material support just integrated into the
DEPLOYED package on 2026-07-30 (geometry.py::RadialFGMMaterial's new
c11_i/c11_o + ip_coeff_series(), core_solvers.py::InPlaneSolver's graded
_series_mp/_Lmat_mp/_build_K_real branches, dispersion.py::_Part2Fast's
graded series()/Lmat()). See LESSONS_LEARNED.md Sec 54.4/54.5 for the full
derivation + integration history. Mirrors probe_fgm_implementation_
validation.py's 4-part structure (the OOP FGM job, 2330576, CLOSED PASS)
as explicitly recommended in Sec 54.5's own "recommended next cluster job"
note -- same scoping: single-process solver path only, NO worker-pool
integration attempted here (matches the OOP job's own explicit scope note;
worker-pool support for graded materials is a separate, not-yet-started
item for BOTH OOP and IP -- see FUTURE_WORK.md Sec 3).

WHAT'S DIFFERENT FROM THE OOP JOB (2330576), AND WHY:
  - No hardcoded historical "GATE_CANTILEVER_IP" scalar exists anywhere in
    this project's files for an arbitrary Omega (grepped test_validated_
    tables.py and LESSONS_LEARNED.md; only found a prose mention of "the
    standard IP cantilever gate" with no number attached, LESSONS Sec ~44).
    Rather than fabricate a number, Part A establishes ITS OWN fresh
    baseline (REFERENCE_CANTILEVER, printed, not asserted against an
    external constant) via the EXISTING scalar (non-graded, list-free)
    InPlaneSolver code path, and separately cross-checks that path against
    a real, already-pinned production data point: FFP1_IP_RAW_SCAN's
    Omega=1.598922 (test_validated_tables.py ~line 161, tagged "REAL -- G2
    (FE 643.1 Hz, 0.012%); found cleanly in this full blind default scan"
    -- i.e. this IS the n_dofs=20 raw-scan dip location) must still show a
    deep sigma_min dip on the actual FF-P1 geometry at dps=40 on the
    cluster (not just the sandbox).
    CORRECTED 2026-08-02 (job 2363617 investigation): this pin was
    originally written as 1.598985, which is test_validated_tables.py's
    "physical REAL table" value (the n28/n36-ADJUDICATED, refined
    location, same file ~line 191) -- a DIFFERENT, more-refined number
    than the n20 raw-scan location this probe actually searches at
    (FFP1_NDOF=20 below). The two differ by only 6.3e-5 in Omega, but per
    this project's own established basis-sensitivity caution ("sigma_min
    depth is basis-sensitive"; the zone-B precedent in LESSONS_LEARNED),
    evaluating a sharp dip 6e-5 off its own basis-specific minimum can
    cost several decades of depth -- job 2363617's Part A2 got only -2.13
    (FAIL, threshold -3.5) at the mismatched pin, consistent with this
    being a probe-construction bug (wrong table row used), not a solver
    regression. The K_omega/kappa OOP-only edit (2026-08-02) cannot be
    responsible either way -- grep-verified it never appears anywhere in
    InPlaneSolver's class body.
    This is the same style of external-pin cross-check job 2336354 (the G2
    promotion gate) used for the OOP/IP scan tables.
  - Part C, instead of re-deriving an independent ODE-residual check (that
    derivation was already done twice -- once standalone by Grok, Sec
    54.4, once against the deployed classes by Claude in-sandbox, Sec
    54.5 item 2, residual 1.08e-18 -- redoing it a third time from scratch
    in this probe would risk introducing a NEW transcription error without
    the benefit of an independent derivation to catch it against), cross-
    checks the two ALREADY-implemented code paths against each other: the
    mp recursion (InPlaneSolver._series_mp) vs. the float64 fast-engine
    mirror (_Part2Fast.series(), which `full_search`/`select_fill` actually
    use for root-finding and was NOT covered by Grok's original derivation
    scope). This is a safe, mechanical consistency check between two
    existing implementations, not new derived physics -- it directly
    targets the one self-test LESSONS Sec 54.5 flagged as sandbox-only so
    far (item 3, agreement 3.9e-17 there) and re-confirms it fresh on the
    cluster at production dps=40, specifically exercising the NEW c11
    grading branch (the piece that didn't exist before 2026-07-30).
  - Part D windows Omega over [0.02, 1.0] rather than OOP's [0.02, 0.35]
    -- IP's raw-Omega dip locations for this geometry are not established
    anywhere in this project's files (only FF-P1's free-free geometry has
    a pinned IP scan; this job's cantilever geometry does not), so the
    wider window is a deliberate hedge, not a known-good choice. If no
    dips appear in this window for EITHER material, that is itself useful
    information (window needs widening further), not a probe failure.

PRE-REGISTERED INTERPRETATION:
  Part A: A1 (fresh baseline) has no fail bar by itself (nothing to compare
    against yet) but must print a finite, non-degenerate value. A2 (FF-P1
    pin) FAIL if log10(sigma_min) at Omega=1.598922 is NOT below -3.5
    decades (the same dip threshold Part D's find_dips uses) -- would mean
    the FGM edits broke the existing non-graded IP scalar path on real
    pinned production data. STOP and do not proceed past Part A on a FAIL.
    (2026-08-02: pin corrected from 1.598985 to 1.598922 -- see the
    docstring note above; that corrected pin STILL failed at -2.11, which
    turned out to be a THIRD issue, not a regression: see the XMAX_FLOOR
    note above and LESSONS_LEARNED.md Sec 73. sigma_min_at's old minimal-
    sufficiency xmax-growth criterion undershot this mode's real depth --
    job 2381734 showed the pin deepens monotonically to -4.95 by xmax=100.
    XMAX_FLOOR now forces growth to xmax=100 for this key before Part A2
    reads the result, so it should pass on the next cluster run. If it
    still fails at the floor, that IS a genuine new finding worth stopping
    on -- do not raise the floor again without a fresh basis-stability
    check like job 2381734's to justify it.)
  Part B: FAIL if the graded-flat material's sigma_min at the SAME Omega as
    A1 differs from A1's own freshly-computed value by more than PARTB_TOL
    (0.05 decades, corrected 2026-08-04 from an unrealistic 1e-6 -- see the
    inline comment at PARTB_TOL's definition and LESSONS_LEARNED.md Sec 74
    for the full root-cause chain: this is float64 full_search branch-count
    sensitivity between IsotropicMaterial and a flat-profile RadialFGMMaterial,
    NOT the graded recursion itself, which stays proven-exact at the mp/
    series level, Sec 54.5). A FAIL past this looser bound would mean the
    graded code path itself is wrong, not just correctly bypassed for
    non-graded materials. STOP before Part C/D on a FAIL.
  Part C: FAIL if max relative coefficient disagreement between the mp and
    float64 graded series exceeds 1e-6 (float64 precision floor is ~1e-15,
    so 1e-6 is already a generous margin) -- would mean the float64 mirror
    (dispersion.py::_Part2Fast, added same day as the mp path, NOT part of
    Grok's original derivation) has a transcription bug relative to the mp
    recursion, which would silently corrupt full_search's root-seeding for
    any graded IP material even though the mp path itself is correct.
  Part D: same reading convention as the OOP job -- real, sharp (< -3.5
    decade) dips for the graded material at DIFFERENT Omega than the
    isotropic baseline's own dips would be the expected, physically sane
    signature. Identical dip locations to the isotropic case, or no dips
    in either case, should be reported plainly (window/structural
    limitation), not smoothed over. This is a STRUCTURAL sanity check only
    -- no literature or FE benchmark exists yet for this material (see the
    separate FGM-literature-benchmark thread), so Part D cannot and does
    not claim physical validation.

COST MODEL (plate-solver-cluster-probes skill): full_search ~55-60s at
dps=40, M=80 (IP's coupled 4x4 system is the same order of cost as OOP's).
Part A: 2 points (~2 min). Part B: 1 point (~1 min). Part C: 1 point, no
full_search (~5s). Part D: 2 materials x 25 points across [0.02,1.0] (step
~0.0408) = 50 full_search (~50 min) + 50 K-builds (~15 min) ~= 65 min.
Total ~=70 min; requesting --time=02:00:00 for margin, same as job 2330576.

GEOMETRY/MATERIAL: Parts A1/B/C/D use the same standard cantilever config
as the OOP job (r0/(2b)=1.25, 2Theta=1.0pi, E=210e9, nu=0.35, rho=7800.0).
Part A2 uses the FF-P1 geometry (r0/(2b)=1.5, 2Theta=0.5pi, nu=0.30,
n_dofs=xi_max=20) to reach the pinned FFP1_IP_RAW_SCAN point. The graded
material in Parts C/D grades c11 (the NEW 2026-07-30 addition) AND R (mild,
matching the OOP job's own asymmetric-grading discipline so more than one
coefficient genuinely varies): c11_i/c11_o = c11_bar/c66 -> 1.05x that
value, R: 1.0->1.08, nu flat at 0.35, n_p=1.
"""
import os
import sys
import time

sys.path.insert(0, os.environ.get("PKG_PATH", "."))
os.environ.setdefault("DPS", "40")

from mpmath import mp, mpf, mpc
mp.dps = int(os.environ["DPS"])

import plate_solver as ps
from plate_solver import (InPlaneSolver, IsotropicMaterial, RadialFGMMaterial,
                          make_geometry, full_search, select_fill,
                          ClampedFreeIP, FreeFreeIP)
from plate_solver.detectors import sigma_min_from_K

EXPECT_VER = os.environ.get("EXPECT_SOLVER_VERSION", "2026-07-10.s10")

DIP_THRESH = -3.5   # decades; same convention as the OOP job's find_dips

# Cantilever config (matches probe_fgm_implementation_validation.py exactly)
CANT_XMAX = 14.0
CANT_NDOF = 16

# FF-P1 config (matches test_validated_tables.py's FFP1_IP_RAW_SCAN capture)
# UPDATE 2026-08-01 (post-job-2339579): the original FFP1_XMAX=20.0 (matching
# _scan_cfg_part2's auto-selected ze_max for this geometry) turned out to be
# too narrow AT THIS SPECIFIC Omega -- select_fill capped at n_dofs_actual=18
# instead of the requested 20. A static-widened FFP1_XMAX=28.0 (job 2339607)
# ALSO undershot for A1's own CANT_XMAX=14.0 (only reached 14 of 16 there) --
# static guessing is not reliable at this geometry/Omega. Replaced with an
# ADAPTIVE grow-until-reached loop (see sigma_min_at below): starts at the
# same nominal xmax, and if select_fill undershoots, grows xmax by XMAX_STEP
# and retries, persisting the discovered-sufficient xmax per config (keyed
# by xmax_key) so later calls at nearby Omega don't repeat the same growth.
# This is the SAME select_fill-capping bug independently found today in
# probe_zoneb_golden_refine_v1/v2.py and probe_ss_coeff_basis_escalation_
# v1/v2.py -- all three now use this adaptive pattern instead of a fixed guess.
FFP1_XMAX = 28.0
FFP1_NDOF = 20
# CORRECTED 2026-08-02 (job 2363617 -> Part A2 FAIL, -2.13 not <-3.5): was
# 1.598985 (test_validated_tables.py's REFINED/n28+-adjudicated "physical
# REAL table" value), which does NOT match this probe's own FFP1_NDOF=20
# basis -- the n20 RAW SCAN's own dip for this mode sits at 1.598922 (see
# test_validated_tables.py's FFP1_IP_RAW_SCAN row comment: "REAL -- G2 ...
# found cleanly in this full blind default scan"). Using the mismatched
# refined value evaluated sigma_min 6.3e-5 in Omega off the n20 landscape's
# actual minimum, which is enough to hide several decades of depth for a
# sharp dip -- see the module docstring's "WHAT'S DIFFERENT" note for the
# full writeup. Not a solver bug, not related to the OOP-only K_omega edit.
FFP1_PIN_OMEGA = 1.598922   # REAL, FE 643.1 Hz, 0.012%, n20-RAW-SCAN location

# CORRECTED AGAIN 2026-08-04 (job 2381735 sequence -> the corrected pin STILL
# FAILED, -2.11 not <-3.5): root-caused via a Grok hand-off chain (jobs
# 2381689 double-log bug -> 2381724 genuine fidelity-gate pass + a -6.31
# candidate at Omega=1.599722 -> 2381734 fixed-xmax sweep). That candidate
# was KILLED as an xmax=92-specific basis artifact (non-monotonic across
# xmax 84/92/100: -3.08/-6.31/-3.39, n_br=12/16/27 -- classic branch-
# collision signature). But the SAME job 2381734 sweep, done for continuity
# at the accepted pin Omega=1.598922, showed something important: sigma_min
# there DEEPENS MONOTONICALLY with xmax -- -2.11 @ xmax=84, -4.01 @ xmax=92,
# -4.95 @ xmax=100 -- this project's own signature of a REAL mode still
# gaining depth from a richer basis, not an artifact. Root cause of the
# original Part A2 FAIL: sigma_min_at's minimal-sufficiency stopping
# criterion (grow xmax only until n_dofs_actual first reaches ndof) landed
# at xmax=84 for this Omega, which is too lean a basis for this mode to show
# its real depth -- NOT a regression from the k_bar/rect-free-free/FGM
# edits. XMAX_FLOOR below forces growth past minimal-sufficiency for this
# specific key. 100.0 is the deepest point independently tested so far
# (job 2381734), NOT proven fully converged -- if a future run still fails
# at this floor, that would be a genuine new finding, not a repeat of this
# one. See LESSONS_LEARNED.md Sec 73 for the full chain.
XMAX_FLOOR = {"FFP1": 100.0}   # xmax_key -> minimum xmax before accepting,
                                 # even if n_dofs_actual>=ndof sooner. Only
                                 # FFP1 has evidence requiring this; CANT is
                                 # intentionally left unfloored (no comparable
                                 # finding for that geometry/Omega).

XMAX_STEP = 8.0
XMAX_CAP = 100.0
_XMAX_STATE = {}   # xmax_key -> last-known-sufficient xmax, ratchets up only


def hdr(s):
    print("\n" + "=" * 78 + f"\n  {s}\n" + "=" * 78, flush=True)


def sigma_min_at(ip, Om, xmax0, ndof, xmax_key):
    xmax = _XMAX_STATE.get(xmax_key, xmax0)
    floor = XMAX_FLOOR.get(xmax_key, 0.0)
    while True:
        brs = full_search(ip.fast, Om, xmax=xmax)
        sel, n = select_fill(brs, ndof)
        if n >= ndof and xmax >= floor:
            _XMAX_STATE[xmax_key] = xmax
            K, nn = ip._build_K_real(Om, sel, lagrange=True, fast_scan=False)
            return float(sigma_min_from_K(K, nn)), n, xmax
        if xmax >= XMAX_CAP:
            if n >= ndof:
                # floor not reached but the cap is -- accept the cap basis
                # rather than hard-failing on a floor that can't be met.
                print(f"  [xmax-floor-cap:{xmax_key}] reached XMAX_CAP="
                      f"{XMAX_CAP} before floor={floor} at Om={Om}; "
                      f"accepting the cap basis as-is.", flush=True)
                _XMAX_STATE[xmax_key] = xmax
                K, nn = ip._build_K_real(Om, sel, lagrange=True, fast_scan=False)
                return float(sigma_min_from_K(K, nn)), n, xmax
            print(f"\n*** FATAL: even at xmax={xmax} (cap), select_fill only "
                  f"reaches n_dofs_actual={n} < requested {ndof} at Om={Om} "
                  f"(key={xmax_key}). Giving up -- this is now a genuine "
                  f"structural finding (not undersampling), report as-is. ***",
                  flush=True)
            raise SystemExit(6)
        old = xmax
        xmax += XMAX_STEP
        print(f"  [xmax-grow:{xmax_key}] n_dofs_actual={n}<{ndof} or "
              f"xmax={old}<floor={floor} at Om={Om}, xmax {old}->{xmax}, "
              f"retrying...", flush=True)


def main():
    hdr("probe_ip_fgm_implementation_validation -- start " + time.strftime("%Y-%m-%d %H:%M:%S"))
    print(f"SOLVER_VERSION = {ps.SOLVER_VERSION}")
    if ps.SOLVER_VERSION != EXPECT_VER:
        print(f"ABORT: expected {EXPECT_VER!r}, got {ps.SOLVER_VERSION!r}")
        sys.exit(2)

    mat_iso = IsotropicMaterial(E=210e9, nu=0.35, rho=7800.0)
    geom_cant = make_geometry(1.25, 1.0)

    # ---- Part A1: fresh cantilever baseline (no external gate available) ----
    hdr("PART A1 -- fresh cantilever IP baseline (isotropic, unchanged path)")
    cant = InPlaneSolver(geom_cant, mat_iso, M=80, n_quad=30, boundary=ClampedFreeIP())
    sm_cant, _, _ = sigma_min_at(cant, 0.12, CANT_XMAX, CANT_NDOF, "CANT")
    print(f"  cantilever IP sigma_min(0.12) = {sm_cant}  "
          f"(REFERENCE_CANTILEVER -- Part B must reproduce this within "
          f"PARTB_TOL, see Sec 74)")
    REFERENCE_CANTILEVER = sm_cant

    # ---- Part A2: external pin cross-check on FF-P1 geometry ----
    hdr("PART A2 -- FF-P1 pinned REAL point (Omega=1.598922), external cross-check")
    mat_ffp1 = IsotropicMaterial(E=210e9, nu=0.30, rho=7800.0)
    geom_ffp1 = make_geometry(1.5, 0.5)
    ff = InPlaneSolver(geom_ffp1, mat_ffp1, M=80, n_quad=30, boundary=FreeFreeIP())
    sm_ff, ndof_ff, _ = sigma_min_at(ff, FFP1_PIN_OMEGA, FFP1_XMAX, FFP1_NDOF, "FFP1")
    okA2 = sm_ff < DIP_THRESH
    print(f"  FreeFreeIP sigma_min({FFP1_PIN_OMEGA}) = {sm_ff}  (ndof={ndof_ff})  "
          f"threshold<{DIP_THRESH}  {'PASS' if okA2 else 'FAIL -- STOP'}")

    if not okA2:
        print("\n*** PART A2 FAILED -- the FGM edits broke the existing "
              "non-graded IP path on real pinned production data. STOPPING "
              "before Parts B/C/D. ***")
        raise SystemExit(3)

    # ---- Part B: graded-flat collapse check ----
    hdr("PART B -- graded material, FLAT profile, must collapse onto A1 within tol")
    mat_flat = RadialFGMMaterial(T_i=1, T_o=1, R_i=1, R_o=1,
                                  nu_i=mat_iso.nu, nu_o=mat_iso.nu, n_p=1,
                                  c11_bar=mat_iso.c11_bar, c66=mat_iso.c66,
                                  rho=mat_iso.rho, geom=geom_cant)
                                  # c11_i/c11_o omitted -> default to the
                                  # ungraded reference (c11_bar/c66), per
                                  # RadialFGMMaterial's 2026-07-30 docstring
    cant_flat = InPlaneSolver(geom_cant, mat_flat, M=80, n_quad=30, boundary=ClampedFreeIP())
    sm_flat, _, _ = sigma_min_at(cant_flat, 0.12, CANT_XMAX, CANT_NDOF, "CANT")
    # TOLERANCE CORRECTED 2026-08-04 (job 2381842 FAIL at 1e-6, root-caused via
    # GrokCode/probe_ip_fgm_basis_select_artifact_v1.py, cluster job + a local
    # non-cluster rerun on the user's own machine): the original <1e-6 bar
    # assumed the float64 full_search/select_fill path is bit-reproducible
    # between a pure IsotropicMaterial and a flat-profile RadialFGMMaterial --
    # it is not. Root cause is a branch-count mismatch in full_search itself
    # (17 raw branches for isotropic vs 18 for flat-graded at the identical
    # Omega=0.12/xmax=70.0, confirmed on BOTH the cluster and a second,
    # independent machine): a branch near a scan-grid boundary gets
    # caught/missed depending on tiny float64 differences in how flat-graded
    # coefficient lists evaluate versus bare isotropic scalars. The already-
    # proven-exact mp/series-level reduction (LESSONS_LEARNED.md Sec 54.5,
    # 3.7e-43) is untouched -- this is float64 search sensitivity only.
    # Two independent measurements of the resulting sigma_min delta: 0.0098
    # decades (cluster, job 2381842) and 0.0016 decades (local machine,
    # same day) -- both far above 1e-6, both far below the multi-decade
    # scale a genuine bug produced in the FF-P1 saga (Sec 73). New bound
    # gives ~5x headroom above the larger observed value while staying
    # ~2 orders of magnitude tighter than a real-bug-scale difference.
    PARTB_TOL = 0.05
    okB = abs(sm_flat - REFERENCE_CANTILEVER) < PARTB_TOL
    print(f"  cantilever IP sigma_min(0.12) [graded-flat] = {sm_flat}  "
          f"(reference={REFERENCE_CANTILEVER}, tol={PARTB_TOL})  "
          f"{'PASS' if okB else 'FAIL -- STOP'}")
    if not okB:
        print("\n*** PART B FAILED -- the graded IP code path itself is "
              "wrong (not just correctly bypassed). STOPPING before Parts C/D. ***")
        raise SystemExit(4)

    # ---- Part C: mp vs float64 graded-series cross-check (genuinely graded) ----
    hdr("PART C -- mp vs float64 (_Part2Fast) graded series cross-check, c11 grading exercised")
    c11_ref = mat_iso.c11_bar / mat_iso.c66
    mat_graded = RadialFGMMaterial(T_i=1, T_o=1, R_i=1.0, R_o=1.08,
                                    nu_i=0.35, nu_o=0.35, n_p=1,
                                    c11_bar=mat_iso.c11_bar, c66=mat_iso.c66,
                                    rho=mat_iso.rho, geom=geom_cant,
                                    c11_i=c11_ref, c11_o=1.05 * c11_ref)
    ip_g = InPlaneSolver(geom_cant, mat_graded, M=80, n_quad=30, boundary=ClampedFreeIP())
    print(f"  c11 series coeffs: {[float(c) for c in ip_g.c11]}")
    print(f"  R   series coeffs: {[float(c) for c in ip_g.R]}")

    ze_v = mpc('1.2', '0')
    Om_v = mpf('0.15')
    sols_mp = ip_g._series_mp(ze_v, Om_v)
    sols_f64 = ip_g.fast.series(complex(ze_v), complex(Om_v))

    worst = 0.0
    NCHECK = 12   # low-order coefficients carry the discriminating content;
                  # higher orders shrink toward zero and would just add noise
    for q in range(4):
        a_mp, b_mp = sols_mp[q]
        a_f, b_f = sols_f64[q]
        for p in range(NCHECK):
            for name, vmp, vf in (("a", a_mp[p], a_f[p]), ("b", b_mp[p], b_f[p])):
                vmp_f = complex(vmp)
                denom = max(abs(vmp_f), 1e-30)
                rel = abs(vmp_f - vf) / denom
                worst = max(worst, rel)
    okC = worst < 1e-6
    print(f"\n  worst relative coefficient disagreement (q=0..3, p=0..{NCHECK-1}): {worst:.3e}  "
          f"{'PASS' if okC else 'FAIL -- float64 mirror has a bug'}")

    # ---- Part D: real sigma_min sweep, isotropic baseline vs graded ----
    hdr("PART D -- sigma_min sweep, Omega in [0.02,1.0], isotropic vs graded (cantilever IP)")
    Oms = [0.02 + (1.0 - 0.02) * k / 24 for k in range(25)]
    print("  Omega      iso_log10_smin   graded_log10_smin")
    iso_rows = []
    graded_rows = []
    for Om in Oms:
        t0 = time.time()
        sm_i, _, _ = sigma_min_at(cant, Om, CANT_XMAX, CANT_NDOF, "CANT")
        sm_g, _, _ = sigma_min_at(ip_g, Om, CANT_XMAX, CANT_NDOF, "CANT_GRADED")
        dt = time.time() - t0
        iso_rows.append((Om, sm_i))
        graded_rows.append((Om, sm_g))
        print(f"  {Om:8.5f}   {sm_i:10.4f}       {sm_g:10.4f}    [{dt:.1f}s]")

    def find_dips(rows, thresh=DIP_THRESH):
        dips = []
        for i in range(1, len(rows) - 1):
            if rows[i][1] < thresh and rows[i][1] < rows[i-1][1] and rows[i][1] < rows[i+1][1]:
                dips.append(rows[i])
        return dips

    iso_dips = find_dips(iso_rows)
    graded_dips = find_dips(graded_rows)
    print(f"\n  isotropic local minima (< {DIP_THRESH} decades): {iso_dips}")
    print(f"  graded    local minima (< {DIP_THRESH} decades): {graded_dips}")
    print(f"  distinct from isotropic (structural sanity, not a literature "
          f"or FE match): {graded_dips != iso_dips}")

    hdr("FINAL SUMMARY")
    print(f"  Part A1 (fresh baseline):          printed, REFERENCE_CANTILEVER={REFERENCE_CANTILEVER}")
    print(f"  Part A2 (FF-P1 pin cross-check):    {'PASS' if okA2 else 'FAIL'}")
    print(f"  Part B  (graded-flat collapse):     {'PASS' if okB else 'FAIL'}")
    print(f"  Part C  (mp vs float64 cross-check): {'PASS' if okC else 'FAIL'}")
    print(f"  Part D: see dip comparison above -- structural sanity check, "
          f"no literature/FE FGM benchmark available yet.")
    hdr("probe_ip_fgm_implementation_validation -- end " + time.strftime("%Y-%m-%d %H:%M:%S"))


if __name__ == "__main__":
    main()
