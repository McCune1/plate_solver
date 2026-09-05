# -*- coding: utf-8 -*-
"""
probe_ff_residual_screen_geomsweep.py -- DIAGNOSTIC/DATA-GENERATION ONLY.
Writes only a plain-text log; no package changes, no SOLVER_VERSION
implications.

CONTEXT: the 2026-07-22 python-vs-independent-FE mode-matching analysis
across the 48 new (r0_2b x 2Theta/pi x part) free-free geometries
(FutureWork/geometry_sweep/, FutureWork/ansys_geomsweep_verification/)
found that, beyond each geometry's first few accepted candidates, 84.0%
of the remaining candidates matched an independent FE mode within 3%
(mean 0.58%, worst 3.00%) -- but that same analysis had to SET ASIDE each
geometry's first min(n_cutoffs, n_accepted) candidates (5 for OOP/part 1,
4 for IP/part 2, per that geometry's own printed cut-off count) because
they carried a much higher mismatch rate (62% unmatched, mean error
23.8%) that tracks the cut-off count exactly -- the same near-cut-off
conditioning degradation this project's residual/parity screen
(detectors.weak_enforcement_residual_oop/_ip, "Section 5b") already
exists to adjudicate at the ONE geometry (FF-P1, r0/2b=1.5, 2Theta=0.5pi,
nu=0.30) where it was originally validated (11/11 correct against Ansys,
LESSONS_LEARNED Sec 10.7/10.8; PAPER1_FREEFREE_DRAFT.tex tab:oop-spurious).

WHAT THIS DOES: for the 215 cut-off-adjacent candidates set aside by that
matching analysis (TARGETS below, one row per (r0_2b, two_T_pi, part,
Omega) -- extracted directly from FutureWork/geometry_sweep/figures_ff_r*/
freefree_manifest.json's own 'diag.accepted' lists, i.e. already-computed,
already-accepted determinant zeros, NOT a new discovery search), this
probe:
  1. Reconstructs the branch set at that Omega via full_search + select_fill
     using the SAME n_dofs and xmax the original geometry-sweep run itself
     used for that (r0_2b, two_T_pi, part) -- read directly from each
     geometry's own manifest entry (diag.n_dofs, diag.max_dim), not
     re-derived or guessed. This is a fixed reconstruction at an
     already-validated root, the same technique probe_ff_oop_spurious_
     table.py already used for the single-geometry supplementary table.
  2. Runs detectors.weak_enforcement_residual_oop (part=1) or
     weak_enforcement_residual_ip (part=2) to get the pointwise edge-
     residual verdict (REAL-like / ARTIFACT-like / AMBIGUOUS for OOP;
     raw rTyy/rTyr/maxDisp for IP, which has no numeric verdict bar yet
     -- see the caveat below).
  3. Cross-checks solver.sigma_min at the same (Omega, sel) as a cheap
     independent depth sanity check (reuses sel, no extra full_search).

HONEST SCOPE CAVEAT -- read before trusting any verdict this probe
prints: `weak_enforcement_residual_oop`'s acceptance thresholds
(real_rM_max=0.1, artifact_rM_min=0.3, etc.) were calibrated ONCE, at ONE
geometry (FF-P1, r0/2b=1.5, 2Theta=0.5pi) and ONE material (nu=0.30).
This probe applies those SAME thresholds, unchanged (per this project's
own "resist tuning, don't move the goalposts" methodology), to 4 NEW
radius ratios, 6 sector angles, and a DIFFERENT material (nu=0.35) --
this is a genuine extrapolation of the thresholds' validated domain, not
just of the search. The physics underlying the screen (a pointwise
Kelvin shear/moment residual on the free edge) has no reason to be
geometry- or material-specific in its qualitative REAL-vs-ARTIFACT
separation (the FF-P1 calibration found a ~2000x gap between the two
populations, not a marginal one), but this is the FIRST time it has
been tested outside that one calibration point. Report the raw
(rV, rM/rTyy, maxW/maxDisp) numbers alongside every verdict so a human
can independently judge whether the same order-of-magnitude separation
holds, rather than trusting the verdict label alone. The in-plane
(part=2) screen has NO validated verdict bar at all yet (module
docstring: "NO numeric calibration exists yet for IP" even at FF-P1) --
its rows are reported as raw ratios only, exactly as
weak_enforcement_residual_ip's own contract promises, not force-labeled.

PRE-REGISTERED INTERPRETATION (write this BEFORE the job runs, honor it
after):
  - A candidate whose OOP verdict is REAL-like, with rM comfortably below
    ~0.1 (ideally in the same 0.002-0.04 range the FF-P1 calibration
    found for its 8 confirmed-real modes) -> genuinely supports promoting
    that specific candidate from "cut-off-adjacent, unadjudicated" to
    "confirmed real" in a future paper table, on the same standard as
    tab:oop/tab:oop-spurious.
  - A candidate whose OOP verdict is ARTIFACT-like, with rM landing in
    the FF-P1 calibration's ~75-350 artifact band -> supports the
    candidate being a weak-enforcement artifact (a spurious determinant
    zero, not a physical mode) at this geometry too, i.e. NOT evidence
    against the detector -- exactly analogous to tab:oop-spurious's 5
    rows at FF-P1.
  - A candidate returning AMBIGUOUS, or a REAL-like/ARTIFACT-like verdict
    whose raw rM sits well outside either FF-P1 calibration band (e.g. an
    rM of 1-10, between the ~0.04 real ceiling and the ~75 artifact
    floor) -> do NOT trust the verdict label at face value; this would be
    the first sign that the FF-P1 thresholds don't transfer cleanly to
    the new geometries/material, and the raw numbers (not the label)
    should go into any writeup, with an explicit note that re-calibration
    may be needed.
  - A residual-screen or branch-reconstruction FAILURE (select_fill short
    of n_dofs, "no null vector", dimension mismatch) -> report plainly;
    do not substitute a neighboring Omega or drop the row silently.
  - IP (part=2) rows: report raw rTyy/rTyr/maxDisp only, no verdict.
    Empirically (per the module docstring) rTyy<0.3 vs >4 has separated
    cleanly on every point checked so far at FF-P1, but this is
    exploratory there too -- doubly so here.

COST: 215 candidates total, split across 4 jobs (one per radius ratio,
matching every other stage of this thread's own r150/r167/r200/r250
convention): r150 has 53, r167/r200/r250 have 54 each. Per-point cost
mirrors probe_ff_oop_spurious_table.py's own estimate (1 full_search,
dps=40, xmax=14-20 depending on geometry, ~55-60s + 1 residual-screen
build, comparable cost to a full K build with per-item series/amplitude
reconstruction, ~60-90s + 1 cheap sigma_min reuse, ~1-2s) ~= ~150s/point
serial. Parallelized 1 point/worker via ProcessPoolExecutor, sized to
the SLURM allocation (RATIO_TAG selects the ~53-54-point subset for that
job; request cpus-per-task >= that count for full parallelism) --
wall time should be a few minutes per job plus queue overhead.

Env var RATIO_TAG (r150/r167/r200/r250) selects which ratio's targets
this run processes -- set by the matching submit_ff_residual_screen_
r*.sh script, not passed on the command line.
"""
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

sys.path.insert(0, os.environ.get("PKG_PATH", "."))
os.environ.setdefault("DPS", "40")

EXPECT_VER = os.environ.get("EXPECT_SOLVER_VERSION", "2026-07-10.s8")
RATIO_TAG = os.environ.get("RATIO_TAG")
if not RATIO_TAG:
    print("FATAL: RATIO_TAG env var not set (expected one of r150/r167/"
          "r200/r250). Set it in the submit script.")
    raise SystemExit(2)

NU_SWEEP = 0.35   # matches geometry_sweep's actual material (see
                   # LESSONS_LEARNED.md Sec 22 item 3 -- IsotropicMaterial
                   # default, confirmed from every sweep run's own
                   # "material this run: nu=0.35" printout), NOT the
                   # FF-P1 tab:oop/tab:ip benchmark's nu=0.30.

# ----------------------------------------------------------------------------
# Auto-extracted 2026-07-22 from FutureWork/geometry_sweep/figures_ff_r*/
# freefree_manifest.json's own 'diag.accepted' lists (already-computed,
# already-accepted candidates from the completed geometry sweep jobs --
# NOT a new discovery search). Each row is one cut-off-adjacent
# candidate: the first min(n_cutoffs, n_accepted) accepted candidates
# per (r0_2b, two_T_pi, part), sorted ascending by Omega, where
# n_cutoffs = len(diag.cut_offs) for that geometry (5 for part=1/OOP,
# 4 for part=2/IP, confirmed uniform across all 24 (ratio,angle)
# combinations this session). These are exactly the candidates the
# python-vs-FE matching analysis (2026-07-22) set aside as needing this
# residual/parity screen before being trusted -- NOT the full 1466-mode
# set (which already has FE-match evidence and doesn't need re-screening).
# Row: (tag, r0_2b, two_T_pi, part, n_dofs, xmax, Omega, Omega_lit)
TARGETS = [
    ('r150', 1.5, 0.25, 1, 20, 20.0, 1.1845658267167267, 46.7648),
    ('r150', 1.5, 0.25, 1, 20, 20.0, 1.2224258528594245, 48.2594),
    ('r150', 1.5, 0.25, 1, 20, 20.0, 1.2970050321098103, 51.2037),
    ('r150', 1.5, 0.25, 1, 20, 20.0, 2.3004974270051886, 90.8200),
    ('r150', 1.5, 0.5, 1, 20, 20.0, 0.017219531036775974, 0.6798),
    ('r150', 1.5, 0.5, 1, 20, 20.0, 0.3168420614324993, 12.5084),
    ('r150', 1.5, 0.5, 1, 20, 20.0, 0.3766154716686498, 14.8682),
    ('r150', 1.5, 0.5, 1, 20, 20.0, 0.5957690495444217, 23.5200),
    ('r150', 1.5, 0.5, 1, 20, 20.0, 0.8304497559844244, 32.7848),
    ('r150', 1.5, 0.75, 1, 16, 14.0, 0.01700729907566782, 0.6714),
    ('r150', 1.5, 0.75, 1, 16, 14.0, 0.13968917668204017, 5.5147),
    ('r150', 1.5, 0.75, 1, 16, 14.0, 0.31964230058270193, 12.6190),
    ('r150', 1.5, 0.75, 1, 16, 14.0, 0.37434067790811687, 14.7784),
    ('r150', 1.5, 0.75, 1, 16, 14.0, 0.6102404500186518, 24.0913),
    ('r150', 1.5, 1.0, 1, 16, 14.0, 0.07713364426395386, 3.0451),
    ('r150', 1.5, 1.0, 1, 16, 14.0, 0.17493079988453508, 6.9060),
    ('r150', 1.5, 1.0, 1, 16, 14.0, 0.18618328846182353, 7.3502),
    ('r150', 1.5, 1.0, 1, 16, 14.0, 0.20520968244526755, 8.1014),
    ('r150', 1.5, 1.0, 1, 16, 14.0, 0.4046136329667267, 15.9735),
    ('r150', 1.5, 1.25, 1, 16, 14.0, 0.04948633958033019, 1.9536),
    ('r150', 1.5, 1.25, 1, 16, 14.0, 0.09509413943045761, 3.7542),
    ('r150', 1.5, 1.25, 1, 16, 14.0, 0.10230644816433546, 4.0389),
    ('r150', 1.5, 1.25, 1, 16, 14.0, 0.12503968131238688, 4.9364),
    ('r150', 1.5, 1.25, 1, 16, 14.0, 0.16164485363967768, 6.3815),
    ('r150', 1.5, 1.5, 1, 16, 14.0, 0.03533152345909243, 1.3948),
    ('r150', 1.5, 1.5, 1, 16, 14.0, 0.0815077229976688, 3.2178),
    ('r150', 1.5, 1.5, 1, 16, 14.0, 0.1273284011188086, 5.0267),
    ('r150', 1.5, 1.5, 1, 16, 14.0, 0.17154373792645178, 6.7723),
    ('r150', 1.5, 1.5, 1, 16, 14.0, 0.2935357602236787, 11.5883),
    ('r150', 1.5, 0.25, 2, 20, 20.0, 0.021597083158334855, 26.7983),
    ('r150', 1.5, 0.25, 2, 20, 20.0, 0.32091508807189967, 398.2010),
    ('r150', 1.5, 0.25, 2, 20, 20.0, 0.8314869268367314, 1031.7337),
    ('r150', 1.5, 0.25, 2, 20, 20.0, 1.0319214020748113, 1280.4387),
    ('r150', 1.5, 0.5, 2, 20, 20.0, 0.32091508807189967, 398.2010),
    ('r150', 1.5, 0.5, 2, 20, 20.0, 0.4078246413998049, 506.0409),
    ('r150', 1.5, 0.5, 2, 20, 20.0, 0.4488250964128827, 556.9155),
    ('r150', 1.5, 0.5, 2, 20, 20.0, 0.5472823500863455, 679.0842),
    ('r150', 1.5, 0.75, 2, 14, 14.0, 0.1996628231474411, 247.7476),
    ('r150', 1.5, 0.75, 2, 14, 14.0, 0.3008667218013114, 373.3244),
    ('r150', 1.5, 0.75, 2, 14, 14.0, 0.32091508807189967, 398.2010),
    ('r150', 1.5, 0.75, 2, 14, 14.0, 0.3368257776302147, 417.9434),
    ('r150', 1.5, 1.0, 2, 14, 14.0, 0.11265364642321388, 139.7840),
    ('r150', 1.5, 1.0, 2, 14, 14.0, 0.18018076496983118, 223.5736),
    ('r150', 1.5, 1.0, 2, 14, 14.0, 0.2737727958568922, 339.7054),
    ('r150', 1.5, 1.0, 2, 14, 14.0, 0.42402386529826985, 526.1414),
    ('r150', 1.5, 1.25, 2, 14, 14.0, 0.07075692311734788, 87.7973),
    ('r150', 1.5, 1.25, 2, 14, 14.0, 0.11298609079896924, 140.1965),
    ('r150', 1.5, 1.25, 2, 14, 14.0, 0.1788117477311018, 221.8749),
    ('r150', 1.5, 1.25, 2, 14, 14.0, 0.21232597372259393, 263.4604),
    ('r150', 1.5, 1.5, 2, 14, 14.0, 0.048597969826797996, 60.3018),
    ('r150', 1.5, 1.5, 2, 14, 14.0, 0.07300772597033979, 90.5902),
    ('r150', 1.5, 1.5, 2, 14, 14.0, 0.12006022296098953, 148.9743),
    ('r150', 1.5, 1.5, 2, 14, 14.0, 0.1571892265509104, 195.0451),
    ('r167', 1.66667, 0.25, 1, 20, 20.0, 1.0069227361972883, 46.6532),
    ('r167', 1.66667, 0.25, 1, 20, 20.0, 1.085743479353516, 50.3052),
    ('r167', 1.66667, 0.25, 1, 20, 20.0, 1.0964088407056487, 50.7993),
    ('r167', 1.66667, 0.25, 1, 20, 20.0, 1.1448003202126054, 53.0414),
    ('r167', 1.66667, 0.25, 1, 20, 20.0, 2.1945734179866747, 101.6800),
    ('r167', 1.66667, 0.5, 1, 20, 20.0, 0.26449634960990875, 12.2548),
    ('r167', 1.66667, 0.5, 1, 20, 20.0, 0.3176510966214744, 14.7176),
    ('r167', 1.66667, 0.5, 1, 20, 20.0, 0.5180450585953696, 24.0023),
    ('r167', 1.66667, 0.5, 1, 20, 20.0, 0.7028142614864562, 32.5631),
    ('r167', 1.66667, 0.5, 1, 20, 20.0, 0.7974079249908315, 36.9459),
    ('r167', 1.66667, 0.75, 1, 16, 14.0, 0.013592256954915671, 0.6298),
    ('r167', 1.66667, 0.75, 1, 16, 14.0, 0.11522951139529942, 5.3389),
    ('r167', 1.66667, 0.75, 1, 16, 14.0, 0.2686788152471268, 12.4485),
    ('r167', 1.66667, 0.75, 1, 16, 14.0, 0.3131999639676649, 14.5113),
    ('r167', 1.66667, 0.75, 1, 16, 14.0, 0.5165171945503693, 23.9315),
    ('r167', 1.66667, 1.0, 1, 16, 14.0, 0.06304060411434029, 2.9208),
    ('r167', 1.66667, 1.0, 1, 16, 14.0, 0.14458580681733468, 6.6990),
    ('r167', 1.66667, 1.0, 1, 16, 14.0, 0.1643198757284942, 7.6133),
    ('r167', 1.66667, 1.0, 1, 16, 14.0, 0.17029563965659605, 7.8902),
    ('r167', 1.66667, 1.0, 1, 16, 14.0, 0.33002159974507583, 15.2907),
    ('r167', 1.66667, 1.25, 1, 16, 14.0, 0.08384304228100697, 3.8847),
    ('r167', 1.66667, 1.25, 1, 16, 14.0, 0.1030800556876784, 4.7760),
    ('r167', 1.66667, 1.25, 1, 16, 14.0, 0.14051151002851686, 6.5102),
    ('r167', 1.66667, 1.25, 1, 16, 14.0, 0.21300090190530507, 9.8689),
    ('r167', 1.66667, 1.25, 1, 16, 14.0, 0.2383758250862774, 11.0445),
    ('r167', 1.66667, 1.5, 1, 16, 14.0, 0.014829728221165663, 0.6871),
    ('r167', 1.66667, 1.5, 1, 16, 14.0, 0.028623270046451314, 1.3262),
    ('r167', 1.66667, 1.5, 1, 16, 14.0, 0.10728603442965776, 4.9708),
    ('r167', 1.66667, 1.5, 1, 16, 14.0, 0.14213231763198653, 6.5853),
    ('r167', 1.66667, 1.5, 1, 16, 14.0, 0.24473182177639718, 11.3390),
    ('r167', 1.66667, 0.25, 2, 20, 20.0, 0.017750528421302515, 25.8493),
    ('r167', 1.66667, 0.25, 2, 20, 20.0, 0.2608198270311902, 379.8204),
    ('r167', 1.66667, 0.25, 2, 20, 20.0, 0.7382544224378805, 1075.0874),
    ('r167', 1.66667, 0.25, 2, 20, 20.0, 0.9255694342382728, 1347.8660),
    ('r167', 1.66667, 0.5, 2, 20, 20.0, 0.017750090003437254, 25.8487),
    ('r167', 1.66667, 0.5, 2, 20, 20.0, 0.26083568785109584, 379.8435),
    ('r167', 1.66667, 0.5, 2, 20, 20.0, 0.3448946745294161, 502.2549),
    ('r167', 1.66667, 0.5, 2, 20, 20.0, 0.3929104819906048, 572.1782),
    ('r167', 1.66667, 0.75, 2, 14, 14.0, 0.017754415862888617, 25.8550),
    ('r167', 1.66667, 0.75, 2, 14, 14.0, 0.16561122630627279, 241.1723),
    ('r167', 1.66667, 0.75, 2, 14, 14.0, 0.2639378608921278, 384.3611),
    ('r167', 1.66667, 0.75, 2, 14, 14.0, 0.29537115845116496, 430.1360),
    ('r167', 1.66667, 1.0, 2, 14, 14.0, 0.09259441629389714, 134.8412),
    ('r167', 1.66667, 1.0, 2, 14, 14.0, 0.15558406247003548, 226.5702),
    ('r167', 1.66667, 1.0, 2, 14, 14.0, 0.23979945291511973, 349.2094),
    ('r167', 1.66667, 1.0, 2, 14, 14.0, 0.26082876735947064, 379.8335),
    ('r167', 1.66667, 1.25, 2, 14, 14.0, 0.017758893231369588, 25.8615),
    ('r167', 1.66667, 1.25, 2, 14, 14.0, 0.05793207954129284, 84.3639),
    ('r167', 1.66667, 1.25, 2, 14, 14.0, 0.09611616255560554, 139.9697),
    ('r167', 1.66667, 1.25, 2, 14, 14.0, 0.148456453754022, 216.1906),
    ('r167', 1.66667, 1.5, 2, 14, 14.0, 0.039644880345483115, 57.7331),
    ('r167', 1.66667, 1.5, 2, 14, 14.0, 0.06150761113747087, 89.5708),
    ('r167', 1.66667, 1.5, 2, 14, 14.0, 0.09890983000317913, 144.0380),
    ('r167', 1.66667, 1.5, 2, 14, 14.0, 0.1346983632802919, 196.1553),
    ('r200', 2.0, 0.25, 1, 20, 20.0, 0.06120654546893885, 3.7755),
    ('r200', 2.0, 0.25, 1, 20, 20.0, 0.7500261101816087, 46.2654),
    ('r200', 2.0, 0.25, 1, 20, 20.0, 0.8033070719914708, 49.5520),
    ('r200', 2.0, 0.25, 1, 20, 20.0, 0.8845658406718172, 54.5645),
    ('r200', 2.0, 0.25, 1, 20, 20.0, 1.468229363603458, 90.5678),
    ('r200', 2.0, 0.5, 1, 20, 20.0, 0.19106982688120716, 11.7861),
    ('r200', 2.0, 0.5, 1, 20, 20.0, 0.23661873242581177, 14.5958),
    ('r200', 2.0, 0.5, 1, 20, 20.0, 0.4055601501261039, 25.0170),
    ('r200', 2.0, 0.5, 1, 20, 20.0, 0.5184990277355656, 31.9836),
    ('r200', 2.0, 0.5, 1, 20, 20.0, 0.5924098753479146, 36.5428),
    ('r200', 2.0, 0.75, 1, 16, 20.0, 0.009435372827814033, 0.5820),
    ('r200', 2.0, 0.75, 1, 16, 20.0, 0.0816439922929659, 5.0362),
    ('r200', 2.0, 0.75, 1, 16, 20.0, 0.19687795767737534, 12.1444),
    ('r200', 2.0, 0.75, 1, 16, 20.0, 0.22686433946188084, 13.9941),
    ('r200', 2.0, 0.75, 1, 16, 20.0, 0.34889196207572015, 21.5214),
    ('r200', 2.0, 1.0, 1, 16, 20.0, 0.04423861421980502, 2.7289),
    ('r200', 2.0, 1.0, 1, 16, 20.0, 0.10335833877305081, 6.3757),
    ('r200', 2.0, 1.0, 1, 16, 20.0, 0.12185626742846399, 7.5167),
    ('r200', 2.0, 1.0, 1, 16, 20.0, 0.13210100624178694, 8.1487),
    ('r200', 2.0, 1.0, 1, 16, 20.0, 0.24739341826347241, 15.2605),
    ('r200', 2.0, 1.25, 1, 16, 20.0, 0.06120654546893885, 3.7755),
    ('r200', 2.0, 1.25, 1, 16, 20.0, 0.07322136402834171, 4.5167),
    ('r200', 2.0, 1.25, 1, 16, 20.0, 0.10876737251655791, 6.7093),
    ('r200', 2.0, 1.25, 1, 16, 20.0, 0.15322163046433945, 9.4515),
    ('r200', 2.0, 1.25, 1, 16, 20.0, 0.20068007732137488, 12.3790),
    ('r200', 2.0, 1.5, 1, 16, 20.0, 0.019612589546385905, 1.2098),
    ('r200', 2.0, 1.5, 1, 16, 20.0, 0.047211165910281204, 2.9122),
    ('r200', 2.0, 1.5, 1, 16, 20.0, 0.06120654546893885, 3.7755),
    ('r200', 2.0, 1.5, 1, 16, 20.0, 0.07880882966523263, 4.8613),
    ('r200', 2.0, 1.5, 1, 16, 20.0, 0.10158447523257594, 6.2662),
    ('r200', 2.0, 0.25, 2, 20, 20.0, 0.012305507200203783, 23.8579),
    ('r200', 2.0, 0.25, 2, 20, 20.0, 0.17780585543876387, 344.7293),
    ('r200', 2.0, 0.25, 2, 20, 20.0, 0.5971621771821134, 1157.7759),
    ('r200', 2.0, 0.25, 2, 20, 20.0, 0.7424737199659655, 1439.5054),
    ('r200', 2.0, 0.5, 2, 20, 20.0, 0.012412470311521203, 24.0653),
    ('r200', 2.0, 0.5, 2, 20, 20.0, 0.31154568921900394, 604.0237),
    ('r200', 2.0, 0.5, 2, 20, 20.0, 0.39585705393249937, 767.4862),
    ('r200', 2.0, 0.5, 2, 20, 20.0, 0.543552314649419, 1053.8373),
    ('r200', 2.0, 0.75, 2, 14, 20.0, 0.11890723145318358, 230.5369),
    ('r200', 2.0, 0.75, 2, 14, 20.0, 0.17786850962695788, 344.8508),
    ('r200', 2.0, 0.75, 2, 14, 20.0, 0.20757593136270208, 402.4475),
    ('r200', 2.0, 0.75, 2, 14, 20.0, 0.2356305010032441, 456.8396),
    ('r200', 2.0, 1.0, 2, 14, 20.0, 0.06549785912176859, 126.9870),
    ('r200', 2.0, 1.0, 2, 14, 20.0, 0.1192285323525378, 231.1598),
    ('r200', 2.0, 1.0, 2, 14, 20.0, 0.17788393154322446, 344.8807),
    ('r200', 2.0, 1.0, 2, 14, 20.0, 0.19033110868435904, 369.0133),
    ('r200', 2.0, 1.25, 2, 14, 20.0, 0.04082314432664935, 79.1478),
    ('r200', 2.0, 1.25, 2, 14, 20.0, 0.0720062552343918, 139.6055),
    ('r200', 2.0, 1.25, 2, 14, 20.0, 0.08672582179969107, 168.1437),
    ('r200', 2.0, 1.25, 2, 14, 20.0, 0.10675006143649576, 206.9666),
    ('r200', 2.0, 1.5, 2, 14, 20.0, 0.04520827555375688, 87.6496),
    ('r200', 2.0, 1.5, 2, 14, 20.0, 0.07043955374328745, 136.5680),
    ('r200', 2.0, 1.5, 2, 14, 20.0, 0.10728838274570021, 208.0103),
    ('r200', 2.0, 1.5, 2, 14, 20.0, 0.14395181849149535, 279.0933),
    ('r250', 2.5, 0.25, 1, 20, 20.0, 0.005985350314739948, 0.5317),
    ('r250', 2.5, 0.25, 1, 20, 20.0, 0.4647442747821423, 41.2816),
    ('r250', 2.5, 0.25, 1, 20, 20.0, 0.5108137667149827, 45.3738),
    ('r250', 2.5, 0.25, 1, 20, 20.0, 0.5358270614236884, 47.5956),
    ('r250', 2.5, 0.25, 1, 20, 20.0, 0.6890791084902717, 61.2084),
    ('r250', 2.5, 0.5, 1, 20, 20.0, 0.0058307845169397006, 0.5179),
    ('r250', 2.5, 0.5, 1, 20, 20.0, 0.04234101070762228, 3.7610),
    ('r250', 2.5, 0.5, 1, 20, 20.0, 0.12600109787775232, 11.1922),
    ('r250', 2.5, 0.5, 1, 20, 20.0, 0.16625654561723807, 14.7680),
    ('r250', 2.5, 0.5, 1, 20, 20.0, 0.29461207501163605, 26.1693),
    ('r250', 2.5, 0.75, 1, 16, 20.0, 0.0058307845169397006, 0.5179),
    ('r250', 2.5, 0.75, 1, 16, 20.0, 0.05307227356945177, 4.7142),
    ('r250', 2.5, 0.75, 1, 16, 20.0, 0.11755436954468551, 10.4419),
    ('r250', 2.5, 0.75, 1, 16, 20.0, 0.13194371720658304, 11.7201),
    ('r250', 2.5, 0.75, 1, 16, 20.0, 0.15006743024235114, 13.3300),
    ('r250', 2.5, 1.0, 1, 16, 20.0, 0.028418233678429515, 2.5243),
    ('r250', 2.5, 1.0, 1, 16, 20.0, 0.0423419810622587, 3.7611),
    ('r250', 2.5, 1.0, 1, 16, 20.0, 0.0677339870019768, 6.0166),
    ('r250', 2.5, 1.0, 1, 16, 20.0, 0.07983201985516605, 7.0912),
    ('r250', 2.5, 1.0, 1, 16, 20.0, 0.10044034357340537, 8.9218),
    ('r250', 2.5, 1.25, 1, 16, 20.0, 0.017724005088368272, 1.5744),
    ('r250', 2.5, 1.25, 1, 16, 20.0, 0.04762541686484556, 4.2304),
    ('r250', 2.5, 1.25, 1, 16, 20.0, 0.07665390243054823, 6.8089),
    ('r250', 2.5, 1.25, 1, 16, 20.0, 0.10064515585516606, 8.9400),
    ('r250', 2.5, 1.25, 1, 16, 20.0, 0.1438164352403259, 12.7747),
    ('r250', 2.5, 1.5, 1, 16, 20.0, 0.05286746128769108, 4.6960),
    ('r250', 2.5, 1.5, 1, 16, 20.0, 0.0664883031475293, 5.9059),
    ('r250', 2.5, 1.5, 1, 16, 20.0, 0.10097654908836827, 8.9694),
    ('r250', 2.5, 1.5, 1, 16, 20.0, 0.11881637994551111, 10.5540),
    ('r250', 2.5, 1.5, 1, 16, 20.0, 0.25468633342890223, 22.6229),
    ('r250', 2.5, 0.25, 2, 20, 20.0, 0.008006359620285534, 22.3527),
    ('r250', 2.5, 0.25, 2, 20, 20.0, 0.02027382939666002, 56.6018),
    ('r250', 2.5, 0.25, 2, 20, 20.0, 0.11024778209624292, 307.7971),
    ('r250', 2.5, 0.25, 2, 20, 20.0, 0.4561929526517344, 1273.6299),
    ('r250', 2.5, 0.5, 2, 20, 20.0, 0.008006359620285534, 22.3527),
    ('r250', 2.5, 0.5, 2, 20, 20.0, 0.02279653844015061, 63.6449),
    ('r250', 2.5, 0.5, 2, 20, 20.0, 0.11022308076663856, 307.7282),
    ('r250', 2.5, 0.5, 2, 20, 20.0, 0.1737057366805692, 484.9633),
    ('r250', 2.5, 0.75, 2, 14, 20.0, 0.00798682813564941, 22.2982),
    ('r250', 2.5, 0.75, 2, 14, 20.0, 0.15135191793823427, 422.5544),
    ('r250', 2.5, 0.75, 2, 14, 20.0, 0.17913990740505595, 500.1348),
    ('r250', 2.5, 0.75, 2, 14, 20.0, 0.20004446799309014, 558.4975),
    ('r250', 2.5, 1.0, 2, 14, 20.0, 0.042837285674276006, 119.5960),
    ('r250', 2.5, 1.0, 2, 14, 20.0, 0.08411476273185123, 234.8372),
    ('r250', 2.5, 1.0, 2, 14, 20.0, 0.11431126963872237, 319.1418),
    ('r250', 2.5, 1.0, 2, 14, 20.0, 0.14239282054336247, 397.5418),
    ('r250', 2.5, 1.25, 2, 14, 20.0, 0.02651951177869291, 74.0389),
    ('r250', 2.5, 1.25, 2, 14, 20.0, 0.04948263187238877, 138.1489),
    ('r250', 2.5, 1.25, 2, 14, 20.0, 0.10236580644555998, 285.7917),
    ('r250', 2.5, 1.25, 2, 14, 20.0, 0.11024778209624292, 307.7971),
    ('r250', 2.5, 1.5, 2, 14, 20.0, 0.030321561461866167, 84.6538),
    ('r250', 2.5, 1.5, 2, 14, 20.0, 0.04616341715013303, 128.8821),
    ('r250', 2.5, 1.5, 2, 14, 20.0, 0.07098884659466874, 198.1914),
    ('r250', 2.5, 1.5, 2, 14, 20.0, 0.12799955242400393, 357.3577),
]


def _worker(args):
    tag, r0_2b, two_T_pi, part, n_dofs, xmax, Om, Om_lit = args
    t0 = time.time()
    try:
        import plate_solver as ps
        from plate_solver.detectors import (full_search, select_fill,
                                             weak_enforcement_residual_oop,
                                             weak_enforcement_residual_ip)

        mat = ps.IsotropicMaterial(E=210e9, nu=NU_SWEEP, rho=7800.0)
        geom = ps.make_geometry(r0_2b, two_T_pi)

        if part == 1:
            solver = ps.OutOfPlaneSolver(geom, mat, M=80, n_quad=30,
                                          boundary=ps.FreeFreeOOP())
        else:
            solver = ps.InPlaneSolver(geom, mat, M=80, n_quad=30,
                                       boundary=ps.FreeFreeIP())

        raw = full_search(solver.fast, Om, xmax=xmax)
        sel, cnt = select_fill(raw, n_dofs)
        if cnt < n_dofs:
            return dict(ok=False, tag=tag, r0_2b=r0_2b, two_T_pi=two_T_pi,
                        part=part, Om=Om, Om_lit=Om_lit,
                        error=f"only {cnt}/{n_dofs} branches filled",
                        dt=time.time() - t0)

        if part == 1:
            res = weak_enforcement_residual_oop(solver, Om, sel, n_dofs=n_dofs)
        else:
            res = weak_enforcement_residual_ip(solver, Om, sel, n_dofs=n_dofs)

        if not res.get("ok", False):
            return dict(ok=False, tag=tag, r0_2b=r0_2b, two_T_pi=two_T_pi,
                        part=part, Om=Om, Om_lit=Om_lit,
                        error=f"residual screen failed: {res.get('reason')}",
                        dt=time.time() - t0)

        s = solver.sigma_min(Om, sel)

        out = dict(ok=True, tag=tag, r0_2b=r0_2b, two_T_pi=two_T_pi,
                   part=part, Om=Om, Om_lit=Om_lit, block=res["block"],
                   log10_sigma_min=float(s), dt=time.time() - t0)
        if part == 1:
            out.update(rV=res["rV"], rM=res["rM"], maxW=res["maxW"],
                       verdict=res["verdict"])
        else:
            # weak_enforcement_residual_ip returns raw ratios only, no
            # verdict bar -- report whatever keys it actually returns
            # rather than assuming a fixed name (rTyy/rTyr/maxDisp per
            # the module docstring, but read from res itself so this
            # doesn't silently drop a differently-named key).
            extra = {k: v for k, v in res.items() if k not in ("ok", "block")}
            out.update(extra)
        return out
    except Exception as exc:
        return dict(ok=False, tag=tag, r0_2b=r0_2b, two_T_pi=two_T_pi,
                    part=part, Om=Om, Om_lit=Om_lit, error=f"{exc}",
                    dt=time.time() - t0)


def main():
    print("=" * 78)
    print(f"  probe_ff_residual_screen_geomsweep -- RATIO_TAG={RATIO_TAG}"
          f" -- start " + time.strftime("%Y-%m-%d %H:%M:%S"))
    print("=" * 78, flush=True)

    import plate_solver as ps
    if ps.SOLVER_VERSION != EXPECT_VER:
        print(f"PREFLIGHT FAIL: SOLVER_VERSION={ps.SOLVER_VERSION!r} != "
              f"{EXPECT_VER!r}")
        raise SystemExit(2)
    print(f"preflight OK: SOLVER_VERSION={ps.SOLVER_VERSION}", flush=True)

    targets = [t for t in TARGETS if t[0] == RATIO_TAG]
    if not targets:
        print(f"FATAL: no targets for RATIO_TAG={RATIO_TAG!r} -- expected "
              f"one of r150/r167/r200/r250.")
        raise SystemExit(2)

    n_workers = int(os.environ.get("N_WORKERS", str(len(targets))))
    print(f"N_WORKERS={n_workers}  NU_SWEEP={NU_SWEEP}")
    n_p1 = sum(1 for t in targets if t[3] == 1)
    n_p2 = sum(1 for t in targets if t[3] == 2)
    print(f"{len(targets)} targets for {RATIO_TAG}: {n_p1} OOP (part 1) + "
          f"{n_p2} IP (part 2)\n", flush=True)

    results = []
    with ProcessPoolExecutor(max_workers=min(n_workers, len(targets))) as ex:
        futs = {ex.submit(_worker, t): t for t in targets}
        done = 0
        for fut in as_completed(futs):
            r = fut.result()
            done += 1
            partlabel = "OOP" if r["part"] == 1 else "IP "
            lit = f"{r['Om_lit']:.3f}" if r["Om_lit"] is not None else "n/a"
            if r["ok"]:
                if r["part"] == 1:
                    print(f"  [{done}/{len(targets)}] {partlabel} "
                          f"r0_2b={r['r0_2b']:g} 2T/pi={r['two_T_pi']:g} "
                          f"Om={r['Om']:.6f} (lit={lit}) block={r['block']:>3s} "
                          f"rV={r['rV']:.4g} rM={r['rM']:.4g} maxW={r['maxW']:.4g} "
                          f"log10smin={r['log10_sigma_min']:.4f} "
                          f"verdict={r['verdict']} [{r['dt']:.0f}s]", flush=True)
                else:
                    extra_str = " ".join(
                        f"{k}={v:.4g}" if isinstance(v, float) else f"{k}={v}"
                        for k, v in r.items()
                        if k not in ("ok", "tag", "r0_2b", "two_T_pi", "part",
                                     "Om", "Om_lit", "dt", "log10_sigma_min",
                                     "block"))
                    print(f"  [{done}/{len(targets)}] {partlabel} "
                          f"r0_2b={r['r0_2b']:g} 2T/pi={r['two_T_pi']:g} "
                          f"Om={r['Om']:.6f} (lit={lit}) block={r['block']:>3s} "
                          f"{extra_str} log10smin={r['log10_sigma_min']:.4f} "
                          f"(NO verdict bar -- IP exploratory, see docstring) "
                          f"[{r['dt']:.0f}s]", flush=True)
            else:
                print(f"  [{done}/{len(targets)}] {partlabel} "
                      f"r0_2b={r['r0_2b']:g} 2T/pi={r['two_T_pi']:g} "
                      f"Om={r['Om']:.6f} FAILED: {r.get('error')} "
                      f"[{r['dt']:.0f}s]", flush=True)
            results.append(r)

    ok_results = [r for r in results if r["ok"]]
    failures = [r for r in results if not r["ok"]]
    oop_ok = [r for r in ok_results if r["part"] == 1]
    real_like = [r for r in oop_ok if r["verdict"] == "REAL-like"]
    artifact_like = [r for r in oop_ok if r["verdict"] == "ARTIFACT-like"]
    ambiguous = [r for r in oop_ok if r["verdict"] == "AMBIGUOUS"]

    print("\n" + "=" * 78)
    print(f"SUMMARY ({RATIO_TAG}): {len(ok_results)}/{len(targets)} ok, "
          f"{len(failures)} failures")
    print(f"  OOP (part 1): {len(oop_ok)} screened -- "
          f"{len(real_like)} REAL-like, {len(artifact_like)} ARTIFACT-like, "
          f"{len(ambiguous)} AMBIGUOUS")
    print(f"  IP  (part 2): {sum(1 for r in ok_results if r['part']==2)} "
          f"screened -- raw ratios only, no verdict bar (see docstring)")
    if ambiguous:
        print("  AMBIGUOUS rows (do not trust either label without a closer look):")
        for r in ambiguous:
            print(f"    r0_2b={r['r0_2b']:g} 2T/pi={r['two_T_pi']:g} "
                  f"Om={r['Om']:.6f} rM={r['rM']:.4g}")
    if failures:
        print("FAILURES:")
        for r in failures:
            print(f"  r0_2b={r['r0_2b']:g} 2T/pi={r['two_T_pi']:g} "
                  f"part={r['part']} Om={r['Om']:.6f}: {r.get('error')}")
    print("=" * 78, flush=True)

    print("\nDiagnostic/data-generation only -- no SOLVER_VERSION action, "
          "no package changes. Read the module docstring's HONEST SCOPE "
          "CAVEAT before using any verdict in a paper claim.", flush=True)
    print("End: " + time.strftime("%Y-%m-%d %H:%M:%S"), flush=True)


if __name__ == "__main__":
    t0 = time.time()
    main()
    print(f"\ntotal wall time: {(time.time() - t0)/60:.1f} min", flush=True)
