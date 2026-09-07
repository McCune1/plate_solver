#!/usr/bin/env python3
"""
Rectangular FFFF release gate -- 2026-09-07.

WHY THIS EXISTS
---------------
Paper 2 review pass #3 found that `github_repo/plate_solver/core_solvers.py`
(the copy that is published, and the copy the repo digest is built from) had
silently fallen a month behind the deployed cluster tree. It still carried the
pre-2026-09-03 free-free corner term

    corner_ff = (T - nu) * (-2) * (TT + TB + WT + WB)

which is the naive same-sign four-term sum. That sum vanishes IDENTICALLY under
the parity identity (Supplementary S.1.3), so it assembles without raising and
silently produces the pre-corner spectrum. The released `RectIPAssembler` also
had no `bc` argument at all, so in-plane FFFF could not even be constructed.
Neither defect throws; both make every tabulated number in Paper 2
unreproducible from the released artifact.

The `.tex`-only currency check in PAPER2_REVIEWER_PROMPT.md cannot see any of
this. This probe is the check that can.

WHAT IT CHECKS (all pre-registered; no tuning, no retuning, read-only)
---------------------------------------------------------------------
  G0  SOLVER_VERSION matches EXPECT_SOLVER_VERSION.
  G1  The deployed free-free corner term is the closed-contour checkerboard
      (TT - TB - WT + WB) and is NOT the naive same-sign sum.
  G2  RectIPAssembler exposes bc='free_free' and rejects an unknown bc.
  G3  Two single-point sigma_min anchors reproduce their stored values to
      RTOL, at every dps in DPS_LIST:
        OOP  l/b=1.5 SYM  Lambda* = 0.964   basis (6, 0)  im_cap 30.0
        IP   l/b=2.0 ANTI Omega*  = 0.5200  basis (5, 3)  im_cap  7.0
      Both are stable to every printed digit at dps = 26, 30 and 40.
  G4  In-plane free-free and clamped-free assemble to materially different
      matrices at the same point (guards a silently-ignored bc).
  G5  (optional, set RELEASE_PKG) every plate_solver/*.py in the release tree
      is byte-identical to the deployed tree.

PRE-REGISTERED INTERPRETATION -- decide BEFORE reading the output
-----------------------------------------------------------------
  RECT_FF_RELEASE_GATE: PASS
        The deployed tree is the checkerboard build and both anchors are
        exactly where Paper 2 left them. If G5 also ran, the release tree is
        the same code. Safe to cut a manuscript-matching tag.

  RECT_FF_RELEASE_GATE: FAIL_CORNER_FORMULA
        The tree being tested is the pre-2026-09-03 assembler. Do not publish
        it, do not build a digest from it, and do not trust any free-free
        rectangular number it produces.

  RECT_FF_RELEASE_GATE: FAIL_IP_BC
        RectIPAssembler has no working free-free option. Same conclusion for
        the 92-row in-plane table.

  RECT_FF_RELEASE_GATE: FAIL_ANCHOR
        The corner formula is right but a stored sigma_min moved. Two very
        different causes, and the printed library versions below tell them
        apart:

          (a) SAME environment as the capture line below -> a code change
              ALTERED A COMPUTED RECTANGULAR FFFF FREQUENCY. That is exactly
              the condition SOLVER_VERSION exists to record: bump it, and
              re-derive Tables 1-2 and the 92-row list before touching the
              paper.

          (b) DIFFERENT scipy/numpy/mpmath -> re-run under the capture
              environment before concluding anything about the code. Note
              that library sensitivity has NOT been observed here: these
              anchors agree to ~15 significant digits between scipy 1.17.1 /
              numpy 2.4.4 and the Mill venv at scipy 1.18.0 / numpy 2.5.1,
              against a 1e-12 relative tolerance. So (b) is the unlikely
              branch; check (a) first.

        In NEITHER case edit a stored anchor merely to make this pass. On
        2026-09-07 three ANNULAR gates in tests/test_solver.py failed on two
        different machines with identical values -- the anchors were stale,
        not the solver -- and the fix was to re-derive them against an
        independent record of what they should be (LESSONS_LEARNED Sec. 1 /
        Sec. 16 / Sec. 18.34), never to overwrite them with whatever came
        out.

  RECT_FF_RELEASE_GATE: FAIL_RELEASE_DRIFT
        Deployed and released trees disagree. Sync release <- deployed (never
        the other way), rebuild the digest, re-run this probe.

Runtime: four single-point probes per dps, ~0.5 s each in a container and a
few seconds on a cluster node. The whole gate is well under a minute; it is
cheap enough to run before every tag.

Env:
  PKG_PATH                 package root to import from (default ".")
  EXPECT_SOLVER_VERSION    default "2026-07-10.s10"
  DPS_LIST                 comma-separated, default "26,30,40"
  RELEASE_PKG              optional path to the release tree's package root
                           (e.g. ".../github_repo") to enable G5
  RTOL                     default 1e-12

Does NOT retune PERSIST_DL or the 0.744 MAC bar, does NOT edit plate_solver,
does NOT bump SOLVER_VERSION, does NOT touch a tabulated frequency.
"""
from __future__ import annotations

import os
import sys
import hashlib

PKG_PATH = os.environ.get("PKG_PATH", ".")
EXPECT = os.environ.get("EXPECT_SOLVER_VERSION", "2026-07-10.s10")
DPS_LIST = [int(x) for x in os.environ.get("DPS_LIST", "26,30,40").split(",")]
RELEASE_PKG = os.environ.get("RELEASE_PKG", "").strip()
RTOL = float(os.environ.get("RTOL", "1e-12"))

sys.path.insert(0, PKG_PATH)

CHECKERBOARD = "p1t * Ut1 * q0t * Vt0 - p1b * Ut1 * q0b * Vt0"
NAIVE_SUM = "p1t * Ut1 * q0t * Vt0 + p1b * Ut1 * q0b * Vt0"
IP_BC_MARKER = "RectIPAssembler: unknown bc"

# Stored anchors. Computed 2026-09-07 on the checkerboard build at
# SOLVER_VERSION 2026-07-10.s10; identical at dps 26, 30 and 40.
# For reference, the pre-2026-09-03 naive-sum build gives 3.950702663724406e-03
# for the OOP anchor -- a factor of 16.6, not a rounding difference.
# Environment the anchors were captured under, 2026-09-07. sigma_min at a
# fixed off-root Lambda is a conditioning number, and the annular gates in
# tests/test_solver.py are known to move between library versions -- so the
# environment is part of the anchor, not incidental to it.
CAPTURE_ENV = "scipy 1.17.1 / numpy 2.4.4 / mpmath 1.3.0 and 1.4.1 (identical)"

ANCHORS = [
    # tag,                lob, sym,   x,      n_real, n_cpair, im_cap, in_plane, n_expect, sigma
    ("OOP l/b=1.5 SYM  Lambda*=0.964",
     1.5, True, 0.964, 6, 0, 30.0, False, 3, 2.3835918992545717e-04),
    ("IP  l/b=2.0 ANTI Omega*=0.5200",
     2.0, False, 0.5200, 5, 3, 7.0, True, 4, 3.900364925307852e-08),
]

FAILURES = []


def fail(code, msg):
    FAILURES.append((code, msg))
    print("  FAIL [%s] %s" % (code, msg), flush=True)


def main():
    print("=" * 72, flush=True)
    print("Rectangular FFFF release gate  (probe_rect_ff_release_gate_2026-09-07)",
          flush=True)
    print("PKG_PATH=%s  DPS_LIST=%s  RTOL=%g" % (PKG_PATH, DPS_LIST, RTOL),
          flush=True)
    try:
        import scipy, numpy, mpmath as _mp
        here = "scipy %s / numpy %s / mpmath %s" % (
            scipy.__version__, numpy.__version__, _mp.__version__)
    except Exception as exc:                       # pragma: no cover
        here = "unavailable (%r)" % (exc,)
    print("this run : %s" % here, flush=True)
    print("captured : %s" % CAPTURE_ENV, flush=True)
    print("=" * 72, flush=True)

    # ---- G0 SOLVER_VERSION -------------------------------------------------
    import plate_solver as ps
    print("\n[G0] SOLVER_VERSION")
    print("  found  : %s" % ps.SOLVER_VERSION)
    print("  expect : %s" % EXPECT)
    if ps.SOLVER_VERSION != EXPECT:
        fail("FAIL_ANCHOR", "SOLVER_VERSION mismatch")
    else:
        print("  ok", flush=True)

    import plate_solver.core_solvers as cs
    with open(cs.__file__, encoding="utf-8") as fh:
        src = fh.read()

    # ---- G1 corner formula -------------------------------------------------
    print("\n[G1] free-free four-corner Kirchhoff jump")
    if CHECKERBOARD not in src:
        fail("FAIL_CORNER_FORMULA",
             "checkerboard (TT - TB - WT + WB) not found in core_solvers.py")
    elif NAIVE_SUM in src:
        fail("FAIL_CORNER_FORMULA",
             "naive same-sign sum (TT + TB + WT + WB) present -- this is the "
             "pre-2026-09-03 assembler, and it vanishes identically")
    else:
        print("  checkerboard present, naive sum absent -- ok", flush=True)

    # ---- G2 in-plane bc switch --------------------------------------------
    print("\n[G2] RectIPAssembler free-free option")
    if IP_BC_MARKER not in src:
        fail("FAIL_IP_BC", "RectIPAssembler has no bc switch in source")
    else:
        try:
            ps.RectIPAssembler(nu_b=0.3, R=1.0, dps=26, bc="free_free")
            print("  bc='free_free' constructs -- ok", flush=True)
        except TypeError as exc:
            fail("FAIL_IP_BC", "bc kwarg rejected: %r" % (exc,))
        try:
            ps.RectIPAssembler(nu_b=0.3, R=1.0, dps=26, bc="nonsense")
            fail("FAIL_IP_BC", "unknown bc was accepted silently")
        except ValueError:
            print("  unknown bc rejected -- ok", flush=True)
        except TypeError:
            pass  # already reported above

    # ---- G3 / G4 numeric anchors ------------------------------------------
    from mpmath import mp, mpf
    from plate_solver.detectors import (
        rect_resolve_branches, rect_select_branches,
        rect_ip_resolve_branches, rect_ip_select_branches,
    )

    def sigma(asm, sym, lob, x, n_real, n_cpair, im_cap, in_plane):
        if in_plane:
            reps = rect_ip_resolve_branches(asm.eng(sym), x, im_cap=im_cap)
            full = rect_ip_select_branches(reps, n_real=n_real, n_cpair=n_cpair)
        else:
            reps = rect_resolve_branches(asm.eng(sym), x, im_cap=im_cap)
            full = rect_select_branches(reps, n_real=n_real, n_cpair=n_cpair)
        K = asm.assemble(full, mpf(str(round(x, 6))), sym, lob)
        return float(asm.equil_sigma(K, full)), len(full)

    print("\n[G3] stored sigma_min anchors")
    for dps in DPS_LIST:
        mp.dps = dps
        for (tag, lob, sym, x, nr, nc, cap, ip, n_exp, ref) in ANCHORS:
            if ip:
                asm = ps.RectIPAssembler(nu_b=0.3, R=1.0, dps=mp.dps,
                                         bc="free_free")
            else:
                asm = ps.RectOOPAssembler(nu_b=0.3, R=1.0, T=1.0, dps=mp.dps,
                                          bc="free_free")
            got, n = sigma(asm, sym, lob, x, nr, nc, cap, ip)
            rel = abs(got - ref) / ref
            ok = (n == n_exp) and (rel <= RTOL)
            print("  dps=%-2d %-34s n=%d (exp %d)  sigma=%.17e  rel=%.3e  %s"
                  % (dps, tag, n, n_exp, got, rel, "ok" if ok else "MOVED"),
                  flush=True)
            if not ok:
                fail("FAIL_ANCHOR",
                     "%s at dps=%d: got %.17e, stored %.17e (rel %.3e), "
                     "branches %d vs %d" % (tag, dps, got, ref, rel, n, n_exp))

    print("\n[G4] in-plane free-free differs from clamped-free")
    mp.dps = 26
    ff = ps.RectIPAssembler(nu_b=0.3, R=1.0, dps=mp.dps, bc="free_free")
    cf = ps.RectIPAssembler(nu_b=0.3, R=1.0, dps=mp.dps, bc="clamped_free")
    s_ff, _ = sigma(ff, False, 2.0, 0.5200, 5, 3, 7.0, True)
    s_cf, _ = sigma(cf, False, 2.0, 0.5200, 5, 3, 7.0, True)
    ratio = abs(s_cf - s_ff) / s_cf
    print("  free_free %.6e   clamped_free %.6e   rel diff %.4f"
          % (s_ff, s_cf, ratio), flush=True)
    if ratio <= 0.5:
        fail("FAIL_IP_BC",
             "bc='free_free' barely changes the in-plane assembly (rel diff "
             "%.4f) -- the switch may be ignored" % ratio)
    else:
        print("  ok", flush=True)

    # ---- G5 release-vs-deployed drift -------------------------------------
    print("\n[G5] release tree vs deployed tree")
    if not RELEASE_PKG:
        print("  skipped (set RELEASE_PKG to the release package root to run)",
              flush=True)
    else:
        dep = os.path.join(PKG_PATH, "plate_solver")
        rel_dir = os.path.join(RELEASE_PKG, "plate_solver")
        if not os.path.isdir(rel_dir):
            fail("FAIL_RELEASE_DRIFT", "no plate_solver/ under %s" % RELEASE_PKG)
        else:
            names = sorted(f for f in os.listdir(dep) if f.endswith(".py"))
            drift = []
            for f in names:
                a = os.path.join(dep, f)
                b = os.path.join(rel_dir, f)
                if not os.path.exists(b):
                    drift.append((f, "missing in release"))
                    continue
                ha = hashlib.md5(open(a, "rb").read()).hexdigest()
                hb = hashlib.md5(open(b, "rb").read()).hexdigest()
                if ha != hb:
                    drift.append((f, "%s != %s" % (ha[:12], hb[:12])))
            print("  compared %d modules" % len(names), flush=True)
            if drift:
                for f, why in drift:
                    print("    DRIFT %-22s %s" % (f, why), flush=True)
                fail("FAIL_RELEASE_DRIFT",
                     "%d module(s) differ between deployed and release"
                     % len(drift))
            else:
                print("  all modules byte-identical -- ok", flush=True)

    # ---- verdict -----------------------------------------------------------
    print("\n" + "=" * 72, flush=True)
    if not FAILURES:
        print("RECT_FF_RELEASE_GATE: PASS", flush=True)
        print("SOLVER_VERSION %s unbumped. PERSIST_DL not retuned. "
              "MAC bar 0.744 not retuned." % ps.SOLVER_VERSION, flush=True)
        return 0
    codes = []
    for code, _ in FAILURES:
        if code not in codes:
            codes.append(code)
    for code in codes:
        print("RECT_FF_RELEASE_GATE: %s" % code, flush=True)
    print("\n%d failure(s):" % len(FAILURES), flush=True)
    for code, msg in FAILURES:
        print("  [%s] %s" % (code, msg), flush=True)
    return 1


if __name__ == "__main__":
    sys.exit(main())
