#!/usr/bin/env python3
"""
Follow-up to job 2455675 (fine reopen, +/-0.02 around l/b=1.0 SYM
FE=1.98249, persist basis (6,3), step 0.0002): the entire scanned window
sat on a nearly FLAT sigma plateau (~1.70-1.704e-6, monotonically
non-increasing toward the right edge Omega=2.0025) with NO 10x rise on
either flank -- not the shape of a genuine isolated notch, and not
obviously the tail of a distant root either (a true approach to a far
root should show sigma falling much faster, not sitting flat).

2455675 could not distinguish between two explanations because its
window was only 0.04 wide:
  (a) this is a real but very broad/shallow near-degenerate trough that
      eventually rises well outside +/-0.02 -- worth flagging as an
      unusual numerical feature even if it is still, in the end, a MISS
      against the sharp FE target, or
  (b) it is the leading edge of the tail approaching the NEXT SYM FE
      root at 2.45439 (FE_LISTS[(1.0, True)][1]), i.e. sigma keeps
      falling smoothly all the way from here to there with no structure
      in between.

This probe scans the full gap between the two FE targets, from just
above 1.98249 out to well short of 2.45439 (stopping at 2.23, comfortably
clear of that neighbor's own persist window so its real dip is not
conflated with this one), at persist basis (6,3), step 0.002 (250 pts,
~5 min at the 2455675 per-point rate -- cheap enough to just look).

Does NOT retune SIGMA_LIST, PERSIST_DL, or Screen B. Does NOT touch the
corner term. Purely descriptive: report the full sigma trace's local
minima (standard 3-point test), the global min, and whether sigma ever
rises 10x above its running minimum anywhere in the window (a "wall"
on either side of SOME point, even if not exactly at 1.98249).

PRE-REGISTERED
--------------
  REOPEN if a local_minima() dip appears anywhere in [1.73, 2.23] within
  3% of FE=1.98249 (looser than 2455675's 1% since we are now hunting
  for a shifted/broadened root, not confirming an exact match).
  FLAT-PLATEAU CONFIRMED (not a hidden root, not yet a wall-approach) if
  sigma stays within one order of magnitude of its 2455675 floor
  (~1.7e-6) across a wide span with no 10x rise anywhere short of 2.23.
  TAIL-OF-NEIGHBOR CONFIRMED if sigma falls monotonically and
  substantially (order(s) of magnitude) as Omega -> 2.23, i.e. it is
  approaching the 2.45439 root rather than sitting on a plateau.
  Either finding leaves the 2455568/2455675 MISS verdict on Table 3
  as-is; this is a diagnostic note for the SM/discussion, not a table
  change.
"""
from __future__ import annotations

import json
import os
import sys
import time

sys.path.insert(0, os.environ.get("PKG_PATH", "."))
os.environ.setdefault("DPS", "30")

from mpmath import mp  # noqa: E402

mp.dps = int(os.environ["DPS"])

import plate_solver as ps  # noqa: E402
import p5_rect_ff_lib as L  # noqa: E402

EXPECT_VER = os.environ.get("EXPECT_SOLVER_VERSION", "2026-07-10.s10")
CHECKPOINT_DIR = os.environ.get("CHECKPOINT_DIR", "./p5_1982_wide_gap_checkpoints")
PKG = os.environ.get("PKG_PATH", ".")
LOB = 1.0
SYM = True
LFE = 1.98249
LFE_NEXT = 2.45439
LO, HI, STEP = 1.73, 2.23, 0.002
IM_CAP = 7.0


def wide_scan(asm):
    path = os.path.join(CHECKPOINT_DIR, "wide_lob1p0_SYM_persist.json")
    if os.path.exists(path):
        print(f"  SKIP {path}", flush=True)
        with open(path) as f:
            return json.load(f)
    n_real, n_cpair = L.n_basis(SYM, persist=True)
    t0 = time.time()
    rows, dips = L.scan_window(asm, LOB, SYM, LO, HI, STEP, IM_CAP,
                                n_real, n_cpair)
    dt = time.time() - t0
    out = {
        "n_real": n_real, "n_cpair": n_cpair,
        "rows": [[Lam, (None if s is None else float(s)), n]
                 for Lam, s, n in rows],
        "dips": [[d[0], float(d[1]), d[2]] for d in dips],
        "elapsed_s": dt,
    }
    os.makedirs(CHECKPOINT_DIR, exist_ok=True)
    with open(path, "w") as f:
        json.dump(out, f)
    return out


def main():
    print("probe_rect_ff_oop_1982_wide_gap_2026-09-05")
    print("start", time.strftime("%Y-%m-%d %H:%M:%S"), f"dps={mp.dps}")
    if not L.deploy_ok(PKG, EXPECT_VER):
        raise SystemExit(2)
    os.makedirs(CHECKPOINT_DIR, exist_ok=True)
    asm = L.make_ff()

    print(f"\n=== wide gap scan l/b={LOB} SYM  window [{LO},{HI}] "
          f"step {STEP}  (FE={LFE:.5f} .. next FE={LFE_NEXT:.5f}) ===",
          flush=True)
    res = wide_scan(asm)
    print(f"  n_real={res['n_real']} n_cpair={res['n_cpair']} "
          f"{len(res['rows'])} pts {res['elapsed_s']:.1f}s", flush=True)

    finite = [(Lam, s, n) for Lam, s, n in res["rows"] if s is not None]
    if not finite:
        print("  ALL NONE -- solver returned no branches in this window",
              flush=True)
        print("=== SUMMARY: no data, inconclusive ===", flush=True)
        return

    Lam_min, s_min, n_min = min(finite, key=lambda t: t[1])
    print(f"  GLOBAL MIN over wide window  Omega={Lam_min:.4f}  "
          f"sigma={s_min:.3e}  n={n_min}", flush=True)

    dips = res["dips"]
    print("  dip detector (local_minima) over wide window: "
          + (", ".join(f"{a:.4f}({b:.2e})" for a, b, _c in dips) or "NONE"),
          flush=True)

    near_1982 = [d for d in dips
                 if abs(d[0] - LFE) / LFE * 100.0 <= 3.0]
    reopen = bool(near_1982)

    # running-min 10x-rise check: walk the ordered trace, track running
    # min so far, flag the first Omega where sigma >= 10x that running
    # min (a "wall" appearing anywhere after some point).
    ordered = sorted(finite, key=lambda t: t[0])
    running_min = None
    wall_at = None
    for Lam, s, n in ordered:
        if running_min is None or s < running_min:
            running_min = s
        elif s >= 10 * running_min and wall_at is None:
            wall_at = (Lam, s, running_min)
    if wall_at:
        print(f"  10x-rise wall first appears at Omega={wall_at[0]:.4f} "
              f"(sigma={wall_at[1]:.3e} vs running min {wall_at[2]:.3e})",
              flush=True)
    else:
        print("  NO 10x-rise wall anywhere in [1.73,2.23] relative to "
              "the running minimum -- sigma trend is monotone/flat "
              "throughout this window", flush=True)

    s_lo = ordered[0][1]
    s_hi = ordered[-1][1]
    print(f"  sigma at left edge (Omega={LO})={s_lo:.3e}   "
          f"sigma at right edge (Omega={HI})={s_hi:.3e}   "
          f"ratio lo/hi={ (s_lo / s_hi) if s_hi else float('nan'):.2f}",
          flush=True)

    print(f"\n  reopen-candidates within 3% of FE={LFE:.5f}: "
          + (", ".join(f"{a:.4f}({b:.2e})" for a, b, _c in near_1982)
             or "NONE"), flush=True)

    if reopen:
        print("=== SUMMARY REOPEN -- a dip within 3% of 1.98249 exists "
              "somewhere in the wider window ===", flush=True)
    elif wall_at is None and s_hi < 10 * s_min:
        print("=== SUMMARY FLAT-PLATEAU CONFIRMED -- sigma stays within "
              "an order of magnitude of its floor across the whole gap, "
              "no wall, no neighbor-tail signature; MISS stands, flag as "
              "numerical curiosity for SM/discussion ===", flush=True)
    elif s_hi < s_lo / 10:
        print("=== SUMMARY TAIL-OF-NEIGHBOR CONFIRMED -- sigma falls "
              "an order of magnitude or more toward Omega=2.23, "
              "consistent with approach to the 2.45439 root rather than "
              "a hidden mode near 1.98249; MISS stands ===", flush=True)
    else:
        print("=== SUMMARY AMBIGUOUS -- neither a clean flat plateau nor "
              "a clean monotone tail; report raw numbers to a human ===",
              flush=True)
    print("end", time.strftime("%Y-%m-%d %H:%M:%S"), flush=True)


if __name__ == "__main__":
    t0 = time.time()
    main()
    print(f"wall time: {time.time() - t0:.1f}s", flush=True)
