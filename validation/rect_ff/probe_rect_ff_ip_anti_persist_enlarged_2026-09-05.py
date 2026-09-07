#!/usr/bin/env python3
"""
Widen the persist criterion for the IP FFFF ANTI branch's near-miss
FAIL_P discovery dips from job 2455570.

Job 2455570's ANTI persist-pass rate is far below SYM's (e.g. l/b=1.5
ANTI: 0/7 dips persisted; l/b=1.0 ANTI 1/6). Most of the FAIL_P dips
below are already within a few percent of a real FE target (the
production-basis discovery estimate O0 is close; it is the n_cpair=3
persist trial, in a window of only O0+/-0.02 at step 0.002, that fails
to converge within PERSIST_DL=0.01 of O0). Two distinct explanations are
live and this probe cannot fully separate them without also widening the
window: (a) n_cpair=3 is not enough basis for ANTI's persist criterion,
so the true minimum near O0 is shallow/displaced at that basis; (b) the
true minimum, once n_cpair is enlarged, sits further from O0 than the
+/-0.02 persist window reaches, so a wider window is required regardless
of basis. This probe changes BOTH at once (n_cpair 3->5, window
+/-0.02->+/-0.05, step unchanged at 0.002) and reports where the new
minimum actually sits relative to the FE target and to the original O0,
so a follow-up can attribute the effect if needed.

Target list is the FAIL_P dips from 2455570 with original nearest-FE
miss <= 6% (near-misses only; the >10% miss FAIL_P dips are almost
certainly unrelated artifacts and are excluded to keep this probe
cheap and focused):
  (lob, ANTI, O0, orig_miss%)
  1.0: 1.40 (5.28), 1.94 (3.15), 2.10 (4.83)
  1.5: 0.80 (1.55), 1.06 (2.58), 1.50 (4.53), 1.60 (1.83), 2.18 (1.02),
       2.36 (5.33)
  2.0: 1.82 (1.34), 1.92 (4.08)
  2.5: 1.08 (1.05), 1.42 (1.56), 1.82 (0.94), 2.02 (1.38), 2.10 (2.21),
       2.18 (1.52), 2.46 (4.74)
  3.0: 0.90 (1.42), 1.28 (0.77), 1.32 (1.96), 1.52 (3.72), 1.70 (2.10),
       1.88 (1.69), 1.94 (1.45), 2.22 (1.35)

PRE-REGISTERED
--------------
  RECOVERED if, at n_cpair=5 and window +/-0.05, a local minimum exists
  within PERSIST_DL=0.01 of O0's NEW location AND within 3% of the FE
  target -- i.e. the point now persists in the enlarged-basis sense.
  Report separately whether the new minimum sits within the OLD
  +/-0.02 window (basis-limited) or only in the wider +/-0.05 to
  +/-0.05 range (window-limited) -- this distinction matters for what
  to fix if the pattern is systematic.
  STILL FAILS if no qualifying minimum exists even at n_cpair=5 and the
  wider window -- these stay named FAIL_P artifacts, not evidence for a
  missing formula term (there is no IP corner jump, confirmed
  2026-09-05).
  Do not retune PERSIST_DL globally, do not call this Screen B, do not
  adopt n_cpair=5/window +/-0.05 as a new default without further
  validation.
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
CHECKPOINT_DIR = os.environ.get("CHECKPOINT_DIR", "./p5_ip_anti_wide_checkpoints")
PKG = os.environ.get("PKG_PATH", ".")
SYM = False  # ANTI only
IM_CAP = 7.0
N_CPAIR_WIDE = 5
HALF = 0.05
STEP = 0.002
MATCH_CUT = 3.0

# (lob, O0, orig_miss_pct) -- FAIL_P dips from 2455570, miss<=6% only.
TARGETS = [
    (1.0, 1.40, 5.28), (1.0, 1.94, 3.15), (1.0, 2.10, 4.83),
    (1.5, 0.80, 1.55), (1.5, 1.06, 2.58), (1.5, 1.50, 4.53),
    (1.5, 1.60, 1.83), (1.5, 2.18, 1.02), (1.5, 2.36, 5.33),
    (2.0, 1.82, 1.34), (2.0, 1.92, 4.08),
    (2.5, 1.08, 1.05), (2.5, 1.42, 1.56), (2.5, 1.82, 0.94),
    (2.5, 2.02, 1.38), (2.5, 2.10, 2.21), (2.5, 2.18, 1.52),
    (2.5, 2.46, 4.74),
    (3.0, 0.90, 1.42), (3.0, 1.28, 0.77), (3.0, 1.32, 1.96),
    (3.0, 1.52, 3.72), (3.0, 1.70, 2.10), (3.0, 1.88, 1.69),
    (3.0, 1.94, 1.45), (3.0, 2.22, 1.35),
]


def wide_persist(asm, lob, O0, ckpt_dir):
    tag = f"antiwide_lob{lob}_O{O0:.4f}".replace(".", "p")
    path = os.path.join(ckpt_dir, f"{tag}.json")
    if os.path.exists(path):
        with open(path) as f:
            return json.load(f)
    lo = max(0.002, round(O0 - HALF, 6))
    hi = round(O0 + HALF, 6)
    rows, dips = L.scan_window_ip(asm, lob, SYM, lo, hi, STEP, IM_CAP,
                                   5, N_CPAIR_WIDE)
    nearest = None
    if dips:
        nearest = min(dips, key=lambda t: abs(t[0] - O0))
    out = {
        "O0": O0, "lo": lo, "hi": hi,
        "nearest": None if nearest is None else [nearest[0], float(nearest[1])],
        "dips": [[a, float(b), c] for a, b, c in dips],
    }
    os.makedirs(ckpt_dir, exist_ok=True)
    with open(path, "w") as f:
        json.dump(out, f)
    return out


def main():
    print("probe_rect_ff_ip_anti_persist_enlarged_2026-09-05")
    print("start", time.strftime("%Y-%m-%d %H:%M:%S"), f"dps={mp.dps}")
    if not L.deploy_ok_ip(PKG, EXPECT_VER):
        raise SystemExit(2)
    os.makedirs(CHECKPOINT_DIR, exist_ok=True)
    asm = L.make_ff_ip()

    print(f"\n=== ANTI wide persist, n_cpair={N_CPAIR_WIDE} window +/-{HALF} "
          f"step {STEP} ({len(TARGETS)} targets) ===", flush=True)
    recovered = []
    for lob, O0, orig_miss in TARGETS:
        t0 = time.time()
        res = wide_persist(asm, lob, O0, CHECKPOINT_DIR)
        dt = time.time() - t0
        if res["nearest"] is None:
            print(f"  l/b={lob} ANTI O0={O0:.4f} (was {orig_miss:.2f}% miss)  "
                  f"NO DIP in +/-{HALF}  ({dt:.1f}s)", flush=True)
            continue
        Ostar, sstar = res["nearest"]
        Ofe, miss = L.nearest_ip_fe(lob, False, Ostar)
        in_old_window = abs(Ostar - O0) <= 0.02
        recovered_here = (miss is not None and miss <= MATCH_CUT)
        print(f"  l/b={lob} ANTI O0={O0:.4f} -> O*={Ostar:.4f} "
              f"sig*={sstar:.2e}  FE={Ofe} miss="
              f"{None if miss is None else '%.2f%%' % miss}  "
              f"in_old_+/-0.02_window={in_old_window}  ({dt:.1f}s)",
              flush=True)
        if recovered_here:
            recovered.append((lob, O0, Ostar, Ofe, miss, in_old_window))

    print(f"\n=== RECOVERED (n_cpair={N_CPAIR_WIDE}, within {MATCH_CUT:.0f}% "
          "of FE) ===", flush=True)
    if not recovered:
        print("  NONE -- ANTI weakness not explained by basis+window alone",
              flush=True)
    for lob, O0, Ostar, Ofe, miss, in_old in recovered:
        loc = "basis-limited (inside old +/-0.02 window)" if in_old \
            else "window-limited (only visible in wider +/-0.05 window)"
        print(f"  l/b={lob}  O0={O0:.4f} -> O*={Ostar:.4f}  FE={Ofe:.5f}  "
              f"miss={miss:.2f}%  {loc}")
    print("end", time.strftime("%Y-%m-%d %H:%M:%S"), flush=True)


if __name__ == "__main__":
    t0 = time.time()
    main()
    print(f"wall time: {time.time() - t0:.1f}s", flush=True)
