#!/usr/bin/env python3
"""
Re-open l/b=1.0 SYM FE=1.98249 (OOP) after job 2455568's own printed
sentinel: at the PERSIST basis (n_real=6, n_cpair=3) sigma AT the exact
FE frequency was 1.7304801007192585e-06 -- five orders of magnitude
below SIGMA_LIST=0.3 and comparable to genuine confirmed matches
elsewhere in this project -- yet the grid-based local-minima detector
(step 0.002, window +/-0.06) reported "dips: NONE" and the job was
scored MISS.

CLAUDE.md's own methodology section: "A standard sigma_min coarse scan
can silently miss real, narrow modes -- treat any 'N modes found' count
as a floor, not a ceiling, until spot-checked at a finer step." That is
exactly the shape of this discrepancy: a near-zero sigma at a single
point with no discrete 3-point dip nearby is the signature of a notch
narrower than the scan step, not evidence of no root.

This probe does NOT retune SIGMA_LIST, PERSIST_DL, or Screen B, and does
NOT touch the corner term. It only re-examines this one Omega window at
a much finer step, at the SAME persist basis 2455568 already used, and
reports the raw sigma profile plus the standard dip detector's verdict
on that finer grid.

PRE-REGISTERED
--------------
  REOPEN (this MISS deserves a second look before it goes in the paper)
  if, in window [1.96249, 2.00249] at persist basis (6,3):
    - the grid GLOBAL minimum sigma is below 1e-3 (two decades under
      SIGMA_LIST), AND
    - it sits within 1% of FE=1.98249, AND
    - it is flanked by sigma at least 10x larger at the sample points
      +/-0.001 away on both sides (a real notch, not a single noisy
      point).
  MISS CONFIRMED if no such point exists, or if the small sigma is not
  flanked by a genuine rise on both sides (i.e. still no isolated dip
  even at 0.0002 step).
  Either way this probe reports both the local_minima()-detector verdict
  (same algorithm as 2455568) AND the raw global-min/profile so a human
  can judge the "confirmed" case even if the automated 3-point test
  still says NONE (a notch narrower than one step will fail that test
  even when it is real).
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
CHECKPOINT_DIR = os.environ.get("CHECKPOINT_DIR", "./p5_1982_reopen_checkpoints")
PKG = os.environ.get("PKG_PATH", ".")
LOB = 1.0
SYM = True
LFE = 1.98249
LO, HI, STEP = 1.96249, 2.00249, 0.0002
IM_CAP = 7.0


def fine_scan(asm):
    path = os.path.join(CHECKPOINT_DIR, "fine_lob1p0_SYM_persist.json")
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
    print("probe_rect_ff_oop_1982_reopen_2026-09-05")
    print("start", time.strftime("%Y-%m-%d %H:%M:%S"), f"dps={mp.dps}")
    if not L.deploy_ok(PKG, EXPECT_VER):
        raise SystemExit(2)
    os.makedirs(CHECKPOINT_DIR, exist_ok=True)
    asm = L.make_ff()

    print(f"\n=== fine reopen l/b={LOB} SYM FE={LFE:.5f} "
          f"window [{LO},{HI}] step {STEP} ===", flush=True)
    res = fine_scan(asm)
    print(f"  n_real={res['n_real']} n_cpair={res['n_cpair']} "
          f"{len(res['rows'])} pts {res['elapsed_s']:.1f}s", flush=True)

    finite = [(Lam, s, n) for Lam, s, n in res["rows"] if s is not None]
    if not finite:
        print("  ALL NONE -- solver returned no branches in this window",
              flush=True)
        print("=== SUMMARY MISS CONFIRMED (no data) ===", flush=True)
        return

    Lam_min, s_min, n_min = min(finite, key=lambda t: t[1])
    print(f"  GLOBAL MIN  Omega={Lam_min:.4f}  sigma={s_min:.3e}  n={n_min}",
          flush=True)

    print("  dip detector (local_minima, same algorithm as 2455568): "
          + (", ".join(f"{a:.4f}({b:.2e})" for a, b, _c in res["dips"])
             or "NONE"), flush=True)

    # profile around the global min, +/- 5 steps
    idx = {round(Lam, 6): (s, n) for Lam, s, n in finite}
    keys = sorted(idx.keys())
    gi = keys.index(round(Lam_min, 6))
    print("  profile around global min:", flush=True)
    for k in keys[max(0, gi - 5):gi + 6]:
        s, n = idx[k]
        marker = " <== min" if abs(k - Lam_min) < 1e-9 else ""
        print(f"    Omega={k:.4f}  sigma={s:.3e}  n={n}{marker}", flush=True)

    # flanked-by-rise check at +/- STEP (one grid step each side)
    left_k = round(Lam_min - STEP, 6)
    right_k = round(Lam_min + STEP, 6)
    left = idx.get(left_k)
    right = idx.get(right_k)
    flanked = (left is not None and right is not None
               and left[0] >= 10 * s_min and right[0] >= 10 * s_min)

    miss_pct = 100.0 * abs(Lam_min - LFE) / LFE
    reopen = (s_min < 1e-3 and miss_pct <= 1.0 and flanked)

    print(f"\n  miss%={miss_pct:.3f}  sigma_min={s_min:.3e}  flanked={flanked}",
          flush=True)
    if reopen:
        print("=== SUMMARY REOPEN -- genuine notch found near FE, "
              "MISS verdict from 2455568 should be revisited ===",
              flush=True)
    else:
        print("=== SUMMARY MISS CONFIRMED (no qualifying notch even at "
              f"step {STEP}) ===", flush=True)
    print("end", time.strftime("%Y-%m-%d %H:%M:%S"), flush=True)


if __name__ == "__main__":
    t0 = time.time()
    main()
    print(f"wall time: {time.time() - t0:.1f}s", flush=True)
