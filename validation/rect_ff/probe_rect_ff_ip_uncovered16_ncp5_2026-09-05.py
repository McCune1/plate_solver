#!/usr/bin/env python3
"""
Targeted mop-up: job 2456002 (full [0.02,2.50] IP re-discovery directly at
n_cpair=3, 3% cut) raised the MATCH count from 38 to 76 unique entries
(40 NEW, 36 of the original 38 reconfirmed) but still left 16 FE targets
UNCOVERED. Two of those 16 are former OLD MATCHes that did NOT reproduce
this time -- (2.0, ANTI, 1.27148) and (3.0, SYM, 1.66346) -- both sitting
right next to a DIFFERENT FE root that DID get claimed instead
((2.0,ANTI,1.28703) NEW; (3.0,SYM,1.67627) OLD, gap only ~0.013-0.015),
consistent with a closely-spaced-root crowding effect (the discovery grid
step / fixed +/-0.02 refine window locks onto the nearer neighbor) rather
than a genuine disappearance -- see LESSONS_LEARNED
`fe-near-degenerate-multiplicity-2026-07-29.md` for the same pattern
elsewhere in this project.

Rather than re-run the discovery step wholesale, this probe goes directly
at each of the 16 UNCOVERED (l/b, parity, FE) triples with the STRONGEST
basis this project has validated so far (n_cpair=5, the exact basis
2455677 already showed recovers ANTI targets that n_cpair=3 alone missed),
scanning a window centered ON THE FE VALUE ITSELF (+/-0.05, step 0.002)
rather than seeded from a possibly-crowded discovery dip. This is the
same basis/window 2455677 validated, applied for the first time to the
SYM branch and to the specific residual list a full re-discovery still
couldn't clear.

Does NOT retune SIGMA_LIST, PERSIST_DL, or MATCH_CUT. Does NOT touch the
corner term (none exists for IP). Does NOT bump SOLVER_VERSION.

PRE-REGISTERED
--------------
  RECOVERED if a local_minima() dip inside [FE-0.05, FE+0.05] lands within
  3% of that FE (the paper's adopted cut).
  STILL UNCOVERED if no such dip exists even at n_cpair=5 in this window
  -- at that point basis/window are no longer the leading suspects (both
  already pushed further than anywhere else in this project's IP work),
  and the honest write-up is "not found," not "must exist, keep looking."
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
CHECKPOINT_DIR = os.environ.get("CHECKPOINT_DIR", "./p5_ip_uncov16_checkpoints")
PKG = os.environ.get("PKG_PATH", ".")
N_REAL, N_CPAIR = 5, 5
WINDOW, STEP, IM_CAP = 0.05, 0.002, 7.0
MATCH_CUT = 3.0

# The 16 FE targets job 2456002 (n_cpair=3 full re-discovery, 3% cut)
# still could not claim. Independently parse-verified at 16 entries.
UNCOVERED_16 = [
    (1.0, False, 1.32984),
    (1.5, True, 1.42432),
    (1.5, False, 1.03335),
    (1.5, False, 1.57124),
    (1.5, False, 2.20253),
    (2.0, True, 1.41422),
    (2.0, False, 1.27148),
    (2.5, True, 1.42658),
    (2.5, True, 2.40468),
    (2.5, False, 1.11248),
    (3.0, True, 1.41422),
    (3.0, True, 1.41791),
    (3.0, True, 1.66346),
    (3.0, False, 1.03153),
    (3.0, False, 1.27024),
    (3.0, False, 2.49371),
]
assert len(UNCOVERED_16) == 16, len(UNCOVERED_16)


def scan_one(asm, lob, sym, fe):
    sl = "SYM" if sym else "ANTI"
    tag = f"uncov_lob{lob}_{sl}_FE{fe:.5f}".replace(".", "p")
    path = os.path.join(CHECKPOINT_DIR, f"{tag}.json")
    if os.path.exists(path):
        print(f"  SKIP {path}", flush=True)
        with open(path) as f:
            return json.load(f)
    lo = max(0.002, round(fe - WINDOW, 6))
    hi = round(fe + WINDOW, 6)
    t0 = time.time()
    rows, dips = L.scan_window_ip(asm, lob, sym, lo, hi, STEP, IM_CAP,
                                  N_REAL, N_CPAIR)
    dt = time.time() - t0
    out = {
        "dips": [[d[0], float(d[1]), d[2]] for d in dips],
        "elapsed_s": dt,
    }
    os.makedirs(CHECKPOINT_DIR, exist_ok=True)
    with open(path, "w") as f:
        json.dump(out, f)
    return out


def main():
    print("probe_rect_ff_ip_uncovered16_ncp5_2026-09-05")
    print("start", time.strftime("%Y-%m-%d %H:%M:%S"), f"dps={mp.dps}")
    if not L.deploy_ok_ip(PKG, EXPECT_VER):
        raise SystemExit(2)
    os.makedirs(CHECKPOINT_DIR, exist_ok=True)
    asm = L.make_ff_ip()

    print(f"\n=== targeted n_cpair={N_CPAIR} scan, window +/-{WINDOW} centered "
          "on each of the 16 UNCOVERED FE targets ===", flush=True)
    n_recovered = 0
    for lob, sym, fe in UNCOVERED_16:
        sl = "SYM" if sym else "ANTI"
        res = scan_one(asm, lob, sym, fe)
        dips = res["dips"]
        print(f"  [l/b={lob} {sl} FE={fe:.5f}] {res['elapsed_s']:.1f}s  "
              "dips: "
              + (", ".join(f"{a:.4f}({b:.2e})" for a, b, _c in dips)
                 or "NONE"), flush=True)
        if not dips:
            print(f"    l/b={lob} {sl} FE={fe:.5f}  STILL UNCOVERED "
                  "(no dip at all)", flush=True)
            continue
        nearest = min(dips, key=lambda d: abs(d[0] - fe))
        miss = 100.0 * abs(nearest[0] - fe) / fe
        if miss <= MATCH_CUT:
            n_recovered += 1
            print(f"    l/b={lob} {sl} FE={fe:.5f}  RECOVERED  "
                  f"O*={nearest[0]:.4f}  sig={nearest[1]:.2e}  "
                  f"miss={miss:.2f}%", flush=True)
        else:
            print(f"    l/b={lob} {sl} FE={fe:.5f}  STILL UNCOVERED  "
                  f"(nearest dip O={nearest[0]:.4f} miss={miss:.2f}%, "
                  f"over the {MATCH_CUT}% cut)", flush=True)

    print(f"\n=== SUMMARY: {n_recovered}/{len(UNCOVERED_16)} recovered at "
          f"n_cpair={N_CPAIR}, window +/-{WINDOW} ===", flush=True)
    print("end", time.strftime("%Y-%m-%d %H:%M:%S"), flush=True)


if __name__ == "__main__":
    t0 = time.time()
    main()
    print(f"wall time: {time.time() - t0:.1f}s", flush=True)
