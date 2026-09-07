#!/usr/bin/env python3
"""
Enlarged-basis rescan of the Omega in [1.20, 1.50] band, all five l/b,
both parities, IP FFFF.

Job 2455570's own discovery scan (production basis n_real=5, n_cpair=0)
found NO dip at all (not even a FAIL_P candidate -- absent from the
printed "dips:" list) anywhere in this band for most (lob,sym)
combinations, even though IP_FE_LISTS has an FE target there for nearly
every one of them:
  (1.0,SYM) 1.32984, 1.41422 | (1.0,ANTI) 1.32984
  (1.5,SYM) 1.41236, 1.42432 | (1.5,ANTI) none in-band
  (2.0,SYM) 1.40011, 1.41422, 1.44333 | (2.0,ANTI) 1.27148 (edge of band)
  (2.5,SYM) 1.22293, 1.41225, 1.42658 | (2.5,ANTI) 1.44246
  (3.0,SYM) 1.41422, 1.41791, 1.41978 | (3.0,ANTI) 1.27024, 1.29465
This band is strikingly aspect-ratio-independent (Omega~1.3-1.45 at
every l/b), unlike most bending/extension branches, and is the
2026-09-05 CLAUDE.md banner's own named hypothesis: "likely basis
(n_real=5, n_cpair=0), not a missing corner" -- there is no IP corner
term to invoke (Seok Part 2 Eq. 1/36, confirmed 2026-09-05).

This probe rescans [1.20, 1.50] at step 0.01 for BOTH the production
basis (n_real=5, n_cpair=0, same as 2455570, for a same-window
sanity-reproduction) AND the enlarged basis (n_real=5, n_cpair=3, the
persist basis) side by side, so any dip that appears only at the larger
basis is directly attributable to n_cpair. Any dip found at n_cpair=3
is confirmed with the standard persist_one_ip (n_cpair=3, +/-0.02
window) exactly as 2455570 did.

PRE-REGISTERED
--------------
  BASIS-STARVATION CONFIRMED for a given (lob,sym) if the n_cpair=3
  rescan shows a dip with sigma < SIGMA_LIST=0.3 near an UNCOVERED FE
  target listed above, where the n_cpair=0 rescan in the same window
  shows no such dip (reproducing 2455570's own silence).
  GENUINE GAP (n_cpair=3 still not enough) if no such dip appears even
  at n_cpair=3 -- escalate to n_cpair=5 or reconsider whether this mode
  family is representable at all, rather than assume basis alone
  explains it.
  Do not call this Screen B. Do not adopt n_cpair=3 as a new discovery
  default without further validation. Do not touch SIGMA_LIST or
  PERSIST_DL. Do not reopen the IP-corner question.
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
CHECKPOINT_DIR = os.environ.get("CHECKPOINT_DIR", "./p5_ip_gap_checkpoints")
PKG = os.environ.get("PKG_PATH", ".")
LO, HI, STEP = 1.20, 1.50, 0.01
IM_CAP = 7.0
GEOMS = (1.0, 1.5, 2.0, 2.5, 3.0)


def scan_gap(asm, lob, sym, n_cpair, tag):
    sl = "SYM" if sym else "ANTI"
    path = os.path.join(CHECKPOINT_DIR, f"gap_lob{lob}_{sl}_{tag}.json")
    if os.path.exists(path):
        print(f"    SKIP {path}", flush=True)
        with open(path) as f:
            return json.load(f)
    t0 = time.time()
    rows, dips = L.scan_window_ip(asm, lob, sym, LO, HI, STEP, IM_CAP,
                                   5, n_cpair)
    dt = time.time() - t0
    listed = [(d[0], float(d[1]), d[2]) for d in dips if d[1] < L.SIGMA_LIST]
    out = {"dips": [[a, b, c] for a, b, c in listed], "elapsed_s": dt}
    os.makedirs(CHECKPOINT_DIR, exist_ok=True)
    with open(path, "w") as f:
        json.dump(out, f)
    return out


def main():
    print("probe_rect_ff_ip_disc_gap_enlarged_2026-09-05")
    print("start", time.strftime("%Y-%m-%d %H:%M:%S"), f"dps={mp.dps}")
    if not L.deploy_ok_ip(PKG, EXPECT_VER):
        raise SystemExit(2)
    os.makedirs(CHECKPOINT_DIR, exist_ok=True)
    asm = L.make_ff_ip()

    print(f"\n=== gap rescan Omega in [{LO},{HI}] step {STEP} ===", flush=True)
    confirmed = []
    for lob in GEOMS:
        for sym in (True, False):
            sl = "SYM" if sym else "ANTI"
            base = scan_gap(asm, lob, sym, 0, "ncp0")
            big = scan_gap(asm, lob, sym, 3, "ncp3")
            print(f"  [l/b={lob} {sl}]  ncp0 {base['elapsed_s']:.1f}s dips: "
                  + (", ".join(f"{a:.4f}({b:.2e})" for a, b, _c in base["dips"])
                     or "NONE"), flush=True)
            print(f"  [l/b={lob} {sl}]  ncp3 {big['elapsed_s']:.1f}s dips: "
                  + (", ".join(f"{a:.4f}({b:.2e})" for a, b, _c in big["dips"])
                     or "NONE"), flush=True)
            for O0, s0, _n0 in big["dips"]:
                pr = L.persist_one_ip(asm, lob, sym, O0, CHECKPOINT_DIR)
                Ofe, miss = L.nearest_ip_fe(lob, sym, O0)
                new_at_ncp3_only = not any(
                    abs(O0 - a) <= 0.02 for a, _b, _c in base["dips"])
                print(f"    ncp3 dip O={O0:.4f} sig={s0:.2e} persist="
                      f"{pr['persist']} FE={Ofe} miss="
                      f"{None if miss is None else '%.2f%%' % miss} "
                      f"NEW_AT_NCP3={new_at_ncp3_only}", flush=True)
                if new_at_ncp3_only and miss is not None and miss <= 5.0:
                    confirmed.append((lob, sl, O0, Ofe, miss, pr["persist"]))

    print("\n=== BASIS-STARVATION CONFIRMED candidates "
          "(new at ncp3, within 5% of an FE target) ===", flush=True)
    if not confirmed:
        print("  NONE -- gap not explained by n_cpair alone in [1.20,1.50]",
              flush=True)
    for lob, sl, O0, Ofe, miss, persist in confirmed:
        print(f"  l/b={lob} {sl}  O={O0:.4f}  FE={Ofe:.5f}  "
              f"miss={miss:.2f}%  persist={persist}")
    print("end", time.strftime("%Y-%m-%d %H:%M:%S"), flush=True)


if __name__ == "__main__":
    t0 = time.time()
    main()
    print(f"wall time: {time.time() - t0:.1f}s", flush=True)
