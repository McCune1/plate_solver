#!/usr/bin/env python3
"""
Overnight job 2/3: Lambda in [2.40, 6.50] at l/b=1.0/1.5/2.5, screen B,
pair B-pass dips to the FULL FE list.

2453559 stopped at 2.40. FE has more modes in (2.4, 6.5): e.g. l/b=1.0
SYM 2.454/3.493/6.160, l/b=1.5 SYM 2.588/3.004/4.404, l/b=2.5 SYM
2.467/3.164/3.195/4.138/...

PRE-REGISTERED
--------------
  Discovery production basis, sigma<0.3 dips, then screen B (same 0.01
  persist cut, do not retune).
  MATCH if B-pass dip is within 5% of some FE Lambda (coarse grid, same
  cut as 2453559).
  SYM 0.41 cannot appear in this window. ANTI 0.04/0.46 cannot either.
  A MATCH here is a higher-mode confirmation, not a formula change.
  Empty windows near a listed FE mode are misses (basis or screen), not
  a reason to revert the corner term.
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
CHECKPOINT_DIR = os.environ.get("CHECKPOINT_DIR", "./p5_highlam_checkpoints")
PKG = os.environ.get("PKG_PATH", ".")
DISC_LO, DISC_HI = 2.40, 6.50
MATCH_CUT = 5.0


def scan_disc(asm, lob, sym):
    sl = "SYM" if sym else "ANTI"
    path = os.path.join(CHECKPOINT_DIR, f"disc_lob{lob}_{sl}.json")
    if os.path.exists(path):
        print(f"  SKIP {path}", flush=True)
        with open(path) as f:
            return json.load(f)
    n_real, n_cpair = L.n_basis(sym, persist=False)
    step = 0.01 if sym else 0.02
    t0 = time.time()
    rows, dips = L.scan_window(asm, lob, sym, DISC_LO, DISC_HI, step, 7.0,
                               n_real, n_cpair)
    listed = [(d[0], float(d[1]), d[2]) for d in dips if d[1] < L.SIGMA_LIST]
    dt = time.time() - t0
    print(f"  [disc l/b={lob} {sl}] {len(rows)} pts {dt:.1f}s  dips: "
          + (", ".join(f"{L0:.4f}({s:.2e})" for L0, s, _n in listed) or "NONE"),
          flush=True)
    out = {"dips": [[a, b, c] for a, b, c in listed], "elapsed_s": dt}
    os.makedirs(CHECKPOINT_DIR, exist_ok=True)
    with open(path, "w") as f:
        json.dump(out, f)
    return out


def main():
    print("probe_rect_ff_oop_highlam_2026-09-04")
    print("start", time.strftime("%Y-%m-%d %H:%M:%S"), f"dps={mp.dps}")
    if not L.deploy_ok(PKG, EXPECT_VER):
        raise SystemExit(2)
    os.makedirs(CHECKPOINT_DIR, exist_ok=True)
    asm = L.make_ff()
    rows_out = []
    print("\n=== High-Lambda discovery + B ===", flush=True)
    for lob in (1.0, 1.5, 2.5):
        for sym in (True, False):
            sl = "SYM" if sym else "ANTI"
            disc = scan_disc(asm, lob, sym)
            for L0, s0, n0 in disc["dips"]:
                pr = L.persist_one(asm, lob, sym, L0, CHECKPOINT_DIR)
                Lfe, miss = L.nearest_fe(lob, sym, L0)
                tag = "FAIL_B"
                if pr["persist"]:
                    tag = "PASS"
                    if miss is not None and miss <= MATCH_CUT:
                        tag = "MATCH"
                print(f"  l/b={lob} {sl} L={L0:.4f} sig={s0:.2e}  "
                      f"B={pr['persist']}  {tag}  FE={Lfe} miss="
                      f"{None if miss is None else '%.2f%%' % miss}",
                      flush=True)
                rows_out.append({
                    "lob": lob, "sym": sl, "L0": L0, "sigma0": s0,
                    "persist": pr["persist"], "tag": tag,
                    "Lfe": Lfe, "miss": miss,
                })

    print("\n=== MATCHES in [2.40, 6.50] ===", flush=True)
    n = 0
    for r in rows_out:
        if r["tag"] != "MATCH":
            continue
        n += 1
        print(f"  l/b={r['lob']} {r['sym']}  L={r['L0']:.4f}  "
              f"FE={r['Lfe']:.5f}  miss={r['miss']:.2f}%")
    print(f"  MATCH count: {n}")
    print("\n=== B-pass but UNMATCHED ===", flush=True)
    for r in rows_out:
        if r["persist"] and r["tag"] != "MATCH":
            print(f"  l/b={r['lob']} {r['sym']}  L={r['L0']:.4f}  "
                  f"nearest FE={r['Lfe']} miss={r['miss']:.1f}%")
    print("end", time.strftime("%Y-%m-%d %H:%M:%S"), flush=True)


if __name__ == "__main__":
    t0 = time.time()
    main()
    print(f"wall time: {time.time() - t0:.1f}s", flush=True)
