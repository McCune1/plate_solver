#!/usr/bin/env python3
"""
Rectangular FFFF in-plane discovery at l/b=1.0/1.5/2.0/2.5/3.0,
Omega in [0.02, 2.50], pair to IP_FE_LISTS from ANSYS job 2455569.

Production basis n_real=5 n_cpair=0 (cantilever IP default).
Exploratory persist trial: n_cpair=3, |dOm|<=0.01. This is NOT Screen B
and is not adopted from this job.

PRE-REGISTERED
--------------
  MATCH if a persist-trial dip is within 5% of a same-parity FE Omega
  (coarse grid, same cut as OOP high-Lambda 2453875).
  UNMATCHED persist dips are trial false positives, not paper modes.
  Empty FE windows are misses (basis or screen), not a reason to add a
  Kirchhoff IP corner (that term does not exist).
  Do not retune |dOm|<=0.01. Do not call this Screen B. Do not bump
  SOLVER_VERSION. MATCH here is a landscape, not table promotion.
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
CHECKPOINT_DIR = os.environ.get("CHECKPOINT_DIR", "./p5_ip_disc_checkpoints")
PKG = os.environ.get("PKG_PATH", ".")
DISC_LO, DISC_HI = 0.02, 2.50
MATCH_CUT = 5.0
GEOMS = (1.0, 1.5, 2.0, 2.5, 3.0)


def scan_disc(asm, lob, sym):
    sl = "SYM" if sym else "ANTI"
    path = os.path.join(CHECKPOINT_DIR, f"disc_ip_lob{lob}_{sl}.json")
    if os.path.exists(path):
        print(f"  SKIP {path}", flush=True)
        with open(path) as f:
            return json.load(f)
    n_real, n_cpair = L.n_basis_ip(persist=False)
    step = 0.01 if sym else 0.02
    t0 = time.time()
    rows, dips = L.scan_window_ip(asm, lob, sym, DISC_LO, DISC_HI, step, 7.0,
                                  n_real, n_cpair)
    listed = [(d[0], float(d[1]), d[2]) for d in dips if d[1] < L.SIGMA_LIST]
    dt = time.time() - t0
    print(f"  [disc IP l/b={lob} {sl}] {len(rows)} pts {dt:.1f}s  dips: "
          + (", ".join(f"{a:.4f}({b:.2e})" for a, b, _c in listed) or "NONE"),
          flush=True)
    out = {"dips": [[a, b, c] for a, b, c in listed], "elapsed_s": dt}
    os.makedirs(CHECKPOINT_DIR, exist_ok=True)
    with open(path, "w") as f:
        json.dump(out, f)
    return out


def main():
    print("probe_rect_ff_ip_disc_2026-09-05")
    print("start", time.strftime("%Y-%m-%d %H:%M:%S"), f"dps={mp.dps}")
    if not L.deploy_ok_ip(PKG, EXPECT_VER):
        raise SystemExit(2)
    os.makedirs(CHECKPOINT_DIR, exist_ok=True)
    asm = L.make_ff_ip()
    rows_out = []
    print("\n=== IP discovery + persist trial (NOT Screen B) ===", flush=True)
    for lob in GEOMS:
        for sym in (True, False):
            sl = "SYM" if sym else "ANTI"
            disc = scan_disc(asm, lob, sym)
            for O0, s0, n0 in disc["dips"]:
                pr = L.persist_one_ip(asm, lob, sym, O0, CHECKPOINT_DIR)
                Ofe, miss = L.nearest_ip_fe(lob, sym, O0)
                tag = "FAIL_P"
                if pr["persist"]:
                    tag = "PASS"
                    if miss is not None and miss <= MATCH_CUT:
                        tag = "MATCH"
                print(f"  l/b={lob} {sl} O={O0:.4f} sig={s0:.2e}  "
                      f"P={pr['persist']}  {tag}  FE={Ofe} miss="
                      f"{None if miss is None else '%.2f%%' % miss}",
                      flush=True)
                rows_out.append({
                    "lob": lob, "sym": sl, "O0": O0, "sigma0": s0,
                    "persist": pr["persist"], "tag": tag,
                    "Ofe": Ofe, "miss": miss,
                })

    print("\n=== MATCHES in [0.02, 2.50] (persist trial, miss<=5%) ===",
          flush=True)
    matches = [r for r in rows_out if r["tag"] == "MATCH"]
    collapsed = {}
    for r in matches:
        key = (r["lob"], r["sym"], round(r["Ofe"], 5) if r["Ofe"] else None)
        if key not in collapsed or r["miss"] < collapsed[key]["miss"]:
            collapsed[key] = r
    for r in collapsed.values():
        print(f"  l/b={r['lob']} {r['sym']}  O={r['O0']:.4f}  "
              f"FE={r['Ofe']:.5f}  miss={r['miss']:.2f}%")
    print(f"  unique MATCH count: {len(collapsed)}")

    print("\n=== persist-PASS but UNMATCHED ===", flush=True)
    for r in rows_out:
        if r["persist"] and r["tag"] != "MATCH":
            print(f"  l/b={r['lob']} {r['sym']}  O={r['O0']:.4f}  "
                  f"nearest FE={r['Ofe']} miss={r['miss']:.1f}%")

    print("\n=== FE in [0.02, 2.50] with no unique MATCH ===", flush=True)
    claimed = set()
    for r in collapsed.values():
        claimed.add((r["lob"], r["sym"], round(r["Ofe"], 5)))
    for (lob, sym), lst in L.IP_FE_LISTS.items():
        sl = "SYM" if sym else "ANTI"
        for Ofe in lst:
            if Ofe < DISC_LO - 1e-6 or Ofe > DISC_HI + 1e-6:
                continue
            if (lob, sl, round(Ofe, 5)) not in claimed:
                print(f"  l/b={lob} {sl}  FE={Ofe:.5f}  UNCOVERED")
    print("end", time.strftime("%Y-%m-%d %H:%M:%S"), flush=True)


if __name__ == "__main__":
    t0 = time.time()
    main()
    print(f"wall time: {time.time() - t0:.1f}s", flush=True)
