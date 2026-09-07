#!/usr/bin/env python3
"""
High-Lambda [2.40, 6.50] at l/b=2.0 and 3.0, screen B, pair B-pass dips
to the FULL FE list from ANSYS job 2454654.

2453876 stopped at 2.40. FE continues: e.g. l/b=2.0 SYM 2.622/2.998/...
l/b=2.0 ANTI 2.552, l/b=3.0 SYM 2.459/2.887/..., l/b=3.0 ANTI 3.143/...

PRE-REGISTERED
--------------
  Discovery production basis, sigma<0.3 dips, then screen B (same 0.01
  persist cut, do not retune).
  MATCH if B-pass dip is within 5% of some FE Lambda (coarse grid, same
  cut as 2453875). Collapse doubles (same FE, keep closer).
  L=4.000 or 6.000 with sigma < 1e-12 is a numerical zero artifact --
  print ZERO_ART, never MATCH, even if an FE sits nearby.
  Empty windows near a listed FE mode in [2.40, 6.50] are misses (basis
  or screen), not a reason to revert the corner term.
  SYM 0.41 / ANTI 0.04 / 0.46 cannot appear in this window.
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
CHECKPOINT_DIR = os.environ.get("CHECKPOINT_DIR", "./p5_highlam_extra_lob_checkpoints")
PKG = os.environ.get("PKG_PATH", ".")
DISC_LO, DISC_HI = 2.40, 6.50
MATCH_CUT = 5.0
GEOMS = (2.0, 3.0)


def is_zero_art(lam, sig):
    if sig is None:
        return False
    near = abs(lam - 4.0) < 0.002 or abs(lam - 6.0) < 0.002
    return near and sig < 1e-12


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
    print("probe_rect_ff_oop_highlam_extra_lob_2026-09-04")
    print("start", time.strftime("%Y-%m-%d %H:%M:%S"), f"dps={mp.dps}")
    if not L.deploy_ok(PKG, EXPECT_VER):
        raise SystemExit(2)
    os.makedirs(CHECKPOINT_DIR, exist_ok=True)
    asm = L.make_ff()
    rows_out = []
    print("\n=== High-Lambda extra-lob discovery + B ===", flush=True)
    for lob in GEOMS:
        for sym in (True, False):
            sl = "SYM" if sym else "ANTI"
            disc = scan_disc(asm, lob, sym)
            for L0, s0, n0 in disc["dips"]:
                pr = L.persist_one(asm, lob, sym, L0, CHECKPOINT_DIR)
                Lfe, miss = L.nearest_fe(lob, sym, L0)
                tag = "FAIL_B"
                if is_zero_art(L0, s0):
                    tag = "ZERO_ART"
                elif pr["persist"]:
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

    print("\n=== MATCHES in [2.40, 6.50] (zeros dropped, doubles collapsed) ===",
          flush=True)
    matches = [r for r in rows_out if r["tag"] == "MATCH"]
    collapsed = {}
    for r in matches:
        key = (r["lob"], r["sym"], round(r["Lfe"], 5) if r["Lfe"] else None)
        if key not in collapsed or r["miss"] < collapsed[key]["miss"]:
            collapsed[key] = r
    for r in collapsed.values():
        print(f"  l/b={r['lob']} {r['sym']}  L={r['L0']:.4f}  "
              f"FE={r['Lfe']:.5f}  miss={r['miss']:.2f}%")
    print(f"  unique MATCH count: {len(collapsed)}")

    print("\n=== ZERO_ART (L=4.000/6.000 machine zeros) ===", flush=True)
    n_z = 0
    for r in rows_out:
        if r["tag"] != "ZERO_ART":
            continue
        n_z += 1
        print(f"  l/b={r['lob']} {r['sym']}  L={r['L0']:.4f}  sig={r['sigma0']:.2e}")
    print(f"  ZERO_ART count: {n_z}")

    print("\n=== B-pass but UNMATCHED ===", flush=True)
    for r in rows_out:
        if r["persist"] and r["tag"] not in ("MATCH", "ZERO_ART"):
            print(f"  l/b={r['lob']} {r['sym']}  L={r['L0']:.4f}  "
                  f"nearest FE={r['Lfe']} miss={r['miss']:.1f}%")

    print("\n=== FE in [2.40, 6.50] with no unique MATCH (empty windows) ===",
          flush=True)
    claimed = set()
    for r in collapsed.values():
        claimed.add((r["lob"], r["sym"], round(r["Lfe"], 5)))
    for (lob, sym), lst in L.FE_LISTS.items():
        if lob not in GEOMS:
            continue
        sl = "SYM" if sym else "ANTI"
        for Lfe in lst:
            if Lfe < DISC_LO - 1e-6 or Lfe > DISC_HI + 1e-6:
                continue
            key = (lob, sl, round(Lfe, 5))
            if key not in claimed:
                print(f"  l/b={lob} {sl}  FE={Lfe:.5f}  UNCOVERED")
    print("end", time.strftime("%Y-%m-%d %H:%M:%S"), flush=True)


if __name__ == "__main__":
    t0 = time.time()
    main()
    print(f"wall time: {time.time() - t0:.1f}s", flush=True)
