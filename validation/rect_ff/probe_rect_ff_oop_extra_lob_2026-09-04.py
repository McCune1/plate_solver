#!/usr/bin/env python3
"""
Overnight job 3/3: extra aspect ratios l/b=2.0 and 3.0, same discovery+B
pipeline as 2453559, Lambda in [0.02, 2.40].

No FE decks exist at these l/b. Exploratory: does screen B still kill
ANTI 0.04/0.46 and leak SYM 0.41 (the family should be l/b-independent
if it is a basis artifact)?

PRE-REGISTERED
--------------
  G2-family: ANTI 0.04 and 0.46 (within 0.015) must FAIL B (same as
  2453559). If they PASS here, B did not travel off the FE geometries.
  G3-family: SYM 0.41 is expected to PASS as KNOWN_LEAK.
  PAPER_CANDIDATES = B-pass, not 0.41, not 0.864 birth artifact.
  Do not treat PAPER_CANDIDATES as validated modes (no FE). They are
  a consistency check of the screen and a seed list if those decks
  are built later.
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
CHECKPOINT_DIR = os.environ.get("CHECKPOINT_DIR", "./p5_extra_lob_checkpoints")
PKG = os.environ.get("PKG_PATH", ".")
DISC_LO, DISC_HI = 0.02, 2.40
GEOMS = (2.0, 3.0)


def is_leak_041(sym, lam):
    return bool(sym) and abs(lam - 0.410) <= 0.015


def is_a04(sym, lam):
    return (not sym) and abs(lam - 0.040) <= 0.015


def is_a46(sym, lam):
    return (not sym) and abs(lam - 0.460) <= 0.015


def is_birth(sym, lam):
    return bool(sym) and 0.862 <= lam <= 0.866


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
          + (", ".join(f"{a:.4f}({b:.2e})" for a, b, _c in listed) or "NONE"),
          flush=True)
    out = {"dips": [[a, b, c] for a, b, c in listed], "elapsed_s": dt}
    os.makedirs(CHECKPOINT_DIR, exist_ok=True)
    with open(path, "w") as f:
        json.dump(out, f)
    return out


def main():
    print("probe_rect_ff_oop_extra_lob_2026-09-04")
    print("start", time.strftime("%Y-%m-%d %H:%M:%S"), f"dps={mp.dps}")
    if not L.deploy_ok(PKG, EXPECT_VER):
        raise SystemExit(2)
    os.makedirs(CHECKPOINT_DIR, exist_ok=True)
    asm = L.make_ff()
    rows_out = []
    print("\n=== l/b=2.0 and 3.0 discovery + B ===", flush=True)
    for lob in GEOMS:
        for sym in (True, False):
            sl = "SYM" if sym else "ANTI"
            disc = scan_disc(asm, lob, sym)
            for L0, s0, n0 in disc["dips"]:
                pr = L.persist_one(asm, lob, sym, L0, CHECKPOINT_DIR)
                leak = is_leak_041(sym, L0)
                birth = is_birth(sym, L0)
                label = "FAIL_B"
                if birth:
                    label = "BIRTH_ARTIFACT"
                elif leak and pr["persist"]:
                    label = "KNOWN_LEAK"
                elif pr["persist"]:
                    label = "PASS"
                print(f"  l/b={lob} {sl} L={L0:.4f} sig={s0:.2e}  "
                      f"B={pr['persist']}  {label}", flush=True)
                rows_out.append({
                    "lob": lob, "sym": sl, "L0": L0, "sigma0": s0,
                    "persist": pr["persist"], "label": label,
                    "a04": is_a04(sym, L0), "a46": is_a46(sym, L0),
                })

    a04p = [r for r in rows_out if r["a04"] and r["persist"]]
    a46p = [r for r in rows_out if r["a46"] and r["persist"]]
    leak = [r for r in rows_out if r["label"] == "KNOWN_LEAK"]
    print("\n=== Family checks ===", flush=True)
    print(f"  ANTI 0.04 B-pass={len(a04p)} (expect 0)")
    print(f"  ANTI 0.46 B-pass={len(a46p)} (expect 0)")
    print(f"  SYM 0.41 KNOWN_LEAK={len(leak)} (expect 2, one per l/b)")
    g2 = len(a04p) == 0 and len(a46p) == 0
    print(f"  G2 travel: {'PASS' if g2 else 'FAIL -- B did not travel'}")

    print("\n=== PAPER_CANDIDATES (no FE -- not validated) ===", flush=True)
    for r in rows_out:
        if r["label"] != "PASS":
            continue
        print(f"  l/b={r['lob']} {r['sym']}  L={r['L0']:.4f}  sig={r['sigma0']:.2e}")
    print("end", time.strftime("%Y-%m-%d %H:%M:%S"), flush=True)


if __name__ == "__main__":
    t0 = time.time()
    main()
    print(f"wall time: {time.time() - t0:.1f}s", flush=True)
