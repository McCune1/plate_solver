#!/usr/bin/env python3
"""
Refine job 2455570's 38 persist-trial MATCHES from a 5% landscape cut
down to a paper-grade 3% cut, item 2 of PAPER2_RECT_FREEFREE_DRAFT.md's
"Remaining work" list ("optional 3% refine of the 38 MATCHES").

2455570's MATCH tag used the exploratory persist trial (n_cpair=3,
NOT Screen B) at whatever Omega the coarse discovery grid (step 0.01
SYM / 0.02 ANTI) happened to land on. This probe re-locates each of
those 38 points with a finer PRODUCTION-basis rescan (n_real=5,
n_cpair=0, same basis IP_FE_LISTS assumes and the same "refine" pattern
probe_rect_ff_oop_refine_pair_2026-09-04.py used for OOP), then
reconfirms the persist criterion (n_cpair=3, standard +/-0.02 window)
at the refined location, and reports the tightened miss% against
IP_FE_LISTS.

PRE-REGISTERED
--------------
  MATCH (tightened) if the refined production-basis Omega* is within 3%
  of the same FE target 2455570 already paired it to.
  DROP if the refined miss% exceeds 3% -- that MATCH was a 5%-cut
  artifact of the coarse discovery grid, not a paper-grade match.
  Persist (n_cpair=3) is reported at Omega* for reference but does NOT
  gate MATCH/DROP here -- persist already gated the original 2455570
  list once; this pass only tightens location and miss%.
  Do not retune Screen B. Do not touch IP_FE_LISTS. Do not call this
  Screen B.
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
CHECKPOINT_DIR = os.environ.get("CHECKPOINT_DIR", "./p5_ip_match_refine_checkpoints")
PKG = os.environ.get("PKG_PATH", ".")
MATCH_CUT = 3.0

# (lob, sym, O0) -- the 38 unique MATCH points from job 2455570's own
# "=== MATCHES in [0.02, 2.50] ===" block, verbatim.
MATCHES_38 = [
    (1.0, True, 1.8500), (1.0, True, 2.0100),
    (1.0, False, 1.2400),
    (1.5, True, 1.6300), (1.5, True, 1.6900), (1.5, True, 2.0100),
    (1.5, True, 2.3000),
    (2.0, True, 1.6500), (2.0, True, 1.7800), (2.0, True, 1.8100),
    (2.0, True, 2.0100), (2.0, True, 2.2800),
    (2.0, False, 0.5200), (2.0, False, 0.8800), (2.0, False, 1.2800),
    (2.5, True, 1.6700), (2.5, True, 1.6900), (2.5, True, 1.7800),
    (2.5, True, 2.0400), (2.5, True, 2.0700), (2.5, True, 2.1300),
    (2.5, True, 2.4000),
    (2.5, False, 0.3800), (2.5, False, 0.7200), (2.5, False, 1.5400),
    (3.0, True, 1.6600), (3.0, True, 1.6900), (3.0, True, 1.8100),
    (3.0, True, 1.8400), (3.0, True, 2.0100), (3.0, True, 2.1500),
    (3.0, True, 2.2100), (3.0, True, 2.2700), (3.0, True, 2.4200),
    (3.0, True, 2.4700),
    (3.0, False, 0.2800), (3.0, False, 0.5800), (3.0, False, 1.5800),
]
assert len(MATCHES_38) == 38, len(MATCHES_38)


def refine_one(asm, lob, sym, O0):
    sl = "SYM" if sym else "ANTI"
    path = os.path.join(
        CHECKPOINT_DIR,
        "ref_lob%s_%s_O%s.json" % (lob, sl, ("%.4f" % O0).replace(".", "p")),
    )
    if os.path.exists(path):
        with open(path) as f:
            return json.load(f)
    n_real, n_cpair = L.n_basis_ip(persist=False)
    lo = max(0.02, round(O0 - 0.02, 6))
    hi = round(O0 + 0.02, 6)
    rows, dips = L.scan_window_ip(asm, lob, sym, lo, hi, 0.001, 7.0,
                                   n_real, n_cpair)
    if dips:
        best = min(dips, key=lambda t: abs(t[0] - O0))
        Ostar, sstar, nstar = best[0], float(best[1]), best[2]
    else:
        finite = [(Om, s, n) for Om, s, n in rows if s is not None]
        if not finite:
            Ostar, sstar, nstar = O0, None, 0
        else:
            Om, s, n = min(finite, key=lambda t: t[1])
            Ostar, sstar, nstar = Om, float(s), n
    out = {
        "O0": O0, "Ostar": Ostar, "sigma": sstar, "n": nstar,
        "dips": [[d[0], float(d[1]), d[2]] for d in dips],
    }
    os.makedirs(CHECKPOINT_DIR, exist_ok=True)
    with open(path, "w") as f:
        json.dump(out, f)
    return out


def main():
    print("probe_rect_ff_ip_match_refine_2026-09-05")
    print("start", time.strftime("%Y-%m-%d %H:%M:%S"), f"dps={mp.dps}")
    if not L.deploy_ok_ip(PKG, EXPECT_VER):
        raise SystemExit(2)
    os.makedirs(CHECKPOINT_DIR, exist_ok=True)
    asm = L.make_ff_ip()
    results = []
    print("\n=== Refine + tighten to 3% ===", flush=True)
    print(f"  {'l/b':>5s} {'cls':>4s} {'O0':>8s} {'O*':>8s} {'sig*':>10s} "
          f"{'FE':>8s} {'miss':>8s}  tag", flush=True)
    for lob, sym, O0 in MATCHES_38:
        sl = "SYM" if sym else "ANTI"
        ref = refine_one(asm, lob, sym, O0)
        pr = L.persist_one_ip(asm, lob, sym, ref["Ostar"], CHECKPOINT_DIR)
        Ofe, miss = L.nearest_ip_fe(lob, sym, ref["Ostar"])
        tag = "DROP"
        if miss is not None and miss <= MATCH_CUT:
            tag = "MATCH"
        rec = {
            "lob": lob, "sym": sl, "O0": O0, "Ostar": ref["Ostar"],
            "sigma": ref["sigma"], "persist": pr["persist"],
            "Ofe": Ofe, "miss": miss, "tag": tag,
        }
        results.append(rec)
        fe = "" if Ofe is None else f"{Ofe:8.5f}"
        ms = "" if miss is None else f"{miss:7.2f}%"
        print(f"  {lob:5.1f} {sl:>4s} {O0:8.4f} {ref['Ostar']:8.4f} "
              f"{(ref['sigma'] or -1):10.2e} {fe:>8s} {ms:>8s}  {tag}  "
              f"persist@O*={pr['persist']}", flush=True)

    matched = [r for r in results if r["tag"] == "MATCH"]
    dropped = [r for r in results if r["tag"] == "DROP"]
    print(f"\n=== SUMMARY: {len(matched)}/38 hold at 3% after refine, "
          f"{len(dropped)}/38 drop ===", flush=True)
    for r in dropped:
        print(f"  DROP l/b={r['lob']} {r['sym']}  O0={r['O0']:.4f} -> "
              f"O*={r['Ostar']:.4f}  FE={r['Ofe']}  miss="
              f"{None if r['miss'] is None else '%.2f%%' % r['miss']}")
    print("end", time.strftime("%Y-%m-%d %H:%M:%S"), flush=True)


if __name__ == "__main__":
    t0 = time.time()
    main()
    print(f"wall time: {time.time() - t0:.1f}s", flush=True)
