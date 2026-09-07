#!/usr/bin/env python3
"""
Overnight job 1/3: refine 2453559 PAPER_CANDIDATES and pair to the FULL
FE lists (not just the 8 fundamentals).

PRE-REGISTERED
--------------
  MATCH if the refined production-basis local min is within 3% of some
  mesh-converged FE Lambda in FE_LISTS (full non-rigid list, KFAC=24.6644).
  The 3% cut is tighter than 2453559's 5% because this job refines.
  UNMATCHED extras are leftover B false-positives, not paper modes.
  Double-dips within 0.04 of the same FE count as ONE match (keep closer).
  l/b=1.0 SYM 2.08 vs FE 1.98249: if refined miss > 3%, drop the 2453559
  "MATCH" and keep the known-miss call from 2453179.
  Do not retune screen B. Do not reopen the corner term.
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
CHECKPOINT_DIR = os.environ.get("CHECKPOINT_DIR", "./p5_refine_pair_checkpoints")
PKG = os.environ.get("PKG_PATH", ".")
MATCH_CUT = 3.0


def refine_one(asm, lob, sym, L0):
    sl = "SYM" if sym else "ANTI"
    path = os.path.join(
        CHECKPOINT_DIR,
        "ref_lob%s_%s_L%s.json" % (lob, sl, ("%.4f" % L0).replace(".", "p")),
    )
    if os.path.exists(path):
        with open(path) as f:
            return json.load(f)
    n_real, n_cpair = L.n_basis(sym, persist=False)
    lo = max(0.02, round(L0 - 0.04, 6))
    hi = round(L0 + 0.04, 6)
    rows, dips = L.scan_window(asm, lob, sym, lo, hi, 0.002, 7.0, n_real, n_cpair)
    if dips:
        best = min(dips, key=lambda t: abs(t[0] - L0))
        Lstar, sstar, nstar = best[0], float(best[1]), best[2]
    else:
        # no interior min; take grid min
        finite = [(Lam, s, n) for Lam, s, n in rows if s is not None]
        if not finite:
            Lstar, sstar, nstar = L0, None, 0
        else:
            Lam, s, n = min(finite, key=lambda t: t[1])
            Lstar, sstar, nstar = Lam, float(s), n
    out = {
        "L0": L0, "Lstar": Lstar, "sigma": sstar, "n": nstar,
        "dips": [[d[0], float(d[1]), d[2]] for d in dips],
    }
    os.makedirs(CHECKPOINT_DIR, exist_ok=True)
    with open(path, "w") as f:
        json.dump(out, f)
    return out


def main():
    print("probe_rect_ff_oop_refine_pair_2026-09-04")
    print("start", time.strftime("%Y-%m-%d %H:%M:%S"), f"dps={mp.dps}")
    if not L.deploy_ok(PKG, EXPECT_VER):
        raise SystemExit(2)
    os.makedirs(CHECKPOINT_DIR, exist_ok=True)
    asm = L.make_ff()
    results = []
    print("\n=== Refine + full-FE pair ===", flush=True)
    print(f"  {'l/b':>5s} {'cls':>4s} {'L0':>8s} {'L*':>8s} {'sig*':>10s} "
          f"{'FE':>8s} {'miss':>8s}  tag", flush=True)
    for lob, sym, L0, _s0 in L.PAPER_2453559:
        sl = "SYM" if sym else "ANTI"
        ref = refine_one(asm, lob, sym, L0)
        pr = L.persist_one(asm, lob, sym, ref["Lstar"], CHECKPOINT_DIR)
        Lfe, miss = L.nearest_fe(lob, sym, ref["Lstar"])
        tag = "UNMATCHED"
        if miss is not None and miss <= MATCH_CUT:
            tag = "MATCH"
        rec = {
            "lob": lob, "sym": sl, "L0": L0, "Lstar": ref["Lstar"],
            "sigma": ref["sigma"], "persist": pr["persist"],
            "Lfe": Lfe, "miss": miss, "tag": tag,
        }
        results.append(rec)
        fe = "" if Lfe is None else f"{Lfe:8.5f}"
        ms = "" if miss is None else f"{miss:7.2f}%"
        print(f"  {lob:5.1f} {sl:>4s} {L0:8.4f} {ref['Lstar']:8.4f} "
              f"{(ref['sigma'] or -1):10.2e} {fe:>8s} {ms:>8s}  {tag}  "
              f"B@L*={pr['persist']}", flush=True)

    print("\n=== MATCHES (refined, miss<=3%, persist still True preferred) ===",
          flush=True)
    matches = [r for r in results if r["tag"] == "MATCH"]
    # collapse doubles: same (lob,sym, Lfe) keep smaller miss
    collapsed = {}
    for r in matches:
        key = (r["lob"], r["sym"], round(r["Lfe"], 5) if r["Lfe"] else None)
        if key not in collapsed or r["miss"] < collapsed[key]["miss"]:
            collapsed[key] = r
    for r in collapsed.values():
        print(f"  l/b={r['lob']} {r['sym']}  L*={r['Lstar']:.4f}  "
              f"FE={r['Lfe']:.5f}  miss={r['miss']:.2f}%  persist={r['persist']}")
    print(f"  unique MATCH count: {len(collapsed)}")

    print("\n=== UNMATCHED extras (B false positives at this FE coverage) ===",
          flush=True)
    for r in results:
        if r["tag"] != "UNMATCHED":
            continue
        print(f"  l/b={r['lob']} {r['sym']}  L*={r['Lstar']:.4f}  "
              f"nearest FE={r['Lfe']}  miss={r['miss']:.1f}%  persist={r['persist']}")

    print("end", time.strftime("%Y-%m-%d %H:%M:%S"), flush=True)


if __name__ == "__main__":
    t0 = time.time()
    main()
    print(f"wall time: {time.time() - t0:.1f}s", flush=True)
