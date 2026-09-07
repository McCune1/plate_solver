#!/usr/bin/env python3
"""
OPTIONAL cluster probe: golden-section polish of remaining Table-1 ANTI
MAC misses, then dump (lob, Lam*, sigma, n_used) for off-cluster MAC.

Does NOT retune Screen B. Does NOT bump SOLVER_VERSION.

These rows still fail MAC locally at dps=30:
  (1.5, ANTI, 2.124)  weakest primary, Leissa (3,2)
  (2.0, ANTI, 1.482)
  (2.5, ANTI, 0.534)  1st ANTI, MAC saturates at 0.655
  (2.5, ANTI, 1.148)
  (3.0, ANTI, 0.936)
  (3.0, ANTI, 1.528)
  (3.0, ANTI, 2.258)

PRE-REGISTERED
  For each seed, persist basis (n_real=5, n_cpair=3, im_cap=30),
  golden-section minimise sigma on [L0-0.02, L0+0.02].
  A. Polished Lam* changes MAC by <0.02 vs the local dps=30 number
     -> coarseness is not the issue; stop.
  B. A row crosses 0.744 at the polished root
     -> report CONFIRMED_POLISH, still do not retune Screen B.
  C. sigma has no interior minimum
     -> say so.

USAGE (cluster)
  python3 -u probe_rect_ff_oop_mac_anti_remain_2026-09-06.py
Off-cluster MAC of the polished Lam* is a separate local rerun of
probe_rect_ff_oop_mac_2026-09-05.py with ONLY= set from this JSON.
"""
from __future__ import annotations
import json, os, sys, time
sys.path.insert(0, os.environ.get("PKG_PATH", "."))
os.environ.setdefault("DPS", "40")
from mpmath import mp, mpf
mp.dps = int(os.environ["DPS"])
import p5_rect_ff_lib as L
import plate_solver as ps

EXPECT_VER = os.environ.get("EXPECT_SOLVER_VERSION", "2026-07-10.s10")
SEEDS = [
    (1.5, False, 2.124),
    (2.0, False, 1.482),
    (2.5, False, 0.534),
    (2.5, False, 1.148),
    (3.0, False, 0.936),
    (3.0, False, 1.528),
    (3.0, False, 2.258),
]


def golden(asm, lob, sym, lo, hi, n_real, n_cpair, im_cap, niter=25):
    phi = (5.0 ** 0.5 - 1.0) / 2.0
    a, b = lo, hi

    def f(x):
        s, n = L.sigma_at(asm, sym, lob, x, n_real, n_cpair, im_cap)
        return s, n

    c = b - phi * (b - a)
    d = a + phi * (b - a)
    fc, nc = f(c)
    fd, nd = f(d)
    hist = [(c, fc, nc), (d, fd, nd)]
    for _ in range(niter):
        if fc is None or fd is None:
            break
        if fc < fd:
            b, d, fd, nd = d, c, fc, nc
            c = b - phi * (b - a)
            fc, nc = f(c)
            hist.append((c, fc, nc))
        else:
            a, c, fc, nc = c, d, fd, nd
            d = a + phi * (b - a)
            fd, nd = f(d)
            hist.append((d, fd, nd))
    valid = [(x, s, n) for x, s, n in hist if s is not None]
    if not valid:
        return None
    x, s, n = min(valid, key=lambda t: t[1])
    return dict(Lam=x, sigma=float(s), n_used=n, lo=lo, hi=hi)


def main():
    print("probe_rect_ff_oop_mac_anti_remain start", time.strftime("%Y-%m-%d %H:%M:%S"), flush=True)
    pkg = os.environ.get("PKG_PATH", ".")
    if not L.deploy_ok(pkg, EXPECT_VER):
        raise SystemExit(2)
    print("preflight OK", ps.SOLVER_VERSION, "dps", mp.dps, flush=True)
    asm = L.make_ff()
    n_real, n_cpair = L.n_basis(False, persist=True)
    out = []
    for lob, sym, L0 in SEEDS:
        sl = "SYM" if sym else "ANTI"
        lo, hi = max(0.02, L0 - 0.02), L0 + 0.02
        print(f"\n{lob} {sl} L0={L0} window [{lo:.4f},{hi:.4f}]", flush=True)
        t0 = time.time()
        rec = golden(asm, lob, sym, lo, hi, n_real, n_cpair, 30.0)
        if rec is None:
            print("  NO MIN", flush=True)
            out.append(dict(lob=lob, sym=sym, L0=L0, ok=False))
            continue
        rec.update(lob=lob, sym=sym, L0=L0, ok=True, dt=time.time() - t0)
        print(f"  Lam*={rec['Lam']:.6f} sigma={rec['sigma']:.3e} n={rec['n_used']} "
              f"[{rec['dt']:.1f}s]", flush=True)
        out.append(rec)
    path = os.environ.get("OUT_JSON", "rect_ff_oop_mac_anti_remain.json")
    json.dump(dict(solver_version=ps.SOLVER_VERSION, dps=mp.dps, results=out),
              open(path, "w"), indent=2)
    print("wrote", path, flush=True)


if __name__ == "__main__":
    main()
