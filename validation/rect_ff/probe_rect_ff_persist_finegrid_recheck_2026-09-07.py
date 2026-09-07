"""
probe_rect_ff_persist_finegrid_recheck_2026-09-07.py

Purpose
-------
Direct, surgical follow-up to job 2457815
(probe_rect_ff_persist_threshold_sensitivity_2026-09-07.py), which scored
all 148 of Paper 2's published production matches against the frozen
Screen-B persistence cut (PERSIST_DL=0.01) and found every one persists,
as required -- but exactly three sit AT the outer edge of what the cut
permits, each needing the full 5-grid-step allowance:

    OOP l/b=1.5 SYM  Lambda*=5.380   (Table tab:oop-high)
    OOP l/b=2.5 ANTI Lambda*=6.260   (Table tab:oop-high)
    IP  l/b=3.0 ANTI Omega*=2.2400   (Table tab:ip92)

Because persist_one/persist_one_ip's own window scan uses a fixed 0.002
step (the frozen production setting, unchanged here), every delta that
scan can report is quantized to a multiple of 0.002 -- 0, 0.002, 0.004,
0.006, 0.008, 0.01, and no value in between. "delta=0.01000" for these
three could mean either (a) the true nearest persist-basis dip really
does sit close to the edge of the +/-0.02 window's outer reach, or
(b) a closer dip exists between two of the coarse grid's sample points
and the 0.002-spaced scan simply stepped over it. Job 2457815's own
threshold sweep cannot distinguish these -- it only ever evaluates
deltas produced by the same fixed-step scan.

This probe resolves that ambiguity for exactly these three candidates,
and ONLY these three, by rescanning the identical +/-0.02 window at a
4x finer step (0.0005 instead of 0.002) using the exact, unmodified,
already-public building blocks p5_rect_ff_lib.sigma_at /.sigma_at_ip,
.local_minima, .frange, and .n_basis/.n_basis_ip(persist=True) -- the
SAME production basis and im_cap that persist_one/persist_one_ip use
internally (im_cap=30.0 OOP, im_cap=7.0 IP; persist-basis n_cpair=3 both
families). This is not a monkeypatch and does not touch persist_one/
persist_one_ip themselves -- it is the identical calculation at a finer
sample spacing, built from the same public functions the project's other
diagnostics already use this way (e.g. the branch-point scan in
probe_rect_ff_ip_symvsanti_2026-09-07.py).

Does NOT retune PERSIST_DL, does NOT edit plate_solver internals, does
NOT touch a tabulated frequency, does NOT bump SOLVER_VERSION, does NOT
edit any paper file.

Fidelity gate: for each of the three candidates, the coarse-grid (0.002
step) local minima recomputed here must reproduce job 2457815's own
result (same nearest-dip location, to within float tolerance) before the
fine-grid (0.0005 step) result is trusted.

Verdict per candidate:
  - fine_delta: nearest persist-basis dip distance at 0.0005 step.
  - resolved: "GENUINE_MARGIN" if fine_delta is still close to 0.01
    (within one coarse grid step, i.e. >= 0.008), meaning the coarse
    scan's classification was correct and this candidate really is
    close to the edge of the frozen cut; "COARSE_ARTIFACT" if fine_delta
    is meaningfully smaller (< 0.008), meaning a closer dip existed
    between coarse sample points and the true margin is more
    comfortable than job 2457815 reported.

Copy/paste on the cluster:
  cd /home/ghmkfh/PythonMill/Plate_Solver_Package
  sbatch submit_rect_ff_persist_finegrid_recheck_2026-09-07.sh

Runtime: 3 candidates x 81-point window (0.0005 step over +/-0.02) = 243
sigma_at/sigma_at_ip calls -- expect well under a minute total.
"""

from __future__ import annotations

import json
import time

from mpmath import mp

import p5_rect_ff_lib as L
import plate_solver as ps

import os

PKG = os.environ.get("PKG_PATH", ".")
EXPECT_VER = os.environ.get("EXPECT_SOLVER_VERSION", "2026-07-10.s10")
mp.dps = int(os.environ.get("DPS", "30"))

FINE_STEP = 0.0005
COARSE_STEP = 0.002
WINDOW = 0.02

CANDIDATES = [
    # (family, lob, sym, L0, coarse_delta_from_2457815)
    ("oop", 1.5, True, 5.380, 0.01000),
    ("oop", 2.5, False, 6.260, 0.01000),
    ("ip", 3.0, False, 2.2400, 0.01000),
]


def scan(asm, family, lob, sym, L0, step):
    if family == "oop":
        n_real, n_cpair = L.n_basis(sym, persist=True)
        im_cap = 30.0
        sigma_fn = L.sigma_at
    else:
        n_real, n_cpair = L.n_basis_ip(persist=True)
        im_cap = 7.0
        sigma_fn = L.sigma_at_ip
    lo = max(0.002, round(L0 - WINDOW, 6))
    hi = round(L0 + WINDOW, 6)
    grid = L.frange(lo, hi, step)
    rows = []
    for x in grid:
        if x <= 0:
            rows.append((x, None, 0))
            continue
        s, n = sigma_fn(asm, sym, lob, x, n_real, n_cpair, im_cap)
        rows.append((x, s, n))
    dips = L.local_minima(rows)
    nearest = None
    delta = None
    if dips:
        nearest = min(dips, key=lambda t: abs(t[0] - L0))
        delta = abs(nearest[0] - L0)
    return {"n_dips": len(dips), "nearest": nearest, "delta": delta,
            "all_dips": dips}


def main():
    print("probe_rect_ff_persist_finegrid_recheck_2026-09-07")
    print("start", time.strftime("%Y-%m-%d %H:%M:%S"), f"dps={mp.dps}",
          flush=True)

    if not L.deploy_ok_ip(PKG, EXPECT_VER):
        raise SystemExit(2)

    asm_oop = L.make_ff()
    asm_ip = L.make_ff_ip()

    results = []
    for family, lob, sym, L0, coarse_delta_2457815 in CANDIDATES:
        asm = asm_oop if family == "oop" else asm_ip
        sl = "SYM" if sym else "ANTI"
        print("=" * 72)
        print(f"{family} l/b={lob} {sl} L0={L0}  "
              f"(job 2457815 coarse delta={coarse_delta_2457815:.5f})")
        print("=" * 72, flush=True)

        t0 = time.time()
        coarse = scan(asm, family, lob, sym, L0, COARSE_STEP)
        print(f"  coarse (0.002 step, reproduction check): "
              f"n_dips={coarse['n_dips']} nearest={coarse['nearest']} "
              f"delta={coarse['delta']}  [{time.time()-t0:.1f}s]", flush=True)
        gate_ok = (coarse["delta"] is not None
                   and abs(coarse["delta"] - coarse_delta_2457815) < 1e-6)
        if not gate_ok:
            print(f"  FIDELITY GATE FAIL: coarse rescan does not reproduce "
                  f"job 2457815's delta ({coarse_delta_2457815}) -- "
                  f"STOP, do not trust the fine-grid result for this "
                  f"candidate.", flush=True)
            results.append({
                "family": family, "lob": lob, "sym": sl, "L0": L0,
                "gate_ok": False, "coarse": coarse, "fine": None,
                "resolved": "GATE_FAIL",
            })
            continue
        print("  gate OK -- coarse rescan reproduces job 2457815.",
              flush=True)

        t1 = time.time()
        fine = scan(asm, family, lob, sym, L0, FINE_STEP)
        print(f"  fine (0.0005 step, 4x finer): n_dips={fine['n_dips']} "
              f"nearest={fine['nearest']} delta={fine['delta']}  "
              f"[{time.time()-t1:.1f}s]", flush=True)

        if fine["delta"] is None:
            resolved = "NO_DIP_AT_FINE_GRID"
        elif fine["delta"] >= 0.008:
            resolved = "GENUINE_MARGIN"
        else:
            resolved = "COARSE_ARTIFACT"
        print(f"  VERDICT: {resolved}  (fine_delta="
              f"{fine['delta']})", flush=True)

        results.append({
            "family": family, "lob": lob, "sym": sl, "L0": L0,
            "gate_ok": True, "coarse": coarse, "fine": fine,
            "resolved": resolved,
        })

    print("\n" + "=" * 72)
    print("SUMMARY")
    print("=" * 72)
    for r in results:
        fd = None if r["fine"] is None else r["fine"]["delta"]
        fd_str = "" if fd is None else f" (fine_delta={fd:.5f})"
        print(f"  {r['family']:3s} lob={r['lob']:<4} {r['sym']:4s} "
              f"L0={r['L0']:<8} -> {r['resolved']}{fd_str}")

    with open("rect_ff_persist_finegrid_recheck_2026-09-07.json", "w") as f:
        json.dump(results, f, indent=1, default=str)
    print("\nWrote rect_ff_persist_finegrid_recheck_2026-09-07.json",
          flush=True)
    print("SOLVER_VERSION unbumped. PERSIST_DL not retuned. "
          "Paper not edited.")


if __name__ == "__main__":
    main()
