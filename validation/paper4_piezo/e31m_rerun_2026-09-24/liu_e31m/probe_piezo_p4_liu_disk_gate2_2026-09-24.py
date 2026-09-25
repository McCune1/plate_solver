"""
probe_piezo_p4_liu_disk_gate2_2026-09-24.py

Paper 4 lever 11: Liu, Wang & Quek (IJSS 39, 2002) solid steel disk with
two PZT-4 skins, both faces of each skin short-circuited. Backs
Supplementary Sec S.10 (tab:liu-t2, tab:liu-t6). Solver:
plate_solver.piezo_disk.PiezoDiskSC (centre-regular 3x3 SC / 2x2 elastic;
NOT the annular 6x6 at r_i=0).

Roots are found FROM SCRATCH: det is sampled on a uniform omega grid
(no published value used to place a bracket), every sign change is
bisected, and the first two roots per (bc, n) are labelled m=1,2 in
frequency order. Published values enter only in the post-hoc score.

PRE-REGISTERED (fixed before this script was run):
  G2a  Elastic limit (h1=0, bare steel): lambda^2 = omega a^2 sqrt(rho H/D)
       within 1e-3 of Paper 3 / Leissa for all 12 (C,S) x n=0..2 x m=1..2.
  G2b  Table 2 (outer C) SC, projection='consistent': |omega/x_kir - 1|
       < 1e-3 for all 6 primary modes.
  G2c  Table 6 (outer S) SC, same bar, all 6 modes.
  G2d  Scan hygiene: exactly the expected roots in the scan window, i.e.
       no (bc, n) has fewer than 2 roots below its window top.
  PASS_ALL needs G2a-G2d.
Diagnostics (reported, not gated):
  D1  |consistent/duan - 1| per mode (expected <~1e-5).
  D2  coupling content: |omega(e31=4.1)/omega(e31=0.01) - 1| per mode
      (expected ~1e-5..1e-4: Liu Tables 2/6 do not test the coupling).
  D3  sign of omega - x_kir (a uniform offset is not print rounding).

Materials: Liu Table 1 (LIU_DISK_KWARGS), e31=+4.1 project sign.
Geometry: r0=0.6 m, host half-thickness h=0.01 m, skins h1=0.002 m.

Run (any machine with the package; ~2-4 min on one core):
    PKG_PATH=/path/to/github_repo python3 probe_piezo_p4_liu_disk_gate2_2026-09-24.py
Writes liu_disk_gate2_2026-09-24_results.json and .csv beside this file.
"""
import csv
import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.environ.get("PKG_PATH", os.path.normpath(os.path.join(HERE, "..", "..")))
sys.path.insert(0, PKG)

from mpmath import mp  # noqa: E402
from plate_solver.piezo_disk import make_liu_disk, LIU_DISK_KWARGS  # noqa: E402

DPS_SCAN = 30
DPS_PROD = 100
N_GRID = 800
WIN = {"C": (50.0, 8000.0), "S": (50.0, 7000.0)}
NS = (0, 1, 2)

# Liu 2002 x_kir (rad/s): Table 2 (C) and Table 6 (S), r0/h = 60.
LIU_KIR = {
    ("C", 0, 1): 902.5, ("C", 1, 1): 1878.2, ("C", 2, 1): 3081.1,
    ("C", 0, 2): 3513.4, ("C", 1, 2): 5373.7, ("C", 2, 2): 7472.1,
    ("S", 0, 1): 435.6, ("S", 1, 1): 1227.5, ("S", 2, 1): 2262.4,
    ("S", 0, 2): 2625.2, ("S", 1, 2): 4282.4, ("S", 2, 2): 6193.9,
}
# Liu 2002 3-D FEA column (ABAQUS), cited only.
LIU_FEM = {
    ("C", 0, 1): 900.1, ("C", 1, 1): 1862.7, ("C", 2, 1): 3050.9,
    ("C", 0, 2): 3475.2, ("C", 1, 2): 5272.1, ("C", 2, 2): 7306.2,
    ("S", 0, 1): 435.2, ("S", 1, 1): 1218.5, ("S", 2, 1): 2242.4,
    ("S", 0, 2): 2606.5, ("S", 1, 2): 4221.0, ("S", 2, 2): 6088.3,
}
# Paper 3 / Leissa solid-disk flexural lambda^2 (nu = 0.3).
LAMBDA2_P3 = {
    ("C", 0, 1): 10.215815, ("C", 1, 1): 21.260405, ("C", 2, 1): 34.877040,
    ("C", 0, 2): 39.771145, ("C", 1, 2): 60.828674, ("C", 2, 2): 84.582670,
    ("S", 0, 1): 4.935127, ("S", 1, 1): 13.898176, ("S", 2, 1): 25.613321,
    ("S", 0, 2): 29.719987, ("S", 1, 2): 48.478936, ("S", 2, 2): 70.117041,
}


def scan_roots(detf, lo, hi, n_grid):
    """Uniform-grid sign-change brackets of a real determinant."""
    ws = [lo + (hi - lo) * i / float(n_grid) for i in range(n_grid + 1)]
    vs = [float(mp.re(detf(w))) for w in ws]
    return [(ws[i], ws[i + 1]) for i in range(n_grid)
            if (vs[i] > 0) != (vs[i + 1] > 0)]


def refine(bisect, br, iters=60):
    return bisect(br[0], br[1], iters)


def lambda2(omega, a=0.6):
    E, nu, rho, H = LIU_DISK_KWARGS["E"], LIU_DISK_KWARGS["nu"], LIU_DISK_KWARGS["rho"], 2 * LIU_DISK_KWARGS["h"]
    D = E * H ** 3 / (12.0 * (1.0 - nu * nu))
    return omega * a * a * (rho * H / D) ** 0.5


def main():
    t0 = time.time()
    rows = []
    ok = {"G2a": True, "G2b": True, "G2c": True, "G2d": True}
    scan = make_liu_disk(dps=DPS_SCAN)
    prod = {p: make_liu_disk(dps=DPS_PROD, projection=p) for p in ("consistent", "duan")}
    weak = make_liu_disk(dps=DPS_PROD, e31=0.01)
    for bc in ("C", "S"):
        lo, hi = WIN[bc]
        for n in NS:
            # --- elastic limit, bare host ---
            br_e = scan_roots(lambda w: scan.elastic_disk_det(w, n, bc, h1=0.0), lo, hi, N_GRID)
            # --- SC coupled ---
            br_s = scan_roots(lambda w: scan.sc_disk_det(w, n, bc), lo, hi, N_GRID)
            if len(br_e) < 2 or len(br_s) < 2:
                ok["G2d"] = False
                print("G2d FAIL %s n=%d: elastic %d / SC %d roots in window" % (bc, n, len(br_e), len(br_s)))
                continue
            for m in (1, 2):
                key = (bc, n, m)
                we = prod["consistent"].elastic_disk_bisect(br_e[m - 1][0], br_e[m - 1][1], n, bc, h1=0.0, iters=60)
                l2 = lambda2(we)
                rel_l2 = l2 / LAMBDA2_P3[key] - 1.0
                wc = prod["consistent"].sc_disk_bisect(br_s[m - 1][0], br_s[m - 1][1], n, bc, iters=60)
                wd = prod["duan"].sc_disk_bisect(br_s[m - 1][0], br_s[m - 1][1], n, bc, iters=60)
                ww = weak.sc_disk_bisect(wc * 0.995, wc * 1.005, n, bc, iters=60)
                rel_kir = wc / LIU_KIR[key] - 1.0
                if abs(rel_l2) >= 1e-3:
                    ok["G2a"] = False
                if abs(rel_kir) >= 1e-3:
                    ok["G2b" if bc == "C" else "G2c"] = False
                rows.append(dict(
                    table=2 if bc == "C" else 6, bc=bc, n=n, m=m,
                    omega_consistent=wc, omega_duan=wd, omega_e31_0p01=ww,
                    x_kir=LIU_KIR[key], x_fem=LIU_FEM[key],
                    rel_kir=rel_kir, rel_fem_cited=wc / LIU_FEM[key] - 1.0,
                    D1_cons_vs_duan=wc / wd - 1.0, D2_coupling=wc / ww - 1.0,
                    omega_elastic_h1_0=we, lambda2=l2, lambda2_P3=LAMBDA2_P3[key],
                    rel_lambda2=rel_l2, f_consistent_Hz=wc / (2 * mp.pi)))
                print("%s n=%d m=%d  omega=%.4f  kir=%.1f  rel=%+.3e  | cons/duan-1=%+.1e  coupling=%+.1e  | lam2 rel=%+.1e"
                      % (bc, n, m, wc, LIU_KIR[key], rel_kir, wc / wd - 1, wc / ww - 1, rel_l2), flush=True)
    for r in rows:
        r["f_consistent_Hz"] = float(r["f_consistent_Hz"])
    verdict = "PASS_ALL" if all(ok.values()) and len(rows) == 12 else "FAIL"
    summ = dict(verdict=verdict, gates=ok, n_rows=len(rows), dps_scan=DPS_SCAN,
                dps_prod=DPS_PROD, n_grid=N_GRID, windows=WIN,
                max_abs_rel_kir=max(abs(r["rel_kir"]) for r in rows) if rows else None,
                max_abs_D1=max(abs(r["D1_cons_vs_duan"]) for r in rows) if rows else None,
                max_abs_D2=max(abs(r["D2_coupling"]) for r in rows) if rows else None,
                all_rel_kir_negative=all(r["rel_kir"] < 0 for r in rows),
                elapsed_s=time.time() - t0)
    base = os.path.join(HERE, "liu_disk_gate2_2026-09-24_results")
    with open(base + ".json", "w") as f:
        json.dump(dict(summary=summ, rows=rows), f, indent=1)
    with open(base + ".csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(json.dumps(summ, indent=1))
    print(verdict)


if __name__ == "__main__":
    main()
