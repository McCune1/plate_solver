#!/usr/bin/env python3
"""Pre-registered Kirchhoff targets for the Paper 4 3-D thin-limit check.

Written 2026-09-23, BEFORE any thin-limit deck ran. Closes (or not) the G3
FAIL of job 2530423 (LESSONS 18.232): is the 0.5-1.2% FE-below-Kirchhoff
offset of the elastic n >= 1 roots thick-plate physics, or an n-dependent
defect in the thin-plate operator?

Bare steel annulus (no skins): r_i = 0.1 m, r_o = 0.6 m, E = 200 GPa,
nu = 0.3, rho = 7800, full thickness t = 20, 10, 5 mm. Roots of the
elastic F-F 4x4 of plate_solver/piezo_solver.py (h1 = 0, h = t/2) --
the same operator, with the Kirchhoff effective-shear row, that job
2530423 scored in G3 (it equals ring_disk to 1e-6, LESSONS 18.226).
All roots below FMAX(t) = 1100 Hz * t/20 mm, n = 0..7, found from scratch
(701-point sign scan at dps 40, bisect at dps 100).

Kirchhoff frequencies scale exactly with t; the script checks that
(f(t)/t equal to 1e-9 across the three t) as a self-test.

Imports by file path, needs only mpmath. Writes
targets_p4_3d_thinlimit_2026-09-23.json next to itself.
"""
import importlib.util
import json
import math
import os

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.normpath(os.path.join(HERE, "..", "..", "plate_solver",
                                    "piezo_solver.py"))
OUT = os.path.join(HERE, "targets_p4_3d_thinlimit_2026-09-23.json")
THICK_MM = (20, 10, 5)
NMAX = 7


def load(path=PKG):
    spec = importlib.util.spec_from_file_location("piezo_solver_p4", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def roots(ps, t_mm):
    h = t_mm / 2000.0
    s40 = ps.PiezoOutOfPlaneSolver(r_i=0.1, r_o=0.6, h=h, E=200e9, nu=0.3,
                                   rho=7800.0, h1=0.0, dps=40)
    s = ps.PiezoOutOfPlaneSolver(r_i=0.1, r_o=0.6, h=h, E=200e9, nu=0.3,
                                 rho=7800.0, h1=0.0, dps=100)
    fmax = 2 * math.pi * 1100.0 * t_mm / 20.0
    lo0 = 20.0 * t_mm / 20.0
    out = {}
    for n in range(NMAX + 1):
        grid = [lo0 + (fmax - lo0) * i / 700 for i in range(701)]
        vals = [s40.elastic_det(w, n).real for w in grid]
        k = 0
        for i in range(700):
            if (vals[i] > 0) == (vals[i + 1] > 0):
                continue
            k += 1
            el = s.elastic_bisect(grid[i], grid[i + 1], n, iters=60)
            out["%d_%d" % (n, k)] = dict(n=n, s=k, el_rad=el,
                                         el_hz=el / 2 / math.pi)
    return out


def main():
    ps = load()
    data = {"meta": dict(written="2026-09-23", solver=os.path.relpath(PKG, HERE),
                         geometry="bare steel r_i=0.1 r_o=0.6 m, E=200 GPa, "
                                  "nu=0.3, rho=7800", thickness_mm=list(THICK_MM))}
    for t in THICK_MM:
        data["t%d" % t] = roots(ps, t)
        print("t=%d mm:" % t, ", ".join("%s %.6f" % (k, v["el_hz"])
              for k, v in sorted(data["t%d" % t].items(),
                                 key=lambda kv: kv[1]["el_hz"])))
    # self-test: Kirchhoff f is exactly proportional to t
    ref = data["t20"]
    worst = 0.0
    for t in THICK_MM[1:]:
        for key, v in data["t%d" % t].items():
            if key in ref:
                r = (v["el_hz"] / t) / (ref[key]["el_hz"] / 20.0) - 1.0
                worst = max(worst, abs(r))
        assert set(data["t%d" % t]) == set(ref), "root sets differ at t=%d" % t
    assert worst < 1e-9, worst
    data["meta"]["scaling_selftest_worst"] = worst
    with open(OUT, "w") as f:
        json.dump(data, f, indent=1, sort_keys=True)
    print("scaling self-test worst %.2e; wrote %s" % (worst, OUT))


if __name__ == "__main__":
    main()
