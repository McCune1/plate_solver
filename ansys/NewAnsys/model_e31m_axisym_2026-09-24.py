# -*- coding: utf-8 -*-
"""
model_e31m_axisym_2026-09-24.py -- model side for scoring the e31 = -4.1
axisymmetric clones (LESSONS Sec 18.244/18.245). Package solvers, default
projection='consistent', dps 30, both e31 signs (the +4.1 column is the
control: it must reproduce the numbers the parent jobs were scored against).
Chunk with MODEL_CASE (p4_12, p4_8, p4_5, p4_cf, p5_ff); appends to
model_e31m_axisym_2026-09-24.json.
"""
import json, math, os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "..")))
from plate_solver.piezo_solver import PiezoOutOfPlaneSolver as P4
from plate_solver.piezo_monolithic import PiezoMonolithicOutOfPlaneSolver as P5
OUT = os.path.join(HERE, "model_e31m_axisym_2026-09-24.json")
H = 0.01


def p4(den, e31):
    return P4(0.1, 0.6, H, 200e9, 0.3, 7800.0, h1=2 * H / den, C11E=132e9, C12E=71e9,
              C13E=73e9, C33E=115e9, e31=e31, e33=14.1, X11=7.124e-9, X33=5.841e-9,
              rho_pzt=7500.0, dps=30)


def ff(den):
    out = {}
    for e31 in (4.1, -4.1):
        s = p4(den, e31)
        rows = []
        for (lo, hi) in ((600, 900), (3300, 4000)):
            el = s.elastic_bisect(lo, hi, 0, iters=50)
            sc = s.coupled_bisect(el * 0.99, el * 1.01, 0, iters=45)
            oc = s.oc_ff_bisect(el * 0.999, el * 1.10, 0, iters=50)
            rows.append(dict(el=el, sc=sc, oc=oc, oc_sc=oc / sc - 1))
        out["%+.1f" % e31] = rows
    return out


def cf():
    out = {}
    for e31 in (4.1, -4.1):
        s = p4(12, e31)
        el = s.elastic_cf_bisect(380, 460, 0, iters=50)
        sc = s.cf_coupled_bisect(el * 0.99, el * 1.01, 0, iters=45)
        oc = s.oc_cf_bisect(el * 0.999, el * 1.10, 0, iters=50)
        out["%+.1f" % e31] = [dict(el=el, sc=sc, oc=oc, oc_sc=oc / sc - 1)]
    return out


def p5ff():
    out = {}
    for e31 in (4.1, -4.1):
        s = P5(0.1, 0.6, 0.01, 132e9, 71e9, 73e9, 115e9, 7500.0, e31=e31, e33=14.1,
               X11=7.124e-9, X33=5.841e-9, dps=30)
        el = s.elastic_bisect(440, 480, 0)
        grid = [el * (1 + 0.3 * i / 60) for i in range(61)]
        v = [s.coupled_det(w, 0).real > 0 for w in grid]
        i = [k for k in range(60) if v[k] != v[k + 1]][0]
        sc = s.coupled_bisect(grid[i], grid[i + 1], 0)
        out["%+.1f" % e31] = [dict(el=el, sc=sc, sc_el=sc / el - 1)]
    return out


CASES = {"p4_12": lambda: ff(12), "p4_8": lambda: ff(8), "p4_5": lambda: ff(5),
         "p4_cf": cf, "p5_ff": p5ff}

if __name__ == "__main__":
    case = os.environ["MODEL_CASE"]
    t0 = time.time()
    res = CASES[case]()
    data = json.load(open(OUT)) if os.path.exists(OUT) else {}
    data[case] = res
    json.dump(data, open(OUT, "w"), indent=1)
    print(case, json.dumps(res), "%.0fs" % (time.time() - t0))
