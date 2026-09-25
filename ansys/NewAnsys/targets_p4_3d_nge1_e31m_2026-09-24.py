#!/usr/bin/env python3
"""e31 = -4.1 re-target (LESSONS Sec 18.244) of targets_p4_3d_nge1_2026-09-23.py.
Identical root search; e31 = -4.1; the SC bracket widened to +-1% and the OC
bracket to -0.1%..+10% (the n = 0 OC split is ~2-4% at -4.1, outside the old
+1.5% bracket). Writes targets_p4_3d_nge1_e31m_2026-09-24.json.

Original docstring:
Pre-registered model targets for the Paper 4 3-D FE n >= 1 check.

Written BEFORE any 3-D deck ran (2026-09-23). The scorer reads the JSON
this writes (targets_p4_3d_nge1_2026-09-23.json); it never recomputes a
target from FE output.

For h1/2h = 1/12 and 1/5, and n = 0..7, every free-free flexural root of
the elastic bilayer 4x4 below 1100 Hz is found from scratch (sign flips of
elastic_det on a 701-point grid, dps 40, then elastic_bisect at dps 100).
Each is then refined in the coupled short-circuit 6x6 (coupled_bisect,
+-0.2% bracket) and the linear open-circuit 4x4 (oc_ff_bisect, -0.1%..
+1.5% bracket), with a sign-flip count on each bracket recorded so a
bracket that caught nothing (or two roots) is visible.

Solver: plate_solver/piezo_solver.py PiezoOutOfPlaneSolver, default
projection='consistent' (LESSONS Sec 18.231), SI units, Duan 2005 steel/
PZT-4 constants -- the same constructor call as tests/test_solver.py.

What the model says, so the FE is read against it and not the reverse:
  * sc_el_split = SC 6x6 / elastic 4x4 - 1: the interior-potential
    (quadratic bubble) stiffening. Tiny: ~1e-6 to 7e-5 relative.
  * oc(n >= 1) is bit-identical to the elastic 4x4 (no net electrode
    charge). The linear OC model omits the interior bubble, which the
    paper (Discussion) shows simply adds, so the physically complete
    statement is OC = SC at n >= 1. That is what the FE is scored on.

Import is by file path so only mpmath is needed (the device and the
cluster venv both have it). Run from Ansys/NewAnsys/ or anywhere:
    python3 targets_p4_3d_nge1_2026-09-23.py
About 10 min for 1/12 and 20 min for 1/5 on one core.
"""
# (no __future__ import: the compute nodes run python3 3.6, job 2530336)

import importlib.util
import json
import math
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.normpath(os.path.join(HERE, "..", "..", "plate_solver",
                                    "piezo_solver.py"))
OUT = os.path.join(HERE, "targets_p4_3d_nge1_e31m_2026-09-24.json")

H = 0.01
FMAX_HZ = 1100.0
NMAX = 7


def load_solver_module(path=PKG):
    spec = importlib.util.spec_from_file_location("piezo_solver_p4", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def make(ps, den, dps):
    return ps.PiezoOutOfPlaneSolver(
        r_i=0.1, r_o=0.6, h=H, E=200e9, nu=0.3, rho=7800.0,
        h1=2.0 * H / den, C11E=132e9, C12E=71e9, C13E=73e9, C33E=115e9,
        e31=-4.1, e33=14.1, X11=7.124e-9, X33=5.841e-9, rho_pzt=7500.0,
        dps=dps)


def flips(fn, lo, hi, m=16):
    g = [lo + (hi - lo) * i / m for i in range(m + 1)]
    v = [fn(w).real > 0 for w in g]
    return sum(1 for i in range(m) if v[i] != v[i + 1])


def roots_for(ps, den, log=print):
    s40 = make(ps, den, 40)
    s = make(ps, den, 100)
    fmax = 2 * math.pi * FMAX_HZ
    out = {}
    for n in range(0, NMAX + 1):
        grid = [20 + (fmax - 20) * i / 700 for i in range(701)]
        vals = [s40.elastic_det(w, n).real for w in grid]
        k = 0
        for i in range(700):
            if (vals[i] > 0) == (vals[i + 1] > 0):
                continue
            k += 1
            t0 = time.time()
            el = s.elastic_bisect(grid[i], grid[i + 1], n, iters=60)
            lo, hi = el * (1 - 1e-2), el * (1 + 1e-2)
            olo, ohi = el * (1 - 1e-3), el * (1 + 1.0e-1)
            nf_sc = flips(lambda w: s40.coupled_det(w, n), lo, hi)
            nf_oc = flips(lambda w: s40.oc_ff_det(w, n), olo, ohi, 128)
            sc = s.coupled_bisect(lo, hi, n, iters=55)
            oc = s.oc_ff_bisect(olo, ohi, n, iters=60)
            rec = dict(
                n=n, s=k, el_rad=el, sc_rad=sc, oc_rad=oc,
                el_hz=el / 2 / math.pi, sc_hz=sc / 2 / math.pi,
                oc_hz=oc / 2 / math.pi,
                sc_el_split=sc / el - 1.0, oc_sc_split=oc / sc - 1.0,
                oc_el_split=oc / el - 1.0,
                bracket_flips_sc=nf_sc, bracket_flips_oc=nf_oc)
            out["%d_%d" % (n, k)] = rec
            log("h1/2h=1/%d n=%d s=%d el=%.9f Hz sc/el-1=%.4e oc/el-1=%.4e "
                "flips sc=%d oc=%d (%.0fs)"
                % (den, n, k, rec["el_hz"], rec["sc_el_split"],
                   rec["oc_el_split"], nf_sc, nf_oc, time.time() - t0))
    return out


def main(argv):
    ps = load_solver_module()
    dens = [int(a) for a in argv[1:]] or [12, 5]
    data = {}
    if os.path.isfile(OUT):
        with open(OUT) as f:
            data = json.load(f)
    data.setdefault("meta", {})
    data["meta"].update(dict(
        written="2026-09-24", e31=-4.1, solver=os.path.relpath(PKG, HERE),
        projection="consistent", fmax_hz=FMAX_HZ, nmax=NMAX,
        geometry="r_i=0.1 r_o=0.6 h=0.01 (SI), Duan 2005 steel/PZT-4"))
    for den in dens:
        data["h1%d" % den] = roots_for(ps, den)
        with open(OUT, "w") as f:
            json.dump(data, f, indent=1, sort_keys=True)
    print("wrote", OUT)


if __name__ == "__main__":
    main(sys.argv)
