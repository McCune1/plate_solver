import os, sys, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))  # repo root (plate_solver/)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mpmath import mp
from plate_solver.geometry import make_geometry, _material_from_env
from plate_solver.core_solvers import InPlaneSolver
from plate_solver.detectors import full_search, select_fill, track
from plate_solver.validation import _scan_cfg_part2
mp.dps = int(os.environ.get('DPS','40'))
mat = _material_from_env()
TWOT=float(os.environ.get('TWOT','0.25'))
geom = make_geometry(1.25, TWOT)
ip = InPlaneSolver(geom, mat, M=80, n_quad=30)
(Om_lo, Om_hi), ns, cut_offs, ze_max = _scan_cfg_part2(r0_over_2b=1.25, two_T_pi=TWOT, mat=mat)
def smin(Om, n, xmax=None, ret_sel=False):
    raw = full_search(ip.fast, Om, xmax=xmax or ze_max)
    sel, cnt = select_fill(raw, n)
    s = float(ip.sigma_min(Om, sel, lagrange=False))
    return (s, cnt, sel) if ret_sel else (s, cnt)

TAG='' if TWOT==0.25 else '_a%g'%TWOT
