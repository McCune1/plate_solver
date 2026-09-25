# Sine-enriched F-F n=0 OC root at 1/8 and 1/5, e31=-4.1 (extends extra_p4 1/12 row);
# model comparator for the FE OC/SC split (FE OC carries the interior potential).
import json, sys
sys.path.insert(0, "..")
from plate_solver.piezo_solver import PiezoOutOfPlaneSolver as P
out = {}
for rat, h1 in (("1/12", 0.02/12), ("1/8", 0.02/8), ("1/5", 0.02/5)):
    so = P(0.1,0.6,0.01,200e9,0.3,7800.0,h1=h1,C11E=132e9,C12E=71e9,C13E=73e9,C33E=115e9,
           e31=-4.1,e33=14.1,X11=7.124e-9,X33=5.841e-9,rho_pzt=7500.0,dps=60)
    el = so.elastic_bisect(600, 900, 0, iters=50)
    sc = so.coupled_bisect(el*0.99, el*1.01, 0, iters=50)
    oc = so.oc_ff_bisect(el*0.999, el*1.10, 0, iters=50)
    ocs = so.oc_ff_sine_bisect(oc*0.999, oc*1.001, 0, iters=50)
    out[rat] = dict(el=el, sc=sc, oc=oc, oc_sine=ocs, oc_over_sc=oc/sc-1, oc_sine_over_sc=ocs/sc-1,
                    oc_sine_over_oc=ocs/oc-1, sc_over_el=sc/el-1)
    print(rat, out[rat], flush=True)
json.dump(out, open("sine_oc_ratios_e31m_2026-09-24.json", "w"), indent=1)
