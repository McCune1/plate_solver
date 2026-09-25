# Paper 4 numbers at e31 = -4.1 not produced by an existing probe (LESSONS Sec 18.246):
#  (a) all 27 Duan Table 4 rows, consistent AND duan projection, and the
#      largest consistent-vs-duan root difference;
#  (b) the sine-enriched open-circuit root at 1/12 (Discussion: the interior
#      and bus stiffenings add).
import json, sys, time
sys.path.insert(0, "..")
from plate_solver.piezo_solver import PiezoOutOfPlaneSolver as P
T = {(0,0,"1/12"):2792,(0,0,"1/8"):2853,(0,0,"1/5"):2989,(0,1,"1/12"):7723,(0,1,"1/8"):7891,(0,1,"1/5"):8268,
     (0,2,"1/12"):15182,(0,2,"1/8"):15512,(0,2,"1/5"):16253,(1,0,"1/12"):2928,(1,0,"1/8"):2992,(1,0,"1/5"):3135,
     (1,1,"1/12"):7965,(1,1,"1/8"):8138,(1,1,"1/5"):8527,(1,2,"1/12"):15482,(1,2,"1/8"):15818,(1,2,"1/5"):16574,
     (2,0,"1/12"):3477,(2,0,"1/8"):3553,(2,0,"1/5"):3723,(2,1,"1/12"):8768,(2,1,"1/8"):8959,(2,1,"1/5"):9387,
     (2,2,"1/12"):16434,(2,2,"1/8"):16791,(2,2,"1/5"):17594}
def s(rat, proj, e31=-4.1, dps=40):
    n, d = map(int, rat.split("/"))
    return P(0.1,0.6,0.01,200e9,0.3,7800.0,h1=0.02*n/d,C11E=132e9,C12E=71e9,C13E=73e9,C33E=115e9,
             e31=e31,e33=14.1,X11=7.124e-9,X33=5.841e-9,rho_pzt=7500.0,dps=dps,projection=proj)
out = {"table4": []}
t0 = time.time()
for (p, nr, rat), tg in sorted(T.items()):
    row = dict(p=p, n_radial=nr, h1=rat, target=tg)
    for proj in ("consistent", "duan"):
        so = s(rat, proj)
        lo, hi = tg*0.985, tg*1.015
        flip = (so.cc_coupled_det(lo, p).real > 0) != (so.cc_coupled_det(hi, p).real > 0)
        w = so.cc_coupled_bisect(lo, hi, p, iters=45)
        row[proj] = dict(omega=w, rel_pct=100*(w/tg-1), flip=bool(flip))
    row["cons_minus_duan_rel"] = row["consistent"]["omega"]/row["duan"]["omega"] - 1
    out["table4"].append(row)
    print(p, nr, rat, tg, "%.2f %+.4f%%  duan %.2f %+.4f%%" % (row["consistent"]["omega"], row["consistent"]["rel_pct"],
          row["duan"]["omega"], row["duan"]["rel_pct"]), flush=True)
c = [r["consistent"]["rel_pct"] for r in out["table4"]]
d = [r["duan"]["rel_pct"] for r in out["table4"]]
out["worst_consistent_pct"] = max(c, key=abs); out["worst_duan_pct"] = max(d, key=abs)
out["max_cons_duan_rel"] = max(abs(r["cons_minus_duan_rel"]) for r in out["table4"])
out["all_flip"] = all(r["consistent"]["flip"] and r["duan"]["flip"] for r in out["table4"])
so = s("1/12", "consistent", dps=60)
el = so.elastic_bisect(600, 900, 0, iters=50)
sc = so.coupled_bisect(el*0.99, el*1.01, 0, iters=50)
oc = so.oc_ff_bisect(el*0.999, el*1.10, 0, iters=50)
ocs = so.oc_ff_sine_bisect(oc*0.999, oc*1.001, 0, iters=50)
out["sine_oc_1_12"] = dict(el=el, sc=sc, oc=oc, oc_sine=ocs, oc_sine_over_oc=ocs/oc-1, sc_over_el=sc/el-1)
print("sine OC:", out["sine_oc_1_12"], flush=True)
print("worst cons %.4f%% duan %.4f%% max cons-duan %.2e all_flip %s (%.0fs)" % (out["worst_consistent_pct"], out["worst_duan_pct"], out["max_cons_duan_rel"], out["all_flip"], time.time()-t0))
json.dump(out, open("extra_p4_e31m_2026-09-24.json", "w"), indent=1)
