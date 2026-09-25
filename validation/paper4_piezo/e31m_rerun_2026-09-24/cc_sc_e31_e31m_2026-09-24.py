# e31 = -4.1 copy (LESSONS Sec 18.246) of consistent_rerun_2026-09-23/cc_sc_e31.py:
# C-C short-circuit e31 sensitivity (the "SC channel cannot see e31" negative
# control of Paper 4). Only e31 changes; the bracket is the model C-C SC root +-15.
import sys, os
sys.path.insert(0, "..")
from plate_solver.piezo_solver import PiezoOutOfPlaneSolver as P
H=0.01
E31=-4.1
k4=dict(r_i=0.1,r_o=0.6,h=H,E=200e9,nu=0.3,rho=7800.0,C11E=132e9,C12E=71e9,C13E=73e9,C33E=115e9,e33=14.1,X11=7.124e-9,X33=5.841e-9,rho_pzt=7500.0,dps=60)
tg={2/12:2791.78,2/8:2852.19,2/5:2987.72}
for pr in ['duan','consistent']:
    for frac,om in tg.items():
        de=1e-3
        f=lambda e: P(h1=frac*H,e31=e,projection=pr,**k4).cc_coupled_bisect(om-15,om+15,0,iters=60)
        wp=f(E31+de); wm=f(E31-de); w0=f(E31)
        d=(wp-wm)/(2*de)
        print(pr,frac/2,w0,'dw/de31=%.6e'%d,'unc@0.01%%=%.2f%%'%(100*1e-4*w0/abs(d)/abs(E31)),flush=True)
