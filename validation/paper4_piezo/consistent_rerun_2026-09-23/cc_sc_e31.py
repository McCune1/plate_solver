from plate_solver.piezo_solver import PiezoOutOfPlaneSolver as P
H=0.01
k4=dict(r_i=0.1,r_o=0.6,h=H,E=200e9,nu=0.3,rho=7800.0,C11E=132e9,C12E=71e9,C13E=73e9,C33E=115e9,e33=14.1,X11=7.124e-9,X33=5.841e-9,rho_pzt=7500.0,dps=60)
tg={2/12:2791.78,2/8:2852.19,2/5:2987.72}
for pr in ['duan','consistent']:
    for frac,om in tg.items():
        de=1e-3
        f=lambda e: P(h1=frac*H,e31=e,projection=pr,**k4).cc_coupled_bisect(om-15,om+15,0,iters=60)
        wp=f(4.1+de); wm=f(4.1-de); w0=f(4.1)
        d=(wp-wm)/(2*de)
        print(pr,frac/2,w0,'dw/de31=%.6e'%d,'unc@0.01%%=%.2f%%'%(100*1e-4*w0/abs(d)/4.1),flush=True)
