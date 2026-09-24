import numpy as np
from setup import *
from plate_solver.detectors import canon
def argcount(Om, X, Y, npts=6000, x0=-0.05, y0=-0.05):
    # counterclockwise rectangle; count zeros inside via winding of det
    pts=[]
    per = 2*(X-x0)+2*(Y-y0)
    def seg(a,b):
        L=abs(b-a); m=max(20,int(npts*L/per))
        return [a+(b-a)*t for t in np.linspace(0,1,m,endpoint=False)]
    c=[complex(x0,y0),complex(X,y0),complex(X,Y),complex(x0,Y)]
    for i in range(4): pts+=seg(c[i],c[(i+1)%4])
    pts.append(pts[0])
    v=np.array([ip.fast.det(z,Om) for z in pts])
    ph=np.angle(v); d=np.diff(ph); d=(d+np.pi)%(2*np.pi)-np.pi
    maxjump=np.max(np.abs(d))
    return np.sum(d)/(2*np.pi), maxjump
def inventory(Om, X=12, Y=45, base=None):
    roots=list(base) if base is not None else list(full_search(ip.fast, Om, xmax=20.0))
    def add(z):
        if z is None or not np.isfinite(z.real) or not np.isfinite(z.imag): return
        z=canon(z)
        if abs(z)<1e-6 or z.real>X or z.imag>Y: return
        if abs(ip.fast.det(z,Om))>1e-6: return
        if any(abs(z-u)<0.05 for u in roots): return
        roots.append(z)
    for a in np.arange(0.5, 8.01, 1.0):
        for b in np.arange(0.5, Y+0.01, 1.0):
            add(ip.fast.newton(complex(a,b),Om,itmax=30))
    for x in np.linspace(0.05,X,300):
        add(ip.fast.newton(complex(x,0),Om,itmax=30))
    roots.sort(key=lambda z:(abs(z.imag),abs(z)))
    return roots
if __name__=="__main__":
    import sys
    Om=float(sys.argv[1])
    t=time.time(); r=inventory(Om); print('inv',round(time.time()-t,1),len(r),[complex(round(z.real,3),round(z.imag,3)) for z in r],flush=True)
    for (X,Y) in ((12,22),(12,45),(30,45)):
        t=time.time(); w,mj=argcount(Om,X,Y); print('argcount',X,Y,round(w,3),'maxjump',round(mj,3),'found_inside',sum(1 for z in r if z.real<X and z.imag<Y),round(time.time()-t,1),flush=True)
