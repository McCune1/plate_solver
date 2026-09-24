# IP cantilever (1.25, 0.25pi) basis convergence with an argument-principle-verified COMPLETE root inventory.
import sys, json, numpy as np
from setup import *
from inv import inventory
Om0=float(sys.argv[1]); ns=[int(x) for x in sys.argv[2].split(',')]; out=sys.argv[3]
half=float(sys.argv[4]) if len(sys.argv)>4 else 0.04
NQ=int(os.environ.get('NQ','30')); MM=int(os.environ.get('MM','80'))
if NQ!=30 or MM!=80:
    import setup; setup.ip = InPlaneSolver(geom, mat, M=MM, n_quad=NQ); ip=setup.ip
inv0 = inventory(Om0)
NINV=len(inv0)
cache={}
def basis(Om):
    key=round(Om,12)
    if key in cache: return cache[key]
    r=track(ip.fast, Om, inv0)
    if len(r)!=NINV: r=inventory(Om)
    cache[key]=r; return r
def f(Om,n):
    r=basis(Om); sel,cnt=select_fill(r,n)
    return float(ip.sigma_min(Om, sel, lagrange=False)), cnt, len(r)
def log(d):
    with open(out,'a') as fh: fh.write(json.dumps(d)+"\n")
log(dict(kind='inv',Om0=Om0,n_inv=NINV,NQ=NQ,MM=MM,dps=mp.dps,roots=[[z.real,z.imag] for z in inv0]))
grid=np.linspace(Om0*(1-half),Om0*(1+half),21)
g=(5**0.5-1)/2
for n in ns:
    t=time.time()
    vals=[]
    for Om in grid:
        s,cnt,nr=f(float(Om),n); vals.append(s)
        log(dict(kind='scan',n=n,Om=float(Om),s=s,cnt=cnt,nr=nr))
    k=int(np.argmin(vals))
    edge = k in (0,len(grid)-1)
    a=grid[max(k-1,0)]; b=grid[min(k+1,len(grid)-1)]
    x1=b-g*(b-a); x2=a+g*(b-a); f1=f(x1,n)[0]; f2=f(x2,n)[0]
    for _ in range(22):
        if f1<f2: b,x2,f2=x2,x1,f1; x1=b-g*(b-a); f1=f(x1,n)[0]
        else: a,x1,f1=x1,x2,f2; x2=a+g*(b-a); f2=f(x2,n)[0]
    om=(a+b)/2; s,cnt,nr=f(om,n)
    shoulders=(vals[0],vals[-1])
    log(dict(kind='min',n=n,cnt=cnt,Om=om,s=s,grid_argmin_edge=edge,shoulder_lo=vals[0],shoulder_hi=vals[-1],t=round(time.time()-t,1)))
