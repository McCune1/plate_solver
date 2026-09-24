# Fine scan of the m2 neighbourhood: unconstrained + constrained sigma_min and det sign on a complete tracked basis
import sys, json, numpy as np
from setup import *
from inv import inventory
from plate_solver.detectors import equilibrated_logdet, sigma_min_from_K
def svals(K,size,passes=2):
    K=K.copy()
    for _ in range(passes):
        for i in range(size):
            s=max(abs(K[i,j]) for j in range(size))
            if s>0:
                for j in range(size): K[i,j]/=s
        for j in range(size):
            s=max(abs(K[i,j]) for i in range(size))
            if s>0:
                for i in range(size): K[i,j]/=s
    A=np.array([[complex(K[i,j]) for j in range(size)] for i in range(size)])
    sv=np.linalg.svd(A,compute_uv=False)
    return [float(np.log10(max(x,1e-300))) for x in sv[-4:][::-1]]
n=int(sys.argv[1]); lo=float(sys.argv[2]); hi=float(sys.argv[3]); npts=int(sys.argv[4]); out=sys.argv[5]
inv0=inventory(0.82); NINV=len(inv0)
prev=inv0
for Om in np.linspace(lo,hi,npts):
    r=track(ip.fast,float(Om),prev)
    if len(r)!=NINV: r=inventory(float(Om))
    prev=r
    sel,cnt=select_fill(r,n)
    Ku,su=ip._build_K_real(float(Om),sel,False,fast_scan=False)
    sg,ld=equilibrated_logdet(Ku,su)
    sv_u=svals(Ku,su); s_u=sv_u[0]
    Kc,sc=ip._build_K_real(float(Om),sel,True,fast_scan=False)
    sv_c=svals(Kc,sc); s_c=sv_c[0]
    sgc,ldc=equilibrated_logdet(Kc,sc)
    with open(out,'a') as f: f.write(json.dumps(dict(n=n,cnt=cnt,Om=float(Om),su=s_u,sc=s_c,svu=sv_u,svc=sv_c,sign=sg,ld=float(ld),signc=sgc,ldc=float(ldc),nr=len(r)))+"\n")
