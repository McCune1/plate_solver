# Sign-bracketed det roots of the unconstrained (u) and Lagrange-constrained (c) IP cantilever K
# on a COMPLETE (argument-principle-verified) tracked basis. Args: n lo hi step invkey out
import sys, json, numpy as np
from setup import *
from plate_solver.detectors import equilibrated_logdet, sigma_min_from_K
n=int(sys.argv[1]); lo=float(sys.argv[2]); hi=float(sys.argv[3]); step=float(sys.argv[4]); key=sys.argv[5]; out=sys.argv[6]
inv0=[complex(a,b) for a,b in json.load(open('inv_cache%s_%s.json'%(TAG,key)))[key]]
NINV=len(inv0)
def basis_at(Om):
    # track from the cached centre inventory in small steps (robust continuation)
    Oc=float(key); m=max(1,int(abs(Om-Oc)/0.002)+1); r=inv0
    for t in np.linspace(Oc,Om,m+1)[1:]:
        r=track(ip.fast,float(t),r)
    if len(r)!=NINV: raise RuntimeError('track lost roots at %g: %d'%(Om,len(r)))
    return r
def ev(Om,lag):
    sel,cnt=select_fill(basis_at(Om),n)
    K,s=ip._build_K_real(Om,sel,lag,fast_scan=False)
    sg,ld=equilibrated_logdet(K,s)
    return sg,float(ld),sigma_min_from_K(K,s),cnt
def log(d):
    with open(out,'a') as f: f.write(json.dumps(d)+"\n")
grid=np.arange(lo,hi+step/2,step)
prev={}
for Om in grid:
    Om=float(Om); row=dict(kind='g',n=n,Om=Om)
    for lag,tag in ((False,'u'),(True,'c')):
        sg,ld,sm,cnt=ev(Om,lag); row[tag]=[sg,ld,sm]; row['cnt']=cnt
        if tag in prev and prev[tag][1]!=sg:
            a,b,sa=prev[tag][0],Om,prev[tag][1]
            for _ in range(11):
                m=(a+b)/2; s2=ev(m,lag)[0]
                if s2==sa: a=m
                else: b=m
            r=(a+b)/2; sgr,ldr,smr,_=ev(r,lag)
            log(dict(kind='root',n=n,mat=tag,Om=r,ld=ldr,smin=smr,bracket=[prev[tag][0],Om]))
        prev[tag]=(Om,sg)
    log(row)
