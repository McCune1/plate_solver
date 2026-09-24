#!/usr/bin/env python3
"""probe_p1_cant_ip_q9fe_2026-09-24.py -- independent plane-stress FE for Paper 1's
in-plane cantilever rows (r0/(2b)=1.25; 2Theta/pi = 0.25, 0.50, 1.00; nu=0.35).

Q2 (9-node Lagrange) quadrilaterals on an (r, theta) tensor grid with the EXACT
polar geometry map (x = r cos t, y = r sin t inside every element), consistent
mass, 4x4 Gauss, plane stress D = E/(1-nu^2)[[1,nu,0],[nu,1,0],[0,0,(1-nu)/2]].
Wall theta=0 clamped (ux=uy=0); wall theta=PHI and both arcs free.
Shift-invert Lanczos (scipy eigsh, sigma=0). Richardson over N = 16/32/64.

Self-check (--validate): the C-C deck geometry of
Ansys/ansys_cc_ip_annular_ext.inp (R_i=4, R_o=8, PHI=90, nu=0.30, both walls
clamped) against its PLANE183 cc_ip_mesh96.txt: 177.736 218.821 299.912 411.818
482.641 527.514 586.336 Hz. 2026-09-24 container run, N=32: 177.747 218.851
299.925 411.861 482.668 527.533 586.385 (<= 0.014%).

Result (2026-09-24 container run; Omega_bar = 2 pi f / wbar, wbar = 2480.106584):
  0.25: 0.348738 0.823991 0.940601 (mode 4 1.551622)
  0.50: 0.121449 0.353054 0.529883 (mode 4 0.943960)
  1.00: 0.039555 0.105826 0.261755 (mode 4 0.458798)
Observed Richardson order p ~ 1.5 (clamped/free corner singularity); the
64-vs-Richardson change is <= 0.01% for every row.
Needs numpy + scipy (cluster venv, or the Cowork container). Runs in ~1 min.
See LESSONS_LEARNED.md Sec 18.233.
"""
import math, numpy as np, scipy.sparse as sp, scipy.sparse.linalg as sla
def run(Ri,Ro,PHI,E,nu,rho,Nr,Nt,nmodes=8,sigma=None,clamp='theta0',order=2):
    p=order; nr=p*Nr+1; nt=p*Nt+1
    rn=np.linspace(Ri,Ro,nr); tn=np.linspace(0,PHI,nt)
    nid=lambda i,j: i*nt+j
    ndof=2*nr*nt
    # Gauss
    g,w=np.polynomial.legendre.leggauss(p+2)
    # 1D Lagrange on equispaced p+1 nodes
    xi_n=np.linspace(-1,1,p+1)
    def lag(x):
        L=np.ones((p+1,)+np.shape(x)); dL=np.zeros((p+1,)+np.shape(x))
        for a in range(p+1):
            for b in range(p+1):
                if b==a: continue
                L[a]*= (x-xi_n[b])/(xi_n[a]-xi_n[b])
            s=0
            for c in range(p+1):
                if c==a: continue
                t=1/(xi_n[a]-xi_n[c])
                for b in range(p+1):
                    if b in (a,c): continue
                    t=t*(x-xi_n[b])/(xi_n[a]-xi_n[b])
                s=s+t
            dL[a]=s
        return L,dL
    L,dL=lag(g)  # (p+1, ng)
    ng=len(g)
    # element-local shape fns N[a,b] at gauss (q1,q2)
    Nsh=np.einsum('aq,bs->abqs',L,L).reshape((p+1)**2,ng,ng)
    dNx=np.einsum('aq,bs->abqs',dL,L).reshape((p+1)**2,ng,ng)  # d/dxi (radial)
    dNe=np.einsum('aq,bs->abqs',L,dL).reshape((p+1)**2,ng,ng)  # d/deta (theta)
    W=np.outer(w,w)
    D=E/(1-nu**2)*np.array([[1,nu,0],[nu,1,0],[0,0,(1-nu)/2]])
    rows=[];cols=[];kv=[];mv=[]
    nen=(p+1)**2
    for ei in range(Nr):
        r1,r2=rn[p*ei],rn[p*ei+p]
        for ej in range(Nt):
            t1,t2=tn[p*ej],tn[p*ej+p]
            R=(r1+r2)/2+(r2-r1)/2*g[:,None]*np.ones((1,ng)); T=(t1+t2)/2+(t2-t1)/2*g[None,:]*np.ones((ng,1))
            # x=r cos t, y= r sin t ; dx/dxi = (dr/dxi) cos t, dx/deta = -r sin t dt/deta
            drx=(r2-r1)/2; dte=(t2-t1)/2
            J11=drx*np.cos(T); J12=drx*np.sin(T); J21=-R*np.sin(T)*dte; J22=R*np.cos(T)*dte
            det=J11*J22-J12*J21
            # inverse: [dxi/dx dxi/dy; deta/dx deta/dy]
            i11=J22/det; i12=-J21/det; i21=-J12/det; i22=J11/det
            dNdx=dNx*i11+dNe*i21; dNdy=dNx*i12+dNe*i22
            B=np.zeros((3,2*nen,ng,ng))
            B[0,0::2]=dNdx; B[1,1::2]=dNdy; B[2,0::2]=dNdy; B[2,1::2]=dNdx
            wd=W*det
            Ke=np.einsum('iaqs,ij,jbqs,qs->ab',B,D,B,wd)
            Me1=np.einsum('aqs,bqs,qs->ab',Nsh,Nsh,wd)*rho
            Me=np.zeros((2*nen,2*nen)); Me[0::2,0::2]=Me1; Me[1::2,1::2]=Me1
            loc=[nid(p*ei+a,p*ej+b) for a in range(p+1) for b in range(p+1)]
            dofs=np.array([[2*n,2*n+1] for n in loc]).ravel()
            rr,cc=np.meshgrid(dofs,dofs,indexing='ij')
            rows.append(rr.ravel()); cols.append(cc.ravel()); kv.append(Ke.ravel()); mv.append(Me.ravel())
    rows=np.concatenate(rows); cols=np.concatenate(cols)
    K=sp.csr_matrix((np.concatenate(kv),(rows,cols)),shape=(ndof,ndof))
    M=sp.csr_matrix((np.concatenate(mv),(rows,cols)),shape=(ndof,ndof))
    fixed=set()
    if clamp in ('theta0','both'):
        for i in range(nr): n=nid(i,0); fixed|={2*n,2*n+1}
    if clamp=='both':
        for i in range(nr):
            for jj in (0,nt-1): n=nid(i,jj); fixed|={2*n,2*n+1}
    free=np.array(sorted(set(range(ndof))-fixed))
    K=K[free][:,free]; M=M[free][:,free]
    if sigma is None: sigma=0.0
    vals,vecs=sla.eigsh(K,k=nmodes,M=M,sigma=sigma,which='LM')
    om=np.sqrt(np.abs(np.sort(vals)))
    return om/(2*np.pi), ndof-len(fixed)

def _main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--validate", action="store_true")
    args = ap.parse_args()
    if args.validate:
        for N in (16, 32):
            f, nd = run(4, 8, np.radians(90), 210e9, 0.30, 7800, N, N, nmodes=7, clamp="both")
            print("C-C check N=%d:" % N, " ".join("%.3f" % x for x in f))
        print("PLANE183 mesh96: 177.736 218.821 299.912 411.818 482.641 527.514 586.336")
    WB = 2480.106584155221
    PR = {0.25: [0.341243, 0.808929, 0.940032], 0.5: [0.120677, 0.347033, 0.529239], 1.0: [0.039466, 0.106529, 0.263049]}
    ST = {0.25: [0.34589, 0.81921, 0.94012], 0.5: [0.12088, 0.35191, 0.52930], 1.0: [0.03948, 0.10556, 0.26129]}
    for a in (0.25, 0.5, 1.0):
        PHI = np.radians(180 * a); fs = []
        for N in (16, 32, 64):
            Nt = N if a <= 0.5 else 2 * N
            f, nd = run(3, 7, PHI, 210e9, 0.35, 7800, N, Nt, nmodes=6, clamp="theta0")
            fs.append(f)
        fs = np.array(fs)
        for k in range(4):
            x, y, z = fs[:, k]
            if (x - y) * (y - z) > 0 and y != z:
                p = math.log(abs(x - y) / abs(y - z)) / math.log(2); finf = z + (z - y) / (2 ** p - 1)
            else:
                p = float("nan"); finf = z
            O = 2 * math.pi * finf / WB
            extra = ("  printed %+.3f%%  Seok-Tiersten %+.3f%%" % ((PR[a][k] - O) / O * 100, (ST[a][k] - O) / O * 100)) if k < 3 else ""
            print("2T/pi=%.2f mode %d  Omega_FE=%.6f  (p=%.2f)%s" % (a, k + 1, O, p, extra))

if __name__ == "__main__":
    _main()
