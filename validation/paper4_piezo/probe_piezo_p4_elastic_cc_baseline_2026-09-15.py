"""
probe_piezo_p4_elastic_cc_baseline_2026-09-15.py

Paper 4 (piezoelectric ring) -- VALIDATED building-block checks, done before
attempting the full 3-branch piezo-coupled forward model (see the sibling
script probe_piezo_p4_cc_forward_model_2026-09-15.py, which is NOT yet
validated). Both checks below use closed-form Bessel functions directly
(mpmath), since the per-branch radial ODE for this isotropic material is
exactly the ordinary/modified Bessel equation -- no need for this project's
Frobenius power-series machinery here.

Geometry: r_i=0.1, r_o=0.6 m (Duan2005 Table 4's own annular target,
R_i/R_o=1/6). h=0.01 is DUAN2005's OWN half-thickness symbol (host layer is
|z|<h, thickness 2h) -- their "r0/h=60" uses this h directly, not "full
thickness" (confirmed numerically: using h=0.01 the fundamental C-C root
lands at ~2712 rad/s, matching Table 4's 2718; using h=0.005 -- i.e.
treating 0.01 as "full thickness" and halving it -- gives no root anywhere
near 2718).

CHECK 1 (elastic-only, h1=0, no piezo layer at all): reproduces Duan2005
Table 4's h1=0 column at p=0 (axisymmetric), n=0,1,2 (2718, 7520, 14783
rad/s) to within ~0.2% (2712.4, 7503.1, 14749.4 computed). This ~0.2%
residual gap is small and CONSISTENT IN SIGN/MAGNITUDE across all three
modes -- consistent with a minor transcription-precision issue in Table 1's
published material constants (more decimal places than transcribed) rather
than a modeling error, but is flagged as an open item, not swept under the
rug: re-derive from Duan2005's OWN more precise constants if this ever
needs to be tighter than 0.2%.

CHECK 2 (pure elastic effect of ADDING the piezo layer's bending stiffness,
e31 forced to 0 -- i.e. an ordinary two-material elastic bilayer, ZERO
electromechanical coupling): reproduces Table 4's h1/2h=1/12 row (2792
rad/s) to within 0.01% (2791.76 computed). This is strong, independent
confirmation that d1, d2, and A2 (LESSONS_LEARNED.md Sec 18.141 C's mass/
stiffness formulas) are transcribed and coded correctly, and independently
confirms Sec 18.141 F's finding: under short-circuit, the piezo layer's
frequency effect is completely dominated by this ordinary added-layer
stiffness, with the true electromechanical coupling contributing only a
tiny (~0.01-0.016%, per Duan2005 Fig. 3) correction on top.

These two checks do NOT exercise the chi/lambda cubic (Eq. 17a/17b) or the
3-branch decoupling at all -- they are the "does the non-piezo skeleton of
this model work" sanity gate, analogous to this project's own "isotropic
no-op" regression discipline for orthotropy/FGM. Both pass. The full 6-dof
coupled model built on top of this skeleton does NOT yet pass its own
equivalent check -- see the sibling script.
"""
import mpmath as mp

mp.mp.dps = 30

E_steel = mp.mpf('200e9'); nu_steel = mp.mpf('0.3'); rho_steel = mp.mpf('7800')
C11E = mp.mpf('132e9'); C13E = mp.mpf('73e9'); C33E = mp.mpf('115e9')
rho_pzt = mp.mpf('7500')
ri = mp.mpf('0.1'); ro = mp.mpf('0.6'); h = mp.mpf('0.01')  # Duan2005's own h

c11_bar = C11E - C13E**2 / C33E


def det_at(omega, n, dsum, A2):
    """4x4 C-C boundary determinant for the ordinary (non-piezo) isotropic
    Kirchhoff annulus: w=0, w'=0 at each edge, 2 branches (+-beta^2) from
    the biharmonic factorization Delta^2 w = (A2*omega^2/dsum) w."""
    mu = mp.sqrt(A2 * omega**2 / dsum)
    k = mp.sqrt(mu)

    def row(r):
        x = k * r
        I = mp.besseli(n, x); K = mp.besselk(n, x)
        J = mp.besselj(n, x); Y = mp.bessely(n, x)
        dI = k * mp.mpf('0.5') * (mp.besseli(n-1, x) + mp.besseli(n+1, x))
        dK = k * mp.mpf('-0.5') * (mp.besselk(n-1, x) + mp.besselk(n+1, x))
        dJ = k * mp.mpf('0.5') * (mp.besselj(n-1, x) - mp.besselj(n+1, x))
        dY = k * mp.mpf('0.5') * (mp.bessely(n-1, x) - mp.bessely(n+1, x))
        return [I, K, J, Y], [dI, dK, dJ, dY]

    w_ri, wp_ri = row(ri)
    w_ro, wp_ro = row(ro)
    M = mp.matrix([w_ri, wp_ri, w_ro, wp_ro])
    for i in range(4):
        mx = max(abs(M[i, j]) for j in range(4)) or mp.mpf(1)
        for j in range(4):
            M[i, j] = M[i, j] / mx
    return mp.det(M)


def bisect(lo, hi, n, dsum, A2, iters=50):
    flo = det_at(mp.mpf(lo), n, dsum, A2)
    for _ in range(iters):
        mid = (lo + hi) / 2
        fm = det_at(mp.mpf(mid), n, dsum, A2)
        if (fm > 0) == (flo > 0):
            lo, flo = mid, fm
        else:
            hi = mid
    return (lo + hi) / 2


if __name__ == "__main__":
    print("=== CHECK 1: bare host, h1=0 (Table 4 h1=0 column, p=0/n=0,1,2) ===")
    d1 = mp.mpf(2) / 3 * E_steel * h**3 / (1 - nu_steel**2)
    A2_bare = 2 * rho_steel * h
    for target, lo, hi in [(2718, 2700, 2730), (7520, 7480, 7560), (14783, 14700, 14860)]:
        root = bisect(lo, hi, 0, d1, A2_bare)
        print(f"  target={target}  root={float(root):.4f}  rel_err={float((root-target)/target)*100:.4f}%")

    print("=== CHECK 2: elastic bilayer, e31=0, h1/2h=1/12 (Table 4 row: 2792) ===")
    h1 = 2 * h / 12
    d2 = mp.mpf(2) / 3 * c11_bar * ((h + h1)**3 - h**3)
    A2_layered = 2 * (rho_steel * h + rho_pzt * h1)
    root = bisect(2700, 2900, 0, d1 + d2, A2_layered)
    print(f"  root={float(root):.4f}  target=2792  rel_err={float((root-2792)/2792)*100:.4f}%")
