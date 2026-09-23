# WARNING 2026-09-22 (LESSONS_LEARNED Sec 18.226): this historical forward-model
# probe uses Q_r alone as the free-edge shear row. The correct row is the
# Kirchhoff effective shear, which adds -2*A1*n^2*(Z'/r^2 - Z/r^3). The two
# agree exactly at n=0 (every published use of this script). Its n>=1
# outputs, e.g. the p=1 root 1573.49 rad/s, are WRONG: the corrected value
# is 1756.01. Use plate_solver.piezo_solver, which is fixed, for n>=1.
"""
probe_piezo_p4_elastic_ff_baseline_2026-09-15.py

Paper 4 (piezoelectric ring) -- VALIDATED elastic-limit building-block
checks for the F-F (free-free) boundary condition, done before trusting the
full 3-branch piezo-coupled F-F model (probe_piezo_p4_ff_forward_model_
2026-09-15.py). Mirrors probe_piezo_p4_elastic_cc_baseline_2026-09-15.py's
role for C-C exactly: prove the non-piezo skeleton first, per this
project's own "no-op on every currently-validated geometry" discipline
(PAPER4_PIEZO_ROADMAP.md Sec 8 step 2) and per this session's own
instructions (elastic limit must be checked BEFORE the piezo-coupled F-F
number is trusted).

Both checks use closed-form Bessel functions directly (mpmath), matching
the ordinary/modified Bessel equation each branch of the isotropic
material's radial ODE reduces to (Delta w = lambda w) -- no need for this
project's Frobenius power-series machinery, exactly as the C-C elastic
baseline script does. The M_rr/Q_r row formulas used here are the
CORRECTED ones derived in probe_piezo_p4_ff_forward_model_2026-09-15.py's
own docstring (read that first -- it documents a real sign bug found this
session in LESSONS_LEARNED.md Sec 18.141 C's own M_rr transcription, fixed
here and cross-checked against this project's own already-validated
Seok-Tiersten M_r row formula, OutOfPlaneSolver._Lmat_mp).

CHECK 1 (bare host, h1=0, no piezo layer at all -- the PURE elastic F-F
ring): cross-checked THIS SESSION, live, against ring_disk.py's own
already-literature-validated (Paper 3, 704/710 vs Vogel/Leissa/Narita/Irie)
F-F ring solver, for the exact geometry this piezo work targets (R_i=0.1,
R_o=0.6, h=0.01 -- Duan2005's own half-thickness symbol -- steel E=200e9,
nu=0.3, rho=7800, n=0/p=0 axisymmetric). ring_disk.py's own machinery
(OutOfPlaneSolver + ring_det_mp, Frobenius series, package venv with scipy
installed to make plate_solver importable) finds a genuine sign-flip root
at native Omega=0.6007252606277491, converting via geometry._omega_lit to
omega = 726.8118579070134 rad/s. THIS script's independent 2-branch
closed-form Bessel construction (below) finds omega = 726.8118575519824
rad/s. Relative difference: -4.88e-8 -- essentially exact agreement (the
tiny residual is bisection/float round-off noise in one or both paths, not
a real discrepancy). This is the single most important check in this
whole F-F thread: it is a LIVE cross-check against Paper 3's own literature
-validated ring solver, not just an internal self-consistency check.

CHECK 2 (pure elastic effect of ADDING the piezo layer's bending stiffness,
e31=0 -- an ordinary two-material elastic bilayer, ZERO electromechanical
coupling -- same role as the C-C elastic baseline's own Check 2): computed
for all three Table-4 thickness ratios (1/12, 1/8, 1/5). No literature
target exists for F-F (Sec 18.144), so this is judged by physical
plausibility against the analogous C-C numbers (Duan2005 Table 4's own
h1=0 -> h1/2h added-stiffness percentages), which it matches in order of
magnitude closely (see results below and the sibling script's own results
table) -- not proof, since F-F genuinely has a different mode shape than
C-C and there is no reason to expect an EXACT match, only a similar-sized
effect from adding the same proportion of stiffer/denser material.

These two checks do NOT exercise the chi/lambda cubic (Eq. 17a/17b) or the
3-branch decoupling at all -- e31_bar=0 is in fact a SINGULAR limit of that
cubic (leading coefficient a, and subleading c,d, are all proportional to
e31_bar and vanish together, degenerating polyroots' cubic solve) -- so
this script deliberately does NOT reuse ff_det()/branches() from the sibling
script; it is a fully independent, from-scratch 2-branch construction,
exactly mirroring how the C-C side keeps ITS elastic baseline in its own
separate script rather than trying to zero out e31_bar inside the 3-branch
model.
"""
import math
import mpmath as mp

mp.mp.dps = 60

E_steel = mp.mpf('200e9'); nu_steel = mp.mpf('0.3'); rho_steel = mp.mpf('7800')
C11E = mp.mpf('132e9'); C12E = mp.mpf('71e9'); C13E = mp.mpf('73e9'); C33E = mp.mpf('115e9')
rho_pzt = mp.mpf('7500')
ri = mp.mpf('0.1'); ro = mp.mpf('0.6'); h = mp.mpf('0.01')  # Duan2005's own h

c11_bar = C11E - C13E**2 / C33E
c12_bar = C12E - C13E**2 / C33E


def det_at(omega, n, dsum, A1, A2):
    """4x4 F-F boundary determinant for the ordinary (non-piezo) isotropic
    Kirchhoff annulus: M_rr=0, Q_r=0 at each edge, 2 branches (lambda=+-mu)
    from the biharmonic factorization Delta^2 w = (A2*omega^2/dsum) w, using
    the CORRECTED row formulas from probe_piezo_p4_ff_forward_model_2026-
    09-15.py's own docstring:
        M_rr-row = [dsum*lambda + 2*A1*n^2/r^2]*Z - (2*A1/r)*dZ/dr
        Q_r-row  = dsum*lambda*dZ/dr
    """
    mu = mp.sqrt(A2 * omega**2 / dsum)
    k = mp.sqrt(mu)

    def rows(r):
        x = k * r
        I = mp.besseli(n, x); K = mp.besselk(n, x)
        J = mp.besselj(n, x); Y = mp.bessely(n, x)
        dI = k * mp.mpf('0.5') * (mp.besseli(n-1, x) + mp.besseli(n+1, x))
        dK = k * mp.mpf('-0.5') * (mp.besselk(n-1, x) + mp.besselk(n+1, x))
        dJ = k * mp.mpf('0.5') * (mp.besselj(n-1, x) - mp.besselj(n+1, x))
        dY = k * mp.mpf('0.5') * (mp.bessely(n-1, x) - mp.bessely(n+1, x))
        cols = [(I, dI, mu), (K, dK, mu), (J, dJ, -mu), (Y, dY, -mu)]
        Mrow, Qrow = [], []
        for Z, dZ, lam in cols:
            Mrow.append((dsum*lam + 2*A1*n**2/r**2)*Z - (2*A1/r)*dZ)
            Qrow.append(dsum*lam*dZ)
        return Mrow, Qrow

    M_ri, Q_ri = rows(ri)
    M_ro, Q_ro = rows(ro)
    M = mp.matrix([M_ri, Q_ri, M_ro, Q_ro])
    for i in range(4):
        s = max(abs(M[i, j]) for j in range(4)) or mp.mpf(1)
        for j in range(4):
            M[i, j] = M[i, j] / s
    return mp.det(M)


def bisect(lo, hi, n, dsum, A1, A2, iters=50):
    flo = det_at(mp.mpf(lo), n, dsum, A1, A2)
    for _ in range(iters):
        mid = (lo + hi) / 2
        fm = det_at(mp.mpf(mid), n, dsum, A1, A2)
        if (fm.real > 0) == (flo.real > 0):
            lo, flo = mid, fm
        else:
            hi = mid
    return (lo + hi) / 2


if __name__ == "__main__":
    print("=== CHECK 1: bare host F-F elastic ring, h1=0 (p=0/n=0, axisymmetric) ===")
    d1 = mp.mpf(2) / 3 * E_steel * h**3 / (1 - nu_steel**2)
    A1_bare = mp.mpf('0.5') * (1 - nu_steel) * d1
    A2_bare = 2 * rho_steel * h
    root_bare = bisect(700, 800, 0, d1, A1_bare, A2_bare)
    ring_disk_target = mp.mpf('726.8118579070134')
    rel = float((root_bare - ring_disk_target) / ring_disk_target) * 100
    print(f"  root={float(root_bare):.10f} rad/s  ring_disk.py target={float(ring_disk_target):.10f}"
          f"  rel_err={rel:.2e}%")

    print()
    print("=== CHECK 2: elastic bilayer, e31=0, all three Table-4 ratios (p=0) ===")
    for frac, label, lo, hi in [
        (mp.mpf(2)/12, '1/12', 700, 800),
        (mp.mpf(2)/8,  '1/8',  700, 800),
        (mp.mpf(2)/5,  '1/5',  700, 850),
    ]:
        h1 = frac * h
        d2 = mp.mpf(2) / 3 * c11_bar * ((h + h1)**3 - h**3)
        dsum = d1 + d2
        A1v = mp.mpf('0.5') * ((1 - nu_steel)*d1 + (1 - c12_bar/c11_bar)*d2)
        A2v = 2 * (rho_steel*h + rho_pzt*h1)
        root = bisect(lo, hi, 0, dsum, A1v, A2v)
        pct = float((root - root_bare) / root_bare) * 100
        print(f"  h1/2h={label}: root={float(root):.4f} rad/s  stiffening vs h1=0: {pct:.3f}%")
