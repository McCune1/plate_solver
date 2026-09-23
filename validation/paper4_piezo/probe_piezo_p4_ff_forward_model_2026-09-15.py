# WARNING 2026-09-22 (LESSONS_LEARNED Sec 18.226): this historical forward-model
# probe uses Q_r alone as the free-edge shear row. The correct row is the
# Kirchhoff effective shear, which adds -2*A1*n^2*(Z'/r^2 - Z/r^3). The two
# agree exactly at n=0 (every published use of this script). Its n>=1
# outputs, e.g. the p=1 root 1573.49 rad/s, are WRONG: the corrected value
# is 1756.01. Use plate_solver.piezo_solver, which is fixed, for n>=1.
"""
probe_piezo_p4_ff_forward_model_2026-09-15.py

Paper 4 (piezoelectric ring) -- standalone, non-package forward-model check
for the F-F (free-free) piezo-coupled annular plate -- the user's actual
primary boundary condition target per PAPER4_PIEZO_ROADMAP.md (C-C was only
ever the validation vehicle, chosen because it is the one BC with a literal
published number table, Duan2005 Table 4). Mirrors
probe_piezo_p4_cc_forward_model_2026-09-15.py's architecture exactly
(same material constants, same geometry, same cubic_coeffs()/branches()/
radial_quad() -- copied byte-for-byte, UNCHANGED) -- only the boundary-row
construction differs, per Duan2005 Eq. (25a,b): free edge is
M_rr = Q = d(phi)/dr = 0 at each edge, replacing C-C's w = w' = phi' = 0.

STATUS 2026-09-15: DERIVED AND SANDBOX-VALIDATED (internally, no literature
target -- NOT yet FE-confirmed). No literature F-F piezo table exists
(confirmed LESSONS_LEARNED.md Sec 18.144 Finding 2 -- neither Duan2005 nor
Liu2002 tabulates a free edge case), so validation here is by (1) exact
algebraic cross-check against Duan2005's own printed M_rr/Q_r formulas
(Eqs. 9a, 10a, re-verified against a 400 DPI rendered page image, not OCR
alone -- see "THE M_rr SIGN BUG" section below for what that caught), (2)
an exact structural match against this project's own already-literature-
validated (Paper 3, 704/710 vs Vogel/Leissa/Narita/Irie) Seok-Tiersten M_r
row formula (OutOfPlaneSolver._Lmat_mp) in the isotropic single-layer
limit, (3) a live numerical elastic-limit regression against ring_disk.py's
own F-F solver for this exact geometry/material, matching to 5e-8 relative,
and (4) physical-plausibility checks on the full coupled model (below). A
real FE (ANSYS SOLID226/PLANE223) cross-check of the FULL piezo-coupled
F-F number is still needed before this is "validated" in the same sense
the C-C model now is (see the sibling ANSYS deck this session wrote,
ansys_p4_ff_piezo_smoke_2026-09-15.inp, not yet run at the time of
writing) -- do not cite these F-F piezo numbers in the paper until that FE
check exists and passes.

RESULTS TABLE (this session, dps=100, p=0 axisymmetric unless noted;
computed via bisect() below, each a genuine sign-flip root, not a dip):

  h1/2h   elastic-only(h1=0)  elastic-bilayer(e31=0)  full 3-branch coupled  coupling-only delta
  1/12    726.8119 rad/s      747.9864 rad/s          747.9929 rad/s        +0.0087%
  1/8     726.8119 rad/s      764.1133 rad/s          764.1336 rad/s        +0.0027%
  1/5     726.8119 rad/s      800.3027 rad/s          800.3738 rad/s        +0.0089%

(elastic-only and elastic-bilayer columns from the sibling script
probe_piezo_p4_elastic_ff_baseline_2026-09-15.py; "coupling-only delta" is
the full-coupled root relative to the elastic-bilayer root, i.e. the pure
electromechanical-coupling contribution with the added-layer stiffness
effect already subtracted out.) Two plausibility checks, not proof, since
there is no literature number: (a) the added-layer stiffening (elastic-
bilayer vs elastic-only column) is 2.9%/5.1%/10.1% across the three
ratios -- close to, and the same order of magnitude as, Duan2005 Table 4's
own C-C added-stiffness percentages for these exact ratios (2.70%/4.93%/
9.95%), which is what physical intuition predicts (adding the same
material in the same proportion should give a similar-magnitude relative
stiffness increase regardless of which mechanical edge condition is
applied); (b) the coupling-only delta (0.003-0.009%) is the same order of
magnitude as Sec 18.141 F's short-circuit prediction (Duan2005 Fig. 3:
0.01-0.016% for C-C) -- coupling is a small correction on top of ordinary
added-layer stiffness under short-circuit for F-F too, as expected from
the same physics (the electrical BC, not the mechanical one, is what
determines this). A p=1 (n=1) run at h1/2h=1/12 also finds a genuine
sign-flip root (1573.49 rad/s, far from p=0's ~748 -- not a near-degenerate
coincidence), confirming the 6x6 determinant behaves sensibly across
circumferential order, not just at the one point tabulated above.

THE M_rr SIGN BUG (found this session, in LESSONS_LEARNED.md Sec 18.141 C's
own transcription, NOT in any deployed code -- C-C never exercises M_rr/Q_r
at all, since its BC is w=w'=phi'=0, so this bug was latent and harmless
until now): Sec 18.141 C transcribed
    M_rr = -(d1+d2)*d2w/dr2 + (d1+d2-2A1)*(...) + (4/pi)*h1*e31_bar*phi
i.e. with the SECOND and THIRD terms carrying a PLUS sign relative to the
first. The actual printed equation (Duan2005 Eq. 9a, p.123, re-rendered at
400 DPI this session and read directly off the image, not OCR) is
    M_rr = -[ (d1+d2)*d2w/dr2 + (d1+d2-2A1)*(dw/r dr + d2w/r2 dtheta2)
              + (4/pi)*h1*e31_bar*phi ]
i.e. ALL THREE terms are inside one outer bracket with a single leading
minus sign -- the second and third terms are NEGATIVE relative to the
first, not positive. (M_theta_theta has the identical bracket structure
and was mis-transcribed the same way; M_r_theta and Q_r, Q_theta were
transcribed correctly -- their brackets were already complete in Sec
18.141 C.) This is a genuine sign error in the derivation notes, not a
tuning issue: caught by directly re-deriving M_rr as it must appear in a
FREE-edge boundary row (this script), where getting the relative sign
between the curvature term and the (w'/r - n^2 w/r^2)/piezo terms wrong
would silently build the WRONG boundary condition while still returning
some frequency (a determinant of six columns is never accidentally
identically zero for a wrong-but-nonsingular matrix, so this class of bug
does not announce itself the way a NaN or an all-same-sign residual does
-- it has to be caught by cross-checking the formula itself, which is
exactly what this session did before writing any root-search code).
Verified two independent ways before trusting the corrected formula (see
"CROSS-CHECKS" below): (a) direct algebraic elimination of w'' via the
branch ODE reduces the corrected M_rr row to EXACTLY this project's own
already-literature-validated Seok-Tiersten M_r row formula
(OutOfPlaneSolver._Lmat_mp) in the isotropic single-layer limit -- an
exact symbolic match, not an approximate one; (b) a live numerical run
(below) reproduces ring_disk.py's own F-F elastic root for this exact
geometry to 5e-8 relative.

DERIVATION -- reducing M_rr/Q_r to closed form in Z, dZ/dr only (no second
or third derivatives needed at all): each branch's radial function Z
satisfies the branch ODE Delta(Z) = lambda_i * Z exactly, i.e.
    Z'' + Z'/r - p^2*Z/r^2 = lambda_i * Z   =>   Z'' = lambda_i*Z - Z'/r + p^2*Z/r^2.
Substituting this into the (corrected) M_rr formula, with phi = chi_i * Z
per branch, and using p for the circumferential Bessel order (Duan2005's
own symbol, this project's "n" argument):

  M_rr-row(i) = [ (d1+d2)*lambda_i + (4/pi)*h1*e31_bar*chi_i + 2*A1*p^2/r^2 ] * Z
                - (2*A1/r) * dZ/dr

  Q_r-row(i)  = [ (d1+d2)*lambda_i + (4/pi)*h1*e31_bar*chi_i ] * dZ/dr

  phi'-row(i) = chi_i * dZ/dr                          (unchanged from C-C;
    Duan2005's own Eq. (22) states d(phi)/dr=0 as the general insulated-
    edge electrical condition BEFORE the C/S/F mechanical split, and it is
    then repeated verbatim as the "b" part of Eqs. 23, 24, AND 25 -- i.e.
    directly confirmed, not inferred, that F-F uses the identical
    electrical row as C-C.)

where A1 = 0.5*[(1-nu_steel)*d1 + (1 - c12_bar/c11_bar)*d2] (Duan2005's own
A1 definition, p.124) -- note A1 needs c12_bar = C12E - C13E**2/C33E, which
the C-C script never needed (C-C's BC never touches A1 at all) and is
added here for the first time; C12E=71e9 per Duan2005 Table 1.

Both M_rr-row and Q_r-row need only Z and dZ/dr, which radial_quad() (below,
copied unchanged from the C-C script) already provides -- no new Bessel-
derivative formula, and no new source of a recurrence-formula bug, is
introduced by this extension.

CROSS-CHECKS (all done before trusting the numbers below):
1. Algebraic: in the single-layer isotropic elastic limit (d2=0, e31_bar=0,
   A1=(1-nu)*d1/2), M_rr-row(lambda) = d1*[lambda*Z + (1-nu)*p^2/r^2*Z -
   (1-nu)/r*dZ], which is EXACTLY (r^2/d1) times OutOfPlaneSolver._Lmat_mp's
   own M_r row (r^2*W'' + nu*(r*W' - p^2*W)) after eliminating W'' the same
   way via the SAME Delta(W)=lambda*W identity -- both reduce to the
   identical closed form. Not an approximate match; the two expressions are
   the same polynomial in Z, dZ, lambda, p, r, nu.
2. Numerical (this session, sandbox, package venv with scipy installed):
   ring_disk.py's own validated F-F elastic ring solver, for R_i=0.1,
   R_o=0.6, h=0.01 (Duan's own half-thickness symbol), steel E=200e9,
   nu=0.3, rho=7800, n=0 (axisymmetric), finds a genuine sign-flip root at
   native Omega=0.6007252606, converting (via geometry._omega_lit) to
   omega = 726.8118579070134 rad/s. This script's OWN independent 2-branch
   elastic (bare host, no piezo layer at all) determinant -- built from the
   corrected M_rr-row/Q_r-row formulas above, using mpmath Bessel functions
   directly, NOT the package's Frobenius series -- finds a sign-flip root
   at omega = 726.8118575520 rad/s. Relative difference: -4.88e-8. This is
   the "elastic-limit regression" PAPER4_PIEZO_ROADMAP.md Sec 8 step 2
   requires before trusting any piezo-coupled number, and it passes at
   essentially machine/bisection precision, not just "close."

STILL OPEN / NOT YET DONE:
- No literature F-F piezo benchmark exists (confirmed twice now, Sec
  18.144); FE (ANSYS SOLID226/PLANE223) is the only external validation
  channel available, and that FE cross-check has not been run yet as of
  this script's writing -- see the sibling ANSYS deck this session also
  wrote. Do not treat the full 3-branch coupled numbers below as more than
  "derived and internally self-consistent" until that FE check exists.
- Only CPT (Kirchhoff) is covered, matching the same scope as the C-C work.
- No worker-path / package-integration regression test -- standalone script
  only, same status as the C-C sibling script.
"""
import mpmath as mp

# dps=100 matches the C-C production cluster gate's own working precision
# (LESSONS_LEARNED.md Sec 18.148, job 2489140: dps=100, extraprec=300,
# PASS_ALL 27/27) -- proven adequate for this exact geometry/material at
# this precision level, not an arbitrary reduction from the C-C sandbox
# script's dps=300 (which was never itself shown to be the minimum needed).
mp.mp.dps = 100

pi = mp.pi

# Material constants (Duan2005 Table 1, PZT4 piezo layer + steel host) --
# identical to the C-C script, PLUS C12E (needed here for c12_bar/A1; C-C's
# BC never touches A1 so C-C never needed it).
C11E = mp.mpf('132e9'); C12E = mp.mpf('71e9'); C13E = mp.mpf('73e9'); C33E = mp.mpf('115e9')
e31 = mp.mpf('4.1'); e33 = mp.mpf('14.1')
X11 = mp.mpf('7.124e-9'); X33 = mp.mpf('5.841e-9')
rho_pzt = mp.mpf('7500')
E_steel = mp.mpf('200e9'); nu_steel = mp.mpf('0.3'); rho_steel = mp.mpf('7800')

# Geometry: same Duan2005 Table 4 target as the C-C script.
ri = mp.mpf('0.1'); ro = mp.mpf('0.6'); h = mp.mpf('0.01')

c11_bar = C11E - C13E**2 / C33E
c12_bar = C12E - C13E**2 / C33E        # new for F-F (needed by A1)
e31_bar = e31 - (C13E / C33E) * e33
Xi33_bar = X33 + e33**2 / C33E
Xi11_bar = X11
d1 = mp.mpf(2) / 3 * E_steel * h**3 / (1 - nu_steel**2)


def cubic_coeffs(omega, h1, d2):
    """Identical to probe_piezo_p4_cc_forward_model_2026-09-15.py's
    cubic_coeffs() -- the chi-cubic (Duan2005 Eq. 17a, h1**4-corrected per
    LESSONS_LEARNED.md Sec 18.147) does not depend on which mechanical edge
    condition is applied; only the boundary rows differ between C-C and F-F."""
    w2 = omega**2
    A2 = 2 * (rho_steel * h + rho_pzt * h1)
    a = -16 * pi * Xi11_bar * Xi33_bar * e31_bar * h1**3
    b = (4*A2*w2*h1**4*Xi11_bar**2 - 8*pi**2*Xi33_bar*e31_bar**2*h1**3
         - 4*pi**4*Xi33_bar**2*(d1 + d2))
    c = 4 * pi * A2 * w2 * h1**4 * Xi11_bar * e31_bar
    d = pi**2 * A2 * w2 * h1**4 * e31_bar**2
    return a, b, c, d


def branches(omega, h1):
    """Identical to the C-C script."""
    d2 = mp.mpf(2) / 3 * c11_bar * ((h + h1)**3 - h**3)
    a, b, c, d = cubic_coeffs(omega, h1, d2)
    roots = mp.polyroots([a, b, c, d], maxsteps=300, extraprec=1500)
    out = []
    for chi in roots:
        lam = 2*pi**2*Xi33_bar*chi / (h1**2*(2*Xi11_bar*chi + e31_bar*pi))
        out.append((chi, lam))
    return out


def radial_quad(n, r, lam):
    """Identical to the C-C script (Z1,Z2,dZ1,dZ2 for Delta w = lam*w)."""
    if lam.real >= 0:
        delta = mp.sqrt(lam); x = delta * r
        Z1 = mp.besseli(n, x); Z2 = mp.besselk(n, x)
        dZ1 = delta * mp.mpf('0.5') * (mp.besseli(n-1, x) + mp.besseli(n+1, x))
        dZ2 = delta * mp.mpf('-0.5') * (mp.besselk(n-1, x) + mp.besselk(n+1, x))
    else:
        delta = mp.sqrt(-lam); x = delta * r
        Z1 = mp.besselj(n, x); Z2 = mp.bessely(n, x)
        dZ1 = delta * mp.mpf('0.5') * (mp.besselj(n-1, x) - mp.besselj(n+1, x))
        dZ2 = delta * mp.mpf('0.5') * (mp.bessely(n-1, x) - mp.bessely(n+1, x))
    return Z1, Z2, dZ1, dZ2


def _A1(d2):
    return mp.mpf('0.5') * ((1 - nu_steel)*d1 + (1 - c12_bar/c11_bar)*d2)


def ff_det(omega, n, h1):
    """6x6 F-F boundary determinant: per edge, M_rr=0, Q_r=0, phi'=0
    (Duan2005 Eq. 25a,b -- free edge; d(phi)/dr=0 confirmed identical to
    C-C's own electrical row, see module docstring). The elastic-bilayer
    (e31_bar=0) sanity check is NOT done by overriding e31_bar in-model here
    -- e31_bar=0 is a genuinely singular limit of the chi-cubic (leading AND
    two subleading coefficients all vanish simultaneously, degenerating
    polyroots' cubic solve to 0=0) -- it is instead a fully independent
    standalone 2-branch construction in the sibling script
    probe_piezo_p4_elastic_ff_baseline_2026-09-15.py, mirroring how the C-C
    side already keeps its elastic baseline in its own separate script."""
    d2 = mp.mpf(2) / 3 * c11_bar * ((h + h1)**3 - h**3)
    lams = branches(omega, h1)
    A1v = _A1(d2)
    K_pref = mp.mpf(4)/pi * h1 * e31_bar  # piezo coefficient on chi in M_rr/Q_r

    def rows_at(r):
        m_row, q_row, phip_row = [], [], []
        for chi, lam in lams:
            Z1, Z2, dZ1, dZ2 = radial_quad(n, r, lam)
            Ki = (d1 + d2)*lam + K_pref*chi
            for Z, dZ in ((Z1, dZ1), (Z2, dZ2)):
                m_row.append((Ki + 2*A1v*n**2/r**2)*Z - (2*A1v/r)*dZ)
                q_row.append(Ki*dZ)
                phip_row.append(chi*dZ)
        return m_row, q_row, phip_row

    m_ri, q_ri, phip_ri = rows_at(ri)
    m_ro, q_ro, phip_ro = rows_at(ro)
    M = mp.matrix([m_ri, q_ri, phip_ri, m_ro, q_ro, phip_ro])
    for j in range(6):
        mx = max(abs(M[i, j]) for i in range(6)) or mp.mpf(1)
        for i in range(6):
            M[i, j] = M[i, j] / mx
    for i in range(6):
        mx = max(abs(M[i, j]) for j in range(6)) or mp.mpf(1)
        for j in range(6):
            M[i, j] = M[i, j] / mx
    return mp.det(M)


def bisect(lo, hi, n, h1, iters=45):
    flo = ff_det(mp.mpf(lo), n, h1)
    for _ in range(iters):
        mid = (lo + hi) / 2
        fm = ff_det(mp.mpf(mid), n, h1)
        if (fm.real > 0) == (flo.real > 0):
            lo, flo = mid, fm
        else:
            hi = mid
    return (lo + hi) / 2


if __name__ == "__main__":
    # Full 3-branch coupled F-F model, p=0, all three Table-4 thickness
    # ratios (no literature target -- see module docstring for what these
    # are checked against instead). Brackets pre-registered from this
    # session's own scan (deep sign-flip dips found by coarse scanning
    # first, standard project discipline).
    targets = [
        (mp.mpf(2)/12, '1/12', 747.9, 748.1),
        (mp.mpf(2)/8,  '1/8',  760, 770),
        (mp.mpf(2)/5,  '1/5',  795, 810),
    ]
    print("=== Full 3-branch coupled F-F model, p=0 (axisymmetric) ===")
    for frac, label, lo, hi in targets:
        h1 = frac * h
        root = bisect(lo, hi, 0, h1, iters=40)
        print(f"  h1/2h={label}: root={float(root):.4f} rad/s")

    print()
    print("=== p=1 (n=1 circumferential), h1/2h=1/12 -- robustness check ===")
    h1 = 2 * h / 12
    root_p1 = bisect(1500, 1650, 1, h1, iters=40)
    print(f"  root={float(root_p1):.4f} rad/s")
