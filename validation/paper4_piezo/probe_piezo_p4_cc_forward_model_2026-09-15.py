"""
probe_piezo_p4_cc_forward_model_2026-09-15.py

Paper 4 (piezoelectric ring) -- standalone, non-package forward-model check
for the C-C (clamped-clamped) piezo-coupled annular plate, using Duan2005's
own closed-form Bessel-family solution directly (mpmath
besselj/bessely/besseli/besselk) instead of this project's Frobenius
power-series machinery -- valid because the per-branch radial ODE
(Delta w = lambda w) is exactly the ordinary/modified Bessel equation of
order n for an ISOTROPIC (non-FGM) material, which is what Duan2005's own
benchmark uses.

STATUS 2026-09-15 (updated, same day, twice): VALIDATED, AND
CLUSTER-CONFIRMED. The h1**2-vs-h1**4 bug described below (found and fixed
this session) is what was causing the ~1% low discrepancy reported earlier
today; six independent Table 4 points matched to <0.05% in the sandbox
(see VALIDATION TABLE below), and the cluster gate built from this script
(Paper4_Piezo/probe_piezo_p4_cc_table4_gate_2026-09-15.py, job 2489140)
then reproduced ALL 27 non-trivial CPT/r0-h=60 Table 4 points at production
precision (dps=100): PASS_ALL, worst rel. err -0.048%, mean abs rel. err
0.021%. See LESSONS_LEARNED.md Sec 18.148 for the full cluster result.
This model is validated for its intended purpose (C-C, as the vehicle for
checking the new electroelastic physics against a literature table) --
folding it into the plate_solver package as a new solver class is a
separate, deliberate architecture decision, not yet done (see
PAPER4_PIEZO_ROADMAP.md's "6x6 not 4x4" framing).

WHAT IS VALIDATED:
- The bare host C-C ring (h1=0, no piezo layer) reproduces Duan2005 Table 4's
  own h1=0 column (p=0 circumferential, n=0,1,2 radial mode index:
  2718, 7520, 14783 rad/s) to within ~0.2% (see companion script
  probe_piezo_p4_elastic_cc_baseline_2026-09-15.py). r_i=0.1, r_o=0.6,
  h=0.01 (Duan2005's own half-thickness symbol -- r0/h=60 uses THEIR h,
  not "full thickness").
- The pure-ELASTIC effect of adding the piezo layer's own bending stiffness
  (d2, with e31_bar forced to 0) reproduces Table 4's h1/2h=1/12 row
  (2792 rad/s) to within 0.01% (companion script). Confirms d1, d2, A2.
- THE FULL 6-DOF THREE-BRANCH COUPLED MODEL (this script), after the fix
  below, reproduces six independent Table 4 points to <0.05%:

  VALIDATION TABLE (2026-09-15, dps=200, extraprec=800, C-C, CPT/r0/h=60):
    p=0 (n_bessel=0), n_radial=0, h1/2h=1/12: computed=2791.79  target=2792   rel_err=-0.0076%
    p=0 (n_bessel=0), n_radial=1, h1/2h=1/12: computed=7722.71  target=7723   rel_err=-0.0038%
    p=0 (n_bessel=0), n_radial=2, h1/2h=1/12: computed=15181.13 target=15182  rel_err=-0.0057%
    p=0 (n_bessel=0), n_radial=0, h1/2h=1/8:  computed=2852.20  target=2853   rel_err=-0.0282%
    p=0 (n_bessel=0), n_radial=0, h1/2h=1/5:  computed=2987.74  target=2989   rel_err=-0.0422%
    p=1 (n_bessel=1), n_radial=0, h1/2h=1/12: computed=2928.22  target=2928   rel_err=+0.0076%

  (n_radial here is Duan2005's Table 4 "n" column, the bisection bracket
  index -- NOT the Bessel order argument passed to cc_det(), which is the
  circumferential index p; the p=0/n_radial=0/h1/2h=1/12 row is this
  script's original __main__ target.) Spans three radial modes, two
  circumferential orders, and three piezo-thickness ratios -- this is a
  multi-point cross-check, not a single-point coincidence.

THE BUG (found and fixed 2026-09-15, same session as the "NOT YET
VALIDATED" status this docstring replaced): cubic_coeffs() had every A2*w2
(inertial / added-mass) term scaled by h1**2 instead of h1**4. Root cause:
Duan2005's own Eq. 17a (p.125, re-OCR'd and hand/sympy cross-multiplied
directly from the printed Eq. 16a numerator over Eq. 16a's own
denominator, set equal to Eq. 16b's lambda(chi), i.e. literally
"mu/chi = lambda", NOT "mu = lambda**2" as this project first assumed) has
h1**2 already present in Eq. 16a's own A2*w2*h1**2*(...) term, and cross-
multiplying by Eq. 16b's own h1**2(...) denominator to clear fractions
introduces a SECOND factor of h1**2 -- i.e. h1**4 overall on every
resulting cubic-in-chi coefficient that traces back to the inertial term,
while the non-inertial (pure piezoelectric/elastic stiffness) terms only
ever picked up h1**3 from a single h1 factor already present in those
terms plus no second multiplication. This was verified two independent
ways: (1) direct hand cross-multiplication of the OCR'd Eq. 16a/16b/17a
text (LESSONS_LEARNED.md's dated entry for today has the full OCR
transcription); (2) sympy symbolic substitution of the (unchanged, already
independently verified against Eq. 16b) branches() lambda(chi) formula
into the literal Eq. 14a relation, which reproduces the corrected cubic's
chi**3 term and the non-inertial part of its chi**2 term bit-for-bit
(sign included), isolating the discrepancy to exactly the A2*w2 terms.
Before the fix, this produced a ~1% LOW root (2765 vs 2792 target) despite
the elastic-only and pure-added-stiffness sub-models each matching to
<0.2% -- the bug was invisible to those simpler checks because neither one
exercises cubic_coeffs() at all.

RULED OUT during this session's diagnostic process (kept here so the next
person doesn't re-walk the same dead ends): wrong C-C boundary condition
choice (Duan2005's own Eq. 23a,b confirms w=w'=phi'=0 at each edge is
exactly what this script already implements); a Bessel-derivative
recurrence-formula bug (cross-checked via mp.diff against all four
formulas, exact match); insufficient mpmath precision as the SOLE cause
(dps=250 and dps=600 pre-fix both converged to the same wrong ~2765 root,
which is consistent with a real algebraic error, not precision noise --
this was suggestive but not conclusive on its own; the sympy cross-check
above is what actually nailed it down).

STILL OPEN / NOT YET DONE:
- The independent from-scratch ODE-shooting cross-check (originally
  next-step #2) was attempted but failed from double-precision overflow
  (branch 3's exponential growth over the ri-to-ro span); not necessary
  now that the sympy algebraic re-derivation has independently confirmed
  the fix, but would still be a nice-to-have belt-and-suspenders check
  before a paper submission, if time permits.
- Only C-C has been checked. F-F (the user's parallel-track BC) still
  needs its own boundary-row derivation (Duan2005 Eq. 25, free-edge:
  M_rr=Q=phi'=0) and its own validation pass before trusting it -- do NOT
  assume the h1**4 fix "obviously" carries over correctly to F-F's
  det-matrix construction without re-deriving that matrix's rows the same
  careful way this session did for C-C.
- Table 4's IPT-based (Mindlin, shear-deformable) rows were not checked
  against (this script is CPT/Kirchhoff only, matching what
  LESSONS_LEARNED already scopes this piezo work to).
- No worker-path / package-integration regression test has been run yet
  -- this remains a standalone script, not yet folded into plate_solver/.
  Per project invariants, do that (and the accompanying worker-path
  re-test) before any cluster submission.

Do NOT hand this off to Grok as a "quick probe" per the project's own
hand-off convention -- this needed the actual boundary-condition/
derivation context (LESSONS_LEARNED Sec 18.141/18.147), which a fresh Grok
session would not have without a much longer hand-off brief than usual.
"""
import mpmath as mp

mp.mp.dps = 300

pi = mp.pi

# Material constants (Duan2005 Table 1, PZT4 piezo layer + steel host)
C11E = mp.mpf('132e9'); C13E = mp.mpf('73e9'); C33E = mp.mpf('115e9')
e31 = mp.mpf('4.1'); e33 = mp.mpf('14.1')
X11 = mp.mpf('7.124e-9'); X33 = mp.mpf('5.841e-9')
rho_pzt = mp.mpf('7500')
E_steel = mp.mpf('200e9'); nu_steel = mp.mpf('0.3'); rho_steel = mp.mpf('7800')

# Geometry: Duan2005 Table 4 target (r1=0.1, r0=0.6, r0/h=60 with h = Duan's
# own half-thickness symbol, i.e. h=0.01 here, NOT full host thickness).
ri = mp.mpf('0.1'); ro = mp.mpf('0.6'); h = mp.mpf('0.01')

c11_bar = C11E - C13E**2 / C33E
e31_bar = e31 - (C13E / C33E) * e33
Xi33_bar = X33 + e33**2 / C33E
Xi11_bar = X11
d1 = mp.mpf(2) / 3 * E_steel * h**3 / (1 - nu_steel**2)


def cubic_coeffs(omega, h1, d2):
    """Duan2005 Eq. 17a cross-multiplied into a*chi^3+b*chi^2+c*chi+d=0.
    FIXED 2026-09-15: every A2*w2 (inertial) term needs h1**4, not h1**2
    -- see this file's module docstring "THE BUG" section for the full
    derivation and the sympy cross-check that pinned it down. The chi**3
    term and the non-inertial part of the chi**2 term were already correct
    and are unchanged."""
    w2 = omega**2
    A2 = 2 * (rho_steel * h + rho_pzt * h1)
    a = -16 * pi * Xi11_bar * Xi33_bar * e31_bar * h1**3
    b = (4*A2*w2*h1**4*Xi11_bar**2 - 8*pi**2*Xi33_bar*e31_bar**2*h1**3
         - 4*pi**4*Xi33_bar**2*(d1 + d2))
    c = 4 * pi * A2 * w2 * h1**4 * Xi11_bar * e31_bar
    d = pi**2 * A2 * w2 * h1**4 * e31_bar**2
    return a, b, c, d


def branches(omega, h1):
    d2 = mp.mpf(2) / 3 * c11_bar * ((h + h1)**3 - h**3)
    a, b, c, d = cubic_coeffs(omega, h1, d2)
    roots = mp.polyroots([a, b, c, d], maxsteps=300, extraprec=1500)
    out = []
    for chi in roots:
        lam = 2*pi**2*Xi33_bar*chi / (h1**2*(2*Xi11_bar*chi + e31_bar*pi))
        out.append((chi, lam))
    return out


def radial_quad(n, r, lam):
    """Z1, Z2, dZ1/dr, dZ2/dr for Delta w = lam*w, per Duan2005 Eq. 20a-d
    (I_p,K_p if lam>0; J_p,Y_p if lam<0), derivatives via the standard
    Bessel recurrences (kept analytic/exact, not numerical differencing;
    cross-checked against mp.diff() this session, exact match)."""
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


def cc_det(omega, n, h1):
    """6x6 C-C boundary determinant: per edge, w=0, w'=0, phi'=0
    (short-circuit ansatz + insulated radial edge, Sec 18.141 F) --
    confirmed 2026-09-15 directly against Duan2005 Eq. 23a,b (own C-C
    definition), exact match. Two-sided (row+column) equilibration before
    mp.det() to combat the extreme dynamic range from the near-pure-
    dielectric branch's large radial decay constant."""
    lams = branches(omega, h1)

    def rows_at(r):
        w_row, wp_row, phip_row = [], [], []
        for chi, lam in lams:
            Z1, Z2, dZ1, dZ2 = radial_quad(n, r, lam)
            w_row += [Z1, Z2]
            wp_row += [dZ1, dZ2]
            phip_row += [chi*dZ1, chi*dZ2]
        return w_row, wp_row, phip_row

    w_ri, wp_ri, phip_ri = rows_at(ri)
    w_ro, wp_ro, phip_ro = rows_at(ro)
    M = mp.matrix([w_ri, wp_ri, phip_ri, w_ro, wp_ro, phip_ro])
    for j in range(6):
        mx = max(abs(M[i, j]) for i in range(6)) or mp.mpf(1)
        for i in range(6):
            M[i, j] = M[i, j] / mx
    for i in range(6):
        mx = max(abs(M[i, j]) for j in range(6)) or mp.mpf(1)
        for j in range(6):
            M[i, j] = M[i, j] / mx
    return mp.det(M)


def bisect(lo, hi, n, h1, iters=40):
    flo = cc_det(mp.mpf(lo), n, h1)
    for _ in range(iters):
        mid = (lo + hi) / 2
        fm = cc_det(mp.mpf(mid), n, h1)
        if (fm.real > 0) == (flo.real > 0):
            lo, flo = mid, fm
        else:
            hi = mid
    return (lo + hi) / 2


if __name__ == "__main__":
    # Duan2005 Table 4, h1/2h=1/12, p=0 (axisymmetric), n=0 (fundamental
    # radial mode): target 2792 rad/s. Post-fix this brackets 2791.79.
    h1 = 2 * h / 12
    for om in (2718, 2765, 2780, 2792, 2810):
        d = cc_det(mp.mpf(om), 0, h1)
        print(om, '+' if d.real > 0 else '-', float(mp.log(abs(d) + mp.mpf('1e-600'))))
    root = bisect(2780, 2800, 0, h1, iters=30)
    print('bisected root =', float(root), ' target=2792')
