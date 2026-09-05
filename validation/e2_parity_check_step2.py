# -*- coding: utf-8 -*-
"""
e2_parity_check_step2.py -- standalone SymPy algebra check, step 2 of 2
(see e2_parity_check_step1.py for the CONTEXT on open item E.2).

Step 1 established whether the raw ssg/csg edge functions are odd or even
under a sign flip. This step carries that through the actual clamped-edge
term as literally coded (`Aj*ai - Bj*bi`, the disp_j*stress_i pairing, no
`_orient` multiplier applied), for all four q_i/q_j parity combinations, to
check whether the naive fix ("just multiply by e.sign like the free edge's
B12 fix") would be correct -- i.e. whether the term already self-flips
correctly (ODD in sg) or needs an explicit correction (EVEN/NEITHER).
Finding: same-angular-parity pairs self-flip correctly already: an added
sign multiplier would double-flip them, so the naive fix is wrong. See
`validation/README.md` and the top-level project notes for the follow-up
(deriving the true sign convention from the source paper's Eq. 43) that
this result motivates.
"""
import sympy as sp

xi, Th, rr = sp.symbols('xi Theta rr', positive=True, real=True)
sg = sp.symbols('sg', real=True)
Wj, Wi, Qyi, Qyj, Tyyi, Tyyj = sp.symbols('Wj Wi Qyi Qyj Tyyi Tyyj')

def ssg_csg(q, sgval):
    ph = q*sp.pi/2
    return (sp.sin(sgval*xi*Th + ph), sp.cos(sgval*xi*Th + ph))

def terms(qi, qj, sgval):
    ssg_i, csg_i = ssg_csg(qi, sgval)
    ssg_j, csg_j = ssg_csg(qj, sgval)
    Aj = ssg_j * Wj
    Bj = xi*csg_j * Wj / rr
    ai = xi*csg_i * Qyi
    bi = ssg_i * Tyyi
    return sp.expand_trig(sp.expand(Aj*ai - Bj*bi))

print("CLAMPED term T(sg) = Aj*ai - Bj*bi  (as literally coded, no _orient)")
print("=" * 78)
for qi in (0, 1):
    for qj in (0, 1):
        Tp = sp.simplify(terms(qi, qj, 1))
        Tm = sp.simplify(terms(qi, qj, -1))
        odd = sp.simplify(Tp + Tm) == 0
        even = sp.simplify(Tp - Tm) == 0
        label = "ODD in sg (auto-flips, no _orient needed)" if odd else \
                "EVEN in sg (sg-invariant, does NOT flip)" if even else \
                "NEITHER (mixed / genuinely needs an explicit fix)"
        print(f"q_i={qi} q_j={qj}:  T(+1)={Tp}")
        print(f"           T(-1)={Tm}")
        print(f"           -> {label}\n")
