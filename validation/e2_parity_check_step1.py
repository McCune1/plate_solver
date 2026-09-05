# -*- coding: utf-8 -*-
"""
e2_parity_check_step1.py -- standalone SymPy algebra check, step 1 of 2.
No plate_solver import; pure symbolic algebra, run with `python3
e2_parity_check_step1.py`.

CONTEXT: backs open item E.2 -- the annular OutOfPlaneSolver/InPlaneSolver's
`_build_K_real` generalizes the FREE edge's orientation via `e.sign` but
hardcodes the CLAMPED edge's orientation to -Theta. Currently harmless
(every existing boundary-condition class places CLAMPED at -Theta), but the
assumption must be understood -- not just patched by analogy -- before any
future boundary-condition layout (e.g. a clamped edge at +Theta, or the
still-open rectangular free-free extension) might violate it. The edge
functions ssg=sin(sign*xi*Theta+phase) / csg=cos(sign*xi*Theta+phase) are
used on both the +Theta and -Theta edges with the same `sign` convention as
coded in plate_solver.core_solvers. This step checks, symbolically, whether
ssg/csg are each odd or even under sign flip (sg -> -sg) for both parity
phases (q=0,1) -- i.e. whether swapping which edge is "positive" flips the
term automatically or needs an explicit orientation correction. Step 2
(e2_parity_check_step2.py) carries this through the actual clamped-edge
boundary term.
"""
import sympy as sp

xi, Th = sp.symbols('xi Theta', positive=True, real=True)
sg = sp.symbols('sg', real=True)  # stands in for e.sign in {+1,-1}

# ssg/csg exactly as coded: ssg=sin(sg*xi*Th+ph), csg=cos(sg*xi*Th+ph), ph=q*pi/2
for q, ph_val in [(0, 0), (1, sp.pi/2)]:
    ssg = sp.sin(sg*xi*Th + ph_val)
    csg = sp.cos(sg*xi*Th + ph_val)
    ssg_p = sp.simplify(ssg.subs(sg, 1))
    ssg_m = sp.simplify(ssg.subs(sg, -1))
    csg_p = sp.simplify(csg.subs(sg, 1))
    csg_m = sp.simplify(csg.subs(sg, -1))
    print(f"q={q} (ph={ph_val}):")
    print(f"  ssg(+1)={ssg_p}   ssg(-1)={ssg_m}   "
          f"{'ODD in sg' if sp.simplify(ssg_p+ssg_m)==0 else 'EVEN in sg' if sp.simplify(ssg_p-ssg_m)==0 else 'NEITHER'}")
    print(f"  csg(+1)={csg_p}   csg(-1)={csg_m}   "
          f"{'ODD in sg' if sp.simplify(csg_p+csg_m)==0 else 'EVEN in sg' if sp.simplify(csg_p-csg_m)==0 else 'NEITHER'}")
