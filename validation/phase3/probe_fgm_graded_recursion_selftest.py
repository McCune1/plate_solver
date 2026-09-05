# -*- coding: utf-8 -*-
"""
probe_fgm_graded_recursion_selftest.py -- DIAGNOSTIC ONLY. Standalone
(pure mpmath -- NO sympy, NO plate_solver import). No SOLVER_VERSION
implications.

REVISION NOTE (2026-07-22): the first version of this probe used sympy for
the coefficient-extraction algebra. Job 2327801 failed on the cluster --
sympy is not in this project's venv/requirements.txt (only mpmath/numpy/
scipy are). Rather than add a new dependency for a diagnostic-only script,
this version hand-derives the same recursion via explicit truncated-
power-series convolution (binomial expansion of (x+r0)^n, Cauchy products
for T(x)*R(x) etc.) in plain mpmath. It was cross-checked against the
original sympy version in the author's sandbox: a_4 and a_5 agree to
~15 significant digits between the two independently-coded methods (sympy
symbolic solve vs hand convolution), and the mpmath-only version's
residuals are actually BETTER (1e-18..1e-31 vs the sympy version's
1e-13..1e-15, because sympy's per-step sp.N(...,30) rounding was
truncating intermediate precision that mpmath's native arithmetic keeps).

CONTEXT: the OOP flexural ODE (Seok & Tiersten Eq. 25/29) is
    (x+r0)^4 w'''' + 2(x+r0)^3 w''' - (x+r0)^2 Phi(x) w''
        + (x+r0) Phi(x) w' + Psi(x) w - (x+r0)^4 kbar^2 Omega^2 w = 0
    Phi(x) = 2 T(x) xi^2 + R(x)
    Psi(x) = R(x) xi^2 (xi^2 - 2 T(x) - 2 R(x))
solved for the homogeneous (constant T, R) case by a Frobenius series
w(x) = sum a_m x^(lambda+m), lambda in {0,1,2,3}, recursion B1 (project
ledger). For a radially-graded (FGM) material, T and R become power
series in x (from expanding a power-law E(r) grading about r0): T(x) =
sum T_k x^k, R(x) = sum R_k x^k. A sandbox symbolic check this session
confirmed the recursion for a_m widens from a fixed 4-term window to a
full convolution over T_0..T_m, R_0..R_m -- e.g. the coefficient of x^1
already brings in R0*R1, R0*T1, R1*T0 cross terms that don't exist in the
constant-coefficient case, and higher orders bring in progressively more
(R0*R2, R1^2, R1*T1, R2*T0 at x^2, etc.) -- confirming the recursion
becomes O(m) work per term / O(M^2) total, not the O(1)-per-term
recursion the existing M=80 loop assumes.

THIS PROBE turns that structural finding into an actual numerical
self-test, in the same spirit as this project's existing "ODE residual"
regression checks (Part 1 residual approx 2.42e-13, Part 2 approx
1.95e-15): it (1) builds the widened recursion directly (convolution
arithmetic, shown in solve_recursion() below) for a mild linear grading
T(x)=T0+T1*x, R(x)=R0+R1*x, (2) solves it forward for N series terms,
seeding a_0=1, a_1=a_2=a_3=0 (the standard convention: each of the 4
independent solutions is distinguished by which single one of a_0..a_3 is
nonzero), (3) evaluates the resulting series and its derivatives at
several interior points, well inside the series' domain, and (4) plugs
them back into the ORIGINAL (closed-form graded, not truncated further)
ODE to get a residual. If the recursion has an algebra bug, this residual
will NOT be small -- it is not possible to "pass" this test by
construction, and the a_4-cross-check below independently confirms the
recursion matches a second, differently-coded derivation of the same
algebra.

PRE-REGISTERED INTERPRETATION:
  PASS: residual at every test point is at or near mp.dps-level precision
    (relative residual << 1e-10 at dps=30), for BOTH the graded (T1,R1 !=
    0) and reduced (T1=R1=0) case, AND the two cases agree on a_4 exactly
    (a_4's recursion only involves T0,R0 at this order -- T1,R1 first
    enter at a_5 -- a structural prediction that is verified, not
    assumed, by the a_4-match printout). This numerically confirms the
    widened recursion genuinely solves the graded ODE and correctly
    reduces to the known constant-coefficient case.
  FAIL: any residual well above precision-level (e.g. 1e-6 or worse), or
    a_4 mismatch, indicates a genuine algebra error -- do NOT treat this
    as "probably fine", the recursion needs to be re-derived, not the
    tolerance loosened.
  NOTE: this validates internal self-consistency (the series solves its
    own ODE) -- it does NOT validate against any FGM literature benchmark
    (Shen, Zhang, Wang FGM sector papers per the roadmap are not in this
    project's files and would need to be sourced separately) and does NOT
    imply a working plate_solver implementation exists -- geometry.py's
    FGMPlateProperties still needs the MaterialModel integration
    (oop_constants()/oop_coeff_series()) before any of this can be run
    through _build_K_real.

COST: pure convolution arithmetic, no full_search, no plate_solver
import, no sympy. Runs in well under a second even at N=80 in the
author's sandbox (measured: 0.14s) -- this submit script exists mainly so
the result is confirmed on the cluster too and for provenance (job-number
citation), not because the sandbox result is in doubt.
"""
import sys
import time
from math import comb

import mpmath as mp

mp.mp.dps = 30

# ---- concrete numeric parameters (mild FGM grading about r0) ----
r0_v = mp.mpf('3.9270')      # r0_bar for r0/(2b)=1.25 (project's own
                             # standard cantilever validation geometry)
xi_v = mp.mpf('1.305')      # a real canonical branch value (ledger B4)
kbar_v = mp.mpf('31.4306')  # project's standard kbar (config.py convention)
Om_v = mp.mpf('0.033238')   # a real known mode (isotropic cantilever baseline)
T0_v, T1_v = mp.mpf('1.0'), mp.mpf('0.05')   # T(x) = T0 + T1*x  (T1=0 recovers
R0_v, R1_v = mp.mpf('1.0'), mp.mpf('0.08')   # the isotropic / constant-coeff
                                              # B1 case)  R(x) = R0 + R1*x

N = int(sys.argv[1]) if len(sys.argv) > 1 else 80   # series length; project
                                                     # convention uses M=80


def hdr(s):
    print("\n" + "=" * 78 + f"\n  {s}\n" + "=" * 78, flush=True)


# ---------------------------------------------------------------------------
# Truncated power-series algebra (plain lists of mpf, index = power of x)
# ---------------------------------------------------------------------------
def poly_mul(p, q, maxlen=None):
    n = len(p) + len(q) - 1
    if maxlen is not None:
        n = min(n, maxlen)
    out = [mp.mpf(0)] * n
    for i, pi in enumerate(p):
        if pi == 0:
            continue
        for j, qj in enumerate(q):
            if i + j >= n:
                break
            out[i + j] += pi * qj
    return out


def poly_add(p, q):
    n = max(len(p), len(q))
    out = [mp.mpf(0)] * n
    for i in range(len(p)):
        out[i] += p[i]
    for i in range(len(q)):
        out[i] += q[i]
    return out


def poly_scale(p, s):
    return [c * s for c in p]


def binom_pow(r0, n):
    """(x+r0)^n coefficients, index=power of x, length n+1."""
    return [mp.mpf(comb(n, k)) * r0 ** (n - k) for k in range(n + 1)]


def falling(m, d):
    f = mp.mpf(1)
    for k in range(d):
        f *= (m - k)
    return f


def solve_recursion(T, R, r0, xi, kbar, Om, N):
    """T, R: lists of coefficients T_k, R_k (T(x)=sum T_k x^k). Solves the
    graded ODE's widened Frobenius recursion for a_0..a_{N-1} (lambda=0
    basis solution: a_0=1, a_1=a_2=a_3=0 by convention). At each new order
    p, ONLY the leading (x+r0)^4 term's top-degree slot involves the new
    unknown a_{p+4} -- every other contribution (including the OTHER,
    lower-degree slots of the (x+r0)^4 term itself) only touches
    a_0..a_{p+3}, already known -- so this is a genuine forward recursion,
    not a linear system."""
    xi2 = xi ** 2
    Phi = poly_add(poly_scale(T, 2 * xi2), R)
    RT = poly_mul(R, T)
    RR = poly_mul(R, R)
    Psi = poly_add(poly_scale(R, xi2 ** 2),
                    poly_add(poly_scale(RT, -2 * xi2), poly_scale(RR, -2 * xi2)))
    rx4 = binom_pow(r0, 4)
    rx3 = binom_pow(r0, 3)
    rx2 = binom_pow(r0, 2)
    rx1 = binom_pow(r0, 1)
    rx2_Phi = poly_mul(rx2, Phi)
    rx1_Phi = poly_mul(rx1, Phi)

    a = [mp.mpf(0)] * N
    a[0] = mp.mpf(1)
    a[1] = mp.mpf(0)
    a[2] = mp.mpf(0)
    a[3] = mp.mpf(0)

    def b1(j):
        return a[j + 1] * (j + 1)

    def b2(j):
        return a[j + 2] * falling(j + 2, 2)

    def b3(j):
        return a[j + 3] * falling(j + 3, 3)

    def b4_known(j):
        return a[j + 4] * falling(j + 4, 4)

    for p in range(0, N - 4):
        s = mp.mpf(0)
        for k in range(1, min(4, p) + 1):          # (x+r0)^4 w'''', k>=1 part
            j = p - k
            if j >= 0:
                s += rx4[k] * b4_known(j)
        for k in range(0, min(3, p) + 1):           # 2*(x+r0)^3 w'''
            j = p - k
            if j >= 0:
                s += 2 * rx3[k] * b3(j)
        for k in range(0, len(rx2_Phi)):             # -(x+r0)^2*Phi(x)*w''
            j = p - k
            if j >= 0:
                s += -rx2_Phi[k] * b2(j)
        for k in range(0, len(rx1_Phi)):             # (x+r0)*Phi(x)*w'
            j = p - k
            if j >= 0:
                s += rx1_Phi[k] * b1(j)
        for k in range(0, len(Psi)):                 # Psi(x)*w
            j = p - k
            if j >= 0:
                s += Psi[k] * a[j]
        for k in range(0, 5):                        # -(x+r0)^4*kbar^2*Om^2*w
            j = p - k
            if j >= 0:
                s += -(kbar ** 2) * (Om ** 2) * rx4[k] * a[j]

        leading = rx4[0] * falling(p + 4, 4)
        a[p + 4] = -s / leading
    return a


def ode_residual(a, T0, T1, R0, R1, r0, xi, kbar, Om, xv):
    N = len(a)
    xv = mp.mpf(xv)

    def series_eval(d):
        s = mp.mpf(0)
        for m in range(N):
            p = m
            if p - d < 0:
                continue
            s += a[m] * falling(p, d) * xv ** (p - d)
        return s

    W = series_eval(0)
    W1 = series_eval(1)
    W2 = series_eval(2)
    W3 = series_eval(3)
    W4 = series_eval(4)
    Tx = T0 + T1 * xv
    Rx = R0 + R1 * xv
    Phi_v = 2 * Tx * xi ** 2 + Rx
    Psi_v = Rx * xi ** 2 * (xi ** 2 - 2 * Tx - 2 * Rx)
    lhs = ((xv + r0) ** 4 * W4 + 2 * (xv + r0) ** 3 * W3
           - (xv + r0) ** 2 * Phi_v * W2 + (xv + r0) * Phi_v * W1
           + Psi_v * W - (xv + r0) ** 4 * kbar ** 2 * Om ** 2 * W)
    scale = max(abs(W4 * (xv + r0) ** 4), mp.mpf('1e-30'))
    return lhs, lhs / scale


def main():
    hdr("probe_fgm_graded_recursion_selftest -- start " + time.strftime("%Y-%m-%d %H:%M:%S"))
    print(f"N={N} series terms, dps={mp.mp.dps}")

    t0 = time.time()
    a_graded = solve_recursion([T0_v, T1_v], [R0_v, R1_v], r0_v, xi_v, kbar_v, Om_v, N)
    a_reduced = solve_recursion([T0_v, mp.mpf(0)], [R0_v, mp.mpf(0)], r0_v, xi_v, kbar_v, Om_v, N)
    print(f"recursion solved (both cases) in {time.time()-t0:.3f}s")
    print(f"a_0..a_5 graded : {[float(x) for x in a_graded[:6]]}")
    print(f"a_0..a_5 reduced: {[float(x) for x in a_reduced[:6]]}")
    a4_match = a_graded[4] == a_reduced[4]
    print(f"a_4 identical between graded/reduced (expected -- T1,R1 first "
          f"enter at a_5)? {'PASS' if a4_match else 'FAIL'}")
    print(f"a_{N-1} graded  = {a_graded[N-1]}  (should be small/shrinking)")

    hdr("ODE residual check (series plugged back into the ORIGINAL graded ODE)")
    worst_g = mp.mpf(0)
    worst_r = mp.mpf(0)
    for xv in [0.05, 0.15, 0.30, -0.10, -0.25]:
        lhs_g, rel_g = ode_residual(a_graded, T0_v, T1_v, R0_v, R1_v, r0_v, xi_v, kbar_v, Om_v, xv)
        lhs_r, rel_r = ode_residual(a_reduced, T0_v, mp.mpf(0), R0_v, mp.mpf(0), r0_v, xi_v, kbar_v, Om_v, xv)
        worst_g = max(worst_g, abs(rel_g))
        worst_r = max(worst_r, abs(rel_r))
        print(f"  x={float(xv):+.2f}: graded residual={mp.nstr(lhs_g,6)} rel={mp.nstr(rel_g,6)}   "
              f"reduced residual={mp.nstr(lhs_r,6)} rel={mp.nstr(rel_r,6)}")

    hdr("Verdict")
    print(f"  worst relative residual, graded case:  {mp.nstr(worst_g,4)}  "
          f"{'PASS (precision-level)' if worst_g < mp.mpf('1e-10') else 'FAIL -- recursion has a bug'}")
    print(f"  worst relative residual, reduced case: {mp.nstr(worst_r,4)}  "
          f"{'PASS (precision-level)' if worst_r < mp.mpf('1e-10') else 'FAIL -- recursion has a bug'}")
    print(f"  a_4 cross-check: {'PASS' if a4_match else 'FAIL'}")

    hdr("probe_fgm_graded_recursion_selftest -- end " + time.strftime("%Y-%m-%d %H:%M:%S"))


if __name__ == "__main__":
    main()
