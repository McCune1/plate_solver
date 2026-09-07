#!/usr/bin/env python3
"""
Sandbox FE/solver check of the 2026-09-03 closed-contour 4-corner derivation
against FIX and the deployed (identically-zero) CURRENT term.

Does NOT touch plate_solver/core_solvers.py. No SOLVER_VERSION bump.

Formulas
--------
CURRENT : deployed corner_ff (4-term same-sign sum; identically zero by
          p1t*q0t + p1b*q0b == 0).  K = tip + left.
FIX     : (T-nu)*(-2)*(p1t*Ut1*q0t*Vt0 - p1t*Uw1*q0t*Vw0)
          structural analogy; best FE miss% so far; never derived.
DERIVED : (T-nu)*(-2)*( TT - TB - WT + WB )
          = 2*FIX under the parity identity.  Closed-contour Eq.(16) trace.

FE anchors (LESSONS_LEARNED.md Sec 18.83; KFAC-decoded ANSYS)
--------------------------------------------------------------
  l/b=1.5 SYM  Lambda_FE=0.96347
  l/b=2.5 SYM  Lambda_FE=0.96559
  l/b=2.5 ANTI Lambda_FE=0.53126
  l/b=2.5 ANTI Lambda_FE=1.14003

PRE-REGISTERED CRITERIA (written before the run; honor disconfirming results)
-----------------------------------------------------------------------------
PART A (algebra; must pass before any FE reading is trusted):
  A1. CURRENT corner contribution max|K_CURRENT - K_parent| ~ 0 on free_free
      (parent IS CURRENT) AND the raw 4-term sum is ~1e-20 or smaller.
  A2. FIX and DERIVED are genuinely nonzero (max|K - K_CURRENT| > 1e-10).
  A3. K_DERIVED - K_CURRENT  ==  2*(K_FIX - K_CURRENT)  to ~1e-10 relative
      (the identity makes DERIVED exactly 2*FIX as a matrix increment).
  A4. Both subclasses' bc='clamped_free' path is byte-identical to parent
      (max|K - K_parent| == 0 at Table-3 SYM Lambda=0.1551).
  FAIL PART A => stop, do not implement, do not interpret FE.

PART B (production-basis dip location vs FE; n_real=6/n_cpair=0 SYM,
        n_real=5/n_cpair=1 ANTI, shipped im_cap=7.0 -- same convention as
        the Sec 18.83 numbers).  For each anchor, take interior local
        minima of sigma_ratio (NOT the global min of the whole bracket)
        and report the one nearest the FE value.
  B-IMPLEMENT DERIVED if, on a majority of anchors that discriminate
      (a formula finds an interior min the others also find, or one finds
      a min within 8% of FE and another does not):
        DERIVED miss% <= FIX miss%   AND   DERIVED miss% <= CURRENT miss%.
      "At least as well as FIX" is the hand-off bar.
  B-DO NOT IMPLEMENT DERIVED if DERIVED miss% is worse than FIX at a
      majority of discriminating anchors (the predicted "2xFIX overshoots"
      failure).  Leave core_solvers.py unchanged.  This is a discrepancy
      between the continuum trace and the discrete residual, not a license
      to silently halve the coefficient.
  B-INCONCLUSIVE if fewer than 2 anchors discriminate, or if every dip
      is shallower than 1e-2 with no isolated structure.

PART C (secondary; Mxy(corner) residual at n_cpair=3, the one trustworthy
        enlarged size per Sec 18.82).  Control: clamped_free Table-3 SYM
        Lambda=0.1551 must stay in the ~2-5% band at n_cpair=3.
  C is NOT sufficient to promote.  Residuals of 17-49% (FIX's own band)
        are "not well satisfied."  Use only as a tie-break / sanity check:
        if DERIVED residual is *worse* than CURRENT at the FE Lambda on
        both free corners, that counts against implementation even if
        PART B is a narrow win.

Decision for core_solvers.py this session:
  Implement DERIVED only if PART A passes AND PART B says B-IMPLEMENT
  AND PART C does not veto.  Otherwise report the numbers and leave
  the deployed (zero) corner_ff in place.
"""
from __future__ import annotations

import os
import sys
import time

sys.path.insert(0, ".")

from mpmath import mp, mpf, mpc  # noqa: E402

DPS = int(os.environ.get("DPS", "30"))
mp.dps = DPS

from plate_solver import RectOOPAssembler  # noqa: E402
from plate_solver.detectors import (  # noqa: E402
    _RECT_AXIS_TOL,
    _rect_absdet_axis_roots,
    _rect_cnewton,
    rect_resolve_branches,
    rect_select_branches,
)


# ---------------------------------------------------------------------------
# Assembler subclass: one copy of assemble(), corner_mode switch on free_free.
# clamped_free path is byte-for-byte the parent.
# ---------------------------------------------------------------------------
class RectOOPAssemblerCornerMode(RectOOPAssembler):
    """corner_mode: 'CURRENT' | 'FIX' | 'DERIVED'."""

    corner_mode = "CURRENT"

    def assemble(self, full, Lam, sym, lob):
        nu, _R, T, mu = self.nu, self.R, self.T, self.mu
        phin = mp.pi / 2 if sym else mpf(0)
        L = mp.pi * mpf(str(lob)) / 2
        Lf = float(Lam)
        eng = self.eng(sym)
        bd = []
        for xi in full:
            e1, e2 = eng._etas(complex(xi), Lf)
            H1, H2 = eng._amp_ratio(complex(xi), Lf)
            e1 = mpc(e1)
            e2 = mpc(e2)
            H1 = mpc(H1)
            H2 = mpc(H2)
            psi = {m: self._pterms(e1, e2, H1, H2, phin, m) for m in (0, 1, 2, 3)}
            bd.append((mpc(xi), psi))

        def U(xi, r, m, x1):
            return xi ** m * mp.sin(xi * x1 + (r - 1) * mp.pi / 2 + m * mp.pi / 2)

        def peval(terms, x2):
            return sum(c * mp.sin(e * x2 + p) for c, e, p in terms)

        P = len(full)
        K = mp.zeros(2 * P, 2 * P)
        mode = self.corner_mode
        for pi_, (xip, psip) in enumerate(bd):
            for ppi, (xipp, psipp) in enumerate(bd):
                I = {(m, k): self._pint(psip[m], psipp[k])
                     for m in (0, 1, 2, 3) for k in (0, 1, 2, 3)}
                p1t = peval(psip[1], mp.pi / 2)
                p1b = peval(psip[1], -mp.pi / 2)
                q0t = peval(psipp[0], mp.pi / 2)
                q0b = peval(psipp[0], -mp.pi / 2)
                for r in (1, 2):
                    for rp in (1, 2):
                        Ut3 = U(xip, r, 3, L)
                        Ut1 = U(xip, r, 1, L)
                        Ut2 = U(xip, r, 2, L)
                        Ut0 = U(xip, r, 0, L)
                        Vt0 = U(xipp, rp, 0, L)
                        Vt1 = U(xipp, rp, 1, L)
                        tip = (Ut3 * Vt0 * I[(0, 0)] + mu * Ut1 * Vt0 * I[(2, 0)]
                               - Ut2 * Vt1 * I[(0, 0)] - nu * Ut0 * Vt1 * I[(2, 0)])
                        Uw0 = U(xip, r, 0, -L)
                        Uw1 = U(xip, r, 1, -L)
                        Vw0 = U(xipp, rp, 0, -L)
                        Vw1 = U(xipp, rp, 1, -L)
                        if self.bc == "clamped_free":
                            Vw2 = U(xipp, rp, 2, -L)
                            wall = (Uw0 * Vw2 * I[(0, 1)] + T * Uw0 * Vw0 * I[(0, 3)]
                                    - Uw1 * Vw2 * I[(0, 0)] - nu * Uw1 * Vw0 * I[(0, 2)]
                                    - (T - nu) * Uw0 * Vw1 * I[(1, 1)])
                            corner = (T - nu) * (p1t * Uw1 * q0t * Vw0
                                                 + p1b * Uw1 * q0b * Vw0
                                                 - 2 * p1t * Ut1 * q0t * Vt0)
                            K[2 * ppi + (rp - 1), 2 * pi_ + (r - 1)] = tip + wall + corner
                        else:
                            Uw3 = U(xip, r, 3, -L)
                            Uw2 = U(xip, r, 2, -L)
                            left = -(Uw3 * Vw0 * I[(0, 0)] + mu * Uw1 * Vw0 * I[(2, 0)]
                                     - Uw2 * Vw1 * I[(0, 0)] - nu * Uw0 * Vw1 * I[(2, 0)])
                            TT = p1t * Ut1 * q0t * Vt0
                            TB = p1b * Ut1 * q0b * Vt0
                            WT = p1t * Uw1 * q0t * Vw0
                            WB = p1b * Uw1 * q0b * Vw0
                            if mode == "CURRENT":
                                corner_ff = (T - nu) * (-2) * (TT + TB + WT + WB)
                            elif mode == "FIX":
                                corner_ff = (T - nu) * (-2) * (TT - WT)
                            elif mode == "DERIVED":
                                corner_ff = (T - nu) * (-2) * (TT - TB - WT + WB)
                            else:
                                raise ValueError(mode)
                            K[2 * ppi + (rp - 1), 2 * pi_ + (r - 1)] = tip + left + corner_ff
        return K


class AsmCURRENT(RectOOPAssemblerCornerMode):
    corner_mode = "CURRENT"


class AsmFIX(RectOOPAssemblerCornerMode):
    corner_mode = "FIX"


class AsmDERIVED(RectOOPAssemblerCornerMode):
    corner_mode = "DERIVED"


def make_asm(cls, bc="free_free"):
    return cls(nu_b=0.3, R=1.0, T=1.0, dps=DPS, bc=bc)


def max_abs_diff(Ka, Kb):
    n = Ka.rows
    m = mpf(0)
    for i in range(n):
        for j in range(n):
            d = abs(Ka[i, j] - Kb[i, j])
            if d > m:
                m = d
    return m


def frange(lo, hi, step):
    n = int(round((hi - lo) / step))
    return [round(lo + step * i, 6) for i in range(n + 1)]


def sigma_at(asm, sym, lob, Lam, n_real, n_cpair, im_cap=7.0):
    eng = asm.eng(sym)
    reps = rect_resolve_branches(eng, Lam, im_cap=im_cap)
    full = rect_select_branches(reps, n_real=n_real, n_cpair=n_cpair)
    if len(full) < 2:
        return None, 0, full
    K = asm.assemble(full, mpf(str(round(Lam, 6))), sym, lob)
    return asm.equil_sigma(K, full), len(full), full


def local_minima(rows):
    """rows: list of (Lam, sigma_or_None, n_used). Interior local mins only."""
    dips = []
    for i in range(1, len(rows) - 1):
        _, s0, _ = rows[i - 1]
        L1, s1, n1 = rows[i]
        _, s2, _ = rows[i + 1]
        if s0 is None or s1 is None or s2 is None:
            continue
        if s1 < s0 and s1 < s2:
            dips.append((L1, s1, n1))
    return dips


def nearest_dip(dips, lam_fe):
    if not dips:
        return None
    return min(dips, key=lambda t: abs(t[0] - lam_fe))


def resolve_wide(eng, Lam, im_cap=30.0, re_max=1.6, cgrid=0.12):
    reps = []
    reps += [complex(x, 0.0) for x in _rect_absdet_axis_roots(eng, Lam, "re", 1e-3, 9.0, 0.005)]
    reps += [complex(0.0, k) for k in _rect_absdet_axis_roots(eng, Lam, "im", 1e-3, im_cap, 0.005)]
    re = cgrid
    while re < re_max:
        im = cgrid
        while im < im_cap:
            z = _rect_cnewton(eng, complex(re, im), Lam, itmax=30)
            if (z is not None and z.real >= _RECT_AXIS_TOL and z.imag >= _RECT_AXIS_TOL
                    and z.imag <= im_cap and not any(abs(z - c) < 5e-3 for c in reps)):
                reps.append(z)
            im += cgrid
        re += cgrid
    uniq = []
    for z in reps:
        if not any(abs(z - u) < 5e-3 for u in uniq):
            uniq.append(z)
    uniq.sort(key=lambda z: (round(abs(z.imag), 4), abs(z.real)))
    return uniq


def nullvec_complex(K):
    K = K.copy()
    n = K.cols
    dc = [mpf(1)] * n
    for j in range(n):
        cn = mp.sqrt(sum(abs(K[i, j]) ** 2 for i in range(n)))
        if cn > 0:
            for i in range(n):
                K[i, j] /= cn
            dc[j] = cn
    _U, S, V = mp.svd(K, compute_uv=True)
    vals = [S[i] for i in range(S.rows)]
    smin = min(vals)
    smax = max(vals)
    idx = vals.index(smin)
    w = [V[j, idx] for j in range(n)]
    c = [w[j] / dc[j] for j in range(n)]
    nrm = max(abs(x) for x in c)
    c = [x / nrm for x in c]
    return c, float(smin / smax)


def corner_residuals(asm, full, Lam, sym, lob, coeffs):
    phin = mp.pi / 2 if sym else mpf(0)
    L = mp.pi * mpf(str(lob)) / 2
    Lf = float(Lam)
    eng = asm.eng(sym)
    bd = []
    for xi in full:
        e1, e2 = eng._etas(complex(xi), Lf)
        H1, H2 = eng._amp_ratio(complex(xi), Lf)
        e1 = mpc(e1)
        e2 = mpc(e2)
        H1 = mpc(H1)
        H2 = mpc(H2)
        psi = {m: asm._pterms(e1, e2, H1, H2, phin, m) for m in (0, 1, 2, 3)}
        bd.append((mpc(xi), psi))

    def U(xi, r, m, x1):
        return xi ** m * mp.sin(xi * x1 + (r - 1) * mp.pi / 2 + m * mp.pi / 2)

    def peval(terms, x2):
        return sum(c * mp.sin(e * x2 + p) for c, e, p in terms)

    corners = {
        "tip_top": (L, mp.pi / 2),
        "tip_bot": (L, -mp.pi / 2),
        "left_top": (-L, mp.pi / 2),
        "left_bot": (-L, -mp.pi / 2),
    }
    out = {}
    for name, (x1c, x2c) in corners.items():
        mxy = mpc(0)
        wv = mpc(0)
        for pi_, (xip, psip) in enumerate(bd):
            for r in (1, 2):
                col = 2 * pi_ + (r - 1)
                cf = coeffs[col]
                mxy += cf * U(xip, r, 1, x1c) * peval(psip[1], x2c)
                wv += cf * U(xip, r, 0, x1c) * peval(psip[0], x2c)
        relw = float(abs(mxy) / abs(wv)) if abs(wv) > 0 else float("nan")
        out[name] = relw
    return out


def raw_corner_scalars(asm, full, Lam, sym, lob):
    """Return (parity_bracket, CURRENT, FIX, DERIVED) at the (0,0,r=1,rp=1) slot."""
    nu, T = asm.nu, asm.T
    phin = mp.pi / 2 if sym else mpf(0)
    L = mp.pi * mpf(str(lob)) / 2
    Lf = float(Lam)
    eng = asm.eng(sym)
    xi = full[0]
    e1, e2 = eng._etas(complex(xi), Lf)
    H1, H2 = eng._amp_ratio(complex(xi), Lf)
    e1, e2, H1, H2 = mpc(e1), mpc(e2), mpc(H1), mpc(H2)
    psip = {m: asm._pterms(e1, e2, H1, H2, phin, m) for m in (0, 1, 2, 3)}
    psipp = psip

    def U(xi, r, m, x1):
        return xi ** m * mp.sin(xi * x1 + (r - 1) * mp.pi / 2 + m * mp.pi / 2)

    def peval(terms, x2):
        return sum(c * mp.sin(e * x2 + p) for c, e, p in terms)

    xip = mpc(xi)
    r = 1
    rp = 1
    p1t = peval(psip[1], mp.pi / 2)
    p1b = peval(psip[1], -mp.pi / 2)
    q0t = peval(psipp[0], mp.pi / 2)
    q0b = peval(psipp[0], -mp.pi / 2)
    Ut1 = U(xip, r, 1, L)
    Uw1 = U(xip, r, 1, -L)
    Vt0 = U(xip, rp, 0, L)
    Vw0 = U(xip, rp, 0, -L)
    TT = p1t * Ut1 * q0t * Vt0
    TB = p1b * Ut1 * q0b * Vt0
    WT = p1t * Uw1 * q0t * Vw0
    WB = p1b * Uw1 * q0b * Vw0
    ident = p1t * q0t + p1b * q0b
    cur = (T - nu) * (-2) * (TT + TB + WT + WB)
    fix = (T - nu) * (-2) * (TT - WT)
    der = (T - nu) * (-2) * (TT - TB - WT + WB)
    return ident, cur, fix, der


FE_ANCHORS = [
    # (name, lob, sym, Lambda_FE, lo, hi, step, n_real, n_cpair)
    ("lb1.5 SYM", 1.5, True, 0.96347, 0.940, 1.020, 0.002, 6, 0),
    ("lb2.5 SYM", 2.5, True, 0.96559, 0.940, 1.000, 0.005, 6, 0),
    ("lb2.5 ANTI-a", 2.5, False, 0.53126, 0.500, 0.580, 0.005, 5, 1),
    ("lb2.5 ANTI-b", 2.5, False, 1.14003, 1.080, 1.180, 0.005, 5, 1),
]


def part_a():
    print("=" * 72)
    print("PART A -- algebra + clamped_free regression")
    print("=" * 72, flush=True)
    parent = make_asm(RectOOPAssembler, bc="free_free")
    cur = make_asm(AsmCURRENT)
    fix = make_asm(AsmFIX)
    der = make_asm(AsmDERIVED)
    lob, Lam, sym = 1.5, 0.96, True
    eng = parent.eng(sym)
    reps = rect_resolve_branches(eng, Lam)
    full = rect_select_branches(reps, n_real=6, n_cpair=0)
    print(f"  free_free SYM Lam={Lam} l/b={lob} n_used={len(full)} dps={mp.dps}", flush=True)

    ident, c_sc, f_sc, d_sc = raw_corner_scalars(parent, full, Lam, sym, lob)
    print(f"  parity p1t*q0t+p1b*q0b = {ident}")
    print(f"  scalar CURRENT = {c_sc}")
    print(f"  scalar FIX     = {f_sc}")
    print(f"  scalar DERIVED = {d_sc}")
    ratio = abs(d_sc / f_sc) if abs(f_sc) > 0 else None
    print(f"  |DERIVED/FIX| (scalar) = {ratio}", flush=True)

    Kp = parent.assemble(full, mpf(str(Lam)), sym, lob)
    Kc = cur.assemble(full, mpf(str(Lam)), sym, lob)
    Kf = fix.assemble(full, mpf(str(Lam)), sym, lob)
    Kd = der.assemble(full, mpf(str(Lam)), sym, lob)
    d_pc = max_abs_diff(Kp, Kc)
    d_fc = max_abs_diff(Kf, Kc)
    d_dc = max_abs_diff(Kd, Kc)
    # Kd - Kc  vs  2*(Kf - Kc)
    two_fix = Kf.copy()
    n = Kf.rows
    for i in range(n):
        for j in range(n):
            two_fix[i, j] = Kc[i, j] + 2 * (Kf[i, j] - Kc[i, j])
    d_2x = max_abs_diff(Kd, two_fix)
    print(f"  max|K_CURRENT - K_parent| = {d_pc}   (expect ~0)")
    print(f"  max|K_FIX     - K_CURRENT| = {d_fc}   (expect >> 0)")
    print(f"  max|K_DERIVED - K_CURRENT| = {d_dc}   (expect >> 0)")
    print(f"  max|K_DERIVED - (K_CURRENT + 2*(K_FIX-K_CURRENT))| = {d_2x}", flush=True)

    print("  clamped_free Table-3 SYM Lambda=0.1551 regression...", flush=True)
    p_cf = make_asm(RectOOPAssembler, bc="clamped_free")
    f_cf = make_asm(AsmFIX, bc="clamped_free")
    d_cf = make_asm(AsmDERIVED, bc="clamped_free")
    eng_cf = p_cf.eng(True)
    reps_cf = rect_resolve_branches(eng_cf, 0.1551)
    full_cf = rect_select_branches(reps_cf, n_real=6, n_cpair=0)
    Kp_cf = p_cf.assemble(full_cf, mpf("0.1551"), True, 1.5)
    Kf_cf = f_cf.assemble(full_cf, mpf("0.1551"), True, 1.5)
    Kd_cf = d_cf.assemble(full_cf, mpf("0.1551"), True, 1.5)
    d_f_cf = max_abs_diff(Kp_cf, Kf_cf)
    d_d_cf = max_abs_diff(Kp_cf, Kd_cf)
    print(f"  max|FIX.clamped_free    - parent| = {d_f_cf}  (must be 0)")
    print(f"  max|DERIVED.clamped_free - parent| = {d_d_cf}  (must be 0)", flush=True)

    a1 = abs(c_sc) < 1e-16 and abs(ident) < 1e-16
    a2 = abs(f_sc) > 1e-12 and abs(d_sc) > 1e-12 and d_fc > 1e-10 and d_dc > 1e-10
    a3 = (d_2x / d_dc < 1e-8) if d_dc > 0 else False
    a4 = d_f_cf == 0 and d_d_cf == 0
    # mpmath 0 comparison: treat tiny as zero
    a4 = abs(d_f_cf) == 0 and abs(d_d_cf) == 0
    ok = a1 and a2 and a3 and a4
    print(f"  A1 parity+CURRENT~0: {a1}")
    print(f"  A2 FIX/DERIVED nonzero: {a2}")
    print(f"  A3 DERIVED == 2*FIX increment: {a3}")
    print(f"  A4 clamped_free byte-identical: {a4}")
    print(f"  PART A: {'PASS' if ok else 'FAIL'}", flush=True)
    return ok


def part_b(asms):
    print("\n" + "=" * 72)
    print("PART B -- production-basis dip location vs FE anchors")
    print("=" * 72, flush=True)
    summary = []
    for name, lob, sym, lam_fe, lo, hi, step, n_real, n_cpair in FE_ANCHORS:
        print(f"\n--- {name}  Lambda_FE={lam_fe}  grid=[{lo},{hi}] step={step} "
              f"n_real={n_real} n_cpair={n_cpair} ---", flush=True)
        grid = frange(lo, hi, step)
        per = {}
        for label, asm in asms:
            t0 = time.time()
            rows = []
            for Lam in grid:
                s, n, _full = sigma_at(asm, sym, lob, Lam, n_real, n_cpair, im_cap=7.0)
                rows.append((Lam, s, n))
            dt = time.time() - t0
            dips = local_minima(rows)
            # sigma at the grid point nearest FE, and at exact FE (rounded 6)
            s_fe, n_fe, _ = sigma_at(asm, sym, lob, lam_fe, n_real, n_cpair, im_cap=7.0)
            nd = nearest_dip(dips, lam_fe)
            dip_str = ", ".join(f"{L:.4f}({s:.3e})" for L, s, _n in dips) or "NONE"
            print(f"  [{label}] {len(grid)} pts in {dt:.1f}s  n_used~{rows[len(rows)//2][2]}", flush=True)
            print(f"         local minima: {dip_str}")
            if nd is None:
                miss = None
                print(f"         nearest dip: NONE   sigma@FE={s_fe}")
            else:
                miss = 100.0 * abs(nd[0] - lam_fe) / lam_fe
                print(f"         nearest dip: Lam={nd[0]:.4f}  sigma={nd[1]:.4e}  "
                      f"miss={miss:.2f}%   sigma@FE={s_fe}")
            per[label] = {"dips": dips, "nearest": nd, "miss": miss, "sigma_fe": s_fe, "n_fe": n_fe}
        summary.append((name, lam_fe, per))
    return summary


def part_c(asms):
    print("\n" + "=" * 72)
    print("PART C -- Mxy(corner) residual at n_cpair=3 (trustworthy size)")
    print("=" * 72, flush=True)

    print("  CONTROL: clamped_free Table-3 SYM l/b=1.5 Lambda=0.1551 n_cpair=3", flush=True)
    asm_cf = make_asm(RectOOPAssembler, bc="clamped_free")
    eng = asm_cf.eng(True)
    reps = resolve_wide(eng, 0.1551, im_cap=30.0)
    full = rect_select_branches(reps, n_real=6, n_cpair=3)
    print(f"    n_used={len(full)} (raw reps={len(reps)})", flush=True)
    K = asm_cf.assemble(full, mpf("0.1551"), True, 1.5)
    coeffs, ratio = nullvec_complex(K)
    res = corner_residuals(asm_cf, full, 0.1551, True, 1.5, coeffs)
    print(f"    sigma_ratio={ratio:.4e}")
    for k, v in res.items():
        print(f"      {k:9s}  |Mxy|/|w| = {v * 100:7.3f}%")
    tip_pct = res["tip_top"] * 100
    control_ok = 0.5 <= tip_pct <= 12.0
    print(f"    control tip residual in ~2-5% band (accept 0.5-12%)? {control_ok}  ({tip_pct:.2f}%)",
          flush=True)

    residual_rows = []
    for name, lob, sym, lam_fe, _lo, _hi, _st, n_real, _ncp in FE_ANCHORS:
        print(f"\n  {name}  Lambda_FE={lam_fe}  n_real={n_real} n_cpair=3 im_cap=30", flush=True)
        # share the branch set across formulas
        parent = asms[0][1]
        eng = parent.eng(sym)
        reps = resolve_wide(eng, lam_fe, im_cap=30.0)
        full = rect_select_branches(reps, n_real=n_real, n_cpair=3)
        print(f"    n_used={len(full)} (raw reps={len(reps)})", flush=True)
        per = {}
        for label, asm in asms:
            K = asm.assemble(full, mpf(str(round(lam_fe, 6))), sym, lob)
            coeffs, ratio = nullvec_complex(K)
            res = corner_residuals(asm, full, lam_fe, sym, lob, coeffs)
            print(f"    [{label}] sigma_ratio={ratio:.4e}")
            for k, v in res.items():
                print(f"       {k:9s}  |Mxy|/|w| = {v * 100:7.3f}%")
            mean_free = 0.25 * sum(res.values())
            per[label] = {"ratio": ratio, "res": res, "mean_pct": mean_free * 100}
        residual_rows.append((name, per))
    return control_ok, residual_rows


def decide(part_a_ok, summary, control_ok, residual_rows):
    print("\n" + "=" * 72)
    print("DECISION")
    print("=" * 72, flush=True)
    if not part_a_ok:
        print("PART A FAIL => do not implement. core_solvers.py unchanged.")
        return "NO_IMPLEMENT_A"

    print("Per-anchor miss% (nearest interior local min to Lambda_FE):")
    print(f"  {'anchor':<16s} {'FE':>8s} {'CURRENT':>12s} {'FIX':>12s} {'DERIVED':>12s}  note")
    n_disc = 0
    n_derived_le_fix = 0
    n_derived_le_cur = 0
    n_derived_worse_fix = 0
    for name, lam_fe, per in summary:
        def fmt(lab):
            m = per[lab]["miss"]
            return "   no-dip" if m is None else f"{m:8.2f}%"

        notes = []
        mc, mf, md = per["CURRENT"]["miss"], per["FIX"]["miss"], per["DERIVED"]["miss"]
        # discriminate if at least two formulas have a dip, or one has a dip
        # within 8% and another has none
        have = [(lab, per[lab]["miss"]) for lab in ("CURRENT", "FIX", "DERIVED")]
        found = [(lab, m) for lab, m in have if m is not None]
        if len(found) >= 2:
            n_disc += 1
            if md is not None and mf is not None:
                if md <= mf + 1e-9:
                    n_derived_le_fix += 1
                else:
                    n_derived_worse_fix += 1
            if md is not None and mc is not None and md <= mc + 1e-9:
                n_derived_le_cur += 1
            elif md is not None and mc is None:
                n_derived_le_cur += 1
        elif len(found) == 1:
            lab, m = found[0]
            if m <= 8.0:
                n_disc += 1
                notes.append(f"only {lab} found a dip")
                if lab == "DERIVED":
                    n_derived_le_fix += 1
                    n_derived_le_cur += 1
                elif lab == "FIX":
                    n_derived_worse_fix += 1
        notes_s = "; ".join(notes)
        print(f"  {name:<16s} {lam_fe:8.5f} {fmt('CURRENT'):>12s} {fmt('FIX'):>12s} "
              f"{fmt('DERIVED'):>12s}  {notes_s}")

    print(f"\n  discriminating anchors: {n_disc}")
    print(f"  DERIVED miss <= FIX miss:     {n_derived_le_fix}/{n_disc}")
    print(f"  DERIVED miss worse than FIX:  {n_derived_worse_fix}/{n_disc}")
    print(f"  DERIVED miss <= CURRENT miss: {n_derived_le_cur}/{n_disc}")

    if n_disc < 2:
        b_verdict = "B-INCONCLUSIVE"
    elif n_derived_le_fix > n_disc / 2 and n_derived_le_cur >= n_derived_worse_fix:
        # majority DERIVED <= FIX, and not systematically worse than CURRENT
        b_verdict = "B-IMPLEMENT"
    elif n_derived_worse_fix > n_disc / 2:
        b_verdict = "B-DO-NOT-IMPLEMENT"
    else:
        b_verdict = "B-INCONCLUSIVE"
    print(f"  PART B verdict: {b_verdict}")

    print("\nPART C mean |Mxy|/|w| at Lambda_FE (n_cpair=3):")
    c_veto = False
    if not control_ok:
        print("  CONTROL residual out of band -- PART C not trusted, no veto from C.")
    else:
        for name, per in residual_rows:
            line = "  " + f"{name:<16s}"
            for lab in ("CURRENT", "FIX", "DERIVED"):
                line += f"  {lab}={per[lab]['mean_pct']:6.1f}%"
            print(line)
            # veto: DERIVED worse than CURRENT on mean residual AND worse than FIX
            if (per["DERIVED"]["mean_pct"] > per["CURRENT"]["mean_pct"]
                    and per["DERIVED"]["mean_pct"] > per["FIX"]["mean_pct"]):
                c_veto = True
                print(f"    veto flag at {name}: DERIVED residual worse than both")
        if c_veto:
            print("  PART C veto: DERIVED residual worse than CURRENT and FIX on at least one anchor.")
        else:
            print("  PART C no veto.")

    if b_verdict == "B-IMPLEMENT" and not c_veto:
        decision = "IMPLEMENT_DERIVED"
    else:
        decision = "NO_IMPLEMENT"
    print(f"\n  FINAL: {decision}")
    print("  core_solvers.py will be left unchanged this run unless FINAL is "
          "IMPLEMENT_DERIVED (applied in a separate explicit edit after this print).")
    return decision


def main():
    print(f"dps={mp.dps}  cwd={os.getcwd()}")
    print("probe_rect_ff_oop_derived_vs_fix_2026-09-03")
    t0 = time.time()
    ok = part_a()
    if not ok:
        print("\nStopping after PART A FAIL.")
        print(f"wall time: {time.time() - t0:.1f}s", flush=True)
        return 1
    asms = [
        ("CURRENT", make_asm(AsmCURRENT)),
        ("FIX", make_asm(AsmFIX)),
        ("DERIVED", make_asm(AsmDERIVED)),
    ]
    summary = part_b(asms)
    control_ok, residual_rows = part_c(asms)
    decide(ok, summary, control_ok, residual_rows)
    print(f"\nwall time: {time.time() - t0:.1f}s", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
