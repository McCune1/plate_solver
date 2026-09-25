# -*- coding: utf-8 -*-
"""
probe_piezo_p5_e15_segment_charge_2026-09-24.py

Paper 5 FUTURE_WORK "Decision needed -- sensing admittance vs e15*gamma_rz".
Decision (2026-09-24): option (a), a segmented-electrode harmonic FE run.
This probe is its PRE-REGISTRATION: it computes, before any FE exists,
what each model predicts for the per-face charge on the inner segment
[r_i, r_*] of the grounded F-F SC ring under a harmonic ring load F at r_F.

Derivation: Project Knowledge/PAPER5_E15_SEGMENT_CHARGE_DERIVATION_2026-09-24.md.
Short form. Gauss on the column r_i<r<r_*, |z|<=H with an insulated inner
rim gives the physical electrode charge on ONE face as

    Q_face = pi * r_* * int_{-H}^{H} D_r dz   (top = bottom by symmetry).

Kirchhoff (what Paper 5 Sec 5.2 uses): D_r = Xi11 E_r, E_r = -phibar' f(z),
    int D_r dz = -I0 Xi11 phibar'          (I0 = 4H/3)
    -> Q_K = -pi r_* I0 Xi11 phibar' = -(Y_sense code value)*F.
Consistent leading order: the transverse shear stress T5 is fixed by
equilibrium (int T5 dz = Q_r), so the shear strain is free, and
    D_r = (e15/c44) T5 + (Xi11 + e15^2/c44) E_r
    int D_r dz = (e15/c44) Q_r - I0 (Xi11 + e15^2/c44) phibar'.
Both terms are O(H^2 k^3 w): same order as the Kirchhoff term.
With phibar ~= c Lap w and Q_r ~= -D_eff (Lap w)' the ratio is a constant,
    kappa = 1 + e15^2/(c44 Xi11) + e15 (c11bar Xi33bar + e31bar^2)
                                     / (c44 Xi11 e31bar),
independent of omega, H, r_F, r_* at leading order.

Gates (sandbox, pre-FE):
  G_A  code q-bracket == -Q_r from the equilibrium integral
       Q_r(r_*) = -(A2 w^2 / r_*) int_{r_i}^{r_*} w r dr  (rel 1e-8)
  G_B  Q_K == -Y_sense code value (rel 1e-12)
  G_C  kappa from the driven solve within 1% of closed-form kappa at
       every test omega, for both c44 values
Outputs the FE prediction table to piezo_p5_e15_segment_charge_predictions.json.

c44: Duan 2005 Table 1 prints C55 = 73 GPa (almost certainly its C13
repeated, LESSONS Sec 18.236); physical PZT-4 is ~26 GPa (Liu 2002). The
Kirchhoff charge does not depend on c44 at all; the corrected charge does.
FE is run at both, which makes the test two-sided.

LF line endings (standalone probe). No package edit. No SOLVER_VERSION use.
"""
import json
import os
import sys
import time

sys.path.insert(0, os.environ.get("PKG_PATH",
                                  os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")))
from plate_solver.piezo_monolithic import PiezoMonolithicOutOfPlaneSolver  # noqa: E402
from mpmath import mp, mpf, quad  # noqa: E402

DPS = int(os.environ.get("PIEZO_DPS", "40"))
mp.dps = DPS

# Same material/geometry as every Paper 5 probe (Duan Table 1 PZT-4).
R_I, R_O, H = 0.1, 0.6, 0.01
C11E, C12E, C13E, C33E = 132e9, 71e9, 73e9, 115e9
RHO = 7500.0
E31, E33, E15 = 4.1, 14.1, 10.5
X11, X33 = 7.124e-9, 5.841e-9
R_F = mpf("0.35")
R_STAR = mpf("0.30")
F = mpf(1)
C44_CASES = {"c44_73": 73e9, "c44_26": 26e9}
# Test omegas (rad/s): Y_sense has poles at 471.36 and ~2220 and a zero
# near 360 rad/s; these sit away from all three.
OMEGAS = [100.0, 200.0, 800.0, 1000.0, 1200.0]


def solver():
    return PiezoMonolithicOutOfPlaneSolver(
        R_I, R_O, H, C11E, C12E, C13E, C33E, RHO,
        e31=E31, e33=E33, X11=X11, X33=X33, dps=DPS)


def q_bracket(s, res, r):
    c = res["c_I"] if r <= res["r_F"] else res["c_II"]
    out = mpf(0)
    idx = 0
    for chi, lam in res["lams"]:
        Z1, Z2, dZ1, dZ2 = s._radial_quad(0, r, lam)
        Ki = res["d"] * lam + res["K_pref"] * chi
        for Z, dZ in ((Z1, dZ1), (Z2, dZ2)):
            out += c[idx] * Ki * dZ
            idx += 1
    return out


def kappa_closed(s, c44):
    c11b = s._c11_bar()
    e31b = s._e31_bar()
    X33b = s._Xi33_bar()
    X11b = s._Xi11_bar()
    e15 = mpf(E15)
    c44 = mpf(c44)
    return (1 + e15 ** 2 / (c44 * X11b)
            + e15 * (c11b * X33b + e31b ** 2) / (c44 * X11b * e31b))


def main():
    t0 = time.time()
    s = solver()
    I0 = 4 * mpf(H) / 3
    A2 = s._A2()
    X11m = mpf(X11)
    e15 = mpf(E15)
    rows = []
    gA = gB = gC = True
    kap_cf = {k: kappa_closed(s, v) for k, v in C44_CASES.items()}
    for om in OMEGAS:
        om_m = mpf(om)
        res = s.driven_ff_force_sc(om_m, R_F, F=F, n=0)
        _, phip = s._driven_eval(res, R_STAR)
        phip = mp.re(phip)
        wint = quad(lambda r: mp.re(s._driven_eval(res, r)[0]) * r,
                    [s.r_i, R_STAR])
        Qr = -(A2 * om_m ** 2 / R_STAR) * wint
        qb = mp.re(q_bracket(s, res, R_STAR))
        relA = abs(qb + Qr) / abs(Qr)
        gA = gA and relA < mpf("1e-8")
        y_code = mp.re(s.Y_sense(om_m, R_F, R_STAR, F=F, n=0))
        QK = -mp.pi * R_STAR * I0 * X11m * phip
        relB = abs(QK + y_code * F) / abs(QK)
        gB = gB and relB < mpf("1e-12")
        row = dict(omega=om, f_Hz=om / (2 * 3.141592653589793),
                   Y_code=float(y_code), Q_K_face=float(QK),
                   Q_r=float(Qr), relA=float(relA), relB=float(relB))
        for key, c44 in C44_CASES.items():
            c44m = mpf(c44)
            shear = mp.pi * R_STAR * (e15 / c44m) * Qr
            perm = -mp.pi * R_STAR * I0 * (e15 ** 2 / c44m) * phip
            QC = QK + shear + perm
            kap = QC / QK
            relC = abs(kap - kap_cf[key]) / abs(kap_cf[key])
            gC = gC and relC < mpf("0.01")
            row[key] = dict(Q_C_face=float(QC), kappa=float(kap),
                            shear_over_QK=float(shear / QK),
                            perm_over_QK=float(perm / QK),
                            kappa_closed=float(kap_cf[key]),
                            relC=float(relC))
        rows.append(row)
        print("omega=%7.1f  Q_K=%+.6e  kappa73=%+.5f  kappa26=%+.5f  relA=%.1e"
              % (om, float(QK), row["c44_73"]["kappa"],
                 row["c44_26"]["kappa"], float(relA)), flush=True)
    verdict = "PASS_ALL" if (gA and gB and gC) else "FAIL"
    out = dict(
        probe="probe_piezo_p5_e15_segment_charge_2026-09-24",
        verdict=verdict, G_A=gA, G_B=gB, G_C=gC, dps=DPS,
        geometry=dict(r_i=R_I, r_o=R_O, H=H, r_F=float(R_F),
                      r_star=float(R_STAR), F=float(F)),
        kappa_closed={k: float(v) for k, v in kap_cf.items()},
        sign_convention=("Q_*_face = physical free charge on ONE face "
                         "electrode segment [r_i,r_*] (top = bottom). "
                         "Paper/code Y_sense = -Q_K_face/F."),
        rows=rows, seconds=time.time() - t0)
    here = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(here, "piezo_p5_e15_segment_charge_predictions.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    print("kappa_closed:", {k: round(float(v), 6) for k, v in kap_cf.items()})
    print("SENTINEL %s  G_A=%s G_B=%s G_C=%s  (%.1f s)"
          % (verdict, gA, gB, gC, time.time() - t0))


if __name__ == "__main__":
    main()
