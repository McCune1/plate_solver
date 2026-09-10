#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Shared Paper 3 overnight helpers. Not imported by plate_solver (no
SOLVER_VERSION implication). Probes at the package root import this.

Literature numbers are transcribed from the source PDFs, not reconstructed:
  SP-160 Table 2.35  -- Leissa_SP160.pdf printed p. 32 (Vogel & Skinner)
  SP-160 Table 2.16  -- printed p. 20, nu=1/3, b/a=0.5, n=0 exact
  SP-160 Table 2.18  -- printed p. 22, outer C inner C (Vogel)
  SP-160 Table 2.20  -- printed p. 22, outer C inner S
  SP-160 Table 2.22  -- printed p. 23, outer C inner F
  SP-160 Table 2.24  -- printed p. 25, outer S inner C
  SP-160 Table 2.26  -- printed p. 26, outer S inner S
  SP-160 Table 2.28  -- printed p. 27, outer S inner F
  SP-160 Table 2.30  -- printed p. 29, outer F inner C
  SP-160 Table 2.33  -- printed p. 30, outer F inner S
  SP-160 Table 2.31  -- printed p. 29, Southwell nu=0.3
  SP-160 Table 2.5   -- printed p. 11, nu=0.33
  Narita 1984 Table 2 -- Narita1984.pdf p. 3
  Narita 1984 Table 3 -- Narita1984.pdf p. 5 (polar orthotropic F-F)
  Irie 1984 Table 2   -- Irie1984.pdf p. 3 (F-F in-plane)

Screen freeze (2026-09-08, after job 2461443). Named, not retuned to
include a published match:
  * Ladder still evaluates delta = 1e-2 .. 1e-6 (unchanged).
  * CLEAN  = inner 5/5.
  * FLOOR  = the only inner miss is delta=1e-6, outer 5/5, shrinking.
    This is the dps=40 floor Sec. 18.38 already named. It does not fail
    a frequency row.
  * NOT    = any inner miss at delta >= 1e-5, or outer/shrinking broken.
  * POS controls must be CLEAN. NEG must be NOT (or at least not CLEAN).
  * A literature MATCH is relative miss vs the printed number. FLOOR
    does not convert a MATCH into a MISS.
"""
from __future__ import annotations

import json
import math
import os

from mpmath import mp, mpf, mpc

from plate_solver.ring_disk import inplane_lambda_irie


# ---------------------------------------------------------------------------
# SP-160 Table 2.35. omega a^2 sqrt(rho/D), free-free annulus.
# Columns b/a = 0.1, 0.3, 0.5, 0.7, 0.9. nu of the printed table is 1/3
# (companion Table 2.34 / Fig. 2.22). (4,0) is not printed.
# ---------------------------------------------------------------------------
VOGEL_BA = (0.1, 0.3, 0.5, 0.7, 0.9)
# (n, s): tuple of 5 values, one per b/a
VOGEL_235 = {
    (2, 0): (5.30, 4.91, 4.28, 3.57, 2.94),
    (3, 0): (12.4, 12.26, 11.4, 9.86, 8.14),
    (0, 1): (8.77, 8.36, 9.32, 13.2, 34.9),
    (1, 1): (20.5, 18.3, 17.2, 22.0, 55.7),
    (2, 1): (34.9, 33.0, 31.1, 37.8, 93.8),
    (3, 1): (53.0, 51.0, 47.4, 55.7, 135.0),
    (0, 2): (38.2, 50.4, 92.3, 251.0, 2238.0),
    (1, 2): (59.0, 58.8, 96.3, 253.0, 2240.0),
}

# Sec. 18.38 five roots at b/a=0.5, nu=0.30 (native Omega). Do not edit.
SEC_18_38_FIVE = (
    (2, 0, 0.1081885366),
    (0, 1, 0.2359132157),
    (3, 0, 0.2894098129),
    (1, 1, 0.4356362253),
    (4, 0, 0.5336384799),
)

NARITA_2A = (  # b/a=0.3, nu=0.3, OOP F-F
    (2, 0, 4.9060),
    (0, 1, 8.3535),
    (3, 0, 12.266),
    (1, 1, 18.292),
)

# SP-160 Table 2.5, completely free circular plate, nu=0.33.
# s=0 starts at n=2 (n=0,1 rigid). s=1 starts at n=0.
LEISSA_25 = {
    (2, 0): 5.253,
    (3, 0): 12.23,
    (4, 0): 21.6,
    (0, 1): 9.084,
    (1, 1): 20.52,
    (2, 1): 35.25,
}

NARITA_2B = (  # solid disk via b/a=0.001, nu=0.33
    (2, 0, 5.2620),
    (0, 1, 9.0689),
    (3, 0, 12.244),
    (1, 1, 20.513),
)

# ---------------------------------------------------------------------------
# SP-160 mixed-edge OOP. Leissa titles list OUTER first. Our API is
# (bc_inner, bc_outer). Transcribed from Leissa_SP160.pdf printed pages
# 23 and 29 (PDF 31 and 37), not reconstructed.
# ---------------------------------------------------------------------------
# Table 2.22 (Vogel ref. 2.46): "Clamped, Free" = outer C, inner F.
# Same b/a columns as Table 2.35. nu of the printed Vogel set is the
# Phase A question (run 1/3 and 0.30). Companion Table 2.21 is Raju
# nu=1/3, s=0 only; Table 2.16 exact 17.51 at b/a=0.5 n=0 is the
# axisymmetric F-inner C-outer check.
VOGEL_222_FC = {
    (0, 0): (10.2, 11.4, 17.7, 43.1, 360),
    (1, 0): (21.1, 19.5, 22.0, 45.3, 362),
    (2, 0): (34.5, 32.5, 32.0, 51.5, 365),
    (3, 0): (51.0, 49.1, 45.8, 61.3, 370),
    (0, 1): (39.5, 51.7, 93.8, 253, 2219),
    (1, 1): (60.0, 59.8, 97.3, 254, 2220),
    (2, 1): (83.4, 79.0, 108.0, 259, 2225),
    (0, 2): (90.4, 132.0, 253.0, 692, 6183),
}
# Table 2.16 (Joga-Rao / Pickett), nu=1/3, b/a=0.5, n=0 axisym exact.
# Leissa r=a is OUTER, r=b is INNER. Keys are (bc_inner, bc_outer).
# F-F is omitted in the source ("all but the free-free cases").
TABLE_216_FC_N0 = 17.51   # outer C, inner F
TABLE_216_CF_N0 = 13.05   # outer F, inner C
TABLE_216_EXACT = {
    ("C", "C"): 89.30,
    ("S", "C"): 64.06,
    ("F", "C"): 17.51,
    ("C", "S"): 59.91,
    ("S", "S"): 40.01,
    ("F", "S"): 5.040,
    ("C", "F"): 13.05,
    ("S", "F"): 4.060,
}

# Table 2.30 (Vogel ref. 2.46): "Free, Clamped" = outer F, inner C.
# n=1 s=0 is below n=0 s=0 at small b/a (Raju Table 2.29).
# The n=0 column at b/a=0.9 is a source defect (Sec. 18.123): n=0 and
# n=1 are near-degenerate there, and Table 2.22 at the same b/a prints
# the expected near-equals (360/362, 2219/2220). Table 2.30 prints
# (0,0)=51.5 vs (1,0)=345 and (0,1)=970 vs (1,1)=2189. Transcribe
# the printed numbers; do not retune freq_gate to absorb them.
VOGEL_230_CF = {
    (1, 0): (3.14, 6.33, 13.3, 37.5, 345),
    (0, 0): (4.23, 6.66, 13.0, 37.0, 51.5),
    (2, 0): (5.62, 7.96, 14.7, 39.3, 347),
    (3, 0): (12.4, 13.27, 18.5, 42.6, 352),
    (0, 1): (25.3, 42.6, 85.1, 239, 970),
    (1, 1): (27.3, 44.6, 86.7, 241, 2189),
    (2, 1): (37.0, 50.9, 91.7, 246, 2194),
    (3, 1): (53.2, 62.1, 100, 253, 2200),
}

# Paper scoring of Table 2.30 (Sec. 18.123). Job 2463845 printed 72/78
# by skipping only (0,0)=51.5 at both nu. The paper reports 72/76 by
# treating the whole n=0 b/a=0.9 column as one named source skip.
# freq_gate is unchanged. 970.0 == 970.
VOGEL_230_NAMED_SKIPS = (
    (0, 0, 0.9, 51.5),
    (0, 1, 0.9, 970.0),
)


def is_vogel_230_named_skip(n, s, ba, printed):
    """True for Table 2.30 n=0 at b/a=0.9 -- source column, not a miss.

    Does not change freq_gate. Paper denominator 76; mill job 2463845
    printed denominator 78.
    """
    return (int(n), int(s), float(ba), float(printed)) in VOGEL_230_NAMED_SKIPS


# ---------------------------------------------------------------------------
# Phase F remaining Vogel combinations (ref. 2.46). Leissa titles list
# OUTER first. API is (bc_inner, bc_outer). None = a blank printed cell;
# skip, do not invent. Transcribed from Leissa_SP160.pdf pages 22, 25,
# 26, 27, 30 (PDF 30, 33, 34, 35, 38).
# ---------------------------------------------------------------------------
# Table 2.18 "Clamped, Clamped" = inner C / outer C.
# (2,0) and (2,1) at b/a=0.9 are blank in the source.
VOGEL_218_CC = {
    (0, 0): (27.3, 45.2, 89.2, 248, 2237),
    (1, 0): (28.4, 46.6, 90.2, 249, 2238),
    (2, 0): (36.7, 51.0, 93.3, 251, None),
    (3, 0): (51.2, 60.0, 99.0, 256, 2243),
    (0, 1): (75.3, 125.0, 246.0, 686, 6167),
    (1, 1): (78.6, 127.0, 248.0, 686, 6167),
    (2, 1): (90.5, 134.0, 253.0, 689, None),
    (3, 1): (112.0, 145.0, 259.0, 694, 6174),
}
# Table 2.20 "Clamped, Simply Supported" = inner S / outer C.
VOGEL_220_SC = {
    (0, 0): (22.6, 33.7, 63.9, 175, 1550),
    (1, 0): (25.1, 35.8, 65.4, 175, 1551),
    (2, 0): (35.4, 42.8, 70.0, 178, 1553),
    (3, 0): (51.0, 54.7, 78.1, 185, 1558),
    (0, 1): (65.6, 104.0, 202.0, 558, 5004),
    (1, 1): (70.5, 107.0, 203.0, 560, 5004),
    (2, 1): (86.7, 116.0, 210.0, 563, 5007),
    (3, 1): (111.0, 130.0, 218.0, 570, 5012),
}
# Table 2.24 "Simply Supported, Clamped" = inner C / outer S.
VOGEL_224_CS = {
    (0, 0): (17.8, 29.9, 59.8, 168, 1535),
    (1, 0): (19.0, 31.4, 61.0, 170, 1536),
    (2, 0): (26.8, 36.2, 64.6, 172, 1538),
    (3, 0): (40.0, 45.4, 71.0, 177, 1541),
    (0, 1): (60.1, 100.0, 198.0, 552, 4989),
    (1, 1): (62.8, 102.0, 200.0, 553, 4989),
    (2, 1): (74.7, 109.0, 205.0, 557, 4992),
    (3, 1): (95.3, 120.0, 211.0, 563, 4997),
}
# Table 2.26 "Simply Supported, Simply Supported" = inner S / outer S.
# (2,1) at b/a=0.3 prints 933; neighbouring s=1 cells are 81.8 / 84.6 /
# 108, so 93.3 is the physically expected value. Transcribe 933; name
# the cell if the solver recovers ~93.3. Do not retune freq_gate.
VOGEL_226_SS = {
    (0, 0): (14.5, 21.1, 40.0, 110, 988),
    (1, 0): (16.7, 23.3, 41.8, 112, 988),
    (2, 0): (25.9, 30.2, 47.1, 116, 993),
    (3, 0): (40.0, 42.0, 56.0, 122, 998),
    (0, 1): (51.7, 81.8, 159.0, 439, 3948),
    (1, 1): (56.5, 84.6, 161.0, 441, 3948),
    (2, 1): (71.7, 933.0, 167.0, 444, 3952),
    (3, 1): (94.7, 108.0, 177.0, 453, 3958),
}
VOGEL_226_NAMED_SKIPS = (
    (2, 1, 0.3, 933.0),
)
# Table 2.28 "Simply Supported, Free" = inner F / outer S.
VOGEL_228_FS = {
    (0, 0): (4.86, 4.66, 5.07, 6.93, 17.7),
    (1, 0): (13.9, 12.8, 11.6, 13.3, 29.7),
    (2, 0): (25.4, 24.1, 22.3, 24.3, 51.2),
    (3, 0): (40.0, 38.8, 35.7, 37.2, 74.5),
    (0, 1): (29.4, 37.0, 65.8, 175, 1550),
    (1, 1): (48.0, 45.8, 69.9, 178, 1553),
    (2, 1): (69.2, 65.1, 81.1, 185, 1558),
    (0, 2): (74.8, 107.0, 203.0, 558, 5004),
}
# Table 2.33 "Free, Simply Supported" = inner S / outer F.
# n=1 s=0 is below n=0 s=0 at small b/a (same pattern as Table 2.30).
VOGEL_233_SF = {
    (1, 0): (2.30, 3.32, 4.86, 8.34, 25.9),
    (0, 0): (3.45, 3.42, 4.11, 6.18, 17.2),
    (2, 0): (5.42, 6.08, 7.98, 13.4, 42.6),
    (3, 0): (12.4, 12.6, 14.0, 20.5, 61.4),
    (0, 1): (20.8, 31.6, 61.0, 170, 1535),
    (1, 1): (24.1, 34.5, 63.3, 172, 1536),
    (2, 1): (35.8, 43.0, 69.7, 177, 1541),
    (3, 1): (53.0, 56.7, 80.3, 185, 1548),
}

# Production list for the nine-combination claim. F-F / F-C / C-F are
# Phases A and D; the rest is Phase F. inner/outer as the local API.
VOGEL_NINE = (
    ("2.35", "F", "F", None),          # already Phase A; not re-swept here
    ("2.22", "F", "C", None),          # already Phase D
    ("2.30", "C", "F", None),          # already Phase D
    ("2.18", "C", "C", VOGEL_218_CC),
    ("2.20", "S", "C", VOGEL_220_SC),
    ("2.24", "C", "S", VOGEL_224_CS),
    ("2.26", "S", "S", VOGEL_226_SS),
    ("2.28", "F", "S", VOGEL_228_FS),
    ("2.33", "S", "F", VOGEL_233_SF),
)
PHASE_F_VOGEL = tuple(t for t in VOGEL_NINE if t[3] is not None)


def is_vogel_226_named_skip(n, s, ba, printed):
    """True for Table 2.26 (2,1) at b/a=0.3 printed 933 -- source typeset."""
    return (int(n), int(s), float(ba), float(printed)) in VOGEL_226_NAMED_SKIPS


# Narita 1984 Table 3. Omega = omega a^2 sqrt(rho h / D_r), F-F.
# Parameters [D_theta/D_r, H/D_r, nu_theta]. Roadmap: one material, two
# b/a. High-modulus graphite epoxy is Table 1's material.
NARITA_3_HMGE = {
    "name": "high modulus graphite epoxy",
    "params": (0.04, 0.055, 0.012),
    "ba": {
        0.5: {
            (2, 0): 0.935,
            (0, 1): 1.918,
            (3, 0): 2.503,
            (4, 0): 4.636,
            (1, 1): 3.961,
        },
        0.3: {
            (2, 0): 1.132,
            (0, 1): 1.693,
            (3, 0): 2.897,
            (4, 0): 5.180,
            (1, 1): 4.519,
        },
    },
}


def make_narita_orthotropic(Dth_Dr, H_Dr, nu_th, E_r=210e9, rho=7800.0):
    """Build OrthotropicMaterial from Narita's (D_theta/D_r, H/D_r, nu_theta).

    R = E_theta/E_r = D_theta/D_r. Reciprocity gives nu_r = nu_theta / R.
    Narita H = D_r nu_theta + 2 D_G, so H/D_r is the ring T:
      T = nu_theta + 2 G (1-nu_r nu_theta)/E_r
    and G follows. flexural_lambda2 then uses D = D_r, matching Narita Omega
    on the isotropic reduction (Table 2).
    """
    from plate_solver.geometry import OrthotropicMaterial
    R = float(Dth_Dr)
    nu_th = float(nu_th)
    H_Dr = float(H_Dr)
    nu_r = nu_th / R
    G = float(E_r) * (H_Dr - nu_th) / (2.0 * (1.0 - nu_r * nu_th))
    return OrthotropicMaterial(
        E_r=E_r, E_theta=R * float(E_r), nu_r=nu_r, G_rtheta=G, rho=rho)


# Table 2.31 Southwell (ref. 2.37), nu=0.3, same BC as 2.30 (inner C,
# outer F). Not on the Vogel b/a grid. Lowest root per n (s=0).
# Hold-out for Phase D C-F.
SOUTHWELL_231 = (
    (0, 0, 0.276, 6.25),
    (0, 0, 0.642, 25.0),
    (0, 0, 0.840, 81.0),
    (1, 0, 0.060, 2.82),
    (1, 0, 0.397, 9.00),
    (1, 0, 0.603, 21.2),
    (1, 0, 0.634, 25.0),
    (1, 0, 0.771, 64.0),
    (1, 0, 0.827, 121.0),
    (2, 0, 0.186, 6.25),
    (2, 0, 0.349, 9.00),
    (2, 0, 0.522, 16.0),
    (2, 0, 0.769, 64.0),
    (2, 0, 0.81, 100),
    (3, 0, 0.43, 16.0),
    (3, 0, 0.59, 25.0),
    (3, 0, 0.71, 49.0),
    (3, 0, 0.82, 100),
)

# Irie 1984 Table 2, F(i)-F(o), nu=0.3, lambda = omega a sqrt(rho(1-nu^2)/E).
# Per beta: n -> (first, second). n=0 radial and torsional stored apart.
# Circular column (beta=0, actually 0.01) is NOT searched as Ri=0 tonight.
IRIE_T2_FF = {
    0.2: {
        "n0_rad": (1.802, 4.748),
        "n0_tor": (3.090, 5.209),
        1: (1.652, 3.842),
        2: (1.110, 2.403),
        3: (2.071, 3.401),
        4: (2.767, 4.389),
    },
    0.4: {
        "n0_rad": (1.452, 5.529),
        "n0_tor": (3.530, 6.445),
        1: (1.683, 4.044),
        2: (0.721, 2.451),
        3: (1.618, 3.346),
        4: (2.482, 4.227),
    },
    0.6: {
        "n0_rad": (1.220, 7.979),
        "n0_tor": (4.867, 9.409),
        1: (1.618, 5.092),
        2: (0.418, 2.474),
        3: (1.043, 3.433),
        4: (1.752, 4.386),
    },
    0.8: {
        "n0_rad": (1.065, 15.754),
        "n0_tor": (9.380, 18.630),
        1: (1.488, 9.458),
        2: (0.178, 2.339),
        3: (0.487, 3.299),
        4: (0.892, 4.290),
    },
}


def ladder_class(ladder):
    """CLEAN / FLOOR / NOT. Does not change SIGNFLIP_DELTAS."""
    if ladder is None:
        return "NOT"
    if ladder.clean:
        return "CLEAN"
    misses = [r.delta for r in ladder.rungs if not r.flip_in]
    if (misses == [1.0e-6]
            and ladder.n_out == len(ladder.rungs)
            and ladder.shrinking):
        return "FLOOR"
    return "NOT"


def rel_miss(solver_val, printed):
    printed = float(printed)
    if printed == 0.0:
        return float("inf")
    return abs(float(solver_val) - printed) / abs(printed)


def freq_gate(rel, printed):
    """MATCH / WEAK / MISS.

    MATCH if rel < 0.01 (covers Narita 5e-4 and 3-figure Vogel 1%).
    WEAK if rel < 0.03. Else MISS. Narita hold-out in Phase A also
    checks rel < 5e-4 explicitly.
    """
    rel = float(rel)
    if rel < 0.01:
        return "MATCH"
    if rel < 0.03:
        return "WEAK"
    return "MISS"


def count_radial_nodes(solver, n, Om, n_r=81, bc_inner="F", bc_outer="F"):
    """Count interior sign changes of Re W(r) along a radius (s-ID).

    Mixed-edge uses the Phase D 4x4, not the F-F `_amp_mp`.
    """
    try:
        from plate_solver.ring_disk import ring_L_mp
        xi = mpc(n)
        Om_m = mpf(Om)
        sols = solver._series_mp(xi, Om_m)
        inner = str(bc_inner).upper()[:1]
        outer = str(bc_outer).upper()[:1]
        if (inner, outer) == ("F", "F"):
            A = solver._amp_mp(xi, Om_m, sols)
        else:
            L = ring_L_mp(solver, n, Om, inner, outer)
            cn = [sum(abs(L[r, q]) ** 2 for r in range(4)) for q in range(4)]
            pc = max(range(4), key=lambda q: cn[q])
            oc = [q for q in range(4) if q != pc]
            rn = [sum(abs(L[r, q]) ** 2 for q in range(4)) for r in range(4)]
            ur = sorted(range(4), key=lambda r: rn[r], reverse=True)[:3]
            sub = mp.matrix(3, 3)
            rhs = mp.matrix(3, 1)
            for i, row in enumerate(ur):
                rhs[i, 0] = -L[row, pc]
                for j, col in enumerate(oc):
                    sub[i, j] = L[row, col]
            try:
                t = mp.lu_solve(sub, rhs)
                A = [mpc(0)] * 4
                A[pc] = mpc(1)
                for j, col in enumerate(oc):
                    A[col] = t[j, 0]
            except Exception:
                A = [mpc(1), mpc(0), mpc(0), mpc(0)]
            W0 = sum(A[q] * solver._ev_mp(sols[q], mpf(0), 0) for q in range(4))
            if abs(W0) > mpf("1e-30"):
                A = [A[q] / W0 for q in range(4)]
        xs = [(-mp.pi / 2) + (mp.pi) * i / (n_r - 1) for i in range(n_r)]
        vals = []
        for x in xs:
            z = complex(solver._W_mp(x, sols, A, 0))
            vals.append(z.real)
        # drop two points at each end (the arcs)
        interior = vals[2:-2]
        signs = [0 if abs(v) < 1e-18 else (1 if v > 0 else -1)
                 for v in interior]
        nchg = 0
        last = 0
        for s in signs:
            if s == 0:
                continue
            if last != 0 and s != last:
                nchg += 1
            last = s
        return nchg
    except Exception as exc:
        return "ERR:%s" % type(exc).__name__


def result_row(res, printed, motion="oop"):
    lam = res.lambda2 if motion == "oop" else None
    rel = rel_miss(lam, printed) if lam is not None else float("nan")
    klass = ladder_class(res.ladder)
    gate = freq_gate(rel, printed) if lam is not None else "NA"
    return {
        "n": res.n,
        "Omega": res.Omega,
        "lambda2": lam,
        "printed": printed,
        "rel": rel,
        "log_abs_det": res.log_abs_det,
        "ladder": klass,
        "n_in": res.ladder.n_in,
        "n_out": res.ladder.n_out,
        "shrinking": res.ladder.shrinking,
        "freq_gate": gate,
        "n_evals": res.n_evals,
    }


# FE n-ID from production PLNSOL,U,Z PNGs (jobs 2461506/2461507).
# PNG 000 = FE mode 7. n = nodal diameters; s = interior radial circles.
# Independent of which solver window found the root. n is None for
# in-plane (UZ ~ 1e-13). Full tables: Project Knowledge/PAPER3_NID_FE_2026-09-08.md
# Each row: (mode_lo, mode_hi, Omega_lit, n, s).
FE_NID_PNG_LO = 7
FE_NID_PNG_HI = 22
FE_NID = {
    "ann_ba01_nu030": (
        (7, 8, 5.297710, 2, 0),
        (9, 9, 8.773326, 0, 1),
        (10, 11, 12.421251, 3, 0),
        (12, 13, 20.386555, 1, 1),
        (14, 15, 21.798514, 4, 0),
        (16, 17, 33.424079, 5, 0),
        (18, 19, 34.878263, 2, 1),
        (20, 20, 38.210428, 0, 2),
        (21, 22, 47.256843, 6, 0),
    ),
    "ann_ba03_nu030": (
        (7, 8, 4.898643, 2, 0),
        (9, 9, 8.352258, 0, 1),
        (10, 11, 12.248222, 3, 0),
        (12, 13, 18.211414, 1, 1),
        (14, 15, 21.745663, 4, 0),
        (16, 17, 32.911320, 2, 1),
        (18, 19, 33.409951, 5, 0),
        (20, 21, 47.253259, 6, 0),
        (22, 22, 50.310374, 0, 2),
    ),
    "ann_ba05_nu030": (
        (7, 8, 4.264461, 2, 0),
        (9, 9, 9.311375, 0, 1),
        (10, 11, 11.408169, 3, 0),
        (12, 13, 17.097298, 1, 1),
        (14, 15, 21.029607, 4, 0),
        (16, 17, 30.934933, 2, 1),
        (18, 19, 32.910338, 5, 0),
        (20, 21, 46.941487, 6, 0),
        (22, 23, 47.245036, 3, 1),
    ),
    "ann_ba07_nu030": (
        (7, 8, 3.565849, 2, 0),
        (9, 10, 9.844975, 3, 0),
        (11, 11, 13.155562, 0, 1),
        (12, 13, 18.668707, 4, 0),
        (14, 15, 21.738907, 1, 1),
        (16, 17, 29.971792, 5, 0),
        (18, 19, 37.447912, 2, 1),
        (20, 21, 43.692941, 6, 0),
        (22, 23, 55.109742, 3, 1),
    ),
    "ann_ba09_nu030": (
        (7, 8, 2.919630, 2, 0),
        (9, 10, 8.123090, 3, 0),
        (11, 12, 15.474280, 4, 0),
        (13, 14, 24.948764, 5, 0),
        (15, 16, 28.217476, None, None),  # IP; |UZ| ~ 1e-13
        (17, 17, 34.659536, 0, 1),
        (18, 19, 36.541315, 6, 0),
        (20, 21, 50.251458, 7, 0),
        (22, 23, 54.353036, 1, 1),
    ),
    "disk_nu033": (
        (7, 8, 5.258423, 2, 0),
        (9, 9, 9.067930, 0, 1),
        (10, 11, 12.228422, 3, 0),
        (12, 13, 20.501026, 1, 1),
        (14, 15, 21.488967, 4, 0),
        (16, 17, 32.987120, 5, 0),
        (18, 19, 35.198791, 2, 1),
        (20, 20, 38.482457, 0, 2),
        (21, 22, 46.680363, 6, 0),
    ),
    # modes 7 and 9 counted; 10-22 same degeneracy ladder as disk_nu033
    "disk_nu030": (
        (7, 8, 5.354847, 2, 0),
        (9, 9, 9.002238, 0, 1),
        (10, 11, 12.424200, 3, 0),
        (12, 13, 20.462445, 1, 1),
        (14, 15, 21.798798, 4, 0),
        (16, 17, 33.424149, 5, 0),
        (18, 19, 35.216048, 2, 1),
        (20, 20, 38.419658, 0, 2),
        (21, 22, 47.256712, 6, 0),
    ),
}


def dump_json(path, payload):
    tmp = path + ".tmp"
    with open(tmp, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2, sort_keys=True)
        fh.write("\n")
    os.replace(tmp, path)


def native_from_irie(lam, geom, mat):
    factor = inplane_lambda_irie(1.0, geom, mat)
    if factor == 0.0:
        raise ZeroDivisionError("Irie conversion factor is zero")
    return float(lam) / factor
