#!/usr/bin/env python3
"""
Off-cluster MAC (Modal Assurance Criterion) cross-check for Paper 2
rectangular FFFF out-of-plane (flexural) modes.

FE side (already on disk, job 2456738, all five QUEUE_OK):
  Ansys/NewAnsys/rect_ff_oop_mac_{symm,asym}_freqs_lob<tag>.txt
  Ansys/NewAnsys/rect_ff_oop_mac_{symm,asym}_m<1..20>_lob<tag>.txt
  nodal columns: X, Y, UZ. Half-model X in [-L,L], Y in [0,B], B=1, L=LOB*B.

Analytical side: this project's deployed RectOOPAssembler(bc='free_free')
null-vector reconstructed at the published Lambda* and evaluated at the
same FE nodes. Shape functions are the same U(xi,r,m,x1)*peval(psi[0],x2)
expansion already used by probe_rect_ff_oop_spurious_screen_v1 (w at
corners) and by archive residual-screen reconstruct_w_M_V. Null vector
via realify + equilibrated_nullvec_mp + invert_realify (the spike-free
path in detectors.py, not a float64 SVD).

MAC (same definition as probe_mac_ambiguity_check / the nu=0.35 in-plane
bake-off, adapted to a scalar UZ field):
  MAC = |dot(uz_fe, w_an)|^2 / (||uz_fe||^2 ||w_an||^2)
on the real part of the phase-aligned analytical field.

KFAC: convert FE Hz in Python with KFAC=24.6644. Do NOT trust the
deck-printed Lambda_FE column (ansys-kfac-runtime-bug Instance 3: the
*VWRITE of LAM(1) after a *DO fill still prints f / (pi*H/8) instead of
f / 24.6644 -- confirmed on these files, e.g. SYM l/b=1 mode 4 is
48.897 Hz -> printed 3112.88, Python 1.98249).

SOLVER_VERSION is not bumped. This is a cross-check, not a new capability.

PRE-REGISTERED (do not retune 0.744 after seeing numbers)
--------------------------------------------------------
Per published Tables 1-2 pair, against the frequency-matched same-parity
FE mode (Hz/KFAC nearest to the paper's Lambda_FE):
  CONFIRMED      MAC >= 0.744  (Paper 1 MAC bar, reused not retuned)
  AMBIGUOUS      0.3 <= MAC < 0.744
  SHAPE_MISMATCH MAC < 0.3
  PAIRING_SWAP   frequency partner MAC < 0.744 but some other same-parity
                 FE mode has MAC >= 0.744
  D_GATE         reconstruction failed / empty FE file / unused-node
                 collapse left too few points
Negative-control leaks (SYM 0.41, ANTI 0.46; not in Tables 1-2) are
reported separately and are NOT used to move the bar.

USAGE
  python -u probe_rect_ff_oop_mac_2026-09-05.py
  SMOKE=1   -> only (1.5, SYM, 0.964) and (1.0, ANTI, 1.366)
  ONLY=...  -> comma list of 'lob:CLASS:Lam' (e.g. 1.5:SYM:0.964)
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.environ.get("PKG_PATH", "."))
os.environ.setdefault("DPS", "30")

from mpmath import mp, mpf, mpc  # noqa: E402

mp.dps = int(os.environ["DPS"])

import plate_solver as ps  # noqa: E402
from plate_solver.detectors import (  # noqa: E402
    equilibrated_nullvec_mp,
    rect_resolve_branches,
    rect_select_branches,
)
import p5_rect_ff_lib as L  # noqa: E402

EXPECT_VER = os.environ.get("EXPECT_SOLVER_VERSION", "2026-07-10.s10")
ANSYS_DIR = os.environ.get(
    "ANSYS_DIR",
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "Ansys", "NewAnsys"),
)
OUT_JSON = os.environ.get(
    "MAC_JSON",
    os.path.join(os.path.dirname(os.path.abspath(__file__)),
                 "rect_ff_oop_mac_2026-09-05.json"),
)
OUT_MD = os.environ.get(
    "MAC_MD",
    os.path.join(os.path.dirname(os.path.abspath(__file__)),
                 "..", "github_repo", "paper", "PAPER2_OOP_MAC_TABLE.md"),
)
MAC_BAR = 0.744
MAC_LO = 0.3
SCALE = 0.5 * np.pi  # x1 = X * pi/(2B), B=1 -> X * pi/2; same for x2/Y

# Frozen from github_repo/paper/PAPER2_RECT_FREEFREE_DRAFT.tex Tables 1-2.
# persist=True rows never appear as production-basis dips.
PAPER_PAIRS = [
    # Table 1 primary (Lambda ~< 2.4)
    dict(table=1, lob=1.0, sym=False, Lam=1.366, Lfe=1.35288, note=""),
    dict(table=1, lob=1.5, sym=True,  Lam=0.964, Lfe=0.96347, note="1st SYM"),
    dict(table=1, lob=1.5, sym=True,  Lam=2.254, Lfe=2.24373, note=""),
    dict(table=1, lob=1.5, sym=False, Lam=0.906, Lfe=0.89838, note="1st ANTI"),
    dict(table=1, lob=1.5, sym=False, Lam=2.124, Lfe=2.07036, note="weakest primary"),
    dict(table=1, lob=2.0, sym=True,  Lam=0.5434, Lfe=0.54338, note="persist basis only", persist=True),
    dict(table=1, lob=2.0, sym=True,  Lam=1.510, Lfe=1.50735, note=""),
    dict(table=1, lob=2.0, sym=True,  Lam=2.234, Lfe=2.22558, note=""),
    dict(table=1, lob=2.0, sym=False, Lam=0.674, Lfe=0.66858, note="1st ANTI"),
    dict(table=1, lob=2.0, sym=False, Lam=1.482, Lfe=1.47069, note=""),
    dict(table=1, lob=2.5, sym=True,  Lam=0.348, Lfe=0.34766, note="1st SYM"),
    dict(table=1, lob=2.5, sym=True,  Lam=0.978, Lfe=0.96559, note=""),
    dict(table=1, lob=2.5, sym=True,  Lam=1.888, Lfe=1.88370, note=""),
    dict(table=1, lob=2.5, sym=True,  Lam=2.278, Lfe=2.27059, note=""),
    dict(table=1, lob=2.5, sym=False, Lam=0.534, Lfe=0.53126, note="1st ANTI"),
    dict(table=1, lob=2.5, sym=False, Lam=1.148, Lfe=1.14003, note=""),
    dict(table=1, lob=2.5, sym=False, Lam=1.918, Lfe=1.90392, note=""),
    dict(table=1, lob=3.0, sym=True,  Lam=0.242, Lfe=0.24127, note="1st SYM"),
    dict(table=1, lob=3.0, sym=True,  Lam=0.6698, Lfe=0.66984, note="persist basis only", persist=True),
    dict(table=1, lob=3.0, sym=True,  Lam=1.320, Lfe=1.31785, note=""),
    dict(table=1, lob=3.0, sym=True,  Lam=2.162, Lfe=2.15417, note=""),
    dict(table=1, lob=3.0, sym=True,  Lam=2.256, Lfe=2.24841, note=""),
    dict(table=1, lob=3.0, sym=False, Lam=0.444, Lfe=0.44044, note="1st ANTI; 0.46-family coincidence"),
    dict(table=1, lob=3.0, sym=False, Lam=0.936, Lfe=0.93044, note=""),
    dict(table=1, lob=3.0, sym=False, Lam=1.528, Lfe=1.51698, note=""),
    dict(table=1, lob=3.0, sym=False, Lam=2.258, Lfe=2.24305, note=""),
    # Table 2 higher
    dict(table=2, lob=1.0, sym=True,  Lam=2.460, Lfe=2.45439, note=""),
    dict(table=2, lob=1.0, sym=True,  Lam=3.530, Lfe=3.49250, note=""),
    dict(table=2, lob=1.0, sym=True,  Lam=6.190, Lfe=6.16026, note=""),
    dict(table=2, lob=1.0, sym=True,  Lam=6.460, Lfe=6.36553, note=""),
    dict(table=2, lob=1.0, sym=False, Lam=3.520, Lfe=3.49250, note=""),
    dict(table=2, lob=1.5, sym=True,  Lam=5.380, Lfe=5.38152, note=""),
    dict(table=2, lob=1.5, sym=False, Lam=3.860, Lfe=3.83354, note=""),
    dict(table=2, lob=2.0, sym=True,  Lam=3.000, Lfe=2.99812, note=""),
    dict(table=2, lob=2.0, sym=True,  Lam=3.650, Lfe=3.62227, note=""),
    dict(table=2, lob=2.0, sym=False, Lam=2.580, Lfe=2.55180, note=""),
    dict(table=2, lob=2.0, sym=False, Lam=4.060, Lfe=4.02613, note=""),
    dict(table=2, lob=2.0, sym=False, Lam=6.200, Lfe=6.17866, note=""),
    dict(table=2, lob=2.5, sym=True,  Lam=3.190, Lfe=3.19522, note=""),
    dict(table=2, lob=2.5, sym=True,  Lam=4.180, Lfe=4.13819, note=""),
    dict(table=2, lob=2.5, sym=True,  Lam=5.410, Lfe=5.34819, note=""),
    dict(table=2, lob=2.5, sym=False, Lam=2.920, Lfe=2.89206, note=""),
    dict(table=2, lob=2.5, sym=False, Lam=4.180, Lfe=4.14962, note=""),
    dict(table=2, lob=2.5, sym=False, Lam=5.720, Lfe=5.67275, note=""),
    dict(table=2, lob=2.5, sym=False, Lam=6.260, Lfe=6.21965, note=""),
    dict(table=2, lob=2.5, sym=False, Lam=6.400, Lfe=6.36131, note=""),
    dict(table=2, lob=3.0, sym=True,  Lam=2.480, Lfe=2.45935, note=""),
    dict(table=2, lob=3.0, sym=True,  Lam=3.340, Lfe=3.32528, note=""),
    dict(table=2, lob=3.0, sym=True,  Lam=3.660, Lfe=3.62467, note=""),
    dict(table=2, lob=3.0, sym=True,  Lam=5.660, Lfe=5.59746, note=""),
    dict(table=2, lob=3.0, sym=True,  Lam=6.180, Lfe=6.14744, note=""),
    dict(table=2, lob=3.0, sym=False, Lam=3.160, Lfe=3.14253, note=""),
    dict(table=2, lob=3.0, sym=False, Lam=4.260, Lfe=4.23529, note=""),
    dict(table=2, lob=3.0, sym=False, Lam=5.560, Lfe=5.51872, note=""),
    dict(table=2, lob=3.0, sym=False, Lam=6.220, Lfe=6.18191, note=""),
    dict(table=2, lob=3.0, sym=False, Lam=6.380, Lfe=6.34339, note=""),
]

# Known Screen-B leaks. Not in Tables 1-2. Negative control only.
LEAK_PAIRS = [
    dict(table=0, lob=1.0, sym=True,  Lam=0.410, Lfe=None, note="SYM 0.41 leak"),
    dict(table=0, lob=1.5, sym=True,  Lam=0.410, Lfe=None, note="SYM 0.41 leak"),
    dict(table=0, lob=2.5, sym=True,  Lam=0.410, Lfe=None, note="SYM 0.41 leak"),
    dict(table=0, lob=1.0, sym=False, Lam=0.460, Lfe=None, note="ANTI 0.46 leak"),
    dict(table=0, lob=1.5, sym=False, Lam=0.460, Lfe=None, note="ANTI 0.46 leak"),
    dict(table=0, lob=2.5, sym=False, Lam=0.460, Lfe=None, note="ANTI 0.46 leak"),
]


def lob_tag(lob):
    return f"{int(round(lob * 100)):03d}"


def parity_token(sym):
    return "symm" if sym else "asym"


def mac(a, b):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    num = abs(float(np.dot(a, b))) ** 2
    den = float(np.dot(a, a)) * float(np.dot(b, b))
    if den <= 0:
        return float("nan")
    return num / den


def _conjugate_pair_cols(full, tol=1e-9):
    pairs = []
    used = set()
    for p, z in enumerate(full):
        if p in used or abs(z.imag) < tol:
            continue
        for q in range(p + 1, len(full)):
            if q in used:
                continue
            w = full[q]
            if abs(w.real - z.real) < tol and abs(w.imag + z.imag) < tol:
                used.add(p)
                used.add(q)
                for r in (0, 1):
                    pairs.append((2 * p + r, 2 * q + r))
                break
    return pairs


def invert_realify(Y, full):
    """X = U Y. Copied from archive residual-screen invert_realify."""
    n = len(Y)
    X = [None] * n
    s2 = mp.sqrt(2)
    I = mpc(0, 1)
    paired = set()
    for (a, b) in _conjugate_pair_cols(full):
        Ya = mpf(Y[a])
        Yb = mpf(Y[b])
        X[a] = (Ya - I * Yb) / s2
        X[b] = (Ya + I * Yb) / s2
        paired.add(a)
        paired.add(b)
    for i in range(n):
        if i not in paired:
            X[i] = mpc(Y[i])
    return X


def load_freq_file(path):
    """Return list of (mode_index, f_Hz, Lambda_python). Ignore printed Lambda."""
    rows = []
    with open(path, encoding="utf-8", errors="replace") as f:
        lines = f.readlines()
    for line in lines[1:]:
        parts = line.split()
        if len(parts) < 2:
            continue
        try:
            midx = int(float(parts[0]))
            fhz = float(parts[1])
        except ValueError:
            continue
        rows.append((midx, fhz, fhz / L.KFAC))
    return rows


def load_mode_file(path):
    """X, Y, UZ arrays, unused ANSYS slots (0,0,0 duplicates) dropped."""
    xs, ys, uz = [], [], []
    with open(path, encoding="utf-8", errors="replace") as f:
        header = True
        for line in f:
            if header:
                header = False
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            xs.append(float(parts[0]))
            ys.append(float(parts[1]))
            uz.append(float(parts[2]))
    x = np.asarray(xs, dtype=float)
    y = np.asarray(ys, dtype=float)
    u = np.asarray(uz, dtype=float)
    key = np.round(x, 6) + 1j * np.round(y, 6)
    _, idx = np.unique(key, return_index=True)
    idx = np.sort(idx)
    return x[idx], y[idx], u[idx], int(x.size), int(idx.size)


def load_fe_block(lob, sym):
    tag = lob_tag(lob)
    tok = parity_token(sym)
    freq_path = os.path.join(ANSYS_DIR, f"rect_ff_oop_mac_{tok}_freqs_lob{tag}.txt")
    if not os.path.isfile(freq_path):
        raise FileNotFoundError(freq_path)
    freqs = load_freq_file(freq_path)
    modes = []
    x0 = y0 = None
    n_raw = n_keep = None
    for midx, fhz, lam in freqs:
        p = os.path.join(ANSYS_DIR, f"rect_ff_oop_mac_{tok}_m{midx}_lob{tag}.txt")
        if not os.path.isfile(p):
            raise FileNotFoundError(p)
        x, y, u, raw, keep = load_mode_file(p)
        if x0 is None:
            x0, y0, n_raw, n_keep = x, y, raw, keep
        else:
            if x.shape != x0.shape or np.max(np.abs(x - x0)) > 1e-8 \
                    or np.max(np.abs(y - y0)) > 1e-8:
                raise RuntimeError(f"node mismatch in {p}")
        modes.append(dict(midx=midx, fhz=fhz, lam=lam, uz=u))
    return dict(x=x0, y=y0, n_raw=n_raw, n_keep=n_keep, modes=modes,
                freq_path=freq_path)


def reconstruct(asm, sym, lob, Lam, persist):
    n_real, n_cpair = L.n_basis(sym, persist=persist)
    im_cap = 30.0 if persist else 7.0
    reps = rect_resolve_branches(asm.eng(sym), Lam, im_cap=im_cap)
    full = rect_select_branches(reps, n_real=n_real, n_cpair=n_cpair)
    if len(full) < 2:
        return dict(ok=False, reason="no branches", n_used=len(full))
    K = asm.assemble(full, mpf(str(round(Lam, 6))), sym, lob)
    sigma = float(asm.equil_sigma(K.copy(), full))
    Kr = K.copy()
    asm._realify_conjugate_pairs(Kr, full)
    n = Kr.rows
    max_im = float(max(abs(mp.im(Kr[i, j])) for i in range(n) for j in range(n)))
    max_re = float(max(abs(mp.re(Kr[i, j])) for i in range(n) for j in range(n)))
    Kr_real = mp.matrix(n, n)
    for i in range(n):
        for j in range(n):
            Kr_real[i, j] = mp.re(Kr[i, j])
    Y = equilibrated_nullvec_mp(Kr_real, n)
    if Y is None:
        return dict(ok=False, reason="no null vector", n_used=len(full),
                    sigma=sigma, max_im=max_im, max_re=max_re)
    X = invert_realify(Y, full)
    phin = mp.pi / 2 if sym else mpf(0)
    Lf = float(Lam)
    eng = asm.eng(sym)
    branches = []
    for xi in full:
        e1, e2 = eng._etas(complex(xi), Lf)
        H1, H2 = eng._amp_ratio(complex(xi), Lf)
        branches.append((
            complex(xi),
            complex(e1), complex(e2),
            complex(H1), complex(H2),
        ))
    Xc = np.array([complex(x) for x in X], dtype=np.complex128)
    return dict(
        ok=True, reason="", n_used=len(full), n_real=n_real, n_cpair=n_cpair,
        im_cap=im_cap, sigma=sigma, max_im=max_im, max_re=max_re,
        X=Xc, branches=branches, phin=float(phin),
    )


def eval_w(rec, x1, x2):
    """Vectorized w(x1,x2) from reconstruct() output."""
    w = np.zeros(x1.shape, dtype=np.complex128)
    phin = rec["phin"]
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        for p, (xi, e1, e2, H1, H2) in enumerate(rec["branches"]):
            p0 = H1 * np.sin(e1 * x2 + phin) + H2 * np.sin(e2 * x2 + phin)
            for r in (1, 2):
                c = rec["X"][2 * p + (r - 1)]
                U0 = np.sin(xi * x1 + (r - 1) * 0.5 * np.pi)
                w += c * p0 * U0
    w = np.nan_to_num(w, nan=0.0, posinf=0.0, neginf=0.0)
    return w


def phase_align_real(w):
    """Global phase so the field is as real as possible; return (w_real, im_ratio)."""
    energy = np.sum(w * w)
    if abs(energy) == 0:
        wr = np.real(w)
        denom = float(np.max(np.abs(w))) if np.max(np.abs(w)) > 0 else 1.0
        return wr, float(np.max(np.abs(np.imag(w)))) / denom
    phi = 0.5 * np.angle(energy)
    wr = w * np.exp(-1j * phi)
    mag = float(np.max(np.abs(wr)))
    if mag <= 0:
        mag = 1.0
    return np.real(wr), float(np.max(np.abs(np.imag(wr)))) / mag


def classify(mac_fe, mac_best, best_is_fe):
    if not np.isfinite(mac_fe) or not np.isfinite(mac_best):
        return "D_GATE"
    if mac_fe >= MAC_BAR:
        return "CONFIRMED"
    if (not best_is_fe) and mac_best >= MAC_BAR:
        return "PAIRING_SWAP"
    if mac_fe >= MAC_LO:
        return "AMBIGUOUS"
    return "SHAPE_MISMATCH"


def pair_key(p):
    sl = "SYM" if p["sym"] else "ANTI"
    return f"{p['lob']:.1f}:{sl}:{p['Lam']}"


def select_pairs():
    pairs = list(PAPER_PAIRS)
    if os.environ.get("SKIP_LEAKS", "0") != "1":
        pairs = pairs + LEAK_PAIRS
    if os.environ.get("SMOKE", "0") == "1":
        keep = {"1.5:SYM:0.964", "1.0:ANTI:1.366"}
        pairs = [p for p in pairs if pair_key(p) in keep]
    only = os.environ.get("ONLY", "").strip()
    if only:
        want = {s.strip() for s in only.split(",") if s.strip()}
        pairs = [p for p in pairs if pair_key(p) in want]
    return pairs


def write_md(results, path):
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    lines = []
    lines.append("# Paper 2 rectangular FFFF OOP MAC table")
    lines.append("")
    lines.append("Off-cluster computation (`probe_rect_ff_oop_mac_2026-09-05.py`).")
    lines.append("FE eigenvectors from job 2456738 (five SHELL281 MAC-extract")
    lines.append("decks). Analytical shapes from `RectOOPAssembler(bc='free_free')`")
    lines.append("at the published $\\Lambda^*$ with the same production / persist")
    lines.append("basis that produced that $\\Lambda^*$. `SOLVER_VERSION` unmodified.")
    lines.append("")
    lines.append("MAC $= |\\langle u_z^{\\mathrm{FE}}, w^{\\mathrm{an}}\\rangle|^2 / "
                 "(\\|u_z^{\\mathrm{FE}}\\|^2 \\|w^{\\mathrm{an}}\\|^2)$ on the")
    lines.append("phase-aligned real part of $w$. Bar 0.744 is Paper 1's, not retuned.")
    lines.append("FE $\\Lambda$ is $f_{\\mathrm{Hz}}/24.6644$ (the deck-printed")
    lines.append("`Lambda_FE` column is the KFAC/*VWRITE bug and is ignored).")
    lines.append("")
    lines.append("## Tables 1–2 matches")
    lines.append("")
    lines.append("| tab | $\\ell/b$ | class | $\\Lambda^*$ | $\\Lambda_{\\mathrm{FE}}$ "
                 "| FE# | MAC(FE) | best FE# | MAC(best) | 2nd MAC | im/re | letter | note |")
    lines.append("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    paper = [r for r in results if r.get("table") in (1, 2)]
    leaks = [r for r in results if r.get("table") == 0]
    for r in paper:
        sl = "SYM" if r["sym"] else "ANTI"
        lines.append(
            "| {table} | {lob:.1f} | {sl} | {Lam} | {Lfe} | {fe_midx} | {mac_fe:.4f} "
            "| {best_midx} | {mac_best:.4f} | {mac_2nd:.4f} | {im_ratio:.3f} | {letter} | {note} |".format(
                table=r["table"], lob=r["lob"], sl=sl, Lam=r["Lam"], Lfe=r["Lfe"],
                fe_midx=r.get("fe_midx", ""), mac_fe=r.get("mac_fe", float("nan")),
                best_midx=r.get("best_midx", ""), mac_best=r.get("mac_best", float("nan")),
                mac_2nd=r.get("mac_2nd", float("nan")), im_ratio=r.get("im_ratio", float("nan")),
                letter=r.get("letter", "D_GATE"), note=r.get("note", ""),
            )
        )
    if leaks:
        lines.append("")
        lines.append("## Negative-control leaks (not in Tables 1–2)")
        lines.append("")
        lines.append("| $\\ell/b$ | class | $\\Lambda^*$ | best FE# | $\\Lambda_{\\mathrm{best}}$ "
                     "| MAC(best) | 2nd MAC | letter | note |")
        lines.append("|---|---|---|---|---|---|---|---|---|")
        for r in leaks:
            sl = "SYM" if r["sym"] else "ANTI"
            lines.append(
                "| {lob:.1f} | {sl} | {Lam} | {best_midx} | {best_lam} | {mac_best:.4f} "
                "| {mac_2nd:.4f} | {letter} | {note} |".format(
                    lob=r["lob"], sl=sl, Lam=r["Lam"],
                    best_midx=r.get("best_midx", ""),
                    best_lam=r.get("best_lam", float("nan")),
                    mac_best=r.get("mac_best", float("nan")),
                    mac_2nd=r.get("mac_2nd", float("nan")),
                    letter=r.get("letter", "D_GATE"), note=r.get("note", ""),
                )
            )
    n_conf = sum(1 for r in paper if r.get("letter") == "CONFIRMED")
    n_cp = sum(1 for r in paper if r.get("letter") == "CONFIRMED_PERSIST")
    n_swap = sum(1 for r in paper if r.get("letter") == "PAIRING_SWAP")
    n_amb = sum(1 for r in paper if r.get("letter") == "AMBIGUOUS")
    n_mis = sum(1 for r in paper if r.get("letter") == "SHAPE_MISMATCH")
    n_gate = sum(1 for r in paper if r.get("letter") == "D_GATE")
    lines.append("")
    lines.append("## Summary")
    lines.append("")
    lines.append(f"- Tables 1–2 pairs scored: {len(paper)}")
    lines.append(f"- CONFIRMED (MAC ≥ {MAC_BAR} vs frequency partner, production basis): {n_conf}")
    lines.append(f"- CONFIRMED_PERSIST (same bar, Screen B n_cpair=3 basis): {n_cp}")
    lines.append(f"- PAIRING_SWAP remaining: {n_swap}")
    lines.append(f"- AMBIGUOUS remaining: {n_amb}")
    lines.append(f"- SHAPE_MISMATCH remaining: {n_mis}")
    lines.append(f"- D_GATE: {n_gate}")
    if paper:
        macs = [r["mac_fe"] for r in paper if np.isfinite(r.get("mac_fe", float("nan")))]
        if macs:
            lines.append(f"- MAC(FE) min / median / max: "
                         f"{min(macs):.4f} / {float(np.median(macs)):.4f} / {max(macs):.4f}")
    lines.append("")
    lines.append("Paper text was not edited. Fold-in proposal lives in the")
    lines.append("session report, not in this file.")
    lines.append("")
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(lines))


def main():
    t0 = time.time()
    print("=" * 78)
    print("  probe_rect_ff_oop_mac_2026-09-05 -- start "
          + time.strftime("%Y-%m-%d %H:%M:%S"))
    print("=" * 78, flush=True)
    pkg = os.environ.get("PKG_PATH", ".")
    if not L.deploy_ok(pkg, EXPECT_VER):
        raise SystemExit(2)
    print(f"preflight OK: SOLVER_VERSION={ps.SOLVER_VERSION}  "
          f"dps={mp.dps}  ANSYS_DIR={ANSYS_DIR}", flush=True)

    pairs = select_pairs()
    print(f"pairs to score: {len(pairs)}", flush=True)
    asm = L.make_ff()
    cache = {}
    results = []

    for i, p in enumerate(pairs, 1):
        sl = "SYM" if p["sym"] else "ANTI"
        persist = bool(p.get("persist", False))
        key = (p["lob"], p["sym"])
        print(f"\n[{i}/{len(pairs)}] l/b={p['lob']:.1f} {sl}  "
              f"Lam*={p['Lam']}  Lfe={p['Lfe']}  persist={persist}",
              flush=True)
        t1 = time.time()
        row = dict(p)
        row.setdefault("persist", persist)
        try:
            if key not in cache:
                cache[key] = load_fe_block(p["lob"], p["sym"])
                blk = cache[key]
                print(f"  FE nodes raw={blk['n_raw']} keep={blk['n_keep']}  "
                      f"n_modes={len(blk['modes'])}", flush=True)
            blk = cache[key]
            rec = reconstruct(asm, p["sym"], p["lob"], p["Lam"], persist)
            if not rec["ok"]:
                row.update(letter="D_GATE", reason=rec.get("reason", ""),
                           n_used=rec.get("n_used", 0),
                           mac_fe=float("nan"), mac_best=float("nan"),
                           mac_2nd=float("nan"), im_ratio=float("nan"))
                print(f"  D_GATE: {rec.get('reason')}", flush=True)
                results.append(row)
                continue
            x1 = blk["x"] * SCALE
            x2 = blk["y"] * SCALE
            w = eval_w(rec, x1, x2)
            w_real, im_ratio = phase_align_real(w)
            macs = []
            for m in blk["modes"]:
                if m["lam"] < 0.2:
                    continue
                mv = mac(m["uz"], w_real)
                macs.append((mv, m["midx"], m["lam"], m["fhz"]))
            macs.sort(key=lambda t: -t[0])
            mac_best = macs[0][0] if macs else float("nan")
            best_midx = macs[0][1] if macs else None
            best_lam = macs[0][2] if macs else float("nan")
            mac_2nd = macs[1][0] if len(macs) > 1 else float("nan")
            fe_midx = None
            mac_fe = float("nan")
            fe_lam = None
            if p["Lfe"] is not None:
                partner = min(blk["modes"], key=lambda m: abs(m["lam"] - p["Lfe"]))
                fe_midx = partner["midx"]
                fe_lam = partner["lam"]
                mac_fe = mac(partner["uz"], w_real)
                dlam = abs(fe_lam - p["Lfe"])
                if dlam > 0.01:
                    print(f"  WARN: nearest FE Lambda {fe_lam:.5f} is "
                          f"{dlam:.4f} from paper Lfe {p['Lfe']}", flush=True)
            letter = classify(mac_fe, mac_best, best_midx == fe_midx) \
                if p["Lfe"] is not None else (
                    "LEAK_HIGH" if (np.isfinite(mac_best) and mac_best >= MAC_BAR)
                    else ("LEAK_AMB" if (np.isfinite(mac_best) and mac_best >= MAC_LO)
                          else "LEAK_LOW")
                )
            row.update(
                letter=letter, reason="", n_used=rec["n_used"],
                n_real=rec["n_real"], n_cpair=rec["n_cpair"],
                im_cap=rec["im_cap"], sigma=rec["sigma"],
                max_im_K=rec["max_im"], im_ratio=im_ratio,
                mac_fe=float(mac_fe), fe_midx=fe_midx, fe_lam=fe_lam,
                mac_best=float(mac_best), best_midx=best_midx,
                best_lam=float(best_lam), mac_2nd=float(mac_2nd),
                n_keep=blk["n_keep"], dt=time.time() - t1,
            )
            print(
                f"  sigma={rec['sigma']:.3e}  n_used={rec['n_used']}  "
                f"im/re={im_ratio:.3f}  MAC(FE#{fe_midx})={mac_fe:.4f}  "
                f"best=FE#{best_midx} {mac_best:.4f}  2nd={mac_2nd:.4f}  "
                f"{letter}  [{time.time()-t1:.1f}s]",
                flush=True,
            )
        except Exception as e:
            row.update(letter="D_GATE", reason=repr(e),
                       mac_fe=float("nan"), mac_best=float("nan"),
                       mac_2nd=float("nan"), im_ratio=float("nan"))
            print(f"  D_GATE exception: {e!r}", flush=True)
        results.append(row)

    # Second pass: persist-basis reconstruction for Table 1-2 rows that
    # did not CONFIRMED at the production basis. Pre-registered: this is
    # the Screen B basis (n_cpair=3, im_cap=30), not a new instrument and
    # not a retune of 0.744. If MAC(FE partner) then exceeds the bar, the
    # frequency match was real and the production-basis eigenvector was
    # just truncated. Leaks are not retried.
    if os.environ.get("PERSIST_RETRY", "1") == "1":
        retry = [r for r in results
                 if r.get("table") in (1, 2)
                 and r.get("letter") != "CONFIRMED"
                 and not r.get("persist")]
        print(f"\npersist-retry of {len(retry)} non-CONFIRMED production rows",
              flush=True)
        for j, row in enumerate(retry, 1):
            sl = "SYM" if row["sym"] else "ANTI"
            print(f"\n[retry {j}/{len(retry)}] l/b={row['lob']:.1f} {sl}  "
                  f"Lam*={row['Lam']}", flush=True)
            t1 = time.time()
            try:
                blk = cache[(row["lob"], row["sym"])]
                rec = reconstruct(asm, row["sym"], row["lob"], row["Lam"], True)
                if not rec["ok"]:
                    row["persist_retry"] = dict(
                        letter="D_GATE", reason=rec.get("reason", ""))
                    print(f"  persist D_GATE: {rec.get('reason')}", flush=True)
                    continue
                x1 = blk["x"] * SCALE
                x2 = blk["y"] * SCALE
                w_real, im_ratio = phase_align_real(eval_w(rec, x1, x2))
                macs = []
                for m in blk["modes"]:
                    if m["lam"] < 0.2:
                        continue
                    macs.append((mac(m["uz"], w_real), m["midx"], m["lam"]))
                macs.sort(key=lambda t: -t[0])
                mac_best = macs[0][0] if macs else float("nan")
                best_midx = macs[0][1] if macs else None
                best_lam = macs[0][2] if macs else float("nan")
                mac_2nd = macs[1][0] if len(macs) > 1 else float("nan")
                partner = None
                mac_fe = float("nan")
                fe_midx = row.get("fe_midx")
                if row.get("Lfe") is not None:
                    partner = min(blk["modes"],
                                  key=lambda m: abs(m["lam"] - row["Lfe"]))
                    fe_midx = partner["midx"]
                    mac_fe = mac(partner["uz"], w_real)
                letter = classify(mac_fe, mac_best, best_midx == fe_midx)
                row["persist_retry"] = dict(
                    letter=letter, mac_fe=float(mac_fe), fe_midx=fe_midx,
                    mac_best=float(mac_best), best_midx=best_midx,
                    best_lam=float(best_lam), mac_2nd=float(mac_2nd),
                    im_ratio=im_ratio, sigma=rec["sigma"],
                    n_used=rec["n_used"], n_cpair=rec["n_cpair"],
                    dt=time.time() - t1,
                )
                # Promote the reported letter if persist recovers the partner.
                if letter == "CONFIRMED":
                    row["letter"] = "CONFIRMED_PERSIST"
                    row["mac_fe"] = float(mac_fe)
                    row["mac_best"] = float(mac_best)
                    row["mac_2nd"] = float(mac_2nd)
                    row["best_midx"] = best_midx
                    row["im_ratio"] = im_ratio
                    row["note"] = (row.get("note") or "") + " persist-basis shape"
                print(
                    f"  persist sigma={rec['sigma']:.3e} n_used={rec['n_used']}  "
                    f"MAC(FE#{fe_midx})={mac_fe:.4f}  best=FE#{best_midx} "
                    f"{mac_best:.4f}  {letter}  [{time.time()-t1:.1f}s]",
                    flush=True,
                )
            except Exception as e:
                row["persist_retry"] = dict(letter="D_GATE", reason=repr(e))
                print(f"  persist D_GATE exception: {e!r}", flush=True)

    payload = dict(
        solver_version=ps.SOLVER_VERSION,
        kfac=L.KFAC, mac_bar=MAC_BAR, ansys_dir=ANSYS_DIR,
        job="2456738", n_pairs=len(results),
        results=results, elapsed_s=time.time() - t0,
    )
    with open(OUT_JSON, "w", encoding="utf-8", newline="\n") as f:
        json.dump(payload, f, indent=2, default=str)
    write_md(results, OUT_MD)
    paper = [r for r in results if r.get("table") in (1, 2)]
    print("\n" + "=" * 78)
    print(f"  wrote {OUT_JSON}")
    print(f"  wrote {OUT_MD}")
    print(f"  Tables 1-2: {len(paper)}  "
          f"CONFIRMED={sum(r.get('letter')=='CONFIRMED' for r in paper)}  "
          f"CONFIRMED_PERSIST={sum(r.get('letter')=='CONFIRMED_PERSIST' for r in paper)}  "
          f"SWAP={sum(r.get('letter')=='PAIRING_SWAP' for r in paper)}  "
          f"AMB={sum(r.get('letter')=='AMBIGUOUS' for r in paper)}  "
          f"MIS={sum(r.get('letter')=='SHAPE_MISMATCH' for r in paper)}  "
          f"GATE={sum(r.get('letter')=='D_GATE' for r in paper)}")
    print(f"  elapsed {time.time()-t0:.1f}s")
    print("=" * 78, flush=True)


if __name__ == "__main__":
    main()
