# -*- coding: utf-8 -*-
"""Production two-sided in-plane classification bar (Rank 18).

This is a POST-PROCESSOR. It does not live inside find_modes_sigmin, does
not change a computed Omega, and does not bump SOLVER_VERSION.

Instrument (LESSONS_LEARNED.md §18.22--§18.24):
  1. Score the candidate on the parity block that OWNS A ROOT at its Omega
     (RULE_LOCALMIN_RESID).
  2. Correlate that reconstruction with FE eigenvectors by the identity-
     weighted MAC below -- never mass-weighted.
  3. Take the maximum MAC over FE modes of the SAME symmetry class as the
     root-owning block. Unrestricted nearest-frequency matching is not
     class-safe: it can hand a candidate an opposite-class FE target, in
     which case MAC is analytically 0.0000 (wrong-class signature, not
     "uncorrelated"). At r200/a100 that accident hit t7 and only t7.
  4. Classify with MAC_THRESHOLD, which sits in the open interval
     (MAC_ART_CEILING, MAC_REAL_FLOOR_SEED) that separates all 24
     calibration candidates at three geometries.

q0 is ANTI and q1 is SYMM about the sector bisector at both FF-P1 and
r200/a100. Do not silently invert that mapping.

R40_BC_KIND: any production call that rebuilds a solver through
_worker_bc() MUST export R40_BC_KIND=free_free into the pool environment
before the pool starts. _worker_bc defaults to clamped_free; leaving it
unset silently builds cantilever workers. Do not change that default --
cantilever runs depend on it.
"""
from __future__ import annotations

import os

import numpy as np

# Seed-based combined window (catalogue Omega). Own-root window is the
# same ceiling with a higher floor; the seed figures are the conservative
# production claim because they need no extra refinement step.
MAC_ART_CEILING = 0.5801
MAC_REAL_FLOOR_SEED = 0.9079
MAC_REAL_FLOOR_OWNROOT = 0.9869
MAC_THRESHOLD = 0.744  # midpoint of (0.5801, 0.9079)

# Block index -> bisector class. Confirmed at FF-P1 and r200/a100.
Q0_CLASS = "ANTI"
Q1_CLASS = "SYMM"

# Local-minimum half-width used by RULE_LOCALMIN_RESID (in Omega).
RULE_DELTA = 0.001


def mac(a, b):
    """Identity-weighted MAC of two stacked nodal displacement vectors."""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    num = abs(float(np.dot(a, b))) ** 2
    den = float(np.dot(a, a)) * float(np.dot(b, b))
    if den <= 0:
        return float("nan")
    return num / den


def stacked_uxuy(ux, uy):
    return np.concatenate([np.asarray(ux, dtype=float),
                           np.asarray(uy, dtype=float)])


def classify_mac(value, threshold=MAC_THRESHOLD):
    """Two-sided call. Returns REAL-like / ARTIFACT-like / INVALID.

    Any threshold in (MAC_ART_CEILING, MAC_REAL_FLOOR_SEED) classifies the
    24-candidate calibration set correctly. Production uses the midpoint.
    A NaN MAC is INVALID, not ARTIFACT-like.
    """
    if value != value:  # NaN
        return "INVALID"
    if value >= threshold:
        return "REAL-like"
    return "ARTIFACT-like"


def block_class(block_name):
    """Map q0/q1 to SYMM/ANTI. Unknown names raise."""
    key = str(block_name).strip().lower()
    if key in ("q0", "0", "anti"):
        return Q0_CLASS
    if key in ("q1", "1", "symm"):
        return Q1_CLASS
    raise ValueError(f"unknown parity block {block_name!r}")


def rule_localmin_resid(sig0, sig1, sig0_m, sig0_p, sig1_m, sig1_p,
                        rTyy0, rTyy1, delta=RULE_DELTA):
    """Pick the root-owning block.

    A block owns a root at Omega only if its log10 sigma_min has a local
    minimum there: sig(Om) < sig(Om-d) and sig(Om) < sig(Om+d). Among
    qualifying blocks take the smaller rTyy. If neither qualifies, fall
    back to argmin sigma.

    Returns (pick, how) with pick in {'q0','q1'}.
    """
    del delta  # documented so callers pass the same d they scanned at
    q0_loc = (sig0 < sig0_m) and (sig0 < sig0_p)
    q1_loc = (sig1 < sig1_m) and (sig1 < sig1_p)
    if q0_loc and not q1_loc:
        return "q0", "single LOCAL_MIN"
    if q1_loc and not q0_loc:
        return "q1", "single LOCAL_MIN"
    if q0_loc and q1_loc:
        pick = "q0" if rTyy0 <= rTyy1 else "q1"
        return pick, "both LOCAL_MIN -> resid"
    pick = "q0" if sig0 < sig1 else "q1"
    return pick, "FALLBACK argmin sigma"


def polar(ux, uy, theta_deg):
    th = np.deg2rad(np.asarray(theta_deg, dtype=float))
    c = np.cos(th)
    s = np.sin(th)
    ux = np.asarray(ux, dtype=float)
    uy = np.asarray(uy, dtype=float)
    ur = ux * c + uy * s
    ut = -ux * s + uy * c
    return ur, ut


def _partner_index(r, th, phi_deg=180.0, ndp=4, th_tol=1.0e-3):
    """Pair each node with its bisector image (r, PHI-theta).

    Exact (round-r, round-theta) keys work at PHI=180 (a100) and fail
    on narrower/wider sectors: those decks are meshed in Cartesian
    coordinates, so midside nodes are chord midpoints whose (r,theta)
    do not land on another printed F10.4 key. Fall back to a nearest-
    neighbour search for the Cartesian image of (r, PHI-theta) within
    half a typical element edge.
    """
    r = np.asarray(r, dtype=float)
    th = np.asarray(th, dtype=float)
    n = len(r)
    r_key = np.round(r, ndp)
    th_key = np.round(th, ndp)
    th_mirr = np.round(float(phi_deg) - th, ndp)
    buckets = {}
    for i in range(n):
        buckets.setdefault((float(r_key[i]), float(th_key[i])), []).append(i)
    partner = np.empty(n, dtype=np.int64)
    missing = 0
    for i in range(n):
        hits = buckets.get((float(r_key[i]), float(th_mirr[i])), [])
        if hits:
            partner[i] = hits[0]
        else:
            partner[i] = -1
            missing += 1
    if missing == 0:
        return partner, missing

    th_rad = np.deg2rad(th)
    xy = np.column_stack((r * np.cos(th_rad), r * np.sin(th_rad)))
    th_img = np.deg2rad(float(phi_deg) - th)
    xy_img = np.column_stack((r * np.cos(th_img), r * np.sin(th_img)))
    span = float(np.ptp(r)) if n else 0.0
    max_dist = max(0.5 * span / 96.0, 1.0e-4)
    try:
        from scipy.spatial import cKDTree
        dist, idx = cKDTree(xy).query(xy_img, k=1)
        dist = np.asarray(dist, dtype=float)
        idx = np.asarray(idx, dtype=np.int64)
    except Exception:
        dist = np.full(n, np.inf)
        idx = np.zeros(n, dtype=np.int64)
        for i in range(n):
            d2 = (xy[:, 0] - xy_img[i, 0]) ** 2 + (xy[:, 1] - xy_img[i, 1]) ** 2
            j = int(np.argmin(d2))
            dist[i] = float(np.sqrt(d2[j]))
            idx[i] = j
    missing = 0
    for i in range(n):
        if partner[i] >= 0:
            continue
        if dist[i] <= max_dist:
            partner[i] = int(idx[i])
        else:
            partner[i] = -1
            missing += 1
    return partner, missing


def parity_ratios(ur, ut, partner):
    """Even/odd residuals for ur and ut about the sector bisector.

    even residual = mean |f(th) - f(mirror)| / mean |f|
    odd  residual = mean |f(th) + f(mirror)| / mean |f|
    SYMM: ur even, ut odd. ANTI: ur odd, ut even.
    """
    ok = partner >= 0
    p = partner[ok]
    ur_i, ur_p = ur[ok], ur[p]
    ut_i, ut_p = ut[ok], ut[p]

    def safe_div(num, den_arr):
        den = float(np.mean(np.abs(den_arr)))
        if den <= 0:
            return float("nan")
        return float(np.mean(np.abs(num))) / den

    return (
        safe_div(ur_i - ur_p, ur_i),
        safe_div(ur_i + ur_p, ur_i),
        safe_div(ut_i - ut_p, ut_i),
        safe_div(ut_i + ut_p, ut_i),
    )


def classify_parity(ur_even, ur_odd, ut_even, ut_odd,
                    clear=0.15, margin=3.0):
    """SYMM / ANTI / AMBIGUOUS from the four polar residuals."""
    symm_res = 0.5 * (ur_even + ut_odd)
    anti_res = 0.5 * (ur_odd + ut_even)
    ur_wants_symm = ur_even * margin < ur_odd
    ur_wants_anti = ur_odd * margin < ur_even
    ut_wants_symm = ut_odd * margin < ut_even
    ut_wants_anti = ut_even * margin < ut_odd
    if (symm_res < clear and anti_res > clear and ur_wants_symm
            and ut_wants_symm and symm_res * margin < anti_res):
        return "SYMM", symm_res, anti_res
    if (anti_res < clear and symm_res > clear and ur_wants_anti
            and ut_wants_anti and anti_res * margin < symm_res):
        return "ANTI", symm_res, anti_res
    return "AMBIGUOUS", symm_res, anti_res


def classify_fe_mode(r, theta_deg, ux, uy, phi_deg=180.0,
                     miss_frac=0.02):
    """Classify one FE eigenvector about the sector bisector.

    A handful of unpaired nodes (F10.4 rounding, crack nodes) does not
    void the class. If more than miss_frac of nodes have no partner the
    mesh is not a bisector pair and the call is AMBIGUOUS.
    """
    r = np.asarray(r, dtype=float)
    theta_deg = np.asarray(theta_deg, dtype=float)
    partner, missing = _partner_index(r, theta_deg, phi_deg=phi_deg)
    n = len(r)
    if n == 0 or (missing / float(n)) > miss_frac:
        return "AMBIGUOUS", float("nan"), float("nan"), missing
    ur, ut = polar(ux, uy, theta_deg)
    ratios = parity_ratios(ur, ut, partner)
    cls, symm_res, anti_res = classify_parity(*ratios)
    return cls, symm_res, anti_res, missing


def class_safe_nearest(target_hz, mode_hz, mode_class, want_class):
    """Nearest FE mode by relative frequency error, restricted to want_class.

    mode_hz and mode_class are dicts keyed by mode index. Returns
    (mode_index, rel_err) or (None, nan) if no mode of that class exists.
    """
    target_hz = float(target_hz)
    best = None
    best_rel = float("nan")
    for m, hz in mode_hz.items():
        if mode_class.get(m) != want_class:
            continue
        if hz <= 0:
            continue
        rel = abs(float(hz) - target_hz) / target_hz
        if best is None or rel < best_rel:
            best, best_rel = m, rel
    return best, best_rel


def unrestricted_nearest(target_hz, mode_hz):
    """Production's historical match: nearest by relative error, any class."""
    target_hz = float(target_hz)
    best = None
    best_rel = float("nan")
    for m, hz in mode_hz.items():
        if hz <= 0:
            continue
        rel = abs(float(hz) - target_hz) / target_hz
        if best is None or rel < best_rel:
            best, best_rel = m, rel
    return best, best_rel


def load_eigvec(path):
    """Load an ANSYS cylindrical dump: r, theta_deg, UX, UY. Skip header."""
    data = np.loadtxt(path, skiprows=1)
    if data.ndim != 2 or data.shape[1] < 4:
        raise RuntimeError(f"unexpected shape {data.shape} in {path}")
    return data[:, 0], data[:, 1], data[:, 2], data[:, 3]


def ensure_free_free_worker_bc():
    """Export R40_BC_KIND=free_free. Call before any ProcessPoolExecutor.

    Does not use setdefault: a leftover clamped_free from an earlier
    solver construction in this process would otherwise win.
    """
    os.environ["R40_BC_KIND"] = "free_free"
    return os.environ["R40_BC_KIND"]
