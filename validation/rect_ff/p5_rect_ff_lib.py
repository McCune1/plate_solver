# Shared helpers for 2026-09-04 overnight P5 probes. Copy with the probes.
from __future__ import annotations

import json
import os

from mpmath import mp, mpf
from plate_solver import RectOOPAssembler, RectIPAssembler
from plate_solver.detectors import (
    rect_resolve_branches,
    rect_select_branches,
    rect_ip_resolve_branches,
    rect_ip_select_branches,
)

CHECKERBOARD = "p1t * Ut1 * q0t * Vt0 - p1b * Ut1 * q0b * Vt0"
OLD_SUM = "p1t * Ut1 * q0t * Vt0 + p1b * Ut1 * q0b * Vt0"
PERSIST_DL = 0.01
SIGMA_LIST = 0.3
KFAC = 24.6644
IP_OFAC = 804.481  # Hz per Omega_bar; job 2455569. Omega = f / IP_OFAC.

# Full non-rigid FE Lambda lists. Hz from first "mode 1." *VWRITE block,
# Lambda = f / 24.6644. mesh1 vs mesh2 <0.02% (lob100/150/250) and
# <0.001% (lob200/300, job 2454654). Modes with Lambda > 8 omitted
# (outside overnight scan windows).
FE_LISTS = {
    (1.0, True): [
        1.98249, 2.45439, 3.49250, 6.16026, 6.36553,
    ],
    (1.0, False): [
        1.35288, 3.49250, 6.16026, 6.93927, 7.74631,
    ],
    (1.5, True): [
        0.96347, 2.24373, 2.58809, 3.00395, 4.40408, 5.38152, 6.57855,
    ],
    (1.5, False): [
        0.89838, 2.07036, 3.83354, 6.04657, 6.52617, 6.70854,
    ],
    (2.0, True): [
        0.54338, 1.50735, 2.22558, 2.62175, 2.99812, 3.62227, 4.87111,
        5.08871, 6.69041, 7.44778,
    ],
    (2.0, False): [
        0.66858, 1.47069, 2.55180, 4.02613, 5.90138, 6.17866, 6.56213,
        7.41437,
    ],
    (2.5, True): [
        0.34766, 0.96559, 1.88370, 2.27059, 2.46702, 3.16376, 3.19522,
        4.13819, 4.75645, 5.34819, 6.61905, 6.80738,
    ],
    (2.5, False): [
        0.53126, 1.14003, 1.90392, 2.89206, 4.14962, 5.67275, 6.21965,
        6.36131,
    ],
    (3.0, True): [
        0.24127, 0.66984, 1.31785, 2.15417, 2.24841, 2.45935, 2.88704,
        3.32528, 3.62467, 4.53021, 4.62253, 5.59746, 6.14744, 6.82968,
        7.88785,
    ],
    (3.0, False): [
        0.44044, 0.93044, 1.51698, 2.24305, 3.14253, 4.23529, 5.51872,
        6.18191, 6.34339, 6.72208, 7.13275, 7.48050,
    ],
}

# Rectangular FFFF IP FE, job 2455569. Omega = f[Hz]/804.481 from mesh2.
# mesh1 vs mesh2 max conv 0.004%. Rigid dropped. Omega>2.50 omitted here
# (first overnight window); full lists in rect_ff_ip_*_lob*.txt.
IP_FE_LISTS = {
    (1.0, True): [
        1.32984, 1.41422, 1.60733, 1.85745, 2.00319,
    ],
    (1.0, False): [
        1.24858, 1.32984, 2.00319, 2.31523,
    ],
    (1.5, True): [
        1.04551, 1.41236, 1.42432, 1.62439, 1.69877, 2.10886, 2.22185, 2.28542,
    ],
    (1.5, False): [
        0.78782, 1.03335, 1.57124, 1.64765, 2.20253, 2.24065,
    ],
    (2.0, True): [
        0.79652, 1.40011, 1.41422, 1.44333, 1.65355, 1.77390, 1.81566, 2.00440,
        2.29856, 2.33104,
    ],
    (2.0, False): [
        0.52557, 0.87891, 1.27148, 1.28703, 1.73452, 1.84467, 2.13652,
    ],
    (2.5, True): [
        0.64049, 1.22293, 1.41225, 1.42658, 1.53449, 1.67003, 1.69350, 1.78140,
        2.04232, 2.07671, 2.13048, 2.38996, 2.40468,
    ],
    (2.5, False): [
        0.37585, 0.71776, 1.06883, 1.11248, 1.44246, 1.53132, 1.80308, 1.99255,
        2.14745, 2.34878,
    ],
    (3.0, True): [
        0.53502, 1.04554, 1.41422, 1.41791, 1.41978, 1.60066, 1.66346, 1.67627,
        1.80401, 1.83870, 2.00478, 2.16287, 2.21145, 2.27554, 2.44758, 2.46034,
    ],
    (3.0, False): [
        0.28169, 0.57789, 0.88739, 1.03153, 1.27024, 1.29465, 1.57876, 1.66509,
        1.91234, 1.99856, 2.25049, 2.27620, 2.49371,
    ],
}

# 2453559 PAPER_CANDIDATES (B-pass, not 0.41 leak).
PAPER_2453559 = [
    (1.0, True, 2.0800, 8.87e-02),
    (1.0, True, 2.2700, 2.95e-02),
    (1.0, False, 1.3600, 1.10e-02),
    (1.0, False, 1.6400, 6.15e-06),
    (1.0, False, 2.1200, 1.13e-04),
    (1.5, True, 0.9600, 3.49e-03),
    (1.5, True, 0.9900, 3.22e-03),
    (1.5, True, 1.9200, 2.03e-01),
    (1.5, True, 2.2500, 3.06e-03),
    (1.5, False, 0.9000, 3.24e-03),
    (1.5, False, 1.3800, 4.31e-05),
    (1.5, False, 1.6400, 8.87e-07),
    (1.5, False, 2.1200, 1.57e-05),
    (2.5, True, 0.3500, 1.04e-02),
    (2.5, True, 0.9800, 1.95e-03),
    (2.5, True, 1.8900, 1.46e-03),
    (2.5, True, 1.9200, 1.11e-02),
    (2.5, True, 2.2800, 2.70e-04),
    (2.5, False, 0.5400, 1.86e-04),
    (2.5, False, 0.6800, 1.20e-04),
    (2.5, False, 1.1400, 1.23e-03),
    (2.5, False, 1.6400, 5.21e-08),
    (2.5, False, 1.9200, 1.05e-07),
    (2.5, False, 2.1200, 8.97e-07),
]

# 2453876 PAPER_CANDIDATES (B-pass, not SYM 0.41 leak). Includes 3.0 ANTI
# 0.440 (0.46-family coincidence with the real 1st ANTI FE 0.44044) and
# the unmatched extras (ANTI 1.64, shallow SYM 1.92/1.93) so refine can
# confirm they stay unmatched.
PAPER_2453876 = [
    (2.0, True, 1.5100, 2.52e-04),
    (2.0, True, 1.5300, 1.23e-02),
    (2.0, True, 1.9300, 4.09e-02),
    (2.0, True, 2.2300, 1.48e-03),
    (2.0, True, 2.2700, 6.98e-03),
    (2.0, False, 0.6800, 9.03e-04),
    (2.0, False, 1.4800, 1.33e-03),
    (2.0, False, 1.6400, 2.47e-07),
    (2.0, False, 2.0800, 5.47e-07),
    (2.0, False, 2.1200, 3.07e-06),
    (3.0, True, 0.2400, 4.38e-03),
    (3.0, True, 1.3200, 1.41e-04),
    (3.0, True, 1.9200, 2.20e-02),
    (3.0, True, 2.1600, 1.94e-04),
    (3.0, True, 2.2600, 1.56e-04),
    (3.0, False, 0.4400, 1.47e-10),
    (3.0, False, 0.5400, 2.50e-05),
    (3.0, False, 0.9400, 9.54e-05),
    (3.0, False, 1.1600, 2.17e-04),
    (3.0, False, 1.5200, 7.30e-04),
    (3.0, False, 1.6400, 1.59e-08),
    (3.0, False, 1.9200, 5.26e-10),
    (3.0, False, 2.1200, 1.97e-07),
    (3.0, False, 2.2600, 7.95e-09),
]

# In-window FE with no 2453876 discovery dip (same class as 1.0 SYM 1.982).
MISS_WINDOWS = [
    (2.0, True, 0.54338),
    (3.0, True, 0.66984),
]

_RECT_TABLE3_SYM = {1.0: 0.3506, 1.5: 0.1551, 2.5: 0.0555}


def deploy_ok(pkg, expect_ver):
    import plate_solver as ps
    if ps.SOLVER_VERSION != expect_ver:
        print("FATAL SOLVER_VERSION", ps.SOLVER_VERSION)
        return False
    src = open(os.path.join(pkg, "plate_solver", "core_solvers.py"),
               encoding="utf-8").read()
    if CHECKERBOARD not in src or OLD_SUM in src:
        print("FATAL: not the checkerboard formula")
        return False
    print("G0 PASS", ps.SOLVER_VERSION, flush=True)
    return True


def deploy_ok_ip(pkg, expect_ver):
    if not deploy_ok(pkg, expect_ver):
        return False
    src = open(os.path.join(pkg, "plate_solver", "core_solvers.py"),
               encoding="utf-8").read()
    if "RectIPAssembler: unknown bc" not in src:
        print("FATAL: RectIPAssembler has no bc='free_free' switch")
        return False
    print("G0 IP bc=free_free present", flush=True)
    return True


def n_basis(sym, persist=False):
    if persist:
        return (6, 3) if sym else (5, 3)
    return (6, 0) if sym else (5, 1)


def n_basis_ip(persist=False):
    # Cantilever IP default is n_real=5 n_cpair=0. Persist trial uses
    # n_cpair=3 -- exploratory, NOT adopted Screen B.
    return (5, 3) if persist else (5, 0)


def frange(lo, hi, step):
    n = int(round((hi - lo) / step))
    return [round(lo + step * i, 6) for i in range(n + 1)]


def local_minima(rows):
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


def nearest_fe(lob, sym, lam):
    lst = FE_LISTS.get((float(lob), bool(sym))) or FE_LISTS.get((lob, sym))
    if not lst:
        return None, None
    Lfe = min(lst, key=lambda x: abs(x - lam))
    miss = 100.0 * abs(lam - Lfe) / Lfe
    return Lfe, miss


def nearest_ip_fe(lob, sym, om):
    lst = IP_FE_LISTS.get((float(lob), bool(sym))) or IP_FE_LISTS.get((lob, sym))
    if not lst:
        return None, None
    Ofe = min(lst, key=lambda x: abs(x - om))
    miss = 100.0 * abs(om - Ofe) / Ofe
    return Ofe, miss


def sigma_at(asm, sym, lob, Lam, n_real, n_cpair, im_cap):
    reps = rect_resolve_branches(asm.eng(sym), Lam, im_cap=im_cap)
    full = rect_select_branches(reps, n_real=n_real, n_cpair=n_cpair)
    if len(full) < 2:
        return None, 0
    K = asm.assemble(full, mpf(str(round(Lam, 6))), sym, lob)
    return asm.equil_sigma(K, full), len(full)


def persist_one(asm, lob, sym, L0, ckpt_dir):
    sl = "SYM" if sym else "ANTI"
    tag = f"persist_lob{lob}_{sl}_L{L0:.4f}".replace(".", "p")
    path = os.path.join(ckpt_dir, f"{tag}.json")
    if os.path.exists(path):
        with open(path) as f:
            return json.load(f)
    n_real, n_cpair = n_basis(sym, persist=True)
    lo = max(0.002, round(L0 - 0.02, 6))
    grid = frange(lo, round(L0 + 0.02, 6), 0.002)
    rows = []
    for Lam in grid:
        s, n = sigma_at(asm, sym, lob, Lam, n_real, n_cpair, 30.0)
        rows.append((Lam, s, n))
    dips = local_minima(rows)
    nearest = None
    persist = False
    if dips:
        nearest = min(dips, key=lambda t: abs(t[0] - L0))
        persist = abs(nearest[0] - L0) <= PERSIST_DL
    out = {
        "L0": L0, "persist": persist,
        "nearest": None if nearest is None else [nearest[0], float(nearest[1])],
        "dips": [[L, float(s), n] for L, s, n in dips],
    }
    os.makedirs(ckpt_dir, exist_ok=True)
    with open(path, "w") as f:
        json.dump(out, f)
    return out


def scan_window(asm, lob, sym, lo, hi, step, im_cap, n_real, n_cpair):
    grid = frange(lo, hi, step)
    rows = []
    for Lam in grid:
        s, n = sigma_at(asm, sym, lob, Lam, n_real, n_cpair, im_cap)
        rows.append((Lam, s, n))
    dips = local_minima(rows)
    return rows, dips


def sigma_at_ip(asm, sym, lob, Om, n_real, n_cpair, im_cap):
    reps = rect_ip_resolve_branches(asm.eng(sym), Om, im_cap=im_cap)
    full = rect_ip_select_branches(reps, n_real=n_real, n_cpair=n_cpair)
    if len(full) < 2:
        return None, 0
    K = asm.assemble(full, mpf(str(round(Om, 6))), sym, lob)
    return asm.equil_sigma(K, full), len(full)


def persist_one_ip(asm, lob, sym, O0, ckpt_dir):
    sl = "SYM" if sym else "ANTI"
    tag = f"persist_ip_lob{lob}_{sl}_O{O0:.4f}".replace(".", "p")
    path = os.path.join(ckpt_dir, f"{tag}.json")
    if os.path.exists(path):
        with open(path) as f:
            return json.load(f)
    n_real, n_cpair = n_basis_ip(persist=True)
    lo = max(0.002, round(O0 - 0.02, 6))
    grid = frange(lo, round(O0 + 0.02, 6), 0.002)
    rows = []
    for Om in grid:
        s, n = sigma_at_ip(asm, sym, lob, Om, n_real, n_cpair, 7.0)
        rows.append((Om, s, n))
    dips = local_minima(rows)
    nearest = None
    persist = False
    if dips:
        nearest = min(dips, key=lambda t: abs(t[0] - O0))
        persist = abs(nearest[0] - O0) <= PERSIST_DL
    out = {
        "O0": O0, "persist": persist,
        "nearest": None if nearest is None else [nearest[0], float(nearest[1])],
        "dips": [[O, float(s), n] for O, s, n in dips],
    }
    os.makedirs(ckpt_dir, exist_ok=True)
    with open(path, "w") as f:
        json.dump(out, f)
    return out


def scan_window_ip(asm, lob, sym, lo, hi, step, im_cap, n_real, n_cpair):
    grid = frange(lo, hi, step)
    rows = []
    for Om in grid:
        s, n = sigma_at_ip(asm, sym, lob, Om, n_real, n_cpair, im_cap)
        rows.append((Om, s, n))
    dips = local_minima(rows)
    return rows, dips


def make_ff():
    return RectOOPAssembler(nu_b=0.3, R=1.0, T=1.0, dps=mp.dps, bc="free_free")


def make_cf():
    return RectOOPAssembler(nu_b=0.3, R=1.0, T=1.0, dps=mp.dps, bc="clamped_free")


def make_ff_ip():
    return RectIPAssembler(nu_b=0.3, R=1.0, dps=mp.dps, bc="free_free")


def make_cf_ip():
    return RectIPAssembler(nu_b=0.3, R=1.0, dps=mp.dps, bc="clamped_free")
