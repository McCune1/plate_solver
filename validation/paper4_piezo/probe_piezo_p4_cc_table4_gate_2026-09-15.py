"""
probe_piezo_p4_cc_table4_gate_2026-09-15.py

Cluster validation gate for the FIXED C-C piezo-coupled annular ring
forward model (see LESSONS_LEARNED.md Sec 18.147 and the sibling script
probe_piezo_p4_cc_forward_model_2026-09-15.py for the full derivation and
the h1**2->h1**4 bug fix this gate exists to confirm at production
precision). Sandbox spot-checks (dps=200-300, extraprec=800-1500) already
matched 6/6 hand-picked points to <0.05%; this gate reproduces ALL 27
non-trivial CPT/r0-h=60 rows of Duan2005 Table 4 (p=0,1,2 circumferential
x n=0,1,2 radial x h1/2h in {1/12,1/8,1/5}) at higher/production precision
as the real validation pass before this model is trusted for the paper.

PRE-REGISTERED PASS/FAIL (decide before reading results, per project
discipline): each of the 27 points must land within 0.15% of its Table 4
target (a generous margin above the <0.05% already seen in the sandbox at
lower dps -- this just guards against a real regression, not a tight
research-grade bound) AND each bisection must see a genuine sign flip
across its bracket (no silent same-sign-both-ends garbage). ALL 27 must
pass for PASS_ALL; any single miss is FAIL_ALL and must be investigated,
not waved through -- per this project's "never accept close enough"
discipline, a single failing point at this stage would mean the sandbox
match was not representative (precision-dependent bias), not that 26/27
is "good enough".

h1=0 (bare host) rows are NOT included here -- they don't exercise the
piezo cubic at all (a=0 identically) and are already independently
validated by probe_piezo_p4_elastic_cc_baseline_2026-09-15.py.

Geometry/material: same as the two sibling scripts (ri=0.1, ro=0.6,
h=0.01 = Duan2005's own half-thickness symbol, PZT4/steel Table 1
constants). Table 4 targets and their brackets are transcribed directly
from the OCR'd table (LESSONS_LEARNED Sec 18.147's re-derivation session);
brackets are +-2% of each target, generous enough for a clean single sign
flip (modes are well-separated in frequency here) but tight enough that a
bracket miss itself is informative.
"""
import json
import os
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

import mpmath as mp

MP_DPS = int(os.environ.get('PIEZO_DPS', '100'))
EXTRAPREC = int(os.environ.get('PIEZO_EXTRAPREC', '300'))
BISECT_ITERS = int(os.environ.get('PIEZO_BISECT_ITERS', '45'))
PASS_TOL_PCT = 0.15

mp.mp.dps = MP_DPS
pi = mp.pi

C11E = mp.mpf('132e9'); C13E = mp.mpf('73e9'); C33E = mp.mpf('115e9')
e31 = mp.mpf('4.1'); e33 = mp.mpf('14.1')
X11 = mp.mpf('7.124e-9'); X33 = mp.mpf('5.841e-9')
rho_pzt = mp.mpf('7500')
E_steel = mp.mpf('200e9'); nu_steel = mp.mpf('0.3'); rho_steel = mp.mpf('7800')

ri = mp.mpf('0.1'); ro = mp.mpf('0.6'); h = mp.mpf('0.01')

c11_bar = C11E - C13E**2 / C33E
e31_bar = e31 - (C13E / C33E) * e33
Xi33_bar = X33 + e33**2 / C33E
Xi11_bar = X11
d1 = mp.mpf(2) / 3 * E_steel * h**3 / (1 - nu_steel**2)


def cubic_coeffs(omega, h1, d2):
    """FIXED per LESSONS_LEARNED Sec 18.147: A2*w2 terms use h1**4."""
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
    roots = mp.polyroots([a, b, c, d], maxsteps=300, extraprec=EXTRAPREC)
    out = []
    for chi in roots:
        lam = 2*pi**2*Xi33_bar*chi / (h1**2*(2*Xi11_bar*chi + e31_bar*pi))
        out.append((chi, lam))
    return out


def radial_quad(n, r, lam):
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


def bisect(lo, hi, n, h1, iters):
    flo = cc_det(mp.mpf(lo), n, h1)
    fhi = cc_det(mp.mpf(hi), n, h1)
    sign_flip = (flo.real > 0) != (fhi.real > 0)
    for _ in range(iters):
        mid = (lo + hi) / 2
        fm = cc_det(mp.mpf(mid), n, h1)
        if (fm.real > 0) == (flo.real > 0):
            lo, flo = mid, fm
        else:
            hi = mid
    return (lo + hi) / 2, sign_flip


# Table 4 (Duan2005 p.137), CPT-based model, r0/h=60. h1_ratio keys are
# h1/2h. Each entry: (p_circumferential, n_radial_label, h1_ratio, target,
# bracket_half_width_fraction).
TABLE4_ROWS = []
_targets_h1_0 = {  # for reference/bracket-centering only, not solved here
    (0, 0): 2718, (0, 1): 7520, (0, 2): 14783,
    (1, 0): 2851, (1, 1): 7755, (1, 2): 15075,
    (2, 0): 3385, (2, 1): 8538, (2, 2): 16002,
}
_targets = {
    ('1/12', (0, 0)): 2792,  ('1/12', (0, 1)): 7723,  ('1/12', (0, 2)): 15182,
    ('1/12', (1, 0)): 2928,  ('1/12', (1, 1)): 7965,  ('1/12', (1, 2)): 15482,
    ('1/12', (2, 0)): 3477,  ('1/12', (2, 1)): 8768,  ('1/12', (2, 2)): 16434,
    ('1/8',  (0, 0)): 2853,  ('1/8',  (0, 1)): 7891,  ('1/8',  (0, 2)): 15512,
    ('1/8',  (1, 0)): 2992,  ('1/8',  (1, 1)): 8138,  ('1/8',  (1, 2)): 15818,
    ('1/8',  (2, 0)): 3553,  ('1/8',  (2, 1)): 8959,  ('1/8',  (2, 2)): 16791,
    ('1/5',  (0, 0)): 2989,  ('1/5',  (0, 1)): 8268,  ('1/5',  (0, 2)): 16253,
    ('1/5',  (1, 0)): 3135,  ('1/5',  (1, 1)): 8527,  ('1/5',  (1, 2)): 16574,
    ('1/5',  (2, 0)): 3723,  ('1/5',  (2, 1)): 9387,  ('1/5',  (2, 2)): 17594,
}
_h1_ratio_val = {'1/12': mp.mpf(1)/12, '1/8': mp.mpf(1)/8, '1/5': mp.mpf(1)/5}

for (ratio_key, (p, n_radial)), target in sorted(_targets.items()):
    TABLE4_ROWS.append((p, n_radial, ratio_key, target))


def _worker(task):
    p, n_radial, ratio_key, target = task
    t0 = time.time()
    h1 = _h1_ratio_val[ratio_key] * 2 * h
    lo = target * 0.98
    hi = target * 1.02
    root, sign_flip = bisect(lo, hi, p, h1, BISECT_ITERS)
    root_f = float(root)
    rel_err_pct = (root_f - target) / target * 100.0
    passed = sign_flip and abs(rel_err_pct) <= PASS_TOL_PCT
    return {
        'p': p, 'n_radial': n_radial, 'h1_ratio': ratio_key,
        'target': target, 'computed': root_f, 'rel_err_pct': rel_err_pct,
        'sign_flip_seen': sign_flip, 'pass': passed,
        'elapsed_s': time.time() - t0,
    }


def main():
    n_workers = int(os.environ.get('N_WORKERS', '4'))
    print(f"job start: dps={MP_DPS} extraprec={EXTRAPREC} "
          f"bisect_iters={BISECT_ITERS} n_workers={n_workers} "
          f"n_points={len(TABLE4_ROWS)}", flush=True)

    results = []
    with ProcessPoolExecutor(max_workers=n_workers) as ex:
        futs = {ex.submit(_worker, task): task for task in TABLE4_ROWS}
        for fut in as_completed(futs):
            r = fut.result()
            results.append(r)
            print(f"  p={r['p']} n={r['n_radial']} h1/2h={r['h1_ratio']:>4s} "
                  f"target={r['target']:>7.1f} computed={r['computed']:>10.3f} "
                  f"rel_err={r['rel_err_pct']:+.4f}% "
                  f"sign_flip={r['sign_flip_seen']} "
                  f"pass={r['pass']} ({r['elapsed_s']:.1f}s)", flush=True)

    results.sort(key=lambda r: (r['p'], r['n_radial'], r['h1_ratio']))
    all_pass = all(r['pass'] for r in results)
    worst = max(results, key=lambda r: abs(r['rel_err_pct']))

    out = {
        'dps': MP_DPS, 'extraprec': EXTRAPREC, 'bisect_iters': BISECT_ITERS,
        'pass_tol_pct': PASS_TOL_PCT, 'n_points': len(results),
        'n_pass': sum(1 for r in results if r['pass']),
        'worst_rel_err_pct': worst['rel_err_pct'],
        'worst_point': {k: worst[k] for k in ('p', 'n_radial', 'h1_ratio')},
        'results': results,
    }
    with open('piezo_p4_cc_table4_gate_results.json', 'w') as f:
        json.dump(out, f, indent=2)

    print("SENTINEL " + ("PASS_ALL" if all_pass else "FAIL_ALL") +
          f" {out['n_pass']}/{out['n_points']} worst_rel_err_pct="
          f"{worst['rel_err_pct']:+.4f}", flush=True)


if __name__ == "__main__":
    main()
