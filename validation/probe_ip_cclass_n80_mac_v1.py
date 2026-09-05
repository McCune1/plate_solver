# -*- coding: utf-8 -*-
"""
probe_ip_cclass_n80_mac_v1.py -- DIAGNOSTIC ONLY.

C-class extract-completeness (Paper 1 §6.6).

Rank 14 scored the 8 keys at 2Theta/pi = 1.25 and 1.50 against a 57-mode
PLANE183 extract (modes 4..60) and every tight FE-frequency match was
MAC-ART. The paper now qualifies that claim as extract-limited. This
probe repeats Rank 14 at ONE of those keys against a DEEPER dump
(modes 4..83 = 80 elastic modes, prefix *_n80) without retuning 0.744
and without touching the Rank 14 57-mode files.

Candidate loop is a ProcessPoolExecutor (C80_WORKERS, default
SLURM_CPUS_PER_TASK). plate_solver N_WORKERS is forced to 1 so inner
full_search does not nest pools.

Does NOT bump SOLVER_VERSION. Does NOT change a computed Omega.

Env:
  RATIO_TAG, ANGLE_TAG, GS_NU (default 0.35)
  MAX_CANDS     optional cap (unset = all)
  C80_WORKERS   candidate-parallelism (default SLURM_CPUS_PER_TASK or 1)

PRE-REGISTERED INTERPRETATION
-----------------------------
C80-0 GATE_FAIL
    Missing n80 dumps / manifest, SOLVER_VERSION != s10, ANGLE_TAG not
    a125/a150, mac() identity fail, or zero candidates reconstructed. VOID.
C80-D_NO_TIGHT
    Gates pass but no candidate is within 0.5% of an FE mode. Matching
    problem. Do not edit the paper.
C80-C_SURVIVES
    Every tight candidate is still MAC-ART at the 80-elastic-mode extract.
    The C-class claim is NOT extract-limited at this key.
C80-REAL_APPEARS
    At least one tight MAC-REAL appears that Rank 14 did not have. The
    abstract's 57-mode qualifier MUST stay. Do not retune 0.744.
C80-A_SPLIT
    Tight REAL and ART both present. Same as REAL_APPEARS for the paper:
    the 57-mode C-class sentence is extract-limited.

Honour the printed letter for THIS key only. All 8 C-class keys are
required before dropping the 57-mode qualifier.
"""
import importlib.util
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

sys.path.insert(0, os.environ.get("PKG_PATH", "."))
os.environ.setdefault("DPS", "40")
# Inner detector/search must not spawn its own pool; we parallelize
# candidates at this layer.
os.environ["N_WORKERS"] = "1"

import numpy as np

EXPECT_VER = os.environ.get("EXPECT_SOLVER_VERSION", "2026-07-10.s10")
PKG = os.environ.get("PKG_PATH", ".")
RATIO_TAG = os.environ.get("RATIO_TAG", "r150").strip().lower()
ANGLE_TAG = os.environ.get("ANGLE_TAG", "a125").strip().lower()
NU = float(os.environ.get("GS_NU", "0.35"))
RULE_DELTA = float(os.environ.get("RULE_DELTA", "0.001"))
MAX_CANDS = os.environ.get("MAX_CANDS", "").strip()
MAX_CANDS = int(MAX_CANDS) if MAX_CANDS else None
TIGHT = 0.005  # 0.5% relative FE error
C80_WORKERS = int(os.environ.get(
    "C80_WORKERS", os.environ.get("SLURM_CPUS_PER_TASK", "1")))

RATIOS = {"r150": 1.5, "r167": 1.66667, "r200": 2.0, "r250": 2.5}
ANGLES = {
    "a025": 0.25, "a050": 0.50, "a075": 0.75,
    "a100": 1.00, "a125": 1.25, "a150": 1.50,
}
if RATIO_TAG not in RATIOS or ANGLE_TAG not in ("a125", "a150"):
    print(f"FATAL: C-class probe takes a125/a150 only, got "
          f"{RATIO_TAG!r}/{ANGLE_TAG!r}")
    raise SystemExit(2)

R0_2B = RATIOS[RATIO_TAG]
TWO_T_PI = ANGLES[ANGLE_TAG]
PHI_DEG = TWO_T_PI * 180.0
THETA = TWO_T_PI * np.pi / 2.0
PREFIX = os.environ.get(
    "FB_PREFIX", f"ip_mac_gs_{RATIO_TAG}_{ANGLE_TAG}_n80")
MODE_LO = int(os.environ.get("FB_MODE_LO", "4"))
MODE_HI = int(os.environ.get("FB_MODE_HI", "83"))  # 80 elastic modes
N_FE = MODE_HI - MODE_LO + 1
TOL = 0.002

# Filled in the parent, inherited/initialized in workers.
FEB = None
MODE_HZ = None
MODE_CLS = None
_ADV = None
_BAR = None
_SOLVER = None
_KBAR = None
_WBAR = None


def hdr(s):
    print("\n" + "=" * 78 + f"\n  {s}\n" + "=" * 78, flush=True)


def find_file(*names):
    roots = [PKG, os.path.join(PKG, "Ansys", "NewAnsys"),
             ".", os.path.join("Ansys", "NewAnsys")]
    for n in names:
        for d in roots:
            p = os.path.join(d, n)
            if os.path.isfile(p):
                return os.path.abspath(p)
    return None


def load_adv():
    path = find_file("probe_ip_mac_adversarial_fullblock_v1.py")
    if path is None:
        print("FATAL: adversarial reconstruction probe not found")
        raise SystemExit(2)
    spec = importlib.util.spec_from_file_location("adv_mac_c80", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    mod.THETA = float(THETA)
    return mod


def load_ip_targets():
    cands = [
        os.environ.get("MANIFEST", ""),
        os.path.join(PKG, f"figures_ff_{RATIO_TAG}_s10v2",
                     "freefree_manifest.json"),
        os.path.join(PKG, "FutureWork", "geometry_sweep",
                     f"figures_ff_{RATIO_TAG}_s10v2",
                     "freefree_manifest.json"),
    ]
    path = next((p for p in cands if p and os.path.isfile(p)), None)
    if path is None:
        print("FATAL: s10v2 manifest not found. On Mill it is")
        print(f"  $PKG_PATH/figures_ff_{RATIO_TAG}_s10v2/freefree_manifest.json")
        raise SystemExit(2)
    print(f"  MANIFEST={path}", flush=True)
    manifest = json.load(open(path))
    out = []
    for key, entry in manifest.items():
        if key.startswith("_"):
            continue
        diag = entry.get("diag") or {}
        if diag.get("part") != 2:
            continue
        try:
            tpi = float(key.rsplit(" ", 1)[1])
        except (ValueError, IndexError):
            continue
        if abs(tpi - TWO_T_PI) > 1e-6:
            continue
        accepted = sorted(diag.get("accepted") or [],
                          key=lambda a: a["Omega"])
        n_dofs = diag.get("n_dofs", 20)
        xmax = diag.get("max_dim", 20.0)
        for a in accepted:
            out.append((key, float(a["Omega"]), int(n_dofs), float(xmax)))
    return out


def _init_worker(feb, mode_hz, mode_cls):
    """Rebuild a solver in each worker. Do not reuse the parent's solver."""
    global FEB, MODE_HZ, MODE_CLS, _ADV, _BAR, _SOLVER, _KBAR, _WBAR
    os.environ["R40_BC_KIND"] = "free_free"
    os.environ["N_WORKERS"] = "1"
    os.environ.setdefault("DPS", "40")
    from mpmath import mp
    mp.dps = int(os.environ.get("DPS", "40"))
    import ip_mac_bar as bar
    import plate_solver as ps
    from plate_solver.geometry import _kbar_wbar
    adv = load_adv()
    mat = ps.IsotropicMaterial(E=210e9, nu=NU, rho=7800.0)
    geom = ps.make_geometry(R0_2B, TWO_T_PI)
    solver = ps.InPlaneSolver(geom, mat, M=80, n_quad=30,
                              boundary=ps.FreeFreeIP())
    k_bar, w_bar = _kbar_wbar(geom, mat)
    FEB = feb
    MODE_HZ = mode_hz
    MODE_CLS = mode_cls
    _ADV = adv
    _BAR = bar
    _SOLVER = solver
    _KBAR = k_bar
    _WBAR = w_bar


def score_one(item):
    """Top-level so ProcessPoolExecutor can pickle it. Returns (i, log, row)."""
    i, ntot, key, Om, n_dofs, xmax = item
    t0 = time.time()
    lines = []

    def log(s):
        lines.append(s)

    try:
        from plate_solver.geometry import _omega_lit
        adv, bar, solver = _ADV, _BAR, _SOLVER
        adv.XMAX0 = xmax
        adv.N_DOFS = n_dofs
        geom = solver.geom
        mat = solver.mat
        f_hz, om_lit = _omega_lit(Om, geom, mat, _WBAR, _KBAR, 2)
        f_hz = float(f_hz)
        log("\n" + "=" * 78 + f"\n  [{i}/{ntot}]  {key}  Omega={Om:.6f}\n"
            + "=" * 78)
        log(f"  f_hz={f_hz:.4f}")
        brs, _x, st, cnt, _n = adv.get_filled_brs(solver, Om, n_dofs,
                                                  xmax=xmax)
        log(f"  fill state={st} cnt={cnt}/{n_dofs}")
        if st == "UNRELIABLE":
            return i, lines, dict(Om=Om, f_hz=f_hz, ok=False, why="UNRELIABLE")
        sh = adv.reconstruct_shared(solver, Om, brs, n_dofs=n_dofs)
        if not sh.get("ok"):
            return i, lines, dict(Om=Om, f_hz=f_hz, ok=False,
                                  why=sh.get("reason"))
        blocks = {}
        for w in ("q0", "q1"):
            r = adv.block_reconstruction(solver, sh, w)
            if r.get("ok"):
                blocks[w] = r
                log(f"  block {w}: rTyy={r['rTyy']:.6g}")
        if len(blocks) < 2:
            return i, lines, dict(Om=Om, f_hz=f_hz, ok=False, why="block fail")
        s0, s1, _c = adv.sigma_pair_at(solver, Om, n_dofs=n_dofs)
        m0, m1, _c0 = adv.sigma_pair_at(solver, Om - RULE_DELTA, n_dofs=n_dofs)
        p0, p1, _c1 = adv.sigma_pair_at(solver, Om + RULE_DELTA, n_dofs=n_dofs)
        if None in (s0, s1, m0, m1, p0, p1):
            return i, lines, dict(Om=Om, f_hz=f_hz, ok=False, why="sigma fail")
        pick, how = bar.rule_localmin_resid(
            s0, s1, m0, p0, m1, p1,
            blocks["q0"]["rTyy"], blocks["q1"]["rTyy"],
        )
        want_cls = bar.block_class(pick)
        log(f"  RULE -> {pick} [{how}] class {want_cls}")
        rp, tp, Ur, Uth = adv.build_grid_from_cfull(
            solver, sh, blocks[pick]["c_full"])
        ux_s, uy_s, bad = adv.solver_field_at_fe(rp, tp, Ur, Uth, FEB[MODE_LO])
        keep = ~bad
        bvec = np.concatenate([ux_s[keep], uy_s[keep]])
        Ablk = np.vstack([np.concatenate([FEB[m][2][keep], FEB[m][3][keep]])
                          for m in range(MODE_LO, MODE_HI + 1)])
        macs = adv.mac_matrix(Ablk, bvec)
        modes = list(range(MODE_LO, MODE_HI + 1))
        mac_by = {m: float(macs[j]) for j, m in enumerate(modes)}
        unres, rel_u = bar.unrestricted_nearest(f_hz, MODE_HZ)
        safe, rel_s = bar.class_safe_nearest(f_hz, MODE_HZ, MODE_CLS, want_cls)
        same = [m for m in modes if MODE_CLS.get(m) == want_cls]
        if not same:
            log(f"  no FE mode of class {want_cls} "
                f"(SYMM/ANTI/AMBIGUOUS split may be empty)")
            mac_same = (float("nan"), -1)
            call = "INVALID"
        else:
            mac_same = max((mac_by[m], m) for m in same)
            call = bar.classify_mac(mac_same[0])
        tight = (rel_u == rel_u and rel_u < TIGHT)
        log(f"  unres #{unres} rel={rel_u*100:.3f}% "
            f"MAC={mac_by.get(unres, float('nan')):.4f}  "
            f"safe #{safe} MAC={mac_by.get(safe, float('nan')):.4f}")
        log(f"  max same-class={mac_same[0]:.4f} at #{mac_same[1]}  "
            f"CALL {call}  tight={tight}")
        dt = time.time() - t0
        log(f"  [{dt:.0f}s]")
        row = dict(
            Om=Om, f_hz=f_hz, ok=True, pick=pick, want_cls=want_cls,
            unres=unres, rel_u=rel_u, mac_near=mac_by.get(unres, float("nan")),
            safe=safe, mac_best_same=mac_same[0], best_same=mac_same[1],
            call=call, tight=tight, rTyy=blocks[pick]["rTyy"], dt=dt,
        )
        return i, lines, row
    except Exception as exc:
        dt = time.time() - t0
        log(f"  WORKER_EXC {type(exc).__name__}: {exc}")
        log(f"  [{dt:.0f}s]")
        return i, lines, dict(Om=Om, f_hz=float("nan"), ok=False,
                              why=f"{type(exc).__name__}: {exc}")


def main():
    t_all = time.time()
    hdr(f"probe_ip_cclass_n80_mac_v1  {RATIO_TAG}/{ANGLE_TAG}  nu={NU}")
    print(f"  r0/2b={R0_2B}  2Theta/pi={TWO_T_PI}  PHI={PHI_DEG} deg",
          flush=True)

    import ip_mac_bar as bar
    bar.ensure_free_free_worker_bc()
    import plate_solver as ps
    print(f"  SOLVER_VERSION={ps.SOLVER_VERSION!r}  "
          f"R40_BC_KIND={os.environ.get('R40_BC_KIND')!r}", flush=True)
    print(f"  N_WORKERS(plate_solver)={os.environ.get('N_WORKERS')!r}  "
          f"C80_WORKERS={C80_WORKERS}", flush=True)
    if ps.SOLVER_VERSION != EXPECT_VER:
        print("READING: C80-0 GATE_FAIL  (SOLVER_VERSION)", flush=True)
        raise SystemExit(3)

    adv = load_adv()
    print(f"  reconstruction imported; adv.THETA set to {adv.THETA:.6f} "
          f"(expect {THETA:.6f})", flush=True)
    if abs(float(adv.THETA) - float(THETA)) > 1e-12:
        print("READING: C80-0 GATE_FAIL  (THETA)", flush=True)
        raise SystemExit(4)

    rng = np.random.default_rng(14)
    a = rng.normal(size=200)
    b = rng.normal(size=200)
    dmac = abs(bar.mac(a, b) - adv.mac(a, b))
    print(f"  mac() identity |d|={dmac:.2e}", flush=True)
    if dmac > 1e-15:
        print("READING: C80-0 GATE_FAIL  (mac)", flush=True)
        raise SystemExit(4)

    freq_path = find_file(f"{PREFIX}_freqs.txt")
    if freq_path is None:
        print(f"READING: C80-0 GATE_FAIL  (missing {PREFIX}_freqs.txt)",
              flush=True)
        raise SystemExit(2)
    FREQ = adv.load_freqs(freq_path)
    feb = {}
    for m in range(MODE_LO, MODE_HI + 1):
        p = find_file(f"{PREFIX}_m{m}_mesh96.txt")
        if p is None:
            print(f"READING: C80-0 GATE_FAIL  (missing mode {m})",
                  flush=True)
            raise SystemExit(2)
        feb[m] = adv.load_fe_eigvec(p)
    print(f"  dumps {len(feb)}/{N_FE}  n_nodes={len(feb[MODE_LO][0])}  "
          f"freqs from {freq_path}", flush=True)

    mode_hz = {m: float(FREQ[m]) for m in range(MODE_LO, MODE_HI + 1)
               if m in FREQ}
    mode_cls = {}
    for m in range(MODE_LO, MODE_HI + 1):
        r, th, ux, uy = feb[m]
        cls, sr, ar, miss = bar.classify_fe_mode(r, th, ux, uy,
                                                 phi_deg=PHI_DEG)
        mode_cls[m] = cls
        if m in (4, 5, 60) or miss:
            print(f"  FE#{m} {cls}  miss={miss}  "
                  f"f={mode_hz.get(m, float('nan')):.4f}", flush=True)
    n_s = sum(1 for c in mode_cls.values() if c == "SYMM")
    n_a = sum(1 for c in mode_cls.values() if c == "ANTI")
    n_u = sum(1 for c in mode_cls.values() if c == "AMBIGUOUS")
    print(f"  live parity SYMM={n_s} ANTI={n_a} AMBIGUOUS={n_u}", flush=True)
    if n_u:
        print("  WARNING: AMBIGUOUS FE modes present; class-safe nearest "
              "will skip them.", flush=True)

    U = np.vstack([np.concatenate([feb[m][2], feb[m][3]])
                   for m in range(MODE_LO, MODE_HI + 1)])
    G = U @ U.T
    norms = np.diag(G).copy()
    den = np.outer(norms, norms)
    MAC = np.where(den > 0, (G * G) / den, np.nan)
    off = MAC.copy()
    np.fill_diagonal(off, -1.0)
    imax = int(np.nanargmax(off))
    i0, j0 = divmod(imax, N_FE)
    max_off = float(off[i0, j0])
    print(f"  FE-vs-FE max off-diag {max_off:.6f} between "
          f"{i0+MODE_LO} and {j0+MODE_LO}  (reported, not gated -- new nu/geom)",
          flush=True)

    targets = load_ip_targets()
    print(f"  IP accepted at this key: {len(targets)}", flush=True)
    if MAX_CANDS is not None:
        targets = targets[:MAX_CANDS]
        print(f"  MAX_CANDS={MAX_CANDS} -> using {len(targets)}", flush=True)
    if not targets:
        print("READING: C80-0 GATE_FAIL  (no IP candidates)", flush=True)
        raise SystemExit(2)

    nw = max(1, min(C80_WORKERS, len(targets)))
    print(f"  scoring {len(targets)} candidates on {nw} processes", flush=True)

    items = [(i, len(targets), key, Om, n_dofs, xmax)
             for i, (key, Om, n_dofs, xmax) in enumerate(targets, 1)]
    by_i = {}
    n_done = 0
    with ProcessPoolExecutor(
            max_workers=nw,
            initializer=_init_worker,
            initargs=(feb, mode_hz, mode_cls)) as ex:
        futs = [ex.submit(score_one, it) for it in items]
        for fut in as_completed(futs):
            i, lines, row = fut.result()
            by_i[i] = (lines, row)
            n_done += 1
            dt = row.get("dt", float("nan"))
            print(f"  progress {n_done}/{len(targets)}  "
                  f"cand {i}  dt={dt:.0f}s  ok={row.get('ok')}",
                  flush=True)

    rows = []
    n_ok = 0
    for i in range(1, len(targets) + 1):
        lines, row = by_i[i]
        print("\n".join(lines), flush=True)
        rows.append(row)
        if row.get("ok"):
            n_ok += 1

    hdr("SUMMARY")
    print(f"  reconstructed {n_ok}/{len(targets)}", flush=True)
    print(f"  {'Om':>10s} {'Hz':>10s} {'blk':3s} {'near':>5s} {'rel%':>7s} "
          f"{'MAC@n':>8s} {'max-sc':>8s} {'at':>4s} tight call",
          flush=True)
    for r in rows:
        if not r.get("ok"):
            print(f"  {r['Om']:10.6f} FAILED {r.get('why')}", flush=True)
            continue
        print(f"  {r['Om']:10.6f} {r['f_hz']:10.3f} {r['pick']:3s} "
              f"{str(r['unres']):>5s} {r['rel_u']*100:7.3f} "
              f"{r['mac_near']:8.4f} {r['mac_best_same']:8.4f} "
              f"#{r['best_same']:<3d} {int(r['tight'])} {r['call']}",
              flush=True)

    hdr("READING (this key only -- honour the letter)")
    if n_ok == 0:
        print("READING: C80-0 GATE_FAIL  (nothing reconstructed)", flush=True)
        raise SystemExit(4)
    tight_rows = [r for r in rows if r.get("ok") and r.get("tight")]
    n_tr = sum(1 for r in tight_rows if r["call"] == "REAL-like")
    n_ta = sum(1 for r in tight_rows if r["call"] == "ARTIFACT-like")
    print(f"  tight (rel<0.5%) : {len(tight_rows)}  "
          f"MAC-REAL={n_tr}  MAC-ART={n_ta}", flush=True)
    if not tight_rows:
        letter = "C80-D_NO_TIGHT"
        note = "No candidate within 0.5% of an FE mode. Do not edit §6.4."
    elif n_tr >= 1 and n_ta >= 1:
        letter = "C80-A_SPLIT"
        note = "Tight REAL appeared at n80. Keep the 57-mode qualifier."
    elif n_tr == len(tight_rows):
        letter = "C80-REAL_APPEARS"
        note = "Every tight match is MAC-REAL at n80. Keep the 57-mode qualifier."
    else:
        letter = "C80-C_SURVIVES"
        note = ("Every tight match is still MAC-ART at 80 elastic modes. "
                "C-class claim is not extract-limited at this key.")
    print(f"\nREADING: {letter}", flush=True)
    print(f"  {note}", flush=True)
    print(f"  key={RATIO_TAG}/{ANGLE_TAG} nu={NU}  n_fe={N_FE}  "
          f"not an 8-key close.", flush=True)
    print(f"elapsed_s={time.time()-t_all:.1f}", flush=True)
    print(f"End: {time.strftime('%Y-%m-%d %H:%M:%S')}", flush=True)


if __name__ == "__main__":
    main()
