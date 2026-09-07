#!/usr/bin/env python3
"""
Local strengthening diagnostics for Paper 2 OOP MAC mismatches.

Does NOT retune 0.744. Does NOT bump SOLVER_VERSION. Off-cluster.

PRE-REGISTERED (per SHAPE_MISMATCH / AMBIGUOUS / PAIRING_SWAP row)
----------------------------------------------------------------
A. Reconstruct at Lambda_FE instead of Lambda*. If MAC(FE partner)
   then >= 0.744: the published Lambda* was too coarse for the
   eigenvector; frequency match stands; report as CONFIRMED_AT_LFE.
B. Keep the production/persist K but take the 2nd and 3rd smallest
   right singular vectors (same equilibration as the first). If SV1
   MAC-matches a lower FE mode and SVk MAC-matches the frequency
   partner at >= 0.744: the kernel is multi-dimensional; report as
   CONFIRMED_SVk. Selecting SVk by MAC against FE is diagnostic, not
   a new published detector -- flag it as such.
C. n_cpair=6, im_cap=30 at Lambda*. If MAC(FE partner) >= 0.744:
   CONFIRMED_NCP6 -- shape was basis-limited, frequency already matched.
D. If none of A/B/C fire: remains a reconstruction miss; do not drop
   the frequency row.

USAGE
  python -u probe_rect_ff_oop_mac_strengthen_2026-09-05.py
  SMOKE=1 -> only 2.0 ANTI 1.482 and 2.5 SYM 3.190
"""
from __future__ import annotations

import importlib.util
import json
import os
import sys
import time

import numpy as np
from mpmath import mp, mpf, mpc, svd_r

sys.path.insert(0, os.environ.get("PKG_PATH", "."))
os.environ.setdefault("DPS", "30")
mp.dps = int(os.environ["DPS"])

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location(
    "macprobe", os.path.join(HERE, "probe_rect_ff_oop_mac_2026-09-05.py"))
M = importlib.util.module_from_spec(spec)
spec.loader.exec_module(M)

import p5_rect_ff_lib as L  # noqa: E402
import plate_solver as ps  # noqa: E402
from plate_solver.detectors import (  # noqa: E402
    rect_resolve_branches, rect_select_branches)

EXPECT_VER = os.environ.get("EXPECT_SOLVER_VERSION", "2026-07-10.s10")
OUT_JSON = os.path.join(HERE, "rect_ff_oop_mac_strengthen_2026-09-05.json")
MAC_BAR = 0.744


def load_mismatch_rows():
    src = os.path.join(HERE, "rect_ff_oop_mac_2026-09-05.json")
    payload = json.load(open(src, encoding="utf-8"))
    keep = []
    for r in payload["results"]:
        if r.get("table") not in (1, 2):
            continue
        if r.get("letter") in ("CONFIRMED", "CONFIRMED_PERSIST"):
            continue
        keep.append(r)
    if os.environ.get("SMOKE", "0") == "1":
        want = {(2.0, False, 1.482), (2.5, True, 3.19)}
        keep = [r for r in keep if (r["lob"], r["sym"], r["Lam"]) in want]
    return keep


def svd_k(asm, sym, lob, Lam, n_cpair, im_cap, k=3):
    n_real, _ = L.n_basis(sym, persist=False)
    if n_cpair >= 3:
        n_real, _ = L.n_basis(sym, persist=True)
    reps = rect_resolve_branches(asm.eng(sym), Lam, im_cap=im_cap)
    full = rect_select_branches(reps, n_real=n_real, n_cpair=n_cpair)
    if len(full) < 2:
        return None
    K = asm.assemble(full, mpf(str(round(float(Lam), 6))), sym, lob)
    Kr = K.copy()
    asm._realify_conjugate_pairs(Kr, full)
    n = Kr.rows
    Kr_real = mp.matrix(n, n)
    for i in range(n):
        for j in range(n):
            Kr_real[i, j] = mp.re(Kr[i, j])
    A = Kr_real.copy()
    dc = [mpf(1)] * n
    for _ in range(3):
        for i in range(n):
            s = max(abs(A[i, j]) for j in range(n))
            if s > 0:
                for j in range(n):
                    A[i, j] /= s
        for j in range(n):
            s = max(abs(A[i, j]) for i in range(n))
            if s > 0:
                for i in range(n):
                    A[i, j] /= s
                dc[j] /= s
    _, S, V = svd_r(A)
    sigmas = [float(S[i]) for i in range(S.rows)]
    smax = max(sigmas) if sigmas else 1.0
    vecs = []
    for kk in range(1, min(k, n) + 1):
        w = [V[n - kk, j] for j in range(n)]
        c = [dc[j] * w[j] for j in range(n)]
        nrm = max(abs(x) for x in c)
        if nrm == 0:
            continue
        Y = np.array([float(x / nrm) for x in c])
        X = M.invert_realify(Y, full)
        phin = float(mp.pi / 2 if sym else 0)
        Lf = float(Lam)
        eng = asm.eng(sym)
        branches = []
        for xi in full:
            e1, e2 = eng._etas(complex(xi), Lf)
            H1, H2 = eng._amp_ratio(complex(xi), Lf)
            branches.append((complex(xi), complex(e1), complex(e2),
                             complex(H1), complex(H2)))
        rec = dict(X=np.array([complex(x) for x in X], dtype=np.complex128),
                   branches=branches, phin=phin)
        vecs.append(dict(
            k=kk,
            sigma_ratio=sigmas[n - kk] / smax if smax else float("nan"),
            rec=rec, n_used=len(full),
        ))
    return dict(full=full, n_used=len(full), n_cpair=n_cpair, vecs=vecs)


def score_rec(rec, blk, Lfe):
    x1 = blk["x"] * M.SCALE
    x2 = blk["y"] * M.SCALE
    w_real, im_ratio = M.phase_align_real(M.eval_w(rec, x1, x2))
    partner = min(blk["modes"], key=lambda m: abs(m["lam"] - Lfe))
    mac_fe = M.mac(partner["uz"], w_real)
    macs = []
    for m in blk["modes"]:
        if m["lam"] < 0.2:
            continue
        macs.append((M.mac(m["uz"], w_real), m["midx"], m["lam"]))
    macs.sort(key=lambda t: -t[0])
    return dict(
        mac_fe=float(mac_fe), fe_midx=partner["midx"], fe_lam=partner["lam"],
        mac_best=float(macs[0][0]) if macs else float("nan"),
        best_midx=macs[0][1] if macs else None,
        best_lam=float(macs[0][2]) if macs else float("nan"),
        mac_2nd=float(macs[1][0]) if len(macs) > 1 else float("nan"),
        im_ratio=im_ratio,
    )


def letter_of(mac_fe, mac_best, best_midx, fe_midx):
    if not np.isfinite(mac_fe):
        return "D_GATE"
    if mac_fe >= MAC_BAR:
        return "CONFIRMED"
    if best_midx != fe_midx and np.isfinite(mac_best) and mac_best >= MAC_BAR:
        return "PAIRING_SWAP"
    if mac_fe >= 0.3:
        return "AMBIGUOUS"
    return "SHAPE_MISMATCH"


def main():
    t0 = time.time()
    print("=" * 78)
    print("  probe_rect_ff_oop_mac_strengthen -- start "
          + time.strftime("%Y-%m-%d %H:%M:%S"))
    print("=" * 78, flush=True)
    if not L.deploy_ok(os.environ.get("PKG_PATH", "."), EXPECT_VER):
        raise SystemExit(2)
    print(f"preflight OK {ps.SOLVER_VERSION}", flush=True)
    rows = load_mismatch_rows()
    print(f"mismatch rows: {len(rows)}", flush=True)
    asm = L.make_ff()
    cache = {}
    out = []
    for i, r in enumerate(rows, 1):
        sl = "SYM" if r["sym"] else "ANTI"
        print(f"\n[{i}/{len(rows)}] l/b={r['lob']:.1f} {sl} L*={r['Lam']} "
              f"Lfe={r['Lfe']} was {r['letter']}", flush=True)
        key = (r["lob"], r["sym"])
        if key not in cache:
            cache[key] = M.load_fe_block(r["lob"], r["sym"])
        blk = cache[key]
        entry = dict(lob=r["lob"], sym=r["sym"], Lam=r["Lam"], Lfe=r["Lfe"],
                     table=r["table"], was=r["letter"])
        trials = [
            ("A_Lfe_prod", r["Lfe"], 1 if not r["sym"] else 0, 7.0),
            ("A_Lfe_persist", r["Lfe"], 3, 30.0),
            ("C_ncp6_Lstar", r["Lam"], 6, 30.0),
            ("C_ncp6_Lfe", r["Lfe"], 6, 30.0),
        ]
        best_hit = None
        for name, Lam, ncp, imc in trials:
            t1 = time.time()
            pack = svd_k(asm, r["sym"], r["lob"], Lam, ncp, imc, k=3)
            if pack is None:
                print(f"  {name}: no branches", flush=True)
                continue
            sv_scores = []
            for v in pack["vecs"]:
                sc = score_rec(v["rec"], blk, r["Lfe"])
                let = letter_of(sc["mac_fe"], sc["mac_best"],
                                sc["best_midx"], sc["fe_midx"])
                sv_scores.append(dict(
                    sv=v["k"], sigma_ratio=v["sigma_ratio"],
                    n_used=v["n_used"], letter=let, **sc,
                ))
                tag = f"SV{v['k']}"
                print(
                    f"  {name} {tag} n={v['n_used']} sig={v['sigma_ratio']:.2e} "
                    f"MAC(FE#{sc['fe_midx']})={sc['mac_fe']:.4f} "
                    f"best=#{sc['best_midx']} {sc['mac_best']:.4f} {let} "
                    f"[{time.time()-t1:.1f}s]",
                    flush=True,
                )
            entry[name] = sv_scores
            sv1 = sv_scores[0] if sv_scores else None
            if sv1 and sv1["letter"] == "CONFIRMED" and best_hit is None:
                best_hit = (name, "SV1", sv1["mac_fe"])
            elif sv_scores:
                for sc in sv_scores[1:]:
                    if sc["letter"] == "CONFIRMED" and best_hit is None:
                        # SV2/3 recovery is diagnostic (FE-informed pick)
                        best_hit = (name, f"SV{sc['sv']}_diag", sc["mac_fe"])
        entry["best_hit"] = best_hit
        out.append(entry)

    n_a = n_b = n_c = n_none = 0
    for e in out:
        h = e.get("best_hit")
        if not h:
            n_none += 1
            continue
        if h[0].startswith("A_") and h[1] == "SV1":
            n_a += 1
        elif h[1].startswith("SV") and h[1] != "SV1":
            n_b += 1
        elif h[0].startswith("C_") and h[1] == "SV1":
            n_c += 1
        else:
            n_none += 1
    payload = dict(
        solver_version=ps.SOLVER_VERSION, n=len(out), results=out,
        n_A_Lfe=n_a, n_B_SVk=n_b, n_C_ncp6=n_c, n_none=n_none,
        elapsed_s=time.time() - t0,
    )
    with open(OUT_JSON, "w", encoding="utf-8", newline="\n") as f:
        json.dump(payload, f, indent=2, default=str)
    print("\n" + "=" * 78)
    print(f"  A (Lfe SV1 CONFIRMED): {n_a}")
    print(f"  B (SVk>=2 CONFIRMED, diagnostic): {n_b}")
    print(f"  C (n_cpair=6 SV1 CONFIRMED): {n_c}")
    print(f"  none: {n_none}")
    print(f"  wrote {OUT_JSON}  elapsed {time.time()-t0:.1f}s")
    print("=" * 78, flush=True)


if __name__ == "__main__":
    main()
