"""
probe_rect_ff_persist_golden_polish_2026-09-07.py

Purpose
-------
Decisive follow-up to the fine-grid recheck (job 2458006), which itself
followed up job 2457815's persist-threshold sensitivity sweep. Job
2457815 found three of Paper 2's 148 published production matches sit at
delta=0.01000 -- the maximum the frozen Screen-B window scan (0.002 step,
+/-0.02 window) can report before a match stops persisting:

    OOP l/b=1.5 SYM  Lambda*=5.380
    OOP l/b=2.5 ANTI Lambda*=6.260
    IP  l/b=3.0 ANTI Omega*=2.2400

Job 2458006 rescanned just these three at 4x finer resolution (0.0005
step) to check whether "delta=0.01000" was a genuine margin or a coarse-
grid artifact. The result was NOT simply reassuring: two of the three
got WORSE, not better, at finer resolution --

    OOP l/b=1.5 SYM:  fine delta = 0.01050  (> PERSIST_DL = 0.01)
    OOP l/b=2.5 ANTI: fine delta = 0.01100  (> PERSIST_DL = 0.01)
    IP  l/b=3.0 ANTI: fine delta = 0.00900  (< PERSIST_DL, and the fine
                       scan also revealed 3 nearby local minima where the
                       coarse scan saw only 1 -- a numerically crowded
                       region, not a single clean dip)

A uniform grid, no matter how fine, still only ever reports the sampled
point closest to sigma's true continuous minimum -- it cannot itself
settle whether the true minimum is inside or outside PERSIST_DL. This
probe settles that with a continuous minimizer instead of a finer grid:
golden-section search (the exact algorithm already used and cluster-run
elsewhere in this project, `probe_rect_ff_oop_mac_anti_remain_2026-09-06
.py`'s `golden()` function, copied verbatim for the OOP family; the IP
family gets an identical copy with only `sigma_at` swapped for
`sigma_at_ip`, per the project's exact-copy/minimal-marked-change
discipline -- not a rewrite, not a monkeypatch), bracketed tightly
around each candidate's own job-2458006 fine-grid dip location so the
search isolates that specific minimum rather than risking convergence
toward a neighboring one (relevant for the IP point, which has close
neighbors).

Does NOT retune PERSIST_DL, does NOT edit plate_solver internals, does
NOT touch a tabulated frequency, does NOT bump SOLVER_VERSION, does NOT
edit any paper file. This is read-only diagnostic work: the verdict
below is input to a human decision about Paper 2's Table tab:oop-high /
tab:ip92, not an automatic edit.

Fidelity gate: for each candidate, evaluate sigma at job 2458006's own
fine-grid nearest-dip location before polishing, and confirm it
reproduces that job's sigma value to relative 1e-3 -- if not, something
about the run environment has drifted and the polish should not be
trusted.

Verdict per candidate: "PERSISTS" if polished delta <= PERSIST_DL=0.01,
"FAILS_SCREEN_B_AT_CONTINUUM" if polished delta > 0.01. The second
outcome does NOT retroactively invalidate the published frequency (the
paper's stated Screen B procedure is the fixed 0.002-step scan, and by
that procedure's own definition every published match legitimately
persists) -- it means the published match's persistence is not robust
to the discretization step used to test it, which is a disclosable
methodological caveat for Paper 2's text.

Copy/paste on the cluster:
  cd /home/ghmkfh/PythonMill/Plate_Solver_Package
  sbatch submit_rect_ff_persist_golden_polish_2026-09-07.sh

Runtime: 3 candidates x 25 golden-section iterations x ~2 evals/iter
= ~150 sigma calls total. Expect well under 2 minutes given job
2458006's fine-grid per-point costs.
"""

from __future__ import annotations

import json
import os
import time

from mpmath import mp

import p5_rect_ff_lib as L
import plate_solver as ps

PKG = os.environ.get("PKG_PATH", ".")
EXPECT_VER = os.environ.get("EXPECT_SOLVER_VERSION", "2026-07-10.s10")
mp.dps = int(os.environ.get("DPS", "30"))

PERSIST_DL = L.PERSIST_DL
assert abs(PERSIST_DL - 0.01) < 1e-12, f"PERSIST_DL drifted: {PERSIST_DL}"

# Job 2458006's own fine-grid results: (family, lob, sym, L0, fine_x,
# fine_sigma, bracket_half_width). Bracket is tight (+/- 4 fine-grid
# steps = +/-0.002) so golden-section isolates the SAME dip found there,
# not a neighbor -- relevant for the IP point (job 2458006 found 3 nearby
# local minima there at fine resolution).
CANDIDATES = [
    dict(family="oop", lob=1.5, sym=True, L0=5.380,
         fine_x=5.3905, fine_sigma=3.1777158097929823e-12, half=0.002),
    dict(family="oop", lob=2.5, sym=False, L0=6.260,
         fine_x=6.249, fine_sigma=7.110450981279162e-20, half=0.002),
    dict(family="ip", lob=3.0, sym=False, L0=2.2400,
         fine_x=2.249, fine_sigma=3.950372707558929e-21, half=0.002),
]


def golden_oop(asm, lob, sym, lo, hi, n_real, n_cpair, im_cap, niter=25):
    """Verbatim copy of golden() from
    probe_rect_ff_oop_mac_anti_remain_2026-09-06.py -- unmodified."""
    phi = (5.0 ** 0.5 - 1.0) / 2.0
    a, b = lo, hi

    def f(x):
        s, n = L.sigma_at(asm, sym, lob, x, n_real, n_cpair, im_cap)
        return s, n

    c = b - phi * (b - a)
    d = a + phi * (b - a)
    fc, nc = f(c)
    fd, nd = f(d)
    hist = [(c, fc, nc), (d, fd, nd)]
    for _ in range(niter):
        if fc is None or fd is None:
            break
        if fc < fd:
            b, d, fd, nd = d, c, fc, nc
            c = b - phi * (b - a)
            fc, nc = f(c)
            hist.append((c, fc, nc))
        else:
            a, c, fc, nc = c, d, fd, nd
            d = a + phi * (b - a)
            fd, nd = f(d)
            hist.append((d, fd, nd))
    valid = [(x, s, n) for x, s, n in hist if s is not None]
    if not valid:
        return None
    x, s, n = min(valid, key=lambda t: t[1])
    return dict(Lam=x, sigma=float(s), n_used=n, lo=lo, hi=hi)


def golden_ip(asm, lob, sym, lo, hi, n_real, n_cpair, im_cap, niter=25):
    """Identical algorithm to golden_oop -- only marked change is
    sigma_at -> sigma_at_ip, per project exact-copy discipline."""
    phi = (5.0 ** 0.5 - 1.0) / 2.0
    a, b = lo, hi

    def f(x):
        s, n = L.sigma_at_ip(asm, sym, lob, x, n_real, n_cpair, im_cap)  # CHANGED
        return s, n

    c = b - phi * (b - a)
    d = a + phi * (b - a)
    fc, nc = f(c)
    fd, nd = f(d)
    hist = [(c, fc, nc), (d, fd, nd)]
    for _ in range(niter):
        if fc is None or fd is None:
            break
        if fc < fd:
            b, d, fd, nd = d, c, fc, nc
            c = b - phi * (b - a)
            fc, nc = f(c)
            hist.append((c, fc, nc))
        else:
            a, c, fc, nc = c, d, fd, nd
            d = a + phi * (b - a)
            fd, nd = f(d)
            hist.append((d, fd, nd))
    valid = [(x, s, n) for x, s, n in hist if s is not None]
    if not valid:
        return None
    x, s, n = min(valid, key=lambda t: t[1])
    return dict(Lam=x, sigma=float(s), n_used=n, lo=lo, hi=hi)


def main():
    print("probe_rect_ff_persist_golden_polish_2026-09-07")
    print("start", time.strftime("%Y-%m-%d %H:%M:%S"), f"dps={mp.dps}",
          flush=True)

    if not L.deploy_ok_ip(PKG, EXPECT_VER):
        raise SystemExit(2)
    print(f"SOLVER_VERSION={ps.SOLVER_VERSION}  PERSIST_DL={PERSIST_DL}",
          flush=True)

    asm_oop = L.make_ff()
    asm_ip = L.make_ff_ip()

    results = []
    for cand in CANDIDATES:
        family = cand["family"]
        lob, sym, L0 = cand["lob"], cand["sym"], cand["L0"]
        sl = "SYM" if sym else "ANTI"
        asm = asm_oop if family == "oop" else asm_ip
        print("=" * 72)
        print(f"{family} l/b={lob} {sl} L0={L0}  "
              f"(job 2458006 fine_x={cand['fine_x']}, "
              f"fine_sigma={cand['fine_sigma']:.6e})")
        print("=" * 72, flush=True)

        if family == "oop":
            n_real, n_cpair = L.n_basis(sym, persist=True)
            im_cap = 30.0
            s_check, _ = L.sigma_at(asm, sym, lob, cand["fine_x"],
                                     n_real, n_cpair, im_cap)
        else:
            n_real, n_cpair = L.n_basis_ip(persist=True)
            im_cap = 7.0
            s_check, _ = L.sigma_at_ip(asm, sym, lob, cand["fine_x"],
                                        n_real, n_cpair, im_cap)

        s_check_f = float(s_check) if s_check is not None else None
        gate_ok = (s_check_f is not None
                   and abs(s_check_f - cand["fine_sigma"])
                       <= 1e-3 * max(abs(cand["fine_sigma"]), 1e-30))
        print(f"  fidelity check: sigma({cand['fine_x']})={s_check_f:.6e} "
              f"vs job 2458006's {cand['fine_sigma']:.6e}  "
              f"{'OK' if gate_ok else 'FAIL'}", flush=True)
        if not gate_ok:
            print("  GATE FAIL -- skipping polish for this candidate, do "
                  "not trust any result for it.", flush=True)
            results.append(dict(family=family, lob=lob, sym=sl, L0=L0,
                                 gate_ok=False, polish=None,
                                 verdict="GATE_FAIL"))
            continue

        lo = cand["fine_x"] - cand["half"]
        hi = cand["fine_x"] + cand["half"]
        t0 = time.time()
        if family == "oop":
            rec = golden_oop(asm, lob, sym, lo, hi, n_real, n_cpair, im_cap)
        else:
            rec = golden_ip(asm, lob, sym, lo, hi, n_real, n_cpair, im_cap)
        dt = time.time() - t0

        if rec is None:
            print(f"  NO MIN in bracket [{lo},{hi}]  [{dt:.1f}s]", flush=True)
            results.append(dict(family=family, lob=lob, sym=sl, L0=L0,
                                 gate_ok=True, polish=None,
                                 verdict="NO_MIN_IN_BRACKET"))
            continue

        polished_delta = abs(rec["Lam"] - L0)
        verdict = ("PERSISTS" if polished_delta <= PERSIST_DL
                    else "FAILS_SCREEN_B_AT_CONTINUUM")
        print(f"  polished: x*={rec['Lam']:.6f}  sigma={rec['sigma']:.3e}  "
              f"n_used={rec['n_used']}  polished_delta={polished_delta:.6f}  "
              f"[{dt:.1f}s]", flush=True)
        print(f"  VERDICT: {verdict}", flush=True)
        results.append(dict(family=family, lob=lob, sym=sl, L0=L0,
                             gate_ok=True, polish=rec,
                             polished_delta=polished_delta,
                             verdict=verdict))

    print("\n" + "=" * 72)
    print("SUMMARY")
    print("=" * 72)
    for r in results:
        pd = r.get("polished_delta")
        pd_str = "" if pd is None else f" (polished_delta={pd:.6f})"
        print(f"  {r['family']:3s} lob={r['lob']:<4} {r['sym']:4s} "
              f"L0={r['L0']:<8} -> {r['verdict']}{pd_str}")

    n_fails = sum(1 for r in results
                  if r["verdict"] == "FAILS_SCREEN_B_AT_CONTINUUM")
    if n_fails:
        print(f"\n{n_fails} candidate(s) exceed PERSIST_DL=0.01 at the "
              f"continuum (golden-section) minimum, despite persisting "
              f"under the frozen 0.002-step production scan. This does "
              f"NOT invalidate the published frequency under Paper 2's "
              f"stated Screen B procedure, but it is a genuine "
              f"discretization-sensitivity caveat worth disclosing "
              f"before this check is cited as closing the rigor gap.")

    with open("rect_ff_persist_golden_polish_2026-09-07.json", "w") as f:
        json.dump(results, f, indent=1, default=str)
    print("\nWrote rect_ff_persist_golden_polish_2026-09-07.json", flush=True)
    print("SOLVER_VERSION unbumped. PERSIST_DL not retuned. "
          "Paper not edited.")


if __name__ == "__main__":
    main()
