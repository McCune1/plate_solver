#!/usr/bin/env python3
"""Score round 2 of the Paper 4 3-D thin-limit check. Gates fixed before the run.

Round 1 (job 2530553, score_p4_3d_thinlimit_2026-09-23.txt) was
pre-registered UNRESOLVED: the coarse -> fine change at t = 5 mm was
0.066-0.101% at n = 1, 3, 4, above its 0.05% bar, so the thin-limit
intercept c was not scored (it was +0.004 to +0.055%, positive).

Round 2 adds the t = 10 coarse deck and an xfine mesh (96 x 288 x 8, one
more uniform halving) at t = 20, 10, 5. Every thickness then has a ratio-2
sequence coarse / fine / xfine, and the mesh error can be removed by
Richardson extrapolation instead of being bounded by a bar.

Inputs: the five round-1 decks and the four round-2 decks, all
ansys_p4_3d_bare_t{t}_{mesh}_2026-09-23*, plus
targets_p4_3d_thinlimit_2026-09-23.json (unchanged, written before round 1).
Loading and mode identification reuse score_p4_3d_thinlimit_2026-09-23.py
and score_p4_3d_nge1_2026-09-23.py unchanged.

Gated set: n = 0..4, first radial root, t = 20, 10, 5 mm.

  R0  All nine transcripts clean; exact element / symmetry counts; three
      rigid modes; every gated (n, 1) identified at purity >= 0.90.
  R1  Monotone ratio-2 convergence at every (t, n): d1 = f_c - f_f and
      d2 = f_f - f_x have the same sign, |d2| < |d1|, and the observed
      order p = log2(d1 / d2) lies in [1, 8].
  R2  |f_x / f_f - 1| <= 0.02% at every (t, n).
  Primary intercept, chosen by rule, not by result:
      if R1 passes at all 15 points: c_inf, from the Richardson limits
          f_inf = f_x + (f_x - f_f) / (2^p - 1);
      else if R2 passes at all 15 points: c_x, from the xfine values;
      else: UNRESOLVED.
      c = (8/3) delta(5) - 2 delta(10) + (1/3) delta(20),
      delta(t) = f / f_Kirchhoff(t) - 1.
  R3  (decisive) |c_primary| <= 0.03% at every gated n.

Reading, fixed now:
  PASS: the n >= 1 thin-plate operator, with the effective-shear row, is
      confirmed in the Kirchhoff limit to 0.03% at n = 0..4. The job 2530423
      G3 offset is thickness physics.
  FAIL with c well below -0.03% (the defect hypothesis predicted -0.2 to
      -0.46%): an n-dependent operator error. Withdraw the Paper 4 sentence.
  FAIL with small positive c: FE discretization not removed by Richardson.
      Report as a FAIL of this bar; do not re-bar.
Always reported: c_f (round-1 fine), c_x, c_inf, p per point, and the
Richardson-limit offsets delta_inf(t).

Usage (from Ansys/NewAnsys/):  python3 score_p4_3d_thinlimit_r2_2026-09-23.py [dir]
"""
import importlib.util
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DATE = "2026-09-23"
TARGETS = os.path.join(HERE, "targets_p4_3d_thinlimit_%s.json" % DATE)
THICK = (20, 10, 5)
MESHES = ("coarse", "fine", "xfine")
GATED_N = (0, 1, 2, 3, 4)
P_RANGE = (1.0, 8.0)
XF_BAR = 0.0002
C_BAR = 0.0003
DEFECT_PRED = "-0.2% to -0.46%"


def _r1():
    p = os.path.join(HERE, "score_p4_3d_thinlimit_%s.py" % DATE)
    spec = importlib.util.spec_from_file_location("thin_r1", p)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


R = _r1()
R.COUNTS["xfine"] = dict(elements=221184, symmetry_nodes=5026)


def intercept(d20, d10, d5):
    return (8.0 / 3.0) * d5 - 2.0 * d10 + (1.0 / 3.0) * d20


def score(directory, targets_path=TARGETS):
    with open(targets_path) as fh:
        T = json.load(fh)
    D = {(t, m): R.load(directory, t, m) for t in THICK for m in MESHES}
    L = []
    a = L.append
    a("Paper 4 3-D thin-limit ROUND 2 score  (%s)" % directory)
    r0 = True
    missing = False
    for t in THICK:
        for m in MESHES:
            d = D[(t, m)]
            a("t=%-2d %-6s raw Hz: %s" % (t, m, ", ".join(
                "%.10g" % x for x in d["freqs"][:12]) or "(none)"))
            for r in d["reasons"]:
                a("    R0: %s" % r)
            r0 &= not d["reasons"]
            missing |= bool(d.get("not_run"))
    a("")
    if missing:
        a("R0: NOT_RUN. RESULT: NOT_RUN")
        return "\n".join(L), "NOT_RUN"
    a("R0: %s" % ("PASS" if r0 else "FAIL"))

    def f(t, m, n):
        rec = D[(t, m)]["table"].get((n, 1))
        return None if rec is None else rec["f"]

    def kir(t, n):
        return T["t%d" % t]["%d_1" % n]["el_hz"]

    a("")
    a("Convergence per (t, n): f_c, f_f, f_x, observed order p, xfine move, Richardson limit")
    r1 = True
    r2 = True
    finf = {}
    for t in THICK:
        for n in GATED_N:
            fc, ff, fx = f(t, "coarse", n), f(t, "fine", n), f(t, "xfine", n)
            if None in (fc, ff, fx):
                a("  t=%d n=%d: missing" % (t, n))
                r1 = r2 = False
                continue
            d1, d2 = fc - ff, ff - fx
            ok1 = d1 != 0 and d2 != 0 and (d1 > 0) == (d2 > 0) and abs(d2) < abs(d1)
            p = math.log(d1 / d2, 2) if ok1 else float("nan")
            ok1 = ok1 and P_RANGE[0] <= p <= P_RANGE[1]
            mv = fx / ff - 1.0
            ok2 = abs(mv) <= XF_BAR
            r1 &= ok1
            r2 &= ok2
            fi = fx + (fx - ff) / (2.0 ** p - 1.0) if ok1 else None
            finf[(t, n)] = fi
            a("  t=%-2d n=%d: %.8f %.8f %.8f  p=%5.2f %s  xfine move %+.4f%% %s  f_inf=%s"
              % (t, n, fc, ff, fx, p, "ok" if ok1 else "R1-FAIL", 100 * mv,
                 "ok" if ok2 else "R2-FAIL",
                 "%.8f" % fi if fi is not None else "n/a"))
    a("R1 (monotone ratio-2 convergence, 1 <= p <= 8): %s" % ("PASS" if r1 else "FAIL"))
    a("R2 (|f_x/f_f - 1| <= %.2f%%): %s" % (100 * XF_BAR, "PASS" if r2 else "FAIL"))

    if r1:
        primary = "c_inf (Richardson)"
    elif r2:
        primary = "c_x (xfine)"
    else:
        primary = None
    a("Primary intercept by rule: %s" % (primary or "none -> UNRESOLVED"))

    a("")
    a("Offsets delta(t) = f/f_K - 1 and thin-limit intercepts (defect hypothesis predicted c ~ %s):"
      % DEFECT_PRED)
    r3 = True
    for n in GATED_N:
        row = {}
        for lab, get in (("f", lambda t: f(t, "fine", n)),
                         ("x", lambda t: f(t, "xfine", n)),
                         ("inf", lambda t: finf.get((t, n)))):
            vals = [get(t) for t in THICK]
            if None in vals:
                row[lab] = None
                continue
            ds = [v / kir(t, n) - 1.0 for v, t in zip(vals, THICK)]
            row[lab] = (ds, intercept(*ds))
        cp = None
        if primary and primary.startswith("c_inf") and row["inf"]:
            cp = row["inf"][1]
        elif primary and primary.startswith("c_x") and row["x"]:
            cp = row["x"][1]
        ok = cp is not None and abs(cp) <= C_BAR
        if primary:
            r3 &= ok
        parts = []
        for lab in ("f", "x", "inf"):
            if row[lab] is None:
                parts.append("c_%s=n/a" % lab)
            else:
                parts.append("c_%s=%+.4f%%" % (lab, 100 * row[lab][1]))
        dinf = row["inf"][0] if row["inf"] else None
        a("  n=%d: %s  %s  delta_inf(20/10/5)=%s"
          % (n, "  ".join(parts),
             ("R3 " + ("PASS" if ok else "FAIL")) if primary else "R3 n/a",
             "/".join("%+.4f%%" % (100 * x) for x in dinf) if dinf else "n/a"))
    if not primary:
        r3s = "UNRESOLVED"
    else:
        r3s = "PASS" if r3 else "FAIL"
    a("R3 (|c_primary| <= %.2f%%): %s" % (100 * C_BAR, r3s))

    if not r0:
        res = "FAIL"
    elif r3s == "UNRESOLVED":
        res = "UNRESOLVED"
    else:
        res = r3s
    a("")
    a("RESULT: %s  (R0 %s, R1 %s, R2 %s, R3 %s; primary %s)" % (
        res, "PASS" if r0 else "FAIL", "PASS" if r1 else "FAIL",
        "PASS" if r2 else "FAIL", r3s, primary or "none"))
    return "\n".join(L), res


def main(argv):
    directory = argv[1] if len(argv) > 1 else HERE
    text, _ = score(directory)
    print(text)
    out = os.path.join(directory, "score_p4_3d_thinlimit_r2_%s.txt" % DATE)
    with open(out, "w", newline="\n") as fh:
        fh.write(text + "\n")
    print("\nwrote", out)


if __name__ == "__main__":
    main(sys.argv)
