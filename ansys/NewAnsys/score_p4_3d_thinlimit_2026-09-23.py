#!/usr/bin/env python3
"""Score the Paper 4 3-D thin-limit check. Gates fixed before any run.

Settles (or not) the G3 FAIL of job 2530423 (LESSONS 18.232): elastic
n >= 1 FE roots of the layered ring 0.5-1.2% below the Kirchhoff 4x4.

Reads the five bare-steel decks' transcripts and outputs and
targets_p4_3d_thinlimit_2026-09-23.json. Parsing and mode identification
reuse score_p4_3d_nge1_2026-09-23.py unchanged (DCT-I of UZ over theta,
purity >= 0.90, rigid < 1 Hz, flexural if RMS UZ > 2x RMS in-plane).

delta_n(t) = f_FE(t) / f_Kirchhoff(t) - 1, fine mesh, first radial root.
Gated set n = 0..4.

  T0  Transcripts clean; element/symmetry counts exact; 3 rigid modes;
      every gated (n,1) identified in all five decks.
  T1  Mesh: |f_coarse / f_fine - 1| <= 0.05% at t = 20 and t = 5.
      A t that fails makes T3 UNRESOLVED.
  T2  Shrinks with thickness: |delta_n(5)| <= 0.5 |delta_n(20)| for every
      gated n. (Pure t^2 gives 0.0625, pure t gives 0.25; a t-independent
      operator error gives 1.)
  T3  Thin limit (decisive): the quadratic through the three fine points,
      delta = c + a t + b t^2, has |c| <= 0.05% at every gated n.
      c = (8/3) delta(5) - 2 delta(10) + (1/3) delta(20).

Reading, fixed now:
  T2 and T3 PASS  -> the n >= 1 thin-plate operator (with the effective-
      shear row) is FE-confirmed in the Kirchhoff limit; the G3 FAIL is
      thickness physics and may be reported as such.
  T3 FAIL at some n -> an n-dependent operator error of about c at that
      n. Stop; do not claim the elastic n >= 1 roots.

Reported, not gated: a and b per n (a is the first-order, edge-type part),
and the t = 20 bare-steel offsets beside the layered-ring G3 offsets.

Usage (from Ansys/NewAnsys/):  python3 score_p4_3d_thinlimit_2026-09-23.py [dir]
"""
import importlib.util
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DATE = "2026-09-23"
TARGETS = os.path.join(HERE, "targets_p4_3d_thinlimit_%s.json" % DATE)

DECKS = [(20, "coarse"), (20, "fine"), (10, "fine"), (5, "fine"), (5, "coarse")]
COUNTS = {"coarse": dict(elements=3456, symmetry_nodes=394),
          "fine": dict(elements=27648, symmetry_nodes=1362)}
GATED_N = (0, 1, 2, 3, 4)
N_RIGID = 3
MESH_BAR = 0.0005
SHRINK_BAR = 0.5
INTERCEPT_BAR = 0.0005
# layered-ring G3 offsets from job 2530423 (score_p4_3d_nge1_2026-09-23.txt)
G3_H112 = {0: -0.00198, 1: -0.01022, 2: -0.00533, 3: -0.00647, 4: -0.00935}


def _nge1():
    p = os.path.join(HERE, "score_p4_3d_nge1_%s.py" % DATE)
    spec = importlib.util.spec_from_file_location("nge1_scorer", p)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


S = _nge1()


def stem(t, mesh):
    return "ansys_p4_3d_bare_t%d_%s" % (t, mesh)


def load(directory, t, mesh):
    st = stem(t, mesh)
    d = dict(t=t, mesh=mesh, reasons=[], freqs=[], modes=[], table={})
    d["reasons"] += S.transcript_problems(
        os.path.join(directory, "%s_%s_out.txt" % (st, DATE)),
        os.path.join(directory, "%s_%s.inp.QUEUE_OK" % (st, DATE)))
    paths = [os.path.join(directory, st + x) for x in
             ("_modes.txt", "_samples.txt", "_samplecoords.txt")]
    miss = [os.path.basename(p) for p in paths if not os.path.isfile(p)]
    if miss:
        d["reasons"].append("missing " + ", ".join(miss))
        d["not_run"] = True
        return d
    counts, freqs = S.parse_modes(paths[0])
    d["freqs"] = freqs
    for k, want in COUNTS[mesh].items():
        if counts.get(k) != want:
            d["reasons"].append("%s = %s, expected %d" % (k, counts.get(k), want))
    coords = S.parse_coords(paths[2])
    if len(coords) != S.NSAMP:
        d["reasons"].append("%d sample coords" % len(coords))
    d["modes"] = S.classify(freqs, S.parse_samples(paths[1]), coords)
    nr = sum(1 for m in d["modes"] if m["kind"] == "rigid")
    if nr != N_RIGID:
        d["reasons"].append("%d rigid modes, expected %d" % (nr, N_RIGID))
    d["table"] = S.assign(d["modes"])
    for n in GATED_N:
        rec = d["table"].get((n, 1))
        if rec is None:
            d["reasons"].append("gated (n=%d, s=1) not identified" % n)
        elif rec["purity"] < S.PURITY:
            d["reasons"].append("(n=%d,1) purity %.3f" % (n, rec["purity"]))
    return d


def score(directory, targets_path=TARGETS):
    with open(targets_path) as f:
        T = json.load(f)
    D = {(t, m): load(directory, t, m) for t, m in DECKS}
    L = []
    a = L.append
    a("Paper 4 3-D thin-limit score  (%s)" % directory)
    t0 = True
    missing = False
    for key in DECKS:
        d = D[key]
        a("t=%-2d %-6s raw Hz: %s" % (key[0], key[1],
          ", ".join("%.10g" % x for x in d["freqs"]) or "(none)"))
        for rec in d["modes"]:
            a("    mode %2d %14.8f Hz  %-8s n=%-3s purity=%s" % (
                rec["mode"], rec["f"], rec["kind"], rec["n"],
                "%.4f" % rec["purity"] if rec["purity"] is not None else "-"))
        for r in d["reasons"]:
            a("    T0: %s" % r)
        t0 &= not d["reasons"]
        missing |= bool(d.get("not_run"))
    a("")
    if missing:
        a("T0: NOT_RUN. RESULT: NOT_RUN")
        return "\n".join(L), "NOT_RUN"
    a("T0: %s" % ("PASS" if t0 else "FAIL"))

    def f(t, mesh, n, s=1):
        rec = D[(t, mesh)]["table"].get((n, s))
        return None if rec is None else rec["f"]

    def delta(t, n, s=1):
        fe = f(t, "fine", n, s)
        k = T["t%d" % t].get("%d_%d" % (n, s))
        return None if fe is None or k is None else fe / k["el_hz"] - 1.0

    a("")
    t1_ok = {20: True, 5: True}
    for t in (20, 5):
        for n in GATED_N:
            fc, ff = f(t, "coarse", n), f(t, "fine", n)
            if fc is None or ff is None:
                t1_ok[t] = False
                a("T1 t=%d n=%d: missing" % (t, n))
                continue
            x = fc / ff - 1.0
            ok = abs(x) <= MESH_BAR
            t1_ok[t] &= ok
            a("T1 t=%-2d n=%d: coarse %.8f fine %.8f  %+.4f%%  %s"
              % (t, n, fc, ff, 100 * x, "PASS" if ok else "FAIL"))
    t1 = all(t1_ok.values())
    a("T1: %s" % ("PASS" if t1 else "FAIL"))

    a("")
    a("delta_n(t) = FE/Kirchhoff - 1 (fine), with the quadratic fit c + a t + b t^2:")
    t2 = True
    t3 = True
    for n in GATED_N:
        d20, d10, d5 = delta(20, n), delta(10, n), delta(5, n)
        if None in (d20, d10, d5):
            a("  n=%d: missing" % n)
            t2 = t3 = False
            continue
        c = (8.0 / 3.0) * d5 - 2.0 * d10 + (1.0 / 3.0) * d20
        # solve for a, b through the same three points (t in mm)
        # d10 - c = 10a + 100b ; d20 - c = 20a + 400b
        b = ((d20 - c) - 2.0 * (d10 - c)) / 200.0
        aa = ((d10 - c) - 100.0 * b) / 10.0
        shrink = abs(d5) / abs(d20) if d20 != 0 else float("inf")
        ok2 = shrink <= SHRINK_BAR
        ok3 = abs(c) <= INTERCEPT_BAR
        t2 &= ok2
        t3 &= ok3
        a("  n=%d: t=20 %+.4f%%  t=10 %+.4f%%  t=5 %+.4f%%  |d5/d20|=%.3f %s  "
          "c=%+.4f%% %s  a*20=%+.4f%%  b*400=%+.4f%%   (layered 1/12 G3: %+.3f%%)"
          % (n, 100 * d20, 100 * d10, 100 * d5, shrink, "PASS" if ok2 else "FAIL",
             100 * c, "PASS" if ok3 else "FAIL", 100 * aa * 20, 100 * b * 400,
             100 * G3_H112[n]))
    a("T2: %s" % ("PASS" if t2 else "FAIL"))
    if not t1:
        t3s = "UNRESOLVED (T1 failed)"
    else:
        t3s = "PASS" if t3 else "FAIL"
    a("T3: %s" % t3s)

    a("")
    a("Information only: other identified modes, fine mesh, delta at t = 20/10/5")
    for key in sorted(T["t20"], key=lambda k: T["t20"][k]["el_hz"]):
        n, s = (int(x) for x in key.split("_"))
        if s == 1 and n in GATED_N:
            continue
        ds = [delta(t, n, s) for t in (20, 10, 5)]
        a("  n=%d s=%d: %s" % (n, s, "  ".join(
            "%+.4f%%" % (100 * x) if x is not None else "n/a" for x in ds)))

    if t0 and t1 and t2 and t3:
        res = "PASS"
    elif t0 and t2 and t3s.startswith("UNRESOLVED"):
        res = "UNRESOLVED"
    else:
        res = "FAIL"
    a("")
    a("RESULT: %s  (T0 %s, T1 %s, T2 %s, T3 %s)" % (
        res, "PASS" if t0 else "FAIL", "PASS" if t1 else "FAIL",
        "PASS" if t2 else "FAIL", t3s))
    return "\n".join(L), res


def main(argv):
    directory = argv[1] if len(argv) > 1 else HERE
    text, _ = score(directory)
    print(text)
    out = os.path.join(directory, "score_p4_3d_thinlimit_%s.txt" % DATE)
    with open(out, "w", newline="\n") as fh:
        fh.write(text + "\n")
    print("\nwrote", out)


if __name__ == "__main__":
    main(sys.argv)
