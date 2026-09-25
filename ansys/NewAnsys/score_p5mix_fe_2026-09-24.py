# -*- coding: utf-8 -*-
"""
score_p5mix_fe_2026-09-24.py -- score the Paper 5 mixed-edge segmented-
electrode harmonic FE (ansys_p5mix_*) against the bars fixed before the run
in targets_p5mix_2026-09-24.json. Python 3.6 compatible, stdlib only.
Run from Ansys/NewAnsys/ after the queue ends:
    python3 score_p5mix_fe_2026-09-24.py targets_p5mix_e31m_2026-09-24.json | tee score_p5mix_e31m_fe_2026-09-24.txt

Pre-registration: PAPER5_MIXED_EDGE_SENSING_2026-09-24.md Sec 4.
R(cut) = Q_FE,face / Q_K,face: physical electrode charge on the top face
segment [r_i, r_cut] (corner node at r_cut split half/half), CHRG sign
fixed by the LS1 calibration (ANSYS 2025R1: CHRG = -charge, job 2531984).

Part A (claims, interior cuts 300 and 450 mm), per BC and c44:
  CORRECTED  finest mesh |R/kappa_closed - 1| < 5% at every omega, both cuts
  KIRCHHOFF  finest mesh |R - 1| < 10% at every omega, both cuts
  else UNRESOLVED (investigate, do not retune).
Part B (descriptive, cuts within 2.5H of a rim): R/kappa_M1 and R/kappa_M2
  at the finest mesh; label M1 / M2 / BOTH / NEITHER at 3%, or
  UNCONVERGED if m160 and m320 differ by more than 2%. Free-rim cuts
  (M1 = M2) are a control and should read BOTH.
"""
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
_TARGETS = next((a for a in sys.argv[1:] if a.endswith(".json")), "targets_p5mix_2026-09-24.json")
T = json.load(open(os.path.join(HERE, _TARGETS)))
KF = T.get("kappa_floor", 0.0)


def rel(x, ref):
    """|x/ref - 1| measured against max(|ref|, kappa_floor) (targets JSON)."""
    return abs(x - ref) / max(abs(ref), KF)
B = T["bars"]
U = T["charge_unit_to_C"]
CUTS = T["cuts_mm"]
RCAL = T["r_star_cal_mm"]
TOL = 1e-3
FINE, MID = 320, 160

_fix = re.compile(r"(\d)([+-]\d{2,3})$")


def num(tok):
    tok = tok.strip()
    if "E" not in tok.upper():
        tok = _fix.sub(r"\1E\2", tok)
    return float(tok)


def read_ls(stem, ls):
    rows = []
    with open(os.path.join(HERE, "%s_ls%d.txt" % (stem, ls))) as fh:
        for line in fh:
            parts = line.split()
            if len(parts) != 7:
                continue
            try:
                rows.append([num(x) for x in parts])
            except ValueError:
                continue
    return rows


def seg(rows, col, rcut):
    """(half-split, inclusive, full-face) sum of column col for r <= rcut."""
    rcol = 1 if col < 4 else 4
    excl = sum(r[col] for r in rows if r[rcol] < rcut - TOL)
    on = [r[col] for r in rows if abs(r[rcol] - rcut) <= TOL]
    if len(on) != 1:
        raise ValueError("expected one node at r=%.4f, got %d" % (rcut, len(on)))
    return excl + 0.5 * on[0], excl + on[0], sum(r[col] for r in rows)


def deck_health(deck):
    out = os.path.join(HERE, deck.replace(".inp", "_out.txt"))
    if not os.path.exists(out):
        return False, "no _out.txt"
    txt = open(out, errors="replace").read()
    m = re.search(r"NUMBER OF ERROR\s+MESSAGES ENCOUNTERED=\s*(\d+)", txt)
    nerr = int(m.group(1)) if m else -1
    return nerr == 0, "errors=%s" % nerr


def main():
    res = {}
    allok = True
    for d in T["decks"]:
        ok, msg = deck_health(d["deck"])
        key = (d["bc"], d["c44_key"], d["ndiv_r"])
        rec = dict(ok=ok, msg=msg, rows=[])
        try:
            rec["qcal"] = seg(read_ls(d["stem"], 1), 2, RCAL)[1]
            for lsd in T["load_steps"][d["bc"]]:
                rows = read_ls(d["stem"], lsd["ls"])
                out = dict(ls=lsd["ls"], omega=lsd["omega"], cut={})
                t300 = seg(rows, 2, CUTS["int_300"])
                b300 = seg(rows, 5, CUTS["int_300"])[0]
                den = max(abs(t300[0]), KF * abs(lsd["cuts"]["int_300"]["Q_K_face"]) / U)
                out["sym"] = abs(t300[0] - b300) / den
                out["full"] = abs(t300[2]) / den
                out["im"] = (max(abs(r[3]) for r in rows)
                             / max(abs(r[2]) for r in rows))
                for ck, rmm in CUTS.items():
                    tgt = lsd["cuts"][ck]
                    q = seg(rows, 2, rmm)[0] * U
                    out["cut"][ck] = dict(
                        Rraw=q / tgt["Q_K_face"],
                        k1=tgt["kappa_model"][d["c44_key"]],
                        k2=tgt.get("kappa_M2", tgt["kappa_model"])[d["c44_key"]])
                rec["rows"].append(out)
        except (IOError, OSError, ValueError, KeyError) as exc:
            rec["ok"] = False
            rec["msg"] = "%s; %s" % (msg, exc)
        res[key] = rec
        allok = allok and rec["ok"]
        print("%-42s %s %s" % (d["deck"], "RUN_OK" if rec["ok"] else "RUN_BAD", rec["msg"]))
    print("G_run:", "PASS" if allok else "FAIL")
    if not allok:
        print("VERDICT: INCOMPLETE -- fix the run before reading any ratio.")
        return 1

    signs = set(1 if v["qcal"] > 0 else -1 for v in res.values())
    gcal = len(signs) == 1
    SGN = signs.pop() if gcal else 1
    lo, hi = B["cal_C_range_F"]
    caps = [abs(v["qcal"]) * U for v in res.values()]
    gcal = gcal and all(lo <= c <= hi for c in caps)
    print("G_cal (one sign, |C_seg| in %.0e-%.0e F): %s  CHRG = %s  |C_seg| range %.4e-%.4e F"
          % (lo, hi, "PASS" if gcal else "FAIL",
             "+charge" if SGN > 0 else "-charge (applied)", min(caps), max(caps)))
    for v in res.values():
        for r in v["rows"]:
            for c in r["cut"].values():
                c["R"] = SGN * c["Rraw"]
    # sym/full are relative to the r* = 300 segment charge; at e31 = -4.1,
    # c44 = 73 that charge is ~0.055 Q_K, so scale the bar by the same floor.
    gsym = all(r["sym"] < B["G_sym_rel"] and r["full"] < B["G_sym_rel"]
               for v in res.values() for r in v["rows"])
    print("G_sym (top=bottom and full-face ~ 0, at r*=300):", "PASS" if gsym else "FAIL")

    gmesh = True
    verdicts = {}
    kc = T["kappa_closed"]
    print("\n==== Part A: interior cuts (claims) ====")
    for bc in ("FF", "CC", "CF", "FC"):
        for ck in ("c44_73", "c44_26"):
            print("\n-- %s  %s  kappa_closed=%+.5f --" % (bc, ck, kc[ck]))
            print("  ndiv   omega   R300       R300/k    R450       R450/k    |Im/Re|  sym     full")
            for nr in (80, MID, FINE):
                for r in res[(bc, ck, nr)]["rows"]:
                    a, b = r["cut"]["int_300"]["R"], r["cut"]["int_450"]["R"]
                    print("  %4d %8.2f  %+.5f  %.5f  %+.5f  %.5f  %.1e  %.1e  %.1e"
                          % (nr, r["omega"], a, a / kc[ck], b, b / kc[ck],
                             r["im"], r["sym"], r["full"]))
            fine = res[(bc, ck, FINE)]["rows"]
            mid = res[(bc, ck, MID)]["rows"]
            for a, b in zip(mid, fine):
                for cut in T["interior_cuts"]:
                    x, y = a["cut"][cut]["R"], b["cut"][cut]["R"]
                    if rel(x, y) >= B["G_mesh_rel"]:
                        gmesh = False
            corr = all(rel(r["cut"][c]["R"], kc[ck]) < B["corrected_rel"]
                       for r in fine for c in T["interior_cuts"])
            kir = all(abs(r["cut"][c]["R"] - 1) < B["kirchhoff_rel"]
                      for r in fine for c in T["interior_cuts"])
            verdicts[(bc, ck)] = "CORRECTED" if corr else ("KIRCHHOFF" if kir else "UNRESOLVED")
            print("  verdict:", verdicts[(bc, ck)])
    print("\nG_mesh (m%d vs m%d within %.0f%%, interior cuts): %s"
          % (MID, FINE, 100 * B["G_mesh_rel"], "PASS" if gmesh else "FAIL"))

    print("\n==== Part B: near-rim cuts (descriptive) ====")
    print("  bc  c44     cut       rim      R/k_M1 range        R/k_M2 range        mesh    label")
    for bc in ("FF", "CC", "CF", "FC"):
        for kind, cuts in (("clamped", T["clamped_rim_cuts"][bc]),
                           ("free", T["free_rim_cuts"][bc])):
            for cut in cuts:
                for ck in ("c44_73", "c44_26"):
                    fine = res[(bc, ck, FINE)]["rows"]
                    mid = res[(bc, ck, MID)]["rows"]
                    r1 = [r["cut"][cut]["R"] / r["cut"][cut]["k1"] for r in fine]
                    r2 = [r["cut"][cut]["R"] / r["cut"][cut]["k2"] for r in fine]
                    conv = all(rel(a["cut"][cut]["R"], b["cut"][cut]["R"]) < B["rim_mesh_rel"]
                               for a, b in zip(mid, fine))
                    m1 = all(rel(r["cut"][cut]["R"], r["cut"][cut]["k1"]) < B["rim_model_rel"]
                             for r in fine)
                    m2 = all(rel(r["cut"][cut]["R"], r["cut"][cut]["k2"]) < B["rim_model_rel"]
                             for r in fine)
                    if not conv:
                        lab = "UNCONVERGED"
                    else:
                        lab = {(True, True): "BOTH", (True, False): "M1",
                               (False, True): "M2", (False, False): "NEITHER"}[(m1, m2)]
                    print("  %s  %s  %-8s  %-7s  %.4f..%.4f   %.4f..%.4f   %-6s  %s"
                          % (bc, ck, cut, kind, min(r1), max(r1), min(r2), max(r2),
                             "ok" if conv else "NO", lab))

    gates = gcal and gsym and gmesh
    allc = all(v == "CORRECTED" for v in verdicts.values())
    allk = all(v == "KIRCHHOFF" for v in verdicts.values())
    if gates and allc:
        final = "CORRECTED -- kappa carries over to C-C, C-F and F-C at both c44 (F-F control reproduced)"
    elif gates and allk:
        final = "KIRCHHOFF -- contradicts job 2531984; investigate before anything else"
    else:
        final = "UNRESOLVED -- investigate, do not retune (verdicts=%s; gates cal=%s sym=%s mesh=%s)" % (
            dict(("%s/%s" % k, v) for k, v in verdicts.items()), gcal, gsym, gmesh)
    print("\nVERDICT:", final)
    return 0


if __name__ == "__main__":
    sys.exit(main())
