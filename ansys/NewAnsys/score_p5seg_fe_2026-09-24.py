# -*- coding: utf-8 -*-
"""
score_p5seg_fe_2026-09-24.py -- score the Paper 5 segmented-electrode
harmonic FE (ansys_p5seg_*) against the pre-registered bars in
targets_p5seg_2026-09-24.json. Python 3.6 compatible, stdlib only.
Run from Ansys/NewAnsys/ after the queue ends:
    python3 score_p5seg_fe_2026-09-24.py | tee score_p5seg_fe_2026-09-24.txt

R = Q_FE,face / Q_K,face (physical electrode charge on the top inner
segment [r_i, r_*], corner node at r_* split half/half).
Pre-registration: PAPER5_E15_SEGMENT_CHARGE_DERIVATION_2026-09-24.md Sec 5.
"""
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
_TARGETS = next((a for a in sys.argv[1:] if a.endswith(".json")), "targets_p5seg_2026-09-24.json")
T = json.load(open(os.path.join(HERE, _TARGETS)))
KF = T.get("kappa_floor", 0.0)
B = T["bars"]
RSTAR = T["r_star_mm"]
U = T["charge_unit_to_C"]
TOL = 1e-3

_fix = re.compile(r"(\d)([+-]\d{2,3})$")


def num(tok):
    tok = tok.strip()
    if "E" not in tok.upper():
        tok = _fix.sub(r"\1E\2", tok)
    return float(tok)


def read_ls(stem, ls):
    p = os.path.join(HERE, "%s_ls%d.txt" % (stem, ls))
    rows = []
    with open(p) as fh:
        for line in fh:
            parts = line.split()
            if len(parts) != 7:
                continue
            try:
                rows.append([num(x) for x in parts])
            except ValueError:
                continue
    return rows


def seg(rows, col):
    """(half-split, inclusive, exclusive, full-face) sums for column col."""
    excl = sum(r[col] for r in rows if r[1 if col < 4 else 4] < RSTAR - TOL)
    rcol = 1 if col < 4 else 4
    on = [r[col] for r in rows if abs(r[rcol] - RSTAR) <= TOL]
    if len(on) != 1:
        raise ValueError("expected exactly one node at r_star, got %d" % len(on))
    full = sum(r[col] for r in rows)
    return excl + 0.5 * on[0], excl + on[0], excl, full


def deck_health(deck):
    out = os.path.join(HERE, deck.replace(".inp", "_out.txt"))
    if not os.path.exists(out):
        return False, "no _out.txt"
    txt = open(out, errors="replace").read()
    m = re.search(r"NUMBER OF ERROR\s+MESSAGES ENCOUNTERED=\s*(\d+)", txt)
    nerr = int(m.group(1)) if m else -1
    ok = (nerr == 0) and ("ERROR: no input file specified" not in txt[:4000])
    return ok, "errors=%s" % nerr


def main():
    res = {}
    allok = True
    for d in T["decks"]:
        ok, msg = deck_health(d["deck"])
        key = (d["c44_key"], d["ndiv_r"])
        try:
            cal = read_ls(d["stem"], 1)
            qcal = seg(cal, 2)[1]      # driven electrode incl. its r_star node
            Rs, extra = [], []
            for lsd in T["load_steps"]:
                rows = read_ls(d["stem"], lsd["ls"])
                th, ti, te, tf = seg(rows, 2)
                bh = seg(rows, 5)[0]
                im = max(abs(r[3]) for r in rows) / max(abs(r[2]) for r in rows)
                qfe = th * U
                R = qfe / lsd["Q_K_face_C"]
                Rs.append(R)
                extra.append(dict(ls=lsd["ls"], omega=lsd["omega"], R=R,
                                  R_incl=ti * U / lsd["Q_K_face_C"],
                                  R_excl=te * U / lsd["Q_K_face_C"],
                                  sym=abs(th - bh) / max(abs(th), KF * abs(lsd["Q_K_face_C"]) / U),
                                  full=abs(tf) / max(abs(th), KF * abs(lsd["Q_K_face_C"]) / U), im=im,
                                  kappa=lsd["kappa"][d["c44_key"]]))
            res[key] = dict(ok=ok, msg=msg, qcal=qcal, rows=extra)
        except (IOError, OSError, ValueError) as exc:
            res[key] = dict(ok=False, msg="%s; %s" % (msg, exc), rows=[])
        allok = allok and res[key]["ok"]
        print("%-40s %s %s" % (d["deck"], "RUN_OK" if res[key]["ok"] else "RUN_BAD",
                               res[key]["msg"]))
    print("G_run:", "PASS" if allok else "FAIL")
    if not allok:
        print("VERDICT: INCOMPLETE -- fix the run before reading any ratio.")
        return 1

    # AMENDED 2026-09-24 after job 2531984 (LESSONS Sec 18.239). LS1 is the
    # sign calibration. The pre-registered text assumed it would show
    # CHRG = +electrode charge. It came out negative in all six decks, as the
    # earlier Paper 4 lever-9 deck already implied (job 2494863: capacitive
    # CHRG < 0 at 80/140 Hz for +1 V). So CHRG = -electrode charge in this
    # ANSYS build. The gate now requires one consistent sign across decks
    # plus |C_seg| inside the units window, and that sign is applied to
    # every R. The raw-sign score is kept in score_p5seg_fe_2026-09-24_raw_precal.txt.
    signs = set(1 if v["qcal"] > 0 else -1 for v in res.values())
    gcal = len(signs) == 1
    SGN = signs.pop() if gcal else 1
    cap = [abs(v["qcal"]) * U for v in res.values()]
    lo, hi = B["cal_C_range_F"]
    gcap = all(lo <= c <= hi for c in cap)
    gcal = gcal and gcap
    print("G_cal (one sign across decks, |C_seg| in %.0e-%.0e F):" % (lo, hi),
          "PASS" if gcal else "FAIL", " CHRG sign convention =",
          "+charge" if SGN > 0 else "-charge (applied)",
          " |C_seg| [F] =", ", ".join("%.4e" % c for c in cap))
    for v in res.values():
        for r in v["rows"]:
            for k in ("R", "R_incl", "R_excl"):
                r[k] = SGN * r[k]
    gsym = all(r["sym"] < B["G_sym_rel"] and r["full"] < B["G_sym_rel"]
               for v in res.values() for r in v["rows"])
    print("G_sym (top=bottom, full-face ~ 0):", "PASS" if gsym else "FAIL")

    verdicts = {}
    gmesh = True
    for ck in ("c44_73", "c44_26"):
        print("\n==== %s ====" % ck)
        print("  ndiv   omega    R(half)     R(incl)     R(excl)     kappa     R/kappa   |Im/Re|  sym      full")
        for nr in (40, 80, 160):
            for r in res[(ck, nr)]["rows"]:
                print("  %4d %7.1f  %+.6f  %+.6f  %+.6f  %+.5f  %.5f  %.1e  %.1e  %.1e"
                      % (nr, r["omega"], r["R"], r["R_incl"], r["R_excl"],
                         r["kappa"], r["R"] / r["kappa"], r["im"], r["sym"], r["full"]))
        fine = res[(ck, 160)]["rows"]
        mid = res[(ck, 80)]["rows"]
        for a, b in zip(mid, fine):
            if abs(a["R"] - b["R"]) / max(abs(b["R"]), KF) >= B["G_mesh_rel"]:
                gmesh = False
        corr = all(abs(r["R"] - r["kappa"]) / max(abs(r["kappa"]), KF) < B["corrected_rel"]
                   for r in fine)
        kirch = all(abs(r["R"] - 1) < B["kirchhoff_rel"] for r in fine)
        verdicts[ck] = "CORRECTED" if corr else ("KIRCHHOFF" if kirch else "UNRESOLVED")
        print("  finest-mesh verdict:", verdicts[ck])
    print("\nG_mesh (m80 vs m160 within %.0f%%):" % (100 * B["G_mesh_rel"]),
          "PASS" if gmesh else "FAIL")
    gates = gcal and gsym and gmesh
    if gates and all(v == "CORRECTED" for v in verdicts.values()):
        final = "CORRECTED -- e15 shear term confirmed at both c44"
    elif gates and all(v == "KIRCHHOFF" for v in verdicts.values()):
        final = "KIRCHHOFF -- Paper 5 Eq. Ysense confirmed as printed"
    else:
        final = "UNRESOLVED -- investigate, do not retune (%s; gates cal=%s sym=%s mesh=%s)" % (
            verdicts, gcal, gsym, gmesh)
    print("VERDICT:", final)
    return 0


if __name__ == "__main__":
    sys.exit(main())
