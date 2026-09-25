# -*- coding: utf-8 -*-
"""
score_parashar_fe_2026-09-24.py -- score the Parashar 2013 n = 0 FE scout
(ansys_parashar_*). Python 3.6, stdlib only. Run from Ansys/NewAnsys/:
    python3 score_parashar_fe_2026-09-24.py | tee score_parashar_fe_2026-09-24.txt
Pre-registration: PARASHAR_PHASEB_SCOUT_2026-09-24.md.
Flexural n = 0 mode at the free outer rim: UY_top, UY_bot same sign and
UX_top, UX_bot opposite sign. The k-th flexural mode above 1 Hz is s = k.
"""
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
T = json.load(open(os.path.join(HERE, "targets_parashar_2026-09-24.json")))
B = T["bars"]
_fix = re.compile(r"(\d)([+-]\d{2,3})$")


def num(t):
    t = t.strip()
    if "E" not in t.upper():
        t = _fix.sub(r"\1E\2", t)
    return float(t)


def modes(stem):
    out = []
    with open(os.path.join(HERE, stem + "_modes.txt")) as fh:
        for line in fh:
            p = line.split()
            if len(p) != 7:
                continue
            try:
                out.append([num(x) for x in p])
            except ValueError:
                continue
    return out


def health(deck):
    out = os.path.join(HERE, deck.replace(".inp", "_out.txt"))
    if not os.path.exists(out):
        return False, "no _out.txt"
    txt = open(out, errors="replace").read()
    m = re.search(r"NUMBER OF ERROR\s+MESSAGES ENCOUNTERED=\s*(\d+)", txt)
    n = int(m.group(1)) if m else -1
    return n == 0, "errors=%s" % n


def main():
    fe = {}
    ok_all = True
    for d in T["decks"]:
        ok, msg = health(d["deck"])
        try:
            fe[(d["bc"], d["piezo"], d["ndiv_r"])] = modes(d["stem"])
        except (IOError, OSError) as exc:
            ok, msg = False, "%s; %s" % (msg, exc)
        ok_all = ok_all and ok
        print("%-42s %s %s" % (d["deck"], "RUN_OK" if ok else "RUN_BAD", msg))
    print("G_run:", "PASS" if ok_all else "FAIL")
    if not ok_all:
        print("VERDICT: INCOMPLETE")
        return 1
    verdict = {}
    gmesh = True
    for bc, g in T["geometry"].items():
        print("\n==== %s  %s ====" % (bc, g["table"]))
        print("  all modes, m120 SC: mode f_kHz class")
        for r in fe[(bc, "sc", 120)]:
            flex = r[1] > 1 and r[3] * r[5] > 0 and r[2] * r[4] < 0
            ext = r[1] > 1 and r[2] * r[4] > 0 and r[3] * r[5] < 0
            print("   %3d %10.4f  %s" % (r[0], r[1] / 1e3,
                                        "flexural" if flex else ("extensional" if ext else "other/rigid")))
        fl = {}
        for key in (("sc", 60), ("sc", 120), ("e0", 120)):
            fl[key] = [r[1] / 1e3 for r in fe[(bc,) + key]
                       if r[1] > 1 and r[3] * r[5] > 0 and r[2] * r[4] < 0]
        print("  s   Parashar_model  exp       FE_sc_m120  FE/model  FE/exp   FE_sc_m60  FE_e0_m120  piezo_shift%")
        good_model = good_exp = True
        for i, (m, e) in enumerate(zip(g["model"], g["exp"])):
            if i >= len(fl[("sc", 120)]):
                print("  %d  (no FE flexural mode found)" % (i + 1))
                good_model = good_exp = False
                continue
            f = fl[("sc", 120)][i]
            fc = fl[("sc", 60)][i] if i < len(fl[("sc", 60)]) else float("nan")
            f0 = fl[("e0", 120)][i] if i < len(fl[("e0", 120)]) else float("nan")
            print("  %d  %10.3f  %10s  %10.4f  %.4f  %8s  %10.4f  %10.4f  %+.3f"
                  % (i + 1, m, "%.3f" % e if e else "-", f, f / m,
                     "%.4f" % (f / e) if e else "-", fc, f0, 100 * (f / f0 - 1)))
            if abs(fc / f - 1) >= B["mesh_rel"]:
                gmesh = False
            if i < 2:
                good_model = good_model and abs(f / m - 1) < B["model_rel_s12"]
                if e:
                    good_exp = good_exp and abs(f / e - 1) < B["model_rel_s12"]
        verdict[bc] = ("MODEL_AND_EXP" if good_model and good_exp else
                       "MODEL" if good_model else "EXP" if good_exp else "NEITHER")
        print("  reading (s = 1, 2 within 3%):", verdict[bc])
    print("\nG_mesh (m60 vs m120 within %.1f%%): %s" % (100 * B["mesh_rel"], "PASS" if gmesh else "FAIL"))
    print("VERDICT:", verdict)
    return 0


if __name__ == "__main__":
    sys.exit(main())
