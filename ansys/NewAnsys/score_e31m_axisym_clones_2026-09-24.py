# -*- coding: utf-8 -*-
"""
score_e31m_axisym_clones_2026-09-24.py -- score the axisymmetric e31 = -4.1
clones of job 2532086 (LESSONS Sec 18.245) against the package model
(model_e31m_axisym_2026-09-24.json, consistent projection, both signs).

Scored quantity: the electrical split, which is what the parent jobs were
scored on (the absolute FE frequency carries a mesh/thickness offset that
is independent of e31):
  P4 F-F  OC/SC - 1 (scboth deck = both faces of each layer grounded), s = 1, 2,
          h1/2h = 1/12, 1/8, 1/5; and the NDIV_ZP = 4 refinement (zp4) decks
  P4 C-F  OC/SC - 1, s = 1
  P5 F-F  SC/e0 - 1 (e0 = coupling-off parent, job 2520476), s = 1
Bar: FE split / model split within 3% (relative) at -4.1, the same bar the
+4.1 parents passed; the +4.1 row is printed from the parent outputs as the
control. Reads the unzipped results_2026-09-24b/ folder for -4.1 outputs and
the parent outputs in this directory for +4.1.
"""
import json, os
HERE = os.path.dirname(os.path.abspath(__file__))
R = os.path.join(HERE, "results_2026-09-24b")
M = json.load(open(os.path.join(HERE, "model_e31m_axisym_2026-09-24.json")))
BAR = 0.03


def modes(path):
    out = {}
    for line in open(path):
        p = line.split()
        if len(p) >= 2:
            try:
                out[int(float(p[0]))] = float(p[1])
            except ValueError:
                pass
    return out


def f(stem, mode, sign):
    p = (os.path.join(R, stem + "_modes_e31m.txt") if sign < 0
         else os.path.join(HERE, stem + "_modes.txt"))
    return modes(p)[mode]


rows = []
ok_all = True
cases = [("P4 FF 1/12", "p4_12", "ansys_p4_epsrel_ff_lanb_oc", "ansys_p4_epsrel_ff_lanb_scboth"),
         ("P4 FF 1/8 ", "p4_8", "ansys_p4_epsrel_ff_lanb_oc_h18", "ansys_p4_epsrel_ff_lanb_scboth_h18"),
         ("P4 FF 1/5 ", "p4_5", "ansys_p4_epsrel_ff_lanb_oc_h15", "ansys_p4_epsrel_ff_lanb_scboth_h15"),
         ("P4 zp4 1/12", "p4_12", "ansys_p4_zp4_ff_oc_h112", "ansys_p4_zp4_ff_sc_h112"),
         ("P4 zp4 1/8 ", "p4_8", "ansys_p4_zp4_ff_oc_h18", "ansys_p4_zp4_ff_sc_h18"),
         ("P4 zp4 1/5 ", "p4_5", "ansys_p4_zp4_ff_oc_h15", "ansys_p4_zp4_ff_sc_h15")]
print("case          s  e31   FE split %   model split %   FE/model   verdict")
for name, key, oc, sc in cases:
    for s_idx, mode in ((0, 2), (1, 3)):
        for sign, tag in ((1, "+4.1"), (-1, "-4.1")):
            fe = f(oc, mode, sign) / f(sc, mode, sign) - 1
            mo = M[key][tag][s_idx]["oc_sc"]
            r = fe / mo
            v = "ok" if abs(r - 1) < BAR else "MISS"
            if sign < 0:
                ok_all = ok_all and v == "ok"
            print("%-12s %d  %s  %10.5f   %10.5f    %.4f    %s%s"
                  % (name, s_idx + 1, tag, 100 * fe, 100 * mo, r, v, "" if sign < 0 else " (control)"))
for sign, tag in ((1, "+4.1"), (-1, "-4.1")):
    fe = f("ansys_p4_epsrel_cf_lanb_oc", 1, sign) / f("ansys_p4_epsrel_cf_lanb_scboth", 1, sign) - 1
    mo = M["p4_cf"][tag][0]["oc_sc"]
    r = fe / mo
    v = "ok" if abs(r - 1) < BAR else "MISS"
    if sign < 0:
        ok_all = ok_all and v == "ok"
    print("%-12s %d  %s  %10.5f   %10.5f    %.4f    %s%s"
          % ("P4 CF 1/12 ", 1, tag, 100 * fe, 100 * mo, r, v, "" if sign < 0 else " (control)"))
e0 = modes(os.path.join(HERE, "ansys_p5_epsrel_ff_lanb_e0_modes.txt"))[2]
for sign, tag in ((1, "+4.1"), (-1, "-4.1")):
    for d in ("scboth", "oc"):
        fe = f("ansys_p5_epsrel_ff_lanb_" + d, 2, sign) / e0 - 1
        mo = M["p5_ff"][tag][0]["sc_el"]
        r = fe / mo
        v = "ok" if abs(r - 1) < BAR else "MISS"
        if sign < 0:
            ok_all = ok_all and v == "ok"
        print("%-12s %d  %s  %10.5f   %10.5f    %.4f    %s%s"
              % ("P5 FF " + d, 1, tag, 100 * fe, 100 * mo, r, v, "" if sign < 0 else " (control)"))
print("\nVERDICT (-4.1 rows, 3%% bar): %s" % ("PASS" if ok_all else "FAIL"))
