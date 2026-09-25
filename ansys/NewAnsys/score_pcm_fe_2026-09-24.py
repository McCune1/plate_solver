# -*- coding: utf-8 -*-
"""
score_pcm_fe_2026-09-24.py -- score the contour-mode FE (ansys_pcm_*)
against targets_pcm_2026-09-24.json (bars fixed before the run).
Python 3.6 compatible, stdlib only. Run from Ansys/NewAnsys/:
    python3 score_pcm_fe_2026-09-24.py | tee score_pcm_fe_2026-09-24.txt
Pre-registration: CONTOUR_MODE_N0_DERIVATION_2026-09-24.md Sec 5.

Extensional (radial) mode = f > 1 Hz, UX_top and UX_bot same sign at
r = RO, and |UX_mid| > |UY_top| (bending has UX_mid ~ 0).
"""
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
T = json.load(open(os.path.join(HERE, "targets_pcm_2026-09-24.json")))
B = T["bars"]
U = T["charge_unit_to_C"]
_fix = re.compile(r"(\d)([+-]\d{2,3})$")


def num(tok):
    tok = tok.strip()
    if "E" not in tok.upper():
        tok = _fix.sub(r"\1E\2", tok)
    return float(tok)


def rows(path, ncol):
    out = []
    with open(path) as fh:
        for line in fh:
            p = line.split()
            if len(p) != ncol:
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


def radial(stem):
    rr = rows(os.path.join(HERE, stem + "_modes.txt"), 7)
    return [r[1] for r in rr
            if r[1] > 1.0 and r[2] * r[4] > 0 and abs(r[6]) > abs(r[3])]


def main():
    ok_all = True
    fe = {}
    for d in T["decks"]:
        ok, msg = health(d["deck"])
        key = (d["e31_key"], d["kind"], d["ndiv_r"])
        try:
            if d["kind"] in ("sc", "oc"):
                fe[key] = radial(d["stem"])
            else:
                qs = []
                for ls in range(1, len(d["harm_freqs_Hz"]) + 1):
                    rr = rows(os.path.join(HERE, "%s_ls%d.txt" % (d["stem"], ls)), 4)
                    qs.append(sum(r[2] for r in rr) * U)
                fe[key] = qs
        except (IOError, OSError, ValueError) as exc:
            ok, msg = False, "%s; %s" % (msg, exc)
        ok_all = ok_all and ok
        print("%-40s %s %s" % (d["deck"], "RUN_OK" if ok else "RUN_BAD", msg))
    print("G_run:", "PASS" if ok_all else "FAIL")
    if not ok_all:
        print("VERDICT: INCOMPLETE -- fix the run before reading any number.")
        return 1

    passes = {}
    for ck, case in sorted(T["cases"].items()):
        print("\n==== %s ====" % ck)
        modes = case["modes"]
        g = dict(cls=True, fr=True, fa=True, split=True, mesh=True, static=True, harm=True)
        for nr in (80, 160):
            if len(fe[(ck, "sc", nr)]) < 6 or len(fe[(ck, "oc", nr)]) < 6:
                g["cls"] = False
        if not g["cls"]:
            print("  G_class FAIL: fewer than 6 extensional modes found",
                  {k: len(v) for k, v in fe.items() if k[0] == ck and k[1] != "harm"})
            passes[ck] = False
            continue
        print("  m   fr_model    fr_FE(160)  dfr%      fa_model    fa_FE(160)  dfa%      "
              "split_mod%  split_FE%  ratio    dfr/(kH)^2  mesh80/160")
        for i, m in enumerate(modes):
            frF, faF = fe[(ck, "sc", 160)][i], fe[(ck, "oc", 160)][i]
            frC = fe[(ck, "sc", 80)][i]
            dfr, dfa = frF / m["fr_Hz"] - 1, faF / m["fa_Hz"] - 1
            spF = 100 * (faF / frF - 1)
            ratio = spF / m["split_pct"]
            dmesh = abs(frC / frF - 1)
            print("  %d  %10.4f  %10.4f  %+.4f  %10.4f  %10.4f  %+.4f  %9.5f  %9.5f  %.4f  %+.4f     %.1e"
                  % (m["m"], m["fr_Hz"], frF, 100 * dfr, m["fa_Hz"], faF, 100 * dfa,
                     m["split_pct"], spF, ratio, dfr / m["kH"] ** 2, dmesh))
            g["mesh"] = g["mesh"] and dmesh < B["mesh_fr_rel"]
            if m["m"] <= 3:
                g["fr"] = g["fr"] and abs(dfr) < B["fr_rel_m123"]
                g["fa"] = g["fa"] and abs(dfa) < B["fr_rel_m123"]
                if m["keff2"] >= B["keff2_min"]:
                    g["split"] = g["split"] and abs(ratio - 1) < B["split_rel_m123"]
        fr2 = modes[1]["fr_Hz"]
        print("  harmonic (top 1 V):  f_Hz   Q_FE[C]   Q_model[C]   ratio   (CHRG sign %+d applied)"
              % T["chrg_sign_expected"])
        for nr in (80, 160):
            qs = [T["chrg_sign_expected"] * q for q in fe[(ck, "harm", nr)]]
            freqs = [d for d in T["decks"] if d["e31_key"] == ck and d["kind"] == "harm"
                     and d["ndiv_r"] == nr][0]["harm_freqs_Hz"]
            cT = case["C_T_F"]
            rs = qs[0] / cT
            print("   m%-3d %9.2f  %+.6e  %+.6e  %.5f  (static: C_T)" % (nr, freqs[0], qs[0], cT, rs))
            if nr == 160:
                g["static"] = abs(rs - 1) < B["static_rel"]
            for f, q, hp in zip(freqs[1:], qs[1:], case["harm_points"]):
                r = q / hp["Q_top_C_per_V"]
                print("   m%-3d %9.2f  %+.6e  %+.6e  %.5f" % (nr, f, q, hp["Q_top_C_per_V"], r))
                if nr == 160 and f < fr2:
                    g["harm"] = g["harm"] and abs(r - 1) < B["harm_rel_below_fr2"]
        print("  gates:", g)
        passes[ck] = all(g.values())
    if all(passes.values()):
        final = "CONFIRMED -- thin-plate contour model matches FE at both e31 signs"
    else:
        final = "UNRESOLVED -- investigate, do not retune (%s)" % passes
    print("\nVERDICT:", final)
    return 0


if __name__ == "__main__":
    sys.exit(main())
