#!/usr/bin/env python3
"""Score the Liu 2002 solid-disk SC FE v2 decks (Paper 4 lever 11).

Reads, for bc in C, S: <prefix>_<bc>_2026-09-24_out.txt (transcript),
<prefix>_<bc>_modes.txt, _samples.txt, _samplecoords.txt, plus the
pre-registered targets_liudisk2_2026-09-24.json. Every ratio is computed
here from the raw Hz MAPDL printed; no deck-computed column is used.

Mode ID (fixed rule, same family as score_p4_3d_nge1_2026-09-23.py):
  rigid      f < 1 Hz (v2 expects NONE: C is clamped; S has UY sym plane,
             UZ midplane ring and a single-node UX pin)
  flexural   RMS UZ over top-face samples > 2 x RMS in-plane |U| there
  n          dominant cos(n theta) energy (13-point DCT-I over theta =
             0..180 step 15, summed over sampled radii); purity reported
  m          frequency rank among flexural modes of that n

Gates F0-F3 are defined in the targets JSON (fixed before the run).

Usage (from Ansys/NewAnsys/ after the queue ran):
    python3 score_liudisk2_fe_2026-09-24.py
    python3 score_liudisk2_fe_2026-09-24.py --selftest   # parses the v1
        (job 2530987) outputs and must reproduce the v1 SCORE table
"""
# (no __future__ import: compute-node python3 is 3.6)

import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))


def floats(line):
    out = []
    for tok in line.split():
        try:
            out.append(float(tok))
        except ValueError:
            return None
    return out


def read_modes(path):
    f = {}
    started = False
    for line in open(path):
        if line.strip().startswith("mode f_Hz"):
            started = True
            continue
        if started:
            v = floats(line)
            if v and len(v) == 2:
                f[int(round(v[0]))] = v[1]
    return f


def read_coords(path):
    c = {}
    for line in open(path):
        v = floats(line)
        if v and len(v) == 5:
            c[int(round(v[0]))] = (v[1], v[2], v[3])
    return c


def read_samples(path):
    s = {}
    for line in open(path):
        v = floats(line)
        if v and len(v) == 6:
            s.setdefault(int(round(v[0])), {})[int(round(v[1]))] = v[2:]
    return s


def transcript_ok(path):
    if not os.path.exists(path):
        return False, "missing"
    txt = open(path, errors="replace").read()
    done = "RUN COMPLETED" in txt
    nerr = None
    for line in txt.splitlines():
        if "NUMBER OF ERROR" in line and "ENCOUNTERED" in line:
            try:
                nerr = int(line.split("=")[-1])
            except ValueError:
                pass
    return (done and nerr == 0), "RUN COMPLETED=%s errors=%s" % (done, nerr)


def classify(freqs, coords, samples):
    top = {k: v for k, v in coords.items() if v[2] > 0}
    radii = sorted(set(round(v[0], 3) for v in top.values()))
    th = [15.0 * j for j in range(13)]
    out = []
    for mode in sorted(freqs):
        f = freqs[mode]
        smp = samples.get(mode, {})
        uz = [smp[k][2] for k in top if k in smp]
        ip = [math.hypot(smp[k][0], smp[k][1]) for k in top if k in smp]
        rms = lambda a: math.sqrt(sum(x * x for x in a) / max(len(a), 1))
        rec = dict(mode=mode, f=f, rigid=f < 1.0, flex=False, n=None, purity=None)
        if not uz:
            out.append(rec)
            continue
        rec["flex"] = (not rec["rigid"]) and rms(uz) > 2.0 * rms(ip)
        E = [0.0] * 13
        for r in radii:
            prof = {}
            for k, (rr, t, z) in top.items():
                if abs(rr - r) < 1e-3 and k in smp:
                    prof[int(round(t / 15.0))] = smp[k][2]
            if len(prof) != 13:
                continue
            for n in range(13):
                w = [0.5 if j in (0, 12) else 1.0 for j in range(13)]
                c = sum(w[j] * prof[j] * math.cos(n * math.radians(th[j])) for j in range(13))
                nrm = sum(w[j] * math.cos(n * math.radians(th[j])) ** 2 for j in range(13))
                E[n] += c * c / nrm
        tot = sum(E)
        if tot > 0:
            n = max(range(13), key=lambda i: E[i])
            rec["n"], rec["purity"] = n, E[n] / tot
        out.append(rec)
    by_n = {}
    for rec in out:
        if rec["flex"] and rec["n"] is not None:
            by_n.setdefault(rec["n"], []).append(rec)
    for n in by_n:
        for i, rec in enumerate(sorted(by_n[n], key=lambda x: x["f"])):
            rec["m"] = i + 1
    return out


def main():
    selftest = "--selftest" in sys.argv
    prefix = "ansys_p4_liu_disk_sc" if selftest else "ansys_p4_liudisk2_sc"
    tg = json.load(open(os.path.join(HERE, "targets_liudisk2_2026-09-24.json")))
    T = {(r["bc"], r["n"], r["m"]): r for r in tg["modes"]}
    fe = {}
    gates = {"F0": True, "F1": True, "F2": True, "F3": True}
    for bc in ("C", "S"):
        ok, why = transcript_ok(os.path.join(HERE, "%s_%s_2026-09-24_out.txt" % (prefix, bc)))
        print("[%s] transcript: %s -> %s" % (bc, why, "OK" if ok else "BAD"))
        gates["F0"] &= ok
        stem = os.path.join(HERE, "%s_%s_" % (prefix, bc))
        try:
            freqs = read_modes(stem + "modes.txt")
            recs = classify(freqs, read_coords(stem + "samplecoords.txt"), read_samples(stem + "samples.txt"))
        except IOError as e:
            print("[%s] NOT_RUN: %s" % (bc, e))
            gates["F0"] = False
            continue
        nrig = sum(1 for r in recs if r["rigid"])
        print("[%s] %d modes, rigid(<1 Hz)=%d" % (bc, len(recs), nrig))
        if nrig and not selftest:
            gates["F0"] = False
        for r in recs:
            if r["flex"] and r.get("m") and r["n"] <= 2 and r["m"] <= 2:
                fe[(bc, r["n"], r["m"])] = (r["f"], r["mode"], r["purity"])
    print()
    print("%-3s %s %s %5s %11s %11s %9s %11s %9s %11s %9s" % (
        "bc", "n", "m", "mode", "FE Hz", "Liu fem", "vs fem%", "CPT(cons)", "vs CPT%", "v1 FE", "vs v1%"))
    for key in sorted(T):
        t = T[key]
        if key not in fe:
            print("%-3s %d %d  MISSING" % key)
            for g in gates:
                gates[g] = False
            continue
        f, mode, pur = fe[key]
        d_fem = 100 * (f / t["f_liu_fem_Hz"] - 1)
        d_cpt = 100 * (f / t["f_cpt_consistent_Hz"] - 1)
        d_v1 = 100 * (f / t["f_v1_fe_2530987_Hz"] - 1)
        print("%-3s %d %d %5d %11.3f %11.3f %+9.3f %11.3f %+9.3f %11.3f %+9.3f  pur=%.2f" % (
            key[0], key[1], key[2], mode, f, t["f_liu_fem_Hz"], d_fem, t["f_cpt_consistent_Hz"], d_cpt,
            t["f_v1_fe_2530987_Hz"], d_v1, pur))
        if abs(d_fem) > 0.6:
            gates["F3"] = False
            if key[1] == 1:
                gates["F1"] = False
        if key[1] in (0, 2) and abs(d_v1) > 0.1:
            gates["F2"] = False
    print()
    for g in sorted(gates):
        print("%s %s  -- %s" % (g, "PASS" if gates[g] else "FAIL", tg["preregistered"][g]))
    verdict = "PASS" if all(gates.values()) else "FAIL"
    print("\nVERDICT %s%s" % (verdict, "  (SELFTEST on v1 outputs: F1/F3 are EXPECTED to fail)" if selftest else ""))


if __name__ == "__main__":
    main()
