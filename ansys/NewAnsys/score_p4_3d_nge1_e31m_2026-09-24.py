#!/usr/bin/env python3
"""e31 = -4.1 re-score (LESSONS Sec 18.245) of score_p4_3d_nge1_2026-09-23.py, for the
job 2532086 clones. Changes: targets from targets_p4_3d_nge1_e31m_2026-09-24.json
(same generator, e31 = -4.1); the G2(b) PLANE223 control splits are the -4.1
PLANE223 clones of the same job (121.12599/118.81638 at 1/12, 132.18292/127.10474
at 1/5); the coupling-off e0 decks are the unchanged parents. Gates and bars are
otherwise identical. Run on a directory holding the clone outputs renamed to the
parent stems (see LESSONS Sec 18.245).

Original docstring: Score the Paper 4 3-D FE check of coupled n >= 1 roots. Gates fixed before any run.

Reads each deck's transcript, _modes.txt, _samples.txt, _samplecoords.txt
and the pre-registered model targets (targets_p4_3d_nge1_2026-09-23.json).
Prints the raw Hz MAPDL wrote and does every ratio here. A missing file is
NOT_RUN. Nothing in this script invents or tunes a number.

Mode identification (fixed rule). Rigid: f < 1 Hz (expect exactly 3: X and
Z translation, rotation about Y). Flexural: RMS UZ over the 65 top-face
samples exceeds 2x the RMS in-plane magnitude there. Circumferential order:
UZ on each of the 5 sampled radii (theta = 0..180 step 15) is split into
cos(m theta), m = 0..12, by the exact 13-point DCT-I; energies are summed
over radii and the dominant m is the order. Purity = its energy share,
>= 0.90 required. Radial index s = rank by frequency among flexural modes
of that order (rigid modes excluded, as in the model scan which starts at
20 rad/s).

Pre-registered gates (fine mesh unless stated; gated set n = 0..4, s = 1
at both h1/2h = 1/12 and 1/5):

  G0  Every transcript clean (RUN COMPLETED, 0 MAPDL errors, no failure
      token); node/element counts equal the builder's; exactly 3 modes
      below 1 Hz; every gated (n, 1) identified with purity >= 0.90.
  G1  Mesh: |f_coarse/f_fine - 1| < 0.5% for e0 and sc, gated set.
  G2  Control against the axisymmetric PLANE223 FE (n = 0 only):
      (a) h1/2h=1/12 e0 within 0.10% of 118.81016 Hz (job 2490706);
      (b) OC/SC split within +-3% (relative) of the PLANE223 split,
          119.14314/118.81103 - 1 at 1/12 (2491761/2491831) and
          127.80931/127.04679 - 1 at 1/5 (2494226).
      This is what makes the n >= 1 electrical results mean anything.
  G3  Elastic n >= 1 operator, including the Kirchhoff effective shear
      (LESSONS Sec 18.226): e0 within 1.0% of the elastic bilayer 4x4
      for n = 1..4. (The retired Q_r-only row put 1/12 n=1 at
      1573.4926 rad/s = 250.43 Hz, 10.4% low; the 1% bar excludes it.)
  G4  No open-circuit split at n >= 1: |f_OC/f_SC - 1| <= 1e-6 for
      n = 1..4, and the outer-electrode bus voltage per unit peak |U|
      is below 1e-3 of the n = 0 value.
  G5  The coupled n >= 1 short-circuit shift (the actual "coupled root"):
      FE split f_SC/f_e0 - 1 against the model's SC 6x6 / elastic 4x4
      split, n = 0..4. PASS within +-10% (relative), PARTIAL within +-25%,
      else FAIL; sign must be positive. A point is UNRESOLVED (not scored)
      if its FE split moves more than 10% from coarse to fine.

Overall PASS needs G0-G5 all PASS. Report everything else (n = 5..7,
s = 2) as information only.

Usage (from Ansys/NewAnsys/ after the queue ran):
    python3 score_p4_3d_nge1_2026-09-23.py [directory]
    python3 score_p4_3d_nge1_2026-09-23.py --selftest
"""
# (no __future__ import: the compute nodes run python3 3.6, job 2530336)

import json
import math
import os
import re
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
DATE = "2026-09-23"
TARGETS = os.path.join(HERE, "targets_p4_3d_nge1_e31m_2026-09-24.json")

RATIOS = (("h112", 12), ("h15", 5))
CASES = ("e0", "sc", "oc")
DECKS = [("h112", "e0", "coarse"), ("h112", "sc", "coarse"),
         ("h112", "e0", "fine"), ("h112", "sc", "fine"), ("h112", "oc", "fine"),
         ("h15", "e0", "coarse"), ("h15", "sc", "coarse"),
         ("h15", "e0", "fine"), ("h15", "sc", "fine"), ("h15", "oc", "fine")]

MESH_COUNTS = {  # from build_p4_3d_nge1_decks_2026-09-23.py
    "coarse": dict(elements=6912, symmetry=690, face=5377),
    "fine": dict(elements=46080, symmetry=2114, face=17649),
}
NSAMP = 78
N_RIGID = 3
RIGID_HZ = 1.0
FLEX_RATIO = 2.0
PURITY = 0.90
GATED_N = (0, 1, 2, 3, 4)

MESH_BAR = 0.005
AXI_E0_H112 = 118.81016
AXI_E0_BAR = 0.0010
AXI_OC_SC = {"h112": 121.12599 / 118.81638 - 1.0,
             "h15": 132.18292 / 127.10474 - 1.0}
AXI_SPLIT_BAR = 0.03
ELASTIC_BAR = 0.010
QR_ONLY_N1_H112_HZ = 1573.4926 / (2 * math.pi)
OC_SC_BAR = 1e-6
VBUS_BAR = 1e-3
G5_PASS, G5_PARTIAL = 0.10, 0.25
SPLIT_MESH_BAR = 0.10

ERROR_TOKEN = "***" + " ERROR " + "***"
ERROR_COUNT = re.compile(r"NUMBER OF ERROR\s+MESSAGES ENCOUNTERED=\s+(-?\d+)")


def stem(tag, case, mesh):
    return "ansys_p4_3d_%s_%s_%s" % (tag, case, mesh)


def read_text(path):
    with open(path, "r", encoding="ascii", errors="replace", newline="") as f:
        return f.read().replace("\r\n", "\n").replace("\r", "\n")


def floats(line):
    out = []
    for tok in line.replace(",", " ").split():
        try:
            out.append(float(tok))
        except ValueError:
            return None
    return out


def transcript_problems(path, sentinel):
    reasons = []
    has_sent = os.path.isfile(sentinel)
    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        reasons.append("transcript missing" + (
            " (false sentinel: .QUEUE_OK exists)" if has_sent else ""))
        return reasons
    text = read_text(path)
    m = ERROR_COUNT.search(text)
    if "RUN COMPLETED" not in text:
        reasons.append("no RUN COMPLETED")
    if m is None:
        reasons.append("no MAPDL error-count line")
    elif int(m.group(1)) != 0:
        reasons.append("MAPDL error count %s" % m.group(1))
    if "ERROR: no input file specified" in text:
        reasons.append("MAPDL read no input file")
    if ERROR_TOKEN in text:
        reasons.append("transcript contains the MAPDL failure token")
    if reasons and has_sent:
        reasons.append("false sentinel: .QUEUE_OK set on an unclean transcript")
    return reasons


def parse_modes(path):
    counts, freqs = {}, []
    for line in read_text(path).splitlines():
        s = line.strip()
        for key in ("elements", "nodes", "symmetry_nodes", "interface_nodes",
                    "outer_nodes"):
            if s.startswith(key + " "):
                counts[key] = int(round(float(s.split()[-1])))
                break
        else:
            v = floats(s)
            if v is not None and len(v) == 2 and v[0] >= 1:
                freqs.append(v[1])
    return counts, freqs


def parse_coords(path):
    rows = {}
    for line in read_text(path).splitlines():
        v = floats(line.strip())
        if v is None or len(v) != 5:
            continue
        rows[int(round(v[0]))] = dict(r=v[1], th=v[2], z=v[3], node=int(round(v[4])))
    return rows


def parse_samples(path):
    by_mode = {}
    for line in read_text(path).splitlines():
        v = floats(line.strip())
        if v is None or len(v) != 6 or v[0] < 1:
            continue
        m, k = int(round(v[0])), int(round(v[1]))
        by_mode.setdefault(m, {})[k] = dict(ux=v[2], uy=v[3], uz=v[4], volt=v[5])
    return by_mode


def dct_energy(vals):
    """Exact DCT-I split of 13 samples at theta = 0..180 step 15 into
    cos(m theta), m = 0..12. Returns per-m energy (weighted Parseval)."""
    N = len(vals) - 1
    w = [0.5 if k in (0, N) else 1.0 for k in range(N + 1)]
    th = [math.pi * k / N for k in range(N + 1)]
    out = []
    for m in range(N + 1):
        c = sum(w[k] * vals[k] * math.cos(m * th[k]) for k in range(N + 1))
        norm = sum(w[k] * math.cos(m * th[k]) ** 2 for k in range(N + 1))
        out.append(c * c / norm)
    return out


def classify(freqs, samples, coords):
    """Per mode: dict(f, kind, n, purity, vbus_rel)."""
    top = [k for k, c in coords.items() if c["z"] > 0]
    radii = sorted({round(coords[k]["r"], 3) for k in top})
    out = []
    for i, f in enumerate(freqs, start=1):
        smp = samples.get(i, {})
        rec = dict(mode=i, f=f, kind=None, n=None, purity=None, vbus_rel=None)
        if len(smp) != len(coords):
            rec["kind"] = "incomplete"
            out.append(rec)
            continue
        uz = [smp[k]["uz"] for k in top]
        uin = [math.hypot(smp[k]["ux"], smp[k]["uy"]) for k in top]
        rms_z = math.sqrt(sum(u * u for u in uz) / len(uz))
        rms_in = math.sqrt(sum(u * u for u in uin) / len(uin))
        peak = max(max(abs(u) for u in uz), max(uin)) or 1.0
        rec["vbus_rel"] = max(abs(smp[k]["volt"]) for k in top) / peak
        if f < RIGID_HZ:
            rec["kind"] = "rigid"
        elif rms_z > FLEX_RATIO * rms_in:
            rec["kind"] = "flex"
        else:
            rec["kind"] = "inplane"
        energy = [0.0] * 13
        for r in radii:
            ks = sorted((k for k in top if round(coords[k]["r"], 3) == r),
                        key=lambda k: coords[k]["th"])
            e = dct_energy([smp[k]["uz"] for k in ks])
            energy = [a + b for a, b in zip(energy, e)]
        tot = sum(energy)
        if tot > 0:
            m = max(range(13), key=lambda j: energy[j])
            rec["n"], rec["purity"] = m, energy[m] / tot
        out.append(rec)
    return out


def assign(modes):
    """(n, s) -> mode record for flexural modes with purity >= PURITY."""
    by_n = {}
    for rec in modes:
        if rec["kind"] == "flex" and rec["purity"] is not None:
            by_n.setdefault(rec["n"], []).append(rec)
    table = {}
    for n, recs in by_n.items():
        for s, rec in enumerate(sorted(recs, key=lambda r: r["f"]), start=1):
            table[(n, s)] = rec
    return table


def load_deck(directory, tag, case, mesh):
    st = stem(tag, case, mesh)
    d = dict(tag=tag, case=case, mesh=mesh, reasons=[], freqs=[], modes=[],
             table={}, counts={})
    d["reasons"] += transcript_problems(
        os.path.join(directory, "%s_%s_out.txt" % (st, DATE)),
        os.path.join(directory, "%s_%s.inp.QUEUE_OK" % (st, DATE)))
    paths = [os.path.join(directory, st + suf) for suf in
             ("_modes.txt", "_samples.txt", "_samplecoords.txt")]
    missing = [os.path.basename(p) for p in paths if not os.path.isfile(p)]
    if missing:
        d["reasons"].append("missing " + ", ".join(missing))
        d["not_run"] = True
        return d
    counts, freqs = parse_modes(paths[0])
    samples = parse_samples(paths[1])
    coords = parse_coords(paths[2])
    d["counts"], d["freqs"] = counts, freqs
    exp = MESH_COUNTS[mesh]
    for key, want in (("elements", exp["elements"]),
                      ("symmetry_nodes", exp["symmetry"]),
                      ("interface_nodes", 2 * exp["face"]),
                      ("outer_nodes", 2 * exp["face"])):
        if counts.get(key) != want:
            d["reasons"].append("%s = %s, expected %d" % (key, counts.get(key), want))
    if len(coords) != NSAMP:
        d["reasons"].append("%d sample coords, expected %d" % (len(coords), NSAMP))
    d["modes"] = classify(freqs, samples, coords)
    n_rigid = sum(1 for m in d["modes"] if m["kind"] == "rigid")
    if n_rigid != N_RIGID:
        d["reasons"].append("%d modes below %.1f Hz, expected %d"
                            % (n_rigid, RIGID_HZ, N_RIGID))
    d["table"] = assign(d["modes"])
    for n in GATED_N:
        rec = d["table"].get((n, 1))
        if rec is None:
            d["reasons"].append("gated (n=%d, s=1) not identified" % n)
        elif rec["purity"] < PURITY:
            d["reasons"].append("(n=%d, s=1) purity %.3f < %.2f"
                                % (n, rec["purity"], PURITY))
    return d


def rel(a, b):
    return a / b - 1.0


def verdict(ok):
    return "PASS" if ok else "FAIL"


def score(directory, targets_path=TARGETS):
    with open(targets_path) as f:
        T = json.load(f)
    decks = {(t, c, m): load_deck(directory, t, c, m) for t, c, m in DECKS}
    L = []
    a = L.append
    a("Paper 4 3-D FE n>=1 score  (%s)" % directory)
    a("targets: %s" % os.path.basename(targets_path))
    a("")
    any_missing = False
    g0 = True
    for key in DECKS:
        d = decks[key]
        a("%-18s raw Hz: %s" % ("_".join(key), ", ".join("%.10g" % x for x in d["freqs"]) or "(none)"))
        for rec in d["modes"]:
            a("    mode %2d  %14.8f Hz  %-8s n=%-3s purity=%s  Vbus/|U|=%s"
              % (rec["mode"], rec["f"], rec["kind"], rec["n"],
                 "%.4f" % rec["purity"] if rec["purity"] is not None else "-",
                 "%.3e" % rec["vbus_rel"] if rec["vbus_rel"] is not None else "-"))
        for r in d["reasons"]:
            a("    G0: %s" % r)
        if d.get("not_run"):
            any_missing = True
        if d["reasons"]:
            g0 = False
    a("")
    if any_missing:
        a("G0: NOT_RUN (a deck's outputs are missing). G1-G5 not scored.")
        a("RESULT: NOT_RUN")
        return "\n".join(L), "NOT_RUN"
    a("G0: %s" % verdict(g0))

    def f(tag, case, mesh, n, s=1):
        rec = decks[(tag, case, mesh)]["table"].get((n, s))
        return None if rec is None else rec["f"]

    def tgt(tag, n, s=1):
        return T[tag].get("%d_%d" % (n, s))

    # G1
    g1 = True
    a("")
    a("G1 mesh (|coarse/fine - 1| < %.1f%%):" % (100 * MESH_BAR))
    for tag, _ in RATIOS:
        for case in ("e0", "sc"):
            for n in GATED_N:
                fc, ff = f(tag, case, "coarse", n), f(tag, case, "fine", n)
                if fc is None or ff is None:
                    a("    %s %s n=%d: missing" % (tag, case, n)); g1 = False; continue
                x = rel(fc, ff)
                ok = abs(x) < MESH_BAR
                g1 &= ok
                a("    %s %s n=%d: coarse %.8f fine %.8f  %+.4f%%  %s"
                  % (tag, case, n, fc, ff, 100 * x, verdict(ok)))
    a("G1: %s" % verdict(g1))

    # G2
    a("")
    g2 = True
    fe0 = f("h112", "e0", "fine", 0)
    if fe0 is None:
        g2 = False; a("G2a: n=0 e0 missing")
    else:
        x = rel(fe0, AXI_E0_H112)
        ok = abs(x) <= AXI_E0_BAR
        g2 &= ok
        a("G2a h112 e0 n=0: 3-D %.8f Hz vs PLANE223 %.5f Hz  %+.4f%%  %s"
          % (fe0, AXI_E0_H112, 100 * x, verdict(ok)))
    for tag, _ in RATIOS:
        fo, fs = f(tag, "oc", "fine", 0), f(tag, "sc", "fine", 0)
        if fo is None or fs is None:
            g2 = False; a("G2b %s: missing" % tag); continue
        sp = rel(fo, fs)
        x = sp / AXI_OC_SC[tag] - 1.0
        ok = abs(x) <= AXI_SPLIT_BAR
        g2 &= ok
        a("G2b %s n=0 OC/SC: 3-D %.5f%% vs PLANE223 %.5f%%  (3-D/axi - 1 = %+.2f%%)  "
          "model oc/sc %.5f%%  %s"
          % (tag, 100 * sp, 100 * AXI_OC_SC[tag], 100 * x,
             100 * tgt(tag, 0)["oc_sc_split"], verdict(ok)))
    a("G2: %s" % verdict(g2))

    # G3
    a("")
    g3 = True
    a("G3 elastic bilayer n>=1 (|FE e0/model - 1| <= %.1f%%); n=0 shown for the offset:"
      % (100 * ELASTIC_BAR))
    for tag, _ in RATIOS:
        for n in GATED_N:
            fe, t = f(tag, "e0", "fine", n), tgt(tag, n)
            if fe is None or t is None:
                a("    %s n=%d: missing" % (tag, n)); g3 &= n == 0; continue
            x = rel(fe, t["el_hz"])
            ok = abs(x) <= ELASTIC_BAR
            if n >= 1:
                g3 &= ok
            a("    %s n=%d s=1: FE %.6f Hz  model %.6f Hz  %+.3f%%  %s"
              % (tag, n, fe, t["el_hz"], 100 * x,
                 verdict(ok) if n >= 1 else "(offset reference)"))
    fe1 = f("h112", "e0", "fine", 1)
    if fe1 is not None:
        a("    retired Q_r-only row, h112 n=1: %.3f Hz -> FE/old - 1 = %+.2f%%"
          % (QR_ONLY_N1_H112_HZ, 100 * rel(fe1, QR_ONLY_N1_H112_HZ)))
    a("G3: %s" % verdict(g3))

    # G4
    a("")
    g4 = True
    a("G4 no OC split at n>=1 (|OC/SC - 1| <= %.0e; Vbus ratio < %.0e of n=0):"
      % (OC_SC_BAR, VBUS_BAR))
    for tag, _ in RATIOS:
        tab = decks[(tag, "oc", "fine")]["table"]
        v0 = tab.get((0, 1), {}).get("vbus_rel") if tab.get((0, 1)) else None
        for n in GATED_N[1:]:
            fo, fs = f(tag, "oc", "fine", n), f(tag, "sc", "fine", n)
            if fo is None or fs is None or v0 in (None, 0.0):
                a("    %s n=%d: missing" % (tag, n)); g4 = False; continue
            x = rel(fo, fs)
            vr = tab[(n, 1)]["vbus_rel"] / v0
            ok = abs(x) <= OC_SC_BAR and vr < VBUS_BAR
            g4 &= ok
            a("    %s n=%d: OC %.10f SC %.10f  OC/SC-1 = %+.3e  Vbus(n)/Vbus(0) = %.3e  %s"
              % (tag, n, fo, fs, x, vr, verdict(ok)))
    a("G4: %s" % verdict(g4))

    # G5
    a("")
    g5_states = []
    a("G5 SC interior-potential shift, FE sc/e0 - 1 vs model SC6x6/el4x4 - 1"
      " (PASS <= %d%%, PARTIAL <= %d%%):" % (100 * G5_PASS, 100 * G5_PARTIAL))
    for tag, _ in RATIOS:
        for n in GATED_N:
            fs, fe = f(tag, "sc", "fine", n), f(tag, "e0", "fine", n)
            cs, ce = f(tag, "sc", "coarse", n), f(tag, "e0", "coarse", n)
            t = tgt(tag, n)
            if None in (fs, fe, cs, ce) or t is None:
                a("    %s n=%d: missing" % (tag, n)); g5_states.append("FAIL"); continue
            sp_f, sp_c, sp_m = rel(fs, fe), rel(cs, ce), t["sc_el_split"]
            mesh_move = abs(sp_c / sp_f - 1.0) if sp_f != 0 else float("inf")
            x = sp_f / sp_m - 1.0
            if mesh_move > SPLIT_MESH_BAR:
                state = "UNRESOLVED"
            elif sp_f <= 0:
                state = "FAIL"
            elif abs(x) <= G5_PASS:
                state = "PASS"
            elif abs(x) <= G5_PARTIAL:
                state = "PARTIAL"
            else:
                state = "FAIL"
            g5_states.append(state)
            a("    %s n=%d: FE %.4e (coarse %.4e, move %.1f%%)  model %.4e  "
              "FE/model-1 = %+.1f%%  %s"
              % (tag, n, sp_f, sp_c, 100 * mesh_move, sp_m, 100 * x, state))
    if "FAIL" in g5_states:
        g5 = "FAIL"
    elif "PARTIAL" in g5_states or "UNRESOLVED" in g5_states:
        g5 = "PARTIAL" if "PARTIAL" in g5_states else "UNRESOLVED"
    else:
        g5 = "PASS"
    a("G5: %s  (%s)" % (g5, ", ".join(g5_states)))

    # information only
    a("")
    a("Information only (not gated): other identified flexural modes, fine mesh")
    for tag, _ in RATIOS:
        tab_e0 = decks[(tag, "e0", "fine")]["table"]
        tab_sc = decks[(tag, "sc", "fine")]["table"]
        for key in sorted(T[tag]):
            n, s = (int(x) for x in key.split("_"))
            if s == 1 and n in GATED_N:
                continue
            t = T[tag][key]
            e0 = tab_e0.get((n, s)); sc = tab_sc.get((n, s))
            if e0 is None:
                a("    %s n=%d s=%d model %.3f Hz: not found in FE" % (tag, n, s, t["el_hz"]))
                continue
            line = ("    %s n=%d s=%d: FE e0 %.4f Hz vs model %.4f  %+.3f%%"
                    % (tag, n, s, e0["f"], t["el_hz"], 100 * rel(e0["f"], t["el_hz"])))
            if sc is not None:
                line += ("   split FE %.3e model %.3e"
                         % (rel(sc["f"], e0["f"]), t["sc_el_split"]))
            a(line)

    overall = "PASS" if (g0 and g1 and g2 and g3 and g4 and g5 == "PASS") else "FAIL"
    if overall == "FAIL" and g0 and g1 and g2 and g3 and g4 and g5 in ("PARTIAL", "UNRESOLVED"):
        overall = "PARTIAL"
    a("")
    a("RESULT: %s  (G0 %s, G1 %s, G2 %s, G3 %s, G4 %s, G5 %s)"
      % (overall, verdict(g0), verdict(g1), verdict(g2), verdict(g3), verdict(g4), g5))
    return "\n".join(L), overall


# ---------------------------------------------------------------- selftest
def _selftest():
    """Write synthetic outputs from the targets (cos(n theta) shapes with
    small FE-like offsets) and check the scorer passes them, then that it
    fails a deliberately broken OC deck. Exercises parsing and the rules,
    not physics."""
    with open(TARGETS) as fh:
        T = json.load(fh)
    tmp = tempfile.mkdtemp(prefix="p4_3d_selftest_")
    try:
        radii = (100.0, 225.0, 350.0, 475.0, 600.0)
        thetas = [15.0 * k for k in range(13)]
        coords = [(r, th, 1.0) for r in radii for th in thetas]
        coords += [(600.0, th, -1.0) for th in thetas]

        def write(tag, case, mesh, broken=False):
            st = stem(tag, case, mesh)
            mc = MESH_COUNTS[mesh]
            modes = [(0.0001, "rigid", 0), (0.0002, "rigid", 1), (0.0003, "rigid", 1)]
            for key, t in T[tag].items():
                if key == "meta":
                    continue
                n, s = (int(x) for x in key.split("_"))
                f = t["el_hz"] * (1 - 0.002 - 0.0005 * n) * (1.0005 if mesh == "coarse" else 1.0)
                if case in ("sc", "oc"):
                    f *= 1 + t["sc_el_split"] * (1.03 if mesh == "fine" else 1.05)
                if case == "oc" and n == 0:
                    f *= 1 + AXI_OC_SC[tag]
                if broken and case == "oc" and n == 1:
                    f *= 1 + 1e-4
                modes.append((f, "flex", n, s))
            modes.sort(key=lambda x: x[0])
            modes = modes[:24]
            with open(os.path.join(tmp, "%s_%s_out.txt" % (st, DATE)), "w") as fh:
                fh.write(" NUMBER OF ERROR   MESSAGES ENCOUNTERED=          0\n RUN COMPLETED\n")
            with open(os.path.join(tmp, st + "_modes.txt"), "w") as fh:
                fh.write("tag x\nelements %12.0f\nnodes %12.0f\nsymmetry_nodes %12.0f\n"
                         "interface_nodes %12.0f\nouter_nodes %12.0f\nmode f_Hz\n"
                         % (mc["elements"], 1, mc["symmetry"], 2 * mc["face"], 2 * mc["face"]))
                for i, m in enumerate(modes, start=1):
                    fh.write("%4.0f   %20.12E\n" % (i, m[0]))
            with open(os.path.join(tmp, st + "_samplecoords.txt"), "w") as fh:
                fh.write("k r_mm theta_deg z_mm node\n")
                for k, (r, th, z) in enumerate(coords, start=1):
                    fh.write("%5.0f   %10.3f   %10.3f   %10.4f   %10.0f\n" % (k, r, th, z, k))
            with open(os.path.join(tmp, st + "_samples.txt"), "w") as fh:
                fh.write("mode k UX UY UZ VOLT\n")
                for i, m in enumerate(modes, start=1):
                    n = m[2]
                    for k, (r, th, z) in enumerate(coords, start=1):
                        c = math.cos(n * math.radians(th))
                        if m[1] == "rigid":
                            ux, uy, uz = (1.0, 0.0, 0.0) if i == 1 else (0.0, 0.0, 1.0)
                        else:
                            uz = c * (0.3 + r / 600.0)
                            ux = 0.05 * c * z
                            uy = 0.02 * z * math.sin(n * math.radians(th))
                        v = 0.0
                        if case == "oc" and z > 0:
                            v = 0.5 if n == 0 else 1e-9
                        fh.write("%4.0f %5.0f %16.8E %16.8E %16.8E %16.8E\n"
                                 % (i, k, ux, uy, uz, v))
        for t, c, m in DECKS:
            write(t, c, m)
        text, res = score(tmp)
        assert res == "PASS", text
        write("h112", "oc", "fine", broken=True)
        text2, res2 = score(tmp)
        assert res2 == "FAIL" and "G4: FAIL" in text2, text2
        print("selftest OK: clean synthetic set PASS, broken OC n=1 -> G4 FAIL")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def main(argv):
    if len(argv) > 1 and argv[1] == "--selftest":
        _selftest()
        return
    directory = argv[1] if len(argv) > 1 else HERE
    text, _ = score(directory)
    print(text)
    out = os.path.join(directory, "score_p4_3d_nge1_e31m_2026-09-24.txt")
    with open(out, "w", newline="\n") as fh:
        fh.write(text + "\n")
    print("\nwrote", out)


if __name__ == "__main__":
    main(sys.argv)
