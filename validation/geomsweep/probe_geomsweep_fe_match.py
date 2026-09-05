# -*- coding: utf-8 -*-
"""
probe_geomsweep_fe_match.py -- RECONSTRUCTION, not a re-run of an existing
script. DIAGNOSTIC/DATA-GENERATION ONLY: reads the checked-in manifests and
ANSYS text results, writes a plain-text log; no package changes, no
SOLVER_VERSION implications. Standalone -- does not import plate_solver.

CONTEXT: default GEOMSWEEP_MANIFEST_DIR is the live Sec 6.4 recapture
(job 2416746): unpack `_s10v2_manifests.zip` and report 1256/1475
(85.2%, mean 0.574%, worst 2.990%). The older figures_ff_r*/ tree in
this directory is the superseded 1231/1474 (83.5%) recapture (job
2328808); set GEOMSWEEP_MANIFEST_DIR to this directory to replay it.
The original 1232/1466 (84.0%) interactive session was never saved;
this script is a from-scratch reconstruction of that analysis from
checked-in ingredients:
  - python candidates: the s10v2 zip (default) or
    validation/geomsweep/figures_ff_r{150,167,200,250}/freefree_manifest.json
  - independent FE ground truth: ansys/geomsweep/geomsweep_{oop,ip}_r*_a*.txt
    (the per-angle *CFOPEN/*VWRITE tables each ANSYS deck writes directly;
    NOT the raw MAPDL solver logs -- see ansys/geomsweep/README.md)

WHAT IS AND ISN'T RECOVERABLE: Sec 6.3's own methodology text says
extensional-family SHELL281 modes are "removed by frequency-matching
against the companion in-plane model" but gives NO numeric tolerance for
that cross-match anywhere in the paper, LESSONS_LEARNED.md, or any script
in either copy of the working directory (confirmed by search, 2026-07-23).
That tolerance was therefore a human judgment call made once, in an
interactive session, and never recorded -- it is NOT recoverable, only
re-decidable. This script makes an explicit, documented choice (a 1.0%
relative Omega_lit tolerance, see CONTAMINATION_TOL_PCT below) and reports
sensitivity to that choice (0.5% / 1.0% / 2.0%) rather than presenting one
number as if it were the original. Do not read this script's output as
proof the original interactive session used the same tolerance -- it
almost certainly did not use exactly this one. What this script DOES
establish is a transparent, rerunnable procedure that any future session
can adjust, rather than an undocumented one-off.

METHOD (mirrors Sec 6.4's own prose description):
  1. Rigid-body floor: FE modes with Omega_lit < RIGID_FLOOR are discarded
     on both the OOP and IP sides (the FE logs' own header comments note
     ~6 near-zero rigid OOP modes, ~3 for IP).
  2. IP Omega_lit is not printed by the deck (its header says "no family
     separation needed" -- PLANE183 has no bending family to interleave).
     It's recovered from f[Hz] using the same linear Omega_lit = C * f_hz
     relationship the OOP deck's own printed (f_hz, Omega_lit) pairs
     already exhibit (Omega_lit/f_hz is constant per radius ratio, since
     it depends only on R_o and the material, not on part or angle) -- C
     is fit per radius ratio from the OOP file's own data, not assumed.
  3. Extensional-family contamination removal (OOP only, see note above):
     an OOP FE mode within CONTAMINATION_TOL_PCT of an IP FE mode at the
     SAME (ratio, angle) is dropped from the OOP list before matching.
  4. Each geometry's first n_cutoffs python candidates (n_cutoffs =
     len(diag.cut_offs) from that geometry's own manifest entry -- 5 for
     every OOP entry, 4 for every IP entry, confirmed uniform across all
     48 geometries) are set aside, exactly as Sec 6.4 describes.
  5. Remaining candidates are nearest-neighbor matched (by Omega_lit)
     against the cleaned FE list for that (ratio, angle, part), restricted
     to the FE-resolved range (candidates above the FE list's own max
     Omega_lit cannot be fairly matched and are excluded from n, not
     counted as failures).
  6. Per-candidate error % = 100 * |cand - FE_match| / FE_match. Tallied:
     n, n matched <=3%, mean of per-geometry mean errors, worst case.

Only the 4 NEW radius ratios (1.5, 1.66667, 2.0, 2.5) x 6 angles x 2 parts
= 48 geometries are included, matching Sec 6.4's stated scope --
figures_ff_ip_anglefill_r125 (r0/2b=1.25) is the ORIGINAL Sec 6.3
validated geometry filling in previously-unrun IP angles, a different
purpose, and is deliberately excluded here.
"""
import csv
import glob
import json
import os
import re
import tempfile
import zipfile

HERE = os.path.dirname(os.path.abspath(__file__))
S10V2_ZIP = os.path.join(HERE, "_s10v2_manifests.zip")


def resolve_manifest_dir():
    """Live paper default: job 2416746 zip. Env or missing zip -> HERE."""
    env = os.environ.get("GEOMSWEEP_MANIFEST_DIR")
    if env:
        return env, "env"
    if os.path.isfile(S10V2_ZIP):
        dest = tempfile.mkdtemp(prefix="s10v2_manifests_")
        with zipfile.ZipFile(S10V2_ZIP) as zf:
            zf.extractall(dest)
        return dest, "s10v2_zip"
    return HERE, "figures_ff_fallback"


# Overridable via env vars so this same script runs unmodified against
# either layout: the curated github_repo/ copy (validation/geomsweep/ +
# ansys/geomsweep/, the defaults below) or the canonical cluster working
# tree (FutureWork/geometry_sweep/ + FutureWork/ansys_geomsweep_verification/,
# where the manifests and FE files actually live -- see
# submit_geomsweep_fe_match.sh, which sets both for the cluster layout).
MANIFEST_DIR, MANIFEST_SOURCE = resolve_manifest_dir()
FE_DIR = os.environ.get(
    "GEOMSWEEP_FE_DIR",
    os.path.normpath(os.path.join(HERE, "..", "..", "ansys", "geomsweep")))
MANIFEST_GLOB = os.path.join(MANIFEST_DIR, "**", "figures_ff_r*", "freefree_manifest.json")
OUT_CSV = os.environ.get("GEOMSWEEP_MATCH_CSV", "geomsweep_fe_match_consolidated.csv")

RIGID_FLOOR = 1.0          # Omega_lit below this is treated as rigid-body
CONTAMINATION_TOL_PCT = float(os.environ.get("GEOMSWEEP_CONTAM_TOL", "1.0"))

RATIO_CODE = {1.5: "r150", 1.66667: "r167", 2.0: "r200", 2.5: "r250"}
ANGLE_CODE = {0.25: "a025", 0.5: "a050", 0.75: "a075",
              1.0: "a100", 1.25: "a125", 1.5: "a150"}

MODE_LINE_OOP = re.compile(
    r"^\s*(\d+)\.\s+([\d.Ee+\-]+)\s+([\d.Ee+\-]+)\s*$")
MODE_LINE_IP = re.compile(r"^\s*(\d+)\.\s+([\d.Ee+\-]+)\s*$")


def load_manifests():
    """Return {(r0_2b, two_T_pi, part): {'omega_lit': [...], 'n_cutoffs': N}}
    for the 4 new-ratio manifests only (figures_ff_ip_anglefill_r125
    excluded -- see module docstring)."""
    out = {}
    for path in sorted(glob.glob(MANIFEST_GLOB, recursive=True)):
        if "anglefill" in path:
            continue
        with open(path) as f:
            data = json.load(f)
        for key, v in data.items():
            if not isinstance(v, dict) or "diag" not in v:
                continue  # e.g. "_cross_geometry_suspects" metadata entries
            parts = key.split()
            if len(parts) != 3:
                continue
            r0_2b = float(parts[1])
            two_T_pi = float(parts[2])
            part = v["diag"]["part"]
            n_cutoffs = len(v["diag"]["cut_offs"])
            omega_lit = list(v["omega_lit"])
            out[(r0_2b, two_T_pi, part)] = {
                "omega_lit": omega_lit,
                "n_cutoffs": n_cutoffs,
            }
    return out


def parse_fe_oop(path):
    """Return list of (mode_idx, f_hz, omega_lit) from an OOP FE file."""
    rows = []
    with open(path) as f:
        for line in f:
            m = MODE_LINE_OOP.match(line)
            if m:
                rows.append((int(m.group(1)), float(m.group(2)), float(m.group(3))))
    return rows


def parse_fe_ip_fhz(path):
    """Return list of (mode_idx, f_hz) from an IP FE file (no Omega_lit
    column printed -- see module docstring)."""
    rows = []
    with open(path) as f:
        for line in f:
            m = MODE_LINE_IP.match(line)
            if m:
                rows.append((int(m.group(1)), float(m.group(2))))
    return rows


def fit_omega_lit_per_fhz(r0_2b):
    """Fit C = Omega_lit / f_hz for this radius ratio from every non-rigid
    row of every angle's OOP file (should be constant; report scatter)."""
    rtag = RATIO_CODE[r0_2b]
    ratios = []
    for atag in ANGLE_CODE.values():
        path = os.path.join(FE_DIR, f"geomsweep_oop_{rtag}_{atag}.txt")
        for _, f_hz, om in parse_fe_oop(path):
            if om >= RIGID_FLOOR and f_hz > 0:
                ratios.append(om / f_hz)
    mean_c = sum(ratios) / len(ratios)
    spread = (max(ratios) - min(ratios)) / mean_c if ratios else float("nan")
    return mean_c, spread, len(ratios)


def load_fe_lists(r0_2b, two_T_pi, c_per_fhz, tol_pct):
    """Return (oop_list, ip_list): cleaned Omega_lit lists (rigid removed,
    OOP extensional-contamination removed) for one (ratio, angle)."""
    rtag = RATIO_CODE[r0_2b]
    atag = ANGLE_CODE[two_T_pi]
    oop_raw = [om for _, f_hz, om in
               parse_fe_oop(os.path.join(FE_DIR, f"geomsweep_oop_{rtag}_{atag}.txt"))
               if om >= RIGID_FLOOR]
    ip_raw = [f_hz * c_per_fhz for _, f_hz in
              parse_fe_ip_fhz(os.path.join(FE_DIR, f"geomsweep_ip_{rtag}_{atag}.txt"))]
    ip_raw = [om for om in ip_raw if om >= RIGID_FLOOR]

    oop_clean = []
    for om in oop_raw:
        contaminated = any(abs(om - ip_om) / ip_om * 100.0 <= tol_pct
                            for ip_om in ip_raw)
        if not contaminated:
            oop_clean.append(om)
    return oop_clean, ip_raw


def nearest_match(cand, fe_list):
    if not fe_list:
        return None, None
    best = min(fe_list, key=lambda x: abs(x - cand))
    err_pct = 100.0 * abs(cand - best) / best
    return best, err_pct


def run(tol_pct, verbose=False, collect_rows=None):
    manifests = load_manifests()
    c_cache = {}
    per_geom_stats = []
    all_errs = []
    n_matched_3pct = 0
    n_total = 0

    for (r0_2b, two_T_pi, part), entry in sorted(manifests.items()):
        if r0_2b not in c_cache:
            c_cache[r0_2b] = fit_omega_lit_per_fhz(r0_2b)
        c_per_fhz, spread, nfit = c_cache[r0_2b]

        oop_list, ip_list = load_fe_lists(r0_2b, two_T_pi, c_per_fhz, tol_pct)
        fe_list = oop_list if part == 1 else ip_list

        all_cands = entry["omega_lit"]
        n_cutoffs = entry["n_cutoffs"]
        fe_max = max(fe_list) if fe_list else float("-inf")

        errs = []
        for idx, cand in enumerate(all_cands):
            cutoff_adjacent = idx < n_cutoffs
            beyond_fe_range = (not cutoff_adjacent) and cand > fe_max
            row = {
                "r0_2b": r0_2b, "two_T_pi": two_T_pi, "part": part,
                "cand_index": idx, "cand_omega_lit": cand,
                "cutoff_adjacent": cutoff_adjacent,
                "beyond_fe_range": beyond_fe_range,
                "fe_match_omega_lit": "", "err_pct": "",
                "matched_within_3pct": "",
            }
            if not cutoff_adjacent and not beyond_fe_range and fe_list:
                fe_match, err_pct = nearest_match(cand, fe_list)
                row["fe_match_omega_lit"] = fe_match
                row["err_pct"] = err_pct
                row["matched_within_3pct"] = err_pct <= 3.0
                errs.append(err_pct)
                all_errs.append(err_pct)
                n_total += 1
                if err_pct <= 3.0:
                    n_matched_3pct += 1
            if collect_rows is not None:
                collect_rows.append(row)

        # Table geomsweep's own mean/max err columns are computed over the
        # MATCHED (<=3%) subset only, not all attempted candidates -- the
        # per-ratio "n" values in that table (276/290/327/339) sum to
        # exactly 1232, the matched numerator, not 1466. Confirmed by
        # reconciling this reconstruction against the printed table.
        matched_errs = [e for e in errs if e <= 3.0]
        if matched_errs:
            per_geom_stats.append({
                "r0_2b": r0_2b, "two_T_pi": two_T_pi, "part": part,
                "n": len(matched_errs),
                "mean": sum(matched_errs) / len(matched_errs),
                "max": max(matched_errs),
            })
        if verbose:
            print(f"  r0_2b={r0_2b:g} 2T/pi={two_T_pi:g} part={part} "
                  f"n_cand={len(entry['omega_lit'])} n_cutoffs={entry['n_cutoffs']} "
                  f"fe_n={len(fe_list)} fe_max={fe_max:.1f} n_matched={len(errs)}")

    mean_of_means = (sum(g["mean"] for g in per_geom_stats) / len(per_geom_stats)
                      if per_geom_stats else float("nan"))
    worst = max((g["max"] for g in per_geom_stats), default=float("nan"))
    pct_matched = 100.0 * n_matched_3pct / n_total if n_total else float("nan")

    return {
        "tol_pct": tol_pct,
        "n_total": n_total,
        "n_matched_3pct": n_matched_3pct,
        "pct_matched": pct_matched,
        "mean_of_means": mean_of_means,
        "worst": worst,
        "n_geoms": len(per_geom_stats),
        "per_geom": per_geom_stats,
    }


def write_csv(rows, path):
    fieldnames = ["r0_2b", "two_T_pi", "part", "cand_index", "cand_omega_lit",
                  "cutoff_adjacent", "beyond_fe_range", "fe_match_omega_lit",
                  "err_pct", "matched_within_3pct"]
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for row in rows:
            w.writerow(row)


def main():
    print(__doc__.split("METHOD")[0].strip()[:600] + "\n...\n")
    print(f"MANIFEST_DIR    = {MANIFEST_DIR}")
    print(f"MANIFEST_SOURCE = {MANIFEST_SOURCE}")
    print(f"FE_DIR          = {FE_DIR}")
    print("=" * 78)
    print("LIVE PAPER (Sec 6.4 / Table S.3): 1256/1475 (85.2%), mean "
          "0.574%, worst 2.990%  -- job 2416746")
    if MANIFEST_SOURCE == "s10v2_zip":
        print("THIS RUN: default _s10v2_manifests.zip (job 2416746).")
    elif MANIFEST_SOURCE == "env":
        print("THIS RUN: GEOMSWEEP_MANIFEST_DIR override.")
    else:
        print("THIS RUN: figures_ff_* fallback (superseded 1231/1474).")
    print("=" * 78)

    primary_rows = []
    r1 = None
    for tol in (0.5, 1.0, 2.0):
        collect = primary_rows if tol == CONTAMINATION_TOL_PCT else None
        r = run(tol, verbose=(tol == CONTAMINATION_TOL_PCT), collect_rows=collect)
        if tol == CONTAMINATION_TOL_PCT:
            r1 = r
        print(f"\n--- contamination tolerance = {tol}% ---")
        print(f"  n_total (attempted) = {r['n_total']}  "
              f"(live paper: 1475)")
        print(f"  n matched <=3%      = {r['n_matched_3pct']} "
              f"({r['pct_matched']:.1f}%)  (live paper: 1256, 85.2%)")
        print(f"  mean of per-geom means = {r['mean_of_means']:.3f}%  "
              f"(live paper: 0.574%)")
        print(f"  worst case = {r['worst']:.3f}%  (live paper: 2.990%)")
        print(f"  n_geometries with >=1 matched candidate = {r['n_geoms']} "
              f"(expected: 48)")

    write_csv(primary_rows, OUT_CSV)
    print(f"\nWrote full per-candidate audit trail (tolerance="
          f"{CONTAMINATION_TOL_PCT}%, {len(primary_rows)} rows, every "
          f"candidate from every geometry -- cutoff-adjacent and "
          f"beyond-FE-range candidates included and flagged, not silently "
          f"dropped) to {OUT_CSV}")

    print(f"\nRESULT (source={MANIFEST_SOURCE}, tol={CONTAMINATION_TOL_PCT}%): "
          f"{r1['n_matched_3pct']}/{r1['n_total']} "
          f"({r1['pct_matched']:.1f}%), mean {r1['mean_of_means']:.3f}%, "
          f"worst {r1['worst']:.3f}%. Live paper: 1256/1475 (85.2%), "
          f"mean 0.574%, worst 2.990%.")
    print("\nSee module docstring: the extensional-contamination tolerance "
          "was never recorded for the original interactive 1232/1466 "
          "session. The live headline is the s10v2 recapture, job 2416746.")


if __name__ == "__main__":
    main()
