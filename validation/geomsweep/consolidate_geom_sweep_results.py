#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
consolidate_geom_sweep_results.py -- flatten the freefree_manifest.json files
written by each geometry-sweep job (submit_ff_geomsweep_*.sh) into one CSV
table covering the full (r0/2b, 2Theta/pi) matrix, both parts.

Standalone: stdlib only (json/csv), does NOT import plate_solver -- run this
anywhere (login node, laptop, this sandbox) after copying the FIGDIR folders
back from the cluster. Does not touch the deployed package; no
SOLVER_VERSION implications.

WHY: each geometry-sweep job writes its own freefree_manifest.json (schema:
top-level key "<TAG> <r0_2b> <two_T_pi>" -> {raw, omega_lit, k_bar, w_bar,
error, diag}, see validation.py::_run_quality_review) into its OWN FIGDIR so
the jobs can run concurrently without clobbering each other. This script is
the thing that puts the pieces back together into one table you can actually
read, and flags the same quality signals _run_quality_review already
computes per-run (marginal/dropped/NaN) so a sparse-looking geometry is
visible immediately instead of silently missing from the table.

Usage:
    python3 consolidate_geom_sweep_results.py [FIGDIR ...]

    With no arguments, looks for the FIGDIRs this project's geometry-sweep
    submit scripts write by convention:
      figures/                             (original FF-P1, r0_2b=1.25)
      figures_ff_r150/  figures_ff_r167/  figures_ff_r200/  figures_ff_r250/
      figures_ff_ip_anglefill_r125/
    (relative to the package root, or pass explicit paths as arguments.)

Output:
    geom_sweep_consolidated.csv  -- one row per (tag, r0_2b, two_T_pi, mode)
    Printed summary: geometry count per ratio, any geometry with 0 accepted
    modes or a NaN flag, and the total mode count found so far.
"""
from __future__ import annotations
import csv
import json
import os
import sys

DEFAULT_FIGDIRS = [
    "figures",
    "figures_ff_r150",
    "figures_ff_r167",
    "figures_ff_r200",
    "figures_ff_r250",
    "figures_ff_ip_anglefill_r125",
]

OUT_CSV = "geom_sweep_consolidated.csv"


def _load_manifest(figdir):
    path = os.path.join(figdir, "freefree_manifest.json")
    if not os.path.isfile(path):
        return None
    with open(path) as f:
        return json.load(f)


def main(argv):
    figdirs = argv[1:] if len(argv) > 1 else DEFAULT_FIGDIRS

    rows = []          # flattened per-mode rows for the CSV
    geom_status = []    # one line per geometry for the printed summary
    n_missing_figdir = 0

    for figdir in figdirs:
        manifest = _load_manifest(figdir)
        if manifest is None:
            n_missing_figdir += 1
            print(f"  [skip] {figdir}: no freefree_manifest.json found "
                  f"(job not run yet, or not copied back from the cluster)")
            continue

        for key, rec in sorted(manifest.items()):
            if key == "_cross_geometry_suspects":
                continue
            parts = key.split(" ")
            if len(parts) != 3:
                continue
            tag, r0_2b_s, two_T_s = parts
            r0_2b = float(r0_2b_s)
            two_T = float(two_T_s)

            if rec.get("error"):
                geom_status.append(
                    (figdir, tag, r0_2b, two_T, 0, f"ERROR: {rec['error']}"))
                continue

            raw = rec.get("raw") or []
            f_hz = rec.get("f_hz") or ([None] * len(raw))
            om_lit = rec.get("omega_lit") or ([None] * len(raw))
            diag = rec.get("diag") or {}
            n_nan = diag.get("n_nan", 0)
            dropped = diag.get("dropped", [])
            marginal_n = sum(
                1 for a in diag.get("accepted", [])
                if a.get("margin_below_floor", 999) < 0.5)

            flags = []
            if not raw:
                flags.append("ZERO MODES FOUND")
            if n_nan:
                flags.append(f"{n_nan} NaN scan pts")
            if dropped:
                flags.append(f"{len(dropped)} dropped candidate(s)")
            if marginal_n:
                flags.append(f"{marginal_n} marginal accepted mode(s)")
            status = "; ".join(flags) if flags else "OK"
            geom_status.append((figdir, tag, r0_2b, two_T, len(raw), status))

            for n, (rw, fz, om) in enumerate(zip(raw, f_hz, om_lit), 1):
                rows.append({
                    "figdir": figdir, "tag": tag,
                    "r0_2b": r0_2b, "two_T_pi": two_T, "mode": n,
                    "raw_omega": rw, "f_hz": fz, "omega_lit": om,
                })

    if not rows and n_missing_figdir == len(figdirs):
        print("\nNo manifests found in any of the given/default FIGDIRs. "
              "Copy the figures_ff_*/ folders back from the cluster first.")
        return 1

    with open(OUT_CSV, "w", newline="") as f:
        w = csv.DictWriter(
            f, fieldnames=["figdir", "tag", "r0_2b", "two_T_pi", "mode",
                           "raw_omega", "f_hz", "omega_lit"])
        w.writeheader()
        w.writerows(rows)
    print(f"\nwrote {len(rows)} rows -> {OUT_CSV}")

    print("\n" + "=" * 78)
    print("GEOMETRY-LEVEL SUMMARY  (figdir, tag, r0/2b, 2Theta/pi, n_modes, flags)")
    print("=" * 78)
    n_bad = 0
    for figdir, tag, r0_2b, two_T, n_modes, status in geom_status:
        if status != "OK":
            n_bad += 1
        print(f"  {figdir:32} {tag:6} r0/2b={r0_2b:<8.4g} 2T/pi={two_T:<6.4g} "
              f"n={n_modes:<3d} {status}")
    print("-" * 78)
    print(f"  {len(geom_status)} geometries consolidated, "
          f"{n_bad} flagged for review, {n_missing_figdir} FIGDIR(s) not "
          f"yet available.")
    if n_bad:
        print("  NOTE: a flag here does NOT mean the result is wrong -- it means "
              "_run_quality_review already thought it was worth a human look. "
              "Cross-check flagged rows against the job's own .out log before "
              "trusting or discarding.")


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
