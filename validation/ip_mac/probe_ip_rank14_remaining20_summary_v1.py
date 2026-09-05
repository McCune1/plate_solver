# -*- coding: utf-8 -*-
"""
probe_ip_rank14_remaining20_summary_v1.py -- DIAGNOSTIC ONLY.

Reads completed Rank 14 IP MAC-bar .out files and prints one table plus
the combined letter. Does not run the solver. Does not retune 0.744.
Does not edit the paper.

Default: the 20 remaining keys (4 ratios x 5 angles, a100 excluded).
Set INCLUDE_A100=1 to add the four already-closed half-annulus keys.

-----------------------------------------------------------------------------
PRE-REGISTERED INTERPRETATION
-----------------------------------------------------------------------------
R14_REM_INCOMPLETE
    Fewer than the expected keys have a READING letter. Report which
    are missing. Not a bar result. Do not retune 0.744. Do not edit
    the paper.
R14_REM_ANY_FAIL
    At least one completed key is R14-0 GATE_FAIL. That key is VOID.
    Do not retune 0.744. Re-run that key only.
R14_REM_LETTERS
    Every expected key has a letter. Print the per-key table, combined
    tight REAL floor / ART ceiling, and whether 0.744 still splits
    every tight candidate. Rank 14 IP closes only when all 24 keys
    (these 20 plus the four a100) have a letter. This script does not
    close it, retune, or edit the paper.
"""
from __future__ import annotations

import glob
import os
import re
import sys

PKG = os.environ.get(
    "PKG_PATH",
    os.path.dirname(os.path.abspath(__file__)) or ".",
)
INCLUDE_A100 = os.environ.get("INCLUDE_A100", "0") == "1"
MAP_PATH = os.environ.get(
    "R14_MAP", os.path.join(PKG, "rank14_remaining20_jobs.txt")
)

RATIOS = ("r150", "r167", "r200", "r250")
REM_ANGLES = ("a025", "a050", "a075", "a125", "a150")
A100 = ("a100",)
EXPECTED = [f"{r}/{a}" for r in RATIOS for a in REM_ANGLES]
if INCLUDE_A100:
    EXPECTED = [f"{r}/{a}" for r in RATIOS for a in A100 + REM_ANGLES]

ROW = re.compile(
    r"^\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(q[01])\s+(\S+)\s+"
    r"([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+#(\d+)\s+([01])\s+"
    r"(REAL-like|ARTIFACT-like)\s*$"
)
KEY_HDR = re.compile(
    r"RATIO_TAG=(\S+)\s+ANGLE_TAG=(\S+)"
)
READING = re.compile(r"^READING:\s+(\S+)")
TIGHT_LINE = re.compile(
    r"tight \(rel<0\.5%\)\s*:\s*(\d+)\s+MAC-REAL=(\d+)\s+MAC-ART=(\d+)"
)
RECON = re.compile(r"reconstructed\s+(\d+)/(\d+)")
THRESH = 0.744


def parse_out(path):
    text = open(path, encoding="utf-8", errors="replace").read()
    kh = KEY_HDR.search(text)
    if not kh:
        return None
    key = f"{kh.group(1)}/{kh.group(2)}"
    rec = RECON.search(text)
    letter = None
    for line in text.splitlines():
        m = READING.match(line)
        if m:
            letter = m.group(1)
    tight_meta = TIGHT_LINE.search(text)
    rows = []
    in_sum = False
    for line in text.splitlines():
        if line.strip().startswith("Om") and "max-sc" in line:
            in_sum = True
            continue
        if in_sum and line.startswith("="):
            break
        if in_sum:
            m = ROW.match(line)
            if m:
                rows.append(
                    dict(
                        maxsc=float(m.group(7)),
                        tight=m.group(9) == "1",
                        call=m.group(10),
                    )
                )
    tight = [r for r in rows if r["tight"]]
    tr = [r for r in tight if r["call"] == "REAL-like"]
    ta = [r for r in tight if r["call"] == "ARTIFACT-like"]
    return dict(
        key=key,
        path=path,
        letter=letter,
        recon=rec.group(0) if rec else "?",
        n_ok=int(rec.group(1)) if rec else 0,
        n_tot=int(rec.group(2)) if rec else 0,
        n_tight=int(tight_meta.group(1)) if tight_meta else len(tight),
        n_tr=int(tight_meta.group(2)) if tight_meta else len(tr),
        n_ta=int(tight_meta.group(3)) if tight_meta else len(ta),
        real_floor=min((r["maxsc"] for r in tr), default=None),
        art_ceil=max((r["maxsc"] for r in ta), default=None),
        mis_r=sum(1 for r in tr if r["maxsc"] < THRESH),
        mis_a=sum(1 for r in ta if r["maxsc"] >= THRESH),
    )


def main():
    print("probe_ip_rank14_remaining20_summary_v1")
    print(f"  PKG={os.path.abspath(PKG)}")
    print(f"  INCLUDE_A100={int(INCLUDE_A100)}  expected={len(EXPECTED)}")
    print(f"  Do not retune {THRESH}. Do not edit the paper from this script.")
    paths = sorted(
        glob.glob(os.path.join(PKG, "probe_ip_rank14_mac_bar_v1_*.out"))
    )
    by_key = {}
    for p in paths:
        rec = parse_out(p)
        if rec is None:
            continue
        prev = by_key.get(rec["key"])
        if prev is None or os.path.getmtime(p) >= os.path.getmtime(prev["path"]):
            by_key[rec["key"]] = rec

    print(f"\n  parsed {len(by_key)} unique keys from {len(paths)} .out files")
    print(
        f"  {'key':<12s} {'letter':<16s} {'recon':<14s} "
        f"{'tight':>5s} {'R/A':>7s} {'floor':>8s} {'ceil':>8s}"
    )
    missing = []
    fails = []
    floors = []
    ceils = []
    mis = 0
    for key in EXPECTED:
        rec = by_key.get(key)
        if rec is None or not rec["letter"]:
            missing.append(key)
            print(f"  {key:<12s} MISSING")
            continue
        if rec["letter"].startswith("R14-0"):
            fails.append(key)
        fl = rec["real_floor"]
        ce = rec["art_ceil"]
        if fl is not None:
            floors.append(fl)
        if ce is not None:
            ceils.append(ce)
        mis += rec["mis_r"] + rec["mis_a"]
        print(
            f"  {key:<12s} {rec['letter']:<16s} {rec['recon']:<14s} "
            f"{rec['n_tight']:5d} {rec['n_tr']:3d}/{rec['n_ta']:<3d} "
            f"{(fl if fl is not None else float('nan')):8.4f} "
            f"{(ce if ce is not None else float('nan')):8.4f}"
        )

    extra = sorted(k for k in by_key if k not in EXPECTED)
    if extra:
        print(f"\n  other keys in .out files (not in this expected set): {extra}")

    print("\n" + "=" * 70)
    print("  READING (honour this letter; do not retune 0.744)")
    print("=" * 70)
    if missing:
        print(f"READING: R14_REM_INCOMPLETE  missing={len(missing)}/{len(EXPECTED)}")
        print("  " + " ".join(missing))
        print("  Wait or re-run those keys. Not a bar result.")
        return 0
    if fails:
        print(f"READING: R14_REM_ANY_FAIL  {fails}")
        print("  Re-run the failed key(s) only. Do not retune 0.744.")
        return 0
    cfl = min(floors) if floors else float("nan")
    cce = max(ceils) if ceils else float("nan")
    print("READING: R14_REM_LETTERS")
    print(f"  {len(EXPECTED)}/{len(EXPECTED)} keys have a letter.")
    print(f"  combined tight REAL floor {cfl:.4f}  ART ceiling {cce:.4f}")
    print(f"  tight misclassifications at {THRESH}: {mis}")
    if mis == 0 and cce < THRESH < cfl:
        print(f"  {THRESH} still splits every tight candidate on this set.")
    else:
        print(f"  {THRESH} does NOT split this set the same way. Report; do not retune.")
    if not INCLUDE_A100:
        print("  a100 (4 keys) not included. Rank 14 IP close is 24/24 keys.")
    else:
        print("  INCLUDE_A100=1: this is the 24-key set.")
    print("  Do not edit the paper from this script.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
