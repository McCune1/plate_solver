"""
probe_rect_ff_persist_threshold_sensitivity_2026-09-07.py

Purpose
-------
Paper 2's frozen "Screen B" persistence cut (p5_rect_ff_lib.PERSIST_DL =
0.01) is the gate that turns a discovery-basis dip into a published
production match: re-scan a +/-0.02 window (step 0.002) around the
published Lambda*/Omega-bar*, at an enlarged persist-basis (n_cpair=3 in
both families), and require the nearest persist-basis dip to sit within
PERSIST_DL of that published value. This is the rectangular companion's
direct analogue of Paper 1's frozen SUBDOM 0.70 threshold and geometry-
sweep tolerance -- both of which Paper 1 already ships an insensitivity/
margin check for (LESSONS_LEARNED.md). Paper 2 does not yet ship one for
its own frozen 0.01 cut. This probe closes exactly that gap.

This probe does NOT retune PERSIST_DL, does NOT edit plate_solver
internals, does NOT monkeypatch anything, and calls the exact, already-
frozen production functions p5_rect_ff_lib.persist_one / .persist_one_ip
VERBATIM (imported, not copy-pasted) once per published candidate. That
single call already returns the full "dips" list AND the "nearest" dip
found -- so delta = |nearest - L0| is known for every one of Paper 2's own
published production matches:

    OOP Table tab:oop-primary   26 rows  (Lambda* <~ 2.4)
    OOP Table tab:oop-high      30 rows  (higher Lambda*, miss <= 3%)
    IP  Table tab:ip92          92 rows  (all 5 aspect ratios, both parity)
    -------------------------------------------------------------
    TOTAL                      148 candidates -- every flexural and
                                extensional frequency match Paper 2
                                publishes anywhere in the main text.

Because delta is known post hoc for every candidate from ONE persist_one/
persist_one_ip call, checking threshold sensitivity costs ZERO extra
plate_solver calls past that: re-classify the same 148 deltas at a whole
sweep of alternate thresholds and report whether the PUBLISHED SET (which
candidates count as "persist") changes anywhere in that band. This is the
project's standalone-algebra-check technique -- only ONE quantity (delta)
needs computing on the cluster; everything past that is arithmetic done
in Python once the run finishes, so it costs nothing to also report it
at a dozen alternate thresholds instead of just one or two.

Verdict fields (see final printed JSON):
  - max_delta: the largest delta among all 148 published matches -- the
    number a reviewer will actually ask for ("how close does the frozen
    cut come to excluding one of your own results?").
  - margin_lo / margin_hi: the widest threshold band [lo, hi], containing
    the frozen 0.01, over which the classified "persist" set is IDENTICAL
    to the set at 0.01. margin_lo > max_delta by construction whenever
    every candidate passes; margin_hi is the smallest swept threshold
    at or above 0.01 that is NOT also safe (or the top of the sweep, if
    none exists) -- i.e. how far the cut could be loosened with no change.
  - flips: any candidate whose persist verdict changes at a swept
    threshold strictly below 0.01 (would have been rejected by a tighter
    cut) or, separately, any candidate that requires threshold > 0.01 to
    persist at all (delta > 0.01 -- would mean the frozen cut is already
    excluding published matches, which must not happen and is a hard
    stop if seen).

Reuses the exact production settings baked into persist_one/persist_one_ip
themselves (im_cap=30.0 OOP, im_cap=7.0 IP; persist-basis n_cpair=3 both
families) -- the same settings that already found every one of these dips
at production time -- so this run introduces no new branch-resolution
risk; only the classification threshold is swept, after the fact, in
plain Python arithmetic. Each persist_one/persist_one_ip call self-
checkpoints to its own JSON file keyed by (lob, sym, L0) under a checkpoint
directory private to this probe -- a re-submitted job resumes for free.

Fidelity gates (must BOTH pass before any candidate is scored):
  G0     SOLVER_VERSION + checkerboard-formula source check (OOP path)
  G0-IP  RectIPAssembler bc='free_free' switch present
  G1     l/b=1.0 IP SYM Omega=1.3300 (Table tab:ip92's first row -- the
         same point job 2456810 independently MAC-confirmed at 0.9997)
         must classify persist=True at threshold 0.01. If this single,
         independently cross-checked point fails, persist_one_ip itself
         is broken and nothing past G1 should be trusted.

Copy/paste on the cluster:
  cd /home/ghmkfh/PythonMill/Plate_Solver_Package
  sbatch submit_rect_ff_persist_threshold_sensitivity_2026-09-07.sh

Fidelity-gate-only smoke test (seconds, a handful of candidates only):
  SMOKE=1 sbatch submit_rect_ff_persist_threshold_sensitivity_2026-09-07.sh

Internal soft wall-clock budget (minutes; default 100):
  MAX_MINUTES=80 sbatch submit_rect_ff_persist_threshold_sensitivity_2026-09-07.sh
"""

from __future__ import annotations

import json
import os
import time

from mpmath import mp

import p5_rect_ff_lib as L
import plate_solver as ps

PKG = os.environ.get("PKG_PATH", ".")
EXPECT_VER = os.environ.get("EXPECT_SOLVER_VERSION", "2026-07-10.s10")
mp.dps = int(os.environ.get("DPS", "30"))
SMOKE = os.environ.get("SMOKE", "0") == "1"
MAX_MINUTES = float(os.environ.get("MAX_MINUTES", "100"))
T_START = time.time()

CHECKPOINT_DIR = os.environ.get(
    "CHECKPOINT_DIR", "./p5_persist_threshold_sensitivity_checkpoints"
)
OUT_JSON = "rect_ff_persist_threshold_sensitivity_2026-09-07.json"

# Frozen production value -- NOT retuned, only referenced.
PERSIST_DL = L.PERSIST_DL
assert abs(PERSIST_DL - 0.01) < 1e-12, f"PERSIST_DL drifted: {PERSIST_DL}"

# Sweep band for the post-hoc classification. Dense near 0.01 (the frozen
# cut), sparse further out -- this is arithmetic only, so it is free to be
# dense.
THRESHOLDS = [
    0.0005, 0.001, 0.0015, 0.002, 0.003, 0.004, 0.005, 0.0075,
    0.01, 0.0125, 0.015, 0.0175, 0.02, 0.025, 0.03, 0.04, 0.05,
]

# ---------------------------------------------------------------------
# Candidate tables, transcribed verbatim from the published paper:
#   github_repo/paper/PAPER2_RECT_FREEFREE_DRAFT.tex Table tab:oop-primary
#   github_repo/paper/PAPER2_RECT_FREEFREE_DRAFT.tex Table tab:oop-high
#   github_repo/paper/PAPER2_IP_MATCH_ROWS.tex        Table tab:ip92
# Each row: (lob, sym, L0_or_O0). sym=True is SYM, sym=False is ANTI.
# ---------------------------------------------------------------------

OOP_PRIMARY = [
    (1.0, False, 1.366),
    (1.5, True, 0.964), (1.5, True, 2.254),
    (1.5, False, 0.906), (1.5, False, 2.124),
    (2.0, True, 0.5434), (2.0, True, 1.510), (2.0, True, 2.234),
    (2.0, False, 0.674), (2.0, False, 1.482),
    (2.5, True, 0.348), (2.5, True, 0.978), (2.5, True, 1.888),
    (2.5, True, 2.278),
    (2.5, False, 0.534), (2.5, False, 1.148), (2.5, False, 1.918),
    (3.0, True, 0.242), (3.0, True, 0.6698), (3.0, True, 1.320),
    (3.0, True, 2.162), (3.0, True, 2.256),
    (3.0, False, 0.444), (3.0, False, 0.936), (3.0, False, 1.528),
    (3.0, False, 2.258),
]
assert len(OOP_PRIMARY) == 26, len(OOP_PRIMARY)

OOP_HIGH = [
    (1.0, True, 2.460), (1.0, True, 3.530), (1.0, True, 6.190),
    (1.0, True, 6.460),
    (1.0, False, 3.520),
    (1.5, True, 5.380),
    (1.5, False, 3.860),
    (2.0, True, 3.000), (2.0, True, 3.650),
    (2.0, False, 2.580), (2.0, False, 4.060), (2.0, False, 6.200),
    (2.5, True, 3.190), (2.5, True, 4.180), (2.5, True, 5.410),
    (2.5, False, 2.920), (2.5, False, 4.180), (2.5, False, 5.720),
    (2.5, False, 6.260), (2.5, False, 6.400),
    (3.0, True, 2.480), (3.0, True, 3.340), (3.0, True, 3.660),
    (3.0, True, 5.660), (3.0, True, 6.180),
    (3.0, False, 3.160), (3.0, False, 4.260), (3.0, False, 5.560),
    (3.0, False, 6.220), (3.0, False, 6.380),
]
assert len(OOP_HIGH) == 30, len(OOP_HIGH)

IP_92 = [
    (1.0, True, 1.3300), (1.0, True, 1.4100), (1.0, True, 1.6100),
    (1.0, True, 1.8600), (1.0, True, 2.0000),
    (1.0, False, 1.2400), (1.0, False, 1.3298), (1.0, False, 2.0000),
    (1.0, False, 2.3200),
    (1.5, True, 1.0500), (1.5, True, 1.4100), (1.5, True, 1.4243),
    (1.5, True, 1.6300), (1.5, True, 1.6900), (1.5, True, 2.1100),
    (1.5, True, 2.2200), (1.5, True, 2.2800),
    (1.5, False, 0.7800), (1.5, False, 1.0333), (1.5, False, 1.5712),
    (1.5, False, 1.6400), (1.5, False, 2.2025), (1.5, False, 2.2400),
    (2.0, True, 0.8100), (2.0, True, 1.4000), (2.0, True, 1.4142),
    (2.0, True, 1.4400), (2.0, True, 1.6500), (2.0, True, 1.7500),
    (2.0, True, 1.8100), (2.0, True, 2.0000), (2.0, True, 2.3000),
    (2.0, True, 2.3300),
    (2.0, False, 0.5200), (2.0, False, 0.8800), (2.0, False, 1.2715),
    (2.0, False, 1.2800), (2.0, False, 1.7400), (2.0, False, 1.8400),
    (2.0, False, 2.1400),
    (2.5, True, 0.6400), (2.5, True, 1.2200), (2.5, True, 1.4100),
    (2.5, True, 1.4266), (2.5, True, 1.5400), (2.5, True, 1.6700),
    (2.5, True, 1.6900), (2.5, True, 1.7800), (2.5, True, 2.0500),
    (2.5, True, 2.0700), (2.5, True, 2.1300), (2.5, True, 2.3900),
    (2.5, True, 2.4047),
    (2.5, False, 0.3800), (2.5, False, 0.7200), (2.5, False, 1.0600),
    (2.5, False, 1.1125), (2.5, False, 1.4400), (2.5, False, 1.5400),
    (2.5, False, 1.8000), (2.5, False, 2.0000), (2.5, False, 2.2000),
    (2.5, False, 2.2800),
    (3.0, True, 0.5400), (3.0, True, 1.0600), (3.0, True, 1.4142),
    (3.0, True, 1.4199), (3.0, True, 1.4200), (3.0, True, 1.6000),
    (3.0, True, 1.6635), (3.0, True, 1.6900), (3.0, True, 1.8000),
    (3.0, True, 1.8400), (3.0, True, 2.0000), (3.0, True, 2.1600),
    (3.0, True, 2.2100), (3.0, True, 2.2700), (3.0, True, 2.4300),
    (3.0, True, 2.4600),
    (3.0, False, 0.2800), (3.0, False, 0.5800), (3.0, False, 0.8800),
    (3.0, False, 1.0315), (3.0, False, 1.2702), (3.0, False, 1.3000),
    (3.0, False, 1.5800), (3.0, False, 1.6600), (3.0, False, 1.9000),
    (3.0, False, 2.0000), (3.0, False, 2.2400), (3.0, False, 2.2800),
    (3.0, False, 2.4937),
]
assert len(IP_92) == 92, len(IP_92)


def budget_left():
    return MAX_MINUTES - (time.time() - T_START) / 60.0


def checkpoint(rows_so_far):
    tmp = OUT_JSON + ".tmp"
    with open(tmp, "w") as f:
        json.dump(rows_so_far, f, indent=1)
    os.replace(tmp, OUT_JSON)


def score_candidate(asm_oop, asm_ip, family, lob, sym, L0):
    if family == "oop":
        pr = L.persist_one(asm_oop, lob, sym, L0, CHECKPOINT_DIR)
    else:
        pr = L.persist_one_ip(asm_ip, lob, sym, L0, CHECKPOINT_DIR)
    nearest = pr.get("nearest")
    delta = None if nearest is None else abs(nearest[0] - L0)
    return {
        "family": family,
        "lob": lob,
        "sym": "SYM" if sym else "ANTI",
        "L0": L0,
        "persist_at_0p01": bool(pr["persist"]),
        "nearest": nearest,
        "delta": delta,
        "n_dips": len(pr.get("dips") or []),
    }


def classify_at(rows, thr):
    """Reclassify every row's persist verdict at an alternate threshold,
    purely from its already-computed delta -- no plate_solver call."""
    out = []
    for r in rows:
        if r["delta"] is None:
            out.append(False)  # no dip found at all -> never persists
        else:
            out.append(r["delta"] <= thr)
    return out


def main():
    print("probe_rect_ff_persist_threshold_sensitivity_2026-09-07")
    print("start", time.strftime("%Y-%m-%d %H:%M:%S"), f"dps={mp.dps}",
          f"smoke={SMOKE}", flush=True)

    if not L.deploy_ok_ip(PKG, EXPECT_VER):
        raise SystemExit(2)
    print(f"SOLVER_VERSION={ps.SOLVER_VERSION}  PERSIST_DL={PERSIST_DL}",
          flush=True)

    os.makedirs(CHECKPOINT_DIR, exist_ok=True)
    asm_oop = L.make_ff()
    asm_ip = L.make_ff_ip()

    print("=" * 72)
    print("G1 -- IP SYM l/b=1.0 Omega=1.3300 fidelity gate")
    print("=" * 72, flush=True)
    g1 = score_candidate(asm_oop, asm_ip, "ip", 1.0, True, 1.3300)
    print(f"  persist={g1['persist_at_0p01']}  delta={g1['delta']}", flush=True)
    if not g1["persist_at_0p01"]:
        print("G1 FAIL -- persist_one_ip did not reproduce a known-published "
              "match as persisting. STOP: do not trust anything past this.")
        raise SystemExit(3)
    print("  G1 PASS", flush=True)

    candidates = []
    for lob, sym, L0 in OOP_PRIMARY:
        candidates.append(("oop", lob, sym, L0))
    for lob, sym, L0 in OOP_HIGH:
        candidates.append(("oop", lob, sym, L0))
    for lob, sym, L0 in IP_92:
        candidates.append(("ip", lob, sym, L0))

    if SMOKE:
        candidates = candidates[:3] + candidates[-3:]
        print(f"SMOKE=1: restricting to {len(candidates)} candidates", flush=True)

    print(f"\n=== Scoring {len(candidates)} published production matches ===",
          flush=True)
    rows = []
    t0 = time.time()
    for i, (family, lob, sym, L0) in enumerate(candidates):
        if budget_left() < 5.0:
            print(f"STOPPING EARLY at {i}/{len(candidates)} -- "
                  f"budget_left={budget_left():.1f} min. Checkpoints are "
                  f"saved per-candidate; resubmit to resume for free.",
                  flush=True)
            break
        row = score_candidate(asm_oop, asm_ip, family, lob, sym, L0)
        rows.append(row)
        sl = row["sym"]
        d = row["delta"]
        print(f"  [{i+1:3d}/{len(candidates)}] {family:3s} lob={lob:<4} "
              f"{sl:4s} L0={L0:<8} persist={row['persist_at_0p01']!s:5s} "
              f"delta={('%.5f' % d) if d is not None else 'NONE':>8}",
              flush=True)
        if (i + 1) % 10 == 0:
            checkpoint(rows)
    checkpoint(rows)
    dt = time.time() - t0
    print(f"\nScoring done: {len(rows)}/{len(candidates)} candidates, "
          f"{dt:.1f}s total, {dt / max(1, len(rows)):.3f}s/candidate avg",
          flush=True)

    complete = len(rows) == len(candidates)

    # ---- Verdict: pure arithmetic on already-computed deltas ----
    deltas = [r["delta"] for r in rows if r["delta"] is not None]
    no_dip = [r for r in rows if r["delta"] is None]
    max_delta = max(deltas) if deltas else None
    baseline = classify_at(rows, PERSIST_DL)
    baseline_all_true = all(baseline) and not no_dip

    sweep = {}
    for thr in THRESHOLDS:
        verdict = classify_at(rows, thr)
        n_persist = sum(verdict)
        same_as_baseline = verdict == baseline
        sweep[f"{thr:.4f}"] = {
            "threshold": thr,
            "n_persist": n_persist,
            "n_total": len(rows),
            "identical_set_to_0p01": same_as_baseline,
        }

    # widest contiguous band around 0.01 with identical classified set
    sorted_thr = sorted(THRESHOLDS)
    idx01 = sorted_thr.index(min(sorted_thr, key=lambda t: abs(t - 0.01)))
    margin_lo = sorted_thr[idx01]
    j = idx01
    while j > 0 and sweep[f"{sorted_thr[j-1]:.4f}"]["identical_set_to_0p01"]:
        j -= 1
        margin_lo = sorted_thr[j]
    margin_hi = sorted_thr[idx01]
    k = idx01
    while (k < len(sorted_thr) - 1
           and sweep[f"{sorted_thr[k+1]:.4f}"]["identical_set_to_0p01"]):
        k += 1
        margin_hi = sorted_thr[k]

    flips_below = [r for r in rows
                   if r["delta"] is not None and PERSIST_DL / 5.0 < r["delta"] <= PERSIST_DL]
    exceeds_frozen = [r for r in rows
                      if r["delta"] is not None and r["delta"] > PERSIST_DL]

    verdict = {
        "n_candidates_total": len(candidates),
        "n_candidates_scored": len(rows),
        "complete": complete,
        "g1_fidelity_gate": g1,
        "baseline_all_persist_at_0p01": baseline_all_true,
        "n_no_dip_found": len(no_dip),
        "no_dip_rows": no_dip,
        "max_delta_confirmed": max_delta,
        "margin_lo": margin_lo,
        "margin_hi": margin_hi,
        "frozen_threshold": PERSIST_DL,
        "candidates_within_5x_of_cut": [
            {"family": r["family"], "lob": r["lob"], "sym": r["sym"],
             "L0": r["L0"], "delta": r["delta"]}
            for r in flips_below
        ],
        "candidates_exceeding_frozen_cut": [
            {"family": r["family"], "lob": r["lob"], "sym": r["sym"],
             "L0": r["L0"], "delta": r["delta"]}
            for r in exceeds_frozen
        ],
        "sweep": sweep,
    }

    print("\n" + "=" * 72)
    print("VERDICT")
    print("=" * 72)
    if not complete:
        print("INCOMPLETE RUN -- resubmit to finish scoring before trusting "
              "any margin/max_delta number below.")
    if exceeds_frozen:
        print(f"HARD STOP CANDIDATE: {len(exceeds_frozen)} published match(es) "
              f"have delta > frozen PERSIST_DL=0.01 -- this should be "
              f"impossible (they are published as persisting) and means "
              f"persist_one/persist_one_ip disagree with production. "
              f"Do not fold into Paper 2 before investigating.")
    elif baseline_all_true and complete:
        print(f"All {len(rows)} published production matches (OOP primary + "
              f"OOP high + IP-92) reproduce persist=True at the frozen "
              f"PERSIST_DL=0.01, as they must. Largest delta among all of "
              f"them is {max_delta:.5f} ({100*max_delta/PERSIST_DL:.1f}% of "
              f"the frozen cut). The classified set is UNCHANGED for any "
              f"threshold in [{margin_lo:.4f}, {margin_hi:.4f}] -- a "
              f"{margin_hi/max(margin_lo,1e-9):.1f}x band around the frozen "
              f"0.01 cut. {'No' if not flips_below else str(len(flips_below))}"
              f" candidate(s) sit within 5x of the cut "
              f"(delta > {PERSIST_DL/5.0:.4f}).")
    else:
        print("Some published matches did NOT reproduce persist=True -- "
              "see no_dip_rows / candidates_exceeding_frozen_cut above.")
    print(json.dumps({k: v for k, v in verdict.items()
                       if k not in ("sweep",)}, indent=1))

    out = {"rows": rows, "verdict": verdict}
    checkpoint(out)
    print(f"\nWrote {OUT_JSON}", flush=True)


if __name__ == "__main__":
    main()
