#!/usr/bin/env python3
"""
Production-candidate IP discovery v2: replace the n_cpair=0 discovery
default with n_cpair=3 (the SAME basis already used, project-wide, as
the persist/confirmation basis) as the discovery scan basis itself.

Motivation (jobs 2455676, 2455677, both 2026-09-05):
  - 2455676 found 14 dips in the Omega in [1.20,1.50] band that exist
    ONLY at n_cpair=3 and are invisible at n_cpair=0 -- the band was
    never "uncovered physics", it was basis starvation in the n_cpair=0
    discovery step.
  - 2455677 found that 20/26 ANTI FAIL_P discovery hits from 2455570
    recover to <=3% miss once persist is re-run at n_cpair=5 with a
    wider +/-0.05 window (10 basis-limited -- already inside the old
    +/-0.02 window, just needed more cpairs; 10 window-limited -- needed
    the wider window too). Only 3/26 ANTI targets (all l/b=1.0: O0=1.40,
    1.94, 2.10) show no dip at all even at n_cpair=5 +/-0.05.
  - 2455678 confirmed 37/38 of the original n_cpair=0-discovered MATCHes
    hold at a tightened 3% cut after refine+persist.

Net effect: n_cpair=0 discovery is demonstrably under-resourced for this
problem. This probe re-runs discovery from scratch over the full
[0.02,2.50] range at n_cpair=3 directly (same scan step convention as
2455570: 0.01 SYM / 0.02 ANTI), then fine-refines every dip found (same
basis, +/-0.02 window, step 0.002, via persist_one_ip -- this sharpens
the dip location and screens out single-point grid noise, it is NOT a
second independent basis check since the basis does not change).

Does NOT retune SIGMA_LIST, PERSIST_DL, or MATCH_CUT semantics beyond
adopting the ALREADY-VALIDATED 3% cut from 2455678 (superseding the
5% trial cut 2455570 used before that refine existed). Does NOT touch
the corner term (none exists for IP). Does NOT bump SOLVER_VERSION.
This is a basis upgrade to the discovery step, not a new physics claim.

PRE-REGISTERED
--------------
  MATCH if a refine-confirmed dip is within 3% of a same-parity IP FE
  Omega (the paper's adopted cut per 2455678).
  Cross-reference every MATCH against the OLD 2455570 n_cpair=0 unique
  MATCH list (38 entries, see p5_rect_ff_lib comments / draft SM): tag
  NEW if this (l/b, parity, FE) triple was not in that list, else OLD.
  UNCOVERED = FE targets in [0.02,2.50] still unmatched at n_cpair=3.
  A non-empty UNCOVERED list after this upgrade is a genuine open item
  for the paper (basis-independent absence), not evidence of a missing
  term -- there is no IP corner to add (Seok Part 2 Eq. 1/36).
"""
from __future__ import annotations

import json
import os
import sys
import time

sys.path.insert(0, os.environ.get("PKG_PATH", "."))
os.environ.setdefault("DPS", "30")

from mpmath import mp  # noqa: E402

mp.dps = int(os.environ["DPS"])

import plate_solver as ps  # noqa: E402
import p5_rect_ff_lib as L  # noqa: E402

EXPECT_VER = os.environ.get("EXPECT_SOLVER_VERSION", "2026-07-10.s10")
CHECKPOINT_DIR = os.environ.get("CHECKPOINT_DIR", "./p5_ip_disc_v2_checkpoints")
PKG = os.environ.get("PKG_PATH", ".")
DISC_LO, DISC_HI = 0.02, 2.50
MATCH_CUT = 3.0
GEOMS = (1.0, 1.5, 2.0, 2.5, 3.0)

# Old 2455570 unique MATCH list (n_cpair=0 discovery, 5% cut), collapsed
# to (lob, sym-as-str, round(Ofe,5)) keys, for NEW-vs-OLD tagging. Copied
# from that job's own printed "MATCHES in [0.02, 2.50]" block.
OLD_MATCHES_38 = {
    (1.0, "SYM", 1.85745), (1.0, "SYM", 2.00319), (1.0, "ANTI", 1.24858),
    (1.5, "SYM", 1.62439), (1.5, "SYM", 1.69877), (1.5, "SYM", 2.10886),
    (1.5, "SYM", 2.28542),
    (2.0, "SYM", 1.65355), (2.0, "SYM", 1.77390), (2.0, "SYM", 1.81566),
    (2.0, "SYM", 2.00440), (2.0, "SYM", 2.29856),
    (2.0, "ANTI", 0.52557), (2.0, "ANTI", 0.87891), (2.0, "ANTI", 1.27148),
    (2.5, "SYM", 1.67003), (2.5, "SYM", 1.69350), (2.5, "SYM", 1.78140),
    (2.5, "SYM", 2.04232), (2.5, "SYM", 2.07671), (2.5, "SYM", 2.13048),
    (2.5, "SYM", 2.38996),
    (2.5, "ANTI", 0.37585), (2.5, "ANTI", 0.71776), (2.5, "ANTI", 1.53132),
    (3.0, "SYM", 1.66346), (3.0, "SYM", 1.67627), (3.0, "SYM", 1.80401),
    (3.0, "SYM", 1.83870), (3.0, "SYM", 2.00478), (3.0, "SYM", 2.16287),
    (3.0, "SYM", 2.21145), (3.0, "SYM", 2.27554), (3.0, "SYM", 2.44758),
    (3.0, "SYM", 2.46034),
    (3.0, "ANTI", 0.28169), (3.0, "ANTI", 0.57789), (3.0, "ANTI", 1.57876),
}
assert len(OLD_MATCHES_38) == 38, len(OLD_MATCHES_38)
# Includes l/b=1.5 SYM 2.10886 (the 2455678 refine DROP at the 3% cut).
# Kept in OLD_MATCHES_38 anyway: this set's only job is NEW-vs-OLD
# provenance tagging against what 2455570 originally reported, not a
# restatement of which entries still pass the 3% cut today -- a v2 hit
# on 2.10886 would correctly tag OLD (it was already known, not a new
# discovery), even though 2455678 separately dropped it from the paper's
# adopted match list at the tighter cut.


def scan_disc_v2(asm, lob, sym):
    sl = "SYM" if sym else "ANTI"
    path = os.path.join(CHECKPOINT_DIR, f"disc2_ip_lob{lob}_{sl}.json")
    if os.path.exists(path):
        print(f"  SKIP {path}", flush=True)
        with open(path) as f:
            return json.load(f)
    n_real, n_cpair = L.n_basis_ip(persist=True)  # (5,3), forced at discovery
    step = 0.01 if sym else 0.02
    t0 = time.time()
    rows, dips = L.scan_window_ip(asm, lob, sym, DISC_LO, DISC_HI, step, 7.0,
                                  n_real, n_cpair)
    listed = [(d[0], float(d[1]), d[2]) for d in dips if d[1] < L.SIGMA_LIST]
    dt = time.time() - t0
    print(f"  [disc-v2 IP l/b={lob} {sl}] n_cpair={n_cpair} {len(rows)} pts "
          f"{dt:.1f}s  dips: "
          + (", ".join(f"{a:.4f}({b:.2e})" for a, b, _c in listed) or "NONE"),
          flush=True)
    out = {"dips": [[a, b, c] for a, b, c in listed], "elapsed_s": dt}
    os.makedirs(CHECKPOINT_DIR, exist_ok=True)
    with open(path, "w") as f:
        json.dump(out, f)
    return out


def main():
    print("probe_rect_ff_ip_disc_v2_ncp3_2026-09-05")
    print("start", time.strftime("%Y-%m-%d %H:%M:%S"), f"dps={mp.dps}")
    if not L.deploy_ok_ip(PKG, EXPECT_VER):
        raise SystemExit(2)
    os.makedirs(CHECKPOINT_DIR, exist_ok=True)
    asm = L.make_ff_ip()
    rows_out = []
    print("\n=== IP discovery v2: scan basis = persist basis (5,3) ===",
          flush=True)
    for lob in GEOMS:
        for sym in (True, False):
            sl = "SYM" if sym else "ANTI"
            disc = scan_disc_v2(asm, lob, sym)
            for O0, s0, n0 in disc["dips"]:
                pr = L.persist_one_ip(asm, lob, sym, O0, CHECKPOINT_DIR)
                Ofe, miss = L.nearest_ip_fe(lob, sym, O0)
                stable = pr["persist"]
                tag = "UNSTABLE"
                if stable:
                    tag = "PASS_NOFE"
                    if miss is not None and miss <= MATCH_CUT:
                        tag = "MATCH"
                key = (lob, sl, round(Ofe, 5) if Ofe is not None else None)
                novelty = ""
                if tag == "MATCH":
                    novelty = "NEW" if key not in OLD_MATCHES_38 else "OLD"
                print(f"  l/b={lob} {sl} O={O0:.4f} sig={s0:.2e}  "
                      f"stable={stable}  {tag} {novelty}  FE={Ofe} miss="
                      f"{None if miss is None else '%.2f%%' % miss}",
                      flush=True)
                rows_out.append({
                    "lob": lob, "sym": sl, "O0": O0, "sigma0": s0,
                    "stable": stable, "tag": tag, "novelty": novelty,
                    "Ofe": Ofe, "miss": miss,
                })

    print(f"\n=== MATCHES in [{DISC_LO},{DISC_HI}] at {MATCH_CUT}% cut, "
          "n_cpair=3 discovery ===", flush=True)
    matches = [r for r in rows_out if r["tag"] == "MATCH"]
    collapsed = {}
    for r in matches:
        key = (r["lob"], r["sym"], round(r["Ofe"], 5) if r["Ofe"] else None)
        if key not in collapsed or r["miss"] < collapsed[key]["miss"]:
            collapsed[key] = r
    n_new = 0
    for r in collapsed.values():
        key = (r["lob"], r["sym"], round(r["Ofe"], 5) if r["Ofe"] else None)
        novelty = "NEW" if key not in OLD_MATCHES_38 else "OLD"
        if novelty == "NEW":
            n_new += 1
        print(f"  l/b={r['lob']} {r['sym']}  O={r['O0']:.4f}  "
              f"FE={r['Ofe']:.5f}  miss={r['miss']:.2f}%  {novelty}")
    print(f"  unique MATCH count: {len(collapsed)}  (NEW vs old-37: {n_new})")

    print("\n=== stable but no FE within cut (PASS_NOFE) ===", flush=True)
    for r in rows_out:
        if r["stable"] and r["tag"] != "MATCH":
            print(f"  l/b={r['lob']} {r['sym']}  O={r['O0']:.4f}  "
                  f"nearest FE={r['Ofe']} miss={r['miss']:.1f}%")

    print(f"\n=== FE in [{DISC_LO},{DISC_HI}] with no unique MATCH "
          "(after n_cpair=3 discovery upgrade) ===", flush=True)
    claimed = set()
    for r in collapsed.values():
        claimed.add((r["lob"], r["sym"], round(r["Ofe"], 5)))
    n_uncovered = 0
    for (lob, sym), lst in L.IP_FE_LISTS.items():
        sl = "SYM" if sym else "ANTI"
        for Ofe in lst:
            if Ofe < DISC_LO - 1e-6 or Ofe > DISC_HI + 1e-6:
                continue
            if (lob, sl, round(Ofe, 5)) not in claimed:
                print(f"  l/b={lob} {sl}  FE={Ofe:.5f}  UNCOVERED")
                n_uncovered += 1
    print(f"  UNCOVERED count: {n_uncovered}")
    print("end", time.strftime("%Y-%m-%d %H:%M:%S"), flush=True)


if __name__ == "__main__":
    t0 = time.time()
    main()
    print(f"wall time: {time.time() - t0:.1f}s", flush=True)
