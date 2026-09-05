# -*- coding: utf-8 -*-
"""
probe_ff_r200_a50_capture_nu030.py -- VALIDATION-RUN CAPTURE, not a
monkeypatch. Writes nothing to the package, no SOLVER_VERSION impact.

Rank 16 of PAPER1_CAVEATS: full one-to-one adjudication exists at
three points that only vary one axis at a time (FF-P1 r150/a50,
r150/a100, r200/a100). This run is the first crossed combination:
  r0/2b=2.0  (Ri/Ro=0.600, same ratio as r200/a100)
  2Θ/π=0.5   (PHI=90deg, same angle as FF-P1)
  ν=0.30     (same material as all three published adjudications)

Plain production detector (no SIGMIN_LOCAL_SPLIT). Job 2410036 is
why. ONE PART PER JOB (PART=1 OOP / PART=2 IP). Resume-skip if
OUT_JSON is already complete=true.

This capture does not close rank 16 by itself. ANSYS decks + merge
+ probe_ff_adjudicate_candidates.py run after both JSONs land.
"""
import json
import os
import sys
import time

sys.path.insert(0, os.environ.get("PKG_PATH", "."))
os.environ.setdefault("DPS", "40")

from mpmath import mp
mp.dps = int(os.environ["DPS"])

import plate_solver as ps
from plate_solver import (OutOfPlaneSolver, InPlaneSolver, IsotropicMaterial,
                          make_geometry, FreeFreeOOP, FreeFreeIP,
                          find_modes_sigmin)
from plate_solver import detectors as det
from plate_solver.geometry import _kbar_wbar
from plate_solver.validation import _scan_cfg_part1, _scan_cfg_part2

EXPECT_VER = os.environ.get("EXPECT_SOLVER_VERSION", "2026-07-10.s10")
JOB_ID = os.environ.get("SLURM_JOB_ID", "local")

R0_2B = 2.0
TWO_T_PI = 0.5     # PHI=90deg -- FF-P1 angle, r200 ratio
NU = 0.30
E_ISO = 210e9
RHO = 7800.0

PART = os.environ.get("PART")
if PART not in ("1", "2"):
    print("FATAL: PART env var must be '1' (OOP) or '2' (IP).")
    raise SystemExit(2)
PART = int(PART)

OUT_JSON = os.environ.get(
    "OUT_JSON",
    f"r200_a50_nu030_capture_{'oop' if PART == 1 else 'ip'}.json")


def hdr(s):
    print("\n" + "=" * 72 + f"\n  {s}\n" + "=" * 72, flush=True)


def preflight():
    ok = True
    if ps.SOLVER_VERSION != EXPECT_VER:
        print(f"PREFLIGHT WARN: SOLVER_VERSION={ps.SOLVER_VERSION!r} != "
              f"{EXPECT_VER!r} (continuing; update EXPECT_SOLVER_VERSION "
              f"if intentional)")
    for fn in ("weak_enforcement_residual_oop", "weak_enforcement_residual_ip"):
        if not hasattr(det, fn):
            print(f"PREFLIGHT FAIL: detectors.{fn} missing")
            ok = False
    if not ok:
        raise SystemExit(2)
    print(f"preflight OK: SOLVER_VERSION={ps.SOLVER_VERSION}", flush=True)


def main():
    hdr("probe_ff_r200_a50_capture_nu030 -- start "
        + time.strftime("%Y-%m-%d %H:%M:%S")
        + f"  PART={PART} ({'OOP' if PART == 1 else 'IP'})")

    if os.path.exists(OUT_JSON):
        with open(OUT_JSON) as fh:
            existing = json.load(fh)
        if existing.get("complete"):
            print(f"  {OUT_JSON} already exists and is marked complete "
                  f"-- RESUME-SKIP, not re-running. Delete the file first "
                  f"if a genuine re-run is intended.", flush=True)
            return
        print(f"  {OUT_JSON} exists but is NOT marked complete -- "
              f"re-running the full scan from scratch.", flush=True)

    preflight()

    mat = IsotropicMaterial(E=E_ISO, nu=NU, rho=RHO)
    geom = make_geometry(R0_2B, TWO_T_PI)
    ri, ro = float(geom.R_i), float(geom.R_o)
    print(f"  Ri={ri} Ro={ro} Ri/Ro={ri/ro:.4f}  (expect 6/10 = 0.6)",
          flush=True)
    if abs(ri / ro - 0.6) > 1e-3:
        print("FATAL: geometry Ri/Ro != 0.6")
        raise SystemExit(2)
    k_bar, w_bar = _kbar_wbar(geom, mat)

    if PART == 1:
        solver = OutOfPlaneSolver(geom, mat, M=80, n_quad=30,
                                   boundary=FreeFreeOOP())
        (lo, hi), ns, cutoffs, xmax = _scan_cfg_part1(R0_2B, TWO_T_PI, mat,
                                                       n_modes_wanted=6)
    else:
        solver = InPlaneSolver(geom, mat, M=80, n_quad=30,
                                boundary=FreeFreeIP())
        (lo, hi), ns, cutoffs, xmax = _scan_cfg_part2(R0_2B, TWO_T_PI, mat,
                                                       n_modes_wanted=6)

    hdr(f"{'OOP' if PART == 1 else 'IP'} free-free production scan "
        f"(default path, n_dofs=20, window auto-derived from cutoffs)")
    print(f"  geometry r0/2b={R0_2B} 2Theta/pi={TWO_T_PI} (PHI="
          f"{TWO_T_PI * 180:.0f}deg)  material nu={NU} E={E_ISO:g} "
          f"rho={RHO}", flush=True)
    print(f"  window ({lo:.4f}, {hi:.4f}) n_scan={ns} xmax={xmax} "
          f"cutoffs={[f'{c:.4f}' for c in cutoffs]}", flush=True)
    print("  SIGMIN_LOCAL_SPLIT is NOT set (plain production default)",
          flush=True)

    t0 = time.time()
    diag = {}
    freqs = find_modes_sigmin(solver, part=PART, Omega_range=(lo, hi),
                              n_scan=ns, n_dofs=20, max_dim=xmax,
                              n_modes_wanted=60, verbose=True, diag=diag)
    dt_min = (time.time() - t0) / 60
    print(f"\n  detection done in {dt_min:.1f} min; {len(freqs)} accepted.",
          flush=True)
    print(f"  diag candidates: {diag.get('candidates', [])}", flush=True)
    print(f"  diag dropped: {diag.get('dropped', [])}", flush=True)

    capture = dict(
        complete=True, job_id=JOB_ID, part=PART,
        solver_version=ps.SOLVER_VERSION, dps=mp.dps,
        r0_2b=R0_2B, two_T_pi=TWO_T_PI, nu=NU, E=E_ISO, rho=RHO,
        window=[float(lo), float(hi)], n_scan=ns, n_dofs=20,
        xmax=float(xmax), k_bar=float(k_bar), w_bar=float(w_bar),
        omegas=[float(f) for f in freqs],
        diag_candidates=[float(c) for c in diag.get("candidates", [])],
        diag_dropped=[[float(d[0]), float(d[1]), str(d[2])]
                      for d in diag.get("dropped", [])],
        wall_min=dt_min,
    )
    tmp_name = OUT_JSON + ".tmp"
    with open(tmp_name, "w") as fh:
        json.dump(capture, fh, indent=1)
    os.replace(tmp_name, OUT_JSON)
    print(f"\n  wrote {OUT_JSON} (complete=true)", flush=True)
    print("\nEnd: " + time.strftime("%Y-%m-%d %H:%M:%S"), flush=True)


if __name__ == "__main__":
    t0 = time.time()
    main()
    print(f"\ntotal wall time: {(time.time() - t0)/60:.1f} min", flush=True)
