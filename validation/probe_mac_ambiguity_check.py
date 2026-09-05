# -*- coding: utf-8 -*-
"""
probe_mac_ambiguity_check.py -- DIAGNOSTIC ONLY. Writes only a small text
report; no package changes, no SOLVER_VERSION implications.

CONTEXT: the FE mode at Omega_lit=111.4627 (f=34.8186 Hz) has TWO
parity-consistent solver candidates (raw Omega=2.828013, 0.16% proximity;
raw Omega=2.838515, 0.52% proximity); five independent techniques (basis
bump, sensitivity map, null-vector gap, parity readout, ODD-block topology
scan -- jobs 2316969/2317010/2317114/2317684/2318031) could not force the
assignment between them, and it was closed as a bounded, documented
ambiguity rather than resolved by fiat (LESSONS_LEARNED Sec 20.4). None of
those five techniques compared mode SHAPE -- each worked with frequency
proximity or internal matrix/basis structure only. This probe is the
identified sixth, genuinely new lever: a Modal Assurance Criterion (MAC)
correlation between each candidate's reconstructed solver mode shape and
the FE eigenvector extracted by ansys_ff_oop_mac_extract.inp (companion
deck, run separately by the user with a licensed ANSYS install -- this
probe is pure Python/SLURM and does not touch ANSYS).

INPUT REQUIRED: mac_eigvec_mesh96.txt (the converged-mesh eigenvector dump
from ansys_ff_oop_mac_extract.inp's Pass 2). As of 2026-07-15, the Ansys
submit script (Ansys/submit_ansys_ff_oop_mac_extract.sh) auto-copies this
file up to the package root right after the Ansys job finishes, so it
should already be sitting next to this script with no manual step -- see
LESSONS_LEARNED Sec 13.5. If FE_EIGVEC_PATH is set explicitly, only that
exact path is tried (no fallback). Otherwise this probe searches, in
order: (1) mac_eigvec_mesh96.txt at the package root (the normal location
after the auto-copy), (2) Ansys/mac_eigvec_mesh96.txt (fallback, in case
the auto-copy didn't run -- older submit script, or the copy step failed).
Columns: r, theta_deg, UX, UY, UZ, one row per FE node, theta_deg in
[0, 90] (the FE mesh's own convention, PHI=0 at one radial edge). This
probe does NOT run if the file can't be found anywhere -- it cannot
substitute for it.

COORDINATE RECONCILIATION: the solver's angular convention is
theta in [-Theta, +Theta] with Theta=pi/4 (FF-P1: 2*Theta/pi=0.5), centered
on the sector bisector; the FE mesh instead spans theta_deg in [0, PHI]
with PHI=90 deg, i.e. the SAME physical sector but referenced from one
radial edge instead of the bisector. Converted here as
theta_solver = theta_deg*pi/180 - Theta -- a rigid re-labeling of the same
points, not a coordinate approximation.

METHOD: for each of the two candidates, reconstruct the solver's OOP mode
shape on its native (r, theta) grid via OutOfPlaneSolver.mode_shape_grid
(the exact function plot_mode_shape/Fig-3 already uses in production --
not a new reconstruction path), build a bilinear interpolator over that
structured grid, evaluate it at every FE node's (r, theta_solver), and
compute the Modal Assurance Criterion
    MAC = |sum_k(w_fe_k * w_solver_k)|^2 / (sum_k(w_fe_k^2) * sum_k(w_solver_k^2))
against the FE UZ vector at those same nodes. MAC is scale- and sign-
invariant by construction (a mode shape's absolute amplitude/sign is
arbitrary), so no separate normalization step is needed beyond the
formula itself.

PRE-REGISTERED INTERPRETATION (before this probe has been run against
real FE data, per the project's standing methodology):
  One candidate's MAC is high (loosely, > 0.8-0.9) and the other's is
    markedly lower -> the high-MAC candidate is the FE mode's true
    solver counterpart; record the assignment and update
    LESSONS_LEARNED Sec 20.4 accordingly. This DOES settle the
    ambiguity, since mode-shape correlation is independent of every one
    of the five previously-exhausted frequency/matrix-structure
    techniques.
  Both candidates' MAC are similarly high, or both similarly low/
    moderate -> MAC does not discriminate here either; report the
    numbers honestly and leave the ambiguity open. Do not round a
    marginal MAC gap up to a decision.
  Either FE UZ-dominance check (also re-verified here, redundantly with
    the APDL deck's own check) fails, or the two mesh passes
    (mesh64/mesh96, if both are supplied) disagree on which mode was
    extracted -> the FE eigenvector itself is not trustworthy yet; stop
    and fix the extraction deck before drawing any conclusion from MAC
    numbers computed against it.
This probe does not, either way, reopen the closed Sec 20.4 ambiguity
record by itself -- it only supplies the one new, previously-untried
diagnostic; a human should read the printed numbers against the
pre-registered criteria above before updating any project document.

COST: cheap. Two mode_shape_grid reconstructions (each one full_search +
one grid build at fixed Omega, ~90-150s, same cost driver as
probe_paper_figures_oop.py) plus a fast numpy interpolation/dot-product
over however many FE nodes there are (thousands at most -- trivial). No
search, no golden section, no discovery. Single CPU is enough; a small
worker count is used only so the two candidates run in parallel.
"""
import os
import sys
import time

sys.path.insert(0, os.environ.get("PKG_PATH", "."))
os.environ.setdefault("DPS", "40")

import numpy as np

EXPECT_VER = os.environ.get("EXPECT_SOLVER_VERSION", "2026-07-10.s10")
# N_DOFS/XI_MAX env-configurable (added 2026-07-15 after job 2324018):
# defaults (20, 20.0) reproduce the original run. LESSONS_LEARNED Sec 20.4's
# next step is a larger-basis re-check (n_dofs=28, maybe 36) of the marginal
# MAC gap (0.9913 vs 0.9998) between the two ODD candidates -- set
# MAC_N_DOFS=28 (or 36) AND MAC_XI_MAX accordingly in the submit script's
# environment to run that without editing this file again. BOTH must move
# together (project's standard "n_dofs+8/xi_max+6" convention, per
# CLAUDE.md's methodology section) -- job 2324023 confirmed why: N_DOFS=28
# alone with the old XI_MAX=20.0 only found enough branches to fill 24/28
# dofs and hard-failed (select_fill needs enough roots inside [0, xi_max] to
# supply the requested dof count; it does not silently proceed short). A
# widening MAC gap at the larger basis would support cand_B (raw 2.838515);
# a flat gap would support genuine near-degeneracy rather than a resolvable
# assignment.
N_DOFS = int(os.environ.get("MAC_N_DOFS", "20"))
XI_MAX = float(os.environ.get("MAC_XI_MAX", "20.0"))
THETA = float(np.pi) / 4.0  # FF-P1: 2*Theta/pi = 0.5 -> Theta = pi/4

CANDIDATES = [
    ("cand_A_2p828013", 2.828013, 0.16),
    ("cand_B_2p838515", 2.838515, 0.52),
]

# FE_EIGVEC_PATH: if the env var is set explicitly, that exact path is the
# ONLY one tried (no fallback -- an explicit override means the caller knows
# better). Otherwise this probe searches FE_EIGVEC_CANDIDATES in order; see
# main() and LESSONS_LEARNED Sec 13.5 for why there's more than one candidate.
FE_EIGVEC_PATH_ENV = os.environ.get("FE_EIGVEC_PATH")
FE_EIGVEC_CANDIDATES = [
    "mac_eigvec_mesh96.txt",                          # package root: normal
                                                       # location after the Ansys
                                                       # submit script's auto-copy
    os.path.join("Ansys", "mac_eigvec_mesh96.txt"),   # fallback: Ansys subfolder,
                                                       # in case the auto-copy
                                                       # didn't run or wasn't done
]
FE_EIGVEC_PATH = FE_EIGVEC_PATH_ENV  # resolved for real in main(); may be None here


def load_fe_eigvec(path):
    """Reads r, theta_deg, UX, UY, UZ from the APDL *VWRITE dump (one
    header line, then whitespace-separated columns). Returns arrays."""
    rows = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or not line[0].isdigit():
                continue  # skip header / blank lines
            parts = line.split()
            if len(parts) < 5:
                continue
            try:
                rows.append([float(x) for x in parts[:5]])
            except ValueError:
                continue
    if not rows:
        raise RuntimeError(f"no numeric data rows parsed from {path}")
    arr = np.array(rows)
    return arr[:, 0], arr[:, 1], arr[:, 2], arr[:, 3], arr[:, 4]


def reconstruct_candidate(label, Om):
    """Runs in-process (not a worker pool -- see module docstring, this is
    cheap and simple enough not to need one; kept single-process for a
    clear, linear log)."""
    t0 = time.time()
    import plate_solver as ps
    from plate_solver.detectors import full_search, select_fill

    mat = ps.IsotropicMaterial(E=210e9, nu=0.30, rho=7800.0)
    geom = ps.make_geometry(1.5, 0.5)
    solver = ps.OutOfPlaneSolver(geom, mat, M=80, n_quad=30,
                                  boundary=ps.FreeFreeOOP())
    raw = full_search(solver.fast, Om, xmax=XI_MAX)
    sel, cnt = select_fill(raw, N_DOFS)
    if cnt < N_DOFS:
        raise RuntimeError(f"{label}: only {cnt}/{N_DOFS} branches filled")
    grid = solver.mode_shape_grid(Om, sel, N_DOFS, n_r=41, n_th=41)
    if grid is None:
        raise RuntimeError(f"{label}: mode_shape_grid reconstruction failed")
    r_phys, th_phys, Wgrid = grid
    dt = time.time() - t0
    print(f"  {label}: reconstructed on ({len(r_phys)}x{len(th_phys)}) grid "
          f"[{dt:.0f}s]", flush=True)
    return np.asarray(r_phys, dtype=float), np.asarray(th_phys, dtype=float), \
        np.asarray(Wgrid, dtype=float)


def mac(a, b):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    num = abs(float(np.dot(a, b))) ** 2
    den = float(np.dot(a, a)) * float(np.dot(b, b))
    if den <= 0:
        return float("nan")
    return num / den


def main():
    print("=" * 78)
    print("  probe_mac_ambiguity_check -- start "
          + time.strftime("%Y-%m-%d %H:%M:%S"))
    print(f"  N_DOFS={N_DOFS} (MAC_N_DOFS env var; default 20)   "
          f"XI_MAX={XI_MAX} (MAC_XI_MAX env var; default 20.0)")
    print("=" * 78, flush=True)

    import plate_solver as ps
    if ps.SOLVER_VERSION != EXPECT_VER:
        print(f"PREFLIGHT FAIL: SOLVER_VERSION={ps.SOLVER_VERSION!r} != "
              f"{EXPECT_VER!r}")
        raise SystemExit(2)
    print(f"preflight OK: SOLVER_VERSION={ps.SOLVER_VERSION}", flush=True)

    global FE_EIGVEC_PATH
    if FE_EIGVEC_PATH_ENV:
        FE_EIGVEC_PATH = FE_EIGVEC_PATH_ENV
        if not os.path.isfile(FE_EIGVEC_PATH):
            print(f"FATAL: FE eigenvector file not found at explicit "
                  f"FE_EIGVEC_PATH={FE_EIGVEC_PATH!r}. Cannot proceed without "
                  f"it -- an explicit FE_EIGVEC_PATH has no fallback.",
                  flush=True)
            raise SystemExit(3)
        print(f"Using explicit FE_EIGVEC_PATH={FE_EIGVEC_PATH!r}", flush=True)
    else:
        found = None
        for cand in FE_EIGVEC_CANDIDATES:
            if os.path.isfile(cand):
                found = cand
                break
        if found is None:
            print(f"FATAL: FE eigenvector file not found at any of "
                  f"{FE_EIGVEC_CANDIDATES}. Run ansys_ff_oop_mac_extract.inp "
                  f"first (from Ansys/, on an ANSYS-licensed allocation) -- "
                  f"its submit script auto-copies mac_eigvec_mesh96.txt up to "
                  f"the package root, so no manual copy step should be needed "
                  f"(see LESSONS_LEARNED Sec 13.5). If it's still missing "
                  f"after that job completed, check its /COM lines and .out "
                  f"log for a mode-selection or UZ-dominance failure -- the "
                  f"file may simply never have been produced. Cannot proceed "
                  f"without it -- there is no further fallback.", flush=True)
            raise SystemExit(3)
        FE_EIGVEC_PATH = found
        print(f"Found FE eigenvector at {FE_EIGVEC_PATH!r} "
              f"(searched: {FE_EIGVEC_CANDIDATES})", flush=True)

    r_fe, th_deg_fe, ux_fe, uy_fe, uz_fe = load_fe_eigvec(FE_EIGVEC_PATH)
    n_fe = len(r_fe)
    print(f"Loaded FE eigenvector: {n_fe} nodes from {FE_EIGVEC_PATH}",
          flush=True)

    uz_max = float(np.max(np.abs(uz_fe)))
    uinplane_max = float(np.max(np.sqrt(ux_fe ** 2 + uy_fe ** 2)))
    print(f"FE UZ-dominance re-check: max|UZ|={uz_max:.4g}  "
          f"max in-plane={uinplane_max:.4g}  "
          f"ratio={uz_max / max(uinplane_max, 1e-30):.2f}x", flush=True)
    if uinplane_max > 0 and uz_max / uinplane_max < 5.0:
        print("  WARNING: UZ does not clearly dominate -- this may not be "
              "a clean OOP mode. Do not trust the MAC numbers below without "
              "checking the extraction deck's mode selection.", flush=True)

    th_solver_fe = th_deg_fe * np.pi / 180.0 - THETA

    results = {}
    for label, Om, prox in CANDIDATES:
        r_phys, th_phys, Wgrid = reconstruct_candidate(label, Om)

        # Self-parity check (added 2026-07-15 after job 2324002 returned
        # MAC=0.0000 for BOTH candidates -- an exact zero, not just "low", is
        # the signature you'd get for free if the reconstructed shape has the
        # WRONG parity: an odd (antisymmetric-about-bisector) FE field dotted
        # against an even (symmetric) solver field integrates to ~0 by
        # symmetry alone, independent of basis size or interpolation. The raw
        # FE dump itself is confirmed odd by hand (r=8: UZ(theta=0)=-1.0,
        # UZ(theta=90)=+1.0 -- an exact antisymmetric pair about the 45-deg
        # bisector). This checks the RECONSTRUCTION's own parity the same
        # way, using only the already-computed Wgrid -- no extra full_search
        # cost -- so a parity mismatch is caught before blaming the MAC math.
        th_mid = Wgrid.shape[1] // 2
        left = Wgrid[:, :th_mid].ravel()
        right_mirror = Wgrid[:, -1:-th_mid - 1:-1].ravel()
        odd_corr = float(np.corrcoef(left, -right_mirror)[0, 1])
        even_corr = float(np.corrcoef(left, right_mirror)[0, 1])
        print(f"  {label}: self-parity check -- corr(W(-th),-W(+th))="
              f"{odd_corr:.4f} [near +1 => reconstruction is ODD, matching "
              f"the FE target]   corr(W(-th),+W(th))={even_corr:.4f} "
              f"[near +1 => reconstruction is EVEN -- WRONG parity for this "
              f"target, would explain an exact-zero MAC below]", flush=True)
        # RESOLVED 2026-07-15 (job 2324010 -> B18, LESSONS_LEARNED bug ledger):
        # both candidates came back even_corr=1.0000 here, confirming the
        # hypothesis above. Root cause found in core_solvers.py's
        # mode_shape_grid: it used cos(xi*theta+ph) where _build_K_real's own
        # validated edge assembly uses sin(...) for the displacement
        # projection (and real_basis_items' docstring says q=0 is the SIN
        # dof) -- a 90-degree phase error, not a sign flip, that silently
        # gave every mode-shape figure the wrong theta-parity. Fixed in
        # mode_shape_grid; this probe needed no change since it only calls
        # that function -- rerun after syncing the fixed core_solvers.py to
        # see corrected self-parity/MAC numbers.

        from scipy.interpolate import RegularGridInterpolator
        interp = RegularGridInterpolator((r_phys, th_phys), Wgrid,
                                          bounds_error=False, fill_value=None)
        pts = np.column_stack([r_fe, th_solver_fe])
        w_solver_at_fe = interp(pts)
        bad = ~np.isfinite(w_solver_at_fe)
        if bad.any():
            print(f"  {label}: {bad.sum()}/{n_fe} FE nodes fell outside the "
                  f"solver grid (extrapolated as NaN) -- excluded from MAC",
                  flush=True)
        m = mac(uz_fe[~bad], w_solver_at_fe[~bad])
        results[label] = dict(Om=Om, proximity_pct=prox, MAC=m,
                               n_used=int((~bad).sum()))
        print(f"  {label}: Om={Om:.6f} (freq. proximity {prox}%)  "
              f"MAC={m:.4f}  (n={int((~bad).sum())}/{n_fe} nodes used)",
              flush=True)

    print("\n" + "=" * 78)
    print("  SUMMARY")
    print("=" * 78, flush=True)
    for label, res in results.items():
        print(f"  {label}: Om={res['Om']:.6f}  MAC={res['MAC']:.4f}",
              flush=True)
    macs = {k: v["MAC"] for k, v in results.items()}
    if all(np.isfinite(list(macs.values()))):
        best = max(macs, key=macs.get)
        worst = min(macs, key=macs.get)
        gap = macs[best] - macs[worst]
        print(f"\n  Highest MAC: {best} ({macs[best]:.4f}); "
              f"gap over the other candidate: {gap:.4f}", flush=True)
        if macs[best] > 0.8 and gap > 0.2:
            print(f"  READING (per pre-registered criteria): {best} is the "
                  f"FE mode's true solver counterpart. Record this in "
                  f"LESSONS_LEARNED Sec 20.4 and the two paper/report docs "
                  f"if a human confirms the UZ-dominance and mesh64/mesh96 "
                  f"consistency checks above look clean.", flush=True)
        else:
            print(f"  READING (per pre-registered criteria): MAC does NOT "
                  f"cleanly discriminate here either (best={macs[best]:.4f}, "
                  f"gap={gap:.4f}). Report honestly; the ambiguity stays "
                  f"open.", flush=True)

    print("\nDiagnostic only -- no SOLVER_VERSION action, no package "
          "changes. Does not by itself reopen or re-decide Sec 20.4; a "
          "human should read this against the pre-registered criteria in "
          "this script's own docstring.", flush=True)
    print("End: " + time.strftime("%Y-%m-%d %H:%M:%S"), flush=True)


if __name__ == "__main__":
    t0 = time.time()
    main()
    print(f"\ntotal wall time: {(time.time() - t0)/60:.1f} min", flush=True)
