#!/usr/bin/env python3
"""
check_deployed_tree.py -- instant, dependency-light preflight: reports what
is ACTUALLY present in the deployed plate_solver/ package tree, as plain
text search (no import, no mpmath, no SLURM job) so it can be run on the
login node in under a second.

WHY: two jobs this session (2315890, 2315891) failed because documented/
assumed-landed work (the 2026-07-09 residual screen, this session's E.3
gap-refinement patch) was not actually present in the tree being run
against. SOLVER_VERSION alone cannot catch this -- both additions are
deliberately additive/non-version-bumping (correctly, per the hard
invariant), which means SOLVER_VERSION is silent about their presence. This
is the SAME failure class as job 2315355 (§13.3/A.12): work confirmed via a
monkeypatched probe copy or documented as landed, but the deploy-to-tree
step silently never ran (or ran against a tree since reset/re-cloned).

Run this BEFORE submitting any job that depends on additive/optional code,
not after a job aborts or ImportErrors.

USAGE:
    python3 check_deployed_tree.py [path/to/plate_solver]
    (defaults to ./plate_solver)
"""
import os
import re
import sys

CHECKS = [
    ("SOLVER_VERSION", "config.py",
     re.compile(
         r'(?:SOLVER_VERSION\s*=\s*["\']([^"\']+)["\']'
         r'|solver_version:\s*str\s*=\s*["\']([^"\']+)["\'])'
     )),
    ("weak_enforcement_residual_oop (2026-07-09 residual screen)",
     "detectors.py", re.compile(r"def\s+weak_enforcement_residual_oop")),
    ("weak_enforcement_residual_ip (2026-07-09 residual screen)",
     "detectors.py", re.compile(r"def\s+weak_enforcement_residual_ip")),
    ("_gap_pass parameter (Addendum 2026-07-09e E.3 gap-refinement patch)",
     "detectors.py", re.compile(r"_gap_pass\s*=\s*False")),
    ("e.sign free-edge orientation fix (B12)", "core_solvers.py",
     re.compile(r"_orient\s*=\s*\[")),
]


def main():
    pkg = sys.argv[1] if len(sys.argv) > 1 else "plate_solver"
    if not os.path.isdir(pkg):
        print(f"ABORT: {pkg} is not a directory (run from the directory "
              f"containing plate_solver/, or pass the path explicitly).")
        sys.exit(1)

    print(f"Deploy-state preflight for {os.path.abspath(pkg)}")
    print("=" * 78)
    any_missing = False
    for label, fname, pattern in CHECKS:
        fpath = os.path.join(pkg, fname)
        if not os.path.exists(fpath):
            print(f"  [MISSING FILE] {fname}  <- cannot check {label}")
            any_missing = True
            continue
        text = open(fpath, "rb").read().decode("utf-8", errors="replace")
        m = pattern.search(text)
        if m:
            ver = next((g for g in m.groups() if g), None)
            extra = f"  = {ver}" if ver else ""
            print(f"  [present] {label}{extra}")
        else:
            print(f"  [ABSENT ] {label}  ({fname})")
            any_missing = True

    print("=" * 78)
    if any_missing:
        print("Some expected additive code is ABSENT from this tree. Before "
              "submitting any job that depends on it: run the matching "
              "apply_*_patch.py / apply_*_deploy.py script against this "
              "tree first (they are idempotent -- safe even if you're not "
              "sure). Do not assume a job will fail loudly and cheaply the "
              "way 2315890/2315891 did; some absences would silently change "
              "results instead of erroring.")
    else:
        print("All checked additive code is present. Safe to proceed.")


if __name__ == "__main__":
    main()
