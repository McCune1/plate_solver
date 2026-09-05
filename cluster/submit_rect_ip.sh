#!/bin/bash
#SBATCH --job-name=rect_ip
#SBATCH --output=rect_ip_%j.out
#SBATCH --error=rect_ip_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --partition=general
#
# ──────────────────────────────────────────────────────────────────────────────
#  Rectangular cantilever IN-PLANE (Part 2) — INTEGRATED validation.
#  Seok, Tiersten & Scarton, JSV 271 (2004) 147–158 (rect Part 2, in-plane).
#
#  2026-07-01: updated for the plate_solver package (module split, see
#  LESSONS_LEARNED.md Sec. 22). Two real changes vs the old submit_rect_ip.sh:
#    1. Invokes `python3 -m plate_solver.cli --rect-ip` instead of
#       `rect_int.py` directly.
#    2. The annular no-op regression gate below was ALREADY known-broken in
#       the old script (LESSONS Sec. 21.1: sorted() over complex dict keys
#       raised TypeError silently swallowed by no `set -e`, AND it compared
#       raw det() magnitude, which B10/Sec.3 established is scale-free and
#       ~0 everywhere -- never a valid invariant). That fix (rebuild K at the
#       canonical fundamental, compare equilibrated log10 sigma_min against
#       the -4.59 anchor) is applied here directly rather than re-shipping
#       the broken version.
#
#  Validates the R40_RECT_IP path (engine + unconstrained assembler + driver)
#  against Table 1 at l/b ∈ {1.0, 1.5, 2.5}, both symmetry classes (SYM/ANTI).
#  Detector: equilibrated sigma_min of the unconstrained Eq.41 variational K=0.
#
#  Single-threaded mpmath (small matrices, minutes). No MPI/ProcessPool.
#  SOLVER_VERSION unchanged; annular path byte-untouched.
#
#  Deploy: copy the plate_solver/ package directory alongside this script (or
#  install it into the job's venv).
# ──────────────────────────────────────────────────────────────────────────────

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export BLIS_NUM_THREADS=1

export DPS=30                       # production precision (sandbox was dps=26)

echo "Job $SLURM_JOB_ID on $(hostname)  start: $(date)"
module load python/3.12.1
source /home/ghmkfh/PythonMill/Plate_Solver_Package/venv/bin/activate
python3 -c "import mpmath, numpy, scipy" || pip install --break-system-packages --user mpmath numpy scipy

# ── 1. RECTANGULAR IN-PLANE validation (R40_RECT_IP=1) ─────────────────────
#    Natural frequencies detected as zeros of the unconstrained det K(Ω̄) via
#    equilibrated sigma_min dips.  Table 1 (Seok 2004 Part 2) accuracy post-hoc.
#    Validates at l/b=1.0, 1.5, 2.5.  Exits after this gate (--rect-ip).
#
echo ""
echo "######################################################################"
echo "# 1) rect-IP validation vs Table 1"
echo "#    l/b in {1.0, 1.5, 2.5}, SYM+ANTI fundamentals, dps=$DPS"
echo "######################################################################"
python3 -u -m plate_solver.cli --rect-ip

# ── 2. ANNULAR regression no-op gate (equilibrated sigma_min, NOT raw det) ──
#    Confirms the rect-IP merge did not perturb the annular OOP path: rebuild
#    K at the canonical fundamental (Om=0.033238) on a FRESH full_search
#    basis (sigma_min depth is basis-sensitive -- never reuse a tracked
#    basis, see plate-solver-sandbox-probe skill Sec. 5), and check
#    log10(sigma_min) against the established anchor (-4.59, tol 0.10).
#
echo ""
echo "######################################################################"
echo "# 2) ANNULAR no-op check (equilibrated sigma_min at the fundamental)"
echo "#    rect-IP merge must not perturb the annular OOP detector"
echo "######################################################################"
python3 - <<'PYEOF'
import os
os.environ.setdefault("DPS", "30")
from mpmath import mp
import plate_solver as ps
from plate_solver.detectors import full_search, select_fill, sigma_min_from_K

geom = ps.make_geometry(1.25, 1.0)
mat = ps.IsotropicMaterial(E=210e9, nu=0.35, rho=7800.0)
oop = ps.OutOfPlaneSolver(geom, mat, M=80, n_quad=30)

Om = 0.033238
brs = full_search(oop.fast, Om, xmax=14.0)     # FRESH basis, never a tracked one
sel, n = select_fill(brs, 16)
K, sz = oop._build_K_real(Om, sel, fast_scan=False)
log_sigma = sigma_min_from_K(K, sz)

ANCHOR, TOL = -4.59, 0.10
ok = abs(log_sigma - ANCHOR) < TOL
print(f"  log10(sigma_min) at Om={Om}: {log_sigma:.4f}  (anchor {ANCHOR}, tol {TOL})")
print(f"  Annular no-op gate: {'PASS' if ok else 'FAIL -- rect-IP merge may have perturbed the annular path'}")
raise SystemExit(0 if ok else 1)
PYEOF

echo ""
echo "######################################################################"
echo "# Job complete: $(date)"
echo "######################################################################"
