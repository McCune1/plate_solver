# plate_solver

Exact free-vibration solver for annular- and rectangular-sector plates, and
the companion code/data for:

> G. McCune, D. Stutts, "Exact Wave-Function Solution and Spurious-Mode
> Discrimination for Completely Free Annular Sector Plates"
> (draft — see [`paper/`](paper/)).
> Mechanical and Aerospace Engineering, Missouri University of Science and
> Technology, Rolla, MO 65409, USA.

Natural frequencies are found as zeros of a characteristic boundary
determinant built from exact radial/transverse wave functions — Frobenius
series for the annular geometry, closed-form for the rectangular one — rather
than an approximate (Rayleigh-Ritz/FE) discretization. The method reproduces
the clamped-free ("cantilever") results of Seok, Tiersten (and Scarton),
*J. Sound Vib.* 271 (2004), and extends it to the completely-free (FFFF)
boundary condition, which those papers did not treat. All computation is
`mpmath` arbitrary-precision arithmetic; frequencies are found from
`det K(Ω) = 0` from scratch — published/FE values are used only for post-hoc
accuracy reporting, never to seed a search.

## Status

- **Cantilever (Phase 1), replicating the source papers — closed.** Isotropic
  annular cantilever: 29/29 tabulated frequencies reproduced. Rectangular
  cantilever, both out-of-plane and in-plane, both symmetry classes.
- **Free-free (Phase 2) — closed for the annular geometry at its primary
  validated benchmark.** Validated against independent finite-element
  (ANSYS) benchmarks: 8/8 isotropic flexural targets to 0.64%, 26 solver
  roots covering 27 extensional FE modes (worst unique-partner 0.052%),
  8/8 orthotropic flexural targets to 1.08% (material adapted from Shi,
  Liang, Wang and Teng 2016; no agreement with their printed row is
  claimed). A broader generalization check (§6.4) then swept 4 radius
  ratios × 6 sector angles × 2 motion types (48 geometries) against a
  second, independent FE benchmark: 1256/1475 candidates (85.2%) matched
  within 3% (mean 0.574%, worst 2.990%). Job 2328808's 1231/1474 (83.5%)
  is a superseded recapture of a churning candidate set; see
  `validation/geomsweep/` and `ansys/geomsweep/`. Rectangular free-free is
  an open item — the
  rectangular assemblers need a new boundary-condition abstraction before
  that work can start.
- **Phase 3** (functionally graded materials, in-plane orthotropy, MIMO
  excitation, shock response) — partially deployed. In-plane orthotropy is
  implemented in `plate_solver/` and externally validated against Q. Wang,
  D. Shi, Q. Liang and F. Ahad, *Appl. Math. Model.* 40 (2016) 9228–9253,
  Table 4 (8/8 targets, max error 0.649%, under that table's own stated
  material `E_r = 20 E_θ`). This is the comparison reported as §6.8 of the
  paper. An earlier run of the same comparison used `E_r = 40 E_θ` — the
  material from a different table in the same source — and reported 0.737%;
  that run is superseded, and the log kept under
  `validation/phase3/probe_ip_orthotropy_wang_table4_postfix_v2_2334289.out`
  is from it, not from the result the paper cites.
  Radially-graded (FGM) material support — out-of-plane T/R/κ grading and,
  newly, in-plane grading — is also implemented and passes this project's
  own regression/self-consistency gates (9/9 clean), but has no external
  literature or finite-element benchmark yet; treat it as self-consistent,
  not externally validated. A general-BC clamped-edge orientation fix is
  implemented for the out-of-plane free-free case. MIMO excitation and
  shock response remain undeployed. None of this is reflected in the
  paper, which scopes only Phases 1–2.

This repository is a **curated snapshot for reproducibility and reference**,
not a live mirror of the research working directory. It contains the
package, its test suite, the scripts and evidence needed to reproduce the
paper's validated tables and figures, and the paper itself. It does not
contain the full exploratory/debugging history (see
[What's not here](#whats-not-here)).

## Repository layout

```
plate_solver/    the package: geometry, boundary conditions, dispersion
                 relations, K-matrix assembly, root detection, plotting,
                 validation, CLI
tests/           unit tests (unittest, no pytest dependency)
validation/      reproduction scripts + evidence backing the paper's
                 validated tables and specific claims (see validation/README.md)
scripts/         figure-generation scripts for the paper's figures
cluster/         SLURM launchers for the full overnight batch (HPC-specific;
                 adapt the SBATCH header to your own scheduler)
data/checkpoints/  checkpoint JSON files consumed by the figure scripts
ansys/           ANSYS APDL decks and text results used as FE benchmarks
paper/           CANONICAL Paper 1 manuscript: LaTeX + compiled PDF of the
                 draft and of its Supplementary Material. See paper/README.md.
```

Each subfolder has its own README with more detail.

### Before rerunning any `submit_*.sh`

All 41 SLURM submit scripts in this repository (`cluster/`, `scripts/`,
`ansys/`, `ansys/geomsweep/`, `validation/`, `validation/geomsweep/`) were
written for one specific cluster deployment and hardcode

```
/home/ghmkfh/PythonMill/Plate_Solver_Package
```

as the package root, plus a matching `venv/bin/activate` path and an
`#SBATCH` header for that scheduler. Change all three before use. The
Python scripts they invoke have no such dependency and run anywhere.

## Installation

```bash
pip install -r requirements.txt
```

Requires `numpy`, `scipy`, `mpmath`. `matplotlib` is needed only for
`plate_solver.plotting`; `openpyxl` only for the Excel export in
`plate_solver.io_export`; `sympy` only for the two standalone
`validation/e2_parity_check_step*.py` algebra checks.

Versions are unpinned because the package uses only long-stable APIs from
each dependency. The results in the paper were produced under Python
3.12.1, with `mpmath` supplying every arbitrary-precision evaluation.
`numpy`/`scipy` are used for the float64 dispersion-root search that
selects which radial branches enter the basis, so a different linear-algebra
build can in principle change a branch count at the margin; the
`tests/test_solver.py` gates exist to catch exactly that, and should be run
before trusting a reproduction against the tabulated values.

## Quickstart

```bash
# Self-tests, cutoff-frequency sanity checks, sigma_min(K) regression gates,
# worker-path reconstruction — fast, no full mode search:
python3 -m unittest tests.test_solver -v

# One-off validation runs against the source papers / FE benchmarks:
python3 -m plate_solver.cli --shi        # Shi/Liang/Wang/Teng 2016 orthotropic FFFF
python3 -m plate_solver.cli --rect       # rectangular cantilever, out-of-plane
python3 -m plate_solver.cli --rect-ip    # rectangular cantilever, in-plane

# Full overnight batch (annular cantilever regression + free-free sweep +
# literature comparison). This is a multi-hour, multi-core job — do not run
# it on a laptop or in a constrained container; see cluster/ for the SLURM
# launcher this was actually run with.
python3 -m plate_solver.cli
```

`SOLVER_VERSION` in `plate_solver/config.py` is a checkpoint-compatibility
key: any checkpoint JSON in `data/checkpoints/` is only valid for the solver
version it was produced with (recorded inside the JSON).

## Reproducing the paper's evidence

- **The validated frequency tables** (§6 of the paper): run
  `validation/test_validated_tables.py` (wrapped by
  `validation/submit_validated_tables_gate.sh` on a SLURM cluster). Its
  clean-pass output is archived at `validation/tables_gate_2318785.out`.
- **The mode-ambiguity / MAC-correlation claim** (§6.3): rerun
  `validation/probe_mac_ambiguity_check.py` against the ANSYS eigenvector
  export in `ansys/`.
- **The supplementary spurious-zero table** (`tab:oop-spurious-sm`,
  Table S.8): rerun `validation/probe_ff_oop_spurious_table.py` for the
  eight physical modes and five in-window spurious zeros
  (`validation/ff_oop_spurious_table_2325529.out`). The sixth Table S.8
  row (raw Ω=2.817613, Fig. 2(e)) is `validation/ff_oop_artifact_2p8176_2414022.out`.
- **The geometry generalization sweep** (§6.4, Table S.3): the python
  candidate side is reproduced by
  `validation/geomsweep/consolidate_geom_sweep_results.py` against the
  manifests in `validation/geomsweep/figures_ff_*/`; the independent FE
  side is in `ansys/geomsweep/`. The **live** §6.4 headline
  (1256/1475, 85.2%, mean 0.574%, worst 2.990%) is job 2416746:
  `validation/geomsweep/geomsweep_fe_match_s10v2_2416746.out` plus
  `_s10v2_consolidated.csv` / `_s10v2_manifests.zip`. The probe
  `probe_geomsweep_fe_match.py` now defaults to that zip. The older
  `figures_ff_*` tree is the superseded 2328808 recapture (1231/1474);
  set `GEOMSWEEP_MANIFEST_DIR` to `validation/geomsweep` to replay it.
  See `validation/geomsweep/README.md`. The
  cut-off-adjacent residual screen is fully reproducible via the three
  `probe_ff_*.py` scripts in `validation/geomsweep/`.
- **The paper's embedded figures**: `scripts/probe_paper_figures_oop.py`
  regenerates the 5 panels in `paper/paper_figures/` from a checkpoint in
  `data/checkpoints/`. `scripts/probe_regen_dispersion_figs.py` and
  `scripts/probe_regen_modeshape_figs.py` regenerate the broader
  dispersion/mode-shape galleries beyond what's embedded in the paper.
- **The finite-element benchmarks**: `ansys/` holds the APDL input decks
  (`*.inp`) and their text results (`*_out.txt`) used as independent
  validation targets throughout.

## References

The exact-wave-function method this package implements and extends:

- J. Seok, H.F. Tiersten, "Free vibrations of annular sector cantilever
  plates, Part 1: out-of-plane motion," *J. Sound Vib.* 271 (2004) 757–772.
- J. Seok, H.F. Tiersten, "Free vibrations of annular sector cantilever
  plates, Part 2: in-plane motion," *J. Sound Vib.* 271 (2004) 773–787.
- J. Seok, H.F. Tiersten, H.A. Scarton, "Free vibrations of rectangular
  cantilever plates, Part 1: out-of-plane motion," *J. Sound Vib.* 271
  (2004) 131–146.
- J. Seok, H.F. Tiersten, H.A. Scarton, "Free vibrations of rectangular
  cantilever plates, Part 2: in-plane motion," *J. Sound Vib.* 271 (2004)
  147–158.
- D. Shi, Q. Liang, Q. Wang, X. Teng, "A unified solution for free vibration
  of orthotropic circular, annular and sector plates with general boundary
  conditions," *J. Vibroengineering* 18 (2016) 3138–3152.

These are copyrighted journal articles and are not redistributed in this
repository. The full bibliography for the extensions and background
literature is in `paper/PAPER1_FREEFREE_DRAFT.tex`.

## What's not here

This repository intentionally excludes the internal research working
directory's day-to-day debugging history: the bug ledger and closure log
(`LESSONS_LEARNED.md`), external peer-review transcripts, and the archive of
retired probe scripts from closed investigations. That material documents
*how* the validated results were reached, not what's needed to run the
solver, check its tests, or reproduce the paper's tables and figures. It
also excludes the reference PDFs (copyright) and ANSYS binary scratch/result
files (`file.db`/`.full`/`.rst`/etc. — large, unreadable without a licensed
ANSYS install, and regenerable by rerunning the included `.inp` decks).

One consequence: a number of in-code comments cite `LESSONS_LEARNED.md`
section numbers, internal review-document items, or SLURM job numbers, and
a few `submit_*.sh` comments mention sibling scripts from that closed-out
work. None of that is required to run or understand the code in this
repository — every substantive finding it points to is either already
reflected in the code's behavior or written up in the paper — but readers
diffing comments against what's actually here shouldn't expect those
references to resolve to anything included.

## Citing

See [`CITATION.cff`](CITATION.cff).

## License

[MIT](LICENSE).
