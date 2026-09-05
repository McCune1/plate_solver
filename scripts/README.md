# scripts/

Figure-regeneration scripts. All read a checkpoint JSON (see
`../data/checkpoints/`, passed via the `CHECKPOINT` environment variable)
and need `plate_solver` plus `matplotlib`.

| File | Produces |
|---|---|
| `probe_paper_figures_oop.py` / `submit_paper_figures_oop.sh` | the 5 figures actually embedded in `paper/PAPER1_FREEFREE_DRAFT.tex` §6 (the σ_min trace and 4 mode-shape panels in `paper/paper_figures/`) — rerun this whenever those panels need to change |
| `probe_regen_dispersion_figs.py` / `submit_regen_dispersion_figs*.sh` | the broader per-geometry dispersion-plot gallery (not embedded in the paper; supplementary) |
| `probe_regen_modeshape_figs.py` / `submit_regen_modeshape_figs*.sh` | the broader per-geometry mode-shape gallery (not embedded in the paper; supplementary) |
| `make_graphical_abstract.py` | JSV graphical abstract (`paper/GRAPHICAL_ABSTRACT.{pdf,png}`); standalone, not embedded in the manuscript. Reads `paper/paper_figures/mac_separation_summary.txt`. |
| `build_repo_digest.py` | rebuilds `PLATE_SOLVER_REPO_DIGEST.md` at the project root from `github_repo/`. Run after any paper or validation-index change that a digest-only reviewer would need. |

The `_p1tail`/`_p2` submit-script variants point at the corresponding
`data/checkpoints/research23_checkpoint_p1tail.json` /
`_p2.json` checkpoints (the overnight batch's checkpointed continuation
runs).

The `submit_*.sh` wrappers hardcode `PKG_PATH` and a venv activation path
from the original cluster deployment (`/home/ghmkfh/PythonMill/Plate_Solver_Package`) —
update both before running them elsewhere.
