# validation/paper4_piezo/

Scripts and evidence behind `paper/PAPER4_PIEZO_DRAFT.tex` and its
Supplementary Material. The solver is `plate_solver/piezo_solver.py`
(`PiezoOutOfPlaneSolver`).

**Which outputs back the paper.** Every analytic number in the current
manuscript comes from `consistent_rerun_2026-09-23/`:
- the `*_consistent_2026-09-23.log` transcripts and `*_results.json` files;
- `p6_table4_consistent.json` (Table 1);
- `cc_sc_e31.log` (clamped–clamped short-circuit e31 column);
- `make_p4_scfix_figs.py` (the admittance zoom and noise-floor figures).

These were re-run on 2026-09-23 with the default
`projection='consistent'`. That is a same-weight Galerkin projection of
the electric enthalpy, with a quadratic potential in each piezo layer and
the reciprocal open-circuit charge arm (h + h1/2).

**The `*.out` files in this folder are older cluster runs.** They used
the earlier model, which the solver still provides as
`projection='duan'` and which reproduces those runs bit-for-bit:
Duan-Quek-Wang's Gauss-law closure of a sinusoidal potential, and a bus
charge taken at the outer electrode.

Four regression probes carry hard-coded open-circuit anchors. Their
anchors were moved to the new values, and the old values are kept in
comments:
- `ff_admittance_sweep`
- `ff_admittance_ratios`
- `h1_ratio_sweep`
- `lever16_freq_noise_floor`

The standalone `*_forward_model_*`, `*_cc_table4_gate_*` and
`*_elastic_*_baseline_*` scripts re-implement Duan's own model. They are
kept as the reproduction of Duan, Quek and Wang's Table 4, not as the
paper's model.

**Running.** Set `PKG_PATH` to the repository root (these scripts sit
two levels below it). The `submit_*.sh` wrappers hardcode the original
cluster paths; change them to your install.

**Finite elements.** The PLANE223 decks, their text outputs, manifests,
submit/build scripts and queue logs for every finite-element job the paper
cites are in `ansys/NewAnsys/` (`ansys_p4_*`). `ansys/README.md` maps each
job to its paper claim.
