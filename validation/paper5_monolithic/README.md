# validation/paper5_monolithic/

Scripts and evidence behind `paper/PAPER5_MONOLITHIC_PIEZO_DRAFT.tex`
and its Supplementary Material. The solver is
`plate_solver/piezo_monolithic.py` (`PiezoMonolithicOutOfPlaneSolver`).

**Which outputs back the paper.** Every number in the current manuscript
comes from `consistent_rerun_2026-09-23/`: the five `*_consistent_2026-09-23.log`
transcripts and the five `*_results.json` files. These were run with the
default `projection='consistent'`: a same-weight Galerkin projection of
the electric enthalpy with the quadratic potential φ̄(1 − z²/H²), which
is exact in the long-wave limit.

**The `*.out` files in this folder are older cluster runs.** They used
the earlier Gauss-law closure of a cosine potential, which is 12/π² high
in the long-wave limit. The solver still provides that model as
`projection='duan'`, and it reproduces those runs bit-for-bit.

The finite-element decks for the free–free check, with their text
outputs and queue log (job 2520476), are in `ansys/NewAnsys/`
(`*_p5_*`); see `ansys/README.md`.

**Segment charge, e15 term (2026-09-24).** `probe_piezo_p5_e15_segment_charge_2026-09-24.py` checks the closed-form factor κ (segment charge = κ × dielectric-flux value) against the driven 12×12 solve, then writes `piezo_p5_e15_segment_charge_predictions.json`. The FE decks in `ansys/NewAnsys/` (`ansys_p5seg_*`) are built from that file and scored against it. The sandbox log is `e15_segment_charge_sandbox_2026-09-24.log`.

**Running.** Set `PKG_PATH` to the repository root. The `submit_*.sh`
wrappers hardcode the original cluster paths; change them to your
install.
