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

**Running.** Set `PKG_PATH` to the repository root. The `submit_*.sh`
wrappers hardcode the original cluster paths; change them to your
install.
