# validation/geomsweep/

Reproduction scripts and evidence for `paper/PAPER1_FREEFREE_DRAFT.tex`
§6.4 ("Geometry generalization: radius ratio and sector angle sweep"), a
48-geometry (4 radius ratios × 6 sector angles × 2 motion types) check
that the free-free detector's accuracy against independent finite-element
benchmarks holds beyond the single geometry (FF-P1) fully adjudicated in
§6/§6.3. Material for this sweep is E=210 GPa, nu=0.35, rho=7800 kg/m^3, H=0.08 m —
a different nu from the nu=0.30 benchmark behind Tables tab:oop/tab:ip
(§sec:conversion), so this is a generalization check
against its own independent FE ground truth, not a re-validation of that
specific geometry (see §6.4's own text for the full caveat).

| File | Backs | Notes |
|---|---|---|
| `submit_ff_geomsweep_r150.sh` / `_r167.sh` / `_r200.sh` / `_r250.sh` | §6.4, Table geomsweep | run `python3 -m plate_solver.cli` with env vars widening `R40_FF_P1`/`R40_FF_P2` to the new radius ratio; no new Python needed since the scan windows are derived purely from the passed-in geometry's own cutoff frequencies |
| `submit_ff_geomsweep_ip_anglefill_r125.sh` | §6.4 | fills in the 3 previously-unrun Part-2 (IP) angles at the original FF-P1 ratio (r0/2b=1.25) so IP has the same 6-angle coverage as OOP |
| `figures_ff_r150/`, `_r167/`, `_r200/`, `_r250/`, `_ip_anglefill_r125/` | §6.4 | each SLURM job's `freefree_manifest.json` (per-mode accept/reject record, schema in `plate_solver/validation.py::_run_quality_review`) and `freefree_results.txt` (human-readable summary); the `freefree_checkpoint.pkl` binary each job also writes is not included — it's a solver resume checkpoint, not needed to verify the published numbers |
| `consolidate_geom_sweep_results.py` | Table geomsweep | flattens the 5 `freefree_manifest.json` files above into `geom_sweep_consolidated.csv`; run with no arguments from this directory. **This is the python-candidate side only** — `raw_omega`/`omega_lit` per accepted mode, no FE data, no match/error columns |
| `geom_sweep_consolidated.csv` | Table geomsweep | output of the above, already generated |
| `probe_geomsweep_fe_match.py` | Table geomsweep, §6.4 footnote | cross-matches the manifests against the independent FE results and reproduces §6.4's summary statistic, writing a full per-candidate audit CSV. Paths are overridable via `GEOMSWEEP_MANIFEST_DIR`/`GEOMSWEEP_FE_DIR`/`GEOMSWEEP_MATCH_CSV`/`GEOMSWEEP_CONTAM_TOL` env vars so the same script runs unmodified against either this curated copy or the canonical `FutureWork/geometry_sweep/` + `FutureWork/ansys_geomsweep_verification/` cluster working tree |
| `submit_geomsweep_fe_match.sh` | — | SLURM wrapper for the above (pure Python stdlib, no `plate_solver` import — the wrapper exists for a job-numbered provenance trail, not because the compute needs a cluster) |
| `geomsweep_fe_match_s10v2_2416746.out` / `.err` / `geomsweep_fe_match_s10v2_consolidated.csv` | §6.4, SM Table S.3 | **live headline.** Job 2416746, tol=1.0%: 1256/1475 (85.2%), mean of per-geometry means 0.574%, worst 2.990%. Per-ratio matched counts 274/294/336/352. `_s10v2_manifests.zip` is the candidate set that job staged. |
| `geomsweep_fe_match_2328808.out` / `.err` | §6.4 footnote (superseded recapture) | Job 2328808: 1231/1474 (83.5%). Kept as the footnote's intermediate figure, not the manuscript headline. |
| `geomsweep_fe_match_reconstruction.out` | — | an earlier same-script run against this repository's curated copy (rather than the cluster working tree), kept for provenance; numbers agree exactly with job 2328808 |
| `geomsweep_fe_match_consolidated.csv` | Table geomsweep | full per-candidate audit trail (2187 rows: every candidate from every geometry, flagged cutoff-adjacent / beyond-FE-range / matched, with its FE match and error% where applicable) — the artifact a referee asked to see archived. **Generated from this repository's curated copy, not pulled from job 2328808 directly** — its aggregate statistics are verified identical to the cluster job's printed output (see below), but if you need the byte-identical cluster-produced CSV, rerun `submit_geomsweep_fe_match.sh` and pull it back |
| `probe_ff_residual_screen_geomsweep.py` | §6.4's cut-off-adjacent paragraph | applies the §5 pointwise edge-residual screen (unmodified, same code as `probe_ff_oop_spurious_table.py`) to the 215 cut-off-adjacent candidates set aside from the FE-matching statistic above (jobs 2327073–2327076) |
| `probe_ff_residual_screen_reconstruction_fix.py` | §6.4's cut-off-adjacent paragraph | resolves the 44/215 candidates (12 OOP, 32 IP) the first pass couldn't reconstruct, via the same `_derive_branches_robust` escalation ladder validated elsewhere in this project (job 2327101) |
| `probe_ff_ip_relabeled_basis.py` | §6.4's cut-off-adjacent paragraph | direct follow-up: evaluates the remaining 25 stuck IP candidates at their achieved (rather than target) basis size, the same basis-relabeling convention §4's convergence tables already use (job 2327405) |
| `submit_ff_residual_screen_r150.sh` … `_r250.sh`, `submit_ff_residual_screen_reconstruction_fix.sh`, `submit_ff_ip_relabeled_basis.sh` | — | SLURM wrappers for the three probes above; adapt the SBATCH header and hardcoded `PKG_PATH`/venv activation to your own scheduler and install location |

**Status: CLOSED (2026-07-23), cluster-confirmed.** An external review
correctly noted that the script producing §6.4's original headline
statistic (1232/1466, 84.0% matched within 3%, mean 0.58%, worst 3.00%)
had been run interactively on 2026-07-22 and never committed — confirmed
absent from this working directory and an independent fresh cluster pull
checked 2026-07-23, so the figure could not be independently audited.
`probe_geomsweep_fe_match.py` is a from-scratch reconstruction built from
the two raw ingredients that *were* checked in (the python candidate
manifests and the independent FE results in `ansys/geomsweep/`), and was
run for real on the cluster as job 2328808 (`geomsweep_fe_match_2328808.out`),
not just in a sandbox.

One necessary step — removing extensional-family contamination from the
flexural FE list by frequency-matching against the companion in-plane
model — has no recorded numeric tolerance anywhere in the paper,
`LESSONS_LEARNED.md`, or any script in either copy of the working
directory; that tolerance was a human judgment call made once during the
original 2026-07-22 session and is not recoverable, only re-decidable.
The script makes an explicit, documented choice (1.0% relative Omega_lit
tolerance) and reports sensitivity across 0.5%/1.0%/2.0% rather than
presenting one number as if it were the original.

**Result (job 2328808, tolerance = 1.0%):** n_total = 1474 (then-published
1466), n matched within 3% = 1231 (83.5%, then-published 1232/84.0%), mean of
per-geometry means = 0.577% (then-published 0.58%), worst case = 3.000%
(then-published 3.00%) — count statistics agree to within 8 candidates out of
~1470 and the error statistics agree to 3 decimal places. That recapture
superseded the original interactive figure (1232/1466, 84.0%). **The live
paper (§6.4 and SM Table S.3) reports a later full recapture: 1256/1475
(85.2%), mean of per-geometry means 0.574%, worst 2.990%**, with a footnote
listing this job's 1231/1474 (83.5%, mean 0.578%, worst 3.000%) and the
original 1232/1466 as superseded recaptures of a churning candidate set.
Do not treat 1231/1474 as the manuscript's current headline.

## `s10` re-closure of the cut-off-adjacent tally (2026-08-13)

The three-script pipeline above (`probe_ff_residual_screen_geomsweep.py` +
`_reconstruction_fix.py` + `probe_ff_ip_relabeled_basis.py`) closed the
cut-off-adjacent tally under the pre-`s10` checkpoint (job 2327101 and
siblings). After the `s10` geomsweep recapture (§6.4's own footnote), all
three stages were re-run under `s10` from `FutureWork/geometry_sweep/`
(`probe_ff_residual_screen_geomsweep_s10v2.py`, jobs 2416866-2416869;
`..._s10v2_reconstruction_fix.py`, job 2417234;
`probe_ff_ip_relabeled_basis_s10v2.py`, job 2419006 — all three now live
in this directory, copied 2026-08-14 from `FutureWork/geometry_sweep/`).
Flexural closed fully at the second stage (119/119:
57 REAL-like, 62 ARTIFACT-like, 0 AMBIGUOUS, now the paper's reported
figure). Extensional needed all three stages: 71/96 resolved directly or
via the escalation ladder, then the remaining 25 (24 stuck at exactly
`n_dofs-2`, 1 shallow-lock-on) closed by the relabeled-basis step with 0
failures, matching the pre-`s10` run's own three-stage 64/7/25 split
exactly. **Combined extensional total: 88 bimodal, 8 middle-band = 96/96,
fully resolved** — the geomsweep cut-off-adjacent tally is now fully
closed for both motion types under `s10`, with the same basis-relabeling
caveat noted for the last 25 (24 at their achieved n=18, 1 at a full n=20
best-effort fill still flagged shallow after re-convergence returned to
the same root) that the pre-`s10` run itself carried. See
`paper/PAPER1_FREEFREE_SUPPLEMENTARY.tex` §S.3 (former main-text
Appendix C.3) and `LESSONS_LEARNED.md` §18.10-§18.12 for the full
derivation. The live manuscript is `paper/`, not a package-root copy.

These scripts import `plate_solver` and were run against the SLURM
cluster build. As with `validation/README.md`, some in-code comments cite
`LESSONS_LEARNED.md` section numbers and specific job numbers — that file
is the private research log and isn't part of this repository, but every
finding it's cited for is either reflected in the code itself or written
up in the paper.
