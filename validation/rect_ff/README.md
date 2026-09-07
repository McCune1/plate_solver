# Rectangular free-free (FFFF) validation — Paper 2

Everything here backs the rectangular completely-free companion paper,
`paper/PAPER2_RECT_FREEFREE_DRAFT.tex` and its Supplementary Material.
Added 2026-09-07; before that date none of it was in the released tree.

## The library

`p5_rect_ff_lib.py` is the shared production library the probes import. It
holds, in one place:

- `FE_LISTS` / `IP_FE_LISTS` — the full non-rigid finite-element target
  lists at ℓ/b = 1.0, 1.5, 2.0, 2.5, 3.0, both parities. Flexural
  Λ = f / 24.6644; extensional Ω̄ = f / 804.481. These are the ground truth
  every match in the paper is measured against.
- `PERSIST_DL = 0.01` and `SIGMA_LIST = 0.3` — the frozen Screen B
  constants. Frozen before the extra-ℓ/b geometries were run; never retuned.
  `persist_one()` and `persist_one_ip()` are the screen itself.
- `n_basis()` / `n_basis_ip()` — the production and persist bases.
- `sigma_at()` / `sigma_at_ip()` — the single-point σ_min primitive: resolve
  branches, select the basis, assemble, equilibrated smallest singular value.
- `deploy_ok()` / `deploy_ok_ip()` — the in-run fidelity gate every probe
  calls before doing anything. It refuses to run unless `SOLVER_VERSION`
  matches and the deployed free-free corner term is the checkerboard.

## Guarding the corner term

Before 2026-09-03 the free-free corner was assembled as the naive same-sign
four-term sum `TT + TB + WT + WB`. That sum is **identically zero** under the
parity identity (Supplementary §S.1.3), so it assembles without raising and
silently produces the pre-corner spectrum. A tree carrying it looks healthy
and reproduces nothing.

Three independent guards now exist, cheapest first:

| guard | where | needs |
|---|---|---|
| `check_release_sync.sh` | this directory | coreutils only |
| `TestRectangularFreeFree` | `tests/test_solver.py` | numpy, scipy, mpmath |
| `probe_rect_ff_release_gate_2026-09-07.py` | this directory | the full package |

`check_release_sync.sh` md5-compares every `plate_solver/*.py` between the
deployed tree and this repository. Run it before cutting any
manuscript-matching tag — the `.tex`-only currency check in
`Project Knowledge/PAPER2_REVIEWER_PROMPT.md` cannot see code drift.

`probe_rect_ff_release_gate_2026-09-07.py` is the full gate: `SOLVER_VERSION`,
the corner formula, the in-plane `bc='free_free'` switch, two stored
single-point σ_min anchors at dps 26/30/40, and (optionally) the
release-vs-deployed md5 comparison. Its verdicts are pre-registered in its
own docstring. `submit_rect_ff_release_gate_2026-09-07.sh` runs it on the
cluster.

The two stored anchors, on the checkerboard build at `SOLVER_VERSION`
`2026-07-10.s10`, identical at dps 26, 30 and 40:

```
OOP  l/b=1.5 SYM  Lambda* = 0.964   basis (6, 0)  im_cap 30.0  ->  2.3835918992545717e-04
IP   l/b=2.0 ANTI Omega*  = 0.5200  basis (5, 3)  im_cap  7.0  ->  3.900364925307852e-08
```

The pre-2026-09-03 build gives `3.950702663724406e-03` for the flexural
anchor — a factor of 16.6, not a rounding difference.

**If an anchor moves, do not edit it.** A moved anchor means a code change
altered a computed rectangular FFFF frequency, which is exactly what
`SOLVER_VERSION` exists to record: bump it and re-derive Tables 1–2 and the
92-row in-plane list before touching the paper.

## The production probes

Grouped by the manuscript section they back. Cluster job numbers are in
Supplementary §S.4.

**Corner term (§3.2, §S.1)**
- `probe_rect_ff_oop_derived_vs_fix_2026-09-03.py` — the algebraic gate
  (naive sum ≡ 0; checkerboard = 2 × single-point) and the four-anchor FE
  cross-check.

**Flexural discovery and refinement (§4.1, §4.2)**
- `probe_rect_ff_oop_extra_lob_2026-09-04.py`,
  `probe_rect_ff_oop_refine_pair_2026-09-04.py`,
  `probe_rect_ff_oop_refine_extra_lob_2026-09-04.py`,
  `probe_rect_ff_oop_highlam_2026-09-04.py`,
  `probe_rect_ff_oop_highlam_extra_lob_2026-09-04.py`.

**The named unmatched square-plate mode (§4.4)**
- `probe_rect_ff_oop_1982_reopen_2026-09-05.py`,
  `probe_rect_ff_oop_1982_wide_gap_2026-09-05.py`.

**Eigenvector cross-check (§4.3, §S.5)**
- `probe_rect_ff_oop_mac_2026-09-05.py`,
  `probe_rect_ff_oop_mac_strengthen_2026-09-05.py`,
  `probe_rect_ff_oop_mac_anti_remain_2026-09-06.py`.
  The MAC bar 0.744 is the annular companion's production call, imported
  unre-tuned; none of these probes retunes it.

**In-plane discovery (§5.2)**
- `probe_rect_ff_ip_disc_2026-09-05.py`,
  `probe_rect_ff_ip_disc_gap_enlarged_2026-09-05.py`,
  `probe_rect_ff_ip_disc_v2_ncp3_2026-09-05.py`,
  `probe_rect_ff_ip_anti_persist_enlarged_2026-09-05.py`,
  `probe_rect_ff_ip_match_refine_2026-09-05.py`,
  `probe_rect_ff_ip_uncovered16_ncp5_2026-09-05.py`.

**Screen B continuum-margin disclosure (§3.3)**
- `probe_rect_ff_persist_threshold_sensitivity_2026-09-07.py` (all 148
  published matches against the frozen fixed-step scan),
  `probe_rect_ff_persist_finegrid_recheck_2026-09-07.py` (4× finer grid on
  the three that sat at the scan's outer edge),
  `probe_rect_ff_persist_golden_polish_2026-09-07.py` (golden-section
  continuum minimisation, which settled them).

Each probe's own docstring carries its pre-registered interpretation
criteria; read those before reading its output.

## Finite-element decks

The APDL decks, queue manifests and submit scripts are in
`ansys/NewAnsys/`, together with the per-aspect-ratio frequency listings the
`FE_LISTS` above were parsed from. The per-mode nodal (X, Y, U_Z) dumps
behind the §4.3 MAC comparison are ~88 MB and are not tracked; regenerate
them from `ansys_rect_ff_oop_mac_extract_lob*.inp`.

ANSYS has one licence seat here: submit decks through
`submit_ansys_queue_rect_ff_*.sh`, never as bare concurrent `sbatch` jobs.
