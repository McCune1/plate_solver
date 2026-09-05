# validation/phase3/

Evidence for the Phase 3 (partially-deployed) capabilities described in the
top-level `README.md`'s Status section: radially-graded (FGM) material
support and in-plane orthotropy. Unlike everything else under `validation/`,
these back a repository-status claim, not a table or figure in
`paper/PAPER1_FREEFREE_DRAFT.tex` — Phase 3 is explicitly out of scope for
that paper (see its Discussion section).

| File | Backs | Notes |
|---|---|---|
| `probe_ip_fgm_implementation_validation.py` | in-plane FGM grading | self-consistency probe: (A1) establishes an isotropic cantilever reference point, (A2) cross-checks a FreeFreeIP point against production `find_modes_sigmin`, (B) confirms a flat-profile graded material collapses onto the isotropic reference within tolerance, (C) confirms the arbitrary-precision and float64 fast paths agree on graded series coefficients to machine precision, (D) an exploratory (non-gating) isotropic-vs-graded σ_min sweep |
| `ip_fgm_implementation_validation_local_20260804_200613.out` | in-plane FGM grading | archived run output for the above — Parts A2/B/C all PASS (Part A1 establishes a reference value rather than a pass/fail; Part D was intentionally stopped early by the researcher partway through, as an exploratory sanity check with no external benchmark, not a gating test) |
| `tables_gate_2364022.out` / `.err` | regression safety | a later re-run of `../test_validated_tables.py` (the same 9-test gate whose original certifying run is `../tables_gate_2318785.out`), confirming the FGM and in-plane-orthotropy code changes introduced no regression to any already-validated table (`Ran 9 tests ... OK`) |
| `probe_ip_orthotropy_wang_table4_postfix_v2.py` | in-plane orthotropy (`mu_theta` fix) | validates the deployed `OrthotropicMaterial.ip_constants()` fix against Wang, Liang, Yao and Zhang, *Applied Mathematical Modelling* 40 (2016) 9228–9253, Table 4, FFFF row — searches for all 8 target modes through the production worker path (`find_modes_sigmin`), no monkeypatching |
| `probe_ip_orthotropy_wang_table4_postfix_v2_2334289.out` | in-plane orthotropy (`mu_theta` fix) | archived cluster run output (job 2334289, resume of job 2334099) — 8/8 targets HIT, max error 0.737% (target 3), confirms both the `ip_constants()` fix and the `nu_bar` worker-tuple wiring are present and load-bearing through the real worker path |
| `probe_fgm_graded_recursion_selftest.py` | out-of-plane FGM grading (T/R/κ) | standalone algebraic self-test: solves the graded and reduced (ungraded-equivalent) Frobenius recursions, confirms the two agree exactly at the order (a_4) where grading first enters, and plugs both series back into the original graded ODE at 5 sample points to check the residual |
| `fgm_graded_recursion_selftest_2327812.out` | out-of-plane FGM grading (T/R/κ) | archived run output (job 2327812) — a_4 cross-check PASS, worst relative ODE residual 7.1e-31 (graded) / 2.1e-31 (reduced), both at machine/arbitrary-precision noise floor |

## What's demonstrated here, and what isn't

**In-plane FGM grading:** self-consistency only (Parts A2/B/C above), plus
the regression-safety re-run. There is no external (literature or
finite-element) benchmark for graded materials anywhere in this project yet.

**In-plane orthotropy (the `mu_theta` fix)** is externally validated: Wang,
Liang, Yao and Zhang, *Applied Mathematical Modelling* 40 (2016) 9228–9253,
Table 4, FFFF row, 8/8 targets matched through the production worker path,
max error 0.737%. Archived raw log: `probe_ip_orthotropy_wang_table4_postfix_v2_2334289.out`.

**Out-of-plane FGM grading (T/R/κ)** has a self-consistency check (the
Frobenius-recursion self-test above, at arbitrary-precision noise floor) but
still has no external literature or finite-element benchmark — same status
as in-plane FGM grading. Archived raw log:
`fgm_graded_recursion_selftest_2327812.out`.
