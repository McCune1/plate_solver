# validation/

Reproduction scripts and their evidence logs, each backing a specific
validated table or claim in `paper/PAPER1_FREEFREE_DRAFT.tex` (canonical
manuscript directory; Supplementary Material is
`paper/PAPER1_FREEFREE_SUPPLEMENTARY.tex`).

| File | Backs | Notes |
|---|---|---|
| `test_validated_tables.py` | §6's pinned frequency tables | unittest suite; also runs the three heavy "discovery" gates (fresh `det K(Ω)=0` search from scratch, not just a check against pinned values) |
| `submit_validated_tables_gate.sh` | — | SLURM wrapper for the above; adapt the SBATCH header to your scheduler |
| `tables_gate_2318785.out` / `.err` | §6 | archived clean-pass output ("Ran 9 tests ... OK") from the run that certified the tables |
| `check_deployed_tree.py` | — | preflight check that a deployed copy of `plate_solver/` matches the expected module set before trusting a batch run |
| `e2_parity_check_step1.py`, `e2_parity_check_step2.py` | — | standalone SymPy algebra check for an open item ("E.2": the annular solver's clamped-edge orientation is hardcoded rather than generalized like the free edge — currently harmless, but must be understood before any new boundary-condition layout, including rectangular free-free); not package-dependent |
| `probe_mac_ambiguity_check.py` / `submit_probe_mac_ambiguity_check.sh` | §6.3 | MAC-correlation check resolving a near-degenerate ODD-parity mode pair against the ANSYS eigenvector export in `ansys/` |
| `probe_ip_orthotropy_wang_table4_Er20_v1.py` | Table `tab:wang` | generating script for job 2408719; copied 2026-08-14 |
| `probe_ip_qin2018_ffff_table2_v7.py` | Table `tab:qin` | generating script for job 2409260; copied 2026-08-14 |
| `probe_ip_qin2018_ffff_pair_prod_n20_2421062.out`, `_n28_2421063.out` | §6.9 pair-split check | n=20/28 at xmax=60: still one production-detector root |
| `probe_ip_qin2018_second_dip_v1.py` / `_2421905.out` | SM §S.3.4 isolated dip | Ω*=1.14758, Ω_Qin=4.4717, depths −5.18 (n=20) / −6.07 (n=28) |
| `probe_ip_qin_pair_mac_v2_2434753.out` | SM §S.3.4 MAC | 0.9996 @ FE#7 ANTI and 0.9871 @ FE#8 SYMM; isolate Hz 459.66/461.60 vs Table 7 FE 459.5/460.6 |
| `probe_oop_shi2014_refine_hits_v1.py` | Table `tab:shi2014oop` | generating script for job 2414020; copied 2026-08-14 |
| `probe_oop_shi2014_t4_unconstrained_v1_2421061.out` | Table `tab:shi2014oop` row 53.649 | unconstrained scan: two isolated stable roots; published value assigned to 53.6553 |
| `probe_oop_shi2014_t3_larger_basis_v1.py` / `submit_oop_shi2014_t3_larger_basis_v1.sh` | Table `tab:shi2014oop` row 38.428 | generating script for job 2421581; copied 2026-08-15 |
| `probe_oop_shi2014_t3_larger_basis_v1_2421581.out` | Table `tab:shi2014oop` row 38.428 | n=20 near-miss clears at n=28/36; refined root 38.4295 / 0.004% / −4.706 |
| `probe_oop_shi2014_table3_t3_local_v1.py` | Table `tab:shi2014oop` n=20 near-miss | generating script for job 2411852; copied 2026-08-14 |
| `probe_ff_adjudicate_candidates.py` | Tables S.4 / S.5 / S.6 | geometry-agnostic residual + FE rematch; copied 2026-08-14. Rank 18 closeout 2026-08-19: exports `R40_BC_KIND=free_free` before the pool; IP `classify_ip()` heuristic is labeled as not the two-sided call. |
| `ip_mac_bar.py` | §6.6 `sec:ipmacbar` | production two-sided IP post-processor (Rank 18 CLOSED 2026-08-19; Rank 14 applied 2026-08-20 at all 24 IP geomsweep keys): identity-weighted MAC, `RULE_LOCALMIN_RESID`, class-safe nearest, threshold 0.744. Cartesian-midside pairing fallback for non-180° sectors. Not inside `find_modes_sigmin`; no `SOLVER_VERSION` bump. |
| `ip_mac/` | §6.6 | Rank 18 calibration + Rank 8 leftover + Rank 14 24-key `.out` logs + C-class 80-mode extract (jobs 2445267--2445274). See `ip_mac/README.md`. |
| `ip_mac_r200a100_fe_parity_mesh96.txt` | §6.6 | SYMM/ANTI classification of the 57-mode r200/a100 mesh96 dump (29/28/0); FE-vs-FE max off-diagonal MAC 0.0115 between modes 4 and 32. |
| `probe_ff_r200_a50_capture_nu030.py` / `submit_ff_r200_a50_capture_nu030_oop.sh` / `_ip.sh` | Table S.6, main-text crossed-geometry subsection | production capture at \(R_i/R_o=0.6\), \(2\Theta=\pi/2\), \(\nu=0.30\) (jobs 2421898/2421899) |
| `merge_r200_a50_nu030_candidates.py` / `submit_ff_adjudicate_r200_a50_nu030.sh` | Table S.6 | merge + adjudicate (jobs 2422897/2422898) |
| `ff_adjudicate_r200_a50_nu030_2422898.json` / `.out` | Table S.6 | 10/10 OOP REAL-like FE-matched; 43/47 IP FE-matched |
| `probe_ff_oop_spurious_table.py` / `submit_ff_oop_spurious_table.sh` | `tab:oop-spurious-sm` (Table S.8) | regenerates the raw Ω / residual values for the 8 physical + 5 in-window spurious out-of-plane points |
| `ff_oop_spurious_table_2325529.out` | `tab:oop-spurious-sm` | archived run output for the above (13/13 confirmed, 0 mismatches). Table S.8's sixth row (Ω=2.817613) is `ff_oop_artifact_2p8176_2414022.out` below |
| `geomsweep/` | §6.4, Table S.3, §S.3.3 | reproduction scripts + evidence for the 48-geometry generalization sweep and its cut-off-adjacent residual screen — see `geomsweep/README.md` |
| `phase3/` | top-level README's Phase 3 status (not a paper claim — Phase 3 is out of scope for the paper) | evidence for in-plane FGM grading self-consistency and the regression-safety re-run after Phase 3 deployment; see `phase3/README.md` for what is and isn't independently reproducible from this repository yet |

Two entries in the table above map to Supplementary Material subsections
that were introduced when the paper's tables were reorganized:
`probe_mac_ambiguity_check.py` backs instrument 6 of the eight-instrument
list in §S.3.2 (`sec:instruments`), and `geomsweep/` backs §S.3.3
(`sec:cutoff-appendix`) as well as §6.4 itself.

## Evidence logs added 2026-08-13, completed 2026-08-13 (second pass)

All seven originally-requested runs are now in this directory, renamed to
drop download-artifact ` (n)` suffixes; duplicates confirmed byte-identical
by `md5sum` and then removed.

| File | Backs | Notes |
|---|---|---|
| `probe_ip_orthotropy_wang_table4_Er20_v1_2408719.out` | Table 8 (`tab:wang`), targets 5–8 | |
| `probe_ip_qin2018_ffff_table2_v7_2409260.out` | Table 9 (`tab:qin`) | |
| `probe_oop_shi2014_table3_r04_v3_2410389.out` | superseded evidence — see note below, **not** the run behind the current Table 10 |
| `ff_omstar_audit_2324817.out` | Table S.1 (`tab:conv-oop`) | |
| `ff_fillfix_omstar_2324822.out` | Table S.2 (`tab:conv-ip`) | |
| `probe_oop_shi2014_refine_hits_v1_2414020.out` | Table 10 (`tab:shi2014oop`), five refined rows | spot-verified line-by-line against the printed table 2026-08-13; exact match on those five refined/window-confirmed rows |
| `probe_oop_shi2014_table3_t3_local_v1_2411852.out` | Table 10 (`tab:shi2014oop`), target-3 n=20 near-miss | spot-verified 2026-08-13: deepest point log₁₀σ_min=-3.297 at Ω_lit=38.428; superseded as the published status by job 2421581 |
| `probe_oop_shi2014_t3_larger_basis_v1_2421581.out` | Table 10 (`tab:shi2014oop`), target-3 refined root | spot-verified 2026-08-15: n=28 refine 38.4295 / 0.004% / −4.706 matches the printed row |
| `probe_ip_ffp1_ip06_n36_branch_threshold_v1_2411851.out` | §S.3.1's IP-06 basis-count footnote | spot-verified 2026-08-13: constant cnt=30 across [1.196,1.199], matching the footnote's negative result for that window |
| `probe_ip_ffp1_ip06_n36_upper_window_v1_2414031.out` | §S.3.1's IP-06 basis-count footnote | spot-verified 2026-08-13: 30↔34 transition band at Ω≈1.2032–1.2040 matches the footnote exactly |
| `ff_adjudicate_r150_a100_nu030_2381555.json` / `.out` | Table S.4 (`tab:r150a100`), ν=0.30 | copied 2026-08-14 from the r150 archive folder |
| `ff_adjudicate_r150_a100_nu035_2347329.json` / `.out` | Table S.4 (`tab:r150a100`), ν=0.35 | copied 2026-08-14 from the r150 archive folder |
| `ff_adjudicate_r200_a100_nu030_2410396.json` / `.out` | Table S.5 (`tab:r200a100`) | present, not independently re-derived against the printed table this pass |
| `probe_ip_r200_a100_flagged8_v1_2421361.out` | §6.5 eight high-rTyy FE matches | residual-direction pass; printed C (REAL control failed depth bar). Superseded as the two-sided call: `ip_mac/probe_ip_rank8_mac_bar_v1_2432456.out` (see `ip_mac/README.md`) scores these eight 6 ARTIFACT-like / 2 REAL-like on the MAC bar; §6.5 and §6.6 report that call. |
| `probe_oop_r150_a100_nu035_om143_v1_2421415.out` | Table S.4 ν=0.35 OOP hole | Ω=1.432743 filled ARTIFACT-like, FE-matched |
| `ff_oop_artifact_2p8176_2414022.out` | Table S.8 last row / Fig. 2(e) | r_V=180.6, r_M=75.52 at raw Ω=2.817613 |

**Table 10 (`tab:shi2014oop`) note.** `..._2410389.out` is the *original*
coarse local-window confirmation, kept because it's still cited in the
paper's own methodology narrative (§6.10 describes superseding it with the
refined result) — but it does not by itself back the table as printed.
The table as printed is backed by `probe_oop_shi2014_refine_hits_v1_2414020.out`
(five of six targets refined), `probe_oop_shi2014_t4_unconstrained_v1_2421061.out`
(row 53.649, two isolated roots), and
`probe_oop_shi2014_t3_larger_basis_v1_2421581.out` (the sixth target,
upgraded from the n=20 near-miss in `..._2411852.out` to a refined root
at n=28/36). All four logs are now present above.

**§S.3.1 note.** The IP-06 basis-count footnote (Tables S.1–S.2's own
`$^\ddagger$` note, "Ω-dependent, not irreproducible") is backed by
`probe_ip_ffp1_ip06_n36_branch_threshold_v1_2411851.out` and
`probe_ip_ffp1_ip06_n36_upper_window_v1_2414031.out`, both now present
above and spot-verified against the footnote text.

## `geomsweep/`'s cut-off-adjacent tally, fully re-closed under `s10` (2026-08-13)

`geomsweep/` (see its own README) backs §6.4 and §S.3.3
(`sec:cutoff-appendix`). As of 2026-08-13, the cut-off-adjacent tally was
fully re-derived under the current `s10` checkpoint in three stages: a
direct-reconstruction pass (jobs 2416866-2416869,
`geomsweep/probe_ff_residual_screen_geomsweep_s10v2.py`), the pre-`s10`
escalation ladder ported to `s10` (job 2417234,
`geomsweep/probe_ff_residual_screen_geomsweep_s10v2_reconstruction_fix.py`),
and the pre-`s10` relabeled-basis best-effort technique also ported to
`s10` (job 2419006, `geomsweep/probe_ff_ip_relabeled_basis_s10v2.py`).
Flexural
closed fully at the escalation-ladder stage: 119/119 candidates resolved,
57 REAL-like / 62 ARTIFACT-like / 0 AMBIGUOUS, now the paper's reported
figure (supersedes the pre-`s10` 53/60/6). Extensional needed all three
stages: 71/96 resolved directly or via the escalation ladder, then the
remaining 25 (24 stuck at exactly `n_dofs-2`, 1 shallow-lock-on) closed by
the relabeled-basis step with 0 failures — 24 evaluated at their achieved
smaller basis, 1 at a full-size best-effort fill still flagged shallow
after an independent re-convergence attempt returned to the same root.
Combined **extensional total: 88 bimodal, 8 middle-band = 96/96, fully
resolved** — the same three-stage structure (64 direct / 7 escalated / 25
relabeled) as the pre-`s10` run's own closure. See §S.3.3 and
`LESSONS_LEARNED.md` §18.11–§18.12 for the full derivation.

`validation/phase3/probe_ip_orthotropy_wang_table4_postfix_v2*` is **not**
the run behind Table 8 (`tab:wang`): it used `E_r = 40 E_θ`, the material from a
different table in the same source paper, and is superseded. It is kept
only as the Phase-3 in-plane-orthotropy self-consistency evidence the
top-level README's Status section refers to.

These scripts import `plate_solver` and were run against the SLURM cluster
build unless noted otherwise (`e2_parity_check_*` are pure SymPy; `sympy`
is listed in `requirements.txt` as an optional dependency for exactly this
reason).

Two things to know before rerunning any `submit_*.sh` here: they hardcode
`PKG_PATH=/home/ghmkfh/PythonMill/Plate_Solver_Package` and a matching venv
activation path from the original cluster deployment — change both to your
own install location. And some in-code comments cite `LESSONS_LEARNED.md`
section numbers and specific job numbers; that file is the private research
log and isn't part of this repository, but every finding it's cited for is
either reflected in the code itself or written up in the paper.
