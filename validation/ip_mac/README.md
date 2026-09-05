# validation/ip_mac/

Run logs backing Paper 1 §6.6 (`sec:ipmacbar`): the identity-weighted
MAC bar (`../ip_mac_bar.py`, threshold 0.744).

These are solver `.out` files from the cluster. They do not change a
computed Ω and do not live inside `find_modes_sigmin`. To re-summarize
the 24-key Rank 14 sweep after adding logs:

```
INCLUDE_A100=1 python3 probe_ip_rank14_remaining20_summary_v1.py
```

(set `PKG_PATH` to this directory, or run from a tree that also has the
a100 logs next to the summary script).

`probe_ip_mac_summary_plot_v1.py` renders the Sec 6.6 summary figure
(`paper/paper_figures/mac_separation_summary.{pdf,png}`, embedded as
Figure 3) from the logs already in this directory — a plotting pass
only, no `plate_solver` import, no cluster submission. Re-run it after
adding new calibration or sweep logs; `LOGDIR`/`OUTDIR` env vars
override the defaults (this directory / `../../paper/paper_figures`).

## Rank 18 — calibration and promotion

| File | What it backs |
|---|---|
| `probe_ip_block_selection_rule_v1_ffp1_2431872.out` | FF-P1 (1.5, 0.5) calibration, job 2431872 |
| `probe_ip_block_selection_rule_v1_r200a100_2431873.out` | r200/a100 (2.0, 1.0), job 2431873 |
| `probe_ip_block_selection_rule_v1_r200a50_2432028.out` | r200/a50 (2.0, 0.5); bit-for-bit gate table, job 2432028 |
| `probe_ip_ownroot_localmin_v1_2431917.out` | own-root floor 0.9869, job 2431917 |
| `probe_ip_mac_adversarial_fullblock_v1_2432150.out` | 57-mode full-block test (t7), job 2432150 |
| `probe_ip_mac_bar_deploy_gate_v1_2432443.out` | production post-processor deploy gate, job 2432443 |
| `probe_ip_rank8_mac_bar_v1_2432456.out` | eight r200/a100 leftover candidates, 6 ART / 2 REAL, job 2432456 |

The three-geometry table in §6.6 (floors 0.9079 / 0.9869 / 0.9994,
ceilings 0.0067 / 0.5801 / 0.3480) is the `READING` block of the three
block-selection logs.

## Rank 14 — 24-key ν=0.35 IP sweep

All 24 keys have a letter. 16 × `R14-A_TRAVELS` for `2Θ/π ≤ 1.0`;
8 × `R14-C_ALL_ART` for `2Θ/π = 1.25` and `1.50`. Combined 1596/1596
reconstructed, 912 tight, 410 REAL / 502 ART, zero misclassifications
at 0.744. Floor 0.8218 (r200/a050), ceiling 0.6803 (r250/a100).

| File | Key |
|---|---|
| `..._2433081.out` | r150/a100 |
| `..._2433083.out` | r167/a100 |
| `..._2432508.out` | r200/a100 |
| `..._2433085.out` | r250/a100 |
| `..._2433497.out` … `..._2433516.out` | remaining 20 keys (resubmit after the Cartesian-midside pairing fix) |

Each log's header has `RATIO_TAG` / `ANGLE_TAG`. Honour the printed
`READING:` letter. Do not retune 0.744.

## C-class n80 extract-completeness (CLOSED, 2026-08-26)

Jobs 2445267--2445274. All eight keys `C80-C_SURVIVES` (0 tight MAC-REAL
at 80 elastic modes, 16-way pool, ~16--21 min/key). The 57-mode qualifier
is dropped from the abstract/§6.6/conclusions. Logs:

| Job | Key | recon | tight REAL/ART | letter |
|---|---|---|---|---|
| 2445267 | r150/a125 | 82/82 | 0/54 | C80-C_SURVIVES |
| 2445268 | r150/a150 | 90/90 | 0/62 | C80-C_SURVIVES |
| 2445269 | r167/a125 | 88/88 | 0/64 | C80-C_SURVIVES |
| 2445270 | r167/a150 | 91/91 | 0/62 | C80-C_SURVIVES |
| 2445271 | r200/a125 | 100/100 | 0/69 | C80-C_SURVIVES |
| 2445272 | r200/a150 | 105/105 | 0/68 | C80-C_SURVIVES |
| 2445273 | r250/a125 | 106/106 | 0/69 | C80-C_SURVIVES |
| 2445274 | r250/a150 | 115/115 | 0/64 | C80-C_SURVIVES |

Do not retune 0.744. Eigenvector dumps stay on the cluster.

Crash logs from the first remaining-20 submit (jobs 2433457–2433495,
`ValueError: max() iterable argument is empty`) are **not** included;
those were pairing failures, superseded by 2433497–2433516.
