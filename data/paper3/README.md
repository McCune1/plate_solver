# data/paper3/

Tabulated solver-versus-literature JSON for Paper 3 (closed ring and
solid disk), plus the transcription module `p3_ring_disk_lib.py` that
holds the printed SP-160 / Narita / Irie / Southwell cells.

These files are the numerical source of the Paper 3 tables. They were
captured under `SOLVER_VERSION` `2026-07-10.s10`. Cluster probe scripts
and mixed-edge APDL decks remain in the working package; they are not
in the Paper 1/2 tags.

| File | Contents |
|---|---|
| `p3_ring_disk_lib.py` | printed literature + `is_vogel_230_named_skip` |
| `p3_phaseA_vogel235.json` | Table 2.35, both ν, Narita 2(a) |
| `p3_phaseB_disk.json` | Narita 2(b), Table 2.5, rank, tiny-hole |
| `p3_phaseD_fc_cf.json` | mixed-edge API table |
| `p3_phaseD_thinring_diag.json` | Tables 2.22 / 2.30 / 2.31 production |
| `p3_phaseD_inverse_ba_close.json` | Southwell inverse b/a (Fig. S.1) |
| `p3_phaseE_irie.json` | Irie Table 2 |
| `p3_phaseF_vogel_ss_narita.json` | mill job 2472653: remaining Vogel 2.18/2.20/2.24/2.26/2.28/2.33, Table 2.16, Narita 3 |
| `p3_mill2473295_irie_ip_mixed.json` | mill job 2473295: Irie 1984 Tables 3-5, in-plane mixed-edge ring (F-C/C-F/C-C), 144/144 |
| `p3_mill2473295_irie_tinyhole.json` | mill job 2473295: Irie's circular column recovered as a tiny-hole ring (beta=0.01), 24/24 |
| `p3_mill2473295_disk_cs.json` | mill job 2473295: solid disk OOP outer clamped and simply supported, SP-160 Tables 2.1, 2.3, 12/12 |
| `p3_mill2473295_guided_ring.json` | mill job 2473295: guided-guided closed ring vs. Bhaskara Rao and Kameswara Rao (J. Solid Mech. 2016), 29/29 |

Note: this batch is named by job number (`mill2473295`), not the next
`phase<letter>` in the A/B/D/E/F sequence above, because `Phase G` is
already taken in this project's convention for the LaTeX-writing
milestone (`Project Knowledge/PAPER3_PHASEG_HANDOFF_PROMPT_2026-09-09.md`),
a different lettering scheme from the data-batch phases here. Four
`p3_phaseG_*.json` deprecation stubs from that naming false start were
deleted 2026-09-11.
