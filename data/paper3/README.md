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
