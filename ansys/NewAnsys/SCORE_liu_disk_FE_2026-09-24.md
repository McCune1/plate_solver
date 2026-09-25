# SCORE — Liu disk FE (lever 11) — job 2530987 — 2026-09-24

> **CORRECTION (Claude review, 2026-09-24): the odd-n "SOFT" diagnosis below is wrong.**
> The +2...+7% n=1 rows are a deck constraint, not half-model/sampling: v1 put
> *global* UX=0 on all 1017 nodes of the r=3 mm hole through the full thickness.
> An n=1 (cos theta) bending mode tilts about the y-axis at the centre, so
> u_x = -z dw/dx is nonzero there; the constraint partially clamps that tilt.
> n=0 and n>=2 have w'(0)=0 and are barely affected (hence even-n sits on Liu fem).
> v1 n=1 rows must not be cited. Fixed decks: `ansys_p4_liudisk2_sc_{C,S}_2026-09-24.inp`
> (hole free; S gets one midplane UX pin), scorer `score_liudisk2_fe_2026-09-24.py`,
> pre-registered gates in `targets_liudisk2_2026-09-24.json`. See
> `Project Knowledge/lever11/CLAUDE_REVIEW_2026-09-24.md`.
> **Confirmed by v2 (job 2531012):** with the hole free, the n=1 rows dropped 3.5–6.7% onto Liu fem;
> see `SCORE_liudisk2_FE_2026-09-24.md`.

**Status:** BOTH decks **PASS** (mesh + modal). Queue `ok=2 fail=0`.  
**Job:** `2530987` on `compute-14-37.mill` (start ~03:29 CDT / America/Chicago).  
**Mesh fix:** `RI_MESH=3 mm` annular topology; **6912** elements; RUN COMPLETED; **0** real MAPDL errors in each `_out.txt`.

Percentages recomputed from **raw FE Hz** vs `targets_liu_disk_2026-09-24.json`. Deck % columns (if any) were **not** used. Mode `(n,m)` from sample UZ Fourier (θ) + radial zero-crossings.

## Primary fundamentals

| Deck | Mode | FE Hz | Gate2 Hz | rel% vs Gate2 | Liu fem Hz | rel% vs Liu fem |
|------|------|------:|---------:|--------------:|-----------:|----------------:|
| C | (0,1) | **143.592** | 143.624 | **−0.022%** | 143.255 | +0.235% |
| S | (0,1) | **69.258** | 69.328 | **−0.102%** | 69.264 | −0.009% |

## Table 2 — clamped (C)

| (n,m) | FE mode # | FE Hz | Gate2 Hz | rel% Gate2 | Liu fem Hz | rel% Liu fem |
|-------|----------:|------:|---------:|-----------:|-----------:|-------------:|
| (0,1) | 1 | 143.592 | 143.624 | **−0.022%** | 143.255 | +0.235% |
| (1,1) | 2 | 315.687 | 298.900 | **+5.616%** | 296.458 | +6.486% |
| (2,1) | 3 | 486.637 | 490.337 | −0.755% | 485.566 | +0.221% |
| (0,2) | 4 | 554.485 | 559.143 | −0.833% | 553.095 | +0.251% |
| (1,2) | 6 | 873.908 | 855.191 | +2.189% | 839.081 | +4.151% |
| (2,2) | 8 | 1165.667 | 1189.148 | −1.975% | 1162.818 | +0.245% |

## Table 6 — simply-supported (S)

| (n,m) | FE mode # | FE Hz | Gate2 Hz | rel% Gate2 | Liu fem Hz | rel% Liu fem |
|-------|----------:|------:|---------:|-----------:|-----------:|-------------:|
| (0,1) | 1 | 69.258 | 69.328 | **−0.102%** | 69.264 | −0.009% |
| (1,1) | 2 | 208.403 | 195.348 | **+6.683%** | 193.930 | +7.463% |
| (2,1) | 3 | 356.969 | 360.053 | −0.857% | 356.889 | +0.022% |
| (0,2) | 4 | 415.014 | 417.791 | −0.665% | 414.837 | +0.043% |
| (1,2) | 7 | 701.071 | 681.524 | +2.868% | 671.793 | +4.358% |
| (2,2) | 9 | 969.317 | 985.735 | −1.666% | 968.983 | +0.034% |

## First ~12 FE frequencies (raw)

| # | C Hz | S Hz |
|--:|-----:|-----:|
| 1 | 143.592 | 69.258 |
| 2 | 315.687 | 208.403 |
| 3 | 486.637 | 356.969 |
| 4 | 554.485 | 415.014 |
| 5 | 708.827 | 554.456 |
| 6 | 873.908 | 639.937† |
| 7 | 962.888 | 701.071 |
| 8 | 1165.667 | 785.059 |
| 9 | 1227.186 | 969.317 |
| 10 | 1247.370 | 1025.405 |
| 11 | 1519.999 | 1047.184 |
| 12 | 1560.910 | 1299.182 |

† S mode 6 ≈ 639.9 Hz is **in-plane / extensional** (UZ ≪ UX); skipped for bending Table 6 matching.

## Verdict

| Check | Result |
|-------|--------|
| MAPDL | **PASS** — both `RUN COMPLETED`, `NUMBER OF ERROR MESSAGES = 0` |
| C (0,1) vs Gate2 143.624 | **PASS** (−0.022%) |
| S (0,1) vs Gate2 69.328 | **PASS** (−0.102%) |
| Even-n bending vs Gate2 | **PASS** — typically |rel| ≲ 2% |
| Even-n vs Liu fem | **PASS** — often |rel| ≲ 0.3% (FE sits on Liu fem) |
| Odd-n `(1,*)` | **SOFT** — systematically **+2%…+7%** vs Gate2 and vs Liu fem (half-model / sampling; not a mesh-fail) |
| SS midplane-ring caveat | **Not triggered** — S fundamentals are slightly **low**, not systematically high |

Honest CPT note: where Gate2 and Liu fem differ, this SOLID226 FE tends to track **Liu fem** (3-D) more tightly than Gate2 kir/CPT — expected on this thin disk (`r0/h=60`).

## Warnings (not failures)

- Shape testing: **288 / 6912** elements violate shape **warning** limits (both decks). Consistent with the tiny `RI_MESH=3 mm` hole corners; mesh still completed at 6912 elems.
- Jobname `.err` files (`p4_liu_disk_sc_*_2026_09.err`) still contain the **appended** mesh ERROR text from failed job **2530974** (same `-j`). Trust the **2530987** `_out.txt` / queue log, not the stale head of `.err`.

## Artifacts

- LIVE: `...\Plate_Solver_Package\Ansys\NewAnsys\` (job outputs + fixed pack)
- Box mirror: `/workspace/liu11/ansys_fe/job2530987/`
- Docs mirror: `/workspace/liu11/ansys/`
