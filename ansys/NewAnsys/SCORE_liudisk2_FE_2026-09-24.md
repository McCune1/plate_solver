# SCORE — Liu disk FE v2 (lever 11) — job 2531012 — 2026-09-24

**Verdict: PASS (F0–F3, all pre-registered in `targets_liudisk2_2026-09-24.json`).**
Queue `ok=2 fail=0` on compute-34-04 with ANSYS 2025R1. Both transcripts show
RUN COMPLETED and 0 errors. The only warnings (2 each) are the known
288/6912 shape warnings near the 3 mm hole.

Claude verified this independently on Makise:
- Rescoring the returned files locally gives output byte-identical to the
  cluster `score_liudisk2_fe_2026-09-24.txt`.
- The transcripts show NIN = 1017 hole nodes with NO constraint applied.
- The C deck has exactly one UX block, the outer clamp.
- The S deck has exactly one UX block, the pin, with NPIN = 1.
- No mode is below 1 Hz. All 12 modes are identified with purity 1.00.

| bc | n | m | FE Hz | Liu fem Hz | vs fem | CPT (consistent) | vs CPT | v1 (2530987) | vs v1 |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|
| C | 0 | 1 | 143.583 | 143.255 | +0.229% | 143.624 | −0.028% | 143.592 | −0.006% |
| C | 1 | 1 | 297.815 | 296.458 | +0.458% | 298.899 | −0.363% | 315.687 | −5.661% |
| C | 2 | 1 | 486.612 | 485.566 | +0.216% | 490.335 | −0.759% | 486.637 | −0.005% |
| C | 0 | 2 | 554.406 | 553.095 | +0.237% | 559.142 | −0.847% | 554.485 | −0.014% |
| C | 1 | 2 | 843.378 | 839.081 | +0.512% | 855.189 | −1.381% | 873.908 | −3.494% |
| C | 2 | 2 | 1165.584 | 1162.818 | +0.238% | 1189.145 | −1.981% | 1165.667 | −0.007% |
| S | 0 | 1 | 69.254 | 69.264 | −0.015% | 69.328 | −0.107% | 69.258 | −0.006% |
| S | 1 | 1 | 194.474 | 193.930 | +0.280% | 195.347 | −0.447% | 208.403 | −6.684% |
| S | 2 | 1 | 356.954 | 356.889 | +0.018% | 360.052 | −0.861% | 356.969 | −0.004% |
| S | 0 | 2 | 414.962 | 414.837 | +0.030% | 417.790 | −0.677% | 415.014 | −0.012% |
| S | 1 | 2 | 673.839 | 671.793 | +0.305% | 681.522 | −1.127% | 701.071 | −3.884% |
| S | 2 | 2 | 969.253 | 968.983 | +0.028% | 985.733 | −1.672% | 969.317 | −0.007% |

## Reading
- **F1 (the fix):** the n=1 rows dropped 3.5–6.7% and now sit on Liu's
  ABAQUS column. This confirms the v1 inner UX=0 diagnosis.
- **F2:** the even-n rows are unchanged to ≤0.014%, so freeing the hole
  cost nothing.
- All 12 FE values are now below CPT, as 3-D vs Kirchhoff must be (v1's
  n=1 rows were above CPT, which was itself a red flag). The gap grows
  with mode order, reaching 2%.
- **Residual, reported but not chased:** the n=1 rows sit about 0.25 pp
  higher against Liu fem than the even-n rows (C +0.46/+0.51 vs +0.22…+0.24;
  S +0.28/+0.31 vs ≤0.03). That is inside the pre-registered 0.6% bar.
  Possible causes are the single mesh, or Liu's ABAQUS clamp and mesh. It
  was not investigated, because it doesn't change any claim.
- This is one mesh, with no convergence study. The claim in SM §S.10 is
  limited to "reproduces Liu's ABAQUS column within 0.6%".

Promoted to `github_repo/ansys/NewAnsys/`: v1 and v2 decks, outputs,
builders, manifests, submits, targets, scorer, score txt and queue logs.
The v1 set is kept for provenance and is marked superseded.
