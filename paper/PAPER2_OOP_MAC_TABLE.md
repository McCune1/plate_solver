# Paper 2 rectangular FFFF OOP MAC table

Off-cluster computation (`probe_rect_ff_oop_mac_2026-09-05.py`).
FE eigenvectors from job 2456738 (five SHELL281 MAC-extract
decks). Analytical shapes from `RectOOPAssembler(bc='free_free')`
at the published $\Lambda^*$ with the same production / persist
basis that produced that $\Lambda^*$. `SOLVER_VERSION` unmodified.

MAC $= |\langle u_z^{\mathrm{FE}}, w^{\mathrm{an}}\rangle|^2 / (\|u_z^{\mathrm{FE}}\|^2 \|w^{\mathrm{an}}\|^2)$ on the
phase-aligned real part of $w$. Bar 0.744 is Paper 1's, not retuned.
FE $\Lambda$ is $f_{\mathrm{Hz}}/24.6644$ (the deck-printed
`Lambda_FE` column is the KFAC/*VWRITE bug and is ignored).

## Tables 1–2 matches

| tab | $\ell/b$ | class | $\Lambda^*$ | $\Lambda_{\mathrm{FE}}$ | FE# | MAC(FE) | best FE# | MAC(best) | 2nd MAC | im/re | letter | note |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 1.0 | ANTI | 1.366 | 1.35288 | 4 | 0.9408 | 4 | 0.9408 | 0.0451 | 0.000 | CONFIRMED |  |
| 1 | 1.5 | SYM | 0.964 | 0.96347 | 4 | 1.0000 | 4 | 1.0000 | 0.0017 | 0.000 | CONFIRMED | 1st SYM |
| 1 | 1.5 | SYM | 2.254 | 2.24373 | 5 | 0.9995 | 5 | 0.9995 | 0.0006 | 0.000 | CONFIRMED |  |
| 1 | 1.5 | ANTI | 0.906 | 0.89838 | 4 | 0.8521 | 4 | 0.8521 | 0.1434 | 0.000 | CONFIRMED | 1st ANTI |
| 1 | 1.5 | ANTI | 2.124 | 2.07036 | 5 | 0.2721 | 5 | 0.2721 | 0.1266 | 0.000 | SHAPE_MISMATCH | weakest primary |
| 1 | 2.0 | SYM | 0.5434 | 0.54338 | 4 | 0.9997 | 4 | 0.9997 | 0.0009 | 0.000 | CONFIRMED | persist basis only |
| 1 | 2.0 | SYM | 1.51 | 1.50735 | 5 | 0.8493 | 5 | 0.8493 | 0.0099 | 0.000 | CONFIRMED_PERSIST |  persist-basis shape |
| 1 | 2.0 | SYM | 2.234 | 2.22558 | 6 | 0.9994 | 6 | 0.9994 | 0.0006 | 0.000 | CONFIRMED |  |
| 1 | 2.0 | ANTI | 0.674 | 0.66858 | 4 | 0.7507 | 4 | 0.7507 | 0.2204 | 0.000 | CONFIRMED | 1st ANTI |
| 1 | 2.0 | ANTI | 1.482 | 1.47069 | 5 | 0.0000 | 4 | 0.6940 | 0.2588 | 0.000 | SHAPE_MISMATCH |  |
| 1 | 2.5 | SYM | 0.348 | 0.34766 | 4 | 0.9991 | 4 | 0.9991 | 0.0004 | 0.000 | CONFIRMED | 1st SYM |
| 1 | 2.5 | SYM | 0.978 | 0.96559 | 5 | 0.8709 | 5 | 0.8709 | 0.0028 | 0.000 | CONFIRMED_PERSIST |  persist-basis shape |
| 1 | 2.5 | SYM | 1.888 | 1.8837 | 6 | 0.9985 | 6 | 0.9985 | 0.0018 | 0.000 | CONFIRMED |  |
| 1 | 2.5 | SYM | 2.278 | 2.27059 | 7 | 0.9991 | 7 | 0.9991 | 0.0008 | 0.000 | CONFIRMED |  |
| 1 | 2.5 | ANTI | 0.534 | 0.53126 | 4 | 0.6548 | 4 | 0.6548 | 0.2739 | 0.000 | AMBIGUOUS | 1st ANTI |
| 1 | 2.5 | ANTI | 1.148 | 1.14003 | 5 | 0.0000 | 4 | 0.6105 | 0.2969 | 0.000 | SHAPE_MISMATCH |  |
| 1 | 2.5 | ANTI | 1.918 | 1.90392 | 6 | 1.0000 | 6 | 1.0000 | 0.0013 | 0.000 | CONFIRMED |  |
| 1 | 3.0 | SYM | 0.242 | 0.24127 | 4 | 0.9591 | 4 | 0.9591 | 0.0012 | 0.000 | CONFIRMED | 1st SYM |
| 1 | 3.0 | SYM | 0.6698 | 0.66984 | 5 | 0.0000 | 8 | 0.4689 | 0.1395 | 0.000 | SHAPE_MISMATCH | persist basis only |
| 1 | 3.0 | SYM | 1.32 | 1.31785 | 6 | 0.9996 | 6 | 0.9996 | 0.0008 | 0.000 | CONFIRMED |  |
| 1 | 3.0 | SYM | 2.162 | 2.15417 | 7 | 0.7788 | 7 | 0.7788 | 0.1187 | 0.000 | CONFIRMED_PERSIST |  persist-basis shape |
| 1 | 3.0 | SYM | 2.256 | 2.24841 | 8 | 0.9994 | 8 | 0.9994 | 0.0002 | 0.000 | CONFIRMED |  |
| 1 | 3.0 | ANTI | 0.444 | 0.44044 | 4 | 1.0000 | 4 | 1.0000 | 0.0014 | 0.000 | CONFIRMED | 1st ANTI; 0.46-family coincidence |
| 1 | 3.0 | ANTI | 0.936 | 0.93044 | 5 | 0.0000 | 4 | 0.5410 | 0.3161 | 0.000 | SHAPE_MISMATCH |  |
| 1 | 3.0 | ANTI | 1.528 | 1.51698 | 6 | 0.3279 | 4 | 0.4980 | 0.3279 | 0.000 | AMBIGUOUS |  |
| 1 | 3.0 | ANTI | 2.258 | 2.24305 | 7 | 0.0000 | 12 | 0.4718 | 0.1496 | 0.000 | SHAPE_MISMATCH |  |
| 2 | 1.0 | SYM | 2.46 | 2.45439 | 5 | 0.9974 | 5 | 0.9974 | 0.0026 | 0.000 | CONFIRMED |  |
| 2 | 1.0 | SYM | 3.53 | 3.4925 | 6 | 0.0010 | 7 | 0.0884 | 0.0142 | 0.000 | SHAPE_MISMATCH |  |
| 2 | 1.0 | SYM | 6.19 | 6.16026 | 7 | 0.1564 | 7 | 0.1564 | 0.0282 | 0.000 | SHAPE_MISMATCH |  |
| 2 | 1.0 | SYM | 6.46 | 6.36553 | 8 | 0.9999 | 8 | 0.9999 | 0.0010 | 0.000 | CONFIRMED |  |
| 2 | 1.0 | ANTI | 3.52 | 3.4925 | 5 | 0.0000 | 7 | 0.4432 | 0.3334 | 0.000 | SHAPE_MISMATCH |  |
| 2 | 1.5 | SYM | 5.38 | 5.38152 | 9 | 0.0000 | 6 | 0.1399 | 0.0962 | 0.000 | SHAPE_MISMATCH |  |
| 2 | 1.5 | ANTI | 3.86 | 3.83354 | 6 | 0.9934 | 6 | 0.9934 | 0.0051 | 0.000 | CONFIRMED |  |
| 2 | 2.0 | SYM | 3.0 | 2.99812 | 8 | 0.8617 | 8 | 0.8617 | 0.0612 | 0.000 | CONFIRMED |  |
| 2 | 2.0 | SYM | 3.65 | 3.62227 | 9 | 0.9998 | 9 | 0.9998 | 0.0004 | 0.000 | CONFIRMED |  |
| 2 | 2.0 | ANTI | 2.58 | 2.5518 | 6 | 0.9992 | 6 | 0.9992 | 0.0018 | 0.000 | CONFIRMED |  |
| 2 | 2.0 | ANTI | 4.06 | 4.02613 | 7 | 0.0000 | 10 | 0.5988 | 0.0905 | 0.000 | SHAPE_MISMATCH |  |
| 2 | 2.0 | ANTI | 6.2 | 6.17866 | 9 | 0.0000 | 8 | 0.4696 | 0.2100 | 0.000 | SHAPE_MISMATCH |  |
| 2 | 2.5 | SYM | 3.19 | 3.19522 | 10 | 0.0000 | 9 | 0.9999 | 0.0001 | 0.000 | PAIRING_SWAP |  |
| 2 | 2.5 | SYM | 4.18 | 4.13819 | 11 | 0.0018 | 5 | 0.3009 | 0.1055 | 0.000 | SHAPE_MISMATCH |  |
| 2 | 2.5 | SYM | 5.41 | 5.34819 | 13 | 0.0000 | 5 | 0.3132 | 0.1252 | 0.000 | SHAPE_MISMATCH |  |
| 2 | 2.5 | ANTI | 2.92 | 2.89206 | 7 | 0.0000 | 11 | 0.7485 | 0.1425 | 0.000 | PAIRING_SWAP |  |
| 2 | 2.5 | ANTI | 4.18 | 4.14962 | 8 | 0.0009 | 11 | 0.6450 | 0.1321 | 0.000 | SHAPE_MISMATCH |  |
| 2 | 2.5 | ANTI | 5.72 | 5.67275 | 9 | 0.0244 | 5 | 0.0587 | 0.0362 | 0.000 | SHAPE_MISMATCH |  |
| 2 | 2.5 | ANTI | 6.26 | 6.21965 | 10 | 0.0234 | 12 | 0.0418 | 0.0238 | 0.000 | SHAPE_MISMATCH |  |
| 2 | 2.5 | ANTI | 6.4 | 6.36131 | 11 | 0.9999 | 11 | 0.9999 | 0.0012 | 0.000 | CONFIRMED |  |
| 2 | 3.0 | SYM | 2.48 | 2.45935 | 9 | 0.0000 | 8 | 0.5501 | 0.1667 | 0.000 | SHAPE_MISMATCH |  |
| 2 | 3.0 | SYM | 3.34 | 3.32528 | 11 | 0.9928 | 11 | 0.9928 | 0.0021 | 0.000 | CONFIRMED |  |
| 2 | 3.0 | SYM | 3.66 | 3.62467 | 12 | 0.0005 | 5 | 0.3116 | 0.0836 | 0.000 | SHAPE_MISMATCH |  |
| 2 | 3.0 | SYM | 5.66 | 5.59746 | 15 | 0.0004 | 5 | 0.3173 | 0.1115 | 0.000 | SHAPE_MISMATCH |  |
| 2 | 3.0 | SYM | 6.18 | 6.14744 | 16 | 0.0000 | 5 | 0.3160 | 0.1163 | 0.000 | SHAPE_MISMATCH |  |
| 2 | 3.0 | ANTI | 3.16 | 3.14253 | 8 | 0.9985 | 8 | 0.9985 | 0.0020 | 0.000 | CONFIRMED |  |
| 2 | 3.0 | ANTI | 4.26 | 4.23529 | 9 | 0.0000 | 12 | 0.6186 | 0.1445 | 0.000 | SHAPE_MISMATCH |  |
| 2 | 3.0 | ANTI | 5.56 | 5.51872 | 10 | 0.9627 | 10 | 0.9627 | 0.0228 | 0.000 | CONFIRMED |  |
| 2 | 3.0 | ANTI | 6.22 | 6.18191 | 11 | 0.0082 | 5 | 0.2237 | 0.1747 | 0.000 | SHAPE_MISMATCH |  |
| 2 | 3.0 | ANTI | 6.38 | 6.34339 | 12 | 0.9999 | 12 | 0.9999 | 0.0015 | 0.000 | CONFIRMED |  |

## Negative-control leaks (not in Tables 1–2)

| $\ell/b$ | class | $\Lambda^*$ | best FE# | $\Lambda_{\mathrm{best}}$ | MAC(best) | 2nd MAC | letter | note |
|---|---|---|---|---|---|---|---|---|
| 1.0 | SYM | 0.41 | 7 | 6.160243914305639 | 0.0059 | 0.0023 | LEAK_LOW | SYM 0.41 leak |
| 1.5 | SYM | 0.41 | 6 | 2.588085256483028 | 0.0098 | 0.0031 | LEAK_LOW | SYM 0.41 leak |
| 2.5 | SYM | 0.41 | 5 | 0.9655865538995475 | 0.0503 | 0.0049 | LEAK_LOW | SYM 0.41 leak |
| 1.0 | ANTI | 0.46 | 5 | 3.4925013785050516 | 0.1853 | 0.0385 | LEAK_LOW | ANTI 0.46 leak |
| 1.5 | ANTI | 0.46 | 5 | 2.0703583707692057 | 0.3327 | 0.0711 | LEAK_AMB | ANTI 0.46 leak |
| 2.5 | ANTI | 0.46 | 5 | 1.1400310974521983 | 0.4341 | 0.1264 | LEAK_AMB | ANTI 0.46 leak |

## Summary

- Tables 1–2 pairs scored: 56
- CONFIRMED (MAC ≥ 0.744 vs frequency partner, production basis): 26
- CONFIRMED_PERSIST (same bar, Screen B n_cpair=3 basis): 3
- PAIRING_SWAP remaining: 2
- AMBIGUOUS remaining: 2
- SHAPE_MISMATCH remaining: 23
- D_GATE: 0
- MAC(FE) min / median / max: 0.0000 / 0.7648 / 1.0000

### Table 1 (primary, $\Lambda\lesssim 2.4$) is the clean subset

18/26 Table 1 rows meet the 0.744 bar against their frequency partner
(15 production + 3 persist-basis recoveries). Every published 1st SYM
is MAC $\ge 0.959$. 1st ANTI: 0.941 / 0.852 / 0.751 / 0.655 / 1.000
at $\ell/b=1.0\to 3.0$ (only $\ell/b=2.5$ sits in the 0.3–0.744 band).
Known Screen-B leaks stay low (SYM 0.41: MAC $\le 0.05$ vs every FE
mode; ANTI 0.46: 0.19–0.43, never $\ge 0.744$).

### One pairing swap MAC actually resolves

$\ell/b=2.5$ SYM $\Lambda^*=3.190$ was frequency-matched to
$\Lambda_{\mathrm{FE}}=3.19522$ (FE#10). MAC $=0.9999$ against FE#9
($\Lambda=3.16376$) and $0$ against FE#10. Those two FE values are
1% apart and both sit in `FE_LISTS`; nearest-frequency matching picked
the slightly closer one, the eigenvector picks the other. This is MAC
doing the job it was introduced for, not a failed reconstruction.

### What MAC does *not* yet confirm

23 Table 1–2 rows remain SHAPE_MISMATCH after a persist-basis retry.
They cluster in (i) ANTI “second” modes whose reconstructed $w$ looks
like the 1st ANTI, (ii) higher-$\Lambda$ Table 2 entries, (iii) the
persist-only $\ell/b=3.0$ SYM $0.6698$ row (frequency miss 0.00%, but
$\sigma\sim 10^{-12}$ and MAC $=0$ vs FE#5). Frequency matches in
Tables 1–2 are independent of this and are not withdrawn.

Paper text was not edited. Fold-in proposal is in the session report,
not in this file.
