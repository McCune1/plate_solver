# Paper 2 — rectangular FFFF in-plane match master table

Built 2026-09-05 from cluster job outputs. Canonical paper
location: `github_repo/paper/`. Not a submitted table — the
main-text Bardell comparison is the 12-row extract in
`PAPER2_RECT_FREEFREE_DRAFT.md` §4. Job numbers belong in the
Supplementary provenance list at submission.

## Conversion (Bardell 1996 → this project)

This project (Seok/Tiersten/Scarton Part 2, Eqs. 22–23):

- Plate: \(-\ell \le x \le \ell\), \(-b \le y \le b\); aspect \(\ell/b\).
- \(\bar\Omega = \omega/\bar\omega\), \(\bar\omega = (\pi/(2b))\sqrt{c_{66}/\rho}\),
  \(c_{66}=E/(2(1+\nu))\).
- FE conversion: \(\bar\Omega = f[\mathrm{Hz}]/804.481\) at
  \(b=1\,\mathrm{m}\), \(E=210\,\mathrm{GPa}\), \(\nu=0.30\),
  \(\rho=7800\,\mathrm{kg\,m^{-3}}\) (independent of \(\ell/b\)).

Bardell, Langley & Dunsdon, *J. Sound Vib.* **191**(3) (1996) 459–467:

- Plate: \(0\le x\le a\), \(0\le y\le b\); aspect \(a/b\).
- \(\Omega_B = \omega a/C_L\), \(C_L^2 = E/[\rho(1-\nu^2)]\).
- Frequencies are printed under Figure 1 mode plots, not tabulated.
- Poisson ratio is not stated numerically. Table 1 S–S–S–S exact
  \((0,1)=(1,0)=1.859\) on the square equals
  \(\pi\sqrt{(1-\nu)/2}=1.8586\) at \(\nu=0.30\), so the
  conversion uses \(\nu=0.30\).

Geometry identification: \(a=2\ell\), \(b_{\mathrm{Bardell}}=2b\), so
\(a/b=\ell/b\). Then

$$\Omega_B = \bar\Omega\cdot\pi\cdot(\ell/b)\cdot\sqrt{(1-\nu)/2}.$$

**Figure 1 (a)/(b) labels are swapped relative to the caption and
to the body text.** Caption: “(a) \(a/b=1\); (b) \(a/b=2\)”.
Panel (a) is drawn 2:1 with unique frequencies
\(\{1.954, 2.961, 3.267, 4.726, 4.784, 5.205\}\); panel (b) is
square with the repeated pair 2.472, 2.472. Body text (p. 462)
assigns repeats to \(a/b=1\) and unique frequencies to \(a/b=2\).
Table 1 (S–S–S–S) does **not** have this swap. Mapping used here
follows the rendered proportions and the repeated-frequency
physics, same caution as Paper 1 §7 / `LESSONS_LEARNED.md` §18.72a.

## Bardell F–F–F–F first six, converted

| \(a/b=\ell/b\) | class | Bardell \(\Omega_B\) | \(\Omega_B\to\bar\Omega\) | \(\bar\Omega^*\) | \(\bar\Omega_{\mathrm{FE}}\) | miss vs Bardell | miss FE vs Bardell | job |
| ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| 1.0 | ANTI | 2.321 | 1.24880 | 1.2400 | 1.24858 | 0.70% | 0.017% | 2456002 |
| 1.0 | SYM | 2.472 | 1.33004 | 1.3300 | 1.32984 | 0.00% | 0.015% | 2456002 |
| 1.0 | ANTI | 2.472 | 1.33004 | 1.3298 | 1.32984 | 0.02% | 0.015% | 2456013 |
| 1.0 | SYM | 2.628 | 1.41397 | 1.4100 | 1.41422 | 0.28% | 0.017% | 2456002 |
| 1.0 | SYM | 2.987 | 1.60713 | 1.6100 | 1.60733 | 0.18% | 0.012% | 2456002 |
| 1.0 | SYM | 3.452 | 1.85732 | 1.8600 | 1.85745 | 0.14% | 0.007% | 2456002 |
| 2.0 | ANTI | 1.954 | 0.52567 | 0.5200 | 0.52557 | 1.08% | 0.018% | 2456002 |
| 2.0 | SYM | 2.961 | 0.79657 | 0.8100 | 0.79652 | 1.69% | 0.006% | 2456002 |
| 2.0 | ANTI | 3.267 | 0.87889 | 0.8800 | 0.87891 | 0.13% | 0.002% | 2456002 |
| 2.0 | ANTI | 4.726 | 1.27139 | 1.2715 | 1.27148 | 0.01% | 0.007% | 2456013 |
| 2.0 | ANTI | 4.784 | 1.28700 | 1.2800 | 1.28703 | 0.54% | 0.003% | 2456002 |
| 2.0 | SYM | 5.205 | 1.40025 | 1.4000 | 1.40011 | 0.02% | 0.010% | 2456002 |

12/12 of Bardell's published F–F–F–F first-six non-zero frequencies at both aspect ratios are matched. Solver vs converted Bardell: max miss 1.69%, median 0.18%. Independent FE vs converted Bardell: max miss 0.018%.

The three largest solver misses (1.69%, 1.08%, 0.70%) are
`n_cpair=3` unique-MATCH roots from job 2456002 that already
passed the 3% FE cut on the refine grid (step 0.002, printed
\(\bar\Omega^*\) to 0.01) and were therefore not sent to the
`n_cpair=5` mop-up. The two mop-up roots that fall inside
Bardell's first six sit at 0.01–0.02%.

## Full 92/92 unique MATCH against `IP_FE_LISTS`

Window \(\bar\Omega\in[0.02,2.50]\), all five aspect ratios,
both parities. 76 from job 2456002 (`n_cpair=3` discovery, 3%
cut) + 16 from job 2456013 (`n_cpair=5` targeted mop-up).
FE from job 2455569 (PLANE183 half-model, mesh2, max mesh-conv
0.004%). `SOLVER_VERSION=2026-07-10.s10`.

| \(\ell/b\) | class | \(\bar\Omega^*\) | \(\bar\Omega_{\mathrm{FE}}\) | miss | job |
| ---: | --- | ---: | ---: | ---: | --- |
| 1.0 | SYM | 1.3300 | 1.32984 | 0.01% | 2456002 |
| 1.0 | SYM | 1.4100 | 1.41422 | 0.30% | 2456002 |
| 1.0 | SYM | 1.6100 | 1.60733 | 0.17% | 2456002 |
| 1.0 | SYM | 1.8600 | 1.85745 | 0.14% | 2456002 |
| 1.0 | SYM | 2.0000 | 2.00319 | 0.16% | 2456002 |
| 1.0 | ANTI | 1.2400 | 1.24858 | 0.69% | 2456002 |
| 1.0 | ANTI | 1.3298 | 1.32984 | 0.00% | 2456013 |
| 1.0 | ANTI | 2.0000 | 2.00319 | 0.16% | 2456002 |
| 1.0 | ANTI | 2.3200 | 2.31523 | 0.21% | 2456002 |
| 1.5 | SYM | 1.0500 | 1.04551 | 0.43% | 2456002 |
| 1.5 | SYM | 1.4100 | 1.41236 | 0.17% | 2456002 |
| 1.5 | SYM | 1.4243 | 1.42432 | 0.00% | 2456013 |
| 1.5 | SYM | 1.6300 | 1.62439 | 0.35% | 2456002 |
| 1.5 | SYM | 1.6900 | 1.69877 | 0.52% | 2456002 |
| 1.5 | SYM | 2.1100 | 2.10886 | 0.05% | 2456002 |
| 1.5 | SYM | 2.2200 | 2.22185 | 0.08% | 2456002 |
| 1.5 | SYM | 2.2800 | 2.28542 | 0.24% | 2456002 |
| 1.5 | ANTI | 0.7800 | 0.78782 | 0.99% | 2456002 |
| 1.5 | ANTI | 1.0333 | 1.03335 | 0.00% | 2456013 |
| 1.5 | ANTI | 1.5712 | 1.57124 | 0.00% | 2456013 |
| 1.5 | ANTI | 1.6400 | 1.64765 | 0.46% | 2456002 |
| 1.5 | ANTI | 2.2025 | 2.20253 | 0.00% | 2456013 |
| 1.5 | ANTI | 2.2400 | 2.24065 | 0.03% | 2456002 |
| 2.0 | SYM | 0.8100 | 0.79652 | 1.69% | 2456002 |
| 2.0 | SYM | 1.4000 | 1.40011 | 0.01% | 2456002 |
| 2.0 | SYM | 1.4142 | 1.41422 | 0.00% | 2456013 |
| 2.0 | SYM | 1.4400 | 1.44333 | 0.23% | 2456002 |
| 2.0 | SYM | 1.6500 | 1.65355 | 0.21% | 2456002 |
| 2.0 | SYM | 1.7500 | 1.77390 | 1.35% | 2456002 |
| 2.0 | SYM | 1.8100 | 1.81566 | 0.31% | 2456002 |
| 2.0 | SYM | 2.0000 | 2.00440 | 0.22% | 2456002 |
| 2.0 | SYM | 2.3000 | 2.29856 | 0.06% | 2456002 |
| 2.0 | SYM | 2.3300 | 2.33104 | 0.04% | 2456002 |
| 2.0 | ANTI | 0.5200 | 0.52557 | 1.06% | 2456002 |
| 2.0 | ANTI | 0.8800 | 0.87891 | 0.12% | 2456002 |
| 2.0 | ANTI | 1.2715 | 1.27148 | 0.00% | 2456013 |
| 2.0 | ANTI | 1.2800 | 1.28703 | 0.55% | 2456002 |
| 2.0 | ANTI | 1.7400 | 1.73452 | 0.32% | 2456002 |
| 2.0 | ANTI | 1.8400 | 1.84467 | 0.25% | 2456002 |
| 2.0 | ANTI | 2.1400 | 2.13652 | 0.16% | 2456002 |
| 2.5 | SYM | 0.6400 | 0.64049 | 0.08% | 2456002 |
| 2.5 | SYM | 1.2200 | 1.22293 | 0.24% | 2456002 |
| 2.5 | SYM | 1.4100 | 1.41225 | 0.16% | 2456002 |
| 2.5 | SYM | 1.4266 | 1.42658 | 0.00% | 2456013 |
| 2.5 | SYM | 1.5400 | 1.53449 | 0.36% | 2456002 |
| 2.5 | SYM | 1.6700 | 1.67003 | 0.00% | 2456002 |
| 2.5 | SYM | 1.6900 | 1.69350 | 0.21% | 2456002 |
| 2.5 | SYM | 1.7800 | 1.78140 | 0.08% | 2456002 |
| 2.5 | SYM | 2.0500 | 2.04232 | 0.38% | 2456002 |
| 2.5 | SYM | 2.0700 | 2.07671 | 0.32% | 2456002 |
| 2.5 | SYM | 2.1300 | 2.13048 | 0.02% | 2456002 |
| 2.5 | SYM | 2.3900 | 2.38996 | 0.00% | 2456002 |
| 2.5 | SYM | 2.4047 | 2.40468 | 0.00% | 2456013 |
| 2.5 | ANTI | 0.3800 | 0.37585 | 1.10% | 2456002 |
| 2.5 | ANTI | 0.7200 | 0.71776 | 0.31% | 2456002 |
| 2.5 | ANTI | 1.0600 | 1.06883 | 0.83% | 2456002 |
| 2.5 | ANTI | 1.1125 | 1.11248 | 0.00% | 2456013 |
| 2.5 | ANTI | 1.4400 | 1.44246 | 0.17% | 2456002 |
| 2.5 | ANTI | 1.5400 | 1.53132 | 0.57% | 2456002 |
| 2.5 | ANTI | 1.8000 | 1.80308 | 0.17% | 2456002 |
| 2.5 | ANTI | 2.0000 | 1.99255 | 0.37% | 2456002 |
| 2.5 | ANTI | 2.2000 | 2.14745 | 2.45% | 2456002 |
| 2.5 | ANTI | 2.2800 | 2.34878 | 2.93% | 2456002 |
| 3.0 | SYM | 0.5400 | 0.53502 | 0.93% | 2456002 |
| 3.0 | SYM | 1.0600 | 1.04554 | 1.38% | 2456002 |
| 3.0 | SYM | 1.4142 | 1.41422 | 0.00% | 2456013 |
| 3.0 | SYM | 1.4199 | 1.41791 | 0.14% | 2456013 |
| 3.0 | SYM | 1.4200 | 1.41978 | 0.02% | 2456002 |
| 3.0 | SYM | 1.6000 | 1.60066 | 0.04% | 2456002 |
| 3.0 | SYM | 1.6635 | 1.66346 | 0.00% | 2456013 |
| 3.0 | SYM | 1.6900 | 1.67627 | 0.82% | 2456002 |
| 3.0 | SYM | 1.8000 | 1.80401 | 0.22% | 2456002 |
| 3.0 | SYM | 1.8400 | 1.83870 | 0.07% | 2456002 |
| 3.0 | SYM | 2.0000 | 2.00478 | 0.24% | 2456002 |
| 3.0 | SYM | 2.1600 | 2.16287 | 0.13% | 2456002 |
| 3.0 | SYM | 2.2100 | 2.21145 | 0.07% | 2456002 |
| 3.0 | SYM | 2.2700 | 2.27554 | 0.24% | 2456002 |
| 3.0 | SYM | 2.4300 | 2.44758 | 0.72% | 2456002 |
| 3.0 | SYM | 2.4600 | 2.46034 | 0.01% | 2456002 |
| 3.0 | ANTI | 0.2800 | 0.28169 | 0.60% | 2456002 |
| 3.0 | ANTI | 0.5800 | 0.57789 | 0.37% | 2456002 |
| 3.0 | ANTI | 0.8800 | 0.88739 | 0.83% | 2456002 |
| 3.0 | ANTI | 1.0315 | 1.03153 | 0.00% | 2456013 |
| 3.0 | ANTI | 1.2702 | 1.27024 | 0.00% | 2456013 |
| 3.0 | ANTI | 1.3000 | 1.29465 | 0.41% | 2456002 |
| 3.0 | ANTI | 1.5800 | 1.57876 | 0.08% | 2456002 |
| 3.0 | ANTI | 1.6600 | 1.66509 | 0.31% | 2456002 |
| 3.0 | ANTI | 1.9000 | 1.91234 | 0.65% | 2456002 |
| 3.0 | ANTI | 2.0000 | 1.99856 | 0.07% | 2456002 |
| 3.0 | ANTI | 2.2400 | 2.25049 | 0.47% | 2456002 |
| 3.0 | ANTI | 2.2800 | 2.27620 | 0.17% | 2456002 |
| 3.0 | ANTI | 2.4937 | 2.49371 | 0.00% | 2456013 |

Row count: 92. `IP_FE_LISTS` entry count: 92.
Zero UNCOVERED.

