# Paper 1 in-plane cantilever rows at 2Θ/π = 0.25: resolution (2026-09-24)

Single-process probes (run 2026-09-24 on a 2-core machine; `import plate_solver` from the repository root). They are not package code. Write-up: Paper 1 Supplementary Material §S.3.8
(Tables S.10 and S.11).

Geometry: r0/(2b) = 1.25 (R_i = 3, R_o = 7). Material: E = 210e9, ν = 0.35, ρ = 7800.
Solver: InPlaneSolver(M=80, n_quad=30), dps = 40, SOLVER_VERSION 2026-07-10.s10.

## Files

- `setup.py`: builds the solver. Set `TWOT=0.5` in the environment for the 2Θ/π = 0.5 geometry.
- `inv.py`: `inventory()` fills a complete ζ-root inventory with ladder-seeded Newton. `argcount()` checks it with an argument-principle count of det(ζ, Ω) on a rectangle. At 0.25: 14 roots in [0,30]×[0,45], all found (Ω = 0.8089 and 0.9400). At 0.5: 13 roots, all found (Ω = 0.35).
- `inv_cache*.json`: cached complete inventories. `roots.py` tracks from them.
- `gapchk.py` / `gapchk.log`: `full_search` completeness at larger `xmax`.
  - At the default `ngrid=10` and xmax = 30, it drops the 4th-lowest root (2.82+4.82i at Ω = 0.339).
  - At `ngrid=20` / `nax=400` and xmax = 36 or 44.6, it is complete.
- `conv.py` → `c_m2.jsonl`, `c_m3.jsonl`: fixed-window σ_min scan plus golden section, on the complete basis, at n = 20 to 32(36).
  - The n = 20 gate reproduces the printed values: 0.808930 vs 0.808929, and 0.940029 vs 0.940032.
- `fine.py`: a σ₁..σ₄ scan. It shows σ₂ does not dip: there is no rank-2 drop.
- `roots.py` → `r_m1.jsonl`, `r_m2.jsonl`, `r_a05m2.jsonl`: det sign scan on the complete tracked basis, plus 11-step bisection.
  - It covers the unconstrained K (u) and the Lagrange-constrained K (c).

## Result: unconstrained det roots A < B, and the constrained root C

The FE value (Q9) is 0.823991 for mode 2 and 0.348738 for mode 1.

| mode | n | A | C | B | pair width | midpoint vs FE |
|---|---|---|---|---|---|---|
| 2 | 20 | 0.808931 (−1.83%) | 0.832999 | 0.835380 | 3.22% | −0.22% |
| 2 | 24 | 0.814646 | 0.828283 | 0.829999 | 1.87% | −0.20% |
| 2 | 28 | 0.817167 | 0.826797 | 0.829142 | 1.45% | −0.10% |
| 2 | 32 | 0.819106 | 0.825618 | 0.827223 | 0.99% | −0.10% |
| 2 | 40 | 0.820851 | – | 0.825988 | 0.62% | −0.07% |
| 1 | 20 | 0.341240 (−2.15%) | 0.347160 | 0.351363 | 2.92% | −0.70% |
| 1 | 32 | 0.346927 | 0.348182 | 0.349034 | 0.61% | −0.22% |

- Mode 3 (FE 0.940601) behaves differently. At n = 20 to 28 the pair is complex: there is no sign change, only a shallow σ_min touch at about 0.9403.
- At n = 36 the pair is real, with A ≈ 0.93975 and C ≈ 0.94075, so the FE value is again bracketed.
- The printed 0.940032 is −0.06% from FE.

