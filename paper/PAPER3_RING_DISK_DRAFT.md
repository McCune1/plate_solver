# Paper 3 draft — Closed ring and solid disk

**Status (2026-09-09).** Markdown draft complete enough to typeset.
LaTeX is the next session. **Contract for that session:**
`Project Knowledge/PAPER3_TEX_WRITING_BRIEF.md`. Paste-able prompt:
`PAPER3_PHASEG_HANDOFF_PROMPT_2026-09-09.md`.

This file is the working manuscript until
`PAPER3_RING_DISK_DRAFT.tex` exists. Live Papers 1–2 stay frozen in
this folder. Do not fold this into either companion. Do not bump
`SOLVER_VERSION` (`2026-07-10.s10`). Job numbers belong in SM
provenance, not main text.

Roadmap: `Project Knowledge/PAPER3_RING_DISK_ROADMAP.md`. Numbers:
`LESSONS_LEARNED.md` §18.112–§18.124. Bans and the number ledger in
the writing brief override this file if they ever disagree.

Working title (not submission-final):

> Exact radial-determinant frequencies of annular and circular plates:
> the closed-ring and solid-disk limits of the sector wave-function method

Venue: JSV, third paper in the Seok–Tiersten exact-solver series.

---

# Exact Radial-Determinant Frequencies of Annular and Circular Plates

**Garrek McCune and Daniel Stutts**
Mechanical and Aerospace Engineering, Missouri University of Science
and Technology, Rolla, MO 65409, USA

*Companion to Paper 1 (annular-sector FFFF) and Paper 2 (rectangular
FFFF). Draft; Phases A–E of the ring/disk programme are closed.*

## Abstract

The Seok–Tiersten–Scarton wave-function method for annular-sector plates
has two classical limits that do not use the two-edge Galerkin
determinant at all. At sector angle \(2\Theta=2\pi\) there are no radial
walls: circumferential order is an integer \(n\), and each frequency is
a zero of a \(4\times4\) radial edge-determinant already present in the
sector assembler. At inner radius \(R_i=0\) the two origin-singular
branches drop out and the same operator becomes a \(2\times2\) at the
outer radius — the classical Kirchhoff disk. Because there are no
\(\theta\)-edges, the Kirchhoff corner jump that blocked simply-supported
sectors never appears, and the nine inner/outer combinations of Vogel
and Skinner are different row pairs of one matrix.

Frequencies are found from \(\det L(n,\Omega)=0\) from scratch. Published
tables are used only in post-hoc reports. A frozen sign-flip ladder
(inner five of five real-part flips, with positive and negative controls
in the same run) replaces the sector residual bar, which does not
travel to this geometry. Independent SHELL281 calculations confirm the
flexural ring and disk, and mixed inner/outer edges, with a one-sided
Kirchhoff-versus-Mindlin offset.

Out-of-plane free–free annuli recover SP-160 Table 2.35 at 80/80
(MATCH+WEAK) across five \(b/a\) and both \(\nu=0.30\) and \(\nu=1/3\),
and Narita’s four five-figure hold-outs at \(b/a=0.3\), \(\nu=0.3\).
The solid-disk \(2\times2\) recovers Narita Table 2(b) 4/4 and Leissa
Table 2.5 6/6; a naive \(4\times4\) at \(R_i=0\) is rank 3 for \(n\ge1\)
and is not the disk. Mixed-edge Tables 2.22 and 2.30 pass at 80/80 and
72/76 (the \(n=0\) column at \(b/a=0.9\) in 2.30 is a named source
defect, internally inconsistent with Vogel’s own \(n=1\) column).
Southwell Table 2.31 scores 13/18 under the same 1%/3% bar; the five
direct misses are a sensitivity to printed \(b/a\), not a solver
defect, and independent FE at \(b/a=0.840\), \(n=0\) gives
\(\lambda^2=132.549\) against the solver 132.990 (+0.333%) and
Southwell’s printed 81.0 (−63.6%). In-plane free–free rings recover
Irie 1984 Table 2 at 48/48. Default clamped-free sector assembly is
unchanged.

**Keywords:** annular plate; circular plate; free vibration; exact
solution; radial determinant; mixed boundary; in-plane vibration.

---

## 1. Introduction

Annular and circular plates are the \(2\pi\) and \(R_i=0\) limits of the
annular-sector geometry treated in Paper 1. Classical Kirchhoff tables
for those limits already exist (Leissa SP-160, compiling Vogel and
Skinner, Southwell, Joga-Rao and Pickett; Narita; Irie). The
contribution here is not another Bessel catalogue. It is the reduction
of a *validated sector solver* onto those tables, both motions, with a
screen that does not use the sector residual bar, plus the mixed-edge
combinations that the sector corner problem blocked.

Seok, Tiersten and Scarton gave exact-wave-function solutions for
annular-sector plates, out-of-plane and in-plane, at clamped-free.
Paper 1 extended the flexural assembler to completely free sectors and
built a residual screen against independent FE. Paper 2 did the
rectangular analogue, where the Kirchhoff corner jump is load-bearing.
Neither paper takes the \(2\pi\) or \(R_i=0\) limit as a production
entry point. That limit is cheaper than the sector: integer \(n\), a
\(4\times4\) (ring) or \(2\times2\) (disk), no \(\theta\)-walls, no
corner jump. It is also the elastic forward model that later
piezoelectric identification will need on the F–F and C–F ring.

The detector remains paper-faithful. Frequencies are zeros of
\(\det L(n,\Omega)=0\) from the same Frobenius radial content as Paper
1. Published values never seed a search. The sign-flip ladder of
§18.38 is frozen before the production tables and is not retuned to
include a published match. Independent cyclic SHELL281 (flexure) and
PLANE183 (in-plane) calculations are a third method, not a fitted
target.

What this paper is not: McGee’s solid *sector* (that is Paper 1’s
assembler with a coordinate singularity); elliptical plates; SS/GUIDED
on sectors (the open corner-jump problem); piezoelectric coupling
(one future-work sentence); a splice into Papers 1 or 2.

---

## 2. Reductions

### 2.1 Closed ring, out-of-plane

\(2\pi\)-periodicity forces integer circumferential order \(n\in
\mathbb{N}_0\). For \(n\ge1\) the \(\cos n\theta/\sin n\theta\) pair is
degenerate (same \(\Omega\)). Rigid-body modes on F–F are \(n=0,s=0\)
translation and \(n=1,s=0\) rotation; they are excluded from elastic
tables, as in Papers 1–2.

The radial operator is the existing \(4\times4\) `_Lmat_mp`. Rows are
\(0=M_r(\mathrm{inner})\), \(1=M_r(\mathrm{outer})\),
\(2=V_r(\mathrm{inner})\), \(3=V_r(\mathrm{outer})\). Searching
\(\det L(n,\Omega)=0\) at fixed integer \(n\) *is* the ring. The four
Frobenius branches span \(\{J_n,Y_n,I_n,K_n\}\) under the project’s
\(r(x)\) map (SM: Bessel span). Boundary-condition row swaps:

- Free: leave \(M_r\) and \(V_r\).
- Clamped: replace both of that edge’s rows with \(W\) and \(W'\).
- Simply supported (not exercised in the production tables): replace
  only the \(V_r\) row with \(W\); leave \(M_r\).

F–C is inner F / outer C (Leissa title “Clamped, Free”). C–F is inner C
/ outer F (Leissa “Free, Clamped”). C–C is admitted by the same API and
was validated against the clamped–clamped strip limit
\(\lambda^2\to 22.3733\,(a/w)^2\) (ratios 0.9999/1.0000/1.0000 at
\(b/a=0.90/0.94/0.97\)); it is Phase F, not a production table here.

### 2.2 Solid disk, out-of-plane

Boundedness at \(r=0\) kills \(Y_n\) and \(K_n\). Two remaining
coefficients and two outer-edge conditions give a \(2\times2\). Regular
solutions behave as \(r^{|n|}\). Inner jet conditions that isolate the
2-D kernel:

- \(n=0\): \(W'=W'''=0\)
- \(n=1\): \(W=W''=0\)
- \(n\ge2\): \(W=W'=0\)

A documented failure: evaluating the *\(4\times4\)* at \(R_i=0\) gives
rank 2 at \(n=0\) and rank 3 (not 2) at \(n\ge1\). Frequencies from
that matrix are not Leissa’s disk. The production path is the
\(2\times2\). A tiny-hole \(4\times4\) at \(b/a=10^{-3}\) is a check,
not the method. That check misses a pre-registered \(10^{-3}\) relative
tolerance on \((2,0)\) and \((0,1)\) (0.94% and 0.55% remainder of a
hole versus a solid); the probe still prints FAIL on that control and
the bar is not loosened.

### 2.3 In-plane ring

`InPlaneSolver._Lmat_mp` is \(N_r\) and \(N_{r\theta}\) at both arcs.
Same role swap: fix \(\xi=n\), search \(\Omega\). No Kirchhoff corner.
Irie’s \(\lambda=\omega a\sqrt{\rho(1-\nu^2)/E}\) is a different symbol
from flexural \(\lambda^2\) and must not reuse the flexural factor.
Irie’s circular column is \(\beta=0.01\), a tiny-hole control, not
\(R_i=0\) unless the 2×2 agrees.

### 2.4 Conversion

\[
\lambda^2 = \omega a^2\sqrt{\rho/D},\qquad a=R_o,\quad
D=\frac{Eh^3}{12(1-\nu^2)}.
\]

The factor is geometry-specific. On the FF-P1 / full-ring steel plate
(\(R_o=8\,\mathrm{m}\), \(H=0.08\,\mathrm{m}\), \(b/a=0.5\)) it is
\(4\pi^2\approx 39.478\ldots\); \(b/a=0.3\) must not reuse it. Unit
test: Paper 1 FF-P1 mode 7 native \(\Omega=1.369611\) converts to Ansys
\(\Omega_\mathrm{lit}=54.0755\) at 0.01%. Never print native \(\Omega\)
labelled as \(\lambda^2\).

Compare \(\nu=1/3\) runs to SP-160 Tables 2.34/2.35; compare
\(\nu=0.30\) to Narita and to FE. Do not mix.

---

## 3. Screen

The sector depth/gap bar failed on five genuine full-ring roots. This
geometry uses, frozen before the production tables:

1. Coarse scan of \(\log_{10}|\det L|\) on a window local to the
   literature target, plus a gap-free coarse scan per \(n\).
2. Golden-section polish of each candidate.
3. Real-part sign-flip ladder, \(\delta=10^{-2}\ldots10^{-6}\). CLEAN =
   inner 5/5. FLOOR = the only inner miss is \(\delta=10^{-6}\), outer
   5/5, shrinking (the dps=40 floor). NOT = any inner miss at
   \(\delta\ge10^{-5}\).
4. Positive controls CLEAN; negative controls NOT, in the same run.

A literature MATCH is relative miss versus the printed number, 1% /
WEAK 3%. FLOOR does not convert a MATCH into a MISS. The ladder is
never retuned to include a published match. This is not Screen B and
not SUBDOM.

A search window is a hypothesis. If \(\log|\det|\) falls monotonically
to an endpoint, the searcher can return that boundary. The value is
meaningless and \(n\)-independent. The ladder catches it (NOT). A
label `edge_flag` (EXACT / NEAR / NONE) is applied after the fact and
does not move the returned \(\Omega\). \(\lambda^2=92.566\) at
Southwell \(n=0\), \(b/a=0.840\) is such a boundary and is never a
frequency.

---

## 4. Out-of-plane free–free annulus

SP-160 Table 2.35, five \(b/a\in\{0.1,0.3,0.5,0.7,0.9\}\), eight
\((n,s)\) rows. **80/80 MATCH+WEAK.** At \(\nu=0.30\) every cell is
MATCH (max miss about 0.4%). At \(\nu=1/3\) several \((2,0)\)/\((3,0)\)
cells are WEAK at 1–2%; those \((2,0)\) numbers sit on Table 2.34
(Raju, \(\nu=1/3\)), not on Vogel 2.35. Paper-facing comparison is
\(\nu=0.30\) against 2.35; 2.35’s \(\nu\) is not the \(1/3\) of 2.34.

Narita 1984 Table 2(a) at \(b/a=0.3\), \(\nu=0.3\) is the five-figure
hold-out: 4.9060, 8.3535, 12.266, 18.292. Printed five figures held.
(3,0) sign-flip was inner 4/5 at \(\delta=10^{-6}\) (precision floor,
not a missing frequency).

The five §18.38 roots at \(b/a=0.5\), \(\nu=0.30\) are a fidelity gate
on the production API, not a new discovery. FE at that geometry
reproduces them bit-for-bit as the elastic listing 4.2645, 9.3114,
11.4082, 17.0973, 21.0296. Mesh convergence (ncirc 64/96/128): mode 7
= 4.264506 / 4.264461 / 4.264453. Contours confirm published
\((n,s)\), including (3,0)/(1,1)/(4,0) that were previously inferred
from frequency alone.

At \(b/a=0.9\), FE modes 15–16 are in-plane (\(|U_Z|\sim10^{-13}\));
do not n-ID them as flexural.

---

## 5. Out-of-plane free disk

Narita Table 2(b), \(\nu=0.33\), four five-figure hold-outs, **4/4
MATCH** (rel \(3\times10^{-6}\) to \(9\times10^{-6}\)). Leissa Table
2.5 six rows **6/6 MATCH** (printed 5.253 versus solver 5.262 is 0.17%,
inside the 1% three-figure bar; Narita is the hold-out).

FE \(\nu=0.33\): 5.2584, 9.0679, 12.2284, 20.5010 versus Narita
(0.07/0.01/0.13/0.06%) and versus the \(2\times2\)
(0.069/0.011/0.13/0.059%). Published \((n,s)\) confirmed on the
contour.

Rank diagnosis (SM): \(4\times4\) at \(R_i=0\) remains rank 3 for
\(n\ge1\). Tiny-hole \(b/a=10^{-3}\) versus \(2\times2\) is a named
FAIL on a pre-registered \(10^{-3}\) bar for (2,0) and (0,1); hole
roots are real and deep, so this is remainder of a hole versus a
solid, not a missed disk frequency.

---

## 6. Mixed inner/outer (F–C and C–F)

This is the stator boundary condition and the content that is not a
FFFF echo of Papers 1–2. No piezo claim.

**Table 2.22** (inner F, outer C): **80/80 MATCH+WEAK**, both \(\nu\),
all five \(b/a\) including 0.9.

**Table 2.30** (inner C, outer F): paper score **72/76**. Mill
production print was 72/78, skipping only the printed 51.5 cell. The
paper treats the entire \(n=0\) column at \(b/a=0.9\) as one named
source skip. Justification is source-internal: at that geometry
\(n=0\) and \(n=1\) are near-degenerate, and Vogel’s own Table 2.22
prints it (360/362, 2219/2220). Table 2.30 prints 51.5 versus 345 and
970 versus 2189. The \(n=1\) entries match the \(4\times4\) to 0.06%;
the \(n=0\) entries match nothing. The 1%/3% bar is not loosened to
absorb the column. Combined 2.22+2.30 is 152/156.

**Table 2.16** axisymmetric exact, \(b/a=0.5\), \(\nu=1/3\), \(n=0\):
F–C printed 17.51, C–F printed 13.05, both MATCH under 1%.

**Table 2.31 (Southwell, \(\nu=0.3\))** is scored **13/18** in the
results table, under the same rule as every other table. Do not invent
a second metric for this table.

The MATCH/MISS split is a sensitivity effect, not an accuracy effect.
Southwell assumed Bessel arguments (round \(\lambda\)) and inverted for
\(b/a\), printed to 2–3 decimals. Holding printed \(\lambda^2\) exact
and root-finding \(b/a\) recovers a CLEAN \(s=0\) root on all 18 rows.
Every row, including MATCH rows, carries a \(b/a\) discrepancy of order
\(10^{-3}\) (control band 0.00016–0.00686). Predicted relative error
\(\mathrm{d}\lambda^2/\mathrm{d}(b/a)\times|\delta(b/a)|/\lambda^2\)
reproduces every row’s forward miss to a factor of 2, 18/18. The five
direct misses are the steep-slope rows.

Load-bearing evidence is FE, not a differential test against Vogel.

Decisive cell, C–F, \(b/a=0.840\), \(n=0\):

| quantity | \(\lambda^2\) | vs FE |
|---|---|---|
| FE | **132.5490** | — |
| solver \(4\times4\) | 132.9899 | **+0.333%** |
| Southwell printed | 81.0 | **−63.6%** |

Nothing in the 40-mode FE spectrum lies below 132.549. Mesh refinement
reproduces 132.548973 versus 132.548990 (relative \(1.3\times10^{-7}\)).
That row is a source error.

Eighteen solver-versus-FE points across three C–F geometries, all
CLEAN, all \(s=0\), all within +0.04% to +0.76%. The solver sits
uniformly *above* FE; the gap grows monotonically in \(n\). That is
Kirchhoff versus Mindlin (SHELL281 carries transverse shear),
one-signed and physical. Mixed-sign scatter would have been the
worrying result.

Deck positive control, F–C, \(b/a=0.5\), \(n=0\): Vogel 17.7, solver
17.715, FE 17.7072 (0.041% / 0.044%).

One figure earns its place: predicted versus observed relative error
on Southwell 2.31, log–log, 18 points, 1:1 line, MATCH/WEAK/MISS
coded. Data: `p3_phaseD_inverse_ba_close.json`.

Do not cite a Vogel-versus-Southwell median ratio of 4.25 without its
distribution (p75 = 8.00, max 7160; 55 of 151 Vogel rows exceed 3) and
the uncontrolled \(b/a\) confound (Vogel’s excess is worst at small
\(b/a\); Southwell’s worst rows sit at large \(b/a\)).

---

## 7. In-plane free–free ring

Irie 1984 Table 2, \(\nu=0.3\), \(\beta=0.2,0.4,0.6,0.8\): **48/48
MATCH**, max relative miss 0.22% (\(\beta=0.8\), \(n=2\) first). Both
\(n=0\) radial and torsional recovered. Frequency-only; no IP MAC is
promised (Paper 2’s in-plane MAC was ugly for the same reason).

IP disk, if published, is Irie’s \(\beta=0.01\) tiny-hole column, not
claimed as \(R_i=0\) unless the 2×2 agrees. OOP disk already exists, so
the paper is not “we refused to do disks in flexure.”

---

## 8. Discussion

What travels from Papers 1–2: Frobenius radial content, conversion
hygiene, FE as an independent method, frozen screens that are not
retuned onto a table.

What does not: the sector residual bar, SUBDOM, Screen B, the
Kirchhoff corner jump (there is no \(\theta\)-edge).

Named limits:

- \(\nu\) on some SP-160 tables is \(1/3\); Narita and FE are 0.30 or
  0.33. Comparisons are not mixed.
- \(n=4\) is not in Table 2.35; FE carries it as (4,0) at \(b/a=0.5\).
- Table 2.30 \(n=0\) at \(b/a=0.9\) is a source defect (72/76).
- Southwell 2.31 is 13/18 direct; FE, not a bespoke inverse metric,
  is what settles the residue.
- Disk tiny-hole \(10^{-3}\) control remains a named FAIL.
- Window-edge values are labelled, not quoted as frequencies.
- C–C, remaining Vogel combinations, and polar-orthotropic F–F
  (Narita Table 3) are reachable with little or no new derivation;
  they are not in this draft’s production tables.

---

## 9. Conclusions

A validated annular-sector solver reduces, without a new ODE, to the
classical ring and disk. The \(4\times4\) recovers Vogel’s free–free
and mixed-edge tables and Irie’s in-plane free–free table. The
\(2\times2\) recovers the solid disk; the \(4\times4\) at \(R_i=0\)
does not. Independent FE confirms the flexural spine, including the
cell where Southwell prints 81 and the plate vibrates at 133.
Piezoelectric identification, using the F–F and C–F ring as a forward
model, is later work.

---

## Supplementary (planned)

- S.1 Bessel span of the four Frobenius branches.
- S.2 Disk regularity, inner jet, and the rank-3 trap.
- S.3 Conversion algebra: one worked FF-P1 number, one worked Irie
  number.
- S.4 Full match lists (2.35, 2.5, 2.22, 2.30, 2.31, Irie 2).
- S.5 Cluster job provenance (not in main text).
- S.6 Southwell inverse-\(b/a\) table and the predicted-versus-observed
  figure.

---

## Graphical abstract (planned)

Ring mode \((2,0)\) plus disk \((2,0)\), one \(\lambda^2(b/a)\) curve
from Table 2.34.

---

## Internal ledger (not for the submitted PDF)

| phase | result | where |
|---|---|---|
| 0 | API + TestRingDisk; mill 8/8 | §18.112, §18.124 |
| A | 2.35 80/80; Narita 2(a) hold-out | §18.114, job 2461503 |
| B | Narita 2(b) 4/4; 2.5 6/6; G3 FAIL named | §18.115, job 2463651 |
| C | FE 11/11 QUEUE_OK; n-ID closed | §18.116, jobs 2461506/07 |
| D | 2.22 80/80; 2.30 72/76 paper / 72/78 mill; Southwell 13/18; mixed-edge FE | §18.117–§18.123, jobs 2463845, 2469111, 2469112 |
| E | Irie T2 48/48 | §18.114, job 2461505 |
| G | this draft | — |

Do not quote \(\lambda^2=92.566\). Do not loosen PHASE_D: FAIL, G3
\(10^{-3}\), or the 1%/3% bar. Do not bump `SOLVER_VERSION`. Do not
edit Papers 1–2.

**Next:** LaTeX in `github_repo/paper/PAPER3_RING_DISK_DRAFT.tex` +
`PAPER3_RING_DISK_SUPPLEMENTARY.tex`, Paper 1/2 preamble. Then the
Southwell predicted-versus-observed figure from
`p3_phaseD_inverse_ba_close.json`. Then a Paper-3-specific review
prompt cloned from `PAPER2_REVIEWER_PROMPT.md`.
