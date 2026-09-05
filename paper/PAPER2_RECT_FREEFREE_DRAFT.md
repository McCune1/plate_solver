# Paper 2 draft — Completely free rectangular plates

**Status (2026-09-04).** Working manuscript, not LaTeX, not submitted.
Out-of-plane (flexure) is validated against independent ANSYS at five
aspect ratios and is ready to draft. In-plane (extension) is **not
started** (no `bc` on `RectIPAssembler`, no FE decks, no Bardell
comparison). The intended paper is **one paper covering both motions**,
matching Paper 1’s annular structure. This file is the OOP half plus an
honest IP placeholder.

Live location: `github_repo/paper/PAPER2_RECT_FREEFREE_DRAFT.md`.
Do not treat this as part of Paper 1 (`PAPER1_FREEFREE_DRAFT.tex`).

Job numbers belong in a later Supplementary provenance list, not in the
submitted main text. They are kept here only while the draft is internal.

---

# Exact Wave-Function Solution of Completely Free Rectangular Plates: Four-Corner Kirchhoff Jump and Mode Screening

**Garrek McCune and Daniel Stutts**
Mechanical and Aerospace Engineering, Missouri University of Science
and Technology, Rolla, MO 65409, USA

*Companion to Paper 1 (annular-sector FFFF). Draft, OOP half complete.*

## Abstract

*(OOP written; IP sentences marked.)*

The Seok–Tiersten–Scarton wave-function boundary-determinant method for
rectangular plates is extended from the cantilever (clamped-free) to
the completely free (FFFF) configuration. Closed-form wave functions
satisfy the governing plate equation exactly; a two-edge variational
statement on the long edges, together with free conditions on the short
edges, produces a finite characteristic determinant whose zeros are
candidate natural frequencies. For out-of-plane (flexural) motion the
Kirchhoff twisting-moment jump at a free–free corner is not optional:
the naive four-corner sum of the cantilever-style terms vanishes
identically by a parity identity, and the physical jump is a
checkerboard combination
\((T-\nu)(-2)(TT-TB-WT+WB)\)
obtained from Seok/Tiersten/Scarton (2004) Eq. (16) on a closed contour
with outward normal. A persistence screen on a modest complex-pair
enlargement of the production basis (Screen B: three conjugate pairs,
dip travel \(|\Delta\Lambda|\le 0.01\)) separates modes that survive
that enlargement from a small family of known artifacts (a symmetric
leak near \(\Lambda=0.41\), antisymmetric extras near \(0.04\), \(0.46\),
and \(1.65\), and machine zeros at \(\Lambda=4\) and \(6\)).

Independent half-model SHELL281 finite-element calculations at five
aspect ratios \(\ell/b=1.0,1.5,2.0,2.5,3.0\), both even and odd
mirror parities, mesh-converged to \(<0.3\%\) (typically \(0.001\%\)),
confirm the screened flexural candidates. After collapsing
double-dips onto the same FE partner, the refined production-basis
matches lie mostly inside \(1\%\), often a few tenths of a percent.
Two low symmetric modes are invisible to the production basis
(no complex pairs) and appear at \(0.00\%\) miss once the persist
basis is used; one symmetric target at \(\ell/b=1.0\)
(\(\Lambda_{\mathrm{FE}}=1.982\)) remains unmatched. Default
clamped-free assembly is unchanged.

*[IP to add: in-plane FFFF assembly has no Kirchhoff corner jump;
validation against ANSYS and against Bardell, Langley and Dunsdon
(1996) at \(a/b=1,2\).]*

**Keywords:** rectangular plate; free vibration; exact solution;
free boundary; Kirchhoff corner; spurious modes.

## 1. Introduction

Seok, Tiersten and Scarton (2004) gave exact-wave-function solutions
for cantilevered rectangular plates, out-of-plane (Part 1) and
in-plane (Part 2). The method is the rectangular counterpart of their
annular-sector papers: analytic wave functions in one direction, a
variational statement on the remaining pair of edges, natural
frequencies as zeros of a small characteristic determinant. Those
papers stop at clamped-free. Completely free plates are the harder
extension, because every edge is natural and Kirchhoff theory
requires a twisting-moment jump at each free–free corner.

A companion paper treats the annular-sector FFFF problem, including
the spurious-zero discrimination instruments that two weakly enforced
circumferential edges demand. The rectangular problem is not a
transplant of that work. The axial solution is closed-form rather
than Frobenius; the contour that produces the corner jump is a
rectangle with four free corners rather than a sector with two; and
the artifact mechanism that must be screened is a mixture of
basis-limited dips and a few \(\ell/b\)-independent leak families,
not the annular residual-screen taxonomy.

This paper (i) re-derives the four-corner flexural jump from
Part 1 Eq. (16) on a closed contour; (ii) shows why a same-sign
four-term sum is identically zero and why the surviving combination
is a checkerboard; (iii) defines a persistence screen that does not
retune a threshold to chase FE; (iv) reports FFFF flexural tables
at five aspect ratios against independent half-model finite
elements; (v) *[IP]* extends the in-plane assembler to four free
edges and compares with FE and with Bardell et al. (1996).

The rectangular cantilever tables of Parts 1 and 2 are already
reproduced by the same package (default `bc='clamped_free'`); that
path is not modified here.

## 2. Flexural formulation

### 2.1 Geometry and determinant

A rectangular plate of half-length \(\ell\) and half-breadth \(b\),
aspect ratio \(\ell/b\), thickness \(H\), isotropic Kirchhoff
flexure, Poisson ratio \(\nu=0.30\). The Seok/Tiersten/Scarton
nondimensional frequency \(\Lambda\) is used throughout. Finite-element
comparisons convert ANSYS cyclic frequencies by the factor
\(K=24.6644\,\mathrm{Hz}\) per unit \(\Lambda\), obtained from the
first non-rigid *VWRITE* Hz block of the half-model decks (the decks’
own \(\Lambda_{\mathrm{FE}}\) column is unusable: the APDL `SQRT` is
missing at runtime).

Even (SYM) and odd (ANTI) classes are the half-model mirror
conditions: SYM \(\leftrightarrow\) \(UY=ROTX=ROTZ=0\) on the cut;
ANTI \(\leftrightarrow\) \(UX=UZ=ROTY=0\). Production bases are
\(n_{\mathrm{real}}=6\), \(n_{\mathrm{cpair}}=0\) (SYM) and
\(5/1\) (ANTI). The persist screen uses three conjugate pairs on
both classes.

### 2.2 Four-corner jump

Part 1 Eq. (16)’s last integral is the twisting-moment contribution
\(\varphi=n_a t_{ab}s_b\) around a contour with outward normal,
oriented counterclockwise. At a free–free corner the integrand jumps.
The four corners of the rectangle are not equivalent under the local
\((n,s)\) frame: TT and WB are even under the global product of edge
tangents; TB and WT pick up a Jacobian \(-1\). The physical jump is
therefore the checkerboard

\[
\mathrm{corner}_{ff}=(T-\nu)(-2)\,(TT-TB-WT+WB).
\]

A naive sum of four cantilever-style terms,
\(p_{1t}q_{0t}+p_{1b}q_{0b}\) and cyclic, vanishes identically
because those two products sum to zero (parity identity). Under that
identity the checkerboard is twice the three-corner “FIX” combination
that had been used as an undervived stand-in. The default cantilever
assembly (three free edges, one clamped) is a different contour and
is not altered; `SOLVER_VERSION` is unchanged.

*[SM to write: line-by-line Eq. (16) trace, corner table TT/TB/WT/WB,
identity \(p_{1t}q_{0t}+p_{1b}q_{0b}\equiv 0\), numerical check that
clamped-free Table 3 is bit-stable.]*

### 2.3 Screen B

Discovery scans \(\sigma_{\min}(\Lambda)\) on the production basis.
Every listed dip with \(\sigma<0.3\) is re-scanned on the persist
basis (\(n_{\mathrm{cpair}}=3\), window \(\pm 0.02\), step \(0.002\)).
A dip **persists** if a persist-basis local minimum remains within
\(0.01\) in \(\Lambda\). That cut was frozen before the extra-lob
geometries were run and is not retuned.

Screen B kills a recurring antisymmetric pair near \(\Lambda=0.04\)
and \(0.46\), and labels a symmetric leak at \(0.41\) that persists
but does not match FE. It does not kill a deep extra near
\(\Lambda=1.65\) (tiny \(\sigma\), no FE partner) or a shallow
symmetric extra near \(1.92\). Those are reported as false
positives, not published modes.

Two genuine symmetric FE modes ( \(\ell/b=2.0\), \(\Lambda=0.543\);
\(\ell/b=3.0\), \(\Lambda=0.670\) ) never appear as production-basis
dips. On the persist basis they sit on the FE values at \(0.00\%\)
miss. Production discovery is therefore a floor, not a census. The
paper will state that limitation rather than enlarge the production
basis by default (adding pairs to every scan is the persist screen,
not a silent change of the published detector).

## 3. Flexural finite-element comparison

Half-model SHELL281, \(E=210\,\mathrm{GPa}\), \(\nu=0.30\),
\(\rho=7800\,\mathrm{kg\,m^{-3}}\), \(H=0.04\,\mathrm{m}\),
\(b=1\,\mathrm{m}\). Mesh 1 uses element size \(b/20\) along the
half-breadth; mesh 2 uses \(b/30\). Divisions along the length scale
as \(\mathrm{NDIVX}=\mathrm{round}(2(\ell/b)\times 20)\) and
\(\times 30\). Mesh 1 vs mesh 2 agrees to \(<0.02\%\) on the older
decks and \(0.001\%\) on \(\ell/b=2.0\) and \(3.0\). Rigid-body
(near-zero Hz) entries are discarded. Cross-parity matches are not
counted. At \(\ell/b=1.0\) (square plate) some SYM and ANTI FE
frequencies coincide; that degeneracy is physical.

A refined candidate is a **MATCH** if the production-basis local
minimum after a \(\pm 0.04\), step-\(0.002\) polish lies within
\(3\%\) of a same-parity FE \(\Lambda\), and Screen B still passes.
Double-dips that share one FE partner count once (closer kept).
High-\(\Lambda\) discovery uses a \(5\%\) coarse cut, then drops
machine zeros at \(4.000\) and \(6.000\) (\(\sigma\sim 10^{-18}\)
to \(10^{-22}\)) and collapses doubles; the tables below keep the
\(3\%\) subset unless noted.

### 3.1 Primary refined matches (\(\Lambda\lesssim 2.4\))

| \(\ell/b\) | class | \(\Lambda^*\) | \(\Lambda_{\mathrm{FE}}\) | miss | note |
| ---: | --- | ---: | ---: | ---: | --- |
| 1.0 | ANTI | 1.366 | 1.35288 | 0.97% | |
| 1.5 | SYM | 0.964 | 0.96347 | 0.06% | 1st SYM |
| 1.5 | SYM | 2.254 | 2.24373 | 0.46% | |
| 1.5 | ANTI | 0.906 | 0.89838 | 0.85% | 1st ANTI |
| 1.5 | ANTI | 2.124 | 2.07036 | 2.59% | weakest primary |
| 2.0 | SYM | 0.5434 | 0.54338 | 0.00% | persist basis only |
| 2.0 | SYM | 1.510 | 1.50735 | 0.18% | |
| 2.0 | SYM | 2.234 | 2.22558 | 0.38% | |
| 2.0 | ANTI | 0.674 | 0.66858 | 0.81% | 1st ANTI |
| 2.0 | ANTI | 1.482 | 1.47069 | 0.77% | |
| 2.5 | SYM | 0.348 | 0.34766 | 0.10% | 1st SYM |
| 2.5 | SYM | 0.978 | 0.96559 | 1.29% | |
| 2.5 | SYM | 1.888 | 1.88370 | 0.23% | |
| 2.5 | SYM | 2.278 | 2.27059 | 0.33% | |
| 2.5 | ANTI | 0.534 | 0.53126 | 0.52% | 1st ANTI |
| 2.5 | ANTI | 1.148 | 1.14003 | 0.70% | |
| 2.5 | ANTI | 1.918 | 1.90392 | 0.74% | |
| 3.0 | SYM | 0.242 | 0.24127 | 0.30% | 1st SYM |
| 3.0 | SYM | 0.6698 | 0.66984 | 0.00% | persist basis only |
| 3.0 | SYM | 1.320 | 1.31785 | 0.16% | |
| 3.0 | SYM | 2.162 | 2.15417 | 0.36% | |
| 3.0 | SYM | 2.256 | 2.24841 | 0.34% | |
| 3.0 | ANTI | 0.444 | 0.44044 | 0.81% | 1st ANTI; 0.46-family coincidence |
| 3.0 | ANTI | 0.936 | 0.93044 | 0.60% | |
| 3.0 | ANTI | 1.528 | 1.51698 | 0.73% | |
| 3.0 | ANTI | 2.258 | 2.24305 | 0.67% | |

First-mode travel is regular: SYM \(0.963\to 0.543\to 0.348\to 0.241\)
as \(\ell/b\) goes \(1.5\to 2.0\to 2.5\to 3.0\); ANTI
\(0.898\to 0.669\to 0.531\to 0.440\). At \(\ell/b=3.0\) the
historical ANTI \(0.46\) leak sits on the real first ANTI; it is
published as that mode with a footnote, not as evidence that the leak
is physical at other aspect ratios.

### 3.2 Higher flexural matches (selected, miss \(\le 3\%\))

| \(\ell/b\) | class | \(\Lambda\) | \(\Lambda_{\mathrm{FE}}\) | miss |
| ---: | --- | ---: | ---: | ---: |
| 1.0 | SYM | 2.460 | 2.45439 | 0.23% |
| 1.0 | SYM | 3.530 | 3.49250 | 1.07% |
| 1.0 | SYM | 6.190 | 6.16026 | 0.48% |
| 1.0 | SYM | 6.460 | 6.36553 | 1.48% |
| 1.0 | ANTI | 3.520 | 3.49250 | 0.79% |
| 1.5 | SYM | 5.380 | 5.38152 | 0.03% |
| 1.5 | ANTI | 3.860 | 3.83354 | 0.69% |
| 2.0 | SYM | 3.000 | 2.99812 | 0.06% |
| 2.0 | SYM | 3.650 | 3.62227 | 0.77% |
| 2.0 | ANTI | 2.580 | 2.55180 | 1.11% |
| 2.0 | ANTI | 4.060 | 4.02613 | 0.84% |
| 2.0 | ANTI | 6.200 | 6.17866 | 0.35% |
| 2.5 | SYM | 3.190 | 3.19522 | 0.16% |
| 2.5 | SYM | 4.180 | 4.13819 | 1.01% |
| 2.5 | SYM | 5.410 | 5.34819 | 1.16% |
| 2.5 | ANTI | 2.920 | 2.89206 | 0.97% |
| 2.5 | ANTI | 4.180 | 4.14962 | 0.73% |
| 2.5 | ANTI | 5.720 | 5.67275 | 0.83% |
| 2.5 | ANTI | 6.260 | 6.21965 | 0.65% |
| 2.5 | ANTI | 6.400 | 6.36131 | 0.61% |
| 3.0 | SYM | 2.480 | 2.45935 | 0.84% |
| 3.0 | SYM | 3.340 | 3.32528 | 0.44% |
| 3.0 | SYM | 3.660 | 3.62467 | 0.97% |
| 3.0 | SYM | 5.660 | 5.59746 | 1.12% |
| 3.0 | SYM | 6.180 | 6.14744 | 0.53% |
| 3.0 | ANTI | 3.160 | 3.14253 | 0.56% |
| 3.0 | ANTI | 4.260 | 4.23529 | 0.58% |
| 3.0 | ANTI | 5.560 | 5.51872 | 0.75% |
| 3.0 | ANTI | 6.220 | 6.18191 | 0.62% |
| 3.0 | ANTI | 6.380 | 6.34339 | 0.58% |

### 3.3 Open flexural items (report, do not hide)

- \(\ell/b=1.0\) SYM \(\Lambda_{\mathrm{FE}}=1.982\): still unmatched.
  A persist dip at \(2.08\) is \(4.9\%\) off; the second FE at \(2.454\)
  is matched separately by \(\Lambda=2.460\).
- Several FE modes near \(2.5\)–\(2.9\) fail Screen B even when a
  production dip is within ~1.5% (e.g. \(2.0\) SYM \(2.622\), \(3.0\)
  SYM \(2.887\), \(1.5\) SYM \(2.588\), \(2.5\) SYM \(2.467\)). The
  screen is not retuned to collect them.
- Not published: SYM \(0.41\) leak; ANTI \(0.04/0.46\) except the
  \(\ell/b=3.0\) coincidence above; ANTI \(\sim 1.65\); shallow SYM
  \(\sim 1.92\); \(\Lambda=4,6\) machine zeros.

## 4. In-plane motion (to be completed)

**Not started.** `RectIPAssembler` has no `bc` parameter; it is
cantilever-only. Expected structure, to be confirmed rather than
assumed:

- In-plane motion is second-order plane elasticity. There is no
  Kirchhoff twisting moment and no corner-jump analogue of §2.2
  (same argument as the annular in-plane case).
- Implementation is an edge-list change: all four edges natural,
  same unconstrained \(\sigma_{\min}\) assembly, plus a `bc`
  switch that leaves clamped-free bit-identical.
- Independent FE: full-plate or half-model PLANE183, same five
  aspect ratios if the OOP set is reused; mesh convergence required.
- Literature: Bardell, Langley and Dunsdon, *J. Sound Vib.* **191**
  (1996) 459–467, first six non-zero in-plane frequencies of
  isotropic F–F–F–F (and C–C–C–C) rectangles at \(a/b=1,2\),
  \(\Omega=\omega a/C_L\) with \(C_L^2=E/[\rho(1-\nu^2)]\). Paper 1
  already warns that the source’s (a)/(b) sub-labels appear swapped
  relative to its captions; any comparison must read the figures,
  not assume the labels.
- Discrimination: do not copy Screen B blindly. In-plane artifacts,
  if they appear, need their own pre-registered test (basis
  enlargement and/or a residual/MAC instrument). Paper 1’s SUBDOM
  is annular-IP-specific until shown otherwise.

Until those items exist, Paper 2 is not submittable. The OOP half
above is the manuscript they land in.

## 5. Discussion

The rectangular FFFF flexural problem is no longer an experimental
assembler. The four-corner jump is a derivation, not a fit; Screen B
is frozen; five aspect ratios of independent FE agree at the level
the annular flexural tables claimed. The remaining flexural issues
are named misses and named artifacts, which is the same honesty
standard Paper 1 used.

What this paper does *not* claim: that production-basis discovery is
complete; that Screen B is a portable residual screen in the annular
sense; that in-plane FFFF is validated; that clamped-free frequencies
changed (`SOLVER_VERSION` unbumped).

The natural place for this work is a companion paper, not an
expansion of the annular manuscript. The annular paper’s contribution
is the two-edge FFFF detector and its discrimination instruments on
a sector. This paper’s contribution is the four-free-edge rectangular
jump and the tables that follow from it, in both motions once IP is
done.

## 6. Conclusions *(flexure only, until §4 is filled)*

The Seok–Tiersten–Scarton rectangular method extends to completely
free flexure once the closed-contour Kirchhoff jump is assembled as
a checkerboard rather than a vanishing four-term sum. A frozen
complex-pair persistence screen, checked against half-model SHELL281
at \(\ell/b=1.0,1.5,2.0,2.5,3.0\), yields a set of even- and
odd-parity frequencies that match FE typically inside \(1\%\). Two
symmetric modes require the persist basis to appear; one square-plate
symmetric target remains open. In-plane FFFF is the other half of
the same paper and is not yet in hand. The elastic FFFF baseline on
both the annular and rectangular geometries is the intended foundation
for piezoelectric constitutive coupling and for using the exact
frequencies to identify those constants; that is later work, not a
claim of this paper.

## Internal provenance (strip before submission)

| item | evidence |
| --- | --- |
| checkerboard corner | `RectOOPAssembler` `bc='free_free'`; identity vs naive sum |
| Screen B | persist \(n_{\mathrm{cpair}}=3\), \(\|\Delta\Lambda\|\le 0.01\) |
| FE \(\ell/b=1.5\) | half-model, job 2398617 lineage |
| FE \(\ell/b=1.0,2.5\) | P5 Phase A queue |
| FE \(\ell/b=2.0,3.0\) | job 2454654 |
| primary refine 1.0/1.5/2.5 | job 2453874, 12 unique MATCH at 3% |
| primary refine 2.0/3.0 | job 2454701, 12 unique MATCH at 3% |
| persist-only SYM hits | job 2454702, both HIT at 0.00% |
| high-\(\Lambda\) 1.0/1.5/2.5 | job 2453875 |
| high-\(\Lambda\) 2.0/3.0 | job 2454703, zeros tagged ZERO_ART |
| G0 | `SOLVER_VERSION=2026-07-10.s10`, checkerboard present, old sum absent |

## Remaining work (order, tools)

Do these in order. OOP FE is done; do not retune Screen B; do not bump
`SOLVER_VERSION` unless a default clamped-free frequency changes.
Piezoelectric coupling is **not** on this list (later paper, after
the elastic FFFF baseline exists on both geometries).

1. **IP no-corner confirmation.** Show rectangular in-plane FFFF has
   no Kirchhoff jump (second-order plane elasticity). Tool:
   `grok-handoff` plus `GENERAL_BC_CORNER_DERIVATION.md` / Seok Part 2.
   Sandbox algebra, not a cluster job.
2. **`RectIPAssembler(bc='free_free')`.** Add `bc` without changing
   default cantilever assembly. Tools: `plate-solver-sandbox-probe`
   (CRLF, worker-path gate), `test_validated_tables.py` / clamped-free
   IP regression.
3. **Rectangular FFFF IP FE decks.** PLANE183, half- or full-model,
   same `ℓ/b` set as OOP if possible. Tools: copy annular IP deck
   pattern in `Ansys/NewAnsys/`; isolated manifest;
   `submit_ansys_queue.sh` discipline (one MAPDL seat). Convert Hz
   with a documented factor; do not trust a broken `Lambda_FE` column.
4. **IP discovery + screen + pair.** Production scan, a
   *pre-registered* IP persist/residual test (do not copy Screen B
   by default), pair to FE, then Bardell et al. (1996) at `a/b=1,2`
   reading the figures not the swapped captions. Tools:
   `plate-solver-cluster-probes`, `p5_rect_ff_lib.py`-style FE lists,
   Mill venv (`module load python/3.12.1`).
5. **Optional OOP leftover.** One persist-window on `ℓ/b=1.0` SYM
   FE 1.982, no Screen B retune. Cluster CPU, not ANSYS.
6. **Supplementary derivation.** Line-by-line Eq. (16) checkerboard
   (TT/TB/WT/WB, parity identity). Can overlap with 1–4. Tool: the
   existing `FFFF_VARIATIONAL_DERIVATION.md` / Paper 1 SM §S.2 style.
7. **Promote this markdown to LaTeX** (`PAPER2_RECT_FREEFREE_DRAFT.tex`
   + SM), strip job numbers from the main text, fill §4 from steps
   1–4. Compile in `github_repo/paper/`.
8. **Review pass** with `PAPER_REVIEWER_PROMPT.md`. Pipeline wiring
   (Phase D) is nice-to-have for IP sweeps, not a publication gate if
   the IP jobs stay single-process like the OOP overnight probes.
