# Paper 2 Supplementary Material draft

**Status (2026-09-05).** Promoted to LaTeX:
`github_repo/paper/PAPER2_RECT_FREEFREE_SUPPLEMENTARY.tex` (7 pp).
This markdown is the working notes, not the manuscript. S.1 corner,
S.2 conversion/Bardell, S.3 92-row list, S.4 job provenance.

Live location (draft): `github_repo/paper/PAPER2_RECT_FREEFREE_SUPPLEMENTARY_DRAFT.md`.

---

## S.1.1 The exact-differential term of Eq. (16)

Seok, Tiersten and Scarton's rectangular Part 1 Eq. (16) is a closed-contour
weak-form identity on the plate boundary $\partial\Omega=c_C\cup c_N$
(clamped and natural portions). Written out in full,

$$
\int_{c_C}\Bigl\{u_3^{(0)}\,\delta\bigl(n_a t_{a3}^{(0)}\bigr)
  - u_{3,b}^{(0)}\,\delta\bigl(n_a t_{ab}^{(1)}\bigr)\Bigr\}\,ds
- \int_{c_N}\Bigl\{\bigl(n_a t_{a3}^{(0)}+(n_a t_{ab}^{(1)} s_b)_{,s}\bigr)\,
  \delta u_3^{(0)} - n_a t_{ab}^{(1)} n_b\,\delta u_{3,n}^{(0)}\Bigr\}\,ds
+ \int_{c_N}\bigl(n_a t_{ab}^{(1)} s_b\,\delta u_3^{(0)}\bigr)_{,s}\,ds = 0 .
$$

Only the third integral is an exact differential in the arclength $s$.
Writing $\varphi := n_a t_{ab}^{(1)} s_b$ (the traction-moment contraction
that reduces to the twisting moment $M_{xy}$ on a smooth edge),

$$
\int_{c_N}(\varphi\,\delta u_3^{(0)})_{,s}\,ds
= \sum_{\text{corners}}\Bigl[\varphi\,\delta u_3^{(0)}\Bigr]_{\text{end}}
  - \Bigl[\varphi\,\delta u_3^{(0)}\Bigr]_{\text{start}},
$$

the sum running over every free–free corner where two arcs of $c_N$ meet;
interior-arc contributions cancel against the $\varphi_{,s}$ term already
present in the second integral (the usual Kirchhoff effective-shear
reduction). A completely free rectangle has $c_C=\varnothing$ and four such
corners.

**Single-corner evaluation and the classical anchor.** Place one corner at
the local origin with both positive local axes pointing into the plate, so
the two free edges meeting there are the local $+x$ and $+y$ half-axes.
Consistent counterclockwise orientation of $\partial\Omega$ gives, on the
vertical edge ($x=0^+$, $s$ increasing in $+y$), $\mathbf n=+\mathbf e_x$,
$\mathbf s=+\mathbf e_y$, so $\varphi_V=t_{xy}^{(1)}$; on the horizontal edge
($y=0^+$, $s$ increasing in $-x$), $\mathbf n=+\mathbf e_y$,
$\mathbf s=-\mathbf e_x$, so $\varphi_H=-t_{yx}^{(1)}$. Because
$t_{ab}^{(1)}$ is symmetric, $\varphi_V-\varphi_H=2\,t_{xy}^{(1)}$, so the
corner's contribution to the first variation is
$+2\,t_{xy}^{(1)}\,\delta u_3^{(0)}$. In the isotropic Kirchhoff limit
$t_{xy}^{(1)}$ is the classical twisting moment $M_{xy}$ (up to the
project's own overall sign convention, fixed as in the annular free–free
derivations by the sign of $\int t_{ab}^{(1)}\delta\kappa_{ab}$), so the
natural condition is $2M_{xy}=0$ — Timoshenko's concentrated corner
reaction $R=2M_{xy}$. This single-corner trace and its classical-limit
check were carried out and PASSED against that pre-registered criterion
before any four-corner assembly was attempted (project file
`GrokCode/archive/by_topic/rect_freefree/RectangularFF.txt`); the
project's existing free–free derivations (annular, variational and Green
forms) all fix the same overall coefficient once $\varphi$'s sign is
fixed, so the same $-2$ is retained here:

$$
\text{corner contribution} = -2\,t_{xy}^{(1)}\,\delta u_3^{(0)}\Big|_{\text{corner}} .
$$

## S.1.2 Four corners, and the discrete trial/test evaluation

In this project's discrete basis, the twisting-moment factor at a corner is
the trial function's mixed second derivative $\partial^2 w/\partial x_1
\partial x_2$ there, and the paired test factor is the test function's
value $w$ there — confirmed directly (`mp.diff`, ~40-digit agreement) that
`U(xi,r,m,x1)` is exactly the $m$-th $x_1$-derivative of `U(xi,r,0,x1)` and
`peval(psi[m],x2)` is exactly the $m$-th $x_2$-derivative of the $x_2$
profile, so every corner term is literally $M_{xy}^{\text{trial}}
(\text{corner})\cdot w^{\text{test}}(\text{corner})$, matching the already-validated
`clamped_free` tip-corner term's own structure. Define, for a given trial
branch $\xi'$ and test branch $\xi''$,

$$
p_{1t}=\left.\frac{\partial\psi_1'}{\partial x_2}\right|_{x_2=+\pi/2},\quad
p_{1b}=\left.\frac{\partial\psi_1'}{\partial x_2}\right|_{x_2=-\pi/2},\qquad
q_{0t}=\left.\psi_0''\right|_{x_2=+\pi/2},\quad
q_{0b}=\left.\psi_0''\right|_{x_2=-\pi/2},
$$

(`p1t/p1b/q0t/q0b` in the code, `core_solvers.py::RectOOPAssembler.assemble`)
and, with $U_{t1},V_{t0}$ the $x_1=+L$ (tip) trial/test factors and
$U_{w1},V_{w0}$ the $x_1=-L$ (wall) factors, the four raw corner terms

$$
\mathrm{TT}=p_{1t}U_{t1}q_{0t}V_{t0},\quad
\mathrm{TB}=p_{1b}U_{t1}q_{0b}V_{t0},\quad
\mathrm{WT}=p_{1t}U_{w1}q_{0t}V_{w0},\quad
\mathrm{WB}=p_{1b}U_{w1}q_{0b}V_{w0}.
$$

TT and WB share the local orientation used in S.1.1 (outward normal along
the trial axis, `+2` sense); TB and WT sit at the diagonal corners, where
the same closed-contour threading picks up the opposite sense, because the
inward-axis Jacobian there is $-1$ (the $F_{\text{in}}-F_{\text{out}}=-2t_{12}$
statement in the deployed code's own comment) — the discrete restatement
of the classical alternating-sign corner-reaction pattern at the four
corners of a free rectangular plate.

## S.1.3 The parity identity, in both symmetry classes

The half-model SYM/ANTI mirror conditions fix the phase of the $m=0$ and
$m=1$ $x_2$-profiles at $x_2=\pm\pi/2$ (`phin`$=\pi/2$ for SYM, $0$ for
ANTI; each profile a sum of $c\sin(e\,x_2+p)$ terms). Direct evaluation of
that phase algebra (`P5_PHASE_A_SIGN_AUDIT_2026-09-02.md` §"Step 1"):

- **SYM** ($m=1$, extra $\pi/2$ phase flips sign under $x_2\to-x_2$;
  $m=0$ does not): $p_{1b}=-p_{1t}$ exactly, $q_{0b}=q_{0t}$ exactly.
- **ANTI** (roles of the two profiles swap): $p_{1b}=p_{1t}$ exactly,
  $q_{0b}=-q_{0t}$ exactly.

Both were confirmed numerically to machine/mp precision (the bracket
$p_{1t}q_{0t}+p_{1b}q_{0b}$ reproduces $0$ to $\sim 10^{-27}$ at
`dps=25`, independent of branch, $\Lambda$, $\ell/b$ or $\nu$). **The same
two relations force $\mathrm{TB}=-\mathrm{TT}$ and $\mathrm{WB}=-\mathrm{WT}$
in *both* classes**, term by term rather than only in the combined bracket
$p_{1t}q_{0t}+p_{1b}q_{0b}\equiv 0$ the audit stated: in SYM,
$\mathrm{TB}=p_{1b}U_{t1}q_{0b}V_{t0}=(-p_{1t})U_{t1}(q_{0t})V_{t0}=-\mathrm{TT}$,
and identically $\mathrm{WB}=-\mathrm{WT}$; in ANTI,
$\mathrm{TB}=p_{1t}U_{t1}(-q_{0t})V_{t0}=-\mathrm{TT}$ and again
$\mathrm{WB}=-\mathrm{WT}$. This stronger, unified statement (both classes
give the identical pairwise relation, not just a vanishing sum) is the
bridge this Supplementary section supplies between the audit's identity
and the final formula below.

**Consequence 1 — why the original four-term same-sign sum is dead.** The
first `free_free` corner term ever coded was the naive
$\mathrm{TT}+\mathrm{TB}+\mathrm{WT}+\mathrm{WB}$. Substituting
$\mathrm{TB}=-\mathrm{TT}$, $\mathrm{WB}=-\mathrm{WT}$ gives
$\mathrm{TT}-\mathrm{TT}+\mathrm{WT}-\mathrm{WT}\equiv 0$ in both classes —
confirmed numerically to $1.6\times10^{-24}$ (SYM) and
$2.8\times10^{-23}$ (ANTI) on real assembled matrices, i.e. exactly zero to
precision, not merely small. This is why `corner_ff` as originally written
contributed nothing to `K` regardless of $\Lambda$, $\ell/b$, or branch —
a vanishing combination, not a wrong sign.

**Consequence 2 — the surviving checkerboard is exactly twice the
single-point shorthand.** The deployed corner term flips the sign of the
diagonal pair instead of summing everything the same way:

$$
\mathrm{corner}_{ff} = (T-\nu)(-2)\,(\mathrm{TT}-\mathrm{TB}-\mathrm{WT}+\mathrm{WB}).
$$

Substituting the same two relations: in SYM,
$\mathrm{TT}-\mathrm{TB}-\mathrm{WT}+\mathrm{WB}
=\mathrm{TT}-(-\mathrm{TT})-\mathrm{WT}+(-\mathrm{WT})
=2(\mathrm{TT}-\mathrm{WT})$; in ANTI the identical substitution gives the
identical result, $2(\mathrm{TT}-\mathrm{WT})$. So the checkerboard bracket
equals **twice** the earlier single-point candidate
$\mathrm{TT}-\mathrm{WT}$ (the tip-top-only, minus wall-top-only pattern
that had been tried as a structural analogy to the validated
`clamped_free` tip term, before this closed-contour trace existed) — in
*both* symmetry classes, by the same two-line substitution, not by
coincidence in one class only. This is the origin of the checkerboard
label: TT and WB (the "diagonal" pair under the local corner frames of
S.1.2) keep one sign, TB and WT (the other diagonal) flip, and the parity
identity collapses the four-term difference to twice a two-term one.

## S.1.4 Numerical checks (not reasoning alone)

- **Regression on the unmodified path.** The `clamped_free` branch of
  `RectOOPAssembler.assemble` is untouched by this term (it uses its own
  `corner`, not `corner_ff`); a subclassed assembler with `corner_ff`
  swapped in for diagnostics reproduces `clamped_free` byte-for-byte
  against the parent (`max|K-K_parent|=0.0` at the validated Table-3 SYM
  benchmark, $\ell/b=1.5$, $\Lambda=0.1551$). `SOLVER_VERSION` is
  unchanged; only the `free_free` corner formula moved.
- **Algebraic gate.** Before any FE comparison was trusted, the
  candidate three-way comparison (`CURRENT` = the dead naive sum,
  `FIX` = the single-point $\mathrm{TT}-\mathrm{WT}$ candidate,
  `DERIVED` = the checkerboard above) required, as a pre-registered
  pass/fail gate: `CURRENT` genuinely $\equiv 0$; `FIX`/`DERIVED`
  genuinely nonzero; and
  $K_{\mathrm{DERIVED}}-K_{\mathrm{CURRENT}} = 2\,(K_{\mathrm{FIX}}-K_{\mathrm{CURRENT}})$
  to $\sim10^{-10}$ relative — the matrix-level restatement of S.1.3's
  algebra. All three passed
  (`probe_rect_ff_oop_derived_vs_fix_2026-09-03.py`, Part A).
- **Independent finite-element cross-check.** At the production basis,
  the checkerboard's nearest-to-FE dip beat both the naive-sum term (which
  contributes nothing, so its dip location is unaffected by any corner
  physics) and the single-point `FIX` candidate at every finite-element
  anchor available at the time (four points across $\ell/b=1.5, 2.5$,
  SYM and ANTI; see `LESSONS_LEARNED.md` §18.83 for the FE decoding and
  the same probe's Part B for the comparison). This is the basis for
  adopting the checkerboard in `core_solvers.py` (2026-09-03) and for
  §3's independent five-aspect-ratio FE table, which is a later,
  larger-sample confirmation of the same formula, not a re-litigation of
  it.

## S.1.5 What this section does not claim

The single-corner Eq. (16) trace (S.1.1) is a first-principles derivation,
checked against the classical $R=2M_{xy}$ limit. The extension from one
corner to the specific four-corner checkerboard sign pattern (S.1.2) is
argued from the closed-contour orientation logic plus the classical
alternating-corner-reaction pattern, and is *validated* rather than
independently re-derived corner-by-corner from Eq. (16) a second time —
the same status Paper 1's annular corner term has relative to its own
Green-identity cross-check. The parity identity (S.1.3) is a fully
first-principles algebraic fact about this project's specific $(m=0,m=1)$
trial-function evaluations, exact in both classes, and is not a numerical
coincidence at one geometry. Nothing here is a claim about the in-plane
formulation, which has no Kirchhoff jump term at all (Part 2 Eq. (1)/(36),
confirmed 2026-09-05; see main draft §4).

---

## S.2 In-plane nondimensionalization, Bardell mapping, and match provenance

Main-text §4 reports in this project's \(\bar\Omega=\omega/\bar\omega\).
Job numbers and the 92-row list stay here, not in the submitted main
text.

### S.2.1 This project's \(\bar\Omega\)

Seok/Tiersten/Scarton Part 2 Eqs. (22)–(23):
\(\bar\omega=(\pi/(2b))\sqrt{c_{66}/\rho}\), \(c_{66}=E/(2(1+\nu))\),
\(\bar\Omega=\omega/\bar\omega\). For the FE decks
(\(b=1\,\mathrm{m}\), \(E=210\,\mathrm{GPa}\), \(\nu=0.30\),
\(\rho=7800\,\mathrm{kg\,m^{-3}}\)) this is
\(\bar\Omega=f[\mathrm{Hz}]/804.481\), independent of \(\ell/b\).
Half-model PLANE183, job 2455569, mesh-conv max \(0.004\%\).

### S.2.2 Bardell, Langley and Dunsdon (1996)

N. S. Bardell, R. S. Langley and J. M. Dunsdon, *J. Sound Vib.*
**191**(3) (1996) 459–467. Source: `Project Knowledge/papers/Bardell.pdf`
(a real PDF; `read-plate-papers`' unzip step does not apply).

Their parameter is \(\Omega_B=\omega a/C_L\) with
\(C_L^2=E/[\rho(1-\nu^2)]\), frequencies printed under Figure 1
(F–F–F–F) rather than tabulated. Poisson ratio is not stated
numerically. Table 1's exact S–S–S–S square pair
\((0,1)=(1,0)=1.859\) equals \(\pi\sqrt{(1-\nu)/2}=1.8586\) at
\(\nu=0.30\).

Plates identify as \(a=2\ell\), so \(a/b=\ell/b\) and
\(\Omega_B=\bar\Omega\cdot\pi\cdot(\ell/b)\cdot\sqrt{(1-\nu)/2}\).

Figure 1 caption: “(a) \(a/b=1\); (b) \(a/b=2\)”. Panel (a) is drawn
2:1 with unique frequencies \(\{1.954, 2.961, 3.267, 4.726, 4.784,
5.205\}\); panel (b) is square with the repeated pair \(2.472, 2.472\).
Body text (p. 462) assigns repeats to \(a/b=1\). Mapping follows the
rendered proportions, not the caption letters (Paper 1 §7 /
`LESSONS_LEARNED.md` §18.72a). Table 1 of the source does not have
this swap.

Printed F–F–F–F first six, assigned that way:

- \(a/b=1\) (square, Fig. 1(b)): \(2.321\), \(2.472\), \(2.472\),
  \(2.628\), \(2.987\), \(3.452\).
- \(a/b=2\) (2:1, Fig. 1(a)): \(1.954\), \(2.961\), \(3.267\),
  \(4.726\), \(4.784\), \(5.205\).

All twelve convert onto matched \(\bar\Omega^*\) in main-text §4.3
(solver max miss \(1.69\%\); independent FE max miss \(0.018\%\)).

### S.2.3 92/92 unique MATCH provenance

Window \(\bar\Omega\in[0.02,2.50]\), five aspect ratios, both
parities. Count \(5+4+8+6+10+7+13+10+16+13=92\) is the
`p5_rect_ff_lib.IP_FE_LISTS` entry total.

| job | role | unique MATCH added |
| --- | --- | ---: |
| 2455569 | FE ground truth (PLANE183 half-model) | — |
| 2455570 | first discovery, \(n_{\mathrm{cpair}}=0\), 5% cut | 38 |
| 2456002 | production re-discovery, \(n_{\mathrm{cpair}}=3\), 3% cut | 76 (of which 36 reconfirm 2455570) |
| 2456013 | mop-up, \(n_{\mathrm{cpair}}=5\), \(\pm 0.05\) on the remaining 16 FE values | +16 → 92/92 |

Full 92-row list: `PAPER2_IP_MATCH_TABLE.md` in this folder.
`SOLVER_VERSION=2026-07-10.s10`. Screen B was not retuned; IP does
not copy it. Two 2456002 “losses” relative to 2455570
(\((2.0,\mathrm{ANTI},1.27148)\), \((3.0,\mathrm{SYM},1.66346)\))
are both present in the 92 (root-crowding, recovered by 2456013).
