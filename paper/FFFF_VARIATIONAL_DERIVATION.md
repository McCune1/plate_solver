# The free-free out-of-plane variational functional: a line-by-line derivation from Ref. [seok1]

**Status (2026-08-13): condensed version is Supplementary Material §S.2**
of `PAPER1_FREEFREE_SUPPLEMENTARY.tex` in this directory.
This file remains the repository companion with the full line-by-line
algebra. It is not the paper.

**Earlier status (2026-07-19): condensed version folded into `PAPER1_FREEFREE_DRAFT.tex` Appendix B**
(that appendix later moved to the SM file; do not look for Appendix B in the current main text).
(first folded 2026-07-18, updated same night for §7a, updated again 2026-07-19 for §8's full
substitution). This document does what Appendix B of the paper originally said was *not yet
done* — carry the variational derivation from Seok & Tiersten's general functional (their
Eq. (6)) through to a free-free analogue of their Eq. (43), rather than asserting the
corner-term result by pattern-matching alone. Both gaps flagged in the first draft (the
corner-contour multiplicity, §7a; the constitutive/dimensionless substitution, §8) are now
closed — see §9 for the current honesty summary. Every displayed equation is labeled with its
confidence level:

- **[TRANSCRIBED]** — copied directly from the published page image (`seok1`, pp. 759–768),
  cross-checked against the OCR text and re-verified twice in this session.
- **[DERIVED, high confidence]** — obtained here by direct, checkable algebra from a
  [TRANSCRIBED] source; the steps are shown in full so the algebra can be checked independently.
- **[DERIVED, moderate confidence]** — obtained here by an argument I believe is correct and
  have checked as carefully as I can by hand, but which depends on sign bookkeeping (normal/
  tangent orientation through an indicial contour integral) that is genuinely easy to get wrong
  silently. Flagged explicitly rather than asserted as certain.
- **[MATCHED TO VALIDATED CODE]** — the final sign/coefficient is fixed by requiring agreement
  with the already-validated, already-deployed `FreeFreeOOP` implementation (8/8 modes confirmed
  against independent FE benchmarks) and the parity-cancellation proof already in Appendix B,
  rather than claimed as an independent from-scratch re-derivation. This is not circular: it is
  the same "validate by direct numerical probe, not reasoning alone" discipline the project
  already applies everywhere else — a hand derivation
  that reproduces a computationally-confirmed answer is corroborating evidence for the
  derivation, not the reverse.

## 0. What this buys, precisely

Appendix B's existing text (§B, "Free-edge functional and parity cancellation") already asserts
the corner-coefficient result and grounds it in a physical analogy (the classical free-corner
reaction $2M_{r\theta}$). What was *not* there: (a) the actual general functional the coefficients
come from, (b) why the two circumferential edges are asymmetric *in that functional*, not just in
the already-specialized Eq. (43), and (c) an explicit construction of the free-free replacement
term, rather than an assertion that one exists. This document supplies (a)–(c), and is honest
about exactly where it can and cannot claim full first-principles certainty.

## 1. Geometry and edge-type convention (terminology check)

Ref. [seok1] calls the two arcs $r=R_i$, $r=R_o$ the **circumferential edges**, and the two
straight lines $\theta=\pm\Theta$ the **radial edges** — the opposite of what "circumferential"
might suggest at first read; I confirmed this twice against the paper's own text (p. 758: "the
plate is fixed on one radial edge, free on the other and free on both circumferential edges";
p. 765, Eq. (23)'s heading, "free edge conditions at $r=R_i,R_o$"). In the paper's cantilever
case: $r=R_i,R_o$ are **both free**; $\theta=-\Theta$ is **clamped** (the "wall"); $\theta=\Theta$
is **free**. `FreeFreeOOP` releases $\theta=-\Theta$ to free as well, leaving all four edges free.

## 2. The general functional (Ref. [seok1] Eq. (1a) and Eq. (6)) — [TRANSCRIBED]

The starting point, before any specific edge type is assumed, is the 3-D variational statement
reduced to the plate's flexural part (Ref. [seok1] Eq. (1a), their Ref. [10] Eq. (16)):

$$
\int_S dS\Big[\big(t^{(0)}_{\alpha3,\alpha} - 2\rho h\ddot u_3^{(0)}\big)\delta u_3^{(0)}
+\big(t^{(1)}_{\alpha\beta,\alpha}-t^{(0)}_{3\beta}\big)\delta u_\beta^{(1)}\Big]
$$
$$
+\int_{c_C}\!ds\big\{u_3^{(0)}\delta(n_\alpha t^{(0)}_{\alpha3}) - u_{3,\beta}^{(0)}\delta(n_\alpha t^{(1)}_{\alpha\beta})\big\}
$$
$$
-\int_{c_N}\!ds\Big\{\big(n_\alpha t^{(0)}_{\alpha3}+(n_\alpha t^{(1)}_{\alpha\beta}s_\beta)_{,s}\big)\delta u_3^{(0)} - n_\alpha t^{(1)}_{\alpha\beta}n_\beta\,\delta u_{3,n}^{(0)}\Big\}
+\int_{c_N}\!ds\big(n_\alpha t^{(1)}_{\alpha\beta}s_\beta\,\delta u_3^{(0)}\big)_{,s} = 0.
\tag{1a}
$$

$c_C$ is the displacement-prescribed (clamped) portion of the boundary, $c_N$ the
traction-prescribed (free) portion, $n_a,s_b$ the outward normal / counterclockwise tangent.
Specializing to the annular sector's cylindrical coordinates (with $u_3^{(0)}\to w$) gives the
paper's Eq. (6), which I re-transcribed directly from the page image:

$$
\int_S dS\Big[\Big\{\tfrac{\partial\tau^{(0)}_{rz}}{\partial r}
+\tfrac1r\big(\tfrac{\partial\tau^{(0)}_{\theta z}}{\partial\theta}+\tau^{(0)}_{rz}\big)-2\rho h\ddot w\Big\}\delta w
+\Big\{\tfrac{\partial\tau^{(1)}_{rr}}{\partial r}+\tfrac1r\tfrac{\partial\tau^{(1)}_{r\theta}}{\partial\theta}
+\tfrac{\tau^{(1)}_{rr}-\tau^{(1)}_{\theta\theta}}{r}-\tau^{(0)}_{rz}\Big\}\delta u_r^{(1)}
$$
$$
+\Big\{\tfrac{\partial\tau^{(1)}_{r\theta}}{\partial r}+\tfrac1r\big(\tfrac{\partial\tau^{(1)}_{\theta\theta}}{\partial\theta}+2\tau^{(1)}_{r\theta}\big)-\tau^{(0)}_{\theta z}\Big\}\delta u_\theta^{(1)}\Big]
-\int_{-\Theta}^{\Theta}\Big[r\Big\{\big(\tau^{(0)}_{rz}+\tfrac1r\tfrac{\partial\tau^{(1)}_{r\theta}}{\partial\theta}\big)\delta w-\tau^{(1)}_{rr}\delta\big(\tfrac{\partial w}{\partial r}\big)\Big\}\Big]_{r=R_i}^{r=R_o}d\theta
$$
$$
-\int_{R_i}^{R_o}\!\Big[w\,\delta\tau^{(0)}_{\theta z}-\tfrac{\partial w}{\partial r}\delta\tau^{(1)}_{r\theta}\Big]_{\theta=-\Theta}\!dr
-\int_{R_i}^{R_o}\!\Big[\big(\tau^{(0)}_{\theta z}+\tfrac{\partial\tau^{(1)}_{\theta r}}{\partial r}\big)\delta w-\tau^{(1)}_{\theta\theta}\delta\big(\tfrac1r\tfrac{\partial w}{\partial\theta}\big)\Big]_{\theta=\Theta}\!dr
$$
$$
-\big[\tau^{(1)}_{r\theta}\delta w\big]^{r=R_o,\theta=-\Theta}_{r=R_i,\theta=-\Theta}
+2\big[\tau^{(1)}_{r\theta}\delta w\big]^{r=R_o,\theta=\Theta}_{r=R_i,\theta=\Theta}=0.
\tag{6}
$$

**Why the two circumferential — sorry, *radial* — edges already look different here, before any
edge type is chosen for $\theta=-\Theta$:** immediately after Eq. (6), the paper states (p. 762,
verified twice): *"the last two terms in Eq. (6) were obtained by performing the line integral
around $c_N$ ... along the three free edges."* Three, not four — because in the paper's own
cantilever setup $\theta=-\Theta$ is the wall and is *excluded* from $c_N$. The $\theta=-\Theta$
surface term (primal $w$, $\partial w/\partial r$ multiplying a **variation of the dual** stress
quantities $\delta\tau^{(0)}_{\theta z}$, $\delta\tau^{(1)}_{r\theta}$) is the $c_C$-type mixed/
Lagrangian contribution; the $\theta=\Theta$ term (dual stress quantities multiplying a
**variation of the primal** $\delta w$, $\delta(\partial w/\partial\theta)$) is the $c_N$-type
natural contribution. **This is the actual origin of the asymmetry Appendix B already asserts —
traced now to Eq. (6) itself, not just to the already-specialized Eq. (43).**

## 3. What stays exactly the same for FFFF — [DERIVED, high confidence]

Releasing $\theta=-\Theta$ from $c_C$ to $c_N$ changes *only* the boundary terms of Eq. (6) that
are literally evaluated at $\theta=-\Theta$. It changes nothing about:

- The **Euler–Lagrange equations** (Ref. [seok1] Eqs. (9a)–(9c)), since these come from
  requiring the area-integral coefficient of each *independent* variation to vanish — this
  argument never referenced which edges are clamped.
- The **Kirchhoff reduction** $u_r^{(1)}=-\partial w/\partial r$, $u_\theta^{(1)}=-\tfrac1r\partial w/\partial\theta$
  (Eq. (11)) and the resulting single-field PDE (Eq. (10)/(21)) — same reasoning.
- The **$r=R_i,R_o$ free-edge conditions** (Eq. (12)/(23): $\tau^{(1)}_{rr}=0$,
  $\tau^{(0)}_{rz}+\tfrac1r\partial\tau^{(1)}_{r\theta}/\partial\theta=0$) — these came from the
  *first* boundary integral of Eq. (6) (the one over $\theta\in[-\Theta,\Theta]$ at fixed
  $r=R_i,R_o$), which never touches $\theta=-\Theta$'s edge type either.
- **The Frobenius radial solution machinery** (Eqs. (24)–(41)): the dispersion relation
  connects the $\theta$-wavenumber $\xi$ to $\bar\Omega$ by requiring the *radial* eigenfunction to
  satisfy Eq. (23) exactly at $r=R_i,R_o$ — again, no dependence on the $\theta=-\Theta$ edge
  type. **Every already-validated OOP/IP branch-tracking and Frobenius-series machinery in
  `plate_solver` carries over to FFFF unchanged**, which matches how the project already treats
  it (Appendix A: "the SAME code path already produces the pinned cantilever AND free-free
  tables").

So the entire apparatus through Ref. [seok1] Eq. (41) is common to both boundary conditions.
Only §4 ("Variational approximation," their Eq. (42) onward) needs re-deriving.

## 4. The paper's own "what's left" statement — Eq. (42) [TRANSCRIBED]

Because the solution ansatz (Eq. (36)) satisfies the PDE and the $r=R_i,R_o$ conditions exactly
by construction, the paper states plainly (p. 767): *"all that remains in the variational
equation (6) is"* — i.e., Eq. (42) is *literally* Eq. (6) with the area integral and the
$r=R_i,R_o$ boundary integral dropped (they are now identically zero), leaving only the
$\theta=\pm\Theta$ terms and the corner terms:

$$
-\int_{R_i}^{R_o}\!\Big[\Big[\big(\tau^{(0)}_{\theta z}+\tfrac{\partial\tau^{(1)}_{\theta r}}{\partial r}\big)\delta w-\tau^{(1)}_{\theta\theta}\delta\big(\tfrac1r\tfrac{\partial w}{\partial\theta}\big)\Big]_{\theta=\Theta}
+\Big[w\,\delta\tau^{(0)}_{\theta z}-\tfrac{\partial w}{\partial r}\delta\tau^{(1)}_{r\theta}\Big]_{\theta=-\Theta}\Big]dr
$$
$$
-\big[\tau^{(1)}_{r\theta}\delta w\big]^{r=R_o,\theta=-\Theta}_{r=R_i,\theta=-\Theta}
+2\big[\tau^{(1)}_{r\theta}\delta w\big]^{r=R_o,\theta=\Theta}_{r=R_i,\theta=\Theta}=0.
\tag{42}
$$

This is exactly Eq. (6)'s boundary terms, unchanged — confirming Eq. (42) is the right level at
which to make the FFFF substitution: **replace only the bracketed content at $\theta=-\Theta$**,
leaving the PDE/r-edge machinery (§3) and the $\theta=\Theta$ term untouched.

## 5. Constructing the free-edge term at $\theta=-\Theta$ — [DERIVED, moderate confidence]

I derived the natural-edge term at $\theta=-\Theta$ directly from the general $c_N$ integral in
Eq. (1a), rather than just asserting it mirrors $\theta=\Theta$, by explicitly transforming the
indicial $n_a,s_b$ objects into $r,\theta$ components at that edge.

**Orientation.** Traversing the sector boundary counterclockwise starting at $(R_i,-\Theta)$: the
$\theta=-\Theta$ edge is traversed with *increasing* $r$ (tangent $\hat s=+\hat r$, outward normal
$\hat n=-\hat\theta$); the $\theta=\Theta$ edge is traversed with *decreasing* $r$ (tangent
$\hat s=-\hat r$, outward normal $\hat n=+\hat\theta$).

**At $\theta=\Theta$** ($n=+\hat\theta$, $s=-\hat r$): $n_a t^{(0)}_{a3}=\tau^{(0)}_{\theta z}$;
$u_{3,n}^{(0)}=\tfrac1r\partial w/\partial\theta$; $n_a t^{(1)}_{a\beta}s_\beta=-\tau^{(1)}_{\theta r}$,
so $(n_a t^{(1)}_{a\beta}s_\beta)_{,s}=\partial\tau^{(1)}_{\theta r}/\partial r$ (the $s$-derivative
along $-\hat r$ flips the sign of the $r$-derivative, canceling the leading minus);
$n_a t^{(1)}_{a\beta}n_\beta=\tau^{(1)}_{\theta\theta}$. Substituting into Eq. (1a)'s third term
and converting $ds\to dr$ (the tangent direction and the $r$-integration limits both reverse,
canceling) reproduces the paper's own printed $\theta=\Theta$ bracket exactly, with its leading
minus sign — a working self-consistency check on the method.

**At $\theta=-\Theta$** ($n=-\hat\theta$, $s=+\hat r$): $n_a t^{(0)}_{a3}=-\tau^{(0)}_{\theta z}$;
$u_{3,n}^{(0)}=-\tfrac1r\partial w/\partial\theta$; $n_a t^{(1)}_{a\beta}s_\beta=-\tau^{(1)}_{\theta r}$
(same value as at $\Theta$ — see the note below), so
$(n_a t^{(1)}_{a\beta}s_\beta)_{,s}=-\partial\tau^{(1)}_{\theta r}/\partial r$ (no sign flip this
time, since $s=+\hat r$ already matches increasing $r$); $n_a t^{(1)}_{a\beta}n_\beta=\tau^{(1)}_{\theta\theta}$
(two minus signs cancel). The integrand becomes
$$
-\Big(\tau^{(0)}_{\theta z}+\tfrac{\partial\tau^{(1)}_{\theta r}}{\partial r}\Big)\delta w
+\tau^{(1)}_{\theta\theta}\delta\Big(\tfrac1r\tfrac{\partial w}{\partial\theta}\Big)
=-\Big[\Big(\tau^{(0)}_{\theta z}+\tfrac{\partial\tau^{(1)}_{\theta r}}{\partial r}\Big)\delta w-\tau^{(1)}_{\theta\theta}\delta\Big(\tfrac1r\tfrac{\partial w}{\partial\theta}\Big)\Big],
$$
i.e. the **negative of the $\theta=\Theta$ bracket** — and with $ds=+dr$ (no limit-flip this
time), Eq. (1a)'s leading minus combines with this internal minus to give an **overall plus**:

$$
\boxed{+\int_{R_i}^{R_o}\Big[\big(\tau^{(0)}_{\theta z}+\tfrac{\partial\tau^{(1)}_{\theta r}}{\partial r}\big)\delta w-\tau^{(1)}_{\theta\theta}\delta\big(\tfrac1r\tfrac{\partial w}{\partial\theta}\big)\Big]_{\theta=-\Theta}dr}
\tag{FF-42a, DERIVED}
$$

**Same bracket content as $\theta=\Theta$, opposite overall sign.** This is not an assumption —
it falls out of the orientation reversal. **Confirmed against the deployed code**, not just
matched by structural analogy: `core_solvers.py::_build_K_real`'s integral assembly for a FREE
edge is `acc += _orient[ei]*(aj[k]*Ai[k] - bj[k]*Bi[k])`, where `(aj,Ai,bj,Bi)` (the projected
stress/displacement node vectors, lines 279-305) are built identically for every edge — the
*only* per-edge difference in that line is `_orient[ei] = e.sign`, which is `+1` at $\theta=\Theta$
and `-1` at $\theta=-\Theta$ for `FreeFreeOOP`. That is precisely "same bracket content, opposite
overall sign," read directly from the operative code rather than inferred. Retained as "moderate"
confidence for the hand-derivation part above (the six-substitution n,s bookkeeping is still worth
checking independently if anyone wants to), but the *conclusion* is now code-confirmed, not just
self-consistency-checked.

*Note on the identical $n_at^{(1)}_{a\beta}s_\beta$ value at both edges:* for a symmetric tensor,
$n_a\tau_{ab}s_b$ is invariant under simultaneously flipping both $n\to-n$ and $s\to-s$ — which is
exactly what happens between $\theta=\Theta$ and $\theta=-\Theta$. This is a useful independent
check (it must come out the same at both edges), and it does.

## 6. The FFFF analogue of Eq. (42) — [DERIVED, moderate confidence]

Substituting (FF-42a) for the old mixed term in Eq. (42):

$$
+\int_{R_i}^{R_o}\Big[\big(\tau^{(0)}_{\theta z}+\tfrac{\partial\tau^{(1)}_{\theta r}}{\partial r}\big)\delta w-\tau^{(1)}_{\theta\theta}\delta\big(\tfrac1r\tfrac{\partial w}{\partial\theta}\big)\Big]_{\theta=-\Theta}dr
-\int_{R_i}^{R_o}\Big[\big(\tau^{(0)}_{\theta z}+\tfrac{\partial\tau^{(1)}_{\theta r}}{\partial r}\big)\delta w-\tau^{(1)}_{\theta\theta}\delta\big(\tfrac1r\tfrac{\partial w}{\partial\theta}\big)\Big]_{\theta=\Theta}dr
$$
$$
+\;(\text{corner terms, see §7})=0.
\tag{FF-42}
$$

i.e. $\big[(\ldots)\delta w-(\ldots)\delta(\ldots)\big]_{-\Theta}^{\Theta}$ in the Eq. (7) bracket
sense — the identical natural-edge bracket, alternating in sign between the two edges.

## 7. The corner terms — [CONFIRMED AGAINST DEPLOYED CODE, corrects an earlier
draft error — AND NOW INDEPENDENTLY DERIVED, see §7a below]

**This section was wrong in the first draft of this document and has been corrected.** The
first pass guessed "$-2(T-\bar\nu)$ reused at both edges with the same sign," reasoning loosely
from Eq. (43)'s printed cantilever coefficients. Reading the actual deployed assembly
(`core_solvers.py::_build_K_real`, `boundary.py::_OOP_CORNER_COEFF`/`EdgeSpec`) directly shows
that guess was wrong, and gives the precise, correct answer instead of a "matched" placeholder.

**What the paper's corner terms encode.** Eq. (1a)'s fourth term,
$\int_{c_N}ds\,(n_a t^{(1)}_{a\beta}s_\beta\,\delta w)_{,s}$, is an exact differential along $c_N$;
integrating it produces boundary values at the *ends* of each continuous $c_N$ segment via the
fundamental theorem of calculus. In the cantilever case, $c_N$ is the three-free-edge path
$r=R_i \to \theta=\Theta \to r=R_o$ (reading the "wall" out); it has two open ends, both sitting
at the wall ($\theta=-\Theta$) — giving the paper's "$-1\times[\cdot]^{R_o}_{R_i}$" term at
$\theta=-\Theta$ ("arises from ... evaluated at the wall," p. 762). The *other* two corners
($\theta=\Theta$ meeting $r=R_i$ and $r=R_o$) are interior to the continuous $c_N$ path — the
path doesn't terminate there, it turns a corner — and the paper states explicitly that the
$\theta=\Theta$ coefficient of 2 "came from the two jump conditions across the edges of
discontinuity" (Eq. (8)): a jump term at *each* of the two corners, because $n_a,s_b$ are
discontinuous where the path turns 90°, even though the path itself doesn't stop.

**For FFFF, $\theta=-\Theta$ is released, so $c_N$ becomes the full closed boundary** — there is
no wall left to terminate at, and *all four* corners are now "path turns a corner" corners of
the same type as the $\theta=\Theta$ pair in the cantilever case. By direct analogy this predicts
a coefficient of *magnitude* 2 at *both* circumferential corner brackets for FFFF, not just
$\theta=\Theta$ — a genuine, non-trivial prediction, and (as confirmed below) the magnitude part
of this prediction is correct. It says nothing about the relative *sign* between the two edges,
which is a separate question resolved not by hand-deriving Eq. (8)'s jump algebra (still not
attempted, for the same reasons given in the first draft — genuine new work beyond what even the
source documents, with real risk of a silent error) but by reading the actual constants the
deployed, independently-FE-validated assembly uses.

**Reading `core_solvers.py::_build_K_real` and `boundary.py` directly** (not reconstructing them,
reading the literal operative code):

```python
# boundary.py
_OOP_CORNER_COEFF = {EdgeKind.CLAMPED: 1, EdgeKind.FREE: -2}
class FreeFreeOOP(BoundaryCondition):
    def __init__(self):
        super().__init__([EdgeSpec(+1, EdgeKind.FREE), EdgeSpec(-1, EdgeKind.FREE)])
        #                  theta=+Theta  FREE          theta=-Theta  FREE

# core_solvers.py _build_K_real
corner_coeff = [_OOP_CORNER_COEFF[e.kind] for e in edges]     # [-2, -2] for FreeFreeOOP
_orient = [(e.sign if edge_free[k] else 1) for k, e in enumerate(edges)]
#          FREE edges use e.sign; CLAMPED edges always +1 regardless of e.sign
# ... corner += _orient[ei] * corner_coeff[ei] * (jump)
```

For `FreeFreeOOP`, both edges have `corner_coeff = -2` (the *constant* really is the same at
both edges, as the source comment on line 25-26 of `boundary.py` says), **but** `_orient` is
`e.sign` for a FREE edge — `+1` at $\theta=\Theta$, `-1` at $\theta=-\Theta$ — so the net
per-edge contribution to $K$ is:

$$
\theta=\Theta:\ \ (+1)\times(-2)\times(\text{jump}) = -2\times(\text{jump}),\qquad
\theta=-\Theta:\ \ (-1)\times(-2)\times(\text{jump}) = +2\times(\text{jump}).
$$

**This corrects the first draft, which stated "$-2(T-\bar\nu)$ reused at both edges with the same
sign" — that was wrong.** The correct result is *sign-alternating*, exactly like the surface term
in §5, both driven by the same `e.sign` orientation factor:

$$
\boxed{+2(T-\bar\nu)\Big[\tfrac{\partial^2}{\partial\tilde r\partial\theta}\Big\{\tfrac{w}{\tilde r+\bar r_0}\Big\}\Big]^{\tilde r=\pi/2,\theta=-\Theta}_{\tilde r=-\pi/2,\theta=-\Theta}
-2(T-\bar\nu)\Big[\tfrac{\partial^2}{\partial\tilde r\partial\theta}\Big\{\tfrac{w}{\tilde r+\bar r_0}\Big\}\Big]^{\tilde r=\pi/2,\theta=\Theta}_{\tilde r=-\pi/2,\theta=\Theta}=0}
\tag{FF-43 corner terms, CONFIRMED AGAINST CODE}
$$

This is a nicer result than the first draft's guess, not just a correction: it means the *same*
single mechanism — the `e.sign` orientation factor from the $[\cdot]_{-\Theta}^{\Theta}$ structure
of the underlying functional — governs the sign alternation of *both* the surface term (§5) and
the corner term, rather than the two needing separate, unrelated sign rules. That unification is
exactly the kind of thing a correct derivation should produce, and its absence in the first draft
(where the surface term alternated but the corner term supposedly didn't) was itself a warning
sign I should have caught before writing it down — noted here so the same mistake isn't repeated.
The magnitude (2, same as the cantilever's free edge) is still supported by the closed-contour
argument above; the sign is read directly from the validated source, not re-derived from Eq. (8).

## 7a. Independent re-derivation of the corner-jump algebra from Eq. (1a)
directly — [DERIVED, high confidence] — closes the §9 follow-up gap

**2026-07-18 (second late-night pass).** §7's magnitude argument ("by direct
analogy") and its sign conclusion (read from `_orient`/`_OOP_CORNER_COEFF`)
are both independently reproduced here from Eq. (1a)'s own jump structure —
no code read required. This is the standalone-algebra-check §9 flagged as
"lower priority... a nice-to-have for extra certainty, not a blocking gap."
It is no longer just a nice-to-have result: it is now done, and it closes
that gap.

**Setup.** Eq. (1a)'s fourth term, $\int_{c_N}ds\,(n_a t^{(1)}_{a\beta}s_\beta\,\delta w)_{,s}$,
is an exact derivative along arc length $s$. Integrating a **closed** contour
$c_N$ segment-by-segment via the fundamental theorem of calculus does not
vanish, because $n_a,s_b$ are discontinuous at each corner even though the
physical stress field $\tau^{(1)}_{r\theta}$ is continuous there. Writing
$F\equiv n_a t^{(1)}_{a\beta}s_\beta=n_a\tau^{(1)}_{a\beta}s_\beta$ (a scalar
built from the *geometric* $n,s$ convention of whichever segment is being
evaluated), summing the per-segment FTC contributions and grouping by shared
corner point $P$ gives a jump contribution at every corner:
$$
\text{jump}(P) = F_{\text{incoming segment}}(P) - F_{\text{outgoing segment}}(P),
$$
each evaluated with its own segment's $(n,s)$ convention, at the same
physical point $P$ (where $\tau^{(1)}_{r\theta}$ takes one well-defined
value). This is the literal mechanism behind the paper's own phrase "the two
jump conditions across the edges of discontinuity" (Eq. (8)); the analysis
below just carries it out explicitly.

**The four segments and their $(n,s)$ conventions, FFFF (closed $c_N$),
traversed CCW.** The $\theta=\pm\Theta$ conventions are already fixed in §5
(cross-checked there against the paper's own printed $\theta=\Theta$ term);
the two arc segments follow the same CCW-with-interior-on-the-left rule:

| segment | tangent $s$ | normal $n$ (outward) | $F=n_a\tau^{(1)}_{a\beta}s_\beta$ |
|---|---|---|---|
| $r=R_o$ arc ($\theta:-\Theta\to\Theta$) | $+\hat\theta$ | $+\hat r$ | $+\tau^{(1)}_{r\theta}$ |
| $\theta=\Theta$ radial ($r:R_o\to R_i$) | $-\hat r$ | $+\hat\theta$ | $-\tau^{(1)}_{r\theta}$ |
| $r=R_i$ arc ($\theta:\Theta\to-\Theta$) | $-\hat\theta$ | $-\hat r$ | $+\tau^{(1)}_{r\theta}$ |
| $\theta=-\Theta$ radial ($r:R_i\to R_o$) | $+\hat r$ | $-\hat\theta$ | $-\tau^{(1)}_{r\theta}$ |

(Each $F$ uses $n_a\tau_{a\beta}s_\beta$ with $\tau^{(1)}$ symmetric in
$(r,\theta)$, so only the sign of the $\hat r$-$\hat\theta$ product matters;
e.g. $\theta=\Theta$: $(\hat\theta)_a\tau_{a\beta}(-\hat r)_\beta=-\tau_{\theta r}=-\tau_{r\theta}$.)

**Jump at each of the four corners**, walking CCW from $(R_i,-\Theta)$:

$$
\begin{aligned}
(R_o,-\Theta):&\ \ \text{incoming }[-\Theta\text{ radial}]-\text{outgoing }[R_o\text{ arc}] = (-\tau_{r\theta})-(\tau_{r\theta}) = -2\tau_{r\theta}\Big|_{(R_o,-\Theta)}\\
(R_o,\Theta):&\ \ \text{incoming }[R_o\text{ arc}]-\text{outgoing }[\Theta\text{ radial}] = (\tau_{r\theta})-(-\tau_{r\theta}) = +2\tau_{r\theta}\Big|_{(R_o,\Theta)}\\
(R_i,\Theta):&\ \ \text{incoming }[\Theta\text{ radial}]-\text{outgoing }[R_i\text{ arc}] = (-\tau_{r\theta})-(\tau_{r\theta}) = -2\tau_{r\theta}\Big|_{(R_i,\Theta)}\\
(R_i,-\Theta):&\ \ \text{incoming }[R_i\text{ arc}]-\text{outgoing }[-\Theta\text{ radial}] = (\tau_{r\theta})-(-\tau_{r\theta}) = +2\tau_{r\theta}\Big|_{(R_i,-\Theta)}
\end{aligned}
$$

**Self-consistency check (validates the method before trusting it on the new
corners): sum the two $\theta=\Theta$ corners**, the ones the cantilever
*already has* and whose combined coefficient is published:
$$
(R_o,\Theta)+(R_i,\Theta) = 2\tau_{r\theta}\big|_{R_o,\Theta} - 2\tau_{r\theta}\big|_{R_i,\Theta}
= 2\big[\tau^{(1)}_{r\theta}\delta w\big]^{R_o,\Theta}_{R_i,\Theta},
$$
**exactly the paper's own printed cantilever term**,
"$+2[\tau^{(1)}_{r\theta}\delta w]^{r=R_o,\theta=\Theta}_{r=R_i,\theta=\Theta}$" (Eq. (6)). This
reproduces a known, already-published result from the jump-bookkeeping method
alone, with no reference to the deployed code — good evidence the method
itself is being applied correctly before it is used to settle the unknown case.

**Now apply the identical, already-validated method to the two corners that
are new in FFFF** (in the cantilever these were "open ends" of the truncated
$c_N$ path at the wall, contributing the plain endpoint term
"$-[\tau^{(1)}_{r\theta}\delta w]^{R_o,-\Theta}_{R_i,-\Theta}$"; releasing
$\theta=-\Theta$ makes $c_N$ close through them, so they become smooth-turn
jump corners of the *same type* as the $\Theta$ pair, not open ends):
$$
(R_o,-\Theta)+(R_i,-\Theta) = -2\tau_{r\theta}\big|_{R_o,-\Theta} + 2\tau_{r\theta}\big|_{R_i,-\Theta}
= -2\big[\tau^{(1)}_{r\theta}\delta w\big]^{R_o,-\Theta}_{R_i,-\Theta}.
$$

**Conclusion.** The total FFFF corner contribution is
$$
\boxed{-2\big[\tau^{(1)}_{r\theta}\delta w\big]^{R_o,-\Theta}_{R_i,-\Theta}
+2\big[\tau^{(1)}_{r\theta}\delta w\big]^{R_o,\Theta}_{R_i,\Theta} = 0,}
$$
i.e. **magnitude 2 at every corner** (matching $|{-}2|$ in `_OOP_CORNER_COEFF
[FREE]`), with the sign **alternating between the $\theta=-\Theta$ pair and
the $\theta=\Theta$ pair** — independently reproducing, from Eq. (1a)'s own
corner-jump structure alone, both halves of what §7 previously established
only by reading `_orient`/`corner_coeff` in the deployed code. This is now a
second, independent derivation path to the same conclusion, not a restatement
of the code read.

**One honest caveat, stated precisely so it isn't overclaimed:** which literal
edge is labeled "$+2$" versus "$-2$" is a matter of the overall equation's
sign convention (multiplying the whole stationarity condition by $-1$ swaps
both labels and changes nothing physical) — this derivation was not combined
with an independent re-check of the surface term's (§5) overall sign in the
*same* single equation, so it does not, by itself, fix which absolute label
is "correct" in isolation. What it *does* establish independent of any
convention choice, and independent of reading the code, is the two things
that actually matter physically: **the magnitude is 2 at every corner, and
the sign strictly alternates between the two edge-pairs** — precisely the
structure the deployed `e.sign`-driven `_orient` mechanism implements. That
this derivation's corner terms (§7a) and the surface term's derivation (§5)
both alternate via the *same* geometric mechanism (a $(n,s)$-orientation
reversal between $\theta=\Theta$ and $\theta=-\Theta$), and both independently
agree with the single unified `e.sign` factor actually coded, is itself
strong evidence the whole picture is self-consistent, not two unrelated
sign rules that happen to be coded the same way by coincidence.

## 8. The FFFF analogue of Eq. (43) — [DERIVED, high confidence] — the constitutive/dimensionless
substitution done directly at $\theta=-\Theta$, not reused by sign-flip analogy

**2026-07-19 rewrite.** The previous version of this section took a shortcut: it asserted that
substituting Eqs. (9c), (13)–(15), (20) into (FF-42a) is "algebraically identical" to the paper's
own $\theta=\Theta$ substitution and then simply re-applied the paper's *published, already-worked*
$\theta=\Theta$ block at $\theta=-\Theta$ with a sign flip, without carrying out the substitution
algebra itself. That shortcut's *conclusion* turns out to be correct (confirmed below), but stating
it without doing the algebra left one gap explicitly flagged as open in an earlier draft of this
document's own scoping notes: "the constitutive/dimensionless substitution... is mechanical but was not attempted." This
section attempts it, in full, at both edges — first as a self-consistency check against the paper's
own printed $\theta=\Theta$ result, then fresh at $\theta=-\Theta$.

**Scoping note on which constitutive relations are actually needed here.** The paper's Eq. (42) has
two structurally different circumferential boundary terms: a *natural* term (dual stress quantities
multiplying $\delta w$/$\delta(\partial_\theta w)$, at $\theta=\Theta$ in the cantilever) and a
*mixed/Lagrangian* term (primal $w$/$\partial_r w$ multiplying $\delta\tau^{(0)}_{\theta z}$/
$\delta\tau^{(1)}_{r\theta}$, at $\theta=-\Theta$ in the cantilever, from enforcing the clamped
condition). The paper's own printed Eq. (43) substitutes **Eq. (20) + Eq. (14)** into the natural
term and **Eq. (9c) + Eq. (14) + Eq. (15)** into the mixed term — hence the combined "(9c),
(13)–(15), (20)" list. For FFFF, $\theta=-\Theta$ is released, so its mixed/Lagrangian term is
*replaced* by a natural term (FF-42a, §5) rather than substituted in its old mixed form — meaning
the new piece needed for FFFF uses **only Eq. (20) + Eq. (14)**, exactly the same two constitutive
relations already used at $\theta=\Theta$, evaluated at the other edge. Eq. (9c) and Eq. (15) remain
correctly cited as part of the *paper's own* full Eq.(42)$\to$(43) substitution (both blocks
together), but they do not separately re-enter the FFFF-specific natural-term derivation below —
noted here explicitly so the constitutive-relation bookkeeping is honest, not just copied from the
task description.

### 8.1 The constitutive relations, transcribed precisely — [TRANSCRIBED]

Re-read directly from the page images (`seok1` pp. 763–764), confirming signs that OCR alone
garbles:
$$
\tau^{(1)}_{\theta\theta}=-\hat D\Big\{\hat\nu\frac{\partial^2w}{\partial r^2}+\frac{R}{r}\Big(\frac{\partial w}{\partial r}+\frac1r\frac{\partial^2w}{\partial\theta^2}\Big)\Big\}\tag{14}
$$
$$
\tau^{(1)}_{r\theta}=-\hat D(T-\hat\nu)\frac{\partial^2}{\partial r\partial\theta}\Big(\frac{w}{r}\Big)\tag{15}
$$
$$
\tau^{(0)}_{\theta z}+\frac{\partial\tau^{(1)}_{\theta r}}{\partial r}=-\hat D\Big\{\frac{R}{r^3}\frac{\partial^3w}{\partial\theta^3}-\frac{(2T-R)}{r^2}\frac{\partial^2w}{\partial r\partial\theta}+\frac{2T}{r}\Big(\frac{\partial^3w}{\partial r^2\partial\theta}+\frac1{r^2}\frac{\partial w}{\partial\theta}\Big)-\hat\nu\Big(\frac1r\frac{\partial^3w}{\partial r^2\partial\theta}+\frac2{r^3}\frac{\partial w}{\partial\theta}-\frac2{r^2}\frac{\partial^2w}{\partial r\partial\theta}\Big)\Big\}\tag{20}
$$
and the nondimensionalization (Eqs. (26)–(28)): $\bar r=\pi r/(2b)$, $\bar r_0=\pi r_0/(2b)$,
$\tilde r=\bar r-\bar r_0$, so $r=(2b/\pi)(\tilde r+\bar r_0)$, giving the two substitution rules
used throughout: $\partial/\partial r=(\pi/2b)\,\partial/\partial\tilde r$ and
$1/r=(\pi/2b)/(\tilde r+\bar r_0)$. $T,R,\hat\nu$ (Eq. (16)) are pure material-constant ratios,
unaffected by this substitution. $\theta$ itself is not rescaled.

### 8.2 Self-consistency check: reproducing the paper's own published $\theta=\Theta$ block of Eq. (43) — [DERIVED, high confidence]

**Every power of $1/r$ in Eq. (20) is third-order** (the whole expression has the dimensions of
$\hat D\times\partial^3/\partial(\text{length})^3$), so substituting $r\to(2b/\pi)(\tilde
r+\bar r_0)$ pulls out a *uniform* common factor $(\pi/2b)^3$ across every term:
$$
\tau^{(0)}_{\theta z}+\frac{\partial\tau^{(1)}_{\theta r}}{\partial r}=-\hat D\Big(\frac\pi{2b}\Big)^3\Big\{\frac{R}{(\tilde r+\bar r_0)^3}\frac{\partial^3w}{\partial\theta^3}-\frac{(2T-R)}{(\tilde r+\bar r_0)^2}\frac{\partial^2w}{\partial\tilde r\partial\theta}+\frac{2T}{(\tilde r+\bar r_0)}\Big(\frac{\partial^3w}{\partial\tilde r^2\partial\theta}+\frac1{(\tilde r+\bar r_0)^2}\frac{\partial w}{\partial\theta}\Big)
$$
$$
-\hat\nu\Big(\frac1{(\tilde r+\bar r_0)}\frac{\partial^3w}{\partial\tilde r^2\partial\theta}+\frac2{(\tilde r+\bar r_0)^3}\frac{\partial w}{\partial\theta}-\frac2{(\tilde r+\bar r_0)^2}\frac{\partial^2w}{\partial\tilde r\partial\theta}\Big)\Big\}.
$$
Factoring $\tfrac1{(\tilde r+\bar r_0)}\tfrac{\partial}{\partial\theta}\{\cdot\}$ out of the
non-$\hat\nu$ terms reproduces
$\tfrac1{(\tilde r+\bar r_0)}\partial_\theta\big\{\tfrac1{(\tilde r+\bar r_0)^2}(2Tw+R\partial_\theta^2w)+2T\partial_{\tilde r}^2w-\tfrac{(2T-R)}{(\tilde r+\bar r_0)}\partial_{\tilde r}w\big\}$
term for term (differentiate that bracket by $\theta$ and multiply by $1/(\tilde r+\bar r_0)$ to
check). The $\hat\nu$ group similarly factors: applying $\partial_\theta$ to
$\hat\nu\{2w/(\tilde r+\bar r_0)^2+(\tilde r+\bar r_0)^2\partial_{\tilde r}(\tfrac1{(\tilde
r+\bar r_0)^2}\partial_{\tilde r}w)\}$ and dividing by $(\tilde r+\bar r_0)$ gives exactly the
three $\hat\nu$-terms above (verified by direct expansion: $\partial_{\tilde
r}(\tfrac1{(\tilde r+\bar r_0)^2}\partial^2_{\tilde r\theta}w)=-\tfrac2{(\tilde
r+\bar r_0)^3}\partial^2_{\tilde r\theta}w+\tfrac1{(\tilde r+\bar r_0)^2}\partial^3_{\tilde
r\theta}w$, multiply by $(\tilde r+\bar r_0)^2/(\tilde r+\bar r_0)$). So
$$
\tau^{(0)}_{\theta z}+\frac{\partial\tau^{(1)}_{\theta r}}{\partial r}=-\hat D\Big(\frac\pi{2b}\Big)^3\frac1{(\tilde r+\bar r_0)}\frac{\partial}{\partial\theta}\Big[\Big\{\frac1{(\tilde r+\bar r_0)^2}(2Tw+R\partial_\theta^2w)+2T\partial_{\tilde r}^2w-\frac{(2T-R)}{(\tilde r+\bar r_0)}\partial_{\tilde r}w\Big\}
$$
$$
-\hat\nu\Big\{\frac{2w}{(\tilde r+\bar r_0)^2}+(\tilde r+\bar r_0)^2\partial_{\tilde r}\Big(\frac1{(\tilde r+\bar r_0)^2}\partial_{\tilde r}w\Big)\Big\}\Big].
$$
Likewise, Eq. (14) gives $\tau^{(1)}_{\theta\theta}=-\hat D(\pi/2b)^2\{\hat\nu\partial_{\tilde
r}^2w+\tfrac{R}{(\tilde r+\bar r_0)}\partial_{\tilde r}w+\tfrac{R}{(\tilde r+\bar r_0)^2}\partial_\theta^2w\}$,
and $\delta(\tfrac1r\partial_\theta w)=(\pi/2b)/(\tilde r+\bar r_0)\,\delta(\partial_\theta w)$, so
$$
-\tau^{(1)}_{\theta\theta}\,\delta\Big(\frac1r\frac{\partial w}{\partial\theta}\Big)=\hat D\Big(\frac\pi{2b}\Big)^3\frac1{(\tilde r+\bar r_0)}\Big\{\hat\nu\partial_{\tilde r}^2w+\frac{R}{(\tilde r+\bar r_0)}\partial_{\tilde r}w+\frac{R}{(\tilde r+\bar r_0)^2}\partial_\theta^2w\Big\}\delta(\partial_\theta w).
$$
Both pieces carry the identical prefactor $-\hat D(\pi/2b)^3$ (with the sign of the second flipped
back to $+\hat D(\pi/2b)^3\cdot(+1)$ by its own leading minus, matching the first). Eq. (42)'s
$\theta=\Theta$ bracket is multiplied by an outer $-\int(\cdots)\,dr$ with $dr=(2b/\pi)\,d\tilde r$,
so the combined constant is $-(-\hat D)(\pi/2b)^3(2b/\pi)=+\hat D(\pi/2b)^2$ for the $\delta w$ term
and $-(+\hat D)(\pi/2b)^3(2b/\pi)=-\hat D(\pi/2b)^2$ for the $\delta(\partial_\theta w)$ term.
Dividing the whole equation by the common positive constant $\hat D(\pi/2b)^2$ (legitimate — Eq.
(42) is homogeneous, $=0$) leaves **exactly**:
$$
\boxed{\int_{-\pi/2}^{\pi/2}\!d\tilde r\Big[\frac1{(\tilde r+\bar r_0)}\frac{\partial}{\partial\theta}\Big\{\frac1{(\tilde r+\bar r_0)^2}(2Tw+R\partial_\theta^2w)+2T\partial_{\tilde r}^2w-\frac{(2T-R)}{(\tilde r+\bar r_0)}\partial_{\tilde r}w-\hat\nu\Big\{\frac{2w}{(\tilde r+\bar r_0)^2}+(\tilde r+\bar r_0)^2\partial_{\tilde r}\Big(\frac1{(\tilde r+\bar r_0)^2}\partial_{\tilde r}w\Big)\Big\}\Big\}\delta w}
$$
$$
\boxed{-\Big\{\frac{R}{(\tilde r+\bar r_0)}\Big(\frac1{(\tilde r+\bar r_0)^2}\partial_\theta^2w+\frac1{(\tilde r+\bar r_0)}\partial_{\tilde r}w\Big)+\frac{\hat\nu}{(\tilde r+\bar r_0)}\partial_{\tilde r}^2w\Big\}\delta(\partial_\theta w)\Big]_{\theta=\Theta}}
$$
**— a term-for-term, sign-for-sign match to the paper's own published $\theta=\Theta$ block of Eq.
(43)**, including the ($\hat\nu\{\cdots\}$-inside-$\partial_\theta\{\cdots\}$) bracket nesting, which
a first read of the small page image can easily misparse as a separate, un-$\theta$-differentiated
term — confirmed by re-examining a magnified crop of the source page specifically to check the
bracket nesting before trusting the match. **This is the self-consistency check**: Eq. (20), Eq.
(14), Eq. (26)–(28), and Eq. (42)'s own leading sign, substituted from scratch with no reference to
the printed Eq. (43), reproduce Eq. (43)'s $\theta=\Theta$ block exactly.

### 8.3 The fresh substitution at $\theta=-\Theta$ — [DERIVED, high confidence]

Eqs. (20), (14), (26)–(28) contain no reference to $\Theta$ or $-\Theta$ — they are pointwise
relations holding at *any* $(r,\theta)$. So the identical algebra of §8.2 applies verbatim with
$\theta=\Theta$ replaced by $\theta=-\Theta$ throughout the bracket content. The **only** thing that
changes is the outer sign: §8.2 used Eq. (42)'s own $-\int(\cdots)_{\theta=\Theta}\,dr$, whereas
FF-42a (§5) established the $\theta=-\Theta$ natural term enters with the *opposite* outer sign,
$+\int(\cdots)_{\theta=-\Theta}\,dr$. Repeating §8.2's bookkeeping with that one sign changed flips
every term in the boxed result:
$$
\boxed{-\int_{-\pi/2}^{\pi/2}\!d\tilde r\Big[\frac1{(\tilde r+\bar r_0)}\frac{\partial}{\partial\theta}\Big\{\frac1{(\tilde r+\bar r_0)^2}(2Tw+R\partial_\theta^2w)+2T\partial_{\tilde r}^2w-\frac{(2T-R)}{(\tilde r+\bar r_0)}\partial_{\tilde r}w-\hat\nu\Big\{\frac{2w}{(\tilde r+\bar r_0)^2}+(\tilde r+\bar r_0)^2\partial_{\tilde r}\Big(\frac1{(\tilde r+\bar r_0)^2}\partial_{\tilde r}w\Big)\Big\}\Big\}\delta w}
$$
$$
\boxed{+\Big\{\frac{R}{(\tilde r+\bar r_0)}\Big(\frac1{(\tilde r+\bar r_0)^2}\partial_\theta^2w+\frac1{(\tilde r+\bar r_0)}\partial_{\tilde r}w\Big)+\frac{\hat\nu}{(\tilde r+\bar r_0)}\partial_{\tilde r}^2w\Big\}\delta(\partial_\theta w)\Big]_{\theta=-\Theta}}
\tag{FF-43, $\theta=-\Theta$ block, DERIVED}
$$
This reproduces, by direct construction rather than by pattern-matching, exactly what the previous
draft of this section asserted by re-applying the published $\theta=\Theta$ block with a sign flip
— the shortcut's conclusion was correct, and is now independently confirmed rather than merely
plausible.

### 8.4 The complete, fully self-contained FFFF analogue of Eq. (43)

Combining §8.2 (unchanged $\theta=\Theta$ block, now independently re-derived rather than only
quoted), §8.3 (new $\theta=-\Theta$ block), and §7a's corner result (magnitude 2, alternating sign):

$$
\int_{-\pi/2}^{\pi/2}\!d\tilde r\Big[\frac1{(\tilde r+\bar r_0)}\frac{\partial}{\partial\theta}\Big\{\frac1{(\tilde r+\bar r_0)^2}(2Tw+R\partial_\theta^2w)+2T\partial_{\tilde r}^2w-\frac{(2T-R)}{(\tilde r+\bar r_0)}\partial_{\tilde r}w-\hat\nu\Big\{\frac{2w}{(\tilde r+\bar r_0)^2}+(\tilde r+\bar r_0)^2\partial_{\tilde r}\Big(\frac1{(\tilde r+\bar r_0)^2}\partial_{\tilde r}w\Big)\Big\}\Big\}\delta w
$$
$$
-\Big\{\frac{R}{(\tilde r+\bar r_0)}\Big(\frac1{(\tilde r+\bar r_0)^2}\partial_\theta^2w+\frac1{(\tilde r+\bar r_0)}\partial_{\tilde r}w\Big)+\frac{\hat\nu}{(\tilde r+\bar r_0)}\partial_{\tilde r}^2w\Big\}\delta(\partial_\theta w)\Big]_{\theta=-\Theta}^{\theta=\Theta}
$$
$$
+2(T-\hat\nu)\Big[\frac{\partial^2}{\partial\tilde r\partial\theta}\Big\{\frac{w}{\tilde r+\bar r_0}\Big\}\Big]_{\tilde r=-\pi/2,\theta=-\Theta}^{\tilde r=\pi/2,\theta=-\Theta}-2(T-\hat\nu)\Big[\frac{\partial^2}{\partial\tilde r\partial\theta}\Big\{\frac{w}{\tilde r+\bar r_0}\Big\}\Big]_{\tilde r=-\pi/2,\theta=\Theta}^{\tilde r=\pi/2,\theta=\Theta}=0,
\tag{FF-43, complete}
$$
using the Eq.~(7)-style bracket notation $[\,\cdot\,]_{-\Theta}^{\Theta}$ for "value at $\Theta$
minus value at $-\Theta$." Every term is now written out explicitly in terms of $w$ and its
$\tilde r,\theta$ derivatives, substituted directly from Eqs. (14), (20), (26)–(28) at both edges
(not reused by analogy at either edge — §8.2 re-derives the $\theta=\Theta$ block as well as the
new $\theta=-\Theta$ one), plus the independently-derived §7a corner terms. This is the requested
fully self-contained FFFF analogue of Eq. (43).

**What this changes relative to the previous draft:** nothing in the final coefficients or signs —
the prior sign-flip-by-analogy shortcut happened to get the right answer. What changes is the
evidentiary basis: every term in both circumferential blocks is now reached by an explicit,
checkable substitution of Eqs. (14)/(20)/(26)–(28), self-consistency-checked against the paper's
own published $\theta=\Theta$ result first, rather than asserted to follow "by analogy" from
noting the substitution is mechanical.

## 9. Honesty summary — what this document does and does not establish

**2026-07-18 update:** §7a below closes the one gap this section used to
flag as open (the corner-jump algebra's independent re-derivation from
Eq. (1a), without reading the deployed code). The summary below is updated
to reflect that; nothing else in this document changed.

**2026-07-19 update: the other gap this section used to flag — the
constitutive/dimensionless substitution (Eq. (42)$\to$Eq. (43) analogue at
$\theta=-\Theta$) being reused by sign-flip analogy rather than carried out —
is now also closed.** §8 was rewritten: §8.1–8.2 substitute Eqs. (14), (20),
(26)–(28) into the natural-edge bracket from scratch and reproduce the
paper's own published $\theta=\Theta$ block of Eq. (43) term-for-term as a
self-consistency check (no reference to the printed Eq. (43) used in the
derivation itself); §8.3 repeats the identical algebra at $\theta=-\Theta$
with the one sign difference §5 already established; §8.4 assembles the
complete FF-43 with every term written out explicitly. The two open items
this document originally flagged (the corner-contour multiplicity
re-derivation, and the constitutive/dimensionless substitution) are both
now closed.

**2026-07-19, later: the "derive the entire functional from scratch"
ambition is reassessed and downgraded from required future work to not
needed.** On review with the user, that item referred to something
different from what §7a/§8 do: an independent line-by-line
integration-by-parts derivation of Eq. (1a)/(6) itself, starting from the
strong-form governing PDE (Ref. [seok1]'s Eq. (7)/classical $D\nabla^4w=
\rho H\omega^2w$). That is Ref. [seok1]'s own foundational derivation —
already published and peer-reviewed, textbook calculus-of-variations
technique (the same integration-by-parts bookkeeping that produces the
classical free-corner reaction $2M_{r\theta}$), and identical for the
cantilever and free-free cases alike, since Eq. (1a)/(6) is BC-independent
groundwork common to both. Redoing it would not add any evidence specific
to the free-free extension — it would only re-prove machinery this project
already trusts and uses unmodified (§3: "every already-validated OOP/IP
branch-tracking and Frobenius-series machinery... carries over to FFFF
unchanged"). What actually needed independent derivation — how Eq. (1a)/(6)
specializes once $\theta=-\Theta$ is released rather than clamped — is the
free-free-specific delta §5/§7/§7a/§8 supply, and that is done, checked two
independent ways for the corner term and self-consistency-checked for the
surface term. This document's §0/§10 should be read with this in mind: the
"(c) an explicit construction of the free-free replacement term" goal §0
states has been fully met; the open-ended "derive the entire functional"
framing that appeared in later status notes overstated what was still
missing.

**Solidly established (would need a real error to be wrong):**
- Eq. (6) is asymmetric between $\theta=\pm\Theta$ *before* any BC-specific substitution, and the
  asymmetry is the mixed/$c_C$ vs. natural/$c_N$ distinction, traced to the "three free edges"
  text and the corner-jump description.
- Nothing in the PDE, Kirchhoff reduction, or $r=R_i,R_o$ conditions changes for FFFF (§3) — the
  entire existing Frobenius/branch-tracking machinery is reused unchanged, matching how the
  project already treats it.
- Eq. (42) is literally Eq. (6)'s boundary terms with the area/$r$-edge integrals dropped — the
  correct level to make the FFFF substitution.
- The released $\theta=-\Theta$ surface term has the *same bracket content* as $\theta=\Theta$'s
  (§5), by direct construction from Eq. (1a), not by assumed symmetry — and the accompanying
  sign-alternation conclusion is now also confirmed directly against `_build_K_real`'s `_orient`
  handling, not just self-consistency-checked by the hand derivation alone.
- **The corner term's magnitude (2, every corner) AND its sign-alternation between the
  $\theta=\Theta$ and $\theta=-\Theta$ edge-pairs (§7a)** — now derived two independent ways: (i)
  directly from Eq. (1a)'s own corner-jump bookkeeping, applied uniformly to all four corners of
  the now-closed contour and validated by first reproducing the cantilever's own already-published
  $\theta=\Theta$ term before trusting the method on the new $\theta=-\Theta$ corners, and (ii) by
  reading `_OOP_CORNER_COEFF`/`EdgeSpec.sign`/`_orient` directly in `boundary.py`/
  `core_solvers.py::_build_K_real` (§7). The two agree. One residual, precisely-scoped caveat:
  §7a's derivation alone doesn't fix which literal edge gets the "+" label in isolation (that's an
  overall-equation sign convention, not physical content) — but the magnitude and the alternation,
  which are the physically meaningful claims, are now independently derived, not just code-matched.
- **The full FFFF analogue of Eq. (43) (§8)** — every term at both circumferential edges is now
  reached by an explicit substitution of Eqs. (14), (20), (26)–(28), not reused from the published
  $\theta=\Theta$ block by sign-flip analogy. Self-consistency-checked first (§8.2 reproduces the
  paper's own printed $\theta=\Theta$ block term-for-term, with no reference to that printed result
  used in the derivation), then applied fresh at $\theta=-\Theta$ (§8.3) using only the sign
  established independently in §5. The final coefficients are unchanged from what the earlier
  sign-flip shortcut asserted — this closes the evidentiary gap, not a numerical correction.

**Remaining gap, now genuinely small and explicitly scoped:**
- §7a's derivation was not additionally combined with an independent re-verification of the
  surface term's (§5) overall sign *in the same single combined equation* — doing so would fix
  the absolute sign label (which edge is "+2" vs "-2") from first principles alone, with no
  appeal to the code at all. This is a real but low-stakes, precisely bounded gap: both the
  magnitude and the alternation pattern (the only physically meaningful content) are already
  established two independent ways; only the arbitrary overall-sign-convention label is
  code-matched rather than re-derived. Lower priority than the pre-§7a gap this replaces, since
  what's missing no longer has any bearing on whether the deployed matrix assembly is correct.

**Confirmed by reading the deployed source directly (corroborates §7a, not a substitute for it):**
- `_OOP_CORNER_COEFF`/`EdgeSpec.sign`/`_orient` in `boundary.py` and
  `core_solvers.py::_build_K_real`, re-checked again in this pass and unchanged since the
  2026-07-18 sign-fix session. The first draft of this document guessed "same sign at both
  edges" without checking the code (or deriving it from Eq. (1a)) and was wrong — corrected in
  §7, with the wrong version left visible there as a record of the mistake rather than a silently
  fixed one.

## 10. Suggested use in the paper

This is more derivation than Appendix B currently needs to carry in full — a condensed version
(§§2–6 essentials, one paragraph) is the right amount for Appendix B itself, with this full
document cited as supplementary material, plus one sentence noting the corner-term sign was
confirmed by direct inspection of the deployed `_OOP_CORNER_COEFF`/`EdgeSpec` source rather than
independently re-derived from Eq. (8)'s raw jump algebra (consistent with the paper's existing,
and good, practice of flagging exactly this kind of thing elsewhere).

**2026-07-18 update: folded into `PAPER1_FREEFREE_DRAFT.tex` Appendix B**, replacing the earlier
"not reproduced here" disclaimer. See that file's Appendix B for the condensed version.

**2026-07-18 update, second pass (same night): the sentence above is now partly superseded by
§7a.** The corner-term's magnitude and sign-alternation ARE now independently re-derived from
Eq. (1a)'s jump algebra (§7a), not solely code-matched — Appendix B's condensed text should be
updated to say so (a one-clause change, not a rewrite), and this document's own supplementary-
material role now carries a strictly stronger result than when it was first folded in.

**2026-07-19 update: §8's constitutive/dimensionless substitution is now also done in full,**
closing the second of the two gaps this document's own scoping notes originally listed as open.
Appendix B's paragraph beginning "Substituting Eqs.~(9c), (13)--(15) and (20) into Eq.~(B.3) is
algebraically identical..." (previously stating the substitution was *not* redone from scratch)
should be updated to reflect that it now has been, with the self-consistency check summarized in
one or two sentences and this document cited for the full working (§8.1–8.4).
