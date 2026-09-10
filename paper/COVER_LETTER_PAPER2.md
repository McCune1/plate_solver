# Cover letter — Paper 2 (rectangular FFFF)

Addressed to the JSV Editor-in-Chief. Keep to one page. Submit together
with Paper 1 (annular FFFF); both share https://github.com/McCune1/plate_solver.
Paper 1 matches tag `v1.0.0` (equivalently `paper1-jsv`); this
manuscript matches `v1.1.0` (equivalently `paper2-jsv`). Do not invent
a second Zenodo DOI story.

---

**Subject:** Submission of "Exact Wave-Function Solution of Completely
Free Rectangular Plates: Four-Corner Kirchhoff Jump and Mode Screening"
(companion to the annular-sector manuscript submitted simultaneously)

Dear Editors,

We submit the enclosed manuscript for consideration as a full-length
article in the *Journal of Sound and Vibration*, with its Supplementary
Material. It is the rectangular companion of a manuscript submitted
simultaneously to this journal, *Exact Wave-Function Solution and
Spurious-Mode Discrimination for Completely Free Annular Sector Plates*.
The two papers share a solver and a public repository; they are not
the same analysis on a second geometry.

**What it does.** Seok, Tiersten and Scarton (*JSV* **271** (2004)
131–146 and 147–158) solved the rectangular plate by an exact
wave-function boundary determinant, but only for the clamped-free
configuration. This paper extends that method to the completely free
rectangle in both flexural and extensional motion. For flexure the
Kirchhoff twisting-moment jump at a free–free corner is not optional:
a naive four-corner sum of the cantilever-style terms vanishes by a
parity identity, and the physical jump is a checkerboard combination
obtained from their Eq. (16) on a closed contour. A persistence screen
(Screen B), frozen before extra aspect ratios were run, separates
modes that survive a modest complex-pair enlargement from a small
family of known artifacts.

**Why it is not a transplant of the annular companion.** The annular
paper screens a two-edge Galerkin projection by pointwise residuals
and by SUBDOM. After the long edges of the rectangle are satisfied
exactly, that residual taxonomy has nothing to act on. The remaining
defect is a mixture of basis-limited dips and a few aspect-ratio-
independent leak families, which is why this paper uses persistence
rather than a residual cut. Completely free in-plane motion has no
Kirchhoff jump at all: the assembler is an edge-list swap.

**Evidence.** Independent half-model SHELL281 calculations at five
aspect ratios, mesh-converged to <0.02%, confirm the screened flexural
candidates (typically inside 1%). Fifteen of Leissa's (1973) F–F–F–F
Ritz values convert onto those matches within 1.10%. An identity-
weighted MAC against the same half-model eigenvectors confirms 36 of
the 56 published flexural rows (every first symmetric mode
MAC ≥ 0.959); known Screen B leaks stay below the bar except where a
leak frequency coincides with a physical mode. The MAC bar 0.744 is
the annular companion's production call, applied unre-tuned as a
conservative shape-confirmation threshold, not recalibrated here.
In-plane, half-model PLANE183 matches every non-rigid target in
Ω̄ ∈ [0.02, 2.50] (92 unique frequencies, 16 recovered at an enlarged
basis in a window centred on the finite-element value). No in-plane
eigenvector cross-check is presented: the in-plane result is a
frequency census. Bardell, Langley and Dunsdon's 1996 first six
F–F–F–F frequencies at a/b = 1 and 2 agree to at most 1.69%
(independent FE versus Bardell at most 0.019%); nineteen of Gorman's
(2004) superposition eigenvalues in the same window agree to at most
1.72% (independent FE versus Gorman at most 0.11%). Default
clamped-free assembly is unchanged (`SOLVER_VERSION` 2026-07-10.s10).

**What we report against ourselves.** Production-basis flexural
discovery is a floor, not a census: two genuine symmetric modes appear
only on the persist basis; one square-plate symmetric target
(Λ_FE = 1.982) is a named unmatched curiosity, not a missing mode.
Several higher-Λ reconstructions do not MAC-confirm; those frequency
rows are kept and the shapes are not claimed. Two Table 2 rows fail
Screen B at the continuum optimum and are retained on independent FE
agreement. Screen B is not a portable residual screen in the annular
sense.

The solver, finite-element decks, validated tables and a from-scratch
regression suite are released under an MIT license as a public GitHub
repository, https://github.com/McCune1/plate_solver. The manuscript is
original, is not under consideration elsewhere except as the companion
named above, and all authors have approved this submission.

Yours sincerely,
Garrek McCune (corresponding author), Daniel Stutts
Mechanical and Aerospace Engineering, Missouri University of Science and
Technology
