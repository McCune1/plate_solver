# Cover letter (draft)

Addressed to the JSV Editor-in-Chief. Keep to one page. The GitHub
repository is public (see the last paragraph); insert the Zenodo DOI there
too once that archive is minted.

---

**Subject:** Submission of "Exact Wave-Function Solution and Spurious-Mode
Discrimination for Completely Free Annular Sector Plates"

Dear Editors,

We submit the enclosed manuscript for consideration as a full-length article
in the *Journal of Sound and Vibration*, with its Supplementary Material.

**What it does.** Seok and Tiersten (*JSV* **271** (2004) 757–772 and
773–787) solved the annular sector plate by an exact wave-function boundary
determinant, but only for the clamped-free configuration, and never released
the numerical machinery. This paper extends that method to the completely
free sector — the hardest of the classical boundary conditions for it,
because releasing both circumferential edges leaves every condition natural —
and reports free-free spectra for both flexural and extensional motion.
The extensional tables have independent published semi-analytical
anchors; the primary flexural tables are FE-validated, with a first
literature check at a third radius ratio.

**Why it is not routine.** Weak enforcement on both edges makes the
determinant acquire zeros that are not natural frequencies. We identify this
failure mode, show it has no counterpart in the cantilever case, and develop
instruments that separate catalogued spurious zeros from physical modes.
The residual screen is robust in flexure and is only one-sided in
extension. A two-sided extensional bar exists separately, on identity-weighted
mode-shape correlation (production threshold 0.744), not on the residual.
One spurious zero has a deeper determinant minimum than any of the eight
physical flexural modes; ranked on that alone it would be the strongest
apparent mode in the spectrum.

**Evidence.** The implementation first reproduces all 29 published
Seok–Tiersten cantilever frequencies from scratch. A fully-default
detector pass is not complete (4 of 26 in-plane modes missed, later
recovered); the published tables use the §5 instruments. Those tables
are matched against mesh-converged finite elements (8 of 8 isotropic
flexural to 0.64%; 8 of 8 orthotropic flexural to 1.08%; 26 extensional
roots accounting for 27 FE frequencies below 1120 Hz, one root covering
the 1114.9/1115.3 Hz pair, worst unique-partner error 0.052%) and, for the
in-plane spectra, against two independently published semi-analytical
solutions with no shared code or mesh. One-to-one flexural adjudication
with zero false negatives among REAL-like candidates holds at all four
corners of {Ri/Ro = 0.5, 0.6} × {2Θ = π/2, π}. The
mode-shape bar classifies the same way at three ν=0.30 geometries and
at all 24 extensional keys of the ν=0.35 geometry sweep (it still
splits for 2Θ/π ≤ 1.0; at 225° and 270° every tight FE-frequency match
is ARTIFACT-like).

**What we report against ourselves.** A fully default single blind pass of
our own detector misses 4 of the 26 in-plane modes; all four are shown to be
genuine roots of the same determinant, which is why we argue the
discrimination instruments are part of the method, not optional
diagnostics. One near-degenerate assignment above the validated window is
left unforced rather than resolved. The ν=0.30 MAC calibration window
is not a universal floor; the production threshold was not retuned.
Flexural tables are FE-checked under a portable residual screen;
extensional tables are currently FE-dependent for the discrimination step.

**Companion paper.** A rectangular FFFF companion, *Exact Wave-Function
Solution of Completely Free Rectangular Plates: Four-Corner Kirchhoff
Jump and Mode Screening*, is submitted simultaneously to this journal.
It is not an expansion of the present manuscript: the rectangular
assembler is closed-form with four free corners, and it uses a
persistence screen rather than the residual/SUBDOM instruments
developed here. Both papers share the public repository below; this
manuscript matches tag `v1.0.0`.

The solver, the finite-element decks, the validated tables and a from-scratch
regression suite are released under an MIT license as a public GitHub
repository, https://github.com/McCune1/plate_solver; a Zenodo DOI will be
added at submission. The manuscript is original,
is not under consideration elsewhere, and all authors have approved this
submission.

Yours sincerely,
Garrek McCune (corresponding author), Daniel Stutts
Mechanical and Aerospace Engineering, Missouri University of Science and
Technology
