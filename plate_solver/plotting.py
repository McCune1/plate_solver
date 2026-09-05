# -*- coding: utf-8 -*-
"""
plate_solver.plotting -- dispersion-branch and mode-shape plotting
(Figures 2 & 3 of the reference papers). Best-effort, exception-safe;
requires matplotlib.

Extracted verbatim (line-range provenance in LESSONS_LEARNED Sec. 22); no
numeric behaviour changed, SOLVER_VERSION not bumped.
"""
from __future__ import annotations
import os
import time
import numpy as np
from mpmath import mp, mpf, mpc, matrix

from .config import PI
from .detectors import full_search, track, select_fill, real_basis_items
from .dispersion import _Part1Fast, _Part2Fast

def compute_dispersion_branches(eng_fast, Om_lo, Om_hi, n_pts=60, xmax=14.0,
                                 max_seconds=180, max_full_search=30,
                                 max_branches=None):
    """Sweep Omega in [Om_lo, Om_hi] and track each canonical root (xi or
    zeta) as a continuous branch.  Returns a list of branches, each a list
    of (Om, Re(z), Im(z)) tuples — the raw data for a Figure 2 style plot.

    eng_fast : a _Part1Fast or _Part2Fast instance (float64 'fast' engine).
    max_branches : if set, keep only the lowest-lying `max_branches` branches
        (ranked by |z| at the smallest Ω where they exist).  The paper's
        Fig. 2 deliberately shows only "the first three imaginary (or complex)
        branches near Ω̃=0"; plotting every branch within `xmax` (including
        far-out high-|Im| branches that are nearly vertical lines) is what
        made the saved PNGs look cluttered/wrong vs. the paper's clean bowls.

    BUGFIX (this pass): the previous version matched branches by their
    *position* in the freshly-sorted root list at each step, not by
    continuity.  Whenever the root count changed (which happens often —
    measured directly: count went 6→4→4→...→3→3 monotonically shrinking
    over a 15-point sweep on the canonical geometry, never recovering) a
    branch's plotted line would silently jump from one physical root to an
    unrelated one, producing the jagged/spiky look in the saved PNGs (very
    unlike the paper's smooth Fig. 2 bowls).  Fixed by greedy nearest-
    neighbor assignment in the (Re, Im) plane between consecutive steps,
    with unmatched new roots starting fresh branches and unmatched old
    roots simply ending (no jump).  A cheaper full_search variant (reduced
    grid density — this is for a qualitative plot, not a determinant-root
    accuracy requirement) keeps the larger max_full_search budget
    affordable; max_seconds is still an absolute backstop.
    """
    t0 = time.time()
    # Sampling plan: a uniform sweep over the FULL range plus extra refinement
    # near Om_lo (where the parabolas turn over toward Re=0 and cross to the
    # imaginary axis — the part the old tracker lost).  We process points in a
    # COARSE-TO-FINE order so that, if the time budget is hit, what we DO have
    # still spans the whole Ω range (coarsely) instead of being truncated at the
    # top — the cause of the "lines too short / cut off near the top" figures.
    # A full root search at every Ω is much costlier than the old Newton
    # continuation, so we also give it a more realistic default budget.
    span = max(Om_hi - Om_lo, 1e-9)
    base = list(np.linspace(Om_lo, Om_hi, n_pts))
    # Dense geometric refinement over the low-Ω region where every branch turns
    # over toward the origin (the real parabola vertex and the imaginary descent
    # both live below the first cut-off).  Capturing this region densely is what
    # lets the real branches reach small Re (so they stitch through the corner)
    # and the imaginary branches bend into the origin — the part the coarse
    # uniform sweep misses.  Cover the first 25% of the range with ~n_pts//2
    # geometrically-clustered points.
    low = list(Om_lo + span * np.linspace(0.0, 0.25, max(8, n_pts // 2)) ** 1.5)
    all_pts = sorted(set([round(float(o), 9) for o in base + low]))

    def _coarse_to_fine(pts):
        # Reorder so endpoints and a coarse spanning subset come first, then the
        # gaps are filled at progressively finer stride.
        n = len(pts); order = []; seen = set()
        stride = max(1, n - 1)
        while stride >= 1:
            for i in range(0, n, stride):
                if i not in seen:
                    seen.add(i); order.append(pts[i])
            if stride == 1:
                break
            stride = max(1, stride // 2)
        for i in range(n):
            if i not in seen:
                order.append(pts[i])
        return order

    Oms = _coarse_to_fine(all_pts)

    def _search(Om):
        # Root-search density (2026-07-16): was ngrid=6/nax=140, a deliberate
        # cheapening vs. full_search's own ngrid=10/nax=240 default, "purely
        # for a continuous qualitative curve". In practice the coarser grid
        # sometimes failed to resolve a branch's root at one Om sample while
        # finding it fine on either side; _value_connect (below) then ends
        # that branch right there and starts a new one when the root
        # reappears, producing visible gaps in the plotted line that aren't
        # in the paper's Fig. 2. Restored to full_search's own default
        # density -- still float64/fast-engine, not the mp-precision
        # discovery path, so this stays cheap; override via
        # FIG2_NGRID/FIG2_NAX if needed.
        ngrid = int(os.environ.get("FIG2_NGRID", "10"))
        nax = int(os.environ.get("FIG2_NAX", "240"))
        return full_search(eng_fast, float(Om), xmax=xmax, ngrid=ngrid, nax=nax)

    # ── Collect every root at every Ω (no tracking yet) ───────────────────────
    snaps = []   # list of (Om, [complex roots])
    for Om in Oms:
        if time.time() - t0 > max_seconds:
            break
        snaps.append((float(Om), _search(Om)))
    if not snaps:
        return []
    # Connect in Ω order regardless of the (coarse-to-fine) processing order.
    snaps.sort(key=lambda s: s[0])

    # ── Connect roots into branches PER AXIS by nearest value ────────────────
    # Real-axis roots (|Im|≈0) live on the Re plane; imaginary-axis roots
    # (|Re|≈0) on the Im plane; everything else is a complex branch.  Roots on a
    # single axis are ORDERED and never cross, so connecting them by nearest
    # value across Ω is far more robust than 2-D Newton continuation (which kept
    # fragmenting the fundamental near Re→0, so the parabolas never reached the
    # origin and so never visually met the imaginary branches — the paper's
    # "right lines attach to the left lines").
    axis_tol = float(os.environ.get("FIG2_AXIS_TOL", "0.25"))   # on-axis cutoff
    gap_tol  = float(os.environ.get("FIG2_GAP_TOL", "0.9"))     # max value jump

    def _value_connect(seq):
        # seq: list of (Om, [values]) sorted by Om.  Returns list of branches,
        # each [(Om, value), ...], via greedy nearest-value matching.
        branches_v = []
        active_v = []   # (branch_idx, last_value)
        for Om, vals in seq:
            vals = sorted(vals)
            new_active = []
            used = [False] * len(vals)
            # match existing branches to nearest unused value
            order = sorted(range(len(active_v)),
                           key=lambda k: active_v[k][1])
            for k in order:
                bi, lastv = active_v[k]
                best = -1; bestd = gap_tol
                for j, v in enumerate(vals):
                    if used[j]:
                        continue
                    d = abs(v - lastv)
                    if d < bestd:
                        bestd = d; best = j
                if best >= 0:
                    used[best] = True
                    branches_v[bi].append((Om, vals[best]))
                    new_active.append((bi, vals[best]))
            for j, v in enumerate(vals):
                if not used[j]:
                    branches_v.append([(Om, v)])
                    new_active.append((len(branches_v) - 1, v))
            active_v = new_active
        return branches_v

    real_seq = [(Om, [z.real for z in rs
                      if abs(z.imag) < axis_tol and z.real > -1e-9])
                for Om, rs in snaps]
    imag_seq = [(Om, [abs(z.imag) for z in rs
                      if abs(z.real) < axis_tol and abs(z.imag) > 1e-9])
                for Om, rs in snaps]

    real_br = [[(Om, v, 0.0) for Om, v in b] for b in _value_connect(real_seq)]
    imag_br = [[(Om, 0.0, v) for Om, v in b] for b in _value_connect(imag_seq)]

    # Complex branches (off both axes): 2-D nearest-neighbor tracking.
    cplx_snaps = [(Om, [z for z in rs
                        if abs(z.imag) >= axis_tol and abs(z.real) >= axis_tol])
                  for Om, rs in snaps]
    cbr = []
    cactive = []   # (branch_idx, last_complex)
    for Om, zs in cplx_snaps:
        new_active = []
        used = [False] * len(zs)
        pairs = sorted(((abs(zs[j] - lv), k, j)
                        for k, (bi, lv) in enumerate(cactive)
                        for j in range(len(zs))), key=lambda t: t[0])
        matched_k = set()
        for d, k, j in pairs:
            if k in matched_k or used[j] or d > 1.0:
                continue
            matched_k.add(k); used[j] = True
            bi = cactive[k][0]
            cbr[bi].append((Om, zs[j].real, zs[j].imag))
            new_active.append((bi, zs[j]))
        for j, z in enumerate(zs):
            if not used[j]:
                cbr.append([(Om, z.real, z.imag)])
                new_active.append((len(cbr) - 1, z))
        cactive = new_active

    real_br = [b for b in real_br if len(b) >= 3]
    imag_br = [b for b in imag_br if len(b) >= 3]
    cbr     = [b for b in cbr if len(b) >= 3]

    # ── PAPER-FAITHFUL BRANCH SELECTION (near-origin only; Fig-2) ─────────────
    # Keep branches that approach the origin (the paper's first few near-Ω̃=0
    # branches) plus real propagating parabolas; drop far-out constant evanescent
    # clutter.  near_tol scales with the plotted extent.  See §Figures.
    def _min_abs(b):
        return min((pt[1] ** 2 + pt[2] ** 2) ** 0.5 for pt in b)

    # Adaptive origin window: a fraction of the largest |Im| we would otherwise
    # plot, floored at a small absolute value so tiny-Ω geometries still work.
    im_extent = max([pt[2] for b in imag_br + cbr for pt in b] + [1.0])
    near_tol = float(os.environ.get("FIG2_NEAR_ORIGIN",
                                    str(max(2.5, 0.45 * im_extent))))
    imag_br = [b for b in imag_br if _min_abs(b) <= near_tol]
    cbr     = [b for b in cbr     if _min_abs(b) <= near_tol]

    # ── Stitch real↔imaginary branch continuations at the cut-offs (Fig-2) ──────
    # A branch is propagating (real ξ) above its cut-off and evanescent (imag ξ)
    # below; they meet at (0,0,Ω_cut).  Join them through that inserted point for
    # the paper's continuous swoosh.  FIG2_STITCH_TOL.
    stitch_tol = float(os.environ.get("FIG2_STITCH_TOL", "0.9"))
    stitched = []
    used_imag = set()
    for rb in real_br:
        # endpoint of the real branch nearest the axis (smallest Re)
        i_lo = min(range(len(rb)), key=lambda i: rb[i][1])
        re_lo, om_re = rb[i_lo][1], rb[i_lo][0]
        # Find the imaginary branch whose lowest point sits at (Im≈0) near the
        # same Ω as this real branch's lowest point — that is the cut-off where
        # the propagating parabola becomes evanescent.  We gate on the IMAGINARY
        # side reaching the origin (im_lo small) and the two meeting Ω's being
        # close, NOT on re_lo being tiny: a sparsely-sampled real parabola may
        # have its lowest captured Re a little above zero, yet it still physically
        # continues to the corner, so we prepend the exact origin point to draw
        # it through (Re=0,Im=0,Ω_cut).
        best = None; bestd = 1e9
        for k, ib in enumerate(imag_br):
            if k in used_imag:
                continue
            j_lo = min(range(len(ib)), key=lambda j: ib[j][2])
            im_lo, om_im = ib[j_lo][2], ib[j_lo][0]
            if im_lo < stitch_tol and abs(om_im - om_re) < 0.20 * span:
                d = im_lo + 2.0 * abs(om_im - om_re) / span + 0.3 * re_lo
                if d < bestd:
                    bestd = d; best = (k, j_lo)
        if best is None or re_lo > 1.6:
            stitched.append(rb)
        else:
            k, j_lo = best; ib = imag_br[k]; used_imag.add(k)
            # imaginary part ordered high→low toward the meet, then the real part
            im_sorted = sorted(ib, key=lambda p: -p[2])
            re_sorted = sorted(rb, key=lambda p: p[1])
            om_meet = 0.5 * (im_sorted[-1][0] + re_sorted[0][0])
            stitched.append(im_sorted + [(om_meet, 0.0, 0.0)] + re_sorted)
    # keep any imaginary branches that were not stitched to a real one
    stitched += [ib for k, ib in enumerate(imag_br) if k not in used_imag]

    branches = stitched + cbr

    # The per-axis builder above already yields a clean, de-fragmented set, and
    # the near-origin filter has removed the far-out vertical clutter; keep all
    # remaining branches by default and only trim if there are an unusually large
    # number — then by length (longest = most physical), preserving the stitched
    # swooshes and the near-origin imaginary curves the paper shows on the left.
    if max_branches is None:
        cap = int(os.environ.get("FIG2_MAX_TOTAL", "10"))
        if len(branches) > cap:
            branches = sorted(branches, key=len, reverse=True)[:cap]
    elif len(branches) > max_branches:
        branches = sorted(branches, key=len, reverse=True)[:max_branches]
    return branches


def plot_dispersion_curves(geom, mat, part, Om_lo, Om_hi, fname,
                            n_pts=48, xmax=14.0, max_seconds=1200, title=None,
                            om_plot_hi=None):
    """Figure 2 style 3D plot: Omega vs Re(branch) and Im(branch).

    part : 1 (out-of-plane, xi) or 2 (in-plane, zeta).
    xmax : branch-search radius — pass the SAME xi_max/ze_max used for the
           mode search on this geometry, so the plotted dispersion curves
           and the modes found by find_modes_sigmin are based on the same
           branch inventory (previously this silently defaulted to 14.0
           regardless of the geometry's actual search radius).
    max_seconds : wall-clock budget for the branch sweep.  The per-axis tracer
           does a full root search at every Ω (robust but ~10x costlier than the
           old Newton continuation), so the default is generous; the sweep is
           coarse-to-fine so a cutoff still spans the full Ω range.  Override
           with FIG2_MAX_SECONDS.
    Saves a PNG to `fname`.  Returns True on success, False (with a printed
    warning) on any failure — never raises.
    """
    try:
        _mxs = os.environ.get("FIG2_MAX_SECONDS")
        if _mxs:
            max_seconds = float(_mxs)
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

        r0b = float(geom.r0_bar)
        nu  = float(mat.nu_bar)
        if part == 1:
            eng = _Part1Fast(r0b, nu, M=80)
            zlabel = r"$\bar{\kappa}\Omega$"
        else:
            eng = _Part2Fast(r0b, nu, float(mat.c11_eff_bar), float(mat.R), M=80)
            zlabel = r"$\bar{\Omega}$"

        # Branch keeping: by default use the per-class caps inside
        # compute_dispersion_branches (FIG2_MAX_REAL / _IMAG / _CPLX), which
        # preserves the imaginary-plane evanescent curves the paper shows on the
        # left.  Only fall back to a single global cap if FIG2_MAX_BRANCHES is
        # explicitly set (legacy behaviour).
        _mbr = os.environ.get("FIG2_MAX_BRANCHES")
        max_br = int(_mbr) if _mbr not in (None, "") else None
        # The paper's Fig. 2 vertical axis stops far below the full mode-scan
        # range (≈3rd cut-off for Part 1, ≈2nd for Part 2).  Plotting to the full
        # Om_hi (here up to 2.5) stretches the curves and adds high-Ω clutter the
        # paper omits.  Clip the sweep to om_plot_hi when provided.
        plot_hi = (Om_hi if om_plot_hi is None
                   else max(Om_lo + 1e-6, min(Om_hi, float(om_plot_hi))))
        branches = compute_dispersion_branches(eng, Om_lo, plot_hi, n_pts=n_pts,
                                                 xmax=xmax, max_seconds=max_seconds,
                                                 max_branches=max_br)
        if not branches:
            print(f"  WARNING: plot_dispersion_curves — no branches found for {fname}")
            return False

        fig = plt.figure(figsize=(7, 6))
        ax = fig.add_subplot(111, projection='3d')
        for b in branches:
            Oms_b  = [p[0] for p in b]
            Re_b   = [p[1] for p in b]
            Im_b   = [p[2] for p in b]
            ax.plot(Re_b, Im_b, Oms_b, color='k', linewidth=0.8)

        re_label = r"$\Re(\xi)$" if part == 1 else r"$\Re(\zeta)$"
        im_label = r"$\Im(\xi)$" if part == 1 else r"$\Im(\zeta)$"
        ax.set_xlabel(re_label)
        ax.set_ylabel(im_label)
        ax.set_zlabel(zlabel)
        if title:
            ax.set_title(title)
        # ORIGIN-AT-CORNER (paper Fig. 2): the paper places (Re=0, Im=0) exactly
        # at the box corner where the two horizontal axes meet, with the vertical
        # axis rising from it.  matplotlib autoscaling instead started each axis
        # at the data minimum, so the "0" ticks floated mid-edge and did not
        # coincide — the readability problem.  Forcing both horizontal limits to
        # start at 0 (and z at 0) pins the origin to the corner.  Axis maxima are
        # rounded up from the plotted data so the box hugs the curves like the
        # paper's.
        re_all = [p[1] for b in branches for p in b] or [1.0]
        im_all = [p[2] for b in branches for p in b] or [1.0]
        xmax_ax = max(1.0, float(np.ceil(max(re_all))))
        ymax_ax = max(1.0, float(np.ceil(max(im_all))))
        ax.set_xlim(0.0, xmax_ax)
        ax.set_ylim(0.0, ymax_ax)
        ax.set_zlim(0.0, plot_hi)
        # Part 1 (paper): Im on the front-left, Re on the front-right.  Part 2
        # (paper): mirrored — Re on the front-left, Im on the front-right.
        ax.view_init(elev=18, azim=(-130 if part == 1 else -50))
        fig.tight_layout()
        fig.savefig(fname, dpi=150)
        plt.close(fig)
        print(f"  Saved dispersion-curve plot: {fname}")
        return True
    except Exception as exc:
        print(f"  WARNING: plot_dispersion_curves failed for {fname}: {exc}")
        return False


def plot_mode_shape(solver, Om_star, brs, n_dofs, mode_num, geom_label, fname,
                     n_r=25, n_th=25):
    """Figure 3 style wireframe plot of one mode shape.

    solver  : OutOfPlaneSolver (mode-shape reconstruction is implemented for
              Part 1 only; Part 2 in-plane shapes are not yet supported and
              this function returns False without error in that case).
    brs     : converged branch-root set at Om_star (float64 canonical roots).
    Saves a PNG to `fname`.  Returns True on success, False (with a printed
    warning) on any failure — never raises.
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

        if not hasattr(solver, 'mode_shape_grid'):
            print(f"  Mode-shape plotting not supported for this solver "
                  f"({type(solver).__name__}); skipping {fname}")
            return False

        grid = solver.mode_shape_grid(Om_star, brs, n_dofs, n_r=n_r, n_th=n_th)
        if grid is None:
            print(f"  WARNING: plot_mode_shape — reconstruction failed for {fname}")
            return False
        r_phys, th_phys, Wgrid = grid

        # r_phys is the PHYSICAL radius in [R_i, R_o]; th_phys the PHYSICAL angle
        # in radians [-Theta,+Theta].  Cartesian conversion then yields the TRUE
        # annular-sector planform (a ring segment with a central hole), matching
        # the paper's Fig. 3 instead of a unit pie wedge.
        R, TH = np.meshgrid(r_phys, th_phys, indexing='ij')
        X = R * np.cos(TH)
        Y = R * np.sin(TH)
        Z = Wgrid.copy()

        # CLAMPED-EDGE DATUM (optional, DEFAULT OFF).  Wall at θ=-Θ (col 0); the weak
        # form (Eq. 43) leaves a nonzero pointwise W there (structural, mp-confirmed).
        # Not forced flat by default (pinning re-creates the free-edge spike).
        # MODE_SHAPE_PIN_CLAMP=1 to pin.  See §Mode-shapes.
        pin = os.environ.get("MODE_SHAPE_PIN_CLAMP", "0") == "1"
        try:
            pin_max = float(os.environ.get("MODE_SHAPE_PIN_MAX", "0.20"))
        except ValueError:
            pin_max = 0.20
        denom = float(np.max(np.abs(Z))) or 1.0
        clamp_resid = float(np.max(np.abs(Z[:, 0]))) / denom
        if pin and clamp_resid <= pin_max:
            Z[:, 0] = 0.0
        elif pin:
            print(f"  note: clamped-edge residual {clamp_resid:.0%} of peak for "
                  f"{os.path.basename(fname)} exceeds {pin_max:.0%}; NOT pinning "
                  f"(structural weak-BC residual — see handoff).")

        # CLAMP-RAMP DISPLAY TRANSFORM (optional, DEFAULT OFF) ─────────────────
        # Purely-visual wall-flattening fade applied to the reconstructed grid (never
        # the coefficients).  MODE_SHAPE_CLAMP_RAMP=1.  Rationale: §Mode-shapes.
        if os.environ.get("MODE_SHAPE_CLAMP_RAMP", "0") == "1" and not (
                pin and clamp_resid <= pin_max):
            th0, th1 = float(th_phys[0]), float(th_phys[-1])
            width = (th1 - th0) or 1.0
            g = (th1 - th_phys) / width            # 1 at wall (col0) → 0 at free
            c_r = Z[:, 0].copy()                   # clamp residual per radius
            Z_ramp = Z - np.outer(c_r, g)
            raw_peak = float(np.max(np.abs(Z))) or 1.0
            new_peak = float(np.max(np.abs(Z_ramp))) or 1.0
            try:
                ramp_max_growth = float(
                    os.environ.get("MODE_SHAPE_CLAMP_RAMP_MAXGROWTH", "1.15"))
            except ValueError:
                ramp_max_growth = 1.15
            if new_peak <= ramp_max_growth * raw_peak:
                Z = Z_ramp
            else:
                print(f"  note: clamp-ramp would inflate peak "
                      f"({new_peak / raw_peak:.2f}×) for "
                      f"{os.path.basename(fname)} (twist-dominated mode); "
                      f"keeping raw shape.")

        fig = plt.figure(figsize=(6, 5))
        ax = fig.add_subplot(111, projection='3d')

        # UNDEFORMED REFERENCE ("shadow"): the flat sector at Z=0, drawn as a
        # light outline beneath the deformed surface so the deflection can be
        # read against the rigid shape — exactly as in the paper's Fig. 3.  Use
        # the ordered perimeter (inner arc -> free radial edge -> outer arc ->
        # clamped radial edge) so it is one clean closed loop, not a full grid.
        px = np.concatenate([X[0, :], X[:, -1], X[-1, ::-1], X[::-1, 0]])
        py = np.concatenate([Y[0, :], Y[:, -1], Y[-1, ::-1], Y[::-1, 0]])
        ax.plot(px, py, np.zeros_like(px), color="0.55", linewidth=1.0)

        # Deformed mode shape.
        ax.plot_wireframe(X, Y, Z, color='k', linewidth=0.6, rstride=1, cstride=1)
        ax.set_title(f"{geom_label}  Mode {mode_num}  "
                      f"$\\Omega$={Om_star:.4f}")
        # Keep the sector to true in-plane scale (so the arcs stay circular, not
        # stretched into ellipses).  Give Z a modest height so the deflection is
        # legible above the flat reference, and view it from a near-plan oblique
        # angle (paper Fig. 3) with the clamped edge reading as the horizon.
        # Z-ratio lowered 0.55->0.18 (2026-07-16): Wgrid is unit-peak normalized
        # above, so 0.55 stretched every mode shape to ~55% of the plate's own
        # planform width regardless of the physical deflection pattern, reading
        # as a twisted ribbon rather than a gentle bend -- visibly steeper than
        # Seok & Tiersten's own Fig. 3 proportions. 0.18 is a starting point
        # calibrated by eye against Fig. 3, not derived from the paper; tune via
        # MODE_SHAPE_Z_RATIO if it still looks off. Purely cosmetic -- Wgrid
        # itself (and every computed frequency) is unchanged, no SOLVER_VERSION
        # implication.
        try:
            z_ratio = float(os.environ.get("MODE_SHAPE_Z_RATIO", "0.18"))
        except ValueError:
            z_ratio = 0.18
        try:
            xs = float(X.max() - X.min()); ys = float(Y.max() - Y.min())
            ax.set_box_aspect((max(xs, 1e-9), max(ys, 1e-9), z_ratio * max(xs, ys)))
        except Exception:
            pass
        ax.view_init(elev=28, azim=-72)
        ax.set_axis_off()
        fig.tight_layout()
        fig.savefig(fname, dpi=150)
        if fname.endswith(".png"):
            pdf_name = fname[:-4] + ".pdf"
            fig.savefig(pdf_name)
            print(f"  Saved mode-shape plot: {fname}")
            print(f"  Saved mode-shape vector: {pdf_name}")
        else:
            print(f"  Saved mode-shape plot: {fname}")
        plt.close(fig)
        return True
    except Exception as exc:
        print(f"  WARNING: plot_mode_shape failed for {fname}: {exc}")
        return False




# ==============================================================================
#  RESEARCH40 OVERNIGHT ALL-STEPS DRIVER
#  (1) cantilever regression spot-check (a few modes; confirms Step-3/6 neutral)
#  (2) parallel FREE-FREE OOP+IP spectra across the geometry sweep (Step 6)
#  (3) literature-comparison tabulation in Omega_lit = w*R_o^2*sqrt(rho*H/D)
#  (4) orthotropic-material scaffold + isotropic-reduction sanity test
#  Every section is best-effort (never aborts the run); checkpoints per geometry.
#  Env knobs (all optional):
#    R40_OVERNIGHT=1     run this driver instead of the legacy cantilever sweep
#    R40_SPOTCHECK=1     cantilever spot-check (default ON)
#    R40_FREEFREE=1      free-free sweep (default ON)
#    R40_ORTHO_TEST=1    orthotropic scaffold sanity (default ON)
#    R40_FF_MODES=6      free-free modes wanted per geometry
#    R40_FF_P1=1.25:0.5,1.25:1.0,1.25:1.5,1.25:0.25,1.25:0.75,1.25:1.25
#    R40_FF_P2=1.25:0.5,1.25:1.0,1.25:0.25     (r0_2b:twoTheta_over_pi, comma-sep)
#    FAST=1              halve scan ranges/points for a quick shake-out
# ==============================================================================

