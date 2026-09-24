# -*- coding: utf-8 -*-
"""
plate_solver.piezo_disk -- SC-coupled Kirchhoff solid-disk 3x3 (and
elastic 2x2) for Paper 4 lever 11 / Liu 2002 CPT.

Architecture (mirrors Paper 3 ring_disk -> disk_L_mp):
  - Annular SC is 6x6 = 3 chi-branches x {Z1, Z2} (regular + singular).
  - Disk SC is 3x3 = 3 chi-branches x {Z1 only} (drop Y_n / K_n).
  - Elastic disk is 2x2 = {I_n, J_n} regular pair (drop K_n / Y_n).

Reuses PiezoOutOfPlaneSolver material reduction, chi-cubic, branches,
K_pref, A1/A2/d1/d2, and _radial_quad via composition. Disk methods
never evaluate at r_i and never call annular 4x4/6x6.

Outer BCs (Kirchhoff + SC electrical):
  C: w = w' = phi' = 0
  S: w = M_r = phi' = 0   (classical simply-supported + SC)

Hard refuse: disk API if the caller tries to pass r_i > 0 as a disk
geometry (constructor locks r_i = 0). Flag DISK_6X6_AT_RI0_IS_NOT_THE_DISK
stays True -- do not evaluate annular coupled_det at r_i=0.

Worker path: this module is main-process only. It is intentionally NOT
wired into plate_solver.workers / ProcessPool; reconstruct from scalars
if a future lever needs a worker (LESSONS_LEARNED.md Sec 18.236).

Projection: default 'consistent' (same default as PiezoOutOfPlaneSolver
since LESSONS Sec 18.231). 'duan' remains available. On Liu Tables 2/6
the two differ by <= 3.4e-6 relative, and the whole electromechanical
contribution is < 4e-5, so Liu CPT checks the disk geometry path and the
laminate reduction, NOT the coupling or the projection.

SOLVER_VERSION is not bumped (isotropic defaults unchanged).
"""
from __future__ import annotations

from mpmath import mp, mpf, matrix

from .piezo_solver import PiezoOutOfPlaneSolver


class DiskPathError(ValueError):
    """Illegal disk / annular geometry path (Paper 3 / lever 11 refuse)."""


class PiezoDiskSC:
    """Solid-disk SC-coupled OOP flexural solver (Liu 2002 CPT column).

    Parameters match PiezoOutOfPlaneSolver materials. Geometry is a solid
    disk of outer radius r_o (no hole). Internally a PiezoOutOfPlaneSolver
    is constructed with a dummy annular r_i used ONLY as a host for
    shared cubic/branch/material helpers -- disk determinants never read
    that r_i and never call annular dets.
    """

    DISK_6X6_AT_RI0_IS_NOT_THE_DISK = True

    def __init__(self, r_o, h, E, nu, rho,
                 h1=0.0,
                 C11E=None, C12E=None, C13E=None, C33E=None,
                 e31=None, e33=None, X11=None, X33=None,
                 rho_pzt=None, dps=100, projection='consistent'):
        r_o = float(r_o)
        if r_o <= 0.0:
            raise DiskPathError("disk needs r_o > 0, got %r" % (r_o,))
        # Dummy annular hole for the shared math host ONLY. Never used by
        # disk_* determinants. Keep r_i clearly away from 0 so accidental
        # annular calls do not silently become "disk at ri=0".
        dummy_ri = 0.5 * r_o
        self._host = PiezoOutOfPlaneSolver(
            r_i=dummy_ri, r_o=r_o, h=h, E=E, nu=nu, rho=rho,
            h1=h1, C11E=C11E, C12E=C12E, C13E=C13E, C33E=C33E,
            e31=e31, e33=e33, X11=X11, X33=X33, rho_pzt=rho_pzt,
            dps=dps, projection=projection)
        self.r_i = 0.0
        self.r_o = r_o
        self.h = float(h)
        self.h1 = float(h1)
        self.dps = int(dps)
        self.projection = projection

    # ---------- refuse helpers ----------

    def _require_disk_geometry(self):
        if not self.DISK_6X6_AT_RI0_IS_NOT_THE_DISK:
            raise DiskPathError(
                "DISK_6X6_AT_RI0_IS_NOT_THE_DISK must stay True")
        if float(self.r_i) > 0.0:
            raise DiskPathError(
                "disk methods refuse r_i > 0 (got r_i=%r); "
                "disk is centre-regularity 2x2/3x3, not a hole"
                % (self.r_i,))

    @staticmethod
    def refuse_annular_at_ri0():
        """Documented hard refuse: never call annular 6x6 at r_i=0."""
        raise DiskPathError(
            "FAIL_DISK_IS_6x6: annular coupled/cc_coupled_det at r_i=0 "
            "is NOT the disk (Paper 3 job 2339453 / lever 11 refuse)")

    # ---------- regular radial column (drop singular partner) ----------

    @staticmethod
    def _regular_radial(n, r, lam):
        """Return (Z, dZ) for the bounded-at-0 member of the Bessel pair.

        lam >= 0 -> modified: I_n (drop K_n)
        lam <  0 -> ordinary: J_n (drop Y_n)

        Same recurrence as PiezoOutOfPlaneSolver._radial_quad Z1/dZ1.
        """
        r = mpf(r)
        if lam.real >= 0:
            delta = mp.sqrt(lam)
            x = delta * r
            Z = mp.besseli(n, x)
            dZ = delta * mpf('0.5') * (
                mp.besseli(n - 1, x) + mp.besseli(n + 1, x))
        else:
            delta = mp.sqrt(-lam)
            x = delta * r
            Z = mp.besselj(n, x)
            dZ = delta * mpf('0.5') * (
                mp.besselj(n - 1, x) - mp.besselj(n + 1, x))
        return Z, dZ

    # ---------- elastic disk 2x2 (Gate 2a elastic limit) ----------

    def elastic_disk_det(self, omega, n, bc_outer="C", h1=None):
        """2x2 elastic disk det on regular pair {I_n, J_n}.

        bc_outer='C': rows w, w'
        bc_outer='S': rows w, M_r

        h1=None uses instance h1; h1=0.0 is bare host (Gate 1 match).

        Warning: do not call PiezoOutOfPlaneSolver annular dets at r_i=0
        as a disk substitute (DISK_6X6_AT_RI0_IS_NOT_THE_DISK).
        """
        self._require_disk_geometry()
        bc = str(bc_outer).upper()
        if bc not in ("C", "S"):
            raise ValueError("elastic disk outer allows C/S; got %r" % bc_outer)
        if h1 is None:
            h1 = self.h1
        H = self._host
        with mp.workdps(self.dps):
            omega = mpf(omega)
            d1 = H._d1()
            d2 = H._d2(h1)
            dsum = d1 + d2
            A1 = H._A1(d2, h1)
            A2 = H._A2(h1)
            mu = mp.sqrt(A2 * omega ** 2 / dsum)
            # Regular columns only: I (lam=+mu), J (lam=-mu)
            r = mpf(self.r_o)
            cols = []
            for lam in (mu, -mu):
                Z, dZ = self._regular_radial(n, r, lam)
                cols.append((Z, dZ, lam))

            w_row = [Z for (Z, dZ, lam) in cols]
            wp_row = [dZ for (Z, dZ, lam) in cols]
            if bc == "C":
                M = matrix([w_row, wp_row])
            else:
                m_row = []
                for (Z, dZ, lam) in cols:
                    m_row.append(
                        (dsum * lam + 2 * A1 * n ** 2 / r ** 2) * Z
                        - (2 * A1 / r) * dZ)
                M = matrix([w_row, m_row])
            if M.rows != 2 or M.cols != 2:
                raise DiskPathError(
                    "FAIL_DISK_IS_NOT_2x2: elastic_disk_det got %sx%s"
                    % (M.rows, M.cols))
            for i in range(2):
                s = max(abs(M[i, j]) for j in range(2)) or mpf(1)
                for j in range(2):
                    M[i, j] = M[i, j] / s
            return mp.det(M)

    def elastic_disk_bisect(self, lo, hi, n, bc_outer="C", h1=None, iters=50):
        if h1 is None:
            h1 = self.h1
        with mp.workdps(self.dps):
            lo = mpf(lo)
            hi = mpf(hi)
            flo = self.elastic_disk_det(lo, n, bc_outer=bc_outer, h1=h1)
            for _ in range(iters):
                mid = (lo + hi) / 2
                fm = self.elastic_disk_det(mid, n, bc_outer=bc_outer, h1=h1)
                if (fm.real > 0) == (flo.real > 0):
                    lo, flo = mid, fm
                else:
                    hi = mid
            return float((lo + hi) / 2)

    # ---------- SC coupled disk 3x3 ----------

    def sc_disk_det(self, omega, n, bc_outer="C", h1=None):
        """3x3 SC-coupled disk det: 3 chi-branches x 1 regular radial each.

        bc_outer='C': rows w, w', phi'
        bc_outer='S': rows w, M_r, phi'

        h1=0 refused (chi-cubic singularity) -- use elastic_disk_det.
        """
        self._require_disk_geometry()
        bc = str(bc_outer).upper()
        if bc not in ("C", "S"):
            raise ValueError("SC disk outer allows C/S; got %r" % bc_outer)
        if h1 is None:
            h1 = self.h1
        if float(h1) == 0.0:
            raise ValueError(
                "sc_disk_det needs h1 > 0 -- h1=0 is a chi-cubic singularity; "
                "use elastic_disk_det(h1=0.0) for the bare-host limit")
        H = self._host
        with mp.workdps(self.dps):
            omega = mpf(omega)
            h1 = mpf(h1)
            d1 = H._d1()
            d2 = H._d2(h1)
            A1v = H._A1(d2, h1)
            K_pref = H._K_pref(h1)
            lams = H._branches(omega, h1)
            if len(lams) != 3:
                raise DiskPathError(
                    "expected 3 chi-branches, got %d" % len(lams))

            r = mpf(self.r_o)
            w_row, wp_row, m_row, phip_row = [], [], [], []
            for chi, lam in lams:
                Z, dZ = self._regular_radial(n, r, lam)
                Ki = (d1 + d2) * lam + K_pref * chi
                w_row.append(Z)
                wp_row.append(dZ)
                m_row.append(
                    (Ki + 2 * A1v * n ** 2 / r ** 2) * Z
                    - (2 * A1v / r) * dZ)
                phip_row.append(chi * dZ)

            if bc == "C":
                M = matrix([w_row, wp_row, phip_row])
            else:
                M = matrix([w_row, m_row, phip_row])
            if M.rows != 3 or M.cols != 3:
                raise DiskPathError(
                    "FAIL_DISK_IS_NOT_3x3: sc_disk_det got %sx%s"
                    % (M.rows, M.cols))
            # Column then row equilibration (same as cc_coupled_det)
            for j in range(3):
                mx = max(abs(M[i, j]) for i in range(3)) or mpf(1)
                for i in range(3):
                    M[i, j] = M[i, j] / mx
            for i in range(3):
                mx = max(abs(M[i, j]) for j in range(3)) or mpf(1)
                for j in range(3):
                    M[i, j] = M[i, j] / mx
            return mp.det(M)

    def sc_disk_bisect(self, lo, hi, n, bc_outer="C", h1=None, iters=45):
        if h1 is None:
            h1 = self.h1
        with mp.workdps(self.dps):
            lo = mpf(lo)
            hi = mpf(hi)
            flo = self.sc_disk_det(lo, n, bc_outer=bc_outer, h1=h1)
            for _ in range(iters):
                mid = (lo + hi) / 2
                fm = self.sc_disk_det(mid, n, bc_outer=bc_outer, h1=h1)
                if (fm.real > 0) == (flo.real > 0):
                    lo, flo = mid, fm
                else:
                    hi = mid
            return float((lo + hi) / 2)

    def sc_disk_logabsdet(self, omega, n, bc_outer="C", h1=None):
        d = self.sc_disk_det(omega, n, bc_outer=bc_outer, h1=h1)
        ad = abs(complex(d))
        if ad <= 0.0:
            return -300.0
        return float(mp.log10(ad))


# Liu / project material bundle (e31 = +4.1 project convention)
LIU_DISK_KWARGS = dict(
    h=0.01,          # host half-thickness (m); full host 2h = 20 mm
    E=200e9,
    nu=0.3,
    rho=7800.0,
    h1=0.002,        # each piezo skin (m)
    C11E=132e9,
    C12E=71e9,
    C13E=73e9,
    C33E=115e9,
    e31=4.1,         # project sign (+); Liu Table 1 prints -4.1
    e33=14.1,
    X11=7.124e-9,
    X33=5.841e-9,
    rho_pzt=7500.0,
)


def make_liu_disk(r_o=0.6, dps=80, projection='consistent', **overrides):
    """Factory: Liu Table 1 steel+PZT-4 solid disk at given r_o."""
    kw = dict(LIU_DISK_KWARGS)
    kw.update(overrides)
    return PiezoDiskSC(r_o=r_o, dps=dps, projection=projection, **kw)
