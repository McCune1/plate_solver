#!/usr/bin/env python3
"""Clone the cluster-confirmed F-F epsrel LANB h1/2h=1/12 deck
(job 2490706, LESSONS_LEARNED.md Sec 18.165) into an OPEN-CIRCUIT
electrical-BC variant. Mesh, units, TB,PIEZ FIX v6, relative permittivity,
LANB, and free-free mechanical BC are unchanged. Only the electrical BC
and the closed-form target change.

Parent SC electrical BC (do not re-open): VOLT=0 on each piezo layer's
OUTER face (z=+- (h+h1)); host/piezo interface unconstrained.

This deck's OC electrical BC (LESSONS_LEARNED.md Sec 18.175): VOLT=0 on
the host/piezo INTERFACE (z=+-h, inner electrodes grounded); outer faces
equipotential and unconstrained (one CP,VOLT set spanning both outers --
parallel bus). Matches the analytic F-F OC model (linear potential, Q=0).

Closed-form target: oc_ff_bisect n=0 h1/2h=1/12
  omega = 750.2394210586 rad/s = 119.4053 Hz
Parent SC FE (job 2490706) mode 2 = 119.242 Hz. OC must land ABOVE that
(stiffening) and near 119.405 Hz. A result that reproduces 119.242 Hz
is a FAIL -- the electrical BC did not take.

Output is LF. Jobname after run_ansys_queue.sh's 24-char truncation is
checked unique against the existing epsrel batch.
"""
from __future__ import annotations

import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
PARENT = os.path.join(HERE, "ansys_p4_epsrel_ff_lanb_2026-09-16.inp")
OUT = os.path.join(HERE, "ansys_p4_epsrel_ff_lanb_oc_2026-09-16.inp")

EXISTING = [
    "ansys_p4_epsrel_ff_e0_2026-09-16.inp",
    "ansys_p4_epsrel_ff_lanb_2026-09-16.inp",
    "ansys_p4_epsrel_ff_subsp_2026-09-16.inp",
    "ansys_p4_epsrel_cc_lanb_2026-09-16.inp",
    "ansys_p4_epsrel_ff_lanb_h18_2026-09-16.inp",
    "ansys_p4_epsrel_ff_lanb_h15_2026-09-16.inp",
    "ansys_p4_epsrel_cc_lanb_rmid_2026-09-16.inp",
]


def jobname(fname: str) -> str:
    stem = os.path.basename(fname)
    if stem.endswith(".inp"):
        stem = stem[:-4]
    if stem.startswith("ansys_"):
        stem = stem[6:]
    stem = re.sub(r"[^A-Za-z0-9]", "_", stem)
    return stem[:24]


def once(text: str, old: str, new: str, label: str) -> str:
    n = text.count(old)
    if n != 1:
        raise SystemExit("FATAL %s: expected 1 occurrence of %r, found %d"
                         % (label, old, n))
    return text.replace(old, new, 1)


BANNER = r"""!
! ============================================================
! OPEN-CIRCUIT CLONE (2026-09-16): cloned from the cluster-confirmed
! h1/2h=1/12 epsrel-fix F-F LANB deck (job 2490706, LESSONS_LEARNED.md
! Sec 18.165). Mesh, mm-N-s-tonne, relative permittivity, FIX v6 TB,PIEZ,
! LANB, free-free mechanical BC: UNCHANGED. Do not re-open those.
!
! THIS DECK changes ONLY the electrical BC, to match the analytic F-F
! open-circuit model (Sec 18.175):
!   inner electrodes (host/piezo interface, z=+-h): VOLT=0
!   outer electrodes (z=+-(h+h1)): CP,VOLT as one parallel bus, unconstrained
! Parent SC deck did the opposite (VOLT=0 on the OUTERS, interface free).
!
! Closed-form target: PiezoOutOfPlaneSolver.oc_ff_bisect n=0 h1/2h=1/12
!   omega = 750.2394210586 rad/s = 119.4053 Hz
! Parent SC FE mode 2 = 119.242 Hz (job 2490706). PASS: a bending mode
! near 119.405 Hz, ABOVE 119.242 Hz (stiffening vs the same mesh), 0 Hz
! rigid-body present, real O(0.1-1) bending discriminator, no
! ill-conditioning warning. A result that reproduces 119.242 Hz is FAIL
! -- the OC electrical BC did not take. Read this deck's own _modes.txt
! / _discriminator.txt / _out.txt, never the .QUEUE_OK sentinel alone.
! ============================================================
"""

SC_BC = """! Electrical BC only -- short-circuit on each piezo layer's OUTER
! (non-interface) face. NO mechanical constraint anywhere (true free-free,
! see header).
NSEL,S,LOC,Y,H+H1
D,ALL,VOLT,0
NSEL,S,LOC,Y,-H-H1
D,ALL,VOLT,0
ALLSEL,ALL
"""

OC_BC = """! Electrical BC -- F-F OPEN-CIRCUIT (Sec 18.175). Inner electrodes
! (host/piezo interface) GROUNDED; outer electrodes EQUIPOTENTIAL and
! unconstrained (one CP set = parallel bus). Do NOT also D,VOLT,0 the
! outers -- that restores short-circuit. No mechanical constraint.
NSEL,S,LOC,Y,H
NSEL,A,LOC,Y,-H
D,ALL,VOLT,0
ALLSEL,ALL
NSEL,S,LOC,Y,H+H1
NSEL,A,LOC,Y,-H-H1
CP,1,VOLT,ALL
ALLSEL,ALL
"""


def main() -> None:
    with open(PARENT, "r", encoding="utf-8", newline="") as f:
        text = f.read()
    text = text.replace("\r\n", "\n").replace("\r", "\n")
    if not text.endswith("\n"):
        text += "\n"

    text = once(
        text,
        "FINISH\n/CLEAR,NOSTART\n/TITLE,",
        BANNER + "FINISH\n/CLEAR,NOSTART\n/TITLE,",
        "banner insert",
    )
    text = once(
        text,
        "/TITLE, Paper 4 F-F piezo ring EPSREL LANB, p=0, h1/2h=1/12, mm-N-s-tonne (2026-09-16)",
        "/TITLE, Paper 4 F-F piezo ring EPSREL LANB OPEN-CIRCUIT, p=0, h1/2h=1/12, mm-N-s-tonne (2026-09-16)",
        "TITLE",
    )
    text = once(
        text,
        "! ELECTRICAL BC: short-circuit, VOLT=0 on the OUTER (free) face of each\n"
        "! piezo layer (z=h+h1 top, z=-h-h1 bottom) -- matching Sec 18.141 F's\n"
        "! adopted electrical BC and Duan2005's own stated \"electrodes...shortly\n"
        "! connected\" setup. The host/piezo INTERFACE (z=+-h) is NOT electrically\n"
        "! constrained (matches the sinusoidal potential ansatz's own zero-at-\n"
        "! electrode-faces-only condition, Sec 18.141 B -- the interface is not an\n"
        "! electrode).",
        "! ELECTRICAL BC: OPEN-CIRCUIT, Sec 18.175. VOLT=0 on the host/piezo\n"
        "! INTERFACE (z=+-h, inner electrodes grounded). Outer faces z=+-(h+h1)\n"
        "! are a conducting electrode (CP,VOLT, one parallel bus) left\n"
        "! unconstrained so Q=0 is the natural electrical condition. This is\n"
        "! the opposite of the parent SC deck. Do not re-open mesh/units/TB,PIEZ.",
        "header electrical BC",
    )
    text = once(text, SC_BC, OC_BC, "SOLU electrical BC")
    text = once(
        text,
        "*CFOPEN,ansys_p4_epsrel_ff_lanb_modes,txt",
        "*CFOPEN,ansys_p4_epsrel_ff_lanb_oc_modes,txt",
        "modes CFOPEN",
    )
    text = once(
        text,
        "('mode f[Hz] tag=epsrel_ff_lanb target=119.058Hz')",
        "('mode f[Hz] tag=epsrel_ff_lanb_oc target=119.405Hz')",
        "modes VWRITE tag",
    )
    text = once(
        text,
        "*CFOPEN,ansys_p4_epsrel_ff_lanb_discriminator,txt",
        "*CFOPEN,ansys_p4_epsrel_ff_lanb_oc_discriminator,txt",
        "discriminator CFOPEN",
    )

    jn = jobname(os.path.basename(OUT))
    seen = {jobname(n): n for n in EXISTING}
    if jn in seen:
        raise SystemExit("FATAL jobname collision %r with %s" % (jn, seen[jn]))
    print("jobname %r <- %s" % (jn, os.path.basename(OUT)))
    for n in EXISTING:
        print("  existing %r <- %s" % (jobname(n), n))

    with open(OUT, "w", encoding="utf-8", newline="\n") as f:
        f.write(text)
    print("wrote %s (%d lines)" % (OUT, text.count("\n")))


if __name__ == "__main__":
    main()
