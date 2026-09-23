#!/usr/bin/env python3
"""Clone the cluster-confirmed F-F epsrel LANB h1/2h=1/12 SC deck
(job 2490706) into a BOTH-FACES-GROUNDED short-circuit control.

Job 2491761's OC FE (inners grounded, outers CP-free) landed at
119.143 Hz -- 0.22% below the OC closed-form 119.405 Hz, and 0.083%
BELOW the parent SC FE 119.242 Hz (outers grounded, interface free).
Those two FE decks do not share an electrode pair, so the missing
stiffening is uninterpretable without this control.

THIS DECK: VOLT=0 on BOTH the interface (z=+-h) AND the outers
(z=+-(h+h1)) -- Duan/analytic SC (both electrode faces of each layer
shorted). Same mesh/units/TB,PIEZ/LANB as parent and as the OC deck.
No CP.

Closed-form target: coupled_bisect / sine-SC = 747.9929 rad/s
= 119.047 Hz. Pre-registered reading (decide before the job returns):
  - If mode 2 ~ 119.05-119.14 Hz and BELOW the OC FE 119.143 Hz, then
    job 2491761 DID show OC stiffening vs true SC; the parent "SC" FE
    (outer-only ground) was a mixed BC, not the analytic SC.
  - If mode 2 ~ 119.14 Hz (matches OC FE), then FE does not resolve
    the 0.30% analytic OC-SC split on this mesh.
  - If mode 2 ~ 119.24 Hz (matches parent), something else is wrong.

Do not re-open mesh/units/TBDATA/LANB.
"""
from __future__ import annotations

import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
PARENT = os.path.join(HERE, "ansys_p4_epsrel_ff_lanb_2026-09-16.inp")
OUT = os.path.join(HERE, "ansys_p4_epsrel_ff_lanb_scboth_2026-09-16.inp")

EXISTING = [
    "ansys_p4_epsrel_ff_e0_2026-09-16.inp",
    "ansys_p4_epsrel_ff_lanb_2026-09-16.inp",
    "ansys_p4_epsrel_ff_subsp_2026-09-16.inp",
    "ansys_p4_epsrel_cc_lanb_2026-09-16.inp",
    "ansys_p4_epsrel_ff_lanb_h18_2026-09-16.inp",
    "ansys_p4_epsrel_ff_lanb_h15_2026-09-16.inp",
    "ansys_p4_epsrel_cc_lanb_rmid_2026-09-16.inp",
    "ansys_p4_epsrel_ff_lanb_oc_2026-09-16.inp",
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
! BOTH-FACES SC CONTROL (2026-09-16): cloned from job 2490706's
! confirmed h1/2h=1/12 F-F LANB. Mesh/units/TB,PIEZ/LANB unchanged.
! VOLT=0 on interface AND outers -- analytic/Duan short-circuit
! (both faces of each piezo layer shorted). Fair SC counterpart to
! job 2491761's OC deck (inners grounded, outers CP-free).
! Target: 119.047 Hz (747.9929 rad/s). See build script docstring
! for the pre-registered reading of 119.05 vs 119.14 vs 119.24.
! ============================================================
"""

SC_OUTER_ONLY = """! Electrical BC only -- short-circuit on each piezo layer's OUTER
! (non-interface) face. NO mechanical constraint anywhere (true free-free,
! see header).
NSEL,S,LOC,Y,H+H1
D,ALL,VOLT,0
NSEL,S,LOC,Y,-H-H1
D,ALL,VOLT,0
ALLSEL,ALL
"""

SC_BOTH = """! Electrical BC -- BOTH-FACES short-circuit control (Sec 18.177).
! VOLT=0 on host/piezo interface AND on both outer faces. No CP.
! Fair SC counterpart of the OC deck (job 2491761).
NSEL,S,LOC,Y,H
NSEL,A,LOC,Y,-H
NSEL,A,LOC,Y,H+H1
NSEL,A,LOC,Y,-H-H1
D,ALL,VOLT,0
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
        "/TITLE, Paper 4 F-F piezo ring EPSREL LANB BOTH-FACES SC, p=0, h1/2h=1/12, mm-N-s-tonne (2026-09-16)",
        "TITLE",
    )
    text = once(text, SC_OUTER_ONLY, SC_BOTH, "SOLU electrical BC")
    text = once(
        text,
        "*CFOPEN,ansys_p4_epsrel_ff_lanb_modes,txt",
        "*CFOPEN,ansys_p4_epsrel_ff_lanb_scboth_modes,txt",
        "modes CFOPEN",
    )
    text = once(
        text,
        "('mode f[Hz] tag=epsrel_ff_lanb target=119.058Hz')",
        "('mode f[Hz] tag=epsrel_ff_lanb_scboth target=119.047Hz')",
        "modes VWRITE tag",
    )
    text = once(
        text,
        "*CFOPEN,ansys_p4_epsrel_ff_lanb_discriminator,txt",
        "*CFOPEN,ansys_p4_epsrel_ff_lanb_scboth_discriminator,txt",
        "discriminator CFOPEN",
    )

    jn = jobname(os.path.basename(OUT))
    seen = {jobname(n): n for n in EXISTING}
    if jn in seen:
        raise SystemExit("FATAL jobname collision %r with %s" % (jn, seen[jn]))
    print("jobname %r <- %s" % (jn, os.path.basename(OUT)))

    with open(OUT, "w", encoding="utf-8", newline="\n") as f:
        f.write(text)
    print("wrote %s (%d lines)" % (OUT, text.count("\n")))


if __name__ == "__main__":
    main()
