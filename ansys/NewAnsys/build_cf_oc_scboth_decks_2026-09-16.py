#!/usr/bin/env python3
"""Optional lever-10 FE: clone the confirmed F-F OC/scboth pair
(jobs 2491761 / 2491831) and change ONLY the mechanical BC to C-F
(inner clamped, outer free -- Paper 3 ring_disk.py convention).

Mesh, mm-N-s-tonne, TB,PIEZ FIX v6, relative permittivity, LANB,
KEYOPT(1)=1001, and each deck's electrical BC are UNCHANGED. Do not
re-open those.

Mechanical clamp matches the C-C deck's inner-radius rows
(ansys_p4_epsrel_cc_lanb_2026-09-16.inp): D,ALL,UX and D,ALL,UY at
X=RI through the stack. Do NOT D,ALL (that would also pin VOLT on the
inner-radius outer-electrode corners and kill the OC bus). Outer
radius stays free.

Closed-form targets, PiezoOutOfPlaneSolver, h1/2h=1/12, n=0,
job 2495528 / in-sandbox SENTINEL (LESSONS_LEARNED.md Sec 18.192):
  scboth / elastic bilayer  421.894681 rad/s = 67.144 Hz
  OC                        423.085352 rad/s = 67.333 Hz
Analytic split +0.282%.

PASS (pre-registered; read each deck's own _modes.txt /
_discriminator.txt / _out.txt; judge MAPDL from NUMBER OF ERROR
MESSAGES ENCOUNTERED and RUN COMPLETED, not the queue script's grep
of the MAPDL error token -- that token must not appear in this
deck's comments):
  - No 0 Hz rigid-body mode (clamp took). FAIL if mode 1 is 0 Hz.
  - scboth: first bending mode near 67.144 Hz.
  - OC: first bending mode ABOVE the same-mesh scboth, near 67.333 Hz.
  - Bending discriminator at r=RO (UY_top ~ UY_bot, UX_top ~ -UX_bot).
  - FAIL if OC reproduces scboth Hz -- electrical BC did not take.
Do not retune mesh / units / TBDATA / LANB if the split is a few
percent off closed form (F-F recovered 93% of +0.301% on this mesh).

Output is LF. Jobnames after run_ansys_queue.sh's 24-char truncation
are asserted unique against the existing p4_epsrel batch.
"""
from __future__ import annotations

import glob
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
OC_PARENT = os.path.join(HERE, "ansys_p4_epsrel_ff_lanb_oc_2026-09-16.inp")
SC_PARENT = os.path.join(HERE, "ansys_p4_epsrel_ff_lanb_scboth_2026-09-16.inp")
OC_OUT = os.path.join(HERE, "ansys_p4_epsrel_cf_lanb_oc_2026-09-16.inp")
SC_OUT = os.path.join(HERE, "ansys_p4_epsrel_cf_lanb_scboth_2026-09-16.inp")

MECH = """! Mechanical BC -- C-F (Paper 3: inner clamped, outer free). UX and UY
! only at r=RI through the stack, same rows as the C-C deck's inner
! radius. Do not D,ALL (that would also pin VOLT on the inner-radius
! outer-electrode corners and kill the OC bus). Outer radius stays free.
NSEL,S,LOC,X,RI
D,ALL,UX,0
D,ALL,UY,0
ALLSEL,ALL
"""

OC_BANNER = r"""!
! ============================================================
! C-F OPEN-CIRCUIT CLONE (2026-09-16): cloned from the confirmed
! F-F OC pair (job 2491761, ansys_p4_epsrel_ff_lanb_oc_2026-09-16.inp).
! Mesh, mm-N-s-tonne, relative permittivity, FIX v6 TB,PIEZ, LANB,
! KEYOPT(1)=1001, and the OC electrical BC: UNCHANGED.
!
! THIS DECK changes ONLY the mechanical BC to C-F (inner clamped,
! outer free). Closed-form target: oc_cf_bisect n=0 h1/2h=1/12
!   omega = 423.085352 rad/s = 67.333 Hz
! Fair SC partner is ansys_p4_epsrel_cf_lanb_scboth_2026-09-16.inp
! (closed-form 67.144 Hz). PASS: first bending mode ABOVE that
! same-mesh scboth, near 67.333 Hz, no 0 Hz rigid body (clamp took),
! bending discriminator at r=RO. Read this deck's own _modes.txt /
! _discriminator.txt / _out.txt, never the .QUEUE_OK sentinel alone.
! ============================================================
"""

SC_BANNER = r"""!
! ============================================================
! C-F BOTH-FACES SC CLONE (2026-09-16): cloned from the confirmed
! F-F scboth pair (job 2491831,
! ansys_p4_epsrel_ff_lanb_scboth_2026-09-16.inp). Mesh/units/TB,PIEZ/
! LANB/KEYOPT(1)=1001 and both-faces VOLT=0: UNCHANGED.
!
! THIS DECK changes ONLY the mechanical BC to C-F (inner clamped,
! outer free). Closed-form target: elastic_cf_bisect / cf_coupled_bisect
! n=0 h1/2h=1/12, 421.895 rad/s = 67.144 Hz. Fair SC partner of the
! C-F OC deck. PASS: first bending mode near 67.144 Hz, no 0 Hz rigid
! body. Read this deck's own files, never the .QUEUE_OK sentinel alone.
! ============================================================
"""


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


def load(path: str) -> str:
    with open(path, "r", encoding="utf-8", newline="") as f:
        text = f.read()
    text = text.replace("\r\n", "\n").replace("\r", "\n")
    if not text.endswith("\n"):
        text += "\n"
    return text


def existing_jobnames():
    seen = {}
    for path in glob.glob(os.path.join(HERE, "ansys_p4_epsrel_*.inp")):
        seen[jobname(path)] = os.path.basename(path)
    return seen


def write_oc(text: str) -> str:
    text = once(
        text,
        "FINISH\n/CLEAR,NOSTART\n/TITLE,",
        OC_BANNER + "FINISH\n/CLEAR,NOSTART\n/TITLE,",
        "OC banner",
    )
    text = once(
        text,
        "/TITLE, Paper 4 F-F piezo ring EPSREL LANB OPEN-CIRCUIT, p=0, h1/2h=1/12, mm-N-s-tonne (2026-09-16)",
        "/TITLE, Paper 4 C-F piezo ring EPSREL LANB OPEN-CIRCUIT, p=0, h1/2h=1/12, mm-N-s-tonne (2026-09-16)",
        "OC TITLE",
    )
    text = once(
        text,
        "! TRUE FREE-FREE: since R_i=0.1 > 0, the ring does NOT touch the axis of\n"
        "! revolution, so (unlike the disk smoke deck, which needed UX=0 at r=0 to\n"
        "! suppress the on-axis singularity) NO mechanical constraint is applied\n"
        "! anywhere -- this is a genuinely unconstrained modal solve. Expect rigid-\n"
        "! body translation (and, if p>0 were modeled, rotation) modes near f=0 Hz\n"
        "! at the bottom of the spectrum; these are correct, not a bug (same\n"
        "! precedent as the disk smoke deck's own mode-1-at-0-Hz note).",
        "! C-F MECHANICAL BC: inner radius r=RI is clamped (UX=UY=0 through the\n"
        "! stack); outer radius is free. No 0 Hz rigid-body mode is expected.\n"
        "! Do not pin VOLT at r=RI (UX and UY only).",
        "OC free-free header",
    )
    text = once(
        text,
        "CP,1,VOLT,ALL\nALLSEL,ALL\nSOLVE\n",
        "CP,1,VOLT,ALL\nALLSEL,ALL\n" + MECH + "SOLVE\n",
        "OC mechanical BC",
    )
    text = once(
        text,
        "*CFOPEN,ansys_p4_epsrel_ff_lanb_oc_modes,txt",
        "*CFOPEN,ansys_p4_epsrel_cf_lanb_oc_modes,txt",
        "OC modes CFOPEN",
    )
    text = once(
        text,
        "('mode f[Hz] tag=epsrel_ff_lanb_oc target=119.405Hz')",
        "('mode f[Hz] tag=epsrel_cf_lanb_oc target=67.333Hz')",
        "OC modes tag",
    )
    text = once(
        text,
        "*CFOPEN,ansys_p4_epsrel_ff_lanb_oc_discriminator,txt",
        "*CFOPEN,ansys_p4_epsrel_cf_lanb_oc_discriminator,txt",
        "OC discriminator CFOPEN",
    )
    text = once(
        text,
        "! unconstrained (one CP set = parallel bus). Do NOT also D,VOLT,0 the\n"
        "! outers -- that restores short-circuit. No mechanical constraint.",
        "! unconstrained (one CP set = parallel bus). Do NOT also D,VOLT,0 the\n"
        "! outers -- that restores short-circuit. Mechanical BC is C-F, below.",
        "OC solu comment",
    )
    return text


def write_sc(text: str) -> str:
    text = once(
        text,
        "FINISH\n/CLEAR,NOSTART\n/TITLE,",
        SC_BANNER + "FINISH\n/CLEAR,NOSTART\n/TITLE,",
        "SC banner",
    )
    text = once(
        text,
        "/TITLE, Paper 4 F-F piezo ring EPSREL LANB BOTH-FACES SC, p=0, h1/2h=1/12, mm-N-s-tonne (2026-09-16)",
        "/TITLE, Paper 4 C-F piezo ring EPSREL LANB BOTH-FACES SC, p=0, h1/2h=1/12, mm-N-s-tonne (2026-09-16)",
        "SC TITLE",
    )
    text = once(
        text,
        "! TRUE FREE-FREE: since R_i=0.1 > 0, the ring does NOT touch the axis of\n"
        "! revolution, so (unlike the disk smoke deck, which needed UX=0 at r=0 to\n"
        "! suppress the on-axis singularity) NO mechanical constraint is applied\n"
        "! anywhere -- this is a genuinely unconstrained modal solve. Expect rigid-\n"
        "! body translation (and, if p>0 were modeled, rotation) modes near f=0 Hz\n"
        "! at the bottom of the spectrum; these are correct, not a bug (same\n"
        "! precedent as the disk smoke deck's own mode-1-at-0-Hz note).",
        "! C-F MECHANICAL BC: inner radius r=RI is clamped (UX=UY=0 through the\n"
        "! stack); outer radius is free. No 0 Hz rigid-body mode is expected.\n"
        "! Do not pin VOLT at r=RI (UX and UY only).",
        "SC free-free header",
    )
    text = once(
        text,
        "NSEL,A,LOC,Y,-H-H1\nD,ALL,VOLT,0\nALLSEL,ALL\nSOLVE\n",
        "NSEL,A,LOC,Y,-H-H1\nD,ALL,VOLT,0\nALLSEL,ALL\n" + MECH + "SOLVE\n",
        "SC mechanical BC",
    )
    text = once(
        text,
        "*CFOPEN,ansys_p4_epsrel_ff_lanb_scboth_modes,txt",
        "*CFOPEN,ansys_p4_epsrel_cf_lanb_scboth_modes,txt",
        "SC modes CFOPEN",
    )
    text = once(
        text,
        "('mode f[Hz] tag=epsrel_ff_lanb_scboth target=119.047Hz')",
        "('mode f[Hz] tag=epsrel_cf_lanb_scboth target=67.144Hz')",
        "SC modes tag",
    )
    text = once(
        text,
        "*CFOPEN,ansys_p4_epsrel_ff_lanb_scboth_discriminator,txt",
        "*CFOPEN,ansys_p4_epsrel_cf_lanb_scboth_discriminator,txt",
        "SC discriminator CFOPEN",
    )
    return text


def main() -> None:
    forbidden = "*** ERROR ***"
    seen = existing_jobnames()
    jobs = []
    for parent, out, fn in (
            (OC_PARENT, OC_OUT, write_oc),
            (SC_PARENT, SC_OUT, write_sc)):
        text = fn(load(parent))
        if forbidden in text:
            raise SystemExit("FATAL %s contains the MAPDL error token"
                             % os.path.basename(out))
        jn = jobname(out)
        if jn in seen:
            raise SystemExit("FATAL jobname collision %r with %s"
                             % (jn, seen[jn]))
        seen[jn] = os.path.basename(out)
        jobs.append((jn, os.path.basename(out)))
        with open(out, "w", encoding="utf-8", newline="\n") as f:
            f.write(text)
        print("wrote %s (%d lines) jobname=%r"
              % (os.path.basename(out), text.count("\n"), jn))
    print("jobnames unique:")
    for jn, name in jobs:
        print("  %r <- %s" % (jn, name))


if __name__ == "__main__":
    main()
