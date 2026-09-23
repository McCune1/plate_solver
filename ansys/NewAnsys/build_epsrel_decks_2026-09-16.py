#!/usr/bin/env python3
"""Build the PLANE223 relative-permittivity diagnostic decks.

Parent decks (already cluster-validated for geometry/mesh/units arithmetic):
  ansys_p4_ff_piezo_ring_UNITS_MM_2026-09-16.inp
  ansys_p4_ff_piezo_coupling_off_UNITS_MM_2026-09-16.inp

Hypothesis being tested (LESSONS_LEARNED.md Sec 18.164 leftover):
  PLANE223 MP,PERx is ALWAYS relative permittivity (Coupled-Field Guide
  2.3.2.1). The parent decks feed Duan2005 Table 1's ABSOLUTE values
  (X11=7.124e-9 F/m, X33=5.841e-9 F/m) into MP,PERX. Coupling-off
  isolation cannot see this (e=0 decouples permittivity from the
  mechanical eigenvalue). This builder changes only that interpretation.

Every substitution is an exact unique-string replace: mismatch or a
second hit aborts. Output is LF. Jobnames after run_ansys_queue.sh's
24-char truncation are asserted unique across the four new decks.
"""
from __future__ import annotations

import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))

EPSREL_BANNER = """!
! ============================================================
! EPSREL FIX (2026-09-16): PLANE223 MP,PERx is RELATIVE permittivity.
! Duan2005 Table 1 X11=7.124e-9 / X33=5.841e-9 are ABSOLUTE (F/m =
! eps_r * eps0). Parent UNITS_MM decks fed those F/m numbers to MP,PERX
! unchanged. Coupled-Field Guide 2.3.2.1: for PLANE223, MP,PERX and
! TB,DPER are ALWAYS relative. The SOLID5/SOLID98 "value<1 means
! absolute" exception does NOT apply. Default EPZRO is left alone
! (8.85e-12 F/m); relative ~805 * EPZRO recovers the Duan absolute
! value, which is also the correct mm-N-s-tonne number (permittivity
! is numerically invariant under that rescale). Host X_INERT is now a
! small relative value (1.0), not 1.0E-9-fed-as-relative. NOTHING else
! in the mesh, TB,PIEZ FIX v6, units, or electrical BC is changed
! unless a deck's own name says so (subsp / e0 / cc).
! ============================================================
"""


def load(fname: str) -> str:
    path = os.path.join(HERE, fname)
    with open(path, "r", encoding="utf-8", newline="") as f:
        text = f.read()
    text = text.replace("\r\n", "\n").replace("\r", "\n")
    if not text.endswith("\n"):
        text += "\n"
    return text


def save(fname: str, text: str) -> None:
    path = os.path.join(HERE, fname)
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write(text)
    print(f"wrote {fname} ({text.count(chr(10))} lines)")


def once(text: str, old: str, new: str, label: str) -> str:
    n = text.count(old)
    if n != 1:
        raise SystemExit(f"FATAL {label}: expected 1 occurrence of {old!r}, found {n}")
    return text.replace(old, new, 1)


def jobname(deck: str) -> str:
    stem = os.path.basename(deck)
    if stem.endswith(".inp"):
        stem = stem[:-4]
    if stem.startswith("ansys_"):
        stem = stem[6:]
    stem = re.sub(r"[^A-Za-z0-9]", "_", stem)
    return stem[:24]


def assert_unique_jobnames(fnames: list[str]) -> None:
    seen: dict[str, str] = {}
    for fname in fnames:
        jn = jobname(fname)
        if jn in seen:
            raise SystemExit(
                f"FATAL jobname collision {jn!r}: {seen[jn]} vs {fname} "
                f"(run_ansys_queue.sh truncates to 24 chars)"
            )
        seen[jn] = fname
        print(f"  jobname {jn!r} <- {fname}")


COUPLED_PARENT = "ansys_p4_ff_piezo_ring_UNITS_MM_2026-09-16.inp"
E0_PARENT = "ansys_p4_ff_piezo_coupling_off_UNITS_MM_2026-09-16.inp"

EPS_BLOCK = (
    "EPS0   = 8.854E-12       ! F/m. PLANE223 MP,PERx is RELATIVE, not absolute.\n"
    "X11ABS = 7.124E-9        ! Duan2005 Table 1 absolute F/m\n"
    "X33ABS = 5.841E-9\n"
    "X11    = X11ABS/EPS0     ! ~804.8 relative\n"
    "X33    = X33ABS/EPS0     ! ~659.7 relative"
)


def patch_permittivity_coupled(text: str) -> str:
    text = once(
        text,
        "X11  = 7.124E-9         ! [F/m] absolute permittivity, numerically UNCHANGED (see units header)\n"
        "X33  = 5.841E-9         ! UNCHANGED",
        EPS_BLOCK,
        "coupled X11/X33 assignment",
    )
    text = once(
        text,
        "X_INERT = 1.0E-9        ! UNCHANGED (FIX v4 value, absolute permittivity, see units header)",
        "X_INERT = 1.0           ! relative (host VOLT well-posed). Was 1.0E-9 fed as if absolute.",
        "coupled X_INERT assignment",
    )
    text = once(
        text,
        "/PREP7\nET,1,PLANE223",
        EPSREL_BANNER + "/PREP7\nET,1,PLANE223",
        "coupled EPSREL banner insert",
    )
    return text


def patch_outputs(text: str, stem: str, modes_tag: str) -> str:
    # Parent coupled CFOPEN names.
    if "ansys_p4_ff_piezo_ring_UNITS_MM_modes" in text:
        text = once(
            text,
            "*CFOPEN,ansys_p4_ff_piezo_ring_UNITS_MM_modes,txt",
            f"*CFOPEN,{stem}_modes,txt",
            f"{stem} modes CFOPEN",
        )
        text = once(
            text,
            "*CFOPEN,ansys_p4_ff_piezo_ring_UNITS_MM_discriminator,txt",
            f"*CFOPEN,{stem}_discriminator,txt",
            f"{stem} discriminator CFOPEN",
        )
        text = once(
            text,
            "('mode      f[Hz]        tag=ff_piezo_layered_ring_p0_h1o2h_1o12 (UNITS_MM rescale, target 119.058 Hz closed-form -- do NOT trust unless stage-1/2 decks confirm first)')",
            f"('{modes_tag}')",
            f"{stem} modes VWRITE tag",
        )
    elif "ansys_p4_ff_piezo_coupling_off_UNITS_MM_modes" in text:
        text = once(
            text,
            "*CFOPEN,ansys_p4_ff_piezo_coupling_off_UNITS_MM_modes,txt",
            f"*CFOPEN,{stem}_modes,txt",
            f"{stem} modes CFOPEN",
        )
        text = once(
            text,
            "*CFOPEN,ansys_p4_ff_piezo_coupling_off_UNITS_MM_discriminator,txt",
            f"*CFOPEN,{stem}_discriminator,txt",
            f"{stem} discriminator CFOPEN",
        )
        text = once(
            text,
            "('mode      f[Hz]        tag=ff_piezo_coupling_off_isolation (UNITS_MM rescale, target 118.81016 Hz = job 2490456's own SI result)')",
            f"('{modes_tag}')",
            f"{stem} modes VWRITE tag",
        )
    else:
        raise SystemExit(f"FATAL {stem}: no parent CFOPEN names found")
    return text


def make_ff_lanb() -> str:
    text = load(COUPLED_PARENT)
    text = patch_permittivity_coupled(text)
    text = once(
        text,
        "/TITLE, Paper 4 F-F piezo layered ring PLANE223 modal, p=0, h1/2h=1/12, mm-N-s-tonne UNITS RESCALE (2026-09-16)",
        "/TITLE, Paper 4 F-F piezo ring EPSREL LANB, p=0, h1/2h=1/12, mm-N-s-tonne (2026-09-16)",
        "ff_lanb TITLE",
    )
    text = patch_outputs(
        text,
        "ansys_p4_epsrel_ff_lanb",
        "mode f[Hz] tag=epsrel_ff_lanb target=119.058Hz",
    )
    return text


def make_ff_subsp() -> str:
    text = make_ff_lanb()
    # TITLE and outputs were already rewritten for lanb; redo those uniquely.
    text = once(
        text,
        "/TITLE, Paper 4 F-F piezo ring EPSREL LANB, p=0, h1/2h=1/12, mm-N-s-tonne (2026-09-16)",
        "/TITLE, Paper 4 F-F piezo ring EPSREL SUBSP, p=0, h1/2h=1/12, mm-N-s-tonne (2026-09-16)",
        "ff_subsp TITLE",
    )
    text = once(
        text,
        "MODOPT,LANB,NMODES,-1.0,,,ON",
        "MODOPT,SUBSP,NMODES,-1.0,,,ON",
        "ff_subsp MODOPT",
    )
    text = once(
        text,
        "*CFOPEN,ansys_p4_epsrel_ff_lanb_modes,txt",
        "*CFOPEN,ansys_p4_epsrel_ff_subsp_modes,txt",
        "ff_subsp modes CFOPEN",
    )
    text = once(
        text,
        "*CFOPEN,ansys_p4_epsrel_ff_lanb_discriminator,txt",
        "*CFOPEN,ansys_p4_epsrel_ff_subsp_discriminator,txt",
        "ff_subsp discriminator CFOPEN",
    )
    text = once(
        text,
        "('mode f[Hz] tag=epsrel_ff_lanb target=119.058Hz')",
        "('mode f[Hz] tag=epsrel_ff_subsp target=119.058Hz')",
        "ff_subsp modes VWRITE tag",
    )
    return text


def make_ff_e0() -> str:
    text = load(E0_PARENT)
    text = once(
        text,
        "X_COMMON = 1.0E-9       ! [F/m]-numerically-UNCHANGED absolute permittivity (see units header)",
        "EPS0   = 8.854E-12       ! F/m. PLANE223 MP,PERx is RELATIVE, not absolute.\n"
        "X11ABS = 7.124E-9        ! Duan2005 Table 1 absolute F/m\n"
        "X33ABS = 5.841E-9\n"
        "X11    = X11ABS/EPS0     ! ~804.8 relative (PZT, e=0 so should NOT move 118.81 Hz)\n"
        "X33    = X33ABS/EPS0     ! ~659.7 relative\n"
        "X_INERT = 1.0            ! relative, host VOLT well-posed",
        "e0 permittivity assignment",
    )
    text = once(
        text,
        "MP,PERX,1,X_COMMON\nMP,PERY,1,X_COMMON\nMP,PERZ,1,X_COMMON",
        "MP,PERX,1,X_INERT\nMP,PERY,1,X_INERT\nMP,PERZ,1,X_INERT",
        "e0 host MP,PER",
    )
    text = once(
        text,
        "MP,PERX,2,X_COMMON\nMP,PERY,2,X_COMMON\nMP,PERZ,2,X_COMMON",
        "MP,PERX,2,X11\nMP,PERY,2,X33\nMP,PERZ,2,X11",
        "e0 pzt MP,PER",
    )
    text = once(
        text,
        "/PREP7\nET,1,PLANE223",
        EPSREL_BANNER + "/PREP7\nET,1,PLANE223",
        "e0 EPSREL banner insert",
    )
    text = once(
        text,
        "/TITLE, Paper 4 F-F piezo-coupling-OFF isolation test, real PZT4 elastic constants, e=0 everywhere, mm-N-s-tonne UNITS RESCALE (2026-09-16)",
        "/TITLE, Paper 4 F-F piezo coupling-OFF EPSREL, e=0, mm-N-s-tonne (2026-09-16)",
        "e0 TITLE",
    )
    text = patch_outputs(
        text,
        "ansys_p4_epsrel_ff_e0",
        "mode f[Hz] tag=epsrel_ff_e0 target=118.81016Hz",
    )
    return text


def make_cc_lanb() -> str:
    text = make_ff_lanb()
    text = once(
        text,
        "/TITLE, Paper 4 F-F piezo ring EPSREL LANB, p=0, h1/2h=1/12, mm-N-s-tonne (2026-09-16)",
        "/TITLE, Paper 4 C-C piezo ring EPSREL LANB, p=0, h1/2h=1/12, mm-N-s-tonne (2026-09-16)",
        "cc_lanb TITLE",
    )
    text = once(
        text,
        "! Electrical BC only -- short-circuit on each piezo layer's OUTER\n"
        "! (non-interface) face. NO mechanical constraint anywhere (true free-free,\n"
        "! see header).\n"
        "NSEL,S,LOC,Y,H+H1\n"
        "D,ALL,VOLT,0\n"
        "NSEL,S,LOC,Y,-H-H1\n"
        "D,ALL,VOLT,0\n"
        "ALLSEL,ALL\n"
        "SOLVE",
        "! Electrical BC: short-circuit on each piezo layer's OUTER face.\n"
        "! Mechanical BC: C-C -- clamp BOTH radial edges (Duan2005 Table 4).\n"
        "! Axisymmetric (r,z): UX=UY=0 at r=RI and r=RO. No 0 Hz rigid-body\n"
        "! mode is expected. Target: 2792 rad/s = 444.40 Hz (p=0,n=0,h1/2h=1/12).\n"
        "NSEL,S,LOC,Y,H+H1\n"
        "D,ALL,VOLT,0\n"
        "NSEL,S,LOC,Y,-H-H1\n"
        "D,ALL,VOLT,0\n"
        "NSEL,S,LOC,X,RI\n"
        "D,ALL,UX,0\n"
        "D,ALL,UY,0\n"
        "NSEL,S,LOC,X,RO\n"
        "D,ALL,UX,0\n"
        "D,ALL,UY,0\n"
        "ALLSEL,ALL\n"
        "SOLVE",
        "cc_lanb C-C BC insert",
    )
    text = once(
        text,
        "*CFOPEN,ansys_p4_epsrel_ff_lanb_modes,txt",
        "*CFOPEN,ansys_p4_epsrel_cc_lanb_modes,txt",
        "cc_lanb modes CFOPEN",
    )
    text = once(
        text,
        "*CFOPEN,ansys_p4_epsrel_ff_lanb_discriminator,txt",
        "*CFOPEN,ansys_p4_epsrel_cc_lanb_discriminator,txt",
        "cc_lanb discriminator CFOPEN",
    )
    text = once(
        text,
        "('mode f[Hz] tag=epsrel_ff_lanb target=119.058Hz')",
        "('mode f[Hz] tag=epsrel_cc_lanb target=444.40Hz')",
        "cc_lanb modes VWRITE tag",
    )
    return text


def main() -> None:
    decks = {
        "ansys_p4_epsrel_ff_lanb_2026-09-16.inp": make_ff_lanb(),
        "ansys_p4_epsrel_ff_subsp_2026-09-16.inp": make_ff_subsp(),
        "ansys_p4_epsrel_ff_e0_2026-09-16.inp": make_ff_e0(),
        "ansys_p4_epsrel_cc_lanb_2026-09-16.inp": make_cc_lanb(),
    }
    names = list(decks)
    print("jobname uniqueness:")
    assert_unique_jobnames(names)
    for fname, text in decks.items():
        if "\r" in text:
            raise SystemExit(f"FATAL {fname}: CR leaked into output")
        if "X11  = 7.124E-9" in text or "X_INERT = 1.0E-9" in text or "X_COMMON = 1.0E-9" in text:
            raise SystemExit(f"FATAL {fname}: leftover absolute-as-relative assignment")
        if "MP,PERX" in text and "X11ABS/EPS0" not in text:
            raise SystemExit(f"FATAL {fname}: permittivity assignment missing")
        save(fname, text)
    print("OK: 4 decks")


if __name__ == "__main__":
    main()
