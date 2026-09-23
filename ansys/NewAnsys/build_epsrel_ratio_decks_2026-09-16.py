#!/usr/bin/env python3
"""Build the F-F epsrel LANB deck for the two remaining Duan2005 Table 4
thickness ratios (h1/2h=1/8 and 1/5), by cloning the cluster-confirmed
h1/2h=1/12 deck (job 2490706, LESSONS_LEARNED.md Sec 18.165) -- same mesh
topology, same mm-N-s-tonne rescale, same relative-permittivity fix, only
H1 (piezo layer thickness) and the closed-form target change.

Closed-form full 3-branch coupled targets (probe_piezo_p4_ff_forward_model_
2026-09-15.py, dps=100, p=0, re-run 2026-09-16 to confirm and get exact
Hz conversion -- see NOTE ON TARGET Hz below):
  h1/2h=1/12: omega=747.9929 rad/s  (already FE-confirmed, job 2490706)
  h1/2h=1/8:  omega=764.1336 rad/s
  h1/2h=1/5:  omega=800.3738 rad/s

NOTE ON TARGET Hz: the h1/2h=1/12 thread has cited "119.058 Hz" as
747.9929 rad/s's closed-form Hz value in ~15 places since Sec 18.153.
Direct recomputation (747.9929/(2*pi) = 119.04677 Hz) does not reproduce
that figure -- a ~0.01% labeling slip, not a physics error (job 2490706's
own PASS verdict, 119.242 vs either label, is unaffected: 0.15% vs
119.058, 0.164% vs 119.04677, still a clean FE match either way). NOT
re-opening that thread or touching any existing file/citation -- flagging
only so the NEW targets below use directly-computed Hz values rather than
propagating the same slip forward:
  h1/2h=1/8:  764.1336 / (2*pi) = 121.61564 Hz
  h1/2h=1/5:  800.3738 / (2*pi) = 127.38345 Hz

Every substitution is an exact unique-string replace: mismatch or a
second hit aborts. Output is LF. Jobnames after run_ansys_queue.sh's
24-char truncation are asserted unique against the whole epsrel batch
(existing four decks + these two new ones).
"""
from __future__ import annotations

import math
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
NEWANSYS_DIR = os.environ.get(
    "NEWANSYS_DIR",
    os.path.expanduser("~/mnt/Plate_Solver_Package/Plate_Solver_Package/Ansys/NewAnsys"),
)

PARENT = "ansys_p4_epsrel_ff_lanb_2026-09-16.inp"

RATIOS = [
    # (denom, stem_tag, omega_rad_s)
    (8, "h18", 764.1336),
    (5, "h15", 800.3738),
]

EXISTING_JOBNAME_SOURCE_DECKS = [
    "ansys_p4_epsrel_ff_e0_2026-09-16.inp",
    "ansys_p4_epsrel_ff_lanb_2026-09-16.inp",
    "ansys_p4_epsrel_ff_subsp_2026-09-16.inp",
    "ansys_p4_epsrel_cc_lanb_2026-09-16.inp",
]


def load(path: str) -> str:
    with open(path, "r", encoding="utf-8", newline="") as f:
        text = f.read()
    text = text.replace("\r\n", "\n").replace("\r", "\n")
    if not text.endswith("\n"):
        text += "\n"
    return text


def save(path: str, text: str) -> None:
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write(text)
    print(f"wrote {path} ({text.count(chr(10))} lines)")


def once(text: str, old: str, new: str, label: str) -> str:
    n = text.count(old)
    if n != 1:
        raise SystemExit(f"FATAL {label}: expected 1 occurrence of {old!r}, found {n}")
    return text.replace(old, new, 1)


def jobname(fname: str) -> str:
    stem = os.path.basename(fname)
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


RATIO_BANNER_TMPL = """!
! ============================================================
! RATIO CLONE (2026-09-16): cloned from the cluster-confirmed h1/2h=1/12
! epsrel-fix deck (job 2490706, LESSONS_LEARNED.md Sec 18.165 -- relative
! permittivity, mm-N-s-tonne rescale, FIX v6 TB,PIEZ; ALL of that is
! unchanged here). Only H1 (piezo layer thickness) and the closed-form
! target below differ. Do not read the "119.058 Hz" / "1/12" references
! in the historical header above this line as applying to this deck --
! they describe the ORIGINAL h1/2h=1/12 derivation session verbatim.
!
! THIS DECK: h1/2h={ratio_label}, closed-form full 3-branch coupled root
! = {omega} rad/s = {hz} Hz (probe_piezo_p4_ff_forward_model_2026-09-15.py,
! dps=100; Hz computed directly as omega/(2*pi), not copied from the
! h1/12 thread's "119.058" label -- see build script docstring for why).
! PASS criterion: a mode within a few percent of {hz} Hz, real O(0.1-1)
! bending discriminator (UY_top~UY_bot, UX_top~-UX_bot at r=RO), and the
! 0 Hz rigid-body mode present. Read this deck's own _modes.txt /
! _discriminator.txt / _out.txt -- do not trust the .QUEUE_OK sentinel
! alone (project standing rule).
! ============================================================
"""


def make_ratio_deck(denom: int, stem_tag: str, omega: float) -> tuple[str, str]:
    hz = omega / (2.0 * math.pi)
    ratio_label = f"1/{denom}"
    parent_path = os.path.join(NEWANSYS_DIR, PARENT)
    text = load(parent_path)

    text = once(
        text,
        "H1     = 2*H/12         ! piezo layer thickness, Table-4 ratio 1/12",
        f"H1     = 2*H/{denom}          ! piezo layer thickness, Table-4 ratio {ratio_label}",
        f"{stem_tag} H1 assignment",
    )
    text = once(
        text,
        "/TITLE, Paper 4 F-F piezo ring EPSREL LANB, p=0, h1/2h=1/12, mm-N-s-tonne (2026-09-16)",
        f"/TITLE, Paper 4 F-F piezo ring EPSREL LANB, p=0, h1/2h={ratio_label}, mm-N-s-tonne (2026-09-16)",
        f"{stem_tag} TITLE",
    )
    text = once(
        text,
        "FINISH\n/CLEAR,NOSTART\n/TITLE,",
        RATIO_BANNER_TMPL.format(ratio_label=ratio_label, omega=f"{omega:.4f}", hz=f"{hz:.5f}")
        + "FINISH\n/CLEAR,NOSTART\n/TITLE,",
        f"{stem_tag} banner insert",
    )

    new_stem = f"ansys_p4_epsrel_ff_lanb_{stem_tag}"
    text = once(
        text,
        "*CFOPEN,ansys_p4_epsrel_ff_lanb_modes,txt",
        f"*CFOPEN,{new_stem}_modes,txt",
        f"{stem_tag} modes CFOPEN",
    )
    text = once(
        text,
        "*CFOPEN,ansys_p4_epsrel_ff_lanb_discriminator,txt",
        f"*CFOPEN,{new_stem}_discriminator,txt",
        f"{stem_tag} discriminator CFOPEN",
    )
    text = once(
        text,
        "('mode f[Hz] tag=epsrel_ff_lanb target=119.058Hz')",
        f"('mode f[Hz] tag=epsrel_ff_lanb_{stem_tag} target={hz:.5f}Hz')",
        f"{stem_tag} modes VWRITE tag",
    )

    out_fname = f"ansys_p4_epsrel_ff_lanb_{stem_tag}_2026-09-16.inp"
    return out_fname, text


def main() -> None:
    outputs: dict[str, str] = {}
    for denom, stem_tag, omega in RATIOS:
        fname, text = make_ratio_deck(denom, stem_tag, omega)
        outputs[fname] = text

    all_jobname_sources = EXISTING_JOBNAME_SOURCE_DECKS + list(outputs)
    print("jobname uniqueness (against full epsrel batch):")
    assert_unique_jobnames(all_jobname_sources)

    for fname, text in outputs.items():
        if "\r" in text:
            raise SystemExit(f"FATAL {fname}: CR leaked into output")
        if "2*H/12" in text:
            raise SystemExit(f"FATAL {fname}: leftover 1/12 H1 assignment")
        if "target=119.058Hz" in text:
            raise SystemExit(f"FATAL {fname}: leftover 1/12 VWRITE target tag")
        save(os.path.join(NEWANSYS_DIR, fname), text)
    print(f"OK: {len(outputs)} decks")


if __name__ == "__main__":
    main()
