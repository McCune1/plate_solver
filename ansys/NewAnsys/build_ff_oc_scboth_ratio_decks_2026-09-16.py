#!/usr/bin/env python3
"""Lever 7 (PAPER4_LEVERS_2026-09-16.md sec B.7): four LF clones giving
OC FE and both-faces-SC FE at h1/2h=1/8 and 1/5, completing the fair
OC-vs-scboth pair that only existed at 1/12 (jobs 2491761/2491831,
LESSONS_LEARNED.md Sec 18.177/18.178).

Clones from the two confirmed 1/12 decks -- NOT from the h18/h15 SC
deck family (those are outer-only-ground, a mixed electrical BC, not a
fair OC or SC partner; see Sec 18.178 and the standing ban in
PAPER4_LEVERS_2026-09-16.md). Mesh, mm-N-s-tonne rescale, relative
permittivity fix, TB,PIEZ, LANB, and free-free mechanical BC are
unchanged -- only H1 (piezo layer thickness) and each deck's own
closed-form target change, exactly as build_epsrel_ratio_decks_2026-09-16.py
already did for the (now-superseded-as-a-partner) plain SC family.

Closed-form targets (PiezoOutOfPlaneSolver, dps=100, iters>=50,
computed directly in this session -- the 1/12 OC number below
reproduces LESSONS_LEARNED.md Sec 18.175/18.177's own
750.2394210586 rad/s bit-for-bit, and the 1/8 / 1/5 sine-SC numbers
reproduce build_epsrel_ratio_decks_2026-09-16.py's own 764.1336 /
800.3738 rad/s -- both cross-checks anchor these new numbers to
already-confirmed ones):

  ratio | OC omega (rad/s)   | OC Hz     | scboth (sine-SC) omega (rad/s) | scboth Hz
  1/12  | 750.2394210586     | 119.40431 | 747.9929429356                 | 119.04677
  1/8   | 767.5496387603     | 122.15932 | 764.1335816243                 | 121.61564
  1/5   | 805.9014661114     | 128.26320 | 800.3737601605                 | 127.38344

PASS bars (pre-registered, same discipline as the 1/12 pair, Sec
18.177/18.178): for each ratio, OC mode 2 Hz ABOVE that ratio's own
scboth mode 2 Hz (stiffening on the SAME mesh), OC mode 2 near its own
closed-form Hz (a few percent, matching the 1/12 pair's -0.220%
mesh-vs-closed-form gap), scboth mode 2 within mesh noise of its own
closed-form sine-SC Hz (matching the 1/12 pair's near-bit-identical
118.811 vs elastic-bilayer coupling-off number), 0 Hz rigid-body mode
present in both, real O(0.1-1) bending discriminator (UY_top~UY_bot,
UX_top~-UX_bot at r=RO) in both, 0 MAPDL "*** ERROR ***" in both. A
result that reproduces the WRONG ratio's Hz, or that reproduces the
1/12 pair's Hz at a different ratio, is a FAIL -- read each deck's own
_modes.txt / _discriminator.txt / _out.txt, never the .QUEUE_OK
sentinel alone.

Every substitution is an exact unique-string replace: mismatch or a
second hit aborts. Output is LF. Jobnames after run_ansys_queue.sh's
24-char truncation are asserted unique against the full existing
p4_epsrel batch AND against each other.
"""
from __future__ import annotations

import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))

OC_PARENT = os.path.join(HERE, "ansys_p4_epsrel_ff_lanb_oc_2026-09-16.inp")
SCBOTH_PARENT = os.path.join(HERE, "ansys_p4_epsrel_ff_lanb_scboth_2026-09-16.inp")

# (denom, ratio_label, stem_tag, oc_omega, oc_hz, scboth_omega, scboth_hz)
RATIOS = [
    (8, "1/8", "h18", 767.5496387603, 122.15932, 764.1335816243, 121.61564),
    (5, "1/5", "h15", 805.9014661114, 128.26320, 800.3737601605, 127.38344),
]

# Full existing p4_epsrel jobname-source set (extends
# build_ff_oc_deck_2026-09-16.py's own EXISTING list with the scboth
# and both ratio-family decks that have landed since).
EXISTING = [
    "ansys_p4_epsrel_ff_e0_2026-09-16.inp",
    "ansys_p4_epsrel_ff_lanb_2026-09-16.inp",
    "ansys_p4_epsrel_ff_subsp_2026-09-16.inp",
    "ansys_p4_epsrel_cc_lanb_2026-09-16.inp",
    "ansys_p4_epsrel_ff_lanb_h18_2026-09-16.inp",
    "ansys_p4_epsrel_ff_lanb_h15_2026-09-16.inp",
    "ansys_p4_epsrel_cc_lanb_rmid_2026-09-16.inp",
    "ansys_p4_epsrel_ff_lanb_oc_2026-09-16.inp",
    "ansys_p4_epsrel_ff_lanb_scboth_2026-09-16.inp",
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
    print("wrote %s (%d lines)" % (path, text.count("\n")))


def once(text: str, old: str, new: str, label: str) -> str:
    n = text.count(old)
    if n != 1:
        raise SystemExit("FATAL %s: expected 1 occurrence of %r, found %d"
                         % (label, old, n))
    return text.replace(old, new, 1)


def jobname(fname: str) -> str:
    stem = os.path.basename(fname)
    if stem.endswith(".inp"):
        stem = stem[:-4]
    if stem.startswith("ansys_"):
        stem = stem[6:]
    stem = re.sub(r"[^A-Za-z0-9]", "_", stem)
    return stem[:24]


def assert_unique_jobnames(fnames: list) -> None:
    seen = {}
    for fname in fnames:
        jn = jobname(fname)
        if jn in seen:
            raise SystemExit(
                "FATAL jobname collision %r: %s vs %s (run_ansys_queue.sh "
                "truncates to 24 chars)" % (jn, seen[jn], fname))
        seen[jn] = fname
        print("  jobname %r <- %s" % (jn, fname))


RATIO_BANNER_TMPL = r"""!
! ============================================================
! LEVER 7 RATIO CLONE (2026-09-16): cloned from the h1/2h=1/12 {kind}
! deck (job {parent_job}, LESSONS_LEARNED.md Sec {parent_sec}). Mesh,
! mm-N-s-tonne rescale, relative permittivity fix, TB,PIEZ, LANB, and
! free-free mechanical BC: UNCHANGED. Do not re-open those. Do NOT
! confuse this with the plain h1/2h={ratio_label} SC deck (outer-only
! ground, job 2490848) -- that is a mixed electrical BC and is NOT the
! fair partner for this pair (PAPER4_LEVERS_2026-09-16.md standing ban).
!
! THIS DECK: h1/2h={ratio_label}, {kind}. Closed-form target
! (PiezoOutOfPlaneSolver, dps=100): {omega} rad/s = {hz} Hz.
! PASS: a bending mode near {hz} Hz, {vs_partner}, 0 Hz rigid-body
! present, real O(0.1-1) bending discriminator (UY_top~UY_bot,
! UX_top~-UX_bot at r=RO), 0 MAPDL "*** ERROR ***". Read this deck's
! own _modes.txt / _discriminator.txt / _out.txt, never the .QUEUE_OK
! sentinel alone.
! ============================================================
"""


def make_deck(parent_text: str, kind: str, parent_job: str, parent_sec: str,
              denom: int, ratio_label: str, stem_tag: str,
              omega: float, hz: float, vs_partner: str,
              cfopen_stem: str, vwrite_tag: str, title_suffix: str) -> str:
    text = parent_text

    text = once(
        text,
        "H1     = 2*H/12         ! piezo layer thickness, Table-4 ratio 1/12",
        "H1     = 2*H/%d          ! piezo layer thickness, Table-4 ratio %s"
        % (denom, ratio_label),
        "%s H1 assignment" % stem_tag,
    )
    text = once(
        text,
        "h1/2h=1/12, mm-N-s-tonne (2026-09-16)",
        "h1/2h=%s, mm-N-s-tonne (2026-09-16)" % ratio_label,
        "%s TITLE" % stem_tag,
    )
    text = once(
        text,
        "FINISH\n/CLEAR,NOSTART\n/TITLE,",
        RATIO_BANNER_TMPL.format(
            kind=kind, parent_job=parent_job, parent_sec=parent_sec,
            ratio_label=ratio_label, omega="%.10f" % omega,
            hz="%.5f" % hz, vs_partner=vs_partner,
        ) + "FINISH\n/CLEAR,NOSTART\n/TITLE,",
        "%s banner insert" % stem_tag,
    )

    # each parent's own CFOPEN/VWRITE stems (oc_ vs scboth_)
    old_modes_cfopen = "*CFOPEN,%s_modes,txt" % cfopen_stem
    new_modes_cfopen = "*CFOPEN,%s_%s_modes,txt" % (cfopen_stem, stem_tag)
    text = once(text, old_modes_cfopen, new_modes_cfopen,
                "%s modes CFOPEN" % stem_tag)

    old_disc_cfopen = "*CFOPEN,%s_discriminator,txt" % cfopen_stem
    new_disc_cfopen = "*CFOPEN,%s_%s_discriminator,txt" % (cfopen_stem, stem_tag)
    text = once(text, old_disc_cfopen, new_disc_cfopen,
                "%s discriminator CFOPEN" % stem_tag)

    text = once(text, vwrite_tag,
                "('mode f[Hz] tag=%s_%s target=%.5fHz')"
                % (cfopen_stem.replace("ansys_p4_epsrel_ff_lanb_", "epsrel_ff_lanb_"),
                   stem_tag, hz),
                "%s modes VWRITE tag" % stem_tag)

    return text


def main() -> None:
    oc_parent_text = load(OC_PARENT)
    scboth_parent_text = load(SCBOTH_PARENT)

    outputs = {}
    for denom, ratio_label, stem_tag, oc_om, oc_hz, sc_om, sc_hz in RATIOS:
        oc_fname = "ansys_p4_epsrel_ff_lanb_%s_oc_2026-09-16.inp" % stem_tag
        outputs[oc_fname] = make_deck(
            oc_parent_text, "OPEN-CIRCUIT", "2491761", "18.176-18.177",
            denom, ratio_label, stem_tag, oc_om, oc_hz,
            "ABOVE this ratio's own scboth mode (stiffening, same mesh)",
            cfopen_stem="ansys_p4_epsrel_ff_lanb_oc",
            vwrite_tag="('mode f[Hz] tag=epsrel_ff_lanb_oc target=119.405Hz')",
            title_suffix="OPEN-CIRCUIT",
        )
        sc_fname = "ansys_p4_epsrel_ff_lanb_%s_scboth_2026-09-16.inp" % stem_tag
        outputs[sc_fname] = make_deck(
            scboth_parent_text, "BOTH-FACES SC CONTROL", "2491831", "18.178",
            denom, ratio_label, stem_tag, sc_om, sc_hz,
            "within mesh noise of this ratio's own sine-SC closed-form Hz",
            cfopen_stem="ansys_p4_epsrel_ff_lanb_scboth",
            vwrite_tag="('mode f[Hz] tag=epsrel_ff_lanb_scboth target=119.047Hz')",
            title_suffix="BOTH-FACES SC",
        )

    all_jobname_sources = EXISTING + list(outputs)
    print("jobname uniqueness (against full epsrel batch):")
    assert_unique_jobnames(all_jobname_sources)

    for fname, text in outputs.items():
        if "\r" in text:
            raise SystemExit("FATAL %s: CR leaked into output" % fname)
        if "2*H/12" in text:
            raise SystemExit("FATAL %s: leftover 1/12 H1 assignment" % fname)
        if "h1/2h=1/12, mm-N-s-tonne (2026-09-16)" in text:
            raise SystemExit("FATAL %s: leftover 1/12 title text" % fname)
        save(os.path.join(HERE, fname), text)
    print("OK: %d decks" % len(outputs))


if __name__ == "__main__":
    main()
