#!/usr/bin/env python3
"""Lever 8 (PAPER4_LEVERS_2026-09-16.md sec B.8): mesh refinement of the
fair OC / both-faces-SC pair. Isolated variable: NDIV_ZP 2 -> 4
(through-thickness divisions of each piezo layer). Units, TB,PIEZ,
LANB, electrical BC, NDIV_R, NDIV_ZH are unchanged.

Why this variable: the FE/analytic OC-split recovery fraction declines
as the layer thickens (93% at 1/12, 90% at 1/8, 86% at 1/5, jobs
2491761/2491831/2494226). NDIV_ZP is currently fixed at 2 and does not
scale with h1. One doubling at 1/12 first (tightest gap); 1/8 and 1/5
are included so the manuscript's "plausible but unconfirmed" trend can
be re-read on the refined mesh.

Clones from the confirmed 1/12 OC (job 2491761) and scboth (job
2491831) parents -- NOT from the lever-7 ratio clones, whose banners
contain the literal MAPDL error-token string that false-positived
job 2494226's queue scan. New decks must not contain that token.

Closed-form targets (unchanged from Sec 18.182):

  ratio | OC Hz     | scboth Hz
  1/12  | 119.40431 | 119.04677
  1/8   | 122.15932 | 121.61564
  1/5   | 128.26320 | 127.38344

Coarse-mesh FE (for comparison, not targets):
  1/12  119.143 vs 118.811  (+0.280% vs analytic +0.301%, 93%)
  1/8   121.839 vs 121.348  (+0.405% vs analytic +0.450%, 90%)
  1/5   127.809 vs 127.047  (+0.600% vs analytic +0.700%, 86%)

PRE-REGISTERED reading of the refined pair, per ratio:
  - 0 Hz rigid-body present, bending discriminator real, 0 MAPDL error
    lines, RUN COMPLETED (read NUMBER OF ERROR MESSAGES ENCOUNTERED,
    never the queue script's own PASS/FAIL label).
  - OC mode 2 ABOVE same-mesh scboth mode 2 (stiffening still there).
  - Scientific question, reported not gated: does the FE/analytic split
    gap shrink (mesh was the cause), stay (ansatz / something else),
    or grow (stop and report)?
"""
from __future__ import annotations

import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))

OC_PARENT = os.path.join(HERE, "ansys_p4_epsrel_ff_lanb_oc_2026-09-16.inp")
SCBOTH_PARENT = os.path.join(HERE, "ansys_p4_epsrel_ff_lanb_scboth_2026-09-16.inp")

# (denom_or_None, ratio_label, stem_tag, oc_hz, sc_hz)
# denom None keeps the parent's H1 = 2*H/12 assignment.
RATIOS = [
    (None, "1/12", "h112", 119.40431, 119.04677),
    (8,    "1/8",  "h18",  122.15932, 121.61564),
    (5,    "1/5",  "h15",  128.26320, 127.38344),
]

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
    "ansys_p4_epsrel_ff_lanb_h18_oc_2026-09-16.inp",
    "ansys_p4_epsrel_ff_lanb_h18_scboth_2026-09-16.inp",
    "ansys_p4_epsrel_ff_lanb_h15_oc_2026-09-16.inp",
    "ansys_p4_epsrel_ff_lanb_h15_scboth_2026-09-16.inp",
]

FORBIDDEN = "*** ERROR ***"


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


BANNER = r"""!
! ============================================================
! LEVER 8 MESH REFINEMENT (2026-09-16): cloned from the h1/2h=1/12
! {kind} parent (job {parent_job}). Isolated change: NDIV_ZP 2 -> 4.
! Units, TB,PIEZ, LANB, electrical BC, NDIV_R=40, NDIV_ZH=6 are
! UNCHANGED. Do not re-open those. Do not confuse this with the
! outer-only-ground h18/h15 SC decks (job 2490848).
!
! THIS DECK: h1/2h={ratio_label}, {kind}, NDIV_ZP=4. Closed-form
! target {hz:.5f} Hz. PASS: bending mode present, 0 Hz rigid-body
! present, real O(0.1-1) bending discriminator, 0 MAPDL error
! lines, RUN COMPLETED. OC sibling must sit ABOVE same-mesh scboth.
! Scientific question (reported, not a run-fail): does the FE vs
! closed-form split gap shrink vs the NDIV_ZP=2 pair?
! ============================================================
"""


def make_deck(parent_text: str, kind: str, parent_job: str,
              denom, ratio_label: str, stem_tag: str, hz: float,
              cfopen_stem: str, vwrite_tag: str) -> str:
    text = parent_text
    text = once(
        text,
        "NDIV_ZP = 2              ! through-thickness divisions, each piezo layer",
        "NDIV_ZP = 4              ! lever 8: doubled from 2; isolate this one mesh variable",
        "%s NDIV_ZP" % stem_tag,
    )
    if denom is not None:
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
            "h1/2h=%s NDIV_ZP=4, mm-N-s-tonne (2026-09-16)" % ratio_label,
            "%s TITLE ratio" % stem_tag,
        )
    else:
        text = once(
            text,
            "h1/2h=1/12, mm-N-s-tonne (2026-09-16)",
            "h1/2h=1/12 NDIV_ZP=4, mm-N-s-tonne (2026-09-16)",
            "%s TITLE zp4" % stem_tag,
        )
    text = once(
        text,
        "FINISH\n/CLEAR,NOSTART\n/TITLE,",
        BANNER.format(kind=kind, parent_job=parent_job,
                      ratio_label=ratio_label, hz=hz)
        + "FINISH\n/CLEAR,NOSTART\n/TITLE,",
        "%s banner insert" % stem_tag,
    )
    new_modes = "*CFOPEN,ansys_p4_zp4_ff_%s_%s_modes,txt" % (
        "oc" if "OPEN" in kind else "sc", stem_tag)
    new_disc = "*CFOPEN,ansys_p4_zp4_ff_%s_%s_discriminator,txt" % (
        "oc" if "OPEN" in kind else "sc", stem_tag)
    text = once(text, "*CFOPEN,%s_modes,txt" % cfopen_stem, new_modes,
                "%s modes CFOPEN" % stem_tag)
    text = once(text, "*CFOPEN,%s_discriminator,txt" % cfopen_stem, new_disc,
                "%s discriminator CFOPEN" % stem_tag)
    text = once(
        text, vwrite_tag,
        "('mode f[Hz] tag=zp4_%s_%s target=%.5fHz')"
        % ("oc" if "OPEN" in kind else "sc", stem_tag, hz),
        "%s modes VWRITE tag" % stem_tag,
    )
    return text


def main() -> None:
    oc_parent = load(OC_PARENT)
    sc_parent = load(SCBOTH_PARENT)
    outputs = {}
    for denom, ratio_label, stem_tag, oc_hz, sc_hz in RATIOS:
        oc_fname = "ansys_p4_zp4_ff_oc_%s_2026-09-16.inp" % stem_tag
        outputs[oc_fname] = make_deck(
            oc_parent, "OPEN-CIRCUIT", "2491761",
            denom, ratio_label, stem_tag, oc_hz,
            cfopen_stem="ansys_p4_epsrel_ff_lanb_oc",
            vwrite_tag="('mode f[Hz] tag=epsrel_ff_lanb_oc target=119.405Hz')",
        )
        sc_fname = "ansys_p4_zp4_ff_sc_%s_2026-09-16.inp" % stem_tag
        outputs[sc_fname] = make_deck(
            sc_parent, "BOTH-FACES SC CONTROL", "2491831",
            denom, ratio_label, stem_tag, sc_hz,
            cfopen_stem="ansys_p4_epsrel_ff_lanb_scboth",
            vwrite_tag="('mode f[Hz] tag=epsrel_ff_lanb_scboth target=119.047Hz')",
        )

    all_jobname_sources = EXISTING + list(outputs)
    print("jobname uniqueness (against full epsrel batch + new decks):")
    assert_unique_jobnames(all_jobname_sources)

    for fname, text in outputs.items():
        if "\r" in text:
            raise SystemExit("FATAL %s: CR leaked into output" % fname)
        if FORBIDDEN in text:
            raise SystemExit(
                "FATAL %s: contains the MAPDL error-token string that "
                "false-positives run_ansys_queue.sh" % fname)
        if "NDIV_ZP = 2              ! through-thickness divisions, each piezo layer" in text:
            raise SystemExit("FATAL %s: leftover NDIV_ZP=2 assignment" % fname)
        save(os.path.join(HERE, fname), text)
    print("OK: %d decks" % len(outputs))


if __name__ == "__main__":
    main()
