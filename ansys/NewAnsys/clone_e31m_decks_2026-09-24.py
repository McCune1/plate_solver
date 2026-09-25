# -*- coding: utf-8 -*-
"""
clone_e31m_decks_2026-09-24.py -- e31 sign migration (LESSONS Sec 18.244).

Every FE deck that backs a COUPLED Paper 4/5/6 claim ran with the project's
old E31 = +4.1 C/m^2. Decision 2026-09-24: adopt the PZT-4 constants of
Liu, Wang & Quek 2002 Table 1 exactly as printed (e31 = -4.1 C/m^2; Duan
2005 Table 1 prints the same). This script clones those decks with ONE
physical change, E31 = -4.1E-3, and renames every *CFOPEN output by
appending _e31m so no parent output is overwritten. Nothing else changes
(mesh, BCs, element, solver, sampling). Coupling-off (e0) and bare-host
decks are not cloned: they do not depend on e31.

Each clone gets a banner stating its parent. The parent's comment header
(which quotes +4.1 model targets) is kept verbatim for traceability; the
-4.1 model values are computed separately and scored after the run.

Run from Ansys/NewAnsys/:  python3 clone_e31m_decks_2026-09-24.py
"""
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
DATE = "2026-09-24"
PARENTS = [
    # Paper 4, axisymmetric PLANE223 (jobs 2490706, 2490848, 2491761, 2491831, 2494226)
    "ansys_p4_epsrel_ff_lanb_2026-09-16.inp",
    "ansys_p4_epsrel_ff_subsp_2026-09-16.inp",
    "ansys_p4_epsrel_cc_lanb_2026-09-16.inp",
    "ansys_p4_epsrel_ff_lanb_h18_2026-09-16.inp",
    "ansys_p4_epsrel_ff_lanb_h15_2026-09-16.inp",
    "ansys_p4_epsrel_ff_lanb_oc_2026-09-16.inp",
    "ansys_p4_epsrel_ff_lanb_scboth_2026-09-16.inp",
    "ansys_p4_epsrel_ff_lanb_h18_oc_2026-09-16.inp",
    "ansys_p4_epsrel_ff_lanb_h18_scboth_2026-09-16.inp",
    "ansys_p4_epsrel_ff_lanb_h15_oc_2026-09-16.inp",
    "ansys_p4_epsrel_ff_lanb_h15_scboth_2026-09-16.inp",
    # Paper 4 mesh refinement (job 2494665)
    "ansys_p4_zp4_ff_oc_h112_2026-09-16.inp", "ansys_p4_zp4_ff_sc_h112_2026-09-16.inp",
    "ansys_p4_zp4_ff_oc_h18_2026-09-16.inp", "ansys_p4_zp4_ff_sc_h18_2026-09-16.inp",
    "ansys_p4_zp4_ff_oc_h15_2026-09-16.inp", "ansys_p4_zp4_ff_sc_h15_2026-09-16.inp",
    # Paper 4 driven admittance (job 2494863) and C-F split (job 2495593)
    "ansys_p4_harm_ff_y_v2_h112_2026-09-16.inp",
    "ansys_p4_epsrel_cf_lanb_oc_2026-09-16.inp",
    "ansys_p4_epsrel_cf_lanb_scboth_2026-09-16.inp",
    # Paper 4 3-D n >= 1 check (job 2530423), coupled decks only
    "ansys_p4_3d_h112_sc_coarse_2026-09-23.inp", "ansys_p4_3d_h112_sc_fine_2026-09-23.inp",
    "ansys_p4_3d_h112_oc_fine_2026-09-23.inp",
    "ansys_p4_3d_h15_sc_coarse_2026-09-23.inp", "ansys_p4_3d_h15_sc_fine_2026-09-23.inp",
    "ansys_p4_3d_h15_oc_fine_2026-09-23.inp",
    # Paper 4 Liu 2002 solid disk (job 2531012)
    "ansys_p4_liudisk2_sc_C_2026-09-24.inp", "ansys_p4_liudisk2_sc_S_2026-09-24.inp",
    # Paper 5 F-F SC/OC (job 2520476)
    "ansys_p5_epsrel_ff_lanb_oc_2026-09-19.inp", "ansys_p5_epsrel_ff_lanb_scboth_2026-09-19.inp",
    # Paper 6 sector SC (job 2528994)
    "ansys_p6_sector_sc_coarse_2026-09-22.inp", "ansys_p6_sector_sc_fine_2026-09-22.inp",
]
E31_LINE = re.compile(r"^(\s*E31\s*=\s*)4\.1E-3(.*)$", re.M)
CFOPEN = re.compile(r"^(\*CFOPEN,)([A-Za-z0-9_]+)(,(?:txt|csv))", re.M | re.I)


def clone_name(parent):
    stem = re.sub(r"_2026-\d\d-\d\d\.inp$", "", parent)
    return "%s_e31m_%s.inp" % (stem, DATE)


def main():
    out = []
    for p in PARENTS:
        src = open(os.path.join(HERE, p)).read()
        n_e = len(E31_LINE.findall(src))
        if n_e != 1:
            raise SystemExit("%s: expected exactly one 'E31 = 4.1E-3' line, found %d" % (p, n_e))
        new = E31_LINE.sub(r"\g<1>-4.1E-3\g<2>    ! e31 SIGN CLONE: was +4.1E-3", src)
        new, n_c = CFOPEN.subn(r"\g<1>\g<2>_e31m\g<3>", new)
        if n_c == 0:
            raise SystemExit("%s: no *CFOPEN outputs found" % p)
        # Four parents (h18/h15 oc/scboth) quote the MAPDL failure token in a
        # comment, which made run_ansys_queue.sh false-FAIL them in job 2494226
        # (Sec 18.185). Neutralise it in comment lines only.
        lines = new.split("\n")
        for i, ln in enumerate(lines):
            if ln.lstrip().startswith("!"):
                lines[i] = ln.replace("*** ERROR ***", "[MAPDL error token]")
        new = "\n".join(lines)
        name = clone_name(p)
        banner = ("! %s\n! e31 SIGN CLONE (%s, clone_e31m_decks_2026-09-24.py, LESSONS Sec 18.244).\n"
                  "! Parent: %s. Only changes: E31 = -4.1E-3 (Liu 2002 / Duan 2005\n"
                  "! Table 1 as printed) and every *CFOPEN output name + '_e31m'.\n"
                  "! Targets quoted in the parent header below are the OLD +4.1 values.\n"
                  "! Do not put the MAPDL failure token in this file.\n" % (name, DATE, p))
        new = banner + new
        assert "*** ERROR" not in new
        with open(os.path.join(HERE, name), "w", newline="\n") as fh:
            fh.write(new)
        out.append((name, n_c))
    man = os.path.join(HERE, "ansys_queue_manifest_e31m_clones_%s.txt" % DATE)
    with open(man, "w", newline="\n") as fh:
        fh.write("# e31 = -4.1 clones of every coupled Paper 4/5/6 FE deck (%s).\n"
                 "# Built by clone_e31m_decks_2026-09-24.py. ONE licence seat.\n#\n" % DATE)
        for n, _ in out:
            fh.write(n + "\n")
    for n, c in out:
        print("%-55s CFOPEN renamed: %d" % (n, c))
    print("wrote %d clones and %s" % (len(out), os.path.basename(man)))


if __name__ == "__main__":
    main()
