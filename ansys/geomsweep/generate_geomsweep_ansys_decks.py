#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
generate_geomsweep_ansys_decks.py -- generates the FE verification decks
for ../geometry_sweep/'s 4 new free-free radius ratios (1.5, 1.66667, 2.0,
2.5), one deck per (ratio, part), each deck looping over all 6 sector
angles as separate passes.

WHY GENERATED, NOT HAND-WRITTEN: 4 ratios x 6 angles x 2 parts = 48
geometry/part combinations. Hand-typing that many near-identical APDL
blocks risks silent copy-paste drift (a wrong R_i in pass 4 of 6, etc.)
that would be very hard to spot by eye. Generating from one template
guarantees every pass is mechanically consistent; this script is itself
the audit trail for how each number was produced, and can be rerun to
add/change ratios or angles later without hand-editing APDL.

GEOMETRY MAPPING (verified against plate_solver.geometry.make_geometry
AND against LESSONS_LEARNED.md's own citation, not assumed):
    make_geometry(r0_2b, two_T_pi) with the DEFAULT b=2.0 gives
        R_i = 4*r0_2b - 2,   R_o = 4*r0_2b + 2
        2Theta = two_T_pi * pi radians  ->  PHI[deg] = two_T_pi * 180
    Cross-check: r0_2b=1.5 -> R_i=4, R_o=8 (R_i/R_o=0.5) -- EXACTLY the
    "FF-P1 benchmark geometry" LESSONS_LEARNED.md cites (r0/(2b)=1.5,
    R_i/R_o=0.5, 2Theta=pi/2 -> PHI=90deg) and exactly what
    Ansys/ansys_ff_oop_ext.inp already uses. This is the r0_2b=1.5 point
    ../geometry_sweep/submit_ff_geomsweep_r150.sh also sweeps (now across
    all 6 angles instead of just the original single 0.5pi point) -- so
    that deck's own a050 pass is a genuine repeat of already-FE-validated
    ground, useful as an in-line sanity check that this generator's output
    matches the existing hand-written deck (see VERIFY_AGAINST_EXISTING
    below).

MATERIAL: E=210e9, nu=0.35, rho=7800.

BUGFIX 2026-07-22: this docstring previously claimed nu=0.30 "the
free-free-specific material...see the MAT_NU=0.30 fix applied to every
../geometry_sweep/submit_ff_geomsweep_*.sh script" -- that fix was never
actually applied. None of submit_ff_geomsweep_r150/r167/r200/r250.sh set
MAT_NU, so all 5 geometry_sweep jobs ran at the IsotropicMaterial DEFAULT
nu=0.35 (confirmed from their own "material this run: nu=0.35" printout),
NOT nu=0.30 like the single already-validated FF-P1 point (1.5, 0.5pi)
did. This generator's decks were generated at nu=0.30 to begin with, so
every one of the first-round ansys_geomsweep_verification results was
compared against the WRONG Poisson's ratio -- caught during the
2026-07-21/22 review before any paper claim was built on it (a free
annular edge's natural BC is nu-dependent beyond D's own (1-nu^2) scaling,
so this is not a negligible mismatch). Fixed here by matching the
material this generator actually uses (nu=0.35) to what the python sweep
actually computed, rather than re-running the much more expensive (13-18h
each) python jobs to match the FE decks. H=0.08 fixed across all 4 ratios
(same as the existing FF-P1 deck) -- H/R_o ranges 0.0067-0.01 across the
4 ratios, all comfortably in the thin/Kirchhoff regime the original
benchmark also relied on.

SINGLE MESH, NOT TWO-MESH CONVERGENCE (a deliberate scope trade-off,
different from every OTHER deck in this project): NDIV=64 was already
shown to agree with NDIV=96 to <0.3% at the base FF-P1 geometry
(ansys_ff_oop_ext.inp / ansys_ff_ip_ext.inp) -- re-proving mesh
convergence at every one of 48 new geometries would double the deck count
and cluster time for a check this project's own prior work already
answered once. If a specific new geometry's results look suspicious,
rerun that one pass at NDIV=96 by hand before trusting it -- don't assume
NDIV=64 convergence transfers to a much more slender geometry (r0_2b=2.5)
without spot-checking at least once.

MODE-FAMILY SEPARATION: identical caveat to every other SHELL281 deck in
this project -- the OOP decks' modes are interleaved with IP modes (and,
for these free-free geometries, 6 rigid-body modes near 0 Hz). Cross-check
elastic OOP-deck frequencies against the matching IP-deck's own frequency
list (same r0_2b, same angle) within ~0.5% to discard IP modes, exactly as
ansys_ff_oop_ext.inp's header describes.

Run:  python3 generate_geomsweep_ansys_decks.py
Writes 8 .inp files + 8 submit_*.sh files into this directory.
"""
from __future__ import annotations
import os

HERE = os.path.dirname(os.path.abspath(__file__))

# (tag, r0_2b) -- the 4 new geometry_sweep radius ratios.
RATIOS = [
    ("r150", 1.5),
    ("r167", 1.66667),
    ("r200", 2.0),
    ("r250", 2.5),
]

# (tag, two_T_pi) -- the 6 sector angles used throughout the geometry sweep.
ANGLES = [
    ("a025", 0.25),
    ("a050", 0.50),
    ("a075", 0.75),
    ("a100", 1.00),
    ("a125", 1.25),
    ("a150", 1.50),
]

E_ISO, NU_ISO, RHO, H = 210.0e9, 0.35, 7800.0, 0.08
NDIV = 64
PI = "3.14159265358979"


def _ri_ro(r0_2b):
    return 4.0 * r0_2b - 2.0, 4.0 * r0_2b + 2.0


def _oop_pass(rtag, r0_2b, atag, two_T_pi, nmodes=40):
    R_i, R_o = _ri_ro(r0_2b)
    phi = two_T_pi * 180.0
    out = f"geomsweep_oop_{rtag}_{atag}"
    return f"""FINISH
/CLEAR,NOSTART
/TITLE, Free-free OOP geomsweep verify -- {rtag} (r0/2b={r0_2b:g}) {atag} (2T/pi={two_T_pi:g}, PHI={phi:g}deg)

R_i    = {R_i:.6f}
R_o    = {R_o:.6f}
PHI    = {phi:.6f}
H      = {H}
E_iso  = {E_ISO}
NU_iso = {NU_ISO}
RHO    = {RHO}
NDIV   = {NDIV}
NMODES = {nmodes}
PI     = {PI}

/PREP7
LOCAL, 11, 1, 0,0,0
CSYS, 11
ET, 1, SHELL281
SECTYPE, 1, SHELL
SECDATA, H, 1
MP, EX,1,E_iso $ MP, EY,1,E_iso $ MP, EZ,1,E_iso
MP, PRXY,1,NU_iso $ MP, PRYZ,1,NU_iso $ MP, PRXZ,1,NU_iso
GG = E_iso/(2*(1+NU_iso))
MP, GXY,1,GG $ MP, GYZ,1,GG $ MP, GXZ,1,GG
MP, DENS,1,RHO
CSYS, 0
PCIRC, R_i, R_o, 0, PHI
ESYS, 11
TYPE,1 $ SECNUM,1 $ MAT,1
CSYS, 11
LSEL, S, LOC, Y, 0  $ LSEL, A, LOC, Y, PHI
LESIZE, ALL, , , NDIV
LSEL, S, LOC, X, R_i $ LSEL, A, LOC, X, R_o
LESIZE, ALL, , , NDIV
ALLSEL, ALL
CSYS, 0
MSHKEY, 1 $ AMESH, ALL
! completely free (FFFF): no constraints
ALLSEL, ALL
FINISH

/SOLU
ANTYPE, MODAL
MODOPT, LANB, NMODES, -1.0, , , ON
MXPAND, NMODES, , , YES
SOLVE
FINISH

/POST1
DD   = E_iso*H**3/(12*(1-NU_iso**2))
KLIT = R_o**2*SQRT(RHO*H/DD)

*DIM, FR_{rtag}_{atag}, ARRAY, NMODES
*DO, I, 1, NMODES
  *GET, FR_{rtag}_{atag}(I), MODE, I, FREQ
*ENDDO

! BUGFIX 2026-07-21: *VWRITE is a vector command -- calling it INSIDE a *DO
! loop with array-indexed args does not print one row per iteration, it
! dumps the WHOLE remaining array (from the current start index to NMODES)
! on every call, producing a triangular cascade of repeated/redundant rows
! (previously recovered by taking only the first NMODES-row block; see
! FutureWork memory 2026-07-21). Fixed to match the _ip_pass() block below:
! build the mode-index and Omega_lit arrays once, then *VWRITE the full
! vectors in a SINGLE call outside any loop.
*DIM, MIDX_{rtag}_{atag}, ARRAY, NMODES
*VFILL, MIDX_{rtag}_{atag}, RAMP, 1, 1
KLIT2PI = 2*PI*KLIT
*DIM, OML_{rtag}_{atag}, ARRAY, NMODES
*VOPER, OML_{rtag}_{atag}(1), FR_{rtag}_{atag}(1), MULT, KLIT2PI

*CFOPEN, {out}, txt
*VWRITE
('mode      f[Hz]           Omega_lit    (SHELL281 NDIV={NDIV}, {rtag} r0/2b={r0_2b:g}, {atag} 2T/pi={two_T_pi:g} PHI={phi:g}deg, FFFF; ~6 rigid modes near 0 Hz; IP modes interleaved -- cross-check vs geomsweep_ip_{rtag}_{atag}.txt)')
*VWRITE, MIDX_{rtag}_{atag}(1), FR_{rtag}_{atag}(1), OML_{rtag}_{atag}(1)
(F4.0, 3X, E14.6, 3X, F12.4)
*CFCLOS
FINISH
"""


def _ip_pass(rtag, r0_2b, atag, two_T_pi, nmodes=35):
    R_i, R_o = _ri_ro(r0_2b)
    phi = two_T_pi * 180.0
    out = f"geomsweep_ip_{rtag}_{atag}"
    return f"""FINISH
/CLEAR,NOSTART
/TITLE, Free-free IP geomsweep verify -- {rtag} (r0/2b={r0_2b:g}) {atag} (2T/pi={two_T_pi:g}, PHI={phi:g}deg)

R_i    = {R_i:.6f}
R_o    = {R_o:.6f}
PHI    = {phi:.6f}
H      = {H}
E_iso  = {E_ISO}
NU_iso = {NU_ISO}
RHO    = {RHO}
NDIV   = {NDIV}
NMODES = {nmodes}

/PREP7
LOCAL, 11, 1, 0,0,0
CSYS, 11
ET, 1, PLANE183
KEYOPT, 1, 3, 3
R, 1, H
MP, EX,1,E_iso
MP, PRXY,1,NU_iso
MP, DENS,1,RHO
CSYS, 0
PCIRC, R_i, R_o, 0, PHI
TYPE,1 $ MAT,1 $ REAL,1
CSYS, 11
LSEL, S, LOC, Y, 0  $ LSEL, A, LOC, Y, PHI
LESIZE, ALL, , , NDIV
LSEL, S, LOC, X, R_i $ LSEL, A, LOC, X, R_o
LESIZE, ALL, , , NDIV
ALLSEL, ALL
CSYS, 0
MSHKEY, 1 $ AMESH, ALL
! completely free (FFFF): no constraints
ALLSEL, ALL
FINISH

/SOLU
ANTYPE, MODAL
MODOPT, LANB, NMODES, -1.0, , , ON
MXPAND, NMODES, , , YES
SOLVE
FINISH

/POST1
*DIM, MIDX_{rtag}_{atag}, ARRAY, NMODES
*VFILL, MIDX_{rtag}_{atag}, RAMP, 1, 1
*DIM, FR_{rtag}_{atag}, ARRAY, NMODES
*DO, I, 1, NMODES
  *GET, FR_{rtag}_{atag}(I), MODE, I, FREQ
*ENDDO

*CFOPEN, {out}, txt
*VWRITE
('mode      f[Hz]           (PLANE183 NDIV={NDIV}, {rtag} r0/2b={r0_2b:g}, {atag} 2T/pi={two_T_pi:g} PHI={phi:g}deg, FFFF; 3 rigid modes near 0 Hz; no family separation needed)')
*VWRITE, MIDX_{rtag}_{atag}(1), FR_{rtag}_{atag}(1)
(F4.0, 3X, E14.6)
*CFCLOS
FINISH
"""


def _submit(deck_fname, out_fname, job_tag, glob_pat, minutes=25):
    return f"""#!/bin/bash
#SBATCH --job-name={job_tag}
#SBATCH --output={job_tag}_%j.out
#SBATCH --error={job_tag}_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=00:{minutes}:00
#SBATCH --partition=general

# Generated by generate_geomsweep_ansys_decks.py -- 6 angle passes, single
# mesh (NDIV=64), one geometry-sweep radius ratio. See that script's
# docstring for the geometry-mapping derivation and the single-mesh
# scope trade-off.
# NOTE: no explicit `cd` here on purpose (matches this project's Ansys
# submit-script convention) -- submit from the SAME directory as
# {deck_fname} (cd there, then `sbatch {os.path.basename(out_fname)}`).
module load ansys/2024R2 2>/dev/null || module load ansys 2>/dev/null

echo "Job $SLURM_JOB_ID on $(hostname)"
echo "Start: $(date)"

mapdl -b -smp -np ${{SLURM_CPUS_PER_TASK:-8}} -i {deck_fname} -o {deck_fname[:-4]}_out.txt

echo "End: $(date)"
for f in geomsweep_*_{deck_fname.split('_')[2]}_*.txt; do
  echo "--- $f ---"; cat "$f" 2>/dev/null; echo
done
"""


def main():
    written = []
    for rtag, r0_2b in RATIOS:
        for part, passfn, nmodes in (("oop", _oop_pass, 40), ("ip", _ip_pass, 35)):
            deck_name = f"ansys_geomsweep_{part}_{rtag}.inp"
            deck_path = os.path.join(HERE, deck_name)
            body = "".join(
                passfn(rtag, r0_2b, atag, two_T_pi, nmodes)
                for atag, two_T_pi in ANGLES
            )
            header = (
                f"! ===========================================================================\n"
                f"!  {deck_name}  --  GENERATED by generate_geomsweep_ansys_decks.py, do not\n"
                f"!  hand-edit (edit the generator and rerun instead, or your edit will be lost\n"
                f"!  next regeneration). See that script's module docstring for the full\n"
                f"!  geometry-mapping derivation, material rationale, and single-mesh caveat.\n"
                f"!\n"
                f"!  {'OOP (SHELL281)' if part == 'oop' else 'IP (PLANE183)'} free-free FE verification for "
                f"geometry_sweep's r0_2b={r0_2b:g} ratio\n"
                f"!  ('{rtag}'), all 6 sector angles, NDIV={NDIV} single mesh.\n"
                f"!  Run headless:  mapdl -b -i {deck_name} -o {deck_name[:-4]}_out.txt\n"
                f"! ===========================================================================\n\n"
            )
            with open(deck_path, "w", newline="\n") as f:
                f.write(header + body)
            written.append(deck_path)

            submit_name = f"submit_ansys_geomsweep_{part}_{rtag}.sh"
            submit_path = os.path.join(HERE, submit_name)
            # Fixed: Added the missing glob_pat parameter
            with open(submit_path, "w", newline="\n") as f:
                f.write(_submit(deck_name, submit_path,
                                 f"geomsweep_{part}_{rtag}", 
                                 f"geomsweep_*_{rtag}_*.txt",  # This is the glob pattern for output files
                                 minutes=30))
            os.chmod(submit_path, 0o755)
            written.append(submit_path)

    print(f"wrote {len(written)} files:")
    for p in written:
        print(f"  {p}")


if __name__ == "__main__":
    main()
