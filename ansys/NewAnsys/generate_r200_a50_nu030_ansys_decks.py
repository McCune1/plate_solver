#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
generate_r200_a50_nu030_ansys_decks.py -- rank 16 FE decks.

r0/2b=2.0, 2Theta/pi=0.5 (PHI=90deg), nu=0.30. Crossed geometry:
r200 ratio × FF-P1 angle. Writes into Ansys/NewAnsys/. Unique MAPDL
-j names (r200a50oop / r200a50ip) so leftover file.lock cannot kill
the job the way 2421046/2421049 died.

IP NMODES=80 (avoid the r150 35-mode ceiling). OOP NMODES=40.
"""
import os

HERE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                    "Ansys", "NewAnsys")

RTAG, R0_2B = "r200", 2.0
ATAG, TWO_T_PI = "a50", 0.5
E_ISO, NU_ISO, RHO, H = 210.0e9, 0.30, 7800.0, 0.08
NDIV = 64
PI = "3.14159265358979"


def _ri_ro(r0_2b):
    return 4.0 * r0_2b - 2.0, 4.0 * r0_2b + 2.0


def _oop_pass(rtag, r0_2b, atag, two_T_pi, nmodes=40):
    R_i, R_o = _ri_ro(r0_2b)
    phi = two_T_pi * 180.0
    out = f"geomsweep_oop_{rtag}_{atag}_nu030"
    return f"""FINISH
/CLEAR,NOSTART
/TITLE, Free-free OOP verify (nu=0.30) -- {rtag} (r0/2b={r0_2b:g}) {atag} (2T/pi={two_T_pi:g}, PHI={phi:g}deg)

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

*DIM, MIDX_{rtag}_{atag}, ARRAY, NMODES
*VFILL, MIDX_{rtag}_{atag}, RAMP, 1, 1
KLIT2PI = 2*PI*KLIT
*DIM, OML_{rtag}_{atag}, ARRAY, NMODES
*VOPER, OML_{rtag}_{atag}(1), FR_{rtag}_{atag}(1), MULT, KLIT2PI

*CFOPEN, {out}, txt
*VWRITE
('mode      f[Hz]           Omega_lit    (SHELL281 NDIV={NDIV}, {rtag} r0/2b={r0_2b:g}, {atag} 2T/pi={two_T_pi:g} PHI={phi:g}deg, nu={NU_ISO:g}, FFFF; ~6 rigid modes near 0 Hz)')
*VWRITE, MIDX_{rtag}_{atag}(1), FR_{rtag}_{atag}(1), OML_{rtag}_{atag}(1)
(F4.0, 3X, E14.6, 3X, F12.4)
*CFCLOS
FINISH
"""


def _ip_pass(rtag, r0_2b, atag, two_T_pi, nmodes=80):
    R_i, R_o = _ri_ro(r0_2b)
    phi = two_T_pi * 180.0
    out = f"geomsweep_ip_{rtag}_{atag}_nu030"
    return f"""FINISH
/CLEAR,NOSTART
/TITLE, Free-free IP verify (nu=0.30) -- {rtag} (r0/2b={r0_2b:g}) {atag} (2T/pi={two_T_pi:g}, PHI={phi:g}deg)

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
('mode      f[Hz]           (PLANE183 NDIV={NDIV}, {rtag} r0/2b={r0_2b:g}, {atag} 2T/pi={two_T_pi:g} PHI={phi:g}deg, nu={NU_ISO:g}, FFFF; 3 rigid modes near 0 Hz)')
*VWRITE, MIDX_{rtag}_{atag}(1), FR_{rtag}_{atag}(1)
(F4.0, 3X, E14.6)
*CFCLOS
FINISH
"""


def _submit(deck_fname, job_tag, mapdl_j, out_glob):
    return f"""#!/bin/bash
#SBATCH --job-name={job_tag}
#SBATCH --output={job_tag}_%j.out
#SBATCH --error={job_tag}_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=01:00:00
#SBATCH --requeue
#SBATCH --partition=requeue

# Rank 16 FE. Submit from Ansys/NewAnsys/.
# -j {mapdl_j} avoids leftover NewAnsys/file.lock (jobs 2421046/2421049).
module load ansys/2024R2 2>/dev/null || module load ansys 2>/dev/null

echo "Job $SLURM_JOB_ID on $(hostname)"
echo "Start: $(date)"

mapdl -b -smp -np ${{SLURM_CPUS_PER_TASK:-8}} -j {mapdl_j} -i {deck_fname} -o {deck_fname[:-4]}_out.txt

echo "End: $(date)"
for f in {out_glob}; do
  echo "--- $f ---"; cat "$f" 2>/dev/null; echo
done
"""


def main():
    os.makedirs(HERE, exist_ok=True)
    written = []
    specs = (
        ("oop", _oop_pass, 40, "r200a50oop",
         "geomsweep_oop_r200_a50_nu030.txt"),
        ("ip", _ip_pass, 80, "r200a50ip",
         "geomsweep_ip_r200_a50_nu030.txt"),
    )
    for part, passfn, nmodes, mapdl_j, out_glob in specs:
        deck_name = f"ansys_{part}_r200_a50_nu030.inp"
        deck_path = os.path.join(HERE, deck_name)
        body = passfn(RTAG, R0_2B, ATAG, TWO_T_PI, nmodes)
        header = (
            f"! Generated by generate_r200_a50_nu030_ansys_decks.py\n"
            f"! {part.upper()} FFFF  r0/2b=2.0  2T/pi=0.5  nu=0.30  "
            f"NMODES={nmodes}\n\n"
        )
        with open(deck_path, "w", newline="\n") as f:
            f.write(header + body)
        written.append(deck_path)

        submit_name = f"submit_ansys_{part}_r200_a50_nu030.sh"
        submit_path = os.path.join(HERE, submit_name)
        with open(submit_path, "w", newline="\n") as f:
            f.write(_submit(deck_name, f"ansys_{part}_r200_a50",
                            mapdl_j, out_glob))
        os.chmod(submit_path, 0o755)
        written.append(submit_path)

    print(f"wrote {len(written)} files:")
    for p in written:
        print(f"  {p}")
    print("\nSubmit from Ansys/NewAnsys:")
    print("  sbatch submit_ansys_oop_r200_a50_nu030.sh")
    print("  sbatch submit_ansys_ip_r200_a50_nu030.sh")


if __name__ == "__main__":
    main()
