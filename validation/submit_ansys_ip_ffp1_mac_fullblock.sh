#!/bin/bash
#SBATCH --job-name=ffp1_ip_mac
#SBATCH --output=ansys_ip_ffp1_mac_fullblock_%j.out
#SBATCH --error=ansys_ip_ffp1_mac_fullblock_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=03:00:00
#SBATCH --partition=requeue
#SBATCH --requeue

# FF-P1 (r0/2b=1.5, 2Theta/pi=0.5, nu=0.30) IP eigenvector dump.
# Distinct prefix ip_mac_ffp1_nu030 so it cannot collide with the
# Rank 14 nu=0.35 dump ip_mac_gs_r150_a050_*.
# Prerequisite for submit_ip_ffp1_1115_pair_mac_v1.sh and
# submit_ip_qin_pair_mac_v1.sh (--dependency=afterok:$this).

export RATIO_TAG=r150
export ANGLE_TAG=a050
export GS_NU=0.30
export FB_PREFIX=ip_mac_ffp1_nu030
export FB_DECK=ansys_ip_ffp1_nu030_mac_fullblock.inp
export FB_JOB=ffp1ipmac
export PKG_PATH=/home/ghmkfh/PythonMill/Plate_Solver_Package

module load ansys/2024R2 2>/dev/null || module load ansys 2>/dev/null
module load python/3.12.1
source /home/ghmkfh/PythonMill/Plate_Solver_Package/venv/bin/activate

cd /home/ghmkfh/PythonMill/Plate_Solver_Package/ || {
  echo "FATAL: cannot cd to package root."
  exit 2
}

if [ ! -f generate_ansys_ip_geomsweep_mac_fullblock.py ]; then
  echo "FATAL: generator not in $(pwd)."
  exit 2
fi

python3 -u generate_ansys_ip_geomsweep_mac_fullblock.py || {
  echo "FATAL: generator failed."
  exit 2
}

cd /home/ghmkfh/PythonMill/Plate_Solver_Package/Ansys/NewAnsys || {
  echo "FATAL: cannot cd to Ansys/NewAnsys."
  exit 2
}

DECK="$FB_DECK"
PREFIX="$FB_PREFIX"
if [ ! -f "$DECK" ]; then
  echo "FATAL: $DECK not in $(pwd) after generate."
  exit 2
fi

echo "Job $SLURM_JOB_ID on $(hostname)"
echo "Working directory: $(pwd)"
echo "GS_NU=$GS_NU PREFIX=$PREFIX JOB=$FB_JOB"
echo "Start: $(date)"

mapdl -b -smp -np ${SLURM_CPUS_PER_TASK:-8} -j "$FB_JOB" -i "$DECK" \
      -o "ansys_ip_ffp1_nu030_mac_fullblock_out.txt" || {
  echo "MAPDL_FAILED (exit $?) -- stopping before the summary greps."
  exit 1
}

OUTTXT="ansys_ip_ffp1_nu030_mac_fullblock_out.txt"
if head -40 "$OUTTXT" 2>/dev/null | grep -q -i "no input file specified"; then
  echo "FATAL: MAPDL never read the deck."
  exit 1
fi

echo "--- frequency table ---"
ls -la "${PREFIX}_freqs.txt" 2>&1
head -8 "${PREFIX}_freqs.txt" 2>/dev/null

echo "--- eigenvector file count (57 expected) ---"
N=$(ls -1 ${PREFIX}_m*_mesh96.txt 2>/dev/null | wc -l)
echo "found $N"
if [ "$N" -lt 57 ]; then
  echo "WARNING: only $N/57 mode files."
fi

echo "--- error/warning scan ---"
grep -n -i -E '\*\*\* ERROR \*\*\*|NUMBER OF ERROR MESSAGES|NUMBER OF WARNING MESSAGES|FATAL' \
  "$OUTTXT" || echo "(no error markers matched)"

echo "End: $(date)"
