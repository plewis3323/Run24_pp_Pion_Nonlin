#!/bin/bash

# Environment Setup
export USER="$(id -u -n)"
export LOGNAME="plewis3323"
export HOME="/sphenix/u/plewis3323"

# Define the absolute path of this script
this_script=$(readlink -f "${BASH_SOURCE}")
this_dir=$(dirname "${this_script}")

# Source sPHENIX environment
source /opt/sphenix/core/bin/sphenix_setup.sh -n
export MYINSTALL="/sphenix/user/plewis3323/install"
source /opt/sphenix/core/bin/setup_local.sh "${MYINSTALL}"

# Update LD_LIBRARY_PATH
export LD_LIBRARY_PATH="${MYINSTALL}/lib:${LD_LIBRARY_PATH}"

# Timestamp for job start
echo " "
echo "START: $(date)"
echo " "

#Check Condor scratch directory
if [[ -n "$_CONDOR_SCRATCH_DIR" && -d "$_CONDOR_SCRATCH_DIR" ]]; then
  cd "$_CONDOR_SCRATCH_DIR"
   # Explicitly copy macro file to Condor scratch
   cp /sphenix/tg/tg01/bulk/plewis3323/Run24_AuAu_New_2024_p007_DST_Calo/D_Calo_run2auau/Condor_process/*.C .
else
   echo "Condor scratch directory NOT set or inaccessible"
   exit 1
fi

# Debugging output
pwd
ls -l *.C

# Run ROOT macro
root.exe -b -q "Fun4All_G4_Pi0_Tbt_3.C(\"$1\")"

# Finish script

echo ' '
echo "END: $(date)"
echo ' '
