#!/bin/bash
# ---------------------------------------------------------------------
#  Per-job Condor executable (runs inside Initialdir = OutDir<q>).
#  EMBED model: this job embeds the pi0/eta gun into ONE pre-simulated
#  background G4Hits DST, selected by the job index (JOBID), and runs
#  over ALL events in that DST (nEvents=0). The macro builds the
#  background path from the index:
#     G4Hits_pythia8_Detroit-0000000028-<JOBID>.root
#  Output ROOT (eff_<species>_run24pp_<JOBID>.root) lands in this OutDir.
# ---------------------------------------------------------------------

source /cvmfs/sphenix.sdcc.bnl.gov/alma9.2-gcc-14.2.0/opt/sphenix/core/bin/sphenix_setup.sh -n new

# Numeric job index from the OutDir name (OutDir<q>). This is BOTH the
# background-DST index the macro embeds into and the output file index.
JOBID=$(basename "$PWD" | tr -dc '0-9')
: "${JOBID:=0}"

# Species: 111 = pi0 (default), 221 = eta. Propagated from run.sh via GetEnv.
PDGID="${PDGID:-111}"

MACRO=/direct/sphenix+tg+tg01/bulk/plewis3323/P+P_NonLin_Run24_Space/Analysis_space/pi0eta_eff_pythia/macro/Fun4All_G4_pythia_eff.C

ls

# nEvents=0 -> process ALL events in this background DST (REPEAT=false
# terminates the loop). USE_EMBED=true; output into this OutDir (".").
# args: fileNum, pdgid, useWeights, nEvents, USE_EMBED, outDir
root -b -q "${MACRO}(${JOBID}, ${PDGID}, false, 0, true, \".\")"

echo "JOB COMPLETE!"
