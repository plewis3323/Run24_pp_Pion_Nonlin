#!/usr/bin/bash
# ---------------------------------------------------------------------
#  Condor submitter for the pi0/eta efficiency analysis (EMBED model).
#
#  Each Condor job embeds the pi0/eta gun into ONE pre-simulated
#  background G4Hits DST and runs over ALL of its events. Job index q
#  (OutDir<q>) selects background file index q, i.e.
#     G4Hits_pythia8_Detroit-0000000028-<q>.root
#  The macro (Fun4All_G4_pythia_eff.C) builds that path from the index,
#  so there is nothing to split -- we just queue one job per index.
#
#  Env knobs:
#    NJOBS=10000        number of Condor jobs = number of background
#                       DSTs to process, indices 0 .. NJOBS-1 (default 10000)
#    PDGID=111          111 = pi0 (default), 221 = eta
#    TEST_NFILES=N      quick test: override NJOBS with a small N
#
#  Usage:
#    ./run.sh                       # production: pi0, NJOBS background DSTs
#    PDGID=221 ./run.sh             # production: eta
#    TEST_NFILES=2 ./run.sh         # tiny test: 2 jobs (indices 0,1)
# ---------------------------------------------------------------------
set -e

export TargetDir="$PWD"/condorout

export PDGID="${PDGID:-111}"
j="${NJOBS:-10000}"
# TEST_NFILES (if set) is a convenience alias for a small NJOBS.
if [ -n "${TEST_NFILES:-}" ]; then
  j="${TEST_NFILES}"
  echo "TEST mode: NJOBS=${j}"
fi

# Cap NJOBS at the number of available background DSTs so we never queue a
# job for a file index that does not exist.
BKG_LIST=/direct/sphenix+tg+tg01/bulk/plewis3323/P+P_NonLin_Run24_Space/Analysis_space/pi0eta_eff_pythia/SIM_DST_SPACE/MB_pythia_Detroit/g4hits.list
if [ -f "${BKG_LIST}" ]; then
  maxfiles=$( wc -l < "${BKG_LIST}" )
  if [ "${j}" -gt "${maxfiles}" ]; then
    echo "NJOBS ${j} > available background DSTs ${maxfiles}; capping to ${maxfiles}"
    j="${maxfiles}"
  fi
fi
njob="${j}"
echo "njob (background DSTs to process): ${njob}"

if [ "${njob}" -eq 0 ]; then
  echo "ERROR: njob=0 (nothing to submit)" >&2
  exit 1
fi

# Clean previous job dirs
mkdir -p "${TargetDir}"
if [ -n "$(ls -A ${TargetDir}/OutDir* 2>/dev/null)" ]; then
  rm -rf ${TargetDir}/OutDir*
fi
ABS_TARGET="$(cd "${TargetDir}" && pwd)"

# One shared copy of the per-job executable.
# rm first: cp onto a still-mapped executable fails with ETXTBSY ("Text file
# busy") if a prior run's copy is momentarily held (stale Lustre exec lock
# after killed jobs). Unlinking drops the old inode so cp writes a fresh one.
rm -f "${ABS_TARGET}/CondorRunJob.sh"
cp "$PWD/CondorRun.sh" "${ABS_TARGET}/CondorRunJob.sh"
chmod +x "${ABS_TARGET}/CondorRunJob.sh"

# Per-job directories (OutDir<q>, q = background/output index)
for ((q=0; q<njob; q++)); do
  mkdir -p "${ABS_TARGET}/OutDir${q}"
  ln -sf "${ABS_TARGET}/CondorRunJob.sh" "${ABS_TARGET}/OutDir${q}/CondorRunJob.sh"
done

# Master submission file
cat > "${ABS_TARGET}/master.sub" <<EOF
+JobFlavour                   = "workday"
Universe                      = vanilla
Notification                  = Never
GetEnv                        = True
request_memory                = 8000MB
Output                        = \$(WorkDir)/job.out
Error                         = \$(WorkDir)/job.err
Log                           = /tmp/condor_pi0eta.\$(Cluster).\$(Process).log
Executable                    = \$(WorkDir)/CondorRunJob.sh
Initialdir                    = \$(WorkDir)
PeriodicHold                  = (NumJobStarts>=1 && JobStatus == 1)

queue WorkDir from (
EOF

for ((q=0; q<njob; q++)); do
  echo "${ABS_TARGET}/OutDir${q}" >> "${ABS_TARGET}/master.sub"
done
echo ")" >> "${ABS_TARGET}/master.sub"

echo "Submitting ${njob} jobs (PDGID=${PDGID})..."
condor_submit "${ABS_TARGET}/master.sub"
