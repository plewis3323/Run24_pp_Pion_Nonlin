#!/usr/bin/env bash
# Robust two-stage local merge of per-job condor outputs.
# Avoids ROOT 6.32 parallel-hadd (-j) flakiness over thousands of inputs:
#   Stage 1: bucket the inputs, sequential `hadd` per bucket (buckets run in parallel via xargs -P).
#   Stage 2: sequential `hadd` of the bucket partials into the final file.
set +u
source /cvmfs/sphenix.sdcc.bnl.gov/alma9.2-gcc-14.2.0/opt/sphenix/core/bin/sphenix_setup.sh -n new >/dev/null 2>&1
set -uo pipefail

SPECIES="${SPECIES:-pi0}"
BASE_DIR="/sphenix/tg/tg01/bulk/plewis3323/P+P_NonLin_Run24_Space/Analysis_space/pi0eta_eff_pythia"
COND_DIR="${BASE_DIR}/macro/condor/condorout"
FINAL_OUT="${BASE_DIR}/output/eff_${SPECIES}_run24pp_MERGED.root"
PART_DIR="${BASE_DIR}/macro/condor/haddParts"
BUCKET_SIZE="${BUCKET_SIZE:-100}"
NPROC="${NPROC:-8}"

rm -rf "${PART_DIR}"; mkdir -p "${PART_DIR}"

mapfile -t ALL < <(ls "${COND_DIR}"/OutDir*/eff_${SPECIES}_run24pp_*.root 2>/dev/null | sort)
N=${#ALL[@]}
echo "[merge] ${N} per-job files, bucket size ${BUCKET_SIZE}, ${NPROC} parallel buckets"
[ "${N}" -eq 0 ] && { echo "[merge] nothing to merge"; exit 1; }

# write bucket file lists
nb=0; i=0
while [ "${i}" -lt "${N}" ]; do
  printf '%s\n' "${ALL[@]:i:BUCKET_SIZE}" > "${PART_DIR}/bucket_${nb}.list"
  i=$((i+BUCKET_SIZE)); nb=$((nb+1))
done
echo "[merge] ${nb} buckets"

# Stage 1: one sequential hadd per bucket, run NPROC at a time
seq 0 $((nb-1)) | xargs -P "${NPROC}" -I{} bash -c '
  b={}
  pd="'"${PART_DIR}"'"
  out="${pd}/part_${b}.root"
  if hadd -k -f "${out}" $(cat "${pd}/bucket_${b}.list") > "${pd}/bucket_${b}.log" 2>&1; then
    echo "[bucket ${b}] OK"
  else
    echo "[bucket ${b}] FAILED (see ${pd}/bucket_${b}.log)"
  fi
'

# verify all parts present
mapfile -t PARTS < <(ls "${PART_DIR}"/part_*.root 2>/dev/null | sort -V)
echo "[merge] stage-1 produced ${#PARTS[@]}/${nb} parts"
if [ "${#PARTS[@]}" -ne "${nb}" ]; then
  echo "[merge] ERROR: missing bucket parts, aborting final merge"; exit 2
fi

# Stage 2: final sequential merge of the bucket partials
echo "[merge] stage-2 final hadd -> ${FINAL_OUT}"
hadd -k -f "${FINAL_OUT}" "${PARTS[@]}"
rc=$?
echo "[merge] final hadd rc=${rc}"
exit "${rc}"
