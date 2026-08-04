#!/usr/bin/bash
# ---------------------------------------------------------------------
#  Per-bucket Condor hadd executable.
#  Merges the per-job eff_<species>_run24pp_*.root files whose OutDir
#  index falls in this bucket (idx % NBUCKETS == BUCKET) into one
#  partial file haddCondorOut/out_<species>_<bucket>.root.
#  A final combine of the bucket files is done by pipeline.sh.
#
#  Args: <bucket 0..NBUCKETS-1> <NBUCKETS> <species>
# ---------------------------------------------------------------------
set -euo pipefail

if [ "$#" -lt 2 ]; then
  echo "Usage: $0 <bucket> <nbuckets> [species]" >&2
  exit 2
fi
BUCKET="$1"
NBK="$2"
SPECIES="${3:-pi0}"

set +u
source /cvmfs/sphenix.sdcc.bnl.gov/alma9.2-gcc-14.2.0/opt/sphenix/core/bin/sphenix_setup.sh -n new
set -u

BASE_DIR="$(pwd)"          # Initialdir = macro/condor
OUT_DIR="${BASE_DIR}/haddCondorOut"
mkdir -p "${OUT_DIR}"

OUT_FILE="${OUT_DIR}/out_${SPECIES}_${BUCKET}.root"
rm -f "${OUT_FILE}"

shopt -s nullglob
CANDIDATES=( "${BASE_DIR}"/condorout/OutDir*/eff_${SPECIES}_run24pp_*.root )
shopt -u nullglob

INPUTS=()
for f in "${CANDIDATES[@]}"; do
  d="$(basename "$(dirname "${f}")")"
  idx="${d#OutDir}"
  [[ "${idx}" =~ ^[0-9]+$ ]] || continue
  if [ $(( idx % NBK )) -eq "${BUCKET}" ]; then
    INPUTS+=( "${f}" )
  fi
done

if [ "${#INPUTS[@]}" -eq 0 ]; then
  echo "[CondorHadd] bucket ${BUCKET}: no inputs, skipping."
  exit 0
fi

echo "[CondorHadd] bucket ${BUCKET}/${NBK} species=${SPECIES} inputs=${#INPUTS[@]} -> ${OUT_FILE}"
# -d: keep hadd's parallel partial files on the condor per-job scratch (sized
# and cleaned by condor) instead of the tiny worker /tmp, which can ENOSPC and
# abort the merge into a header-only file.
HADD_TMP="${_CONDOR_SCRATCH_DIR:-${OUT_DIR}}"
hadd -j 8 -d "${HADD_TMP}" -k "${OUT_FILE}" "${INPUTS[@]}"
echo "[CondorHadd] bucket ${BUCKET} done."
