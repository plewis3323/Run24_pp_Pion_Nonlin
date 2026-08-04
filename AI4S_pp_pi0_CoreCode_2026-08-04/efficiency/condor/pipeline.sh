#!/usr/bin/bash
# =====================================================================
#  pipeline.sh  --  universal control script for the pi0/eta efficiency
#                   analysis over the MB pythia8_Detroit G4Hits DSTs.
#
#  ONE script you edit. Flip the STAGE TOGGLES on/off and tune CONFIG,
#  then just run:  ./pipeline.sh
#
#  Stages (run in order, each independently toggleable):
#    1. SUBMIT  -> condor_submit one job per background DST index (EMBED
#                  model): each job embeds the gun into one background
#                  G4Hits DST and runs ALL of its events.
#    2. WAIT    -> block until the submitted condor cluster drains.
#    3. HADD    -> merge the per-job ROOT files into one final ROOT,
#                  either locally or via a condorized bucketed hadd.
#
#  Typical uses:
#    - First time / production:  DO_SUBMIT=1 DO_WAIT=1 DO_HADD=1
#    - Re-merge only:            DO_SUBMIT=0 DO_WAIT=0 DO_HADD=1
#    - Small test:               set TEST_NFILES + NJOBS small below
# =====================================================================
set -euo pipefail
cd "$(cd "$(dirname "$0")" && pwd)"   # always run from macro/condor

# ---------------------------------------------------------------------
#  STAGE TOGGLES  (1 = run this stage, 0 = skip)
# ---------------------------------------------------------------------
# (each may be overridden on the command line, e.g.  DO_HADD=0 ./pipeline.sh)
DO_SUBMIT="${DO_SUBMIT:-1}"   # submit the analysis condor jobs
DO_WAIT="${DO_WAIT:-1}"       # wait for those jobs to finish before hadd
DO_HADD="${DO_HADD:-1}"       # merge the outputs

# ---------------------------------------------------------------------
#  CONFIG  (edit these; env overrides win, e.g.  PDGID=221 ./pipeline.sh)
# ---------------------------------------------------------------------
PDGID="${PDGID:-111}"                       # 111 = pi0, 221 = eta
NJOBS="${NJOBS:-10000}"                      # analysis jobs = background DSTs to process
TEST_NFILES="${TEST_NFILES:-}"              # "" = use NJOBS; or e.g. 4 for a quick test

# hadd settings
HADD_MODE="${HADD_MODE:-local}"             # "local" or "condor"
HADD_BUCKETS="${HADD_BUCKETS:-100}"         # condor hadd: number of parallel bucket jobs
HADD_THREADS="${HADD_THREADS:-8}"           # ROOT hadd -j threads (per hadd call)
FINAL_OUT_DIR="${FINAL_OUT_DIR:-/direct/sphenix+tg+tg01/bulk/plewis3323/P+P_NonLin_Run24_Space/Analysis_space/pi0eta_eff_pythia/output}"
# Scratch dir for hadd's parallel (-j) partial files. MUST live on a roomy
# filesystem: with -j N, hadd writes N large partial ROOT files here. The
# default (/tmp) is the ~31GB node root LV and overflows on big merges,
# which aborts hadd leaving a header-only output. Point it at Lustre.
HADD_TMP="${HADD_TMP:-${FINAL_OUT_DIR}/.haddtmp}"

# ---------------------------------------------------------------------
#  Derived
# ---------------------------------------------------------------------
SPECIES=$([ "${PDGID}" -eq 221 ] && echo eta || echo pi0)
FINAL_OUT="${FINAL_OUT_DIR}/eff_${SPECIES}_run24pp_MERGED.root"
SUBMIT_CLUSTER=""

# =====================================================================
#  ENV  (sphenix_setup references unset vars, so relax -u just for it)
# =====================================================================
set +u
source /cvmfs/sphenix.sdcc.bnl.gov/alma9.2-gcc-14.2.0/opt/sphenix/core/bin/sphenix_setup.sh -n new
set -u

# ---------------------------------------------------------------------
#  wait_for_cluster <clusterid> [poll_interval_secs]
#  Robustly block until a condor cluster fully drains.
#
#  Guards against the transient-empty-condor_q trap: under a big (10k-job)
#  load the schedd can momentarily time out / return no rows, which naively
#  looks like "drained" and would let the downstream hadd run on unfinished
#  jobs (that once produced an empty 883-byte merge). We therefore
#    (a) only trust a poll whose condor_q EXIT STATUS is 0, and
#    (b) require several CONSECUTIVE confirmed-empty polls before declaring
#        the cluster done.
#  We also condor_release each cycle so any job auto-held by PeriodicHold
#  ((NumJobStarts>=1 && JobStatus==1), i.e. evicted->idle) resumes instead
#  of sitting held forever and stalling the drain.
# ---------------------------------------------------------------------
wait_for_cluster() {
  local cl="$1"
  local interval="${2:-30}"
  local need="${WAIT_CONFIRMATIONS:-5}"   # consecutive empty polls required
  local zeros=0
  local out rc n
  while true; do
    rc=0
    out="$(condor_q "${cl}" -af ClusterId 2>/dev/null)" || rc=$?
    if [ "${rc}" -ne 0 ]; then
      # schedd hiccup: do NOT trust as empty; reset the streak
      zeros=0
    else
      n="$(printf '%s' "${out}" | grep -c . || true)"
      # recover evicted->held jobs so they finish rather than block the drain
      condor_release "${cl}" >/dev/null 2>&1 || true
      if [ "${n}" -eq 0 ]; then
        zeros=$((zeros + 1))
        if [ "${zeros}" -ge "${need}" ]; then break; fi
      else
        zeros=0
      fi
    fi
    sleep "${interval}"
  done
}

echo "=== pipeline.sh  species=${SPECIES} (pdgid=${PDGID})  NJOBS=${NJOBS}  hadd=${HADD_MODE} ==="

# =====================================================================
#  STAGE 1: SUBMIT
#  Delegates the per-index job queueing + condor_submit to run.sh
#  (passing config via env). Captures the cluster id so STAGE 2 can wait.
# =====================================================================
if [ "${DO_SUBMIT}" -eq 1 ]; then
  echo "--- [SUBMIT] ---"
  SUB_LOG="$(mktemp)"
  NJOBS="${NJOBS}" PDGID="${PDGID}" TEST_NFILES="${TEST_NFILES}" \
    ./run.sh 2>&1 | tee "${SUB_LOG}"
  SUBMIT_CLUSTER="$(grep -oE 'submitted to cluster [0-9]+' "${SUB_LOG}" | grep -oE '[0-9]+' | tail -1 || true)"
  rm -f "${SUB_LOG}"
  echo "[SUBMIT] cluster = ${SUBMIT_CLUSTER:-<none>}"
fi

# =====================================================================
#  STAGE 2: WAIT  (poll condor until the submitted cluster drains)
# =====================================================================
if [ "${DO_WAIT}" -eq 1 ]; then
  echo "--- [WAIT] ---"
  if [ -n "${SUBMIT_CLUSTER}" ]; then
    echo "[WAIT] waiting on cluster ${SUBMIT_CLUSTER} ..."
    wait_for_cluster "${SUBMIT_CLUSTER}" 30
    echo "[WAIT] cluster ${SUBMIT_CLUSTER} finished."
  else
    echo "[WAIT] no cluster id captured (nothing submitted this run); skipping wait."
  fi
fi

# =====================================================================
#  STAGE 3: HADD
#  Merges condorout/OutDir*/eff_<species>_run24pp_*.root -> FINAL_OUT.
#  Two interchangeable blocks below; pick with HADD_MODE.
# =====================================================================
if [ "${DO_HADD}" -eq 1 ]; then
  echo "--- [HADD] (${HADD_MODE}) ---"
  mkdir -p "${FINAL_OUT_DIR}"

  HADD_JOPT=()
  if [ "${HADD_THREADS}" -gt 0 ]; then
    mkdir -p "${HADD_TMP}"
    # -j: parallel merge; -d: keep the (large) partial files off /tmp
    HADD_JOPT=(-j "${HADD_THREADS}" -d "${HADD_TMP}")
  fi

  shopt -s nullglob
  ALL_OUT=( condorout/OutDir*/eff_${SPECIES}_run24pp_*.root )
  shopt -u nullglob
  echo "[HADD] found ${#ALL_OUT[@]} per-job files for species=${SPECIES}"
  if [ "${#ALL_OUT[@]}" -eq 0 ]; then
    echo "[HADD] nothing to merge; exiting." >&2
    exit 1
  fi

  # ---------------- HADD BLOCK A: LOCAL ----------------
  if [ "${HADD_MODE}" = "local" ]; then
    rm -f "${FINAL_OUT}"
    hadd "${HADD_JOPT[@]}" -k "${FINAL_OUT}" "${ALL_OUT[@]}"
    echo "[HADD] local merge -> ${FINAL_OUT}"

  # ---------------- HADD BLOCK B: CONDOR (bucketed) ----------------
  elif [ "${HADD_MODE}" = "condor" ]; then
    BASE_DIR="$(pwd)"
    HOUT="${BASE_DIR}/haddCondorOut"
    mkdir -p "${HOUT}"
    rm -f "${HOUT}"/out_${SPECIES}_*.root "${HOUT}"/master_hadd.sub
    chmod +x "${BASE_DIR}/CondorHadd.sh"

    cat > "${HOUT}/master_hadd.sub" <<EOF
+JobFlavour      = "workday"
Universe         = vanilla
Notification     = Never
GetEnv           = True
request_memory   = 2000MB
Executable       = ${BASE_DIR}/CondorHadd.sh
Initialdir       = ${BASE_DIR}
Output           = ${HOUT}/hadd_\$(Process).out
Error            = ${HOUT}/hadd_\$(Process).err
Log              = /tmp/condor_pi0eta_hadd.\$(Cluster).log
Arguments        = \$(Bucket) ${HADD_BUCKETS} ${SPECIES}
PeriodicHold     = (NumJobStarts>=1 && JobStatus == 1)

queue Bucket from (
EOF
    for ((b=0; b<HADD_BUCKETS; b++)); do echo "${b}" >> "${HOUT}/master_hadd.sub"; done
    echo ")" >> "${HOUT}/master_hadd.sub"

    HSUB_LOG="$(mktemp)"
    condor_submit "${HOUT}/master_hadd.sub" | tee "${HSUB_LOG}"
    HCLUSTER="$(grep -oE 'submitted to cluster [0-9]+' "${HSUB_LOG}" | grep -oE '[0-9]+' | tail -1 || true)"
    rm -f "${HSUB_LOG}"

    echo "[HADD] waiting on hadd cluster ${HCLUSTER} ..."
    wait_for_cluster "${HCLUSTER}" 20

    # Final combine of the per-bucket partials
    shopt -s nullglob
    BUCKETS=( "${HOUT}"/out_${SPECIES}_*.root )
    shopt -u nullglob
    rm -f "${FINAL_OUT}"
    hadd "${HADD_JOPT[@]}" -k "${FINAL_OUT}" "${BUCKETS[@]}"
    echo "[HADD] condor merge (${#BUCKETS[@]} buckets) -> ${FINAL_OUT}"
  else
    echo "[HADD] unknown HADD_MODE='${HADD_MODE}' (use 'local' or 'condor')" >&2
    exit 1
  fi
fi

echo "=== pipeline.sh done ==="
