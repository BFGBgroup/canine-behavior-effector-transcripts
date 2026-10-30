#!/usr/bin/bash
set -euo pipefail


PARALLEL=2
THREADS=20
SAMTOOLS_EXTRA_THREADS=19
EXPECTED_TASKS=19

TSV="/path/to/project/legacy/five-dogs/data/statistics/five_dogs_fastq_pair_stats.tsv"
MANIFEST="/path/to/project/analysis_20260817/5-dogs/10.finemaping/preparation/1.1star-index/20260820_005233/five_dogs_star_index_manifest.tsv"

SCRIPT_DIR="/path/to/repository/scripts/08_expression_ase_integration/preparation/1.2.1star"
SINGLE_SCRIPT="${SCRIPT_DIR}/1.2.1star.single-bulk-5.sh"
LOG_DIR="${SCRIPT_DIR}/log"
OUTPUT_ROOT="/path/to/project/analysis_20260817/5-dogs/10.finemaping/preparation/1.2.1star"

BASH_BIN="/usr/bin/bash"
AWK_BIN="/usr/bin/awk"
XARGS_BIN="/usr/bin/xargs"
FLOCK_BIN="/usr/bin/flock"
CUT_BIN="/usr/bin/cut"
HEAD_BIN="/usr/bin/head"
TR_BIN="/usr/bin/tr"
MKTEMP_BIN="/usr/bin/mktemp"
MKDIR_BIN="/usr/bin/mkdir"
RM_BIN="/usr/bin/rm"
WC_BIN="/usr/bin/wc"
DATE_BIN="/usr/bin/date"

RUN_ID="${RUN_ID:-$(${DATE_BIN} '+%Y%m%d_%H%M%S')}"
if [[ ! "${RUN_ID}" =~ ^[0-9]{8}_[0-9]{6}$ ]]; then
  echo "[ERROR] RUN_ID must use YYYYMMDD_HHMMSS format: ${RUN_ID}" >&2
  exit 1
fi

RUN_DIR="${OUTPUT_ROOT}/${RUN_ID}"
MASTER_LOG="${LOG_DIR}/1.2.1star.batch-bulk-5.${RUN_ID}.log"
SUCCESS_LIST="${LOG_DIR}/1.2.1star.batch-bulk-5.${RUN_ID}.success.tsv"
FAILED_LIST="${LOG_DIR}/1.2.1star.batch-bulk-5.${RUN_ID}.failed.tsv"
COMPLETION_MARKER="${RUN_DIR}/execution_completion.marker"

${MKDIR_BIN} -p "${OUTPUT_ROOT}" "${LOG_DIR}"
for reserved_path in "${MASTER_LOG}" "${SUCCESS_LIST}" "${FAILED_LIST}"; do
  if [[ -e "${reserved_path}" ]]; then
    echo "[ERROR] this run run of log or status file already exists, refuse to overwrite: ${reserved_path}" >&2
    exit 1
  fi
done
if ! ${MKDIR_BIN} "${RUN_DIR}"; then
  echo "[ERROR] The run directory already exists or cannot be created; refusing to overwrite: ${RUN_DIR}" >&2
  exit 1
fi

: > "${MASTER_LOG}"
exec >> "${MASTER_LOG}" 2>&1
printf 'run_id\tanimal_id\tterm\tstatus\texit_code\tstarted_at\tfinished_at\toutput_path\ttask_log\treason\n' > "${SUCCESS_LIST}"
printf 'run_id\tanimal_id\tterm\tstatus\texit_code\tstarted_at\tfinished_at\toutput_path\ttask_log\treason\n' > "${FAILED_LIST}"

timestamp() {
  "${DATE_BIN}" '+%F %T'
}

log() {
  printf '[%s ] %s\n' "$(timestamp)" "$*"
}

STATE_DIR=""
cleanup() {
  if [[ -n "${STATE_DIR}" && -d "${STATE_DIR}" ]]; then
    "${RM_BIN}" -rf "${STATE_DIR}"
  fi
}

on_error() {
  local rc=$?
  trap - ERR
  log "error: batch processing at step ${BASH_LINENO[0]} row exited abnormally, rc=${rc}"
  log "execution_status=failed"
  exit "${rc}"
}

on_signal() {
  local signal="$1"
  trap - ERR
  log "received ${signal} signal, batch processing terminated"
  log "execution_status=failed"
  exit 1
}

trap cleanup EXIT
trap on_error ERR
trap 'on_signal TERM' TERM
trap 'on_signal INT' INT

log "start: batch run five dogs five-only paired-end STAR alignment"
log "RUN_ID=${RUN_ID}"
log "TSV=${TSV}"
log "MANIFEST=${MANIFEST}"
log "SINGLE_SCRIPT=${SINGLE_SCRIPT}"
log "RUN_DIR=${RUN_DIR}"
log "PARALLEL=${PARALLEL}"
log "STAR_THREADS_PER_TASK=${THREADS}"
log "SAMTOOLS_ADDITIONAL_THREADS=${SAMTOOLS_EXTRA_THREADS}"
log "MASTER_LOG=${MASTER_LOG}"
log "SUCCESS_LIST=${SUCCESS_LIST}"
log "FAILED_LIST=${FAILED_LIST}"

for required_file in "${TSV}" "${MANIFEST}" "${SINGLE_SCRIPT}"; do
  [[ -s "${required_file}" ]] || { log "error: required input is missing or empty ${required_file}"; exit 1; }
done

EXPECTED_TASK_HEADER=$'animal_id\tterm\traw_data_name'
ACTUAL_TASK_HEADER="$(${HEAD_BIN} -n 1 "${TSV}" | ${TR_BIN} -d '\r')"
[[ "${ACTUAL_TASK_HEADER}" == "${EXPECTED_TASK_HEADER}" ]] || {
  log "error: task table unexpected header: ${ACTUAL_TASK_HEADER}"
  exit 1
}

EXPECTED_MANIFEST_HEADER=$'run_id\tanimal_id\tterm\tmax_read_length\tsjdbOverhang\tfq1_path\tfq2_path\tstar_index_dir\tgenome_fasta\tgtf_file\tannotation_scope\treadlen_source\thistorical_readlen_gtf'
ACTUAL_MANIFEST_HEADER="$(${HEAD_BIN} -n 1 "${MANIFEST}" | ${TR_BIN} -d '\r')"
[[ "${ACTUAL_MANIFEST_HEADER}" == "${EXPECTED_MANIFEST_HEADER}" ]] || {
  log "error: shared index manifest unexpected header: ${ACTUAL_MANIFEST_HEADER}"
  exit 1
}

RAW_TASK_ROWS="$(( $(${WC_BIN} -l < "${TSV}") - 1 ))"
RAW_MANIFEST_ROWS="$(( $(${WC_BIN} -l < "${MANIFEST}") - 1 ))"
[[ "${RAW_TASK_ROWS}" -eq "${EXPECTED_TASKS}" ]] || {
  log "error: task-table data row count is ${RAW_TASK_ROWS}, expected ${EXPECTED_TASKS}"
  exit 1
}
[[ "${RAW_MANIFEST_ROWS}" -eq "${EXPECTED_TASKS}" ]] || {
  log "error: manifest data row count is ${RAW_MANIFEST_ROWS}, expected ${EXPECTED_TASKS}"
  exit 1
}

mapfile -t TASKS < <(
  "${AWK_BIN}" -F'\t' '
    NR > 1 && $1 != "" && $2 != "" {
      gsub(/\r/, "", $1)
      gsub(/\r/, "", $2)
      print $1 "\t" $2
    }
  ' "${TSV}" | "${AWK_BIN}" '!seen[$0]++'
)

TOTAL="${#TASKS[@]}"
[[ "${TOTAL}" -eq "${EXPECTED_TASKS}" ]] || {
  log "error: unique animal_id+term task count is ${TOTAL}, expected ${EXPECTED_TASKS}"
  exit 1
}


if ! "${AWK_BIN}" -F'\t' '
  NR == FNR {
    if (FNR > 1) {
      gsub(/\r/, "", $1)
      gsub(/\r/, "", $2)
      task[$1 FS $2] = 1
    }
    next
  }
  FNR > 1 {
    gsub(/\r/, "", $2)
    gsub(/\r/, "", $3)
    manifest[$2 FS $3]++
  }
  END {
    for (key in task) {
      if (manifest[key] != 1) exit 1
    }
    for (key in manifest) {
      if (!(key in task)) exit 1
    }
  }
' "${TSV}" "${MANIFEST}"; then
  log "error: task table and 13 column manifest of animal_id/term is not strict one-to-one corresponding"
  exit 1
fi

STATE_DIR="$(${MKTEMP_BIN} -d -t finemaping_star_align_batch_XXXXXX)"
COMPLETED_FILE="${STATE_DIR}/completed.count"
SUCCESS_COUNT_FILE="${STATE_DIR}/success.count"
FAIL_COUNT_FILE="${STATE_DIR}/fail.count"
LOCK_FILE="${STATE_DIR}/lock"
printf '0\n' > "${COMPLETED_FILE}"
printf '0\n' > "${SUCCESS_COUNT_FILE}"
printf '0\n' > "${FAIL_COUNT_FILE}"

run_one() {
  local animal="$1"
  local term="$2"
  local task_log="${LOG_DIR}/1.2.1star.${animal}.${term}.${RUN_ID}.log"
  local output_path="${RUN_DIR}/${animal}/${term}/align_star"
  local started_at finished_at rc

  started_at="$(${DATE_BIN} '+%F %T')"
  log "start library: animal_id=${animal} term=${term} TASK_LOG=${task_log}"
  set +e
  THREADS="${THREADS}" SAMTOOLS_EXTRA_THREADS="${SAMTOOLS_EXTRA_THREADS}" \
    "${BASH_BIN}" "${SINGLE_SCRIPT}" "${animal}" "${term}" "${RUN_DIR}" \
    > "${task_log}" 2>&1
  rc=$?
  set -e
  finished_at="$(${DATE_BIN} '+%F %T')"

  (
    "${FLOCK_BIN}" -x 200
    printf '%s\n' "$(( $(<"${COMPLETED_FILE}") + 1 ))" > "${COMPLETED_FILE}"
    if [[ "${rc}" -eq 0 ]]; then
      printf '%s\n' "$(( $(<"${SUCCESS_COUNT_FILE}") + 1 ))" > "${SUCCESS_COUNT_FILE}"
      printf '%s\t%s\t%s\tsuccess\t0\t%s\t%s\t%s\t%s\tNA\n' \
        "${RUN_ID}" "${animal}" "${term}" "${started_at}" "${finished_at}" \
        "${output_path}" "${task_log}" >> "${SUCCESS_LIST}"
    else
      printf '%s\n' "$(( $(<"${FAIL_COUNT_FILE}") + 1 ))" > "${FAIL_COUNT_FILE}"
      printf '%s\t%s\t%s\tfailed\t%s\t%s\t%s\t%s\t%s\tworker_exit_nonzero\n' \
        "${RUN_ID}" "${animal}" "${term}" "${rc}" "${started_at}" "${finished_at}" \
        "${output_path}" "${task_log}" >> "${FAILED_LIST}"
    fi
  ) 200> "${LOCK_FILE}"

  log "library finished: completed $(<"${COMPLETED_FILE}")/${TOTAL}, success $(<"${SUCCESS_COUNT_FILE}"), failed $(<"${FAIL_COUNT_FILE}"), last=${animal} ${term} rc=${rc}"
  return 0
}

export -f timestamp log run_one
export RUN_ID RUN_DIR LOG_DIR THREADS SAMTOOLS_EXTRA_THREADS SINGLE_SCRIPT
export BASH_BIN FLOCK_BIN CUT_BIN DATE_BIN COMPLETED_FILE SUCCESS_COUNT_FILE FAIL_COUNT_FILE LOCK_FILE TOTAL SUCCESS_LIST FAILED_LIST

log "stage 1/2: parallel run 19 libraries"
set +e
printf '%s\n' "${TASKS[@]}" | "${XARGS_BIN}" -I {} -P "${PARALLEL}" "${BASH_BIN}" -lc '
  set -euo pipefail
  task="$1"
  animal=$(printf "%s\n" "${task}" | "${CUT_BIN}" -f1)
  term=$(printf "%s\n" "${task}" | "${CUT_BIN}" -f2)
  run_one "${animal}" "${term}"
' _ {}
XARGS_RC=$?
set -e

COMPLETED="$(<"${COMPLETED_FILE}")"
SUCCESS="$(<"${SUCCESS_COUNT_FILE}")"
FAIL="$(<"${FAIL_COUNT_FILE}")"
log "stage 2/2: summary execution status"
log "completed ${COMPLETED}/${TOTAL}, success ${SUCCESS}, failed ${FAIL}, xargs_rc=${XARGS_RC}"

if [[ "${COMPLETED}" -ne "${TOTAL}" || "${SUCCESS}" -ne "${TOTAL}" || "${FAIL}" -ne 0 || "${XARGS_RC}" -ne 0 ]]; then
  log "error: at least one libraries task not successful completed"
  log "execution_status=failed"
  exit 1
fi

printf 'run_id\t%s\nexecution_status\tsuccess\ncompleted_at\t%s\n' \
  "${RUN_ID}" "$(${DATE_BIN} '+%F %T')" > "${COMPLETION_MARKER}"
log "COMPLETION_MARKER=${COMPLETION_MARKER}"
log "task completed / done"
log "execution_status=success"
log "reminder: execution completed does not imply BAM, BAI, GeneCounts, alignment rate or biological quality already passed; please in run after separately check result."
