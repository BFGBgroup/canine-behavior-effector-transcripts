#!/usr/bin/bash
set -euo pipefail


PARALLEL=2
THREADS=30
EXPECTED_TASKS=19
EXPECTED_GROUPS=$'101\t100\n151\t150'

TSV="/path/to/project/legacy/five-dogs/data/statistics/five_dogs_fastq_pair_stats.tsv"
READLEN_SOURCE_DIR="/path/to/project/analysis_20260622/5-dogs/10.finemaping/preparation/1.1star-index"
GENOME_FASTA="/path/to/project/legacy/reference/GCF_000002285.3_CanFam3.1_genomic_chr.fna"
GTF_FILE="/path/to/project/analysis_20260817/reference/Manual-annotation/1.1_build_five_only_new_gene_gtf/20260818_204722/GCF_000002285.3_CanFam3.1_genomic_chr.plus_five_only_new_gene.gtf"

SCRIPT_DIR="/path/to/repository/scripts/08_expression_ase_integration/preparation/1.1star-index"
SINGLE_SCRIPT="${SCRIPT_DIR}/1.1.build_star_index_single-bulk-5.sh"
LOG_DIR="${SCRIPT_DIR}/log"
OUTPUT_ROOT="/path/to/project/analysis_20260817/5-dogs/10.finemaping/preparation/1.1star-index"

BASH_BIN="/usr/bin/bash"
AWK_BIN="/usr/bin/awk"
SORT_BIN="/usr/bin/sort"
XARGS_BIN="/usr/bin/xargs"
FLOCK_BIN="/usr/bin/flock"
CUT_BIN="/usr/bin/cut"
HEAD_BIN="/usr/bin/head"
TR_BIN="/usr/bin/tr"
MKTEMP_BIN="/usr/bin/mktemp"
MV_BIN="/usr/bin/mv"
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
MANIFEST="${RUN_DIR}/five_dogs_star_index_manifest.tsv"
MASTER_LOG="${LOG_DIR}/1.1.batch_build_star_index_five_only.${RUN_ID}.log"
SUCCESS_LIST="${LOG_DIR}/1.1.batch_build_star_index_five_only.${RUN_ID}.success.tsv"
FAILED_LIST="${LOG_DIR}/1.1.batch_build_star_index_five_only.${RUN_ID}.failed.tsv"
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
printf 'run_id\tmax_read_length\tsjdbOverhang\tstatus\texit_code\tstarted_at\tfinished_at\toutput_path\ttask_log\treason\n' > "${SUCCESS_LIST}"
printf 'run_id\tmax_read_length\tsjdbOverhang\tstatus\texit_code\tstarted_at\tfinished_at\toutput_path\ttask_log\treason\n' > "${FAILED_LIST}"

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

log "start: build five dogs five-only shared STAR index"
log "RUN_ID=${RUN_ID}"
log "TSV=${TSV}"
log "READLEN_SOURCE_DIR=${READLEN_SOURCE_DIR}"
log "GENOME_FASTA=${GENOME_FASTA}"
log "GTF_FILE=${GTF_FILE}"
log "SINGLE_SCRIPT=${SINGLE_SCRIPT}"
log "RUN_DIR=${RUN_DIR}"
log "PARALLEL=${PARALLEL}"
log "THREADS_PER_STAR_TASK=${THREADS}"
log "MASTER_LOG=${MASTER_LOG}"
log "SUCCESS_LIST=${SUCCESS_LIST}"
log "FAILED_LIST=${FAILED_LIST}"

for required_file in "${TSV}" "${GENOME_FASTA}" "${GTF_FILE}" "${SINGLE_SCRIPT}"; do
  [[ -s "${required_file}" ]] || { log "error: required input is missing or empty ${required_file}"; exit 1; }
done
[[ -d "${READLEN_SOURCE_DIR}" ]] || { log "error: cannot find read record directory ${READLEN_SOURCE_DIR}"; exit 1; }

EXPECTED_TASK_HEADER=$'animal_id\tterm\traw_data_name'
ACTUAL_TASK_HEADER="$(${HEAD_BIN} -n 1 "${TSV}" | ${TR_BIN} -d '\r')"
[[ "${ACTUAL_TASK_HEADER}" == "${EXPECTED_TASK_HEADER}" ]] || {
  log "error: task table unexpected header: ${ACTUAL_TASK_HEADER}"
  exit 1
}

RAW_TASK_ROWS="$(( $(${WC_BIN} -l < "${TSV}") - 1 ))"
[[ "${RAW_TASK_ROWS}" -eq "${EXPECTED_TASKS}" ]] || {
  log "error: task-table data row count is ${RAW_TASK_ROWS}, expected ${EXPECTED_TASKS}"
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

TOTAL_TASKS="${#TASKS[@]}"
[[ "${TOTAL_TASKS}" -eq "${EXPECTED_TASKS}" ]] || {
  log "error: unique animal_id+term task count is ${TOTAL_TASKS}, expected ${EXPECTED_TASKS}"
  exit 1
}

STATE_DIR="$(${MKTEMP_BIN} -d -t finemaping_star_index_batch_XXXXXX)"
MANIFEST_TMP="${STATE_DIR}/five_dogs_star_index_manifest.tsv"
EXPECTED_READLEN_HEADER=$'animal_id\tterm\traw_dir\tfq1_path\tfq2_path\treads1\treads2\tmin_len1\tmax_len1\tuniq_len_count1\tmin_len2\tmax_len2\tuniq_len_count2\tmax_read_length\tsjdbOverhang\tgenome_fasta\tgtf_file'

printf 'run_id\tanimal_id\tterm\tmax_read_length\tsjdbOverhang\tfq1_path\tfq2_path\tstar_index_dir\tgenome_fasta\tgtf_file\tannotation_scope\treadlen_source\thistorical_readlen_gtf\n' > "${MANIFEST_TMP}"

log "stage 1/3: check 19 tasks, historical read record and clean FASTQ"
for task in "${TASKS[@]}"; do
  IFS=$'\t' read -r animal term <<< "${task}"
  readlen_file="${READLEN_SOURCE_DIR}/${animal}/${term}/${animal}_${term}.readlen.tsv"
  [[ -s "${readlen_file}" ]] || { log "error: read record missing or is empty ${readlen_file}"; exit 1; }

  readlen_rows="$(${WC_BIN} -l < "${readlen_file}")"
  [[ "${readlen_rows}" -eq 2 ]] || {
    log "error: read record expected column header plus 1 data records, actual ${readlen_rows} row ${readlen_file}"
    exit 1
  }

  actual_readlen_header="$(${HEAD_BIN} -n 1 "${readlen_file}" | ${TR_BIN} -d '\r')"
  [[ "${actual_readlen_header}" == "${EXPECTED_READLEN_HEADER}" ]] || {
    log "error: read record unexpected header ${readlen_file}"
    exit 1
  }

  IFS=$'\t' read -r recorded_animal recorded_term fq1_path fq2_path max_read_length sjdb_overhang recorded_genome_fasta historical_gtf < <(
    "${AWK_BIN}" -F'\t' 'NR == 2 {gsub(/\r/, "", $0); print $1 "\t" $2 "\t" $4 "\t" $5 "\t" $14 "\t" $15 "\t" $16 "\t" $17}' "${readlen_file}"
  )

  [[ "${recorded_animal}" == "${animal}" && "${recorded_term}" == "${term}" ]] || {
    log "error: recorded animal_id/term does not match the task table: ${readlen_file}"
    exit 1
  }
  [[ "${recorded_genome_fasta}" == "${GENOME_FASTA}" ]] || {
    log "error: reference FASTA in the historical read record differs from the current FASTA: ${readlen_file}"
    exit 1
  }
  [[ -s "${fq1_path}" ]] || { log "error: R1 missing or is empty ${fq1_path}"; exit 1; }
  [[ -s "${fq2_path}" ]] || { log "error: R2 missing or is empty ${fq2_path}"; exit 1; }
  [[ "${max_read_length}" =~ ^[0-9]+$ && "${sjdb_overhang}" =~ ^[0-9]+$ ]] || {
    log "error: read length or sjdbOverhang is not nonnegative integer ${readlen_file}"
    exit 1
  }
  [[ "$((max_read_length - 1))" -eq "${sjdb_overhang}" ]] || {
    log "error: sjdbOverhang does not equal max_read_length-1 ${readlen_file}"
    exit 1
  }

  star_index_dir="${RUN_DIR}/shared/sjdbOverhang_${sjdb_overhang}/STARindex"
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "${RUN_ID}" "${animal}" "${term}" "${max_read_length}" "${sjdb_overhang}" \
    "${fq1_path}" "${fq2_path}" "${star_index_dir}" "${GENOME_FASTA}" "${GTF_FILE}" \
    "five-dogs-only" "${readlen_file}" "${historical_gtf}" >> "${MANIFEST_TMP}"
done

if ! GROUPS_TEXT="$(${AWK_BIN} -F'\t' 'NR > 1 {print $4 "\t" $5}' "${MANIFEST_TMP}" | ${SORT_BIN} -n -k2,2 -u)"; then
  log "error: Cannot extract shared index groups from the temporary manifest ${MANIFEST_TMP}"
  exit 1
fi
[[ -n "${GROUPS_TEXT}" ]] || { log "error: no shared index groups were extracted"; exit 1; }
[[ "${GROUPS_TEXT}" == "${EXPECTED_GROUPS}" ]] || {
  log "error: actual read group is not expected of 101/100 and 151/150"
  printf '%s\n' "${GROUPS_TEXT}"
  exit 1
}

"${MV_BIN}" "${MANIFEST_TMP}" "${MANIFEST}"
TOTAL_GROUPS="$(printf '%s\n' "${GROUPS_TEXT}" | "${AWK_BIN}" 'END {print NR}')"
log "input check completed: ${TOTAL_TASKS} libraries, ${TOTAL_GROUPS} items shared index group; MANIFEST=${MANIFEST}"

COMPLETED_FILE="${STATE_DIR}/completed.count"
SUCCESS_COUNT_FILE="${STATE_DIR}/success.count"
FAIL_COUNT_FILE="${STATE_DIR}/fail.count"
LOCK_FILE="${STATE_DIR}/lock"
printf '0\n' > "${COMPLETED_FILE}"
printf '0\n' > "${SUCCESS_COUNT_FILE}"
printf '0\n' > "${FAIL_COUNT_FILE}"

run_one() {
  local max_read_length="$1"
  local sjdb_overhang="$2"
  local task_log="${LOG_DIR}/1.1.build_star_index_five_only.sjdbOverhang_${sjdb_overhang}.${RUN_ID}.log"
  local output_path="${RUN_DIR}/shared/sjdbOverhang_${sjdb_overhang}/STARindex"
  local started_at finished_at rc

  started_at="$(${DATE_BIN} '+%F %T')"
  log "start index group: max_read_length=${max_read_length} sjdbOverhang=${sjdb_overhang} TASK_LOG=${task_log}"
  set +e
  "${BASH_BIN}" "${SINGLE_SCRIPT}" \
    "${max_read_length}" "${sjdb_overhang}" "${RUN_DIR}" "${GENOME_FASTA}" "${GTF_FILE}" \
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
        "${RUN_ID}" "${max_read_length}" "${sjdb_overhang}" "${started_at}" "${finished_at}" \
        "${output_path}" "${task_log}" >> "${SUCCESS_LIST}"
    else
      printf '%s\n' "$(( $(<"${FAIL_COUNT_FILE}") + 1 ))" > "${FAIL_COUNT_FILE}"
      printf '%s\t%s\t%s\tfailed\t%s\t%s\t%s\t%s\t%s\tworker_exit_nonzero\n' \
        "${RUN_ID}" "${max_read_length}" "${sjdb_overhang}" "${rc}" "${started_at}" "${finished_at}" \
        "${output_path}" "${task_log}" >> "${FAILED_LIST}"
    fi
  ) 200> "${LOCK_FILE}"

  log "index group finished: completed $(<"${COMPLETED_FILE}")/${TOTAL_GROUPS}, success $(<"${SUCCESS_COUNT_FILE}"), failed $(<"${FAIL_COUNT_FILE}"), last sjdbOverhang=${sjdb_overhang} rc=${rc}"
  return "${rc}"
}

export -f timestamp log run_one
export RUN_ID RUN_DIR LOG_DIR THREADS SINGLE_SCRIPT GENOME_FASTA GTF_FILE
export BASH_BIN FLOCK_BIN CUT_BIN DATE_BIN COMPLETED_FILE SUCCESS_COUNT_FILE FAIL_COUNT_FILE LOCK_FILE TOTAL_GROUPS SUCCESS_LIST FAILED_LIST

log "stage 2/3: parallel build two sets shared index"
set +e
printf '%s\n' "${GROUPS_TEXT}" | "${XARGS_BIN}" -I {} -P "${PARALLEL}" "${BASH_BIN}" -lc '
  group="$1"
  max_read_length=$(printf "%s\n" "${group}" | "${CUT_BIN}" -f1)
  sjdb_overhang=$(printf "%s\n" "${group}" | "${CUT_BIN}" -f2)
  run_one "${max_read_length}" "${sjdb_overhang}"
' _ {}
XARGS_RC=$?
set -e

COMPLETED="$(<"${COMPLETED_FILE}")"
SUCCESS="$(<"${SUCCESS_COUNT_FILE}")"
FAIL="$(<"${FAIL_COUNT_FILE}")"
log "stage 3/3: summary execution status"
log "completed ${COMPLETED}/${TOTAL_GROUPS}, success ${SUCCESS}, failed ${FAIL}, xargs_rc=${XARGS_RC}"

if [[ "${COMPLETED}" -ne "${TOTAL_GROUPS}" || "${SUCCESS}" -ne "${TOTAL_GROUPS}" || "${FAIL}" -ne 0 || "${XARGS_RC}" -ne 0 ]]; then
  log "error: at least one items shared index task not successful completed"
  log "execution_status=failed"
  exit 1
fi

printf 'run_id\t%s\nexecution_status\tsuccess\ncompleted_at\t%s\n' \
  "${RUN_ID}" "$(${DATE_BIN} '+%F %T')" > "${COMPLETION_MARKER}"
log "COMPLETION_MARKER=${COMPLETION_MARKER}"
log "task completed / done"
log "execution_status=success"
log "reminder: execution completed does not imply index completeness, format, count or biological quality already passed; please in run after separately check result."
