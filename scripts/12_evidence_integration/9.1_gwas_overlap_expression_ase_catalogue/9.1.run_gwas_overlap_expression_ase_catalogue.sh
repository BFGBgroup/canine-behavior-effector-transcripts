#!/usr/bin/bash
set -euo pipefail


SCRIPT_DIR="/path/to/repository/scripts/12_evidence_integration/9.1_gwas_overlap_expression_ase_catalogue"
SCRIPT="${SCRIPT_DIR}/9.1.overlap_gwas_with_no_gwas_expression_ase_catalogue.py"
OUTPUT_ROOT="/path/to/project/analysis_20260817/5-dogs/10.finemaping/formal/9.1_gwas_overlap_expression_ase_catalogue"
LOG_DIR="${SCRIPT_DIR}/log"
RUN_ID="${RUN_ID:-$(/usr/bin/date '+%Y%m%d_%H%M%S')}"
RUN_DIR="${OUTPUT_ROOT}/${RUN_ID}"
MASTER_LOG="${LOG_DIR}/9.1.gwas_overlap_expression_ase_catalogue.${RUN_ID}.log"
COMPLETION_MARKER="${RUN_DIR}/execution_completion.marker"
SUCCESS_FILE="${LOG_DIR}/9.1.gwas_overlap_expression_ase_catalogue.${RUN_ID}.success.tsv"
FAILED_FILE="${LOG_DIR}/9.1.gwas_overlap_expression_ase_catalogue.${RUN_ID}.failed.tsv"
CONCURRENT_TASKS=1
THREADS_PER_TASK=1

log() {
    printf '[%s ] %s\n' "$(/usr/bin/date '+%F %T')" "$*"
}

if [[ ! "${RUN_ID}" =~ ^[0-9]{8}_[0-9]{6}$ ]]; then
    log "[error] RUN_ID must use YYYYMMDD_HHMMSS format: ${RUN_ID}"
    exit 1
fi

if [[ ! -f "${SCRIPT}" ]]; then
    log "[error] Python main program does not exist: ${SCRIPT}"
    exit 1
fi

if [[ -e "${RUN_DIR}" || -e "${MASTER_LOG}" || -e "${SUCCESS_FILE}" || -e "${FAILED_FILE}" ]]; then
    log "[error] run directory or main log already exists, refuse to overwrite: ${RUN_DIR}"
    exit 1
fi

/usr/bin/mkdir -p "${LOG_DIR}" "${OUTPUT_ROOT}"

/usr/bin/mkdir "${RUN_DIR}"
exec >>"${MASTER_LOG}" 2>&1
export RUN_ID
printf 'run_id\ttask_id\tstatus\texit_code\toutput_path\ttask_log\treason\n' > "${SUCCESS_FILE}"
printf 'run_id\ttask_id\tstatus\texit_code\toutput_path\ttask_log\treason\n' > "${FAILED_FILE}"
CHILD_PID=""
COMPLETED=0
finish_execution() {
    local code=$?
    if [[ "${COMPLETED}" -ne 1 ]]; then
        [[ "${code}" -ne 0 ]] || code=1
        printf '%s\t9.1\tfailed\t%s\t%s\t%s\tcommand_error_or_interruption\n' \
            "${RUN_ID}" "${code}" "${RUN_DIR}" "${MASTER_LOG}" >> "${FAILED_FILE}"
        log "execution_status=failed; exit_code=${code}"
    fi
}
stop_execution() {
    trap '' INT TERM HUP
    if [[ -n "${CHILD_PID}" ]]; then
        kill -TERM "${CHILD_PID}" 2>/dev/null || true
        wait "${CHILD_PID}" 2>/dev/null || true
    fi
    exit 143
}
trap finish_execution EXIT
trap stop_execution INT TERM HUP

log "[9.1] start"
log "[9.1] script=${SCRIPT}"
log "[9.1] run_id=${RUN_ID}"
log "[9.1] output=${RUN_DIR}"
log "[9.1] master_log=${MASTER_LOG}"
log "[9.1] concurrent_tasks=${CONCURRENT_TASKS}; threads_per_task=${THREADS_PER_TASK}; resource_baseline=preserved_serial_script"

set +e
/usr/bin/python3.10 -u "${SCRIPT}" &
CHILD_PID=$!
wait "${CHILD_PID}"
PYTHON_EXIT_CODE=$?
CHILD_PID=""
set -e

if (( PYTHON_EXIT_CODE != 0 )); then
    log "[error] Python main program run failed, exit_code=${PYTHON_EXIT_CODE}"
    log "execution_status=failed"
    exit "${PYTHON_EXIT_CODE}"
fi

printf 'run_id=%s\nexecution_status=success\nfinished_at=%s\n' \
    "${RUN_ID}" \
    "$(/usr/bin/date '+%F %T')" \
    >"${COMPLETION_MARKER}"

printf '%s\t9.1\tsuccess\t0\t%s\t%s\tnormal_completion\n' \
    "${RUN_ID}" "${RUN_DIR}" "${MASTER_LOG}" >> "${SUCCESS_FILE}"
COMPLETED=1
log "[done] [9.1] finish; task completed / done"
log "execution_status=success"
log "[reminder] Execution completion does not confirm output completeness, correct formatting, expected counts, or biological quality; check the results separately."
