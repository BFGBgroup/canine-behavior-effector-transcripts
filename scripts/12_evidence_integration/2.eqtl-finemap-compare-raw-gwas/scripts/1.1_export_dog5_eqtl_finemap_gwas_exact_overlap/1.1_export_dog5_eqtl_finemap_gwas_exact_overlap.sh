#!/usr/bin/bash
set -euo pipefail


SCRIPT_DIR="/path/to/repository/scripts/12_evidence_integration/2.eqtl-finemap-compare-raw-gwas/scripts/1.1_export_dog5_eqtl_finemap_gwas_exact_overlap"
PYTHON_SCRIPT="${SCRIPT_DIR}/1.1.export_dog5_eqtl_finemap_gwas_exact_overlap.py"
INPUT_RUN_DIR="/path/to/project/analysis_20260817/5-dogs/10.finemaping/formal/9.1_gwas_overlap_expression_ase_catalogue/20260917_191751"
INPUT_EXACT_OVERLAPS="${INPUT_RUN_DIR}/gwas_catalogue_exact_variant_overlaps_p_le_1e-4.tsv"
INPUT_RUN_SUMMARY="${INPUT_RUN_DIR}/gwas_catalogue_overlap_run_summary.tsv"
INPUT_COMPLETION_MARKER="${INPUT_RUN_DIR}/execution_completion.marker"
OUTPUT_ROOT="/path/to/project/analysis_20260817/results/5-dogs/2.eqtl-finemap-compare-raw-gwas/1.1_export_dog5_eqtl_finemap_gwas_exact_overlap"
LOG_DIR="${SCRIPT_DIR}/log"
RUN_ID="${RUN_ID:-$(/usr/bin/date '+%Y%m%d_%H%M%S')}"
RUN_DIR="${OUTPUT_ROOT}/${RUN_ID}"
MASTER_LOG="${LOG_DIR}/1.1.export_dog5_eqtl_finemap_gwas_exact_overlap.${RUN_ID}.log"
COMPLETION_MARKER="${RUN_DIR}/execution_completion.marker"
CONCURRENT_TASKS=1
THREADS_PER_TASK=1

log() {
    printf '[%s ] %s\n' "$(/usr/bin/date '+%F %T')" "$*"
}

if [[ ! "${RUN_ID}" =~ ^[0-9]{8}_[0-9]{6}$ ]]; then
    log "[error] RUN_ID must use YYYYMMDD_HHMMSS format: ${RUN_ID}"
    exit 1
fi

/usr/bin/mkdir -p "${LOG_DIR}" "${OUTPUT_ROOT}"


if ! (set -o noclobber; : >"${MASTER_LOG}") 2>/dev/null; then
    log "[error] this run main log already exists, refuse to overwrite: ${MASTER_LOG}"
    exit 1
fi

exec >>"${MASTER_LOG}" 2>&1


SUCCESS_TSV="${LOG_DIR}/1.1.export_dog5_eqtl_finemap_gwas_exact_overlap.${RUN_ID}.success.tsv"
FAILED_TSV="${LOG_DIR}/1.1.export_dog5_eqtl_finemap_gwas_exact_overlap.${RUN_ID}.failed.tsv"
STATUS_HEADER='run_id\texecution_status\texit_code\toutput_path'
(set -o noclobber; printf '%b\n' "${STATUS_HEADER}" >"${SUCCESS_TSV}")
(set -o noclobber; printf '%b\n' "${STATUS_HEADER}" >"${FAILED_TSV}")
CHILD_PID=""
MARKER_OWNED=0
NORMAL_END=0
finish() {
    code=$?
    trap - EXIT
    if (( code == 0 && NORMAL_END == 1 )); then
        if printf '%s\tsuccess\t0\t%s\n' "${RUN_ID}" "${RUN_DIR}" >>"${SUCCESS_TSV}"; then
            log "task completed / done; execution_status=success"
            log "[reminder] execution completed does not imply output format, counts or biological quality have passed review; please check separately."
            exit 0
        fi
        code=1
    fi
    if (( code != 0 || NORMAL_END != 1 )); then
        (( code != 0 )) || code=1
        if (( MARKER_OWNED == 1 )); then
            /usr/bin/rm -f -- "${COMPLETION_MARKER}"
        fi
        printf '%s\tfailed\t%s\t%s\n' "${RUN_ID}" "${code}" "${RUN_DIR}" >>"${FAILED_TSV}"
        log "execution_status=failed; exit_code=${code}"
    fi
    exit "${code}"
}
stop_child() {
    trap '' INT TERM
    if [[ -n "${CHILD_PID}" ]]; then
        kill -TERM "${CHILD_PID}" 2>/dev/null || true
        wait "${CHILD_PID}" 2>/dev/null || true
    fi
    exit "$1"
}
trap finish EXIT
trap 'stop_child 130' INT
trap 'stop_child 143' TERM


log "[1.1] start"
log "[1.1] python_script=${PYTHON_SCRIPT}"
log "[1.1] input_run_dir=${INPUT_RUN_DIR}"
log "[1.1] run_id=${RUN_ID}"
log "[1.1] output=${RUN_DIR}"
log "[1.1] master_log=${MASTER_LOG}"
log "[1.1] concurrent_tasks=${CONCURRENT_TASKS}; threads_per_task=${THREADS_PER_TASK}; resource_baseline=preserved_serial_script"

for required_path in \
    "${PYTHON_SCRIPT}" \
    "${INPUT_EXACT_OVERLAPS}" \
    "${INPUT_RUN_SUMMARY}" \
    "${INPUT_COMPLETION_MARKER}"
do
    if [[ ! -f "${required_path}" ]]; then
        log "[error] required input does not exist: ${required_path}"
        exit 1
    fi
done

if [[ -e "${RUN_DIR}" ]]; then
    log "[error] run directory already exists, refuse to overwrite: ${RUN_DIR}"
    exit 1
fi

/usr/bin/mkdir "${RUN_DIR}"

export RUN_ID RUN_DIR INPUT_RUN_DIR
export OMP_NUM_THREADS="${THREADS_PER_TASK}"
export OPENBLAS_NUM_THREADS="${THREADS_PER_TASK}"
export MKL_NUM_THREADS="${THREADS_PER_TASK}"
export NUMEXPR_NUM_THREADS="${THREADS_PER_TASK}"
export VECLIB_MAXIMUM_THREADS="${THREADS_PER_TASK}"

set +e
/usr/bin/python3.10 -B -u "${PYTHON_SCRIPT}" &
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

MARKER_OWNED=1
printf 'run_id=%s\nexecution_status=success\nfinished_at=%s\n' \
    "${RUN_ID}" \
    "$(/usr/bin/date '+%F %T')" \
    >"${COMPLETION_MARKER}"

NORMAL_END=1
