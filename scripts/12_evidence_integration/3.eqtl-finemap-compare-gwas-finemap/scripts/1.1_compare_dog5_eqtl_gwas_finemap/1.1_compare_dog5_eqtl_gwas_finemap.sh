#!/usr/bin/bash
set -euo pipefail


PYTHON_BIN="/usr/bin/python3.10"
SCRIPT_DIR="/path/to/repository/scripts/12_evidence_integration/3.eqtl-finemap-compare-gwas-finemap/scripts/1.1_compare_dog5_eqtl_gwas_finemap"
RUNNER="${SCRIPT_DIR}/1.1.compare_dog5_eqtl_gwas_finemap.py"
OUTPUT_ROOT="/path/to/project/analysis_20260817/results/5-dogs/3.eqtl-finemap-compare-gwas-finemap/1.1_compare_dog5_eqtl_gwas_finemap"
EQTL_EXPORT_DIR="/path/to/project/analysis_20260817/5-dogs/10.finemaping/formal/8.6_no_gwas_export_expression_ase_catalogue_deliverables/20260827_224912"
EQTL_FINAL_QC="/path/to/project/analysis_20260817/5-dogs/10.finemaping/formal/8.7_no_gwas_final_catalogue_qc_report/20260828_122551/no_gwas_expression_ase_final_qc_summary.tsv"
EQTL_FINAL_QC_COMPLETION="/path/to/project/analysis_20260817/5-dogs/10.finemaping/formal/8.7_no_gwas_final_catalogue_qc_report/20260828_122551/execution_completion.marker"
GWAS_RUN_DIR="/path/to/project/analysis_20260709/5-6dogs/2.gwas-finemap/work/1.10_run_dog_behavior_susie_rss_finemapping/20260802_205352"

AUDIT_DIR="${SCRIPT_DIR}/log"
LOCK_FILE="${AUDIT_DIR}/.1.1_compare_dog5_eqtl_gwas_finemap.lock"
CONCURRENT_TASKS=1
THREADS_PER_TASK=1

log() {
    printf '[%s ] %s\n' "$(/usr/bin/date '+%F %T')" "$*"
}

/usr/bin/mkdir -p "${AUDIT_DIR}" "${OUTPUT_ROOT}"
exec 9>"${LOCK_FILE}"
if ! /usr/bin/flock -n 9; then
    log "[error] existing the same step launcher run, refusing concurrent execution"
    exit 1
fi

RUN_ID="${RUN_ID:-$(/usr/bin/date '+%Y%m%d_%H%M%S')}"
RUN_DIR="${OUTPUT_ROOT}/${RUN_ID}"
COMPLETION_MARKER="${RUN_DIR}/execution_completion.marker"

if [[ ! "${RUN_ID}" =~ ^[0-9]{8}_[0-9]{6}$ ]]; then
    log "[error] RUN_ID must use YYYYMMDD_HHMMSS format: ${RUN_ID}"
    exit 1
fi
if [[ -e "${RUN_DIR}" ]]; then
    log "[error] run directory already exists, refuse to overwrite: ${RUN_DIR}"
    exit 1
fi

MASTER_LOG="${AUDIT_DIR}/1.1_compare_dog5_eqtl_gwas_finemap.${RUN_ID}.log"
SUCCESS_FILE="${AUDIT_DIR}/1.1_compare_dog5_eqtl_gwas_finemap.${RUN_ID}.success.tsv"
FAILED_FILE="${AUDIT_DIR}/1.1_compare_dog5_eqtl_gwas_finemap.${RUN_ID}.failed.tsv"

(set -o noclobber; : >"${MASTER_LOG}")
exec >>"${MASTER_LOG}" 2>&1
(set -o noclobber; printf 'run_id\tstatus\texit_code\tstarted_at\tfinished_at\toutput_path\ttask_log\treason\n' >"${SUCCESS_FILE}")
(set -o noclobber; printf 'run_id\tstatus\texit_code\tstarted_at\tfinished_at\toutput_path\ttask_log\treason\n' >"${FAILED_FILE}")

started_at="$(/usr/bin/date '+%F %T')"
CHILD_PID=""
FAILURE_REASON=""
NORMAL_END=0
MARKER_OWNED=0

finish() {
    local code=$?
    local finished_at
    trap - EXIT

    trap '' INT TERM
    finished_at="$(/usr/bin/date '+%F %T')"
    if (( code == 0 && NORMAL_END == 1 )); then
        MARKER_OWNED=1
        if printf 'run_id=%s\nexecution_status=success\nfinished_at=%s\n' \
            "${RUN_ID}" "${finished_at}" >"${COMPLETION_MARKER}"; then
            if printf '%s\tsuccess\t0\t%s\t%s\t%s\t%s\t%s\n' \
                "${RUN_ID}" "${started_at}" "${finished_at}" "${RUN_DIR}" \
                "${MASTER_LOG}" "WORKER_EXIT_0" >>"${SUCCESS_FILE}"; then
                log "[done] five dogs expression+ASE with GWAS fine-map exact compare task completed / done"
                log "execution_status=success"
                log "[reminder] Execution completion does not confirm output completeness, correct formatting, expected counts, or biological quality; check the results separately."
                exit 0
            fi
        fi
        code=1
        FAILURE_REASON="execution_status_write_failed"
    fi
    (( code != 0 )) || code=1
    if (( MARKER_OWNED == 1 )); then
        /usr/bin/rm -f -- "${COMPLETION_MARKER}"
    fi
    printf '%s\tfailed\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "${RUN_ID}" "${code}" "${started_at}" "${finished_at}" "${RUN_DIR}" \
        "${MASTER_LOG}" "${FAILURE_REASON:-launcher_exit_before_normal_end}" >>"${FAILED_FILE}"
    log "execution_status=failed; exit_code=${code}; reason=${FAILURE_REASON:-launcher_exit_before_normal_end}"
    exit "${code}"
}

stop_child() {
    trap '' INT TERM
    FAILURE_REASON="received_signal=$1"
    if [[ -n "${CHILD_PID}" ]]; then
        kill -TERM "${CHILD_PID}" 2>/dev/null || true
        wait "${CHILD_PID}" 2>/dev/null || true
    fi
    exit "$2"
}
trap finish EXIT
trap 'stop_child INT 130' INT
trap 'stop_child TERM 143' TERM

log "[start] run_id=${RUN_ID}"
log "[input] eqtl_export_dir=${EQTL_EXPORT_DIR}"
log "[input] eqtl_final_qc=${EQTL_FINAL_QC}"
log "[input] eqtl_final_qc_completion=${EQTL_FINAL_QC_COMPLETION}"
log "[input] gwas_run_dir=${GWAS_RUN_DIR}"
log "[output] run_dir=${RUN_DIR}"
log "[log] master_log=${MASTER_LOG}; success_file=${SUCCESS_FILE}; failed_file=${FAILED_FILE}"
log "[resource] concurrent_tasks=${CONCURRENT_TASKS}; threads_per_task=${THREADS_PER_TASK}; baseline=preserved_serial"

required_inputs=(
    "${RUNNER}"
    "${EQTL_EXPORT_DIR}/no_gwas_expression_ase_priority_variants_ppa_ge_0.1.tsv.gz"
    "${EQTL_EXPORT_DIR}/no_gwas_expression_ase_credible_set_variants.tsv.gz"
    "${EQTL_EXPORT_DIR}/no_gwas_expression_ase_export_summary.tsv"
    "${EQTL_EXPORT_DIR}/execution_completion.marker"
    "${EQTL_FINAL_QC}"
    "${EQTL_FINAL_QC_COMPLETION}"
    "${GWAS_RUN_DIR}/locus_run_manifest.tsv"
    "${GWAS_RUN_DIR}/overall_summary.tsv"
    "${GWAS_RUN_DIR}/branch_summary.tsv"
    "${GWAS_RUN_DIR}/model_parameters.tsv"
    "${GWAS_RUN_DIR}/COMPLETED.ok"
)
for required_path in "${required_inputs[@]}"; do
    if [[ ! -f "${required_path}" ]]; then
        FAILURE_REASON="missing_required_input=${required_path}"
        log "[error] required input does not exist: ${required_path}"
        exit 1
    fi
done

export OMP_NUM_THREADS="${THREADS_PER_TASK}"
export OPENBLAS_NUM_THREADS="${THREADS_PER_TASK}"
export MKL_NUM_THREADS="${THREADS_PER_TASK}"
export NUMEXPR_NUM_THREADS="${THREADS_PER_TASK}"
export VECLIB_MAXIMUM_THREADS="${THREADS_PER_TASK}"

set +e
"${PYTHON_BIN}" -B -u "${RUNNER}" \
    --run-dir "${RUN_DIR}" \
    --eqtl-export-dir "${EQTL_EXPORT_DIR}" \
    --eqtl-final-qc "${EQTL_FINAL_QC}" \
    --gwas-run-dir "${GWAS_RUN_DIR}" &
CHILD_PID=$!
wait "${CHILD_PID}"
exit_code=$?
CHILD_PID=""
set -e

if (( exit_code != 0 )); then
    FAILURE_REASON="worker_exit_code=${exit_code}"
    log "[failed] Python worker run failed, exit_code=${exit_code}"
    exit "${exit_code}"
fi

NORMAL_END=1
