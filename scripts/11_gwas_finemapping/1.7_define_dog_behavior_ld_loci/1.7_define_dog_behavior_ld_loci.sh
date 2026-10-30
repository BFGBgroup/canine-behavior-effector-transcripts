#!/usr/bin/env bash


set -Eeuo pipefail

SCRIPT_DIR="/path/to/repository/scripts/11_gwas_finemapping/1.7_define_dog_behavior_ld_loci"
PYTHON_BIN="/path/to/software/miniconda3/envs/dogs-harmony-latest/bin/python"
PLINK2_BIN="/path/to/software/miniconda3/envs/dogs-harmony-latest/bin/plink2"
WORKER="${SCRIPT_DIR}/1.7_define_dog_behavior_ld_loci.py"

STEP12_RUN_DIR="/path/to/project/analysis_20260709/5-6dogs/2.gwas-finemap/work/1.2_build_pnas_dog_behavior_locus_manifest/20260801_202130"
TRAIT_METADATA="${STEP12_RUN_DIR}/trait_metadata.tsv"
PNAS_TABLE1_LOCI="${STEP12_RUN_DIR}/pnas_table1_loci.tsv"
CANDIDATE_PRODUCER="/path/to/repository/scripts/11_gwas_finemapping/1.2_build_pnas_dog_behavior_locus_manifest/1.2_build_pnas_dog_behavior_locus_manifest.py"

STEP15_RUN_DIR="/path/to/project/analysis_20260709/5-6dogs/2.gwas-finemap/work/1.5_prepare_grls_pilot_ld_reference_plink2_king_primary/20260802_134457"
FINAL_LD_PREFIX="${STEP15_RUN_DIR}/ld_reference_qc_plink2_king_primary"
STEP15_COMPLETION="${STEP15_RUN_DIR}/COMPLETED.ok"
LD_PRODUCER="/path/to/repository/scripts/11_gwas_finemapping/1.5_prepare_grls_pilot_ld_reference_plink2_king_primary/1.5_prepare_grls_pilot_ld_reference_plink2_king_primary.py"

STEP16_RUN_DIR="/path/to/project/analysis_20260709/5-6dogs/2.gwas-finemap/work/1.6_audit_dog_behavior_candidate_ld_compatibility/20260802_154520"
COMPATIBILITY_TABLE="${STEP16_RUN_DIR}/candidate_ld_compatibility.tsv"
STEP16_COMPLETION="${STEP16_RUN_DIR}/COMPLETED.ok"
AUDIT_PRODUCER="/path/to/repository/scripts/11_gwas_finemapping/1.6_audit_dog_behavior_candidate_ld_compatibility/1.6_audit_dog_behavior_candidate_ld_compatibility.py"

OUTPUT_ROOT="/path/to/project/analysis_20260709/5-6dogs/2.gwas-finemap/work/1.7_define_dog_behavior_ld_loci"
THREADS="${THREADS:-16}"
LD_R2_THRESHOLD="0.7"
LD_WINDOW_KB="9999999"
LD_WINDOW_VARIANT_COUNT="99999999"
FALLBACK_WINDOW_BP="1000000"
RUN_ID="${RUN_ID:-$(date '+%Y%m%d_%H%M%S')}"
RUN_DIR="${OUTPUT_ROOT}/${RUN_ID}"

AUDIT_DIR="${SCRIPT_DIR}/log"
LOCK_FILE="${AUDIT_DIR}/1.7_define_dog_behavior_ld_loci.lock"
MASTER_LOG="${AUDIT_DIR}/1.7_define_dog_behavior_ld_loci.${RUN_ID}.log"
SUCCESS_FILE="${AUDIT_DIR}/1.7_define_dog_behavior_ld_loci.${RUN_ID}.success.tsv"
FAILED_FILE="${AUDIT_DIR}/1.7_define_dog_behavior_ld_loci.${RUN_ID}.failed.tsv"

TASK_ID="define_dog_behavior_ld_loci"
STARTED_AT="$(date '+%Y-%m-%d %H:%M:%S')"
STATUS_RECORDED=0

timestamp() {
    date '+[%Y-%m-%d %H:%M:%S ]'
}

log() {
    printf '%s %s\n' "$(timestamp)" "$*"
}

record_status() {
    local status="$1"
    local exit_code="$2"
    local reason="$3"
    local destination="${FAILED_FILE}"
    if [[ "${status}" == "success" ]]; then
        destination="${SUCCESS_FILE}"
    fi
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "${TASK_ID}" "${RUN_ID}" "${status}" "${exit_code}" \
        "${STARTED_AT}" "$(date '+%Y-%m-%d %H:%M:%S')" \
        "${RUN_DIR}" "${MASTER_LOG}" "${reason}" >> "${destination}"
    STATUS_RECORDED=1
}

on_error() {
    local exit_code=$?
    trap - ERR
    if [[ "${STATUS_RECORDED}" -eq 0 ]]; then
        record_status "failed" "${exit_code}" "launcher_or_worker_failed"
    fi
    log "task failed: exit_code=${exit_code}"
    log "failure list: ${FAILED_FILE}"
    exit "${exit_code}"
}

mkdir -p "${AUDIT_DIR}"
exec 9>"${LOCK_FILE}"
if ! flock -n 9; then
    printf '%s another task for this step is already running: %s\n' "$(timestamp)" "${LOCK_FILE}" >&2
    exit 1
fi

if [[ ! "${RUN_ID}" =~ ^[0-9]{8}_[0-9]{6}$ ]]; then
    printf '%s RUN_ID format invalid, expected YYYYMMDD_HHMMSS: %s\n' "$(timestamp)" "${RUN_ID}" >&2
    exit 1
fi

if [[ -e "${RUN_DIR}" || -e "${MASTER_LOG}" || -e "${SUCCESS_FILE}" || -e "${FAILED_FILE}" ]]; then
    printf '%s refuse to overwrite existing run: %s\n' "$(timestamp)" "${RUN_ID}" >&2
    exit 1
fi

: > "${MASTER_LOG}"
printf 'task_id\trun_id\tstatus\texit_code\tstarted_at\tfinished_at\toutput_path\ttask_log\treason\n' > "${SUCCESS_FILE}"
printf 'task_id\trun_id\tstatus\texit_code\tstarted_at\tfinished_at\toutput_path\ttask_log\treason\n' > "${FAILED_FILE}"
exec >> "${MASTER_LOG}" 2>&1
trap on_error ERR

log "task started"
log "candidate compatibility table: ${COMPATIBILITY_TABLE}"
log "final main LD prefix: ${FINAL_LD_PREFIX}"
log "PLINK2：${PLINK2_BIN}"
log "parameter: r2>=${LD_R2_THRESHOLD}; ld_window_kb=${LD_WINDOW_KB}; ld_window_variant_count=${LD_WINDOW_VARIANT_COUNT}; fallback_window_bp=${FALLBACK_WINDOW_BP}"
log "output directory: ${RUN_DIR}"
log "standard log: ${MASTER_LOG}"

required_nonempty=(
    "${PYTHON_BIN}"
    "${PLINK2_BIN}"
    "${WORKER}"
    "${COMPATIBILITY_TABLE}"
    "${TRAIT_METADATA}"
    "${PNAS_TABLE1_LOCI}"
    "${STEP16_COMPLETION}"
    "${CANDIDATE_PRODUCER}"
    "${AUDIT_PRODUCER}"
    "${FINAL_LD_PREFIX}.bed"
    "${FINAL_LD_PREFIX}.bim"
    "${FINAL_LD_PREFIX}.fam"
    "${STEP15_COMPLETION}"
    "${LD_PRODUCER}"
)
for required_path in "${required_nonempty[@]}"; do
    if [[ ! -s "${required_path}" ]]; then
        log "required input does not exist or is empty: ${required_path}"
        false
    fi
done

set +e
"${PYTHON_BIN}" "${WORKER}" \
    --compatibility-table "${COMPATIBILITY_TABLE}" \
    --trait-metadata "${TRAIT_METADATA}" \
    --pnas-table1-loci "${PNAS_TABLE1_LOCI}" \
    --step16-completion "${STEP16_COMPLETION}" \
    --candidate-producer "${CANDIDATE_PRODUCER}" \
    --audit-producer "${AUDIT_PRODUCER}" \
    --ld-prefix "${FINAL_LD_PREFIX}" \
    --step15-completion "${STEP15_COMPLETION}" \
    --ld-producer "${LD_PRODUCER}" \
    --plink2-bin "${PLINK2_BIN}" \
    --output-dir "${RUN_DIR}" \
    --threads "${THREADS}" \
    --ld-r2-threshold "${LD_R2_THRESHOLD}" \
    --ld-window-kb "${LD_WINDOW_KB}" \
    --ld-window-variant-count "${LD_WINDOW_VARIANT_COUNT}" \
    --fallback-window-bp "${FALLBACK_WINDOW_BP}"
WORKER_EXIT=$?
set -e

SUMMARY_FILE="${RUN_DIR}/locus_definition_summary.tsv"
LD_SUMMARY_FILE="${RUN_DIR}/ld_locus_definition_summary.tsv"
REGION_FILE="${RUN_DIR}/all_locus_regions.tsv"
ASSIGNMENT_FILE="${RUN_DIR}/candidate_locus_assignment.tsv"
NON_FINEMAPPABLE_FILE="${RUN_DIR}/non_finemappable_candidate_regions.tsv"
MERGE_AUDIT_FILE="${RUN_DIR}/locus_merge_audit.tsv"
REGION_SOURCE_AUDIT_FILE="${RUN_DIR}/locus_region_source_audit.tsv"
PNAS_RECONCILIATION_FILE="${RUN_DIR}/pnas_table1_locus_reconciliation.tsv"
COMPLETION_MARKER="${RUN_DIR}/COMPLETED.ok"
SUMMARY_OK=0
if [[ -s "${SUMMARY_FILE}" ]]; then
    if awk -F '\t' '
        $1 == "overall_status" && ($2 == "PASS" || $2 == "PASS_WITH_WARNINGS") {found=1}
        END {exit(found ? 0 : 1)}
    ' "${SUMMARY_FILE}"; then
        SUMMARY_OK=1
    fi
fi

if [[ "${WORKER_EXIT}" -eq 0 && "${SUMMARY_OK}" -eq 1 \
    && -s "${LD_SUMMARY_FILE}" && -s "${REGION_FILE}" \
    && -s "${ASSIGNMENT_FILE}" && -s "${NON_FINEMAPPABLE_FILE}" \
    && -s "${MERGE_AUDIT_FILE}" && -s "${REGION_SOURCE_AUDIT_FILE}" \
    && -s "${PNAS_RECONCILIATION_FILE}" \
    && -s "${COMPLETION_MARKER}" ]]; then
    record_status "success" "0" "dog_behavior_ld_locus_definition_completed"
    log "total tasks: 1"
    log "success: 1; failed: 0"
    log "success list: ${SUCCESS_FILE}"
    log "failure list: ${FAILED_FILE}"
    log "task completed / done"
    trap - ERR
    exit 0
fi

record_status "failed" "${WORKER_EXIT}" "worker_exit_or_output_validation_failed"
log "total tasks: 1"
log "success: 0; failed: 1"
log "failure list: ${FAILED_FILE}"
log "the task did not reach a valid completion state"
trap - ERR
exit 1


