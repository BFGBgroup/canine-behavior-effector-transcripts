#!/usr/bin/env bash


set -Eeuo pipefail

SCRIPT_DIR="/path/to/repository/scripts/11_gwas_finemapping/1.6_audit_dog_behavior_candidate_ld_compatibility"
PYTHON_BIN="/path/to/software/miniconda3/envs/dogs-harmony-latest/bin/python"
WORKER="${SCRIPT_DIR}/1.6_audit_dog_behavior_candidate_ld_compatibility.py"

STEP12_RUN_DIR="/path/to/project/analysis_20260709/5-6dogs/2.gwas-finemap/work/1.2_build_pnas_dog_behavior_locus_manifest/20260801_202130"
CANDIDATE_TABLE="${STEP12_RUN_DIR}/gwas_p_lrt_candidate_variants.tsv"
TRAIT_METADATA="${STEP12_RUN_DIR}/trait_metadata.tsv"
STEP12_COMPLETION="${STEP12_RUN_DIR}/COMPLETED.ok"
CANDIDATE_PRODUCER="/path/to/repository/scripts/11_gwas_finemapping/1.2_build_pnas_dog_behavior_locus_manifest/1.2_build_pnas_dog_behavior_locus_manifest.py"

ORIGINAL_LD_PREFIX="/path/to/project/legacy/reference/gwas/LD/GRLS_Genotyping_Axiom_Canine_Merged_2023_08_02/AxiomGT1v2.merged"
STEP15_RUN_DIR="/path/to/project/analysis_20260709/5-6dogs/2.gwas-finemap/work/1.5_prepare_grls_pilot_ld_reference_plink2_king_primary/20260802_134457"
FINAL_LD_PREFIX="${STEP15_RUN_DIR}/ld_reference_qc_plink2_king_primary"
PRELIMINARY_MAF_REMOVE="${STEP15_RUN_DIR}/ld_variant_remove.preliminary_maf.txt"
PRELIMINARY_MISSINGNESS_REMOVE="${STEP15_RUN_DIR}/ld_variant_remove.preliminary_missingness.txt"
FINAL_QC_REASONS="${STEP15_RUN_DIR}/ld_plink2_king_primary_variant_remove.final_qc_reasons.tsv"
STEP15_COMPLETION="${STEP15_RUN_DIR}/COMPLETED.ok"
LD_PRODUCER="/path/to/repository/scripts/11_gwas_finemapping/1.5_prepare_grls_pilot_ld_reference_plink2_king_primary/1.5_prepare_grls_pilot_ld_reference_plink2_king_primary.py"

OUTPUT_ROOT="/path/to/project/analysis_20260709/5-6dogs/2.gwas-finemap/work/1.6_audit_dog_behavior_candidate_ld_compatibility"
RUN_ID="${RUN_ID:-$(date '+%Y%m%d_%H%M%S')}"
RUN_DIR="${OUTPUT_ROOT}/${RUN_ID}"

AUDIT_DIR="${SCRIPT_DIR}/log"
LOCK_FILE="${AUDIT_DIR}/1.6_audit_dog_behavior_candidate_ld_compatibility.lock"
MASTER_LOG="${AUDIT_DIR}/1.6_audit_dog_behavior_candidate_ld_compatibility.${RUN_ID}.log"
SUCCESS_FILE="${AUDIT_DIR}/1.6_audit_dog_behavior_candidate_ld_compatibility.${RUN_ID}.success.tsv"
FAILED_FILE="${AUDIT_DIR}/1.6_audit_dog_behavior_candidate_ld_compatibility.${RUN_ID}.failed.tsv"

TASK_ID="audit_dog_behavior_candidate_ld_compatibility"
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
    printf '%s RUN_ID format invalid, expected YYYYMMDD_HHMMSS: %s\n' \
        "$(timestamp)" "${RUN_ID}" >&2
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
log "candidate site table: ${CANDIDATE_TABLE}"
log "original LD prefix: ${ORIGINAL_LD_PREFIX}"
log "final main LD prefix: ${FINAL_LD_PREFIX}"
log "output directory: ${RUN_DIR}"
log "standard log: ${MASTER_LOG}"

required_nonempty=(
    "${PYTHON_BIN}"
    "${WORKER}"
    "${CANDIDATE_TABLE}"
    "${TRAIT_METADATA}"
    "${STEP12_COMPLETION}"
    "${CANDIDATE_PRODUCER}"
    "${ORIGINAL_LD_PREFIX}.bed"
    "${ORIGINAL_LD_PREFIX}.bim"
    "${ORIGINAL_LD_PREFIX}.fam"
    "${FINAL_LD_PREFIX}.bed"
    "${FINAL_LD_PREFIX}.bim"
    "${FINAL_LD_PREFIX}.fam"
    "${PRELIMINARY_MAF_REMOVE}"
    "${FINAL_QC_REASONS}"
    "${STEP15_COMPLETION}"
    "${LD_PRODUCER}"
)
for required_path in "${required_nonempty[@]}"; do
    if [[ ! -s "${required_path}" ]]; then
        log "required input does not exist or is empty: ${required_path}"
        false
    fi
done
if [[ ! -f "${PRELIMINARY_MISSINGNESS_REMOVE}" ]]; then
    log "initial screening missingness remove list does not exist: ${PRELIMINARY_MISSINGNESS_REMOVE}"
    false
fi

set +e
"${PYTHON_BIN}" "${WORKER}" \
    --candidate-table "${CANDIDATE_TABLE}" \
    --trait-metadata "${TRAIT_METADATA}" \
    --step12-completion "${STEP12_COMPLETION}" \
    --candidate-producer "${CANDIDATE_PRODUCER}" \
    --original-ld-prefix "${ORIGINAL_LD_PREFIX}" \
    --final-ld-prefix "${FINAL_LD_PREFIX}" \
    --preliminary-maf-remove "${PRELIMINARY_MAF_REMOVE}" \
    --preliminary-missingness-remove "${PRELIMINARY_MISSINGNESS_REMOVE}" \
    --final-qc-reasons "${FINAL_QC_REASONS}" \
    --step15-completion "${STEP15_COMPLETION}" \
    --ld-producer "${LD_PRODUCER}" \
    --output-dir "${RUN_DIR}"
WORKER_EXIT=$?
set -e

SUMMARY_FILE="${RUN_DIR}/candidate_ld_compatibility_summary.tsv"
COMPATIBILITY_FILE="${RUN_DIR}/candidate_ld_compatibility.tsv"
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

if [[ "${WORKER_EXIT}" -eq 0 && "${SUMMARY_OK}" -eq 1 && -s "${COMPATIBILITY_FILE}" && -s "${COMPLETION_MARKER}" ]]; then
    record_status "success" "0" "candidate_ld_compatibility_audit_completed"
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


