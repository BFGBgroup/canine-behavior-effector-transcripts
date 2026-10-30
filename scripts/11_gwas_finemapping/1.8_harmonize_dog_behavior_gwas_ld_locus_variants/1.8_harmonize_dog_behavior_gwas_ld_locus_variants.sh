#!/usr/bin/env bash


set -Eeuo pipefail

SCRIPT_DIR="/path/to/repository/scripts/11_gwas_finemapping/1.8_harmonize_dog_behavior_gwas_ld_locus_variants"
PYTHON_BIN="/path/to/software/miniconda3/envs/dogs-harmony-latest/bin/python"
WORKER="${SCRIPT_DIR}/1.8_harmonize_dog_behavior_gwas_ld_locus_variants.py"

STEP11_RUN_DIR="/path/to/project/analysis_20260709/5-6dogs/2.gwas-finemap/work/1.1_validate_dog_behavior_gwas_finemap_inputs/20260801_182359"
TRAIT_INPUT_CONTRACT="${STEP11_RUN_DIR}/trait_input_contract.tsv"
STEP11_PRODUCER="/path/to/repository/scripts/11_gwas_finemapping/1.1_validate_dog_behavior_gwas_finemap_inputs/1.1_validate_dog_behavior_gwas_finemap_inputs.py"

STEP12_RUN_DIR="/path/to/project/analysis_20260709/5-6dogs/2.gwas-finemap/work/1.2_build_pnas_dog_behavior_locus_manifest/20260801_202130"
TRAIT_METADATA="${STEP12_RUN_DIR}/trait_metadata.tsv"
STEP12_PRODUCER="/path/to/repository/scripts/11_gwas_finemapping/1.2_build_pnas_dog_behavior_locus_manifest/1.2_build_pnas_dog_behavior_locus_manifest.py"

STEP15_RUN_DIR="/path/to/project/analysis_20260709/5-6dogs/2.gwas-finemap/work/1.5_prepare_grls_pilot_ld_reference_plink2_king_primary/20260802_134457"
FINAL_LD_PREFIX="${STEP15_RUN_DIR}/ld_reference_qc_plink2_king_primary"
FINAL_LD_FREQUENCY="${STEP15_RUN_DIR}/intermediate/11_joint_round_02_final_validation.afreq"
STEP15_COMPLETION="${STEP15_RUN_DIR}/COMPLETED.ok"
STEP15_PRODUCER="/path/to/repository/scripts/11_gwas_finemapping/1.5_prepare_grls_pilot_ld_reference_plink2_king_primary/1.5_prepare_grls_pilot_ld_reference_plink2_king_primary.py"

STEP17_RUN_DIR="/path/to/project/analysis_20260709/5-6dogs/2.gwas-finemap/work/1.7_define_dog_behavior_ld_loci/20260802_171634"
LOCUS_REGIONS="${STEP17_RUN_DIR}/all_locus_regions.tsv"
NON_FINEMAPPABLE_CANDIDATES="${STEP17_RUN_DIR}/non_finemappable_candidate_regions.tsv"
STEP17_COMPLETION="${STEP17_RUN_DIR}/COMPLETED.ok"
STEP17_PRODUCER="/path/to/repository/scripts/11_gwas_finemapping/1.7_define_dog_behavior_ld_loci/1.7_define_dog_behavior_ld_loci.py"

OUTPUT_ROOT="/path/to/project/analysis_20260709/5-6dogs/2.gwas-finemap/work/1.8_harmonize_dog_behavior_gwas_ld_locus_variants"
MAX_MISSING_RATE="0.05"
AMBIGUOUS_FREQUENCY_TOLERANCE="0.10"
AMBIGUOUS_FREQUENCY_SEPARATION="0.10"
MIN_VARIANTS_FOR_FINEMAPPING="2"
PROGRESS_EVERY="100000"
RUN_ID="${RUN_ID:-$(date '+%Y%m%d_%H%M%S')}"
RUN_DIR="${OUTPUT_ROOT}/${RUN_ID}"

AUDIT_DIR="${SCRIPT_DIR}/log"
LOCK_FILE="${AUDIT_DIR}/1.8_harmonize_dog_behavior_gwas_ld_locus_variants.lock"
MASTER_LOG="${AUDIT_DIR}/1.8_harmonize_dog_behavior_gwas_ld_locus_variants.${RUN_ID}.log"
SUCCESS_FILE="${AUDIT_DIR}/1.8_harmonize_dog_behavior_gwas_ld_locus_variants.${RUN_ID}.success.tsv"
FAILED_FILE="${AUDIT_DIR}/1.8_harmonize_dog_behavior_gwas_ld_locus_variants.${RUN_ID}.failed.tsv"

TASK_ID="harmonize_dog_behavior_gwas_ld_locus_variants"
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
log "1.7 all locus: ${LOCUS_REGIONS}"
log "complete GWAS path contract: ${TRAIT_INPUT_CONTRACT}"
log "final LD prefix: ${FINAL_LD_PREFIX}"
log "final LD A1 frequency: ${FINAL_LD_FREQUENCY}"
log "parameter: max_missing_rate=${MAX_MISSING_RATE}; ambiguous_tolerance=${AMBIGUOUS_FREQUENCY_TOLERANCE}; ambiguous_separation=${AMBIGUOUS_FREQUENCY_SEPARATION}"
log "output directory: ${RUN_DIR}"
log "standard log: ${MASTER_LOG}"

required_nonempty=(
    "${PYTHON_BIN}"
    "${WORKER}"
    "${TRAIT_INPUT_CONTRACT}"
    "${TRAIT_METADATA}"
    "${FINAL_LD_PREFIX}.bed"
    "${FINAL_LD_PREFIX}.bim"
    "${FINAL_LD_PREFIX}.fam"
    "${FINAL_LD_FREQUENCY}"
    "${STEP15_COMPLETION}"
    "${LOCUS_REGIONS}"
    "${NON_FINEMAPPABLE_CANDIDATES}"
    "${STEP17_COMPLETION}"
    "${STEP11_PRODUCER}"
    "${STEP12_PRODUCER}"
    "${STEP15_PRODUCER}"
    "${STEP17_PRODUCER}"
)
for required_path in "${required_nonempty[@]}"; do
    if [[ ! -s "${required_path}" ]]; then
        log "required input does not exist or is empty: ${required_path}"
        false
    fi
done

set +e
"${PYTHON_BIN}" "${WORKER}" \
    --locus-regions "${LOCUS_REGIONS}" \
    --non-finemappable-candidates "${NON_FINEMAPPABLE_CANDIDATES}" \
    --step17-completion "${STEP17_COMPLETION}" \
    --trait-input-contract "${TRAIT_INPUT_CONTRACT}" \
    --trait-metadata "${TRAIT_METADATA}" \
    --ld-prefix "${FINAL_LD_PREFIX}" \
    --ld-frequency "${FINAL_LD_FREQUENCY}" \
    --step15-completion "${STEP15_COMPLETION}" \
    --step11-producer "${STEP11_PRODUCER}" \
    --step12-producer "${STEP12_PRODUCER}" \
    --step15-producer "${STEP15_PRODUCER}" \
    --step17-producer "${STEP17_PRODUCER}" \
    --output-dir "${RUN_DIR}" \
    --max-missing-rate "${MAX_MISSING_RATE}" \
    --ambiguous-frequency-tolerance "${AMBIGUOUS_FREQUENCY_TOLERANCE}" \
    --ambiguous-frequency-separation "${AMBIGUOUS_FREQUENCY_SEPARATION}" \
    --min-variants-for-finemapping "${MIN_VARIANTS_FOR_FINEMAPPING}" \
    --progress-every "${PROGRESS_EVERY}"
WORKER_EXIT=$?
set -e

OVERALL_SUMMARY="${RUN_DIR}/harmonization_overall_summary.tsv"
LOCUS_SUMMARY="${RUN_DIR}/harmonization_summary.tsv"
LOCUS_MANIFEST="${RUN_DIR}/locus_input_manifest.tsv"
NON_FINEMAPPABLE_AUDIT="${RUN_DIR}/non_finemappable_candidate_carry_forward.tsv"
PARAMETERS="${RUN_DIR}/pipeline_parameters.tsv"
LINEAGE="${RUN_DIR}/input_lineage.tsv"
COMPLETION_MARKER="${RUN_DIR}/COMPLETED.ok"
SUMMARY_OK=0
if [[ -s "${OVERALL_SUMMARY}" ]]; then
    if awk -F '\t' '
        $1 == "overall_status" && ($2 == "PASS" || $2 == "PASS_WITH_WARNINGS") {found=1}
        END {exit(found ? 0 : 1)}
    ' "${OVERALL_SUMMARY}"; then
        SUMMARY_OK=1
    fi
fi

if [[ "${WORKER_EXIT}" -eq 0 && "${SUMMARY_OK}" -eq 1 \
    && -s "${LOCUS_SUMMARY}" && -s "${LOCUS_MANIFEST}" \
    && -s "${NON_FINEMAPPABLE_AUDIT}" && -s "${PARAMETERS}" \
    && -s "${LINEAGE}" && -s "${COMPLETION_MARKER}" ]]; then
    record_status "success" "0" "dog_behavior_gwas_ld_locus_harmonization_completed"
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


