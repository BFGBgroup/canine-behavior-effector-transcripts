#!/usr/bin/env bash


set -Eeuo pipefail

SCRIPT_DIR="/path/to/repository/scripts/11_gwas_finemapping/1.5_prepare_grls_pilot_ld_reference_plink2_king_primary"
WORKER="${SCRIPT_DIR}/1.5_prepare_grls_pilot_ld_reference_plink2_king_primary.py"
INPUT_PREFIX="/path/to/project/legacy/reference/gwas/LD/GRLS_Genotyping_Axiom_Canine_Merged_2023_08_02/AxiomGT1v2.merged"
STEP11_DIR="/path/to/project/analysis_20260709/5-6dogs/2.gwas-finemap/work/1.1_validate_dog_behavior_gwas_finemap_inputs/20260801_182359"
STEP12_DIR="/path/to/project/analysis_20260709/5-6dogs/2.gwas-finemap/work/1.2_build_pnas_dog_behavior_locus_manifest/20260801_202130"
OUTPUT_ROOT="/path/to/project/analysis_20260709/5-6dogs/2.gwas-finemap/work/1.5_prepare_grls_pilot_ld_reference_plink2_king_primary"
PLINK2_BIN="/path/to/software/miniconda3/envs/dogs-harmony-latest/bin/plink2"

THREADS="${THREADS:-16}"
MIND="0.05"
GENO="0.05"
MAF="0.05"
HWE="1e-6"
PRUNE_WINDOW="500kb"
PRUNE_STEP="1"
PRUNE_R2="0.2"
KING_RELATEDNESS="0.177"
MAX_FINAL_QC_ITERATIONS="10"
RUN_ID="${RUN_ID:-$(date +%Y%m%d_%H%M%S)}"
RUN_DIR="${OUTPUT_ROOT}/${RUN_ID}"
AUDIT_DIR="${SCRIPT_DIR}/log"
MASTER_LOG="${AUDIT_DIR}/1.5_prepare_grls_pilot_ld_reference_plink2_king_primary.${RUN_ID}.log"
SUCCESS_TSV="${AUDIT_DIR}/1.5_prepare_grls_pilot_ld_reference_plink2_king_primary.${RUN_ID}.success.tsv"
FAILED_TSV="${AUDIT_DIR}/1.5_prepare_grls_pilot_ld_reference_plink2_king_primary.${RUN_ID}.failed.tsv"
LOCK_FILE="${AUDIT_DIR}/1.5_prepare_grls_pilot_ld_reference_plink2_king_primary.lock"
STARTED_AT="$(date '+%Y-%m-%d %H:%M:%S')"

mkdir -p "${AUDIT_DIR}"
exec 9>"${LOCK_FILE}"
if ! flock -n 9; then
    printf 'another task for this step is already running: %s\n' "${LOCK_FILE}" >&2
    exit 1
fi

if [[ -e "${RUN_DIR}" ]]; then
    printf 'refuse to overwrite existing run directory: %s\n' "${RUN_DIR}" >&2
    exit 1
fi

printf 'run_id\tstatus\texit_code\tstarted_at\tfinished_at\toutput_path\ttask_log\treason\n' >"${SUCCESS_TSV}"
printf 'run_id\tstatus\texit_code\tstarted_at\tfinished_at\toutput_path\ttask_log\treason\n' >"${FAILED_TSV}"
exec >>"${MASTER_LOG}" 2>&1

record_status() {
    local status="$1"
    local exit_code="$2"
    local message="$3"
    local status_file="${FAILED_TSV}"
    if [[ "${status}" == "SUCCESS" ]]; then
        status_file="${SUCCESS_TSV}"
    fi
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "${RUN_ID}" "${status}" "${exit_code}" "${STARTED_AT}" \
        "$(date '+%Y-%m-%d %H:%M:%S')" "${RUN_DIR}" "${MASTER_LOG}" "${message}" \
        >>"${status_file}"
}

on_error() {
    local exit_code=$?
    record_status "FAILED" "${exit_code}" "worker_or_launcher_failed"
    printf '[%s ] FAILED exit=%s run_dir=%s\n' "$(date '+%F %T')" "${exit_code}" "${RUN_DIR}"
    exit "${exit_code}"
}
trap on_error ERR

printf '[%s ] START run_id=%s\n' "$(date '+%F %T')" "${RUN_ID}"
printf 'worker=%s\ninput_prefix=%s\nstep11_dir=%s\nstep12_dir=%s\noutput_dir=%s\nplink2=%s\nthreads=%s\nmind=%s\ngeno=%s\nmaf=%s\nhwe=%s\nprune_window=%s\nprune_step=%s\nprune_r2=%s\nking_relatedness=%s\nmax_final_qc_iterations=%s\n' \
    "${WORKER}" "${INPUT_PREFIX}" "${STEP11_DIR}" "${STEP12_DIR}" "${RUN_DIR}" \
    "${PLINK2_BIN}" "${THREADS}" "${MIND}" "${GENO}" "${MAF}" "${HWE}" \
    "${PRUNE_WINDOW}" "${PRUNE_STEP}" "${PRUNE_R2}" "${KING_RELATEDNESS}" \
    "${MAX_FINAL_QC_ITERATIONS}"

python3 "${WORKER}" \
    --input-prefix "${INPUT_PREFIX}" \
    --step11-dir "${STEP11_DIR}" \
    --step12-dir "${STEP12_DIR}" \
    --plink2-bin "${PLINK2_BIN}" \
    --output-dir "${RUN_DIR}" \
    --threads "${THREADS}" \
    --mind "${MIND}" \
    --geno "${GENO}" \
    --maf "${MAF}" \
    --hwe "${HWE}" \
    --prune-window "${PRUNE_WINDOW}" \
    --prune-step "${PRUNE_STEP}" \
    --prune-r2 "${PRUNE_R2}" \
    --king-relatedness "${KING_RELATEDNESS}" \
    --max-final-qc-iterations "${MAX_FINAL_QC_ITERATIONS}"

SUMMARY="${RUN_DIR}/ld_reference_qc_summary.tsv"
COMPLETION="${RUN_DIR}/COMPLETED.ok"
[[ -s "${SUMMARY}" ]]
[[ -s "${COMPLETION}" ]]
awk -F '\t' '$1 == "overall_status" && $2 == "PASS" {found=1} END {exit(found ? 0 : 1)}' "${SUMMARY}"

record_status "SUCCESS" "0" "plink2_king_primary_ld_reference_qc_completed"
trap - ERR
printf '[%s ] done run_dir=%s\n' "$(date '+%F %T')" "${RUN_DIR}"


