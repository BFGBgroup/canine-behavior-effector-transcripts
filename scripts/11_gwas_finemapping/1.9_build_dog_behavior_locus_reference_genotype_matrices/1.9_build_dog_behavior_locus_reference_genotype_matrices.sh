#!/usr/bin/env bash


set -Eeuo pipefail

SCRIPT_DIR="/path/to/repository/scripts/11_gwas_finemapping/1.9_build_dog_behavior_locus_reference_genotype_matrices"
WORKER_R="${SCRIPT_DIR}/1.9_build_dog_behavior_locus_reference_genotype_matrices.R"
LOG_DIR="${SCRIPT_DIR}/log"

RSCRIPT_BIN="/path/to/software/environments/plasma_py27/bin/Rscript"
PLINK2_BIN="/path/to/software/miniconda3/envs/dogs-harmony-latest/bin/plink2"

STEP18_RUN="/path/to/project/analysis_20260709/5-6dogs/2.gwas-finemap/work/1.8_harmonize_dog_behavior_gwas_ld_locus_variants/20260802_175127"
LOCUS_MANIFEST="${STEP18_RUN}/locus_input_manifest.tsv"
STEP18_COMPLETION="${STEP18_RUN}/COMPLETED.ok"
STEP18_PRODUCER="/path/to/repository/scripts/11_gwas_finemapping/1.8_harmonize_dog_behavior_gwas_ld_locus_variants/1.8_harmonize_dog_behavior_gwas_ld_locus_variants.py"

STEP15_RUN="/path/to/project/analysis_20260709/5-6dogs/2.gwas-finemap/work/1.5_prepare_grls_pilot_ld_reference_plink2_king_primary/20260802_134457"
LD_PREFIX="${STEP15_RUN}/ld_reference_qc_plink2_king_primary"
STEP15_COMPLETION="${STEP15_RUN}/COMPLETED.ok"
STEP15_PRODUCER="/path/to/repository/scripts/11_gwas_finemapping/1.5_prepare_grls_pilot_ld_reference_plink2_king_primary/1.5_prepare_grls_pilot_ld_reference_plink2_king_primary.py"

OUTPUT_ROOT="/path/to/project/analysis_20260709/5-6dogs/2.gwas-finemap/work/1.9_build_dog_behavior_locus_reference_genotype_matrices"
MAX_PARALLEL="${MAX_PARALLEL:-4}"
PLINK2_THREADS_PER_TASK="${PLINK2_THREADS_PER_TASK:-1}"
FREQUENCY_TOLERANCE="${FREQUENCY_TOLERANCE:-1e-5}"

RUN_ID="$(date '+%Y%m%d_%H%M%S')"
RUN_DIR="${OUTPUT_ROOT}/${RUN_ID}"
MASTER_LOG="${LOG_DIR}/1.9_build_dog_behavior_locus_reference_genotype_matrices.${RUN_ID}.log"
SUCCESS_TSV="${LOG_DIR}/1.9_build_dog_behavior_locus_reference_genotype_matrices.${RUN_ID}.success.tsv"
FAILED_TSV="${LOG_DIR}/1.9_build_dog_behavior_locus_reference_genotype_matrices.${RUN_ID}.failed.tsv"
SKIPPED_TSV="${LOG_DIR}/1.9_build_dog_behavior_locus_reference_genotype_matrices.${RUN_ID}.skipped.tsv"
LOCK_FILE="${LOG_DIR}/1.9_build_dog_behavior_locus_reference_genotype_matrices.lock"
STATUS_LOCK="${LOG_DIR}/1.9_build_dog_behavior_locus_reference_genotype_matrices.${RUN_ID}.status.lock"

mkdir -p "${LOG_DIR}"
exec 9>"${LOCK_FILE}"
if ! flock -n 9; then
  printf '[%s] [ERROR] another task for this step is already running, refuse duplicate start\n' "$(date '+%F %T')" >&2
  exit 1
fi

if [[ -e "${RUN_DIR}" ]]; then
  printf '[%s] [ERROR] output directory already exists, refuse to overwrite: %s\n' "$(date '+%F %T')" "${RUN_DIR}" >&2
  exit 1
fi

mkdir -p "${RUN_DIR}/task_logs"
exec >"${MASTER_LOG}" 2>&1

log_message() {
  printf '[%s ] [%s] %s\n' "$(date '+%F %T')" "$1" "$2"
}

on_error() {
  local exit_code=$?
  local line_no=$1
  log_message "ERROR" "Launcher failed at line ${line_no} with exit code ${exit_code}."
  exit "${exit_code}"
}
trap 'on_error ${LINENO}' ERR

STATUS_HEADER=$'run_id\tlocus_id\ttrait\tanalysis_scope\tstatus\texit_code\tstarted_at\tfinished_at\toutput_path\ttask_log\treason'
printf '%s\n' "${STATUS_HEADER}" >"${SUCCESS_TSV}"
printf '%s\n' "${STATUS_HEADER}" >"${FAILED_TSV}"
printf '%s\n' "${STATUS_HEADER}" >"${SKIPPED_TSV}"

log_message "START" "Step 1.9 start, run_id=${RUN_ID}"
log_message "CONFIG" "MAX_PARALLEL=${MAX_PARALLEL}, PLINK2_THREADS_PER_TASK=${PLINK2_THREADS_PER_TASK}"

[[ "${MAX_PARALLEL}" =~ ^[1-9][0-9]*$ ]] || { log_message "ERROR" "MAX_PARALLEL must be positive integer"; exit 1; }
[[ "${PLINK2_THREADS_PER_TASK}" =~ ^[1-9][0-9]*$ ]] || { log_message "ERROR" "PLINK2_THREADS_PER_TASK must be positive integer"; exit 1; }

for required_file in \
  "${WORKER_R}" \
  "${RSCRIPT_BIN}" \
  "${PLINK2_BIN}" \
  "${LOCUS_MANIFEST}" \
  "${STEP18_COMPLETION}" \
  "${STEP18_PRODUCER}" \
  "${LD_PREFIX}.bed" \
  "${LD_PREFIX}.bim" \
  "${LD_PREFIX}.fam" \
  "${STEP15_COMPLETION}" \
  "${STEP15_PRODUCER}"; do
  if [[ ! -e "${required_file}" ]]; then
    log_message "ERROR" "required input does not exist: ${required_file}"
    exit 1
  fi
done

STEP15_RECORDED_PREFIX="$(awk -F'\t' '$1 == "primary_prefix" {sub(/\r$/, "", $2); print $2}' "${STEP15_COMPLETION}")"
STEP15_RECORDED_SAMPLES="$(awk -F'\t' '$1 == "primary_samples" {sub(/\r$/, "", $2); print $2}' "${STEP15_COMPLETION}")"
STEP15_RECORDED_VARIANTS="$(awk -F'\t' '$1 == "primary_variants" {sub(/\r$/, "", $2); print $2}' "${STEP15_COMPLETION}")"
if [[ "${STEP15_RECORDED_PREFIX}" != "${LD_PREFIX}" ]] \
  || [[ ! "${STEP15_RECORDED_SAMPLES}" =~ ^[1-9][0-9]*$ ]] \
  || [[ ! "${STEP15_RECORDED_VARIANTS}" =~ ^[1-9][0-9]*$ ]]; then
  log_message "ERROR" "Step 1.5 completion marker does not match the current LD reference panel"
  exit 1
fi
BED_VARIANT_COUNT="$(wc -l <"${LD_PREFIX}.bim")"
BED_SAMPLE_COUNT="$(wc -l <"${LD_PREFIX}.fam")"
BED_BYTE_COUNT="$(stat -c '%s' "${LD_PREFIX}.bed")"
BED_BYTES_PER_VARIANT=$(( (BED_SAMPLE_COUNT + 3) / 4 ))
BED_EXPECTED_BYTE_COUNT=$(( 3 + BED_BYTES_PER_VARIANT * BED_VARIANT_COUNT ))
BED_MAGIC="$(od -An -tx1 -N3 "${LD_PREFIX}.bed" | tr -d '[:space:]')"
if [[ "${BED_MAGIC}" != "6c1b01" ]] \
  || (( BED_BYTE_COUNT != BED_EXPECTED_BYTE_COUNT )) \
  || (( BED_SAMPLE_COUNT != STEP15_RECORDED_SAMPLES )) \
  || (( BED_VARIANT_COUNT != STEP15_RECORDED_VARIANTS )); then
  log_message "ERROR" "BED header, expected file size or BIM/FAM count differs from the Step 1.5 completion record"
  exit 1
fi
log_message "INPUT" "BED validate passed: ${BED_SAMPLE_COUNT} samples, ${BED_VARIANT_COUNT} variants, ${BED_BYTE_COUNT} bytes"
grep -Eq $'^overall_status\t(PASS|PASS_WITH_WARNINGS)$' "${STEP18_COMPLETION}" || {
  log_message "ERROR" "Step 1.8 completion marker is not PASS/PASS_WITH_WARNINGS"
  exit 1
}

EXPECTED_HEADER=$'locus_id\ttrait\tanalysis_scope\tlocus_input_status\tharmonized_variant_count\tn_locus_median\tlocus_directory\tharmonized_variants\tvariant_order\tallele_harmonization_audit\tdropped_variants\tlocus_metadata'
ACTUAL_HEADER="$(head -n 1 "${LOCUS_MANIFEST}" | tr -d '\r')"
if [[ "${ACTUAL_HEADER}" != "${EXPECTED_HEADER}" ]]; then
  log_message "ERROR" "locus_input_manifest.tsv field order differs from the expected order"
  exit 1
fi

TOTAL_TASKS="$(( $(wc -l <"${LOCUS_MANIFEST}") - 1 ))"
if (( TOTAL_TASKS < 1 )); then
  log_message "ERROR" "locus_input_manifest.tsv no task"
  exit 1
fi
log_message "INPUT" "discovery ${TOTAL_TASKS} items locus; LD reference panel prefix=${LD_PREFIX}"

{
  printf 'input_role\tabsolute_path\tproducer_script\n'
  printf 'locus_manifest\t%s\t%s\n' "${LOCUS_MANIFEST}" "${STEP18_PRODUCER}"
  printf 'step18_completion\t%s\t%s\n' "${STEP18_COMPLETION}" "${STEP18_PRODUCER}"
  printf 'ld_reference_bed\t%s\t%s\n' "${LD_PREFIX}.bed" "${STEP15_PRODUCER}"
  printf 'ld_reference_bim\t%s\t%s\n' "${LD_PREFIX}.bim" "${STEP15_PRODUCER}"
  printf 'ld_reference_fam\t%s\t%s\n' "${LD_PREFIX}.fam" "${STEP15_PRODUCER}"
  printf 'step15_completion\t%s\t%s\n' "${STEP15_COMPLETION}" "${STEP15_PRODUCER}"
} >"${RUN_DIR}/input_lineage.tsv"

{
  printf 'parameter\tvalue\n'
  printf 'run_id\t%s\n' "${RUN_ID}"
  printf 'max_parallel\t%s\n' "${MAX_PARALLEL}"
  printf 'plink2_threads_per_task\t%s\n' "${PLINK2_THREADS_PER_TASK}"
  printf 'frequency_tolerance\t%s\n' "${FREQUENCY_TOLERANCE}"
  printf 'plink2_species_mode\tdog_38_autosomes_no_xy\n'
  printf 'plink2_allow_extra_chr\tyes\n'
  printf 'plink2_export\tA include-alt\n'
  printf 'plink2_export_counted_allele\tREF_equals_BIM_A2\n'
  printf 'matrix_counted_allele\tBIM_A1_after_2_minus_REF_dosage\n'
  printf 'missing_genotype_policy\tper_variant_mean_imputation\n'
  printf 'zero_variance_policy\tdrop_before_LD\n'
  printf 'main_branch\tall_harmonized_variants\n'
  printf 'strict_branch\texclude_ambiguous_AT_CG\n'
} >"${RUN_DIR}/pipeline_parameters.tsv"

append_status() {
  local destination=$1
  local row=$2
  {
    flock -x 200
    printf '%b\n' "${row}" >>"${destination}"
  } 200>"${STATUS_LOCK}"
}

run_locus() {
  local task_index=$1
  local locus_id=$2
  local trait=$3
  local analysis_scope=$4
  local started_at
  local finished_at
  local exit_code
  local locus_output="${RUN_DIR}/${analysis_scope}/${trait}/${locus_id}"
  local task_log="${RUN_DIR}/task_logs/${analysis_scope}/${trait}/${locus_id}.log"

  started_at="$(date '+%F %T')"
  mkdir -p "$(dirname "${locus_output}")" "$(dirname "${task_log}")"


  if "${RSCRIPT_BIN}" "${WORKER_R}" \
      --locus-id "${locus_id}" \
      --manifest "${LOCUS_MANIFEST}" \
      --ld-prefix "${LD_PREFIX}" \
      --plink2 "${PLINK2_BIN}" \
      --plink2-threads "${PLINK2_THREADS_PER_TASK}" \
      --frequency-tolerance "${FREQUENCY_TOLERANCE}" \
      --output-dir "${locus_output}" \
      >"${task_log}" 2>&1; then
    exit_code=0
  else
    exit_code=$?
  fi
  finished_at="$(date '+%F %T')"

  if (( exit_code == 0 )) \
    && [[ -s "${locus_output}/COMPLETED.ok" ]] \
    && grep -q $'^overall_status\tPASS' "${locus_output}/COMPLETED.ok" \
    && [[ -s "${locus_output}/main/X_ref.rds" ]] \
    && [[ -s "${locus_output}/main/R_signed.rds" ]] \
    && [[ -s "${locus_output}/strict_no_ambiguous/X_ref.rds" ]] \
    && [[ -s "${locus_output}/strict_no_ambiguous/R_signed.rds" ]] \
    && [[ -s "${locus_output}/matrix_qc.tsv" ]]; then
    append_status "${SUCCESS_TSV}" "${RUN_ID}\t${locus_id}\t${trait}\t${analysis_scope}\tSUCCESS\t0\t${started_at}\t${finished_at}\t${locus_output}\t${task_log}\tPASS"
    log_message "PROGRESS" "completed ${task_index}/${TOTAL_TASKS}: ${locus_id}"
  else
    append_status "${FAILED_TSV}" "${RUN_ID}\t${locus_id}\t${trait}\t${analysis_scope}\tFAILED\t${exit_code}\t${started_at}\t${finished_at}\t${locus_output}\t${task_log}\tWORKER_OR_OUTPUT_VALIDATION_FAILED"
    log_message "PROGRESS" "failed ${task_index}/${TOTAL_TASKS}: ${locus_id}, see ${task_log}"
  fi
  return 0
}

active_jobs=0
dispatched_tasks=0
skipped_tasks=0
seen_tasks=0

while IFS=$'\t' read -r locus_id trait analysis_scope locus_input_status _remaining; do
  seen_tasks=$((seen_tasks + 1))
  if [[ "${locus_input_status}" != "READY_FOR_GENOTYPE_MATRIX" ]]; then
    now="$(date '+%F %T')"
    append_status "${SKIPPED_TSV}" "${RUN_ID}\t${locus_id}\t${trait}\t${analysis_scope}\tSKIPPED\t0\t${now}\t${now}\tNA\tNA\t${locus_input_status}"
    skipped_tasks=$((skipped_tasks + 1))
    log_message "PROGRESS" "skip ${seen_tasks}/${TOTAL_TASKS}: ${locus_id}, reason=${locus_input_status}"
    continue
  fi

  run_locus "${seen_tasks}" "${locus_id}" "${trait}" "${analysis_scope}" &
  dispatched_tasks=$((dispatched_tasks + 1))
  active_jobs=$((active_jobs + 1))
  if (( active_jobs >= MAX_PARALLEL )); then
    wait -n
    active_jobs=$((active_jobs - 1))
  fi
done < <(tail -n +2 "${LOCUS_MANIFEST}")

while (( active_jobs > 0 )); do
  wait -n
  active_jobs=$((active_jobs - 1))
done

success_count="$(( $(wc -l <"${SUCCESS_TSV}") - 1 ))"
failed_count="$(( $(wc -l <"${FAILED_TSV}") - 1 ))"
skipped_count="$(( $(wc -l <"${SKIPPED_TSV}") - 1 ))"
terminal_count=$((success_count + failed_count + skipped_count))

if (( terminal_count != TOTAL_TASKS )); then
  log_message "ERROR" "Terminal task count ${terminal_count} differs from the input task count ${TOTAL_TASKS}"
  exit 1
fi

LOCUS_SUMMARY="${RUN_DIR}/locus_matrix_manifest.tsv"
QC_SUMMARY="${RUN_DIR}/matrix_qc_summary.tsv"
printf 'locus_id\ttrait\tanalysis_scope\tlocus_status\tsample_count\tmain_variant_count\tstrict_variant_count\tlocus_output\tmain_x_ref\tmain_r_signed\tmain_z_beta_se\tstrict_x_ref\tstrict_r_signed\tstrict_z_beta_se\n' >"${LOCUS_SUMMARY}"
printf 'locus_id\ttrait\tanalysis_scope\tbranch\tstatus\tsample_count\tinput_variant_count\tzero_variance_dropped_count\toutput_variant_count\tmissing_genotype_count_before_imputation\tmissing_genotype_fraction_before_imputation\tmissing_genotype_count_after_imputation\tmax_a1_frequency_difference\tr_symmetry_max_abs\tr_diagonal_max_abs_difference\tr_min_eigenvalue\tn_locus_median\tn_obs_min\tn_obs_max\n' >"${QC_SUMMARY}"

while IFS=$'\t' read -r _run_id _locus_id _trait _scope _status _exit_code _started _finished locus_output _task_log _reason; do
  tail -n +2 "${locus_output}/locus_matrix_manifest.tsv" >>"${LOCUS_SUMMARY}"
  tail -n +2 "${locus_output}/matrix_qc.tsv" >>"${QC_SUMMARY}"
done < <(tail -n +2 "${SUCCESS_TSV}")

main_total="$(awk -F'\t' 'NR>1 {sum += $6} END {print sum+0}' "${LOCUS_SUMMARY}")"
strict_total="$(awk -F'\t' 'NR>1 {sum += $7} END {print sum+0}' "${LOCUS_SUMMARY}")"

{
  printf 'metric\tvalue\n'
  printf 'input_locus_count\t%s\n' "${TOTAL_TASKS}"
  printf 'dispatched_locus_count\t%s\n' "${dispatched_tasks}"
  printf 'success_locus_count\t%s\n' "${success_count}"
  printf 'failed_locus_count\t%s\n' "${failed_count}"
  printf 'skipped_locus_count\t%s\n' "${skipped_count}"
  printf 'main_total_locus_variant_count\t%s\n' "${main_total}"
  printf 'strict_total_locus_variant_count\t%s\n' "${strict_total}"
} >"${RUN_DIR}/overall_summary.tsv"

if (( failed_count > 0 )); then
  log_message "FAILED" "${failed_count} items locus failed; see ${FAILED_TSV}"
  exit 1
fi

overall_status="PASS"
if (( skipped_count > 0 )); then
  overall_status="PASS_WITH_WARNINGS"
fi

{
  printf 'overall_status\t%s\n' "${overall_status}"
  printf 'run_id\t%s\n' "${RUN_ID}"
  printf 'input_locus_count\t%s\n' "${TOTAL_TASKS}"
  printf 'success_locus_count\t%s\n' "${success_count}"
  printf 'failed_locus_count\t%s\n' "${failed_count}"
  printf 'skipped_locus_count\t%s\n' "${skipped_count}"
  printf 'completed_at\t%s\n' "$(date '+%F %T')"
} >"${RUN_DIR}/COMPLETED.ok"

log_message "DONE" "Step 1.9 completed; status=${overall_status}, success=${success_count}, skip=${skipped_count}"
log_message "OUTPUT" "run directory: ${RUN_DIR}"
