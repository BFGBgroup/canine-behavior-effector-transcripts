#!/usr/bin/env bash


set -Eeuo pipefail

SCRIPT_DIR="/path/to/repository/scripts/11_gwas_finemapping/1.10_run_dog_behavior_susie_rss_finemapping"
WORKER_R="${SCRIPT_DIR}/1.10_run_dog_behavior_susie_rss_finemapping.R"
LOG_DIR="${SCRIPT_DIR}/log"

RSCRIPT_BIN="/path/to/software/environments/plasma_py27/bin/Rscript"
EXPECTED_SUSIER_VERSION="0.16.5"

STEP19_RUN="/path/to/project/analysis_20260709/5-6dogs/2.gwas-finemap/work/1.9_build_dog_behavior_locus_reference_genotype_matrices/20260802_194132"
LOCUS_MANIFEST="${STEP19_RUN}/locus_matrix_manifest.tsv"
QC_SUMMARY="${STEP19_RUN}/matrix_qc_summary.tsv"
STEP19_COMPLETION="${STEP19_RUN}/COMPLETED.ok"
STEP19_LAUNCHER="/path/to/repository/scripts/11_gwas_finemapping/1.9_build_dog_behavior_locus_reference_genotype_matrices/1.9_build_dog_behavior_locus_reference_genotype_matrices.sh"
STEP19_WORKER="/path/to/repository/scripts/11_gwas_finemapping/1.9_build_dog_behavior_locus_reference_genotype_matrices/1.9_build_dog_behavior_locus_reference_genotype_matrices.R"
STEP18_RUN="/path/to/project/analysis_20260709/5-6dogs/2.gwas-finemap/work/1.8_harmonize_dog_behavior_gwas_ld_locus_variants/20260802_175127"
STEP18_LOCUS_MANIFEST="${STEP18_RUN}/locus_input_manifest.tsv"
STEP18_COMPLETION="${STEP18_RUN}/COMPLETED.ok"
STEP18_PRODUCER="/path/to/repository/scripts/11_gwas_finemapping/1.8_harmonize_dog_behavior_gwas_ld_locus_variants/1.8_harmonize_dog_behavior_gwas_ld_locus_variants.py"

OUTPUT_ROOT="/path/to/project/analysis_20260709/5-6dogs/2.gwas-finemap/work/1.10_run_dog_behavior_susie_rss_finemapping"
RUN_ID="${RUN_ID:-$(date '+%Y%m%d_%H%M%S')}"
RUN_DIR="${OUTPUT_ROOT}/${RUN_ID}"

BRANCHES="${BRANCHES:-main,strict_no_ambiguous}"
ANALYSIS_SCOPE="${ANALYSIS_SCOPE:-all}"
LOCUS_ID_FILTER="${LOCUS_ID_FILTER:-}"
MAX_PARALLEL="${MAX_PARALLEL:-4}"

MASTER_LOG="${LOG_DIR}/1.10_run_dog_behavior_susie_rss_finemapping.${RUN_ID}.log"
SUCCESS_TSV="${LOG_DIR}/1.10_run_dog_behavior_susie_rss_finemapping.${RUN_ID}.success.tsv"
FAILED_TSV="${LOG_DIR}/1.10_run_dog_behavior_susie_rss_finemapping.${RUN_ID}.failed.tsv"
SKIPPED_TSV="${LOG_DIR}/1.10_run_dog_behavior_susie_rss_finemapping.${RUN_ID}.skipped.tsv"
LOCK_FILE="${LOG_DIR}/1.10_run_dog_behavior_susie_rss_finemapping.lock"
STATUS_LOCK="${LOG_DIR}/1.10_run_dog_behavior_susie_rss_finemapping.${RUN_ID}.status.lock"

mkdir -p "${LOG_DIR}"
exec 9>"${LOCK_FILE}"
if ! flock -n 9; then
  printf '[%s ] [ERROR] another task for this step is already running, refuse duplicate start\n' "$(date '+%F %T')" >&2
  exit 1
fi

if [[ -e "${RUN_DIR}" ]]; then
  printf '[%s ] [ERROR] output directory already exists, refuse to overwrite: %s\n' "$(date '+%F %T')" "${RUN_DIR}" >&2
  exit 1
fi

mkdir -p "${RUN_DIR}/task_logs"
: >"${MASTER_LOG}"
exec >>"${MASTER_LOG}" 2>&1

log_message() {
  printf '[%s ] [%s] %s\n' "$(date '+%F %T')" "$1" "$2"
}

on_error() {
  local exit_code=$?
  local line_number=$1
  log_message "ERROR" "Launcher failed at line ${line_number} with exit code ${exit_code}"
  exit "${exit_code}"
}
trap 'on_error ${LINENO}' ERR

STATUS_HEADER=$'run_id\tlocus_id\ttrait\tanalysis_scope\tbranch\tstatus\texit_code\tstarted_at\tfinished_at\toutput_path\ttask_log\treason'
printf '%s\n' "${STATUS_HEADER}" >"${SUCCESS_TSV}"
printf '%s\n' "${STATUS_HEADER}" >"${FAILED_TSV}"
printf '%s\n' "${STATUS_HEADER}" >"${SKIPPED_TSV}"

log_message "START" "Step 1.10 start, run_id=${RUN_ID}"
log_message "CONFIG" "BRANCHES=${BRANCHES}, ANALYSIS_SCOPE=${ANALYSIS_SCOPE}, LOCUS_ID_FILTER=${LOCUS_ID_FILTER:-ALL}, MAX_PARALLEL=${MAX_PARALLEL}"

case "${BRANCHES}" in
  main)
    BRANCH_ARRAY=("main")
    ;;
  strict_no_ambiguous)
    BRANCH_ARRAY=("strict_no_ambiguous")
    ;;
  main,strict_no_ambiguous)
    BRANCH_ARRAY=("main" "strict_no_ambiguous")
    ;;
  strict_no_ambiguous,main)
    BRANCH_ARRAY=("strict_no_ambiguous" "main")
    ;;
  *)
    log_message "ERROR" "BRANCHES must be one of main, strict_no_ambiguous or both of comma combination"
    exit 1
    ;;
esac
[[ "${ANALYSIS_SCOPE}" == "all" || "${ANALYSIS_SCOPE}" == "primary" || "${ANALYSIS_SCOPE}" == "exploratory" ]] || {
  log_message "ERROR" "ANALYSIS_SCOPE must be one of all, primary or exploratory"
  exit 1
}
[[ "${MAX_PARALLEL}" =~ ^[1-9][0-9]*$ ]] || {
  log_message "ERROR" "MAX_PARALLEL must be positive integer"
  exit 1
}
[[ "${RUN_ID}" =~ ^[A-Za-z0-9._-]+$ ]] || {
  log_message "ERROR" "RUN_ID can only contains letter, number, point, underscore and hyphen"
  exit 1
}

for required_path in \
  "${WORKER_R}" \
  "${RSCRIPT_BIN}" \
  "${LOCUS_MANIFEST}" \
  "${QC_SUMMARY}" \
  "${STEP19_COMPLETION}" \
  "${STEP19_LAUNCHER}" \
  "${STEP19_WORKER}" \
  "${STEP18_LOCUS_MANIFEST}" \
  "${STEP18_COMPLETION}" \
  "${STEP18_PRODUCER}"; do
  if [[ ! -s "${required_path}" ]]; then
    log_message "ERROR" "required input does not exist or is empty: ${required_path}"
    exit 1
  fi
done

grep -Eq $'^overall_status\t(PASS|PASS_WITH_WARNINGS)$' "${STEP19_COMPLETION}" || {
  log_message "ERROR" "Step 1.9 completion marker is not PASS/PASS_WITH_WARNINGS"
  exit 1
}

grep -Eq $'^overall_status\t(PASS|PASS_WITH_WARNINGS)$' "${STEP18_COMPLETION}" || {
  log_message "ERROR" "Step 1.8 completion marker is not PASS/PASS_WITH_WARNINGS"
  exit 1
}

EXPECTED_MANIFEST_HEADER=$'locus_id\ttrait\tanalysis_scope\tlocus_status\tsample_count\tmain_variant_count\tstrict_variant_count\tlocus_output\tmain_x_ref\tmain_r_signed\tmain_z_beta_se\tstrict_x_ref\tstrict_r_signed\tstrict_z_beta_se'
ACTUAL_MANIFEST_HEADER="$(head -n 1 "${LOCUS_MANIFEST}" | tr -d '\r')"
if [[ "${ACTUAL_MANIFEST_HEADER}" != "${EXPECTED_MANIFEST_HEADER}" ]]; then
  log_message "ERROR" "Step 1.9 locus manifest field order differs from the expected order"
  exit 1
fi

EXPECTED_QC_HEADER=$'locus_id\ttrait\tanalysis_scope\tbranch\tstatus\tsample_count\tinput_variant_count\tzero_variance_dropped_count\toutput_variant_count\tmissing_genotype_count_before_imputation\tmissing_genotype_fraction_before_imputation\tmissing_genotype_count_after_imputation\tmax_a1_frequency_difference\tr_symmetry_max_abs\tr_diagonal_max_abs_difference\tr_min_eigenvalue\tn_locus_median\tn_obs_min\tn_obs_max'
ACTUAL_QC_HEADER="$(head -n 1 "${QC_SUMMARY}" | tr -d '\r')"
if [[ "${ACTUAL_QC_HEADER}" != "${EXPECTED_QC_HEADER}" ]]; then
  log_message "ERROR" "Step 1.9 matrix QC field order differs from the expected order"
  exit 1
fi

EXPECTED_STEP18_HEADER=$'locus_id\ttrait\tanalysis_scope\tlocus_input_status\tharmonized_variant_count\tn_locus_median\tlocus_directory\tharmonized_variants\tvariant_order\tallele_harmonization_audit\tdropped_variants\tlocus_metadata'
ACTUAL_STEP18_HEADER="$(head -n 1 "${STEP18_LOCUS_MANIFEST}" | tr -d '\r')"
if [[ "${ACTUAL_STEP18_HEADER}" != "${EXPECTED_STEP18_HEADER}" ]]; then
  log_message "ERROR" "Step 1.8 locus manifest field order differs from the expected order"
  exit 1
fi

INPUT_LOCUS_COUNT="$(( $(wc -l <"${LOCUS_MANIFEST}") - 1 ))"
QC_ROW_COUNT="$(( $(wc -l <"${QC_SUMMARY}") - 1 ))"
if (( INPUT_LOCUS_COUNT < 1 || QC_ROW_COUNT != INPUT_LOCUS_COUNT * 2 )); then
  log_message "ERROR" "locus manifest with main/strict QC row count cannot one-to-one corresponding"
  exit 1
fi

INSTALLED_SUSIER_VERSION="$("${RSCRIPT_BIN}" --vanilla -e 'cat(as.character(utils::packageVersion("susieR")))')"
if [[ "${INSTALLED_SUSIER_VERSION}" != "${EXPECTED_SUSIER_VERSION}" ]]; then
  log_message "ERROR" "susieR version mismatch: actual=${INSTALLED_SUSIER_VERSION}; expected=${EXPECTED_SUSIER_VERSION}"
  exit 1
fi


export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

is_selected_task() {
  local locus_id=$1
  local analysis_scope=$2
  if [[ "${ANALYSIS_SCOPE}" != "all" && "${analysis_scope}" != "${ANALYSIS_SCOPE}" ]]; then
    return 1
  fi
  if [[ -n "${LOCUS_ID_FILTER}" && "${locus_id}" != "${LOCUS_ID_FILTER}" ]]; then
    return 1
  fi
  return 0
}

SELECTED_LOCUS_COUNT=0
while IFS=$'\t' read -r locus_id _trait analysis_scope _remaining; do
  if is_selected_task "${locus_id}" "${analysis_scope}"; then
    SELECTED_LOCUS_COUNT=$((SELECTED_LOCUS_COUNT + 1))
  fi
done < <(tail -n +2 "${LOCUS_MANIFEST}")

if (( SELECTED_LOCUS_COUNT < 1 )); then
  log_message "ERROR" "current ANALYSIS_SCOPE/LOCUS_ID_FILTER no matching task"
  exit 1
fi
BRANCH_COUNT=${#BRANCH_ARRAY[@]}
TOTAL_TASKS=$((SELECTED_LOCUS_COUNT * BRANCH_COUNT))
log_message "INPUT" "Step 1.9 total ${INPUT_LOCUS_COUNT} items locus; this run select ${SELECTED_LOCUS_COUNT} items locus, ${BRANCH_COUNT} branches, total ${TOTAL_TASKS} tasks"
log_message "INPUT" "locus manifest=${LOCUS_MANIFEST}"
log_message "INPUT" "matrix QC=${QC_SUMMARY}"
log_message "INPUT" "Step 1.8 locus metadata manifest=${STEP18_LOCUS_MANIFEST}"

{
  printf 'input_role\tabsolute_path\tproducer_script\n'
  printf 'locus_matrix_manifest\t%s\t%s;%s\n' "${LOCUS_MANIFEST}" "${STEP19_LAUNCHER}" "${STEP19_WORKER}"
  printf 'matrix_qc_summary\t%s\t%s;%s\n' "${QC_SUMMARY}" "${STEP19_LAUNCHER}" "${STEP19_WORKER}"
  printf 'step19_completion\t%s\t%s\n' "${STEP19_COMPLETION}" "${STEP19_LAUNCHER}"
  printf 'step18_locus_manifest\t%s\t%s\n' "${STEP18_LOCUS_MANIFEST}" "${STEP18_PRODUCER}"
  printf 'step18_completion\t%s\t%s\n' "${STEP18_COMPLETION}" "${STEP18_PRODUCER}"
} >"${RUN_DIR}/input_lineage.tsv"

{
  printf 'parameter\tvalue\n'
  printf 'run_id\t%s\n' "${RUN_ID}"
  printf 'branches\t%s\n' "${BRANCHES}"
  printf 'output_layout\tRUN_ID/branch/analysis_scope/trait/locus_id\n'
  printf 'analysis_scope\t%s\n' "${ANALYSIS_SCOPE}"
  printf 'locus_id_filter\t%s\n' "${LOCUS_ID_FILTER:-ALL}"
  printf 'max_parallel\t%s\n' "${MAX_PARALLEL}"
  printf 'susieR_version\t%s\n' "${INSTALLED_SUSIER_VERSION}"
  printf 'z_input\tz_for_ld_a1\n'
  printf 'n_input\tn_locus_median\n'
  printf 'reference_sample_count_source\tnrow_X_ref\n'
  printf 'z_method\twald\n'
  printf 'L\tmin(10,variant_count)\n'
  printf 'coverage\t0.95\n'
  printf 'min_abs_corr\t0.5\n'
  printf 'estimate_residual_variance\tFALSE\n'
  printf 'refine\tTRUE\n'
  printf 'R_finite\tTRUE\n'
  printf 'R_mismatch\teb\n'
  printf 'max_iter\t50\n'
  printf 'check_input\tTRUE\n'
  printf 'official_reference\thttps://stephenslab.github.io/susieR/reference/susie_rss.html\n'
  printf 'official_mismatch_guidance\thttps://stephenslab.github.io/susieR/articles/rss_mismatch.html\n'
} >"${RUN_DIR}/model_parameters.tsv"

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
  local branch=$5
  local started_at
  local finished_at
  local exit_code
  local completion_status
  local locus_output="${RUN_DIR}/${branch}/${analysis_scope}/${trait}/${locus_id}"
  local task_log="${RUN_DIR}/task_logs/${branch}/${analysis_scope}/${trait}/${locus_id}.log"

  started_at="$(date '+%F %T')"
  mkdir -p "$(dirname "${locus_output}")" "$(dirname "${task_log}")"

  if "${RSCRIPT_BIN}" --vanilla "${WORKER_R}" \
      --locus-id "${locus_id}" \
      --manifest "${LOCUS_MANIFEST}" \
      --qc-summary "${QC_SUMMARY}" \
      --metadata-manifest "${STEP18_LOCUS_MANIFEST}" \
      --branch "${branch}" \
      --output-dir "${locus_output}" \
      >"${task_log}" 2>&1; then
    exit_code=0
  else
    exit_code=$?
  fi
  finished_at="$(date '+%F %T')"

  completion_status=""
  if [[ -s "${locus_output}/COMPLETED.ok" ]]; then
    completion_status="$(awk -F'\t' 'NR==2 {print $1}' "${locus_output}/COMPLETED.ok")"
  fi

  if (( exit_code == 0 )) \
    && [[ "${completion_status}" == "PASS" || "${completion_status}" == "PASS_WITH_WARNINGS" ]] \
    && [[ -s "${locus_output}/susie_fit.rds" ]] \
    && [[ -s "${locus_output}/variant_pip.tsv.gz" ]] \
    && [[ -s "${locus_output}/credible_sets.tsv" ]] \
    && [[ -s "${locus_output}/credible_set_members.tsv.gz" ]] \
    && [[ -s "${locus_output}/rss_diagnostics.tsv" ]] \
    && [[ -s "${locus_output}/kriging_rss.tsv.gz" ]] \
    && [[ -s "${locus_output}/locus_summary.tsv" ]] \
    && grep -q $'\tSUCCESS\t' "${locus_output}/run_status.tsv"; then
    append_status "${SUCCESS_TSV}" "${RUN_ID}\t${locus_id}\t${trait}\t${analysis_scope}\t${branch}\tSUCCESS\t0\t${started_at}\t${finished_at}\t${locus_output}\t${task_log}\t${completion_status}"
    log_message "PROGRESS" "completed ${task_index}/${TOTAL_TASKS}: ${locus_id}/${branch}, status=${completion_status}"
  else
    append_status "${FAILED_TSV}" "${RUN_ID}\t${locus_id}\t${trait}\t${analysis_scope}\t${branch}\tFAILED\t${exit_code}\t${started_at}\t${finished_at}\t${locus_output}\t${task_log}\tWORKER_OR_OUTPUT_VALIDATION_FAILED"
    log_message "PROGRESS" "failed ${task_index}/${TOTAL_TASKS}: ${locus_id}/${branch}, see ${task_log}"
  fi
  return 0
}

active_jobs=0
selected_index=0
dispatched_tasks=0
skipped_tasks=0

while IFS=$'\t' read -r locus_id trait analysis_scope locus_status _remaining; do
  if ! is_selected_task "${locus_id}" "${analysis_scope}"; then
    continue
  fi
  for branch in "${BRANCH_ARRAY[@]}"; do
    selected_index=$((selected_index + 1))
    if [[ "${locus_status}" != "PASS" ]]; then
      current_time="$(date '+%F %T')"
      append_status "${SKIPPED_TSV}" "${RUN_ID}\t${locus_id}\t${trait}\t${analysis_scope}\t${branch}\tSKIPPED\t0\t${current_time}\t${current_time}\tNA\tNA\tSTEP19_LOCUS_STATUS_${locus_status}"
      skipped_tasks=$((skipped_tasks + 1))
      log_message "PROGRESS" "skip ${selected_index}/${TOTAL_TASKS}: ${locus_id}/${branch}, reason=${locus_status}"
      continue
    fi

    run_locus "${selected_index}" "${locus_id}" "${trait}" "${analysis_scope}" "${branch}" &
    dispatched_tasks=$((dispatched_tasks + 1))
    active_jobs=$((active_jobs + 1))
    if (( active_jobs >= MAX_PARALLEL )); then
      wait -n
      active_jobs=$((active_jobs - 1))
    fi
  done
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
  log_message "ERROR" "Terminal task count ${terminal_count} differs from the selected task count ${TOTAL_TASKS}"
  exit 1
fi

LOCUS_RUN_MANIFEST="${RUN_DIR}/locus_run_manifest.tsv"
printf 'locus_id\ttrait\tanalysis_scope\tbranch\tstatus\tmodel_status\tgwas_lead_variant_id\tgwas_lead_pos\tgwas_lead_p_lrt\tgwas_lead_p_lrt_tier\tgwas_lead_in_harmonized_variants\tvariant_count\treference_sample_count\tn_locus_median\tmin_observed_p_lrt\tconverged\tniter\tL\tcoverage\testimate_s_rss\tcredible_set_count\tmax_pip\ttop_variant_id\ttop_variant_pos\tartifact_flag\tR_sensitivity_flag\tR_reliability_flag\tkriging_outlier_count\toutput_dir\tvariant_pip\tcredible_sets\tcredible_set_members\trss_diagnostics\tkriging_diagnostics\tsusie_fit\n' >"${LOCUS_RUN_MANIFEST}"

while IFS=$'\t' read -r _run_id _locus_id _trait _scope _branch _status _exit_code _started _finished locus_output _task_log _reason; do
  tail -n +2 "${locus_output}/locus_summary.tsv" >>"${LOCUS_RUN_MANIFEST}"
done < <(tail -n +2 "${SUCCESS_TSV}")

total_variant_count="$(awk -F'\t' 'NR>1 {sum += $12} END {print sum+0}' "${LOCUS_RUN_MANIFEST}")"
total_credible_set_count="$(awk -F'\t' 'NR>1 {sum += $21} END {print sum+0}' "${LOCUS_RUN_MANIFEST}")"
warning_locus_count="$(awk -F'\t' 'NR>1 && $5 == "PASS_WITH_WARNINGS" {count++} END {print count+0}' "${LOCUS_RUN_MANIFEST}")"

BRANCH_SUMMARY="${RUN_DIR}/branch_summary.tsv"
printf 'branch\tselected_task_count\tsuccess_task_count\tfailed_task_count\tskipped_task_count\tpass_with_warnings_count\tmodeled_variant_count\tcredible_set_count\n' >"${BRANCH_SUMMARY}"
for branch in "${BRANCH_ARRAY[@]}"; do
  branch_success="$(awk -F'\t' -v branch="${branch}" 'NR>1 && $5 == branch {count++} END {print count+0}' "${SUCCESS_TSV}")"
  branch_failed="$(awk -F'\t' -v branch="${branch}" 'NR>1 && $5 == branch {count++} END {print count+0}' "${FAILED_TSV}")"
  branch_skipped="$(awk -F'\t' -v branch="${branch}" 'NR>1 && $5 == branch {count++} END {print count+0}' "${SKIPPED_TSV}")"
  branch_warning="$(awk -F'\t' -v branch="${branch}" 'NR>1 && $4 == branch && $5 == "PASS_WITH_WARNINGS" {count++} END {print count+0}' "${LOCUS_RUN_MANIFEST}")"
  branch_variants="$(awk -F'\t' -v branch="${branch}" 'NR>1 && $4 == branch {sum += $12} END {print sum+0}' "${LOCUS_RUN_MANIFEST}")"
  branch_credible_sets="$(awk -F'\t' -v branch="${branch}" 'NR>1 && $4 == branch {sum += $21} END {print sum+0}' "${LOCUS_RUN_MANIFEST}")"
  branch_terminal=$((branch_success + branch_failed + branch_skipped))
  if (( branch_terminal != SELECTED_LOCUS_COUNT )); then
    log_message "ERROR" "Branch ${branch}: terminal task count ${branch_terminal} differs from the selected locus count ${SELECTED_LOCUS_COUNT}"
    exit 1
  fi
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "${branch}" "${SELECTED_LOCUS_COUNT}" "${branch_success}" \
    "${branch_failed}" "${branch_skipped}" "${branch_warning}" \
    "${branch_variants}" "${branch_credible_sets}" >>"${BRANCH_SUMMARY}"
done

{
  printf 'metric\tvalue\n'
  printf 'input_locus_count\t%s\n' "${INPUT_LOCUS_COUNT}"
  printf 'selected_locus_count\t%s\n' "${SELECTED_LOCUS_COUNT}"
  printf 'selected_branch_count\t%s\n' "${BRANCH_COUNT}"
  printf 'selected_task_count\t%s\n' "${TOTAL_TASKS}"
  printf 'dispatched_task_count\t%s\n' "${dispatched_tasks}"
  printf 'success_task_count\t%s\n' "${success_count}"
  printf 'failed_task_count\t%s\n' "${failed_count}"
  printf 'skipped_task_count\t%s\n' "${skipped_count}"
  printf 'pass_with_warnings_task_count\t%s\n' "${warning_locus_count}"
  printf 'modeled_variant_count\t%s\n' "${total_variant_count}"
  printf 'credible_set_count\t%s\n' "${total_credible_set_count}"
} >"${RUN_DIR}/overall_summary.tsv"

if (( failed_count > 0 )); then
  {
    printf 'overall_status\tFAILED\n'
    printf 'run_id\t%s\n' "${RUN_ID}"
    printf 'branches\t%s\n' "${BRANCHES}"
    printf 'selected_locus_count\t%s\n' "${SELECTED_LOCUS_COUNT}"
    printf 'selected_task_count\t%s\n' "${TOTAL_TASKS}"
    printf 'success_task_count\t%s\n' "${success_count}"
    printf 'failed_task_count\t%s\n' "${failed_count}"
    printf 'skipped_task_count\t%s\n' "${skipped_count}"
    printf 'completed_at\t%s\n' "$(date '+%F %T')"
  } >"${RUN_DIR}/COMPLETED.failed"
  log_message "FAILED" "${failed_count} items locus-branch task failed; see ${FAILED_TSV}"
  exit 1
fi

overall_status="PASS"
if (( skipped_count > 0 || warning_locus_count > 0 )); then
  overall_status="PASS_WITH_WARNINGS"
fi

{
  printf 'overall_status\t%s\n' "${overall_status}"
  printf 'run_id\t%s\n' "${RUN_ID}"
  printf 'branches\t%s\n' "${BRANCHES}"
  printf 'selected_locus_count\t%s\n' "${SELECTED_LOCUS_COUNT}"
  printf 'selected_task_count\t%s\n' "${TOTAL_TASKS}"
  printf 'success_task_count\t%s\n' "${success_count}"
  printf 'failed_task_count\t%s\n' "${failed_count}"
  printf 'skipped_task_count\t%s\n' "${skipped_count}"
  printf 'completed_at\t%s\n' "$(date '+%F %T')"
} >"${RUN_DIR}/COMPLETED.ok"

log_message "DONE" "Step 1.10 task completed / done; status=${overall_status}, success task=${success_count}, skip task=${skipped_count}"
log_message "OUTPUT" "run directory: ${RUN_DIR}"
log_message "OUTPUT" "branch summary: ${BRANCH_SUMMARY}"
log_message "OUTPUT" "success list: ${SUCCESS_TSV}"
log_message "OUTPUT" "failure list: ${FAILED_TSV}"
log_message "OUTPUT" "skip list: ${SKIPPED_TSV}"
