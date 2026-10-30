#!/usr/bin/env bash
set -euo pipefail


PYTHON_BIN="/usr/bin/python3"
SCRIPT_DIR="/path/to/repository/scripts/12_evidence_integration/8.new-gene-allelome-eqtl-gwas-finemap-integration/scripts/1.1_integrate_dog5_newgene_allelome_eqtl_gwas_finemap"
PYTHON_SCRIPT="/path/to/repository/scripts/12_evidence_integration/8.new-gene-allelome-eqtl-gwas-finemap-integration/scripts/1.1_integrate_dog5_newgene_allelome_eqtl_gwas_finemap/1.1.integrate_dog5_newgene_allelome_eqtl_gwas_finemap.py"

STEP6_LINKS="/path/to/project/analysis_20260817/results/5-dogs/6.allelome-compare-eqtl-finemap/1.1_compare_dog5_allelome_eqtl_finemap/20260910_183740/dog5_allelome_eqtl_finemap_expression_gene_link_overlap.tsv.gz"
STEP6_SUMMARY="/path/to/project/analysis_20260817/results/5-dogs/6.allelome-compare-eqtl-finemap/1.1_compare_dog5_allelome_eqtl_finemap/20260910_183740/dog5_allelome_eqtl_finemap_run_summary.tsv"

STEP3_RUN_DIR="/path/to/project/analysis_20260817/results/5-dogs/3.eqtl-finemap-compare-gwas-finemap/1.1_compare_dog5_eqtl_gwas_finemap/20260910_131929"
STEP3_MAIN_DIR="${STEP3_RUN_DIR}/main"
STEP3_BOTH_CREDIBLE="${STEP3_MAIN_DIR}/dog5_expression_ase_ppa_ge_0.1_gwas_both_credible_sets_exact_overlap.tsv"
STEP3_NONCREDIBLE="${STEP3_MAIN_DIR}/dog5_expression_ase_ppa_ge_0.1_gwas_noncredible_locus_variant_overlap.tsv"
STEP3_SUMMARY="${STEP3_MAIN_DIR}/dog5_expression_ase_gwas_finemap_overlap_summary.tsv"
GWAS_RUN_DIR="/path/to/project/analysis_20260709/5-6dogs/2.gwas-finemap/work/1.10_run_dog_behavior_susie_rss_finemapping/20260802_205352"

OUTPUT_ROOT="/path/to/project/analysis_20260817/results/5-dogs/8.new-gene-allelome-eqtl-gwas-finemap-integration/1.1_integrate_dog5_newgene_allelome_eqtl_gwas_finemap"
RUN_ID="${RUN_ID:-$(/usr/bin/date '+%Y%m%d_%H%M%S')}"
RUN_DIR="${OUTPUT_ROOT}/${RUN_ID}"
TASK_LOG="${RUN_DIR}/1.1_integrate_dog5_newgene_allelome_eqtl_gwas_finemap.log"

AUDIT_DIR="${SCRIPT_DIR}/log"
MASTER_LOG="${AUDIT_DIR}/1.1_integrate_dog5_newgene_allelome_eqtl_gwas_finemap.log"
SUCCESS_FILE="${AUDIT_DIR}/1.1_integrate_dog5_newgene_allelome_eqtl_gwas_finemap.success.tsv"
FAILED_FILE="${AUDIT_DIR}/1.1_integrate_dog5_newgene_allelome_eqtl_gwas_finemap.failed.tsv"
LOCK_FILE="${AUDIT_DIR}/1.1_integrate_dog5_newgene_allelome_eqtl_gwas_finemap.lock"

/usr/bin/mkdir -p "${AUDIT_DIR}"
exec 9>"${LOCK_FILE}"
if ! /usr/bin/flock -n 9; then
  printf '[%s ] another task for this step is already running, no new task was started: %s\n' \
    "$(/usr/bin/date '+%Y-%m-%d %H:%M:%S')" "${LOCK_FILE}" >&2
  exit 75
fi

: > "${MASTER_LOG}"
printf 'run_id\tstatus\texit_code\tstarted_at\tfinished_at\toutput_path\ttask_log\treason\n' > "${SUCCESS_FILE}"
printf 'run_id\tstatus\texit_code\tstarted_at\tfinished_at\toutput_path\ttask_log\treason\n' > "${FAILED_FILE}"
exec >> "${MASTER_LOG}" 2>&1

STARTED_AT="$(/usr/bin/date '+%Y-%m-%d %H:%M:%S')"

log_message() {
  printf '[%s ] %s\n' "$(/usr/bin/date '+%Y-%m-%d %H:%M:%S')" "$1"
}

log_message "task started; total_tasks=1; run_id=${RUN_ID}"
log_message "Step 3 main input directory: ${STEP3_MAIN_DIR}"
log_message "Step 6 relationship input: ${STEP6_LINKS}"
log_message "GWAS fine-map input directory: ${GWAS_RUN_DIR}"
log_message "planned output directory: ${RUN_DIR}"
log_message "resource configuration: concurrent_tasks=1; threads_per_task=1; preserved_serial_script"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1

set +e
"${PYTHON_BIN}" "${PYTHON_SCRIPT}" \
  --step6-links "${STEP6_LINKS}" \
  --step6-summary "${STEP6_SUMMARY}" \
  --step3-both-credible "${STEP3_BOTH_CREDIBLE}" \
  --step3-noncredible "${STEP3_NONCREDIBLE}" \
  --step3-summary "${STEP3_SUMMARY}" \
  --gwas-run-dir "${GWAS_RUN_DIR}" \
  --output-root "${OUTPUT_ROOT}" \
  --run-id "${RUN_ID}"
EXIT_CODE=$?
set -e

STATUS="success"
REASON="python_exit_code_0_and_normal_end"
if (( EXIT_CODE != 0 )); then
  STATUS="failed"
  REASON="python_exit_code_${EXIT_CODE}"
fi

FINISHED_AT="$(/usr/bin/date '+%Y-%m-%d %H:%M:%S')"
if [[ "${STATUS}" == "success" ]]; then
  printf 'execution_status=success\nrun_id=%s\ncompleted_at=%s\n' \
    "${RUN_ID}" "${FINISHED_AT}" > "${RUN_DIR}/execution_completion.marker"
  printf '%s\tsuccess\t0\t%s\t%s\t%s\t%s\t%s\n' \
    "${RUN_ID}" "${STARTED_AT}" "${FINISHED_AT}" "${RUN_DIR}" \
    "${TASK_LOG}" "${REASON}" >> "${SUCCESS_FILE}"
  log_message "task summary: total=1; success=1; failed=0"
  log_message "task completed / done; output=${RUN_DIR}"
  log_message "execution_status=success"
  log_message "completion only indicates that the main program exited normally; output completeness, format, counts and biological quality require separate review"
  exit 0
fi

printf '%s\tfailed\t%s\t%s\t%s\t%s\t%s\t%s\n' \
  "${RUN_ID}" "${EXIT_CODE}" "${STARTED_AT}" "${FINISHED_AT}" \
  "${RUN_DIR}" "${TASK_LOG}" "${REASON}" >> "${FAILED_FILE}"
log_message "task summary: total=1; success=0; failed=1; reason=${REASON}"
log_message "failure list: ${FAILED_FILE}"
log_message "execution_status=failed"
exit 1
