#!/usr/bin/env bash
set -Eeuo pipefail


SCRIPT_DIR="/path/to/repository/scripts/11_gwas_finemapping/1.2_build_pnas_dog_behavior_locus_manifest"
PYTHON_BIN="/path/to/software/miniconda3/envs/dogs-harmony-latest/bin/python"
WORKER="${SCRIPT_DIR}/1.2_build_pnas_dog_behavior_locus_manifest.py"

STEP11_RUN_ROOT="/path/to/project/analysis_20260709/5-6dogs/2.gwas-finemap/work/1.1_validate_dog_behavior_gwas_finemap_inputs"
STEP11_RUN_ID="${STEP11_RUN_ID:-20260801_182359}"
STEP11_RUN_DIR="${STEP11_RUN_ROOT}/${STEP11_RUN_ID}"
GWAS_DIR="/path/to/project/legacy/reference/gwas/summary_statistics"
ARTICLE_XLSX="/path/to/project/analysis_20260709/5-6dogs/2.gwas-finemap/article/pnas.2421757122.sd01 (1).xlsx"
ARTICLE_PDF="/path/to/project/analysis_20260709/5-6dogs/2.gwas-finemap/article/GWAS for behavioral traits in golden retrievers identifiesgenes implicated in human temperament, mental health,and cognition.pdf"

OUTPUT_ROOT="/path/to/project/analysis_20260709/5-6dogs/2.gwas-finemap/work/1.2_build_pnas_dog_behavior_locus_manifest"
RUN_ID="${RUN_ID:-$(date '+%Y%m%d_%H%M%S')}"
RUN_DIR="${OUTPUT_ROOT}/${RUN_ID}"

AUDIT_DIR="${SCRIPT_DIR}/log"
LOCK_FILE="${AUDIT_DIR}/1.2_build_pnas_dog_behavior_locus_manifest.lock"
MASTER_LOG="${AUDIT_DIR}/1.2_build_pnas_dog_behavior_locus_manifest.${RUN_ID}.log"
SUCCESS_FILE="${AUDIT_DIR}/1.2_build_pnas_dog_behavior_locus_manifest.${RUN_ID}.success.tsv"
FAILED_FILE="${AUDIT_DIR}/1.2_build_pnas_dog_behavior_locus_manifest.${RUN_ID}.failed.tsv"

timestamp() {
  date '+[%Y-%m-%d %H:%M:%S ]'
}

log() {
  printf '%s %s\n' "$(timestamp)" "$*"
}

mkdir -p "${AUDIT_DIR}"
exec 9>"${LOCK_FILE}"
if ! flock -n 9; then
  printf '%s Another step 1.2 task is already running.\n' "$(timestamp)" >&2
  exit 1
fi

if [[ ! "${RUN_ID}" =~ ^[0-9]{8}_[0-9]{6}$ ]]; then
  printf '%s Invalid RUN_ID; expected YYYYMMDD_HHMMSS: %s\n' \
    "$(timestamp)" "${RUN_ID}" >&2
  exit 1
fi

if [[ -e "${RUN_DIR}" || -e "${MASTER_LOG}" || -e "${SUCCESS_FILE}" || -e "${FAILED_FILE}" ]]; then
  printf '%s Refusing to overwrite an existing run: %s\n' \
    "$(timestamp)" "${RUN_ID}" >&2
  exit 1
fi

: > "${MASTER_LOG}"
printf 'task_id\trun_id\tstatus\texit_code\tstarted_at\tfinished_at\toutput_path\ttask_log\treason\n' > "${SUCCESS_FILE}"
printf 'task_id\trun_id\tstatus\texit_code\tstarted_at\tfinished_at\toutput_path\ttask_log\treason\n' > "${FAILED_FILE}"
exec >> "${MASTER_LOG}" 2>&1

STARTED_AT="$(date '+%Y-%m-%d %H:%M:%S')"
TASK_ID="build_dog_behavior_gwas_summary_candidate_manifest"

log "Task started"
log "Run ID: ${RUN_ID}"
log "Validated step 1.1 run: ${STEP11_RUN_DIR}"
log "Primary GWAS summary directory: ${GWAS_DIR}"
log "Optional PNAS Dataset S1 audit input: ${ARTICLE_XLSX}"
log "Optional PNAS Table 1 audit source: ${ARTICLE_PDF}"
log "Output directory: ${RUN_DIR}"
log "Master log: ${MASTER_LOG}"

set +e
"${PYTHON_BIN}" "${WORKER}" \
  --validation-dir "${STEP11_RUN_DIR}" \
  --gwas-dir "${GWAS_DIR}" \
  --article-xlsx "${ARTICLE_XLSX}" \
  --article-pdf "${ARTICLE_PDF}" \
  --output-dir "${RUN_DIR}"
WORKER_EXIT=$?
set -e

FINISHED_AT="$(date '+%Y-%m-%d %H:%M:%S')"
SUMMARY_FILE="${RUN_DIR}/manifest_summary.tsv"
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

if [[ "${WORKER_EXIT}" -eq 0 && "${SUMMARY_OK}" -eq 1 && -s "${COMPLETION_MARKER}" ]]; then
  printf '%s\t%s\tsuccess\t0\t%s\t%s\t%s\t%s\t%s\n' \
    "${TASK_ID}" "${RUN_ID}" "${STARTED_AT}" "${FINISHED_AT}" \
    "${RUN_DIR}" "${MASTER_LOG}" "summary_candidate_manifest_built" >> "${SUCCESS_FILE}"
  log "Successful tasks: 1"
  log "Failed tasks: 0"
  log "Success table: ${SUCCESS_FILE}"
  log "Failure table: ${FAILED_FILE}"
  log "Task completed / done"
  exit 0
fi

REASON="worker_exit_or_output_validation_failed"
printf '%s\t%s\tfailed\t%s\t%s\t%s\t%s\t%s\t%s\n' \
  "${TASK_ID}" "${RUN_ID}" "${WORKER_EXIT}" "${STARTED_AT}" "${FINISHED_AT}" \
  "${RUN_DIR}" "${MASTER_LOG}" "${REASON}" >> "${FAILED_FILE}"
log "Successful tasks: 0"
log "Failed tasks: 1"
log "Failure table: ${FAILED_FILE}"
log "Task failed; no valid completion marker was accepted"
exit 1
