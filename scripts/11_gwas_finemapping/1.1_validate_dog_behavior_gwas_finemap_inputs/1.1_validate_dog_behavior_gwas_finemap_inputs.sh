#!/usr/bin/env bash
set -Eeuo pipefail


SCRIPT_DIR="/path/to/repository/scripts/11_gwas_finemapping/1.1_validate_dog_behavior_gwas_finemap_inputs"
PYTHON_BIN="/path/to/software/miniconda3/envs/dogs-harmony-latest/bin/python"
WORKER="${SCRIPT_DIR}/1.1_validate_dog_behavior_gwas_finemap_inputs.py"

GWAS_ROOT="/path/to/project/legacy/reference/gwas"
GWAS_DIR="${GWAS_ROOT}/summary_statistics"
GWAS_README="${GWAS_ROOT}/README.md"
ARTICLE_XLSX="/path/to/project/analysis_20260709/5-6dogs/2.gwas-finemap/article/pnas.2421757122.sd01 (1).xlsx"
LD_PREFIX="${GWAS_ROOT}/LD/GRLS_Genotyping_Axiom_Canine_Merged_2023_08_02/AxiomGT1v2.merged"
LD_README="${GWAS_ROOT}/LD/README.txt"
PLINK_BIN="/path/to/software/miniconda3/envs/dogs-harmony-latest/bin/plink"
RSCRIPT_BIN="/path/to/software/environments/plasma_py27/bin/Rscript"

OUTPUT_ROOT="/path/to/project/analysis_20260709/5-6dogs/2.gwas-finemap/work/1.1_validate_dog_behavior_gwas_finemap_inputs"
RUN_ID="${RUN_ID:-$(date '+%Y%m%d_%H%M%S')}"
RUN_DIR="${OUTPUT_ROOT}/${RUN_ID}"

AUDIT_DIR="${SCRIPT_DIR}/log"
LOCK_FILE="${AUDIT_DIR}/1.1_validate_dog_behavior_gwas_finemap_inputs.lock"
MASTER_LOG="${AUDIT_DIR}/1.1_validate_dog_behavior_gwas_finemap_inputs.${RUN_ID}.log"
SUCCESS_FILE="${AUDIT_DIR}/1.1_validate_dog_behavior_gwas_finemap_inputs.${RUN_ID}.success.tsv"
FAILED_FILE="${AUDIT_DIR}/1.1_validate_dog_behavior_gwas_finemap_inputs.${RUN_ID}.failed.tsv"

timestamp() {
  date '+[%Y-%m-%d %H:%M:%S ]'
}

log() {
  printf '%s %s\n' "$(timestamp)" "$*"
}

mkdir -p "${AUDIT_DIR}"
exec 9>"${LOCK_FILE}"
if ! flock -n 9; then
  printf '%s another task for this step is already running, no new task was started.\n' "$(timestamp)" >&2
  exit 1
fi

if [[ ! "${RUN_ID}" =~ ^[0-9]{8}_[0-9]{6}$ ]]; then
  printf '%s RUN_ID format error, expected YYYYMMDD_HHMMSS: %s\n' "$(timestamp)" "${RUN_ID}" >&2
  exit 1
fi

if [[ -e "${RUN_DIR}" || -e "${MASTER_LOG}" || -e "${SUCCESS_FILE}" || -e "${FAILED_FILE}" ]]; then
  printf '%s RUN_ID corresponding output or audit file already exists, refuse to overwrite: %s\n' \
    "$(timestamp)" "${RUN_ID}" >&2
  exit 1
fi

: > "${MASTER_LOG}"
printf 'task_id\trun_id\tstatus\texit_code\tstarted_at\tfinished_at\toutput_path\ttask_log\treason\n' > "${SUCCESS_FILE}"
printf 'task_id\trun_id\tstatus\texit_code\tstarted_at\tfinished_at\toutput_path\ttask_log\treason\n' > "${FAILED_FILE}"
exec >> "${MASTER_LOG}" 2>&1

STARTED_AT="$(date '+%Y-%m-%d %H:%M:%S')"
TASK_ID="validate_dog_behavior_gwas_finemap_inputs"

log "task started"
log "run identifier: ${RUN_ID}"
log "GWAS input directory: ${GWAS_DIR}"
log "PNAS supplement table: ${ARTICLE_XLSX}"
log "LD reference prefix: ${LD_PREFIX}"
log "output directory: ${RUN_DIR}"
log "main log: ${MASTER_LOG}"

set +e
"${PYTHON_BIN}" "${WORKER}" \
  --gwas-dir "${GWAS_DIR}" \
  --gwas-readme "${GWAS_README}" \
  --article-xlsx "${ARTICLE_XLSX}" \
  --ld-prefix "${LD_PREFIX}" \
  --ld-readme "${LD_README}" \
  --plink-bin "${PLINK_BIN}" \
  --rscript-bin "${RSCRIPT_BIN}" \
  --output-dir "${RUN_DIR}"
WORKER_EXIT=$?
set -e

FINISHED_AT="$(date '+%Y-%m-%d %H:%M:%S')"
SUMMARY_FILE="${RUN_DIR}/validation_summary.tsv"
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
    "${RUN_DIR}" "${MASTER_LOG}" "validation_passed" >> "${SUCCESS_FILE}"
  log "success task count: 1"
  log "failed task count: 0"
  log "success list: ${SUCCESS_FILE}"
  log "failure list: ${FAILED_FILE}"
  log "task completed / done"
  exit 0
fi

REASON="worker_exit_or_output_validation_failed"
printf '%s\t%s\tfailed\t%s\t%s\t%s\t%s\t%s\t%s\n' \
  "${TASK_ID}" "${RUN_ID}" "${WORKER_EXIT}" "${STARTED_AT}" "${FINISHED_AT}" \
  "${RUN_DIR}" "${MASTER_LOG}" "${REASON}" >> "${FAILED_FILE}"
log "success task count: 0"
log "failed task count: 1"
log "failure list: ${FAILED_FILE}"
log "task failed, not written done marker"
exit 1
