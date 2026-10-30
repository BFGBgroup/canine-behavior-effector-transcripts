#!/usr/bin/env bash
set -euo pipefail


PYTHON_BIN="/usr/bin/python3"
SCRIPT_DIR="/path/to/repository/scripts/12_evidence_integration/8.new-gene-allelome-eqtl-gwas-finemap-integration/summary/scripts"
PYTHON_SCRIPT="${SCRIPT_DIR}/1.1_build_dog5_newgene_allelome_eqtl_gwas_finemap_summary.py"

SOURCE_RUN_DIR="/path/to/project/analysis_20260817/results/5-dogs/8.new-gene-allelome-eqtl-gwas-finemap-integration/1.1_integrate_dog5_newgene_allelome_eqtl_gwas_finemap/20260910_195332"
CANDIDATES="${SOURCE_RUN_DIR}/dog5_newgene_eqtl_gwas_finemap_all_candidates.tsv"
PCGENE_LINKS="${SOURCE_RUN_DIR}/dog5_integrated_newgene_pcGene_links.tsv.gz"
CANDIDATE_PAIRS="${SOURCE_RUN_DIR}/dog5_candidate_newgene_pcGene_mechanism_pairs.tsv"
SOURCE_SUMMARY="${SOURCE_RUN_DIR}/dog5_newgene_allelome_eqtl_gwas_finemap_summary.tsv"
SOURCE_MANIFEST="${SOURCE_RUN_DIR}/dog5_newgene_allelome_eqtl_gwas_finemap_manifest.tsv"
FIVE_ONLY_NEW_GENE_MAP="/path/to/project/analysis_20260817/reference/Manual-annotation/1.1_build_five_only_new_gene_gtf/20260818_204722/five_only_new_gene_loci.transcript_map.tsv"

OUTPUT_ROOT="/path/to/project/analysis_20260817/results/5-dogs/8.new-gene-allelome-eqtl-gwas-finemap-integration/summary/results"
RUN_ID="${RUN_ID:-$(/usr/bin/date '+%Y%m%d_%H%M%S')}"
RUN_DIR="${OUTPUT_ROOT}/${RUN_ID}"

LOG_DIR="${SCRIPT_DIR}/log"
MASTER_LOG="${LOG_DIR}/1.1_build_dog5_newgene_allelome_eqtl_gwas_finemap_summary.log"
SUCCESS_FILE="${LOG_DIR}/1.1_build_dog5_newgene_allelome_eqtl_gwas_finemap_summary.success.tsv"
FAILED_FILE="${LOG_DIR}/1.1_build_dog5_newgene_allelome_eqtl_gwas_finemap_summary.failed.tsv"
LOCK_FILE="${LOG_DIR}/1.1_build_dog5_newgene_allelome_eqtl_gwas_finemap_summary.lock"

/usr/bin/mkdir -p "${LOG_DIR}"
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

log_message() {
  printf '[%s ] %s\n' "$(/usr/bin/date '+%Y-%m-%d %H:%M:%S')" "$1"
}

STARTED_AT="$(/usr/bin/date '+%Y-%m-%d %H:%M:%S')"
log_message "task started; total_tasks=1; run_id=${RUN_ID}"
log_message "Step 8 input directory: ${SOURCE_RUN_DIR}"
log_message "five-dog-specific new gene map: ${FIVE_ONLY_NEW_GENE_MAP}"
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
  --candidates "${CANDIDATES}" \
  --pcgene-links "${PCGENE_LINKS}" \
  --candidate-pairs "${CANDIDATE_PAIRS}" \
  --source-summary "${SOURCE_SUMMARY}" \
  --source-manifest "${SOURCE_MANIFEST}" \
  --five-only-new-gene-map "${FIVE_ONLY_NEW_GENE_MAP}" \
  --output-root "${OUTPUT_ROOT}" \
  --run-id "${RUN_ID}"
EXIT_CODE=$?
set -e

FINISHED_AT="$(/usr/bin/date '+%Y-%m-%d %H:%M:%S')"
REASON="python_exit_code_${EXIT_CODE}"
if (( EXIT_CODE == 0 )); then
  if printf 'execution_status=success\nrun_id=%s\ncompleted_at=%s\n' \
    "${RUN_ID}" "${FINISHED_AT}" > "${RUN_DIR}/execution_completion.marker"; then
    REASON="python_exit_code_0_and_normal_end"
    printf '%s\tsuccess\t0\t%s\t%s\t%s\t%s\t%s\n' \
      "${RUN_ID}" "${STARTED_AT}" "${FINISHED_AT}" "${RUN_DIR}" \
      "${MASTER_LOG}" "${REASON}" >> "${SUCCESS_FILE}"
    log_message "task summary: total=1; success=1; failed=0"
    log_message "task completed / done; output=${RUN_DIR}"
    log_message "execution_status=success"
    log_message "completion only indicates that the main program exited normally; output completeness, format, counts and biological quality require separate review"
    exit 0
  fi
  EXIT_CODE=74
  REASON="execution_completion_marker_write_failed"
fi

printf '%s\tfailed\t%s\t%s\t%s\t%s\t%s\t%s\n' \
  "${RUN_ID}" "${EXIT_CODE}" "${STARTED_AT}" "${FINISHED_AT}" \
  "${RUN_DIR}" "${MASTER_LOG}" "${REASON}" >> "${FAILED_FILE}"
log_message "task summary: total=1; success=0; failed=1; reason=${REASON}"
log_message "failure list: ${FAILED_FILE}"
log_message "execution_status=failed"
exit "${EXIT_CODE}"
