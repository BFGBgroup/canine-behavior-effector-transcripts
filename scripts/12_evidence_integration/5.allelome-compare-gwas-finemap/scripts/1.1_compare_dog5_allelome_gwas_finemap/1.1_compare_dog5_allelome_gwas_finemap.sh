#!/usr/bin/bash
set -euo pipefail


PYTHON_BIN="/usr/bin/python3.10"
SCRIPT_DIR="/path/to/repository/scripts/12_evidence_integration/5.allelome-compare-gwas-finemap/scripts/1.1_compare_dog5_allelome_gwas_finemap"
PYTHON_SCRIPT="/path/to/repository/scripts/12_evidence_integration/5.allelome-compare-gwas-finemap/scripts/1.1_compare_dog5_allelome_gwas_finemap/1.1.compare_dog5_allelome_gwas_finemap.py"
ALLELOME_LINKS="/path/to/project/analysis_20260817/5-dogs/7.Allelome.LINK/1.4_annotate_and_merge_dog5_allelome_links/20260907_213429/dog5_allelome_links.annotated.tsv"
ALLELOME_SUMMARY="/path/to/project/analysis_20260817/5-dogs/7.Allelome.LINK/1.4_annotate_and_merge_dog5_allelome_links/20260907_213429/dog5_allelome_annotation.summary.tsv"
NEW_GENE_GTF="/path/to/project/analysis_20260817/reference/Manual-annotation/1.1_build_five_only_new_gene_gtf/20260818_204722/five_only_new_gene_loci.new_only.gtf"
NEW_GENE_METADATA="/path/to/project/analysis_20260817/reference/Manual-annotation/1.1_build_five_only_new_gene_gtf/20260818_204722/five_only_new_gene_loci.transcript_map.tsv"
NEW_GENE_BUILD_SUMMARY="/path/to/project/analysis_20260817/reference/Manual-annotation/1.1_build_five_only_new_gene_gtf/20260818_204722/five_only_new_gene_loci.summary.tsv"
GWAS_RUN_DIR="/path/to/project/analysis_20260709/5-6dogs/2.gwas-finemap/work/1.10_run_dog_behavior_susie_rss_finemapping/20260802_205352"
OUTPUT_ROOT="/path/to/project/analysis_20260817/results/5-dogs/5.allelome-compare-gwas-finemap/1.1_compare_dog5_allelome_gwas_finemap"
AUDIT_DIR="${SCRIPT_DIR}/log"
MASTER_LOG="${AUDIT_DIR}/1.1_compare_dog5_allelome_gwas_finemap.log"
SUCCESS_FILE="${AUDIT_DIR}/1.1_compare_dog5_allelome_gwas_finemap.success.tsv"
FAILED_FILE="${AUDIT_DIR}/1.1_compare_dog5_allelome_gwas_finemap.failed.tsv"
LOCK_FILE="${AUDIT_DIR}/.1.1_compare_dog5_allelome_gwas_finemap.lock"

/usr/bin/mkdir -p "${AUDIT_DIR}"
exec 9>"${LOCK_FILE}"
if ! /usr/bin/flock -n 9; then
  echo "[$(/usr/bin/date '+%Y-%m-%d %H:%M:%S') ] existing the same launcher in run, refusing a concurrent launch" >&2
  exit 1
fi

RUN_ID="${RUN_ID:-$(/usr/bin/date '+%Y%m%d_%H%M%S')}"
[[ "${RUN_ID}" =~ ^[0-9]{8}_[0-9]{6}$ ]] || { echo "RUN_ID must be YYYYMMDD_HHMMSS" >&2; exit 2; }
RUN_DIR="${OUTPUT_ROOT}/${RUN_ID}"
if [[ -e "${RUN_DIR}" ]]; then
  echo "[$(/usr/bin/date '+%Y-%m-%d %H:%M:%S') ] run directory already exists, refuse to overwrite: ${RUN_DIR}" >&2
  exit 1
fi

: >"${MASTER_LOG}"
printf 'run_id\tstatus\texit_code\tstarted_at\tfinished_at\toutput_path\ttask_log\treason\n' >"${SUCCESS_FILE}"
printf 'run_id\tstatus\texit_code\tstarted_at\tfinished_at\toutput_path\ttask_log\treason\n' >"${FAILED_FILE}"
exec >>"${MASTER_LOG}" 2>&1

started_at="$(/usr/bin/date '+%Y-%m-%d %H:%M:%S')"
printf '[%s ] [start] run_id=%s\n' "${started_at}" "${RUN_ID}"
printf '[%s ] [input] allelome_links=%s\n' "${started_at}" "${ALLELOME_LINKS}"
printf '[%s ] [input] gwas_run_dir=%s\n' "${started_at}" "${GWAS_RUN_DIR}"
printf '[%s ] [output] run_dir=%s\n' "${started_at}" "${RUN_DIR}"
printf '[%s ] [audit] master_log=%s success=%s failed=%s\n' \
  "${started_at}" "${MASTER_LOG}" "${SUCCESS_FILE}" "${FAILED_FILE}"

export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
printf '[%s ] resources: concurrent_tasks=1 threads_per_task=1(serial settings retained)\n' "${started_at}"
child_pid=""
stop_worker() {
  local code="$1"
  trap '' TERM INT
  if [[ -n "${child_pid}" ]]; then
    kill -TERM "${child_pid}" 2>/dev/null || true
    wait "${child_pid}" 2>/dev/null || true
  fi
  printf '%s\tfailed\t%s\t%s\t%s\t%s\t%s\tsignal_termination\n' "${RUN_ID}" "${code}" "${started_at}" "$(/usr/bin/date '+%Y-%m-%d %H:%M:%S')" "${RUN_DIR}" "${MASTER_LOG}" >>"${FAILED_FILE}"
  exit "${code}"
}
trap 'stop_worker 143' TERM
trap 'stop_worker 130' INT
exit_code=0
"${PYTHON_BIN}" -B "${PYTHON_SCRIPT}" \
  --allelome-links "${ALLELOME_LINKS}" \
  --allelome-summary "${ALLELOME_SUMMARY}" \
  --new-gene-gtf "${NEW_GENE_GTF}" \
  --new-gene-metadata "${NEW_GENE_METADATA}" \
  --new-gene-build-summary "${NEW_GENE_BUILD_SUMMARY}" \
  --gwas-run-dir "${GWAS_RUN_DIR}" \
  --output-root "${OUTPUT_ROOT}" \
  --run-id "${RUN_ID}" &
child_pid=$!
if wait "${child_pid}"; then
  exit_code=0
else
  exit_code=$?
fi

finished_at="$(/usr/bin/date '+%Y-%m-%d %H:%M:%S')"
task_log="${RUN_DIR}/1.1.compare_dog5_allelome_gwas_finemap.log"
if [[ "${exit_code}" -eq 0 ]]; then
  printf '%s\tsuccess\t0\t%s\t%s\t%s\t%s\t%s\n' \
    "${RUN_ID}" "${started_at}" "${finished_at}" "${RUN_DIR}" "${task_log}" \
    "worker_exit_0" >>"${SUCCESS_FILE}"
  printf '[%s ] [audit] total=1 success=1 failed=0\n' "${finished_at}"
  printf 'execution_status=success\n' >"${RUN_DIR}/execution_completion.marker"
  printf '[%s ] task completed / done; execution_status=success: %s\n' "${finished_at}" "${RUN_DIR}"
  printf '[%s ] please separately check result completeness, format with biological interpretation; execution success does not equal QC passed.\n' "${finished_at}"
  exit 0
fi

reason="worker_exit_code=${exit_code}"
printf '%s\tfailed\t%s\t%s\t%s\t%s\t%s\t%s\n' \
  "${RUN_ID}" "${exit_code}" "${started_at}" "${finished_at}" "${RUN_DIR}" \
  "${task_log}" "${reason}" >>"${FAILED_FILE}"
printf '[%s ] [audit] total=1 success=0 failed=1\n' "${finished_at}"
printf '[%s ] [failed] %s\n' "${finished_at}" "${reason}"
exit 1
