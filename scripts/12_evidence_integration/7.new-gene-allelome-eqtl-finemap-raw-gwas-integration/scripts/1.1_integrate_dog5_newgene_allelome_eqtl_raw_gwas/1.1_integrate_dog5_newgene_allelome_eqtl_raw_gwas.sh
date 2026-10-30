#!/usr/bin/bash
set -euo pipefail


SCRIPT_DIR="/path/to/repository/scripts/12_evidence_integration/7.new-gene-allelome-eqtl-finemap-raw-gwas-integration/scripts/1.1_integrate_dog5_newgene_allelome_eqtl_raw_gwas"
LOG_DIR="${SCRIPT_DIR}/log"
LOCK_FILE="${LOG_DIR}/1.1_integrate_dog5_newgene_allelome_eqtl_raw_gwas.lock"
PYTHON_BIN="/usr/bin/python3.10"
PYTHON_SCRIPT="${SCRIPT_DIR}/1.1.integrate_dog5_newgene_allelome_eqtl_raw_gwas.py"

EQTL_GWAS_RUN_DIR="/path/to/project/analysis_20260817/results/5-dogs/2.eqtl-finemap-compare-raw-gwas/1.1_export_dog5_eqtl_finemap_gwas_exact_overlap/20260917_195542"
EQTL_GWAS_DETAIL="${EQTL_GWAS_RUN_DIR}/five_dogs_eqtl_finemap_gwas_exact_overlaps.csv"
EQTL_GWAS_SUMMARY="${EQTL_GWAS_RUN_DIR}/five_dogs_eqtl_finemap_gwas_exact_overlap_overall_summary.csv"
ALLELOME_EQTL_RUN_DIR="/path/to/project/analysis_20260817/results/5-dogs/6.allelome-compare-eqtl-finemap/1.1_compare_dog5_allelome_eqtl_finemap/20260910_183740"
ALLELOME_EQTL_LINKS="${ALLELOME_EQTL_RUN_DIR}/dog5_allelome_eqtl_finemap_expression_gene_link_overlap.tsv.gz"
ALLELOME_EQTL_SUMMARY="${ALLELOME_EQTL_RUN_DIR}/dog5_allelome_eqtl_finemap_run_summary.tsv"
NEW_GENE_RUN_DIR="/path/to/project/analysis_20260817/reference/Manual-annotation/1.1_build_five_only_new_gene_gtf/20260818_204722"
NEW_GENE_TRANSCRIPT_MAP="${NEW_GENE_RUN_DIR}/five_only_new_gene_loci.transcript_map.tsv"
NEW_GENE_BUILD_SUMMARY="${NEW_GENE_RUN_DIR}/five_only_new_gene_loci.summary.tsv"

OUTPUT_ROOT="/path/to/project/analysis_20260817/results/5-dogs/7.new-gene-allelome-eqtl-finemap-raw-gwas-integration/1.1_integrate_dog5_newgene_allelome_eqtl_raw_gwas"
RUN_ID="${RUN_ID:-$(/usr/bin/date '+%Y%m%d_%H%M%S')}"
RUN_DIR="${OUTPUT_ROOT}/${RUN_ID}"
COMPLETION_MARKER="${RUN_DIR}/execution_completion.marker"
SUCCESS_TSV="${LOG_DIR}/1.1_integrate_dog5_newgene_allelome_eqtl_raw_gwas.${RUN_ID}.success.tsv"
FAILED_TSV="${LOG_DIR}/1.1_integrate_dog5_newgene_allelome_eqtl_raw_gwas.${RUN_ID}.failed.tsv"
CONCURRENT_TASKS=1
THREADS_PER_TASK=1

log() {
    printf '[%s ] %s\n' "$(/usr/bin/date '+%Y-%m-%d %H:%M:%S')" "$*"
}

if [[ ! "${RUN_ID}" =~ ^[0-9]{8}_[0-9]{6}$ ]]; then
    printf '[%s ] execution_status=failed; RUN_ID must use YYYYMMDD_HHMMSS format: %s\n' \
        "$(/usr/bin/date '+%Y-%m-%d %H:%M:%S')" "${RUN_ID}" >&2
    exit 1
fi

/usr/bin/mkdir -p "${LOG_DIR}"
exec 9>"${LOCK_FILE}"
if ! /usr/bin/flock -n 9; then
    printf '[%s ] execution_status=failed; existing the same step task holds run lock: %s\n' \
        "$(/usr/bin/date '+%Y-%m-%d %H:%M:%S')" "${LOCK_FILE}" >&2
    exit 75
fi


MASTER_LOG="${LOG_DIR}/1.1_integrate_dog5_newgene_allelome_eqtl_raw_gwas.${RUN_ID}.log"
if [[ -e "${MASTER_LOG}" || -e "${SUCCESS_TSV}" || -e "${FAILED_TSV}" || -e "${RUN_DIR}" ]]; then
    log "[error] batch already exists, refuse to overwrite: ${RUN_ID}"
    exit 1
fi
: >"${MASTER_LOG}"
exec >>"${MASTER_LOG}" 2>&1

printf 'run_id\texecution_status\texit_code\toutput_path\n' >"${SUCCESS_TSV}"
printf 'run_id\texecution_status\texit_code\toutput_path\n' >"${FAILED_TSV}"
NORMAL_END=0
MARKER_OWNED=0
CHILD_PID=""
finish() {
    local code=$?
    trap - EXIT
    if (( code == 0 && NORMAL_END == 1 )); then
        if printf '%s\tsuccess\t0\t%s\n' "${RUN_ID}" "${RUN_DIR}" >>"${SUCCESS_TSV}"; then
            log "task completed / done; execution_status=success"
            log "completion only indicates that the main program exited normally; output completeness, format, counts and biological quality require separate review"
            exit 0
        fi
        code=1
    fi
    (( code != 0 )) || code=1

    if (( MARKER_OWNED == 1 )); then
        /usr/bin/rm -f -- "${COMPLETION_MARKER}"
    fi
    printf '%s\tfailed\t%s\t%s\n' "${RUN_ID}" "${code}" "${RUN_DIR}" >>"${FAILED_TSV}"
    log "execution_status=failed; exit_code=${code}"
    exit "${code}"
}
stop_child() {
    trap '' INT TERM HUP
    if [[ -n "${CHILD_PID}" ]]; then
        kill -TERM "${CHILD_PID}" 2>/dev/null || true
        wait "${CHILD_PID}" 2>/dev/null || true
    fi
    exit 143
}
trap finish EXIT
trap stop_child INT TERM HUP

log "[stage 1/3] start five-dog novel gene multiple level evidence integrate"
log "python_script=${PYTHON_SCRIPT}"
log "eqtl_gwas_run_dir=${EQTL_GWAS_RUN_DIR}"
log "allelome_eqtl_run_dir=${ALLELOME_EQTL_RUN_DIR}"
log "new_gene_run_dir=${NEW_GENE_RUN_DIR}"
log "output=${RUN_DIR}"
log "concurrent_tasks=${CONCURRENT_TASKS}; threads_per_task=${THREADS_PER_TASK}; resource_baseline=preserved_serial_script"

export OMP_NUM_THREADS="${THREADS_PER_TASK}"
export OPENBLAS_NUM_THREADS="${THREADS_PER_TASK}"
export MKL_NUM_THREADS="${THREADS_PER_TASK}"
export NUMEXPR_NUM_THREADS="${THREADS_PER_TASK}"
export VECLIB_MAXIMUM_THREADS="${THREADS_PER_TASK}"

log "[stage 2/3] invoke Python main program"
set +e
"${PYTHON_BIN}" -B -u "${PYTHON_SCRIPT}" \
    --eqtl-gwas-detail "${EQTL_GWAS_DETAIL}" \
    --eqtl-gwas-summary "${EQTL_GWAS_SUMMARY}" \
    --allelome-eqtl-links "${ALLELOME_EQTL_LINKS}" \
    --allelome-eqtl-summary "${ALLELOME_EQTL_SUMMARY}" \
    --new-gene-transcript-map "${NEW_GENE_TRANSCRIPT_MAP}" \
    --new-gene-build-summary "${NEW_GENE_BUILD_SUMMARY}" \
    --output-root "${OUTPUT_ROOT}" \
    --run-id "${RUN_ID}" &
CHILD_PID=$!
wait "${CHILD_PID}"
exit_code=$?
CHILD_PID=""
set -e

if (( exit_code != 0 )); then
    log "execution_status=failed; Python main program exit code=${exit_code}"
    exit "${exit_code}"
fi

MARKER_OWNED=1
printf 'execution_status=success\nrun_id=%s\ncompleted_at=%s\n' \
    "${RUN_ID}" \
    "$(/usr/bin/date '+%Y-%m-%d %H:%M:%S')" \
    >"${COMPLETION_MARKER}"

log "[stage 3/3] run entry point completed normally"
log "result directory: ${RUN_DIR}"
NORMAL_END=1
