#!/usr/bin/env bash


set -euo pipefail

SCRIPT_PATH="/path/to/repository/scripts/05_joint_genotyping_phasing/dog11/1.3_run_dog11_gvcf_calling/1.3_run_dog11_gvcf_calling.batch.sh"
MAIN_SCRIPT="/path/to/repository/scripts/05_joint_genotyping_phasing/dog11/1.3_run_dog11_gvcf_calling/1.3_run_dog11_gvcf_calling.sh"
TARGETS_TSV="/path/to/project/analysis_20260709/5-6dogs/1.ASE-processing/1.2_prepare_dog11_gvcf_calling/dog11_gvcf_input_targets.tsv"
OUT_DIR="/path/to/project/analysis_20260709/5-6dogs/1.ASE-processing/1.3_run_dog11_gvcf_calling"
BATCH_LOG="${OUT_DIR}/1.3_run_dog11_gvcf_calling.batch.log"

log_stage() {
  printf '[%(%F %T)T ] %s\n' -1 "$*" | tee -a "$BATCH_LOG"
}

fail() {
  log_stage "ERROR: $*"
  exit 1
}

mkdir -p "$OUT_DIR"
: > "$BATCH_LOG"

log_stage "start 1.3 dog11 gVCF calling batch wrapper."
log_stage "wrapper script path: ${SCRIPT_PATH}"
log_stage "main script path: ${MAIN_SCRIPT}"
log_stage "input targets: ${TARGETS_TSV}"
log_stage "output directory: ${OUT_DIR}"

if [[ $# -ne 0 ]]; then
  fail "this batch wrapper does not accept arguments; default direct run."
fi

if [[ ! -s "$MAIN_SCRIPT" ]]; then
  fail "main script does not exist or is empty: ${MAIN_SCRIPT}"
fi

if [[ ! -s "$TARGETS_TSV" ]]; then
  fail "Step 2 targets does not exist or is empty: ${TARGETS_TSV}"
fi

log_stage "invoke main script: bash ${MAIN_SCRIPT}"
bash "$MAIN_SCRIPT" 2>&1 | tee -a "$BATCH_LOG"

log_stage "main script return completed."
log_stage "main script run summary: ${OUT_DIR}/run_summary.tsv"
log_stage "main script batch status: ${OUT_DIR}/dog11_gvcf_batch_status.tsv"
log_stage "task completed / done"
