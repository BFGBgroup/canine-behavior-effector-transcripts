#!/usr/bin/env bash


set -euo pipefail
shopt -s nullglob

SCRIPT_PATH="/path/to/repository/scripts/05_joint_genotyping_phasing/dog11/1.3_run_dog11_gvcf_calling/1.3_run_dog11_gvcf_calling.sh"
STEP2_OUT_DIR="/path/to/project/analysis_20260709/5-6dogs/1.ASE-processing/1.2_prepare_dog11_gvcf_calling"
TARGETS_TSV="${STEP2_OUT_DIR}/dog11_gvcf_input_targets.tsv"

OUT_DIR="/path/to/project/analysis_20260709/5-6dogs/1.ASE-processing/1.3_run_dog11_gvcf_calling"
GVCF_ROOT="${OUT_DIR}/gvcf"
LOG_ROOT="${OUT_DIR}/logs"
TMP_ROOT="${OUT_DIR}/tmp"
STATE_DIR="${OUT_DIR}/state"
STATUS_TSV="${OUT_DIR}/dog11_gvcf_batch_status.tsv"
SUCCESS_TSV="${OUT_DIR}/dog11_gvcf_batch_success.tsv"
FAILED_TSV="${OUT_DIR}/dog11_gvcf_batch_failed.tsv"
RUN_SUMMARY_TSV="${OUT_DIR}/run_summary.tsv"
LOG_FILE="${OUT_DIR}/1.3_run_dog11_gvcf_calling.log"
LOCK_FILE="${OUT_DIR}/1.3_run_dog11_gvcf_calling.lock"
GATK_ENV_BIN="/path/to/software/miniconda3/envs/dogs-harmony-latest/bin"
GATK_BIN="${GATK_ENV_BIN}/gatk"
JAVA_BIN="${GATK_ENV_BIN}/java"


export PATH="${GATK_ENV_BIN}:${PATH}"

PARALLEL_JOBS=12
THREADS=5
JAVA_MEM="48g"

log_stage() {
  printf '[%(%F %T)T ] %s\n' -1 "$*" | tee -a "$LOG_FILE"
}

fail() {
  log_stage "ERROR: $*"
  exit 1
}

require_file() {
  local path="$1"
  local label="$2"

  if [[ ! -s "$path" ]]; then
    fail "${label} does not exist or is empty: ${path}"
  fi
}

file_exists_flag() {
  local path="$1"

  if [[ -s "$path" ]]; then
    printf 'yes'
  else
    printf 'no'
  fi
}

wait_for_slot() {
  while [[ "$(jobs -pr | wc -l)" -ge "$PARALLEL_JOBS" ]]; do
    sleep 30
  done
}

run_one_sample() {
  local project="$1"
  local sample_id="$2"
  local reference_fasta="$3"
  local reference_fai="$4"
  local reference_dict="$5"
  local input_bam="$6"
  local input_bai="$7"
  local ready_for_gvcf="$8"
  local preflight_status="$9"
  local target_status="${10}"

  local output_gvcf="${GVCF_ROOT}/${sample_id}/${sample_id}.g.vcf.gz"
  local output_gvcf_tbi="${output_gvcf}.tbi"
  local sample_log="${LOG_ROOT}/${sample_id}.gvcf.log"
  local sample_tmp_dir="${TMP_ROOT}/${sample_id}"
  local sample_status="${STATE_DIR}/${sample_id}.status.tsv"
  local output_gvcf_exists
  local output_gvcf_tbi_exists
  local batch_status="pending"

  mkdir -p "$(dirname "$output_gvcf")" "$LOG_ROOT" "$sample_tmp_dir" "$STATE_DIR"
  touch "$sample_log"

  {
    printf '[%(%F %T)T ] start sample: %s/%s\n' -1 "$project" "$sample_id"
    printf '[%(%F %T)T ] input BAM: %s\n' -1 "$input_bam"
    printf '[%(%F %T)T ] output GVCF: %s\n' -1 "$output_gvcf"
    printf '[%(%F %T)T ] preflight check status: %s; target status: %s\n' -1 "$preflight_status" "$target_status"
  } >> "$sample_log"

  output_gvcf_exists="$(file_exists_flag "$output_gvcf")"
  output_gvcf_tbi_exists="$(file_exists_flag "$output_gvcf_tbi")"

  if [[ "$target_status" != "ready" || "$ready_for_gvcf" != "yes" || "$preflight_status" != "pass" ]]; then
    batch_status="skipped_not_ready"
  elif [[ ! -s "$reference_fasta" || ! -s "$reference_fai" || ! -s "$reference_dict" || ! -s "$input_bam" || ! -s "$input_bai" ]]; then
    batch_status="failed_missing_required_input"
  elif [[ "$output_gvcf_exists" == "yes" && "$output_gvcf_tbi_exists" == "yes" ]]; then
    batch_status="already_done"
  elif [[ "$output_gvcf_exists" == "yes" || "$output_gvcf_tbi_exists" == "yes" ]]; then
    batch_status="failed_incomplete_existing_output"
    printf '[%(%F %T)T ] ERROR: incomplete existing output detected, please check manually: %s / %s\n' \
      -1 "$output_gvcf" "$output_gvcf_tbi" >> "$sample_log"
  else
    if "$GATK_BIN" --java-options "-Xmx${JAVA_MEM} -Djava.io.tmpdir=${sample_tmp_dir}" HaplotypeCaller \
      -R "$reference_fasta" \
      -I "$input_bam" \
      -O "$output_gvcf" \
      -ERC GVCF \
      --create-output-variant-index true \
      --native-pair-hmm-threads "$THREADS" >> "$sample_log" 2>&1; then
      output_gvcf_exists="$(file_exists_flag "$output_gvcf")"
      output_gvcf_tbi_exists="$(file_exists_flag "$output_gvcf_tbi")"
      if [[ "$output_gvcf_exists" == "yes" && "$output_gvcf_tbi_exists" == "yes" ]]; then
        batch_status="run_success"
      else
        batch_status="failed_missing_output_after_gatk"
      fi
    else
      batch_status="failed_gatk_exit"
    fi
  fi

  {
    printf '[%(%F %T)T ] sample status: %s\n' -1 "$batch_status"
    printf '[%(%F %T)T ] single sample task completed / done\n' -1
  } >> "$sample_log"

  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$project" "$sample_id" "$input_bam" "$output_gvcf" "$output_gvcf_tbi" \
    "$preflight_status" "$target_status" "$output_gvcf_exists" \
    "$output_gvcf_tbi_exists" "$batch_status" > "$sample_status"
}

mkdir -p "$OUT_DIR" "$GVCF_ROOT" "$LOG_ROOT" "$TMP_ROOT" "$STATE_DIR"
: > "$LOG_FILE"

exec 9> "$LOCK_FILE"
if ! /usr/bin/flock -n 9; then
  fail "existing another items 1.3 task holds run lock: ${LOCK_FILE}"
fi

rm -f "${STATE_DIR}"/*.status.tsv

log_stage "start dog11 gVCF calling batch execution."
log_stage "script path: ${SCRIPT_PATH}"
log_stage "input targets: ${TARGETS_TSV}"
log_stage "output directory: ${OUT_DIR}"
log_stage "gVCF output root directory: ${GVCF_ROOT}"
log_stage "GATK path: ${GATK_BIN}"
log_stage "Java path: ${JAVA_BIN}"
log_stage "parallel sample count: ${PARALLEL_JOBS}"
log_stage "single sample thread count: ${THREADS}; single sample Java memory: ${JAVA_MEM}"

if [[ $# -ne 0 ]]; then
  fail "this script does not accept arguments; default direct run, parallel count fixed is ${PARALLEL_JOBS}."
fi

require_file "$TARGETS_TSV" "Step 2 targets"

if [[ ! -x "$GATK_BIN" ]]; then
  fail "GATK does not exist or is not executable: ${GATK_BIN}"
fi
if [[ ! -x "$JAVA_BIN" ]]; then
  fail "Java does not exist or is not executable: ${JAVA_BIN}"
fi

expected_targets_header=$'project\tsample_id\treference_fasta\treference_fai\treference_dict\tinput_bam\tinput_bai\tready_for_gvcf\tpreflight_status\tstatus'
read -r targets_header < "$TARGETS_TSV"

if [[ "$targets_header" != "$expected_targets_header" ]]; then
  fail "targets unexpected header."
fi

TOTAL_SAMPLE_COUNT=0

while IFS=$'\t' read -r project sample_id reference_fasta reference_fai reference_dict input_bam input_bai ready_for_gvcf preflight_status target_status; do
  TOTAL_SAMPLE_COUNT=$((TOTAL_SAMPLE_COUNT + 1))
  wait_for_slot
  log_stage "submit sample ${TOTAL_SAMPLE_COUNT}: ${project}/${sample_id}"
  run_one_sample "$project" "$sample_id" "$reference_fasta" "$reference_fai" "$reference_dict" \
    "$input_bam" "$input_bai" "$ready_for_gvcf" "$preflight_status" "$target_status" &
done < <(tail -n +2 "$TARGETS_TSV")

log_stage "submitted all sample, waiting for background tasks to finish."
wait
log_stage "all background tasks have returned, starting the summary status."

printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
  "project" "sample_id" "input_bam" "output_gvcf" "output_gvcf_tbi" \
  "preflight_status" "target_status" "output_gvcf_exists" \
  "output_gvcf_tbi_exists" "batch_status" > "$STATUS_TSV"

printf '%s\t%s\t%s\t%s\n' "project" "sample_id" "output_gvcf" "batch_status" > "$SUCCESS_TSV"
printf '%s\t%s\t%s\t%s\n' "project" "sample_id" "input_bam" "batch_status" > "$FAILED_TSV"

ALREADY_DONE_COUNT=0
RUN_SUCCESS_COUNT=0
FAILED_COUNT=0
SKIPPED_COUNT=0
STATUS_ROW_COUNT=0

for sample_status in "${STATE_DIR}"/*.status.tsv; do
  IFS=$'\t' read -r project sample_id input_bam output_gvcf output_gvcf_tbi preflight_status target_status output_gvcf_exists output_gvcf_tbi_exists batch_status < "$sample_status"
  STATUS_ROW_COUNT=$((STATUS_ROW_COUNT + 1))
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$project" "$sample_id" "$input_bam" "$output_gvcf" "$output_gvcf_tbi" \
    "$preflight_status" "$target_status" "$output_gvcf_exists" \
    "$output_gvcf_tbi_exists" "$batch_status" >> "$STATUS_TSV"

  case "$batch_status" in
    already_done)
      ALREADY_DONE_COUNT=$((ALREADY_DONE_COUNT + 1))
      printf '%s\t%s\t%s\t%s\n' "$project" "$sample_id" "$output_gvcf" "$batch_status" >> "$SUCCESS_TSV"
      ;;
    run_success)
      RUN_SUCCESS_COUNT=$((RUN_SUCCESS_COUNT + 1))
      printf '%s\t%s\t%s\t%s\n' "$project" "$sample_id" "$output_gvcf" "$batch_status" >> "$SUCCESS_TSV"
      ;;
    skipped_not_ready)
      SKIPPED_COUNT=$((SKIPPED_COUNT + 1))
      printf '%s\t%s\t%s\t%s\n' "$project" "$sample_id" "$input_bam" "$batch_status" >> "$FAILED_TSV"
      ;;
    *)
      FAILED_COUNT=$((FAILED_COUNT + 1))
      printf '%s\t%s\t%s\t%s\n' "$project" "$sample_id" "$input_bam" "$batch_status" >> "$FAILED_TSV"
      ;;
  esac
done

if [[ "$STATUS_ROW_COUNT" -ne "$TOTAL_SAMPLE_COUNT" ]]; then
  FAILED_COUNT=$((FAILED_COUNT + 1))
  log_stage "The status-record count (${STATUS_ROW_COUNT}) differs from the submitted-sample count (${TOTAL_SAMPLE_COUNT})."
fi

{
  printf 'key\tvalue\n'
  printf 'step\t%s\n' "1.3_run_dog11_gvcf_calling"
  printf 'script_path\t%s\n' "$SCRIPT_PATH"
  printf 'input_targets_tsv\t%s\n' "$TARGETS_TSV"
  printf 'output_dir\t%s\n' "$OUT_DIR"
  printf 'gvcf_root\t%s\n' "$GVCF_ROOT"
  printf 'log_root\t%s\n' "$LOG_ROOT"
  printf 'tmp_root\t%s\n' "$TMP_ROOT"
  printf 'state_dir\t%s\n' "$STATE_DIR"
  printf 'status_tsv\t%s\n' "$STATUS_TSV"
  printf 'success_tsv\t%s\n' "$SUCCESS_TSV"
  printf 'failed_tsv\t%s\n' "$FAILED_TSV"
  printf 'total_sample_count\t%s\n' "$TOTAL_SAMPLE_COUNT"
  printf 'status_row_count\t%s\n' "$STATUS_ROW_COUNT"
  printf 'already_done_count\t%s\n' "$ALREADY_DONE_COUNT"
  printf 'run_success_count\t%s\n' "$RUN_SUCCESS_COUNT"
  printf 'skipped_count\t%s\n' "$SKIPPED_COUNT"
  printf 'failed_count\t%s\n' "$FAILED_COUNT"
  printf 'parallel_jobs\t%s\n' "$PARALLEL_JOBS"
  printf 'threads_per_sample\t%s\n' "$THREADS"
  printf 'java_mem_per_sample\t%s\n' "$JAVA_MEM"
  printf 'gatk_bin\t%s\n' "$GATK_BIN"
  printf 'java_bin\t%s\n' "$JAVA_BIN"
  printf 'lock_file\t%s\n' "$LOCK_FILE"
} > "$RUN_SUMMARY_TSV"

log_stage "write status: ${STATUS_TSV}"
log_stage "write success: ${SUCCESS_TSV}"
log_stage "write failed: ${FAILED_TSV}"
log_stage "total samples=${TOTAL_SAMPLE_COUNT}; status_rows=${STATUS_ROW_COUNT}; already_done=${ALREADY_DONE_COUNT}; run_success=${RUN_SUCCESS_COUNT}; skipped=${SKIPPED_COUNT}; failed=${FAILED_COUNT}"

if [[ "$FAILED_COUNT" -gt 0 ]]; then
  log_stage "Some samples failed; see ${FAILED_TSV}"
  exit 1
fi

log_stage "task completed / done"
