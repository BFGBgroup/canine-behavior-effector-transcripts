#!/usr/bin/env bash


set -euo pipefail
shopt -s nullglob

SCRIPT_PATH="/path/to/repository/scripts/05_joint_genotyping_phasing/dog11/1.6_rerun_dog11_whatsHap/1.6_rerun_dog11_whatsHap.sh"
STEP2_TARGETS="/path/to/project/analysis_20260709/5-6dogs/1.ASE-processing/1.2_prepare_dog11_gvcf_calling/dog11_gvcf_input_targets.tsv"
STEP5_SUCCESS="/path/to/project/analysis_20260709/5-6dogs/1.ASE-processing/1.5_prepare_dog11_phase_input/dog11_phase_input_success.tsv"

REFERENCE_FASTA="/path/to/project/legacy/reference/GCF_000002285.3_CanFam3.1_genomic_chr.fna"
REFERENCE_FAI="/path/to/project/legacy/reference/GCF_000002285.3_CanFam3.1_genomic_chr.fna.fai"

OUT_DIR="/path/to/project/analysis_20260709/5-6dogs/1.ASE-processing/1.6_rerun_dog11_whatsHap"
PHASED_ROOT="${OUT_DIR}/phased_vcf"
LOG_ROOT="${OUT_DIR}/logs"
TMP_ROOT="${OUT_DIR}/tmp"
STATE_DIR="${OUT_DIR}/state"
STATUS_TSV="${OUT_DIR}/dog11_whatsHap_status.tsv"
SUCCESS_TASK_TSV="${OUT_DIR}/dog11_whatsHap_success.tsv"
FAILED_TASK_TSV="${OUT_DIR}/dog11_whatsHap_failed.tsv"
RUN_SUMMARY_TSV="${OUT_DIR}/run_summary.tsv"
LOG_FILE="${OUT_DIR}/1.6_rerun_dog11_whatsHap.log"
LOCK_FILE="${OUT_DIR}/1.6_rerun_dog11_whatsHap.lock"

ENV_BIN="/path/to/software/miniconda3/envs/dogs-harmony-latest/bin"
WHATSHAP_BIN="${ENV_BIN}/whatshap"
BCFTOOLS_BIN="${ENV_BIN}/bcftools"
SAMTOOLS_BIN="${ENV_BIN}/samtools"
export PATH="${ENV_BIN}:${PATH}"

EXPECTED_SAMPLE_COUNT=11
EXPECTED_CHROMOSOME_COUNT=38
EXPECTED_TASK_COUNT=418
PARALLEL_TASKS=60
BCFTOOLS_THREADS=1

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
  [[ -s "$path" ]] || fail "${label} does not exist or is empty: ${path}"
}

wait_for_slot() {
  while [[ "$(jobs -pr | wc -l)" -ge "$PARALLEL_TASKS" ]]; do
    sleep 20
  done
}

run_one_task() {
  local project="$1"
  local sample_id="$2"
  local chromosome="$3"
  local input_bam="$4"
  local input_vcf="$5"
  local task_id="${sample_id}.${chromosome}"
  local sample_output_dir="${PHASED_ROOT}/${sample_id}"
  local output_vcf="${sample_output_dir}/${sample_id}.${chromosome}.whatshap.phased.vcf.gz"
  local output_tbi="${output_vcf}.tbi"
  local task_log="${LOG_ROOT}/${task_id}.whatshap.log"
  local task_tmp="${TMP_ROOT}/${task_id}"
  local single_sample_vcf="${task_tmp}/${sample_id}.${chromosome}.input.vcf.gz"
  local temporary_output="${task_tmp}/${sample_id}.${chromosome}.whatshap.phased.vcf.gz"
  local task_state="${STATE_DIR}/${task_id}.status.tsv"
  local output_samples=0
  local output_records=0
  local phased_het_count=0
  local ps_count=0
  local status="pending"

  mkdir -p "$sample_output_dir" "$task_tmp"
  : > "$task_log"

  {
    printf '[%(%F %T)T ] start task: %s/%s/%s\n' -1 "$project" "$sample_id" "$chromosome"
    printf '[%(%F %T)T ] input phase-input VCF: %s\n' -1 "$input_vcf"
    printf '[%(%F %T)T ] input WGS BAM: %s\n' -1 "$input_bam"
    printf '[%(%F %T)T ] output WhatsHap VCF: %s\n' -1 "$output_vcf"
  } >> "$task_log"

  if [[ -s "$output_vcf" && -s "$output_tbi" ]]; then
    output_samples="$($BCFTOOLS_BIN query -l "$output_vcf" | wc -l)"
    output_records="$($BCFTOOLS_BIN view -H "$output_vcf" | wc -l)"
    read -r phased_het_count ps_count < <(
      "$BCFTOOLS_BIN" query -f '[%GT\t%PS\n]' "$output_vcf" \
        | awk '$1=="0|1" || $1=="1|0" {n++; if($2!="." && $2!="") ps[$2]=1} END {print n+0, length(ps)}'
    )
    status="already_done"
  elif [[ -e "$output_vcf" || -e "$output_tbi" ]]; then
    status="failed_incomplete_existing_output"
  elif [[ -e "$single_sample_vcf" || -e "${single_sample_vcf}.tbi" || -e "$temporary_output" || -e "${temporary_output}.tbi" ]]; then
    status="failed_stale_temporary_output"
  else
    printf '[%(%F %T)T ] stage 1/4: from 11 sample VCF extract target sample.\n' -1 >> "$task_log"
    if "$BCFTOOLS_BIN" view --threads "$BCFTOOLS_THREADS" -s "$sample_id" -Oz -o "$single_sample_vcf" "$input_vcf" >> "$task_log" 2>&1 \
      && "$BCFTOOLS_BIN" index --threads "$BCFTOOLS_THREADS" -t "$single_sample_vcf" >> "$task_log" 2>&1; then

      printf '[%(%F %T)T ] stage 2/4: use WGS reads run WhatsHap local phase.\n' -1 >> "$task_log"
      if "$WHATSHAP_BIN" phase \
          --reference "$REFERENCE_FASTA" \
          --sample "$sample_id" \
          --chromosome "$chromosome" \
          --tag PS \
          --output "$temporary_output" \
          "$single_sample_vcf" "$input_bam" >> "$task_log" 2>&1; then

        printf '[%(%F %T)T ] stage 3/4: create output index.\n' -1 >> "$task_log"
        if "$BCFTOOLS_BIN" index --threads "$BCFTOOLS_THREADS" -t "$temporary_output" >> "$task_log" 2>&1; then
          printf '[%(%F %T)T ] stage 4/4: check sample, record count, phased het and PS count.\n' -1 >> "$task_log"
          output_samples="$($BCFTOOLS_BIN query -l "$temporary_output" | wc -l)"
          output_records="$($BCFTOOLS_BIN view -H "$temporary_output" | wc -l)"
          read -r phased_het_count ps_count < <(
            "$BCFTOOLS_BIN" query -f '[%GT\t%PS\n]' "$temporary_output" \
              | awk '$1=="0|1" || $1=="1|0" {n++; if($2!="." && $2!="") ps[$2]=1} END {print n+0, length(ps)}'
          )

          if [[ "$output_samples" -ne 1 ]]; then
            status="failed_output_sample_count"
          elif [[ "$output_records" -le 0 ]]; then
            status="failed_empty_output"
          elif [[ "$phased_het_count" -le 0 || "$ps_count" -le 0 ]]; then
            status="failed_no_phased_heterozygote"
          else
            mv "$temporary_output" "$output_vcf"
            mv "${temporary_output}.tbi" "$output_tbi"
            rm -f "$single_sample_vcf" "${single_sample_vcf}.tbi"
            status="run_success"
          fi
        else
          status="failed_output_index_exit"
        fi
      else
        status="failed_whatshap_exit"
      fi
    else
      status="failed_sample_extract_exit"
    fi
  fi

  {
    printf '[%(%F %T)T ] task status: %s; output_samples=%s; output_records=%s; phased_het=%s; ps_count=%s\n' \
      -1 "$status" "$output_samples" "$output_records" "$phased_het_count" "$ps_count"
    printf '[%(%F %T)T ] single task completed / done\n' -1
  } >> "$task_log"

  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$project" "$sample_id" "$chromosome" "$input_bam" "$input_vcf" "$output_vcf" \
    "$output_records" "$phased_het_count" "$ps_count" "$status" > "$task_state"
}

if [[ $# -ne 0 ]]; then
  printf 'ERROR: this script does not accept arguments.\n' >&2
  exit 1
fi

mkdir -p "$OUT_DIR" "$PHASED_ROOT" "$LOG_ROOT" "$TMP_ROOT" "$STATE_DIR"
: > "$LOG_FILE"

exec 9> "$LOCK_FILE"
if ! /usr/bin/flock -n 9; then
  fail "existing another items 1.6 task holds run lock: ${LOCK_FILE}"
fi

log_stage "start rerunning dog11 WhatsHap local phase."
log_stage "script path: ${SCRIPT_PATH}"
log_stage "input BAM list: ${STEP2_TARGETS}"
log_stage "input phase-input success table: ${STEP5_SUCCESS}"
log_stage "output directory: ${OUT_DIR}"
log_stage "parallel task count=${PARALLEL_TASKS}; bcftools threads=${BCFTOOLS_THREADS}; expected task count=${EXPECTED_TASK_COUNT}"

require_file "$STEP2_TARGETS" "Step 2 BAM list"
require_file "$STEP5_SUCCESS" "Step 5 phase-input success table"
require_file "$REFERENCE_FASTA" "CanFam3.1 FASTA"
require_file "$REFERENCE_FAI" "CanFam3.1 FASTA index"
[[ -x "$WHATSHAP_BIN" ]] || fail "WhatsHap does not exist or is not executable: ${WHATSHAP_BIN}"
[[ -x "$BCFTOOLS_BIN" ]] || fail "bcftools does not exist or is not executable: ${BCFTOOLS_BIN}"
[[ -x "$SAMTOOLS_BIN" ]] || fail "samtools does not exist or is not executable: ${SAMTOOLS_BIN}"

expected_targets_header=$'project\tsample_id\treference_fasta\treference_fai\treference_dict\tinput_bam\tinput_bai\tready_for_gvcf\tpreflight_status\tstatus'
read -r targets_header < "$STEP2_TARGETS"
[[ "$targets_header" == "$expected_targets_header" ]] || fail "Step 2 BAM list unexpected header."

expected_phase_header=$'chromosome\tphase_input_vcf\toutput_records\tstatus'
read -r phase_header < "$STEP5_SUCCESS"
[[ "$phase_header" == "$expected_phase_header" ]] || fail "Step 5 success table has an unexpected header."

rm -f "${STATE_DIR}"/*.status.tsv

declare -A PROJECT_BY_SAMPLE=()
declare -A BAM_BY_SAMPLE=()
declare -A BAI_BY_SAMPLE=()
SAMPLE_COUNT=0
while IFS=$'\t' read -r project sample_id reference_fasta reference_fai reference_dict input_bam input_bai ready_for_gvcf preflight_status target_status; do
  [[ "$ready_for_gvcf" == "yes" && "$preflight_status" == "pass" && "$target_status" == "ready" ]] || \
    fail "BAM the list contains samples that are not ready: ${sample_id}"
  [[ -z "${BAM_BY_SAMPLE[$sample_id]+x}" ]] || fail "BAM list contains duplicate sample: ${sample_id}"
  require_file "$input_bam" "${sample_id} WGS BAM"
  require_file "$input_bai" "${sample_id} WGS BAM index"

  mapfile -t bam_sm < <("$SAMTOOLS_BIN" view -H "$input_bam" | awk -F '\t' '$1=="@RG" {for(i=1;i<=NF;i++) if($i~/^SM:/) print substr($i,4)}' | sort -u)
  [[ "${#bam_sm[@]}" -eq 1 && "${bam_sm[0]}" == "$sample_id" ]] || \
    fail "BAM @RG/SM does not match the sample ID: ${sample_id}; BAM=${input_bam}"

  PROJECT_BY_SAMPLE["$sample_id"]="$project"
  BAM_BY_SAMPLE["$sample_id"]="$input_bam"
  BAI_BY_SAMPLE["$sample_id"]="$input_bai"
  SAMPLE_COUNT=$((SAMPLE_COUNT + 1))
done < <(tail -n +2 "$STEP2_TARGETS")
[[ "$SAMPLE_COUNT" -eq "$EXPECTED_SAMPLE_COUNT" ]] || fail "unexpected sample count: expected=${EXPECTED_SAMPLE_COUNT}; observed=${SAMPLE_COUNT}"

declare -A VCF_BY_CHROM=()
CHROMOSOME_COUNT=0
while IFS=$'\t' read -r chromosome phase_input_vcf output_records status; do
  [[ "$status" == "run_success" || "$status" == "already_done" ]] || fail "Step 5 success table contains a non-success status: ${chromosome}/${status}"
  [[ "$chromosome" =~ ^chr([1-9]|[12][0-9]|3[0-8])$ ]] || fail "chromosome name does not match chr1-chr38: ${chromosome}"
  [[ -z "${VCF_BY_CHROM[$chromosome]+x}" ]] || fail "Step 5 success table contains duplicate chromosomes: ${chromosome}"
  require_file "$phase_input_vcf" "${chromosome} phase-input VCF"
  require_file "${phase_input_vcf}.tbi" "${chromosome} phase-input VCF index"
  input_sample_count="$($BCFTOOLS_BIN query -l "$phase_input_vcf" | wc -l)"
  [[ "$input_sample_count" -eq "$EXPECTED_SAMPLE_COUNT" ]] || fail "${chromosome} phase-input unexpected sample count: ${input_sample_count}"
  VCF_BY_CHROM["$chromosome"]="$phase_input_vcf"
  CHROMOSOME_COUNT=$((CHROMOSOME_COUNT + 1))
done < <(tail -n +2 "$STEP5_SUCCESS")
[[ "$CHROMOSOME_COUNT" -eq "$EXPECTED_CHROMOSOME_COUNT" ]] || fail "unexpected chromosome count: expected=${EXPECTED_CHROMOSOME_COUNT}; observed=${CHROMOSOME_COUNT}"

mapfile -t CHROMOSOMES < <(awk -F '\t' '$1 ~ /^chr([1-9]|[12][0-9]|3[0-8])$/ {print $1}' "$REFERENCE_FAI")
mapfile -t SAMPLES < <(printf '%s\n' "${!BAM_BY_SAMPLE[@]}" | sort)

SUBMITTED_COUNT=0
for sample_id in "${SAMPLES[@]}"; do
  for chromosome in "${CHROMOSOMES[@]}"; do
    wait_for_slot
    SUBMITTED_COUNT=$((SUBMITTED_COUNT + 1))
    if (( SUBMITTED_COUNT == 1 || SUBMITTED_COUNT % 20 == 0 || SUBMITTED_COUNT == EXPECTED_TASK_COUNT )); then
      log_stage "submit task ${SUBMITTED_COUNT}/${EXPECTED_TASK_COUNT}: ${sample_id}/${chromosome}"
    fi
    run_one_task "${PROJECT_BY_SAMPLE[$sample_id]}" "$sample_id" "$chromosome" \
      "${BAM_BY_SAMPLE[$sample_id]}" "${VCF_BY_CHROM[$chromosome]}" &
  done
done

[[ "$SUBMITTED_COUNT" -eq "$EXPECTED_TASK_COUNT" ]] || fail "unexpected number of submitted tasks: expected=${EXPECTED_TASK_COUNT}; observed=${SUBMITTED_COUNT}"
log_stage "submitted all task, waiting for background tasks to finish."
wait
log_stage "all task already return, starting the summary status."

printf 'project\tsample_id\tchromosome\tinput_bam\tinput_phase_vcf\toutput_phased_vcf\toutput_records\tphased_het_count\tps_count\tstatus\n' > "$STATUS_TSV"
printf 'sample_id\tchromosome\toutput_phased_vcf\tphased_het_count\tps_count\tstatus\n' > "$SUCCESS_TASK_TSV"
printf 'sample_id\tchromosome\toutput_phased_vcf\tstatus\n' > "$FAILED_TASK_TSV"

STATUS_ROW_COUNT=0
RUN_SUCCESS_COUNT=0
ALREADY_DONE_COUNT=0
FAILED_COUNT=0
TOTAL_PHASED_HET=0
TOTAL_PS=0

for sample_id in "${SAMPLES[@]}"; do
  for chromosome in "${CHROMOSOMES[@]}"; do
    task_state="${STATE_DIR}/${sample_id}.${chromosome}.status.tsv"
    if [[ ! -s "$task_state" ]]; then
      FAILED_COUNT=$((FAILED_COUNT + 1))
      printf '%s\t%s\t%s\t%s\n' "$sample_id" "$chromosome" "${PHASED_ROOT}/${sample_id}/${sample_id}.${chromosome}.whatshap.phased.vcf.gz" "failed_missing_state" >> "$FAILED_TASK_TSV"
      continue
    fi

    IFS=$'\t' read -r project state_sample state_chromosome input_bam input_vcf output_vcf output_records phased_het_count ps_count status < "$task_state"
    STATUS_ROW_COUNT=$((STATUS_ROW_COUNT + 1))
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      "$project" "$state_sample" "$state_chromosome" "$input_bam" "$input_vcf" "$output_vcf" \
      "$output_records" "$phased_het_count" "$ps_count" "$status" >> "$STATUS_TSV"

    case "$status" in
      run_success)
        RUN_SUCCESS_COUNT=$((RUN_SUCCESS_COUNT + 1))
        TOTAL_PHASED_HET=$((TOTAL_PHASED_HET + phased_het_count))
        TOTAL_PS=$((TOTAL_PS + ps_count))
        printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$state_sample" "$state_chromosome" "$output_vcf" "$phased_het_count" "$ps_count" "$status" >> "$SUCCESS_TASK_TSV"
        ;;
      already_done)
        ALREADY_DONE_COUNT=$((ALREADY_DONE_COUNT + 1))
        TOTAL_PHASED_HET=$((TOTAL_PHASED_HET + phased_het_count))
        TOTAL_PS=$((TOTAL_PS + ps_count))
        printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$state_sample" "$state_chromosome" "$output_vcf" "$phased_het_count" "$ps_count" "$status" >> "$SUCCESS_TASK_TSV"
        ;;
      *)
        FAILED_COUNT=$((FAILED_COUNT + 1))
        printf '%s\t%s\t%s\t%s\n' "$state_sample" "$state_chromosome" "$output_vcf" "$status" >> "$FAILED_TASK_TSV"
        ;;
    esac
  done
done

{
  printf 'key\tvalue\n'
  printf 'step\t%s\n' "1.6_rerun_dog11_whatsHap"
  printf 'script_path\t%s\n' "$SCRIPT_PATH"
  printf 'input_bam_targets\t%s\n' "$STEP2_TARGETS"
  printf 'input_phase_success\t%s\n' "$STEP5_SUCCESS"
  printf 'reference_fasta\t%s\n' "$REFERENCE_FASTA"
  printf 'output_dir\t%s\n' "$OUT_DIR"
  printf 'phased_vcf_root\t%s\n' "$PHASED_ROOT"
  printf 'sample_count\t%s\n' "$SAMPLE_COUNT"
  printf 'chromosome_count\t%s\n' "$CHROMOSOME_COUNT"
  printf 'expected_task_count\t%s\n' "$EXPECTED_TASK_COUNT"
  printf 'status_row_count\t%s\n' "$STATUS_ROW_COUNT"
  printf 'run_success_count\t%s\n' "$RUN_SUCCESS_COUNT"
  printf 'already_done_count\t%s\n' "$ALREADY_DONE_COUNT"
  printf 'failed_count\t%s\n' "$FAILED_COUNT"
  printf 'total_phased_het_count\t%s\n' "$TOTAL_PHASED_HET"
  printf 'total_ps_count\t%s\n' "$TOTAL_PS"
  printf 'parallel_tasks\t%s\n' "$PARALLEL_TASKS"
  printf 'bcftools_threads\t%s\n' "$BCFTOOLS_THREADS"
  printf 'whatshap_bin\t%s\n' "$WHATSHAP_BIN"
  printf 'whatshap_version\t%s\n' "$($WHATSHAP_BIN --version)"
  printf 'lock_file\t%s\n' "$LOCK_FILE"
} > "$RUN_SUMMARY_TSV"

log_stage "write the status table: ${STATUS_TSV}"
log_stage "write success table: ${SUCCESS_TASK_TSV}"
log_stage "write the failure table: ${FAILED_TASK_TSV}"
log_stage "total tasks=${EXPECTED_TASK_COUNT}; status_rows=${STATUS_ROW_COUNT}; run_success=${RUN_SUCCESS_COUNT}; already_done=${ALREADY_DONE_COUNT}; failed=${FAILED_COUNT}; phased_het=${TOTAL_PHASED_HET}; PS=${TOTAL_PS}"

if [[ "$FAILED_COUNT" -gt 0 || "$STATUS_ROW_COUNT" -ne "$EXPECTED_TASK_COUNT" ]]; then
  log_stage "Some tasks failed or are missing; see ${FAILED_TASK_TSV}"
  exit 1
fi

log_stage "task completed / done"
