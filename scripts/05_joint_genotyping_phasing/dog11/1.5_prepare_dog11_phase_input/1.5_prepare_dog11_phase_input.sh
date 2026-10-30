#!/usr/bin/env bash


set -euo pipefail
shopt -s nullglob

SCRIPT_PATH="/path/to/repository/scripts/05_joint_genotyping_phasing/dog11/1.5_prepare_dog11_phase_input/1.5_prepare_dog11_phase_input.sh"
STEP4_OUT_DIR="/path/to/project/analysis_20260709/5-6dogs/1.ASE-processing/1.4_joint_genotype_dog11"
SUCCESS_TSV="${STEP4_OUT_DIR}/dog11_joint_genotype_success.tsv"

REFERENCE_FASTA="/path/to/project/legacy/reference/GCF_000002285.3_CanFam3.1_genomic_chr.fna"
REFERENCE_FAI="/path/to/project/legacy/reference/GCF_000002285.3_CanFam3.1_genomic_chr.fna.fai"

OUT_DIR="/path/to/project/analysis_20260709/5-6dogs/1.ASE-processing/1.5_prepare_dog11_phase_input"
HARD_FILTERED_ROOT="${OUT_DIR}/hard_filtered_vcf"
PHASE_INPUT_ROOT="${OUT_DIR}/phase_input_vcf"
STATS_ROOT="${OUT_DIR}/stats"
LOG_ROOT="${OUT_DIR}/logs"
TMP_ROOT="${OUT_DIR}/tmp"
STATE_DIR="${OUT_DIR}/state"
STATUS_TSV="${OUT_DIR}/dog11_phase_input_status.tsv"
SUCCESS_CHROM_TSV="${OUT_DIR}/dog11_phase_input_success.tsv"
FAILED_CHROM_TSV="${OUT_DIR}/dog11_phase_input_failed.tsv"
RUN_SUMMARY_TSV="${OUT_DIR}/run_summary.tsv"
LOG_FILE="${OUT_DIR}/1.5_prepare_dog11_phase_input.log"
LOCK_FILE="${OUT_DIR}/1.5_prepare_dog11_phase_input.lock"

BCFTOOLS_BIN="/path/to/software/miniconda3/envs/dogs-harmony-latest/bin/bcftools"

EXPECTED_SAMPLE_COUNT=11
EXPECTED_CHROMOSOME_COUNT=38
PARALLEL_CHROMOSOMES=40
THREADS=2
MIN_DP=5
MIN_GQ=20
MAX_MISSING_RATE="0.10"


SITE_HARD_FILTER='QUAL<30 || INFO/QD<2.0 || INFO/FS>60.0 || INFO/SOR>3.0 || INFO/MQ<40.0 || INFO/MQRankSum<-12.5 || INFO/ReadPosRankSum<-8.0'

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
  while [[ "$(jobs -pr | wc -l)" -ge "$PARALLEL_CHROMOSOMES" ]]; do
    sleep 20
  done
}

run_one_chromosome() {
  local chromosome="$1"
  local input_vcf="$2"
  local input_tbi="${input_vcf}.tbi"
  local hard_filtered_vcf="${HARD_FILTERED_ROOT}/${chromosome}.dog11.biallelic_snps.hard_filtered.vcf.gz"
  local hard_filtered_tbi="${hard_filtered_vcf}.tbi"
  local phase_input_vcf="${PHASE_INPUT_ROOT}/${chromosome}.dog11.biallelic_snps.pass.callrate90.vcf.gz"
  local phase_input_tbi="${phase_input_vcf}.tbi"
  local stats_file="${STATS_ROOT}/${chromosome}.phase_input.bcftools.stats.txt"
  local chromosome_log="${LOG_ROOT}/${chromosome}.prepare_phase_input.log"
  local chromosome_tmp="${TMP_ROOT}/${chromosome}"
  local temporary_hard_filtered="${chromosome_tmp}/${chromosome}.hard_filtered.vcf.gz"
  local temporary_phase_input="${chromosome_tmp}/${chromosome}.phase_input.vcf.gz"
  local temporary_stats="${chromosome_tmp}/${chromosome}.phase_input.stats.txt"
  local chromosome_state="${STATE_DIR}/${chromosome}.status.tsv"
  local input_samples=0
  local output_samples=0
  local output_records=0
  local status="pending"

  mkdir -p "$chromosome_tmp"
  : > "$chromosome_log"

  {
    printf '[%(%F %T)T ] start chromosome: %s\n' -1 "$chromosome"
    printf '[%(%F %T)T ] input joint VCF: %s\n' -1 "$input_vcf"
    printf '[%(%F %T)T ] output phase-input VCF: %s\n' -1 "$phase_input_vcf"
  } >> "$chromosome_log"

  if [[ -s "$hard_filtered_vcf" && -s "$hard_filtered_tbi" && -s "$phase_input_vcf" && -s "$phase_input_tbi" && -s "$stats_file" ]]; then
    input_samples="$($BCFTOOLS_BIN query -l "$input_vcf" | wc -l)"
    output_samples="$($BCFTOOLS_BIN query -l "$phase_input_vcf" | wc -l)"
    output_records="$($BCFTOOLS_BIN view -H "$phase_input_vcf" | wc -l)"
    status="already_done"
  elif [[ -e "$hard_filtered_vcf" || -e "$hard_filtered_tbi" || -e "$phase_input_vcf" || -e "$phase_input_tbi" || -e "$stats_file" ]]; then
    status="failed_incomplete_existing_output"
  elif [[ -e "$temporary_hard_filtered" || -e "${temporary_hard_filtered}.tbi" || -e "$temporary_phase_input" || -e "${temporary_phase_input}.tbi" || -e "$temporary_stats" ]]; then
    status="failed_stale_temporary_output"
  elif [[ ! -s "$input_vcf" || ! -s "$input_tbi" ]]; then
    status="failed_missing_input"
  else
    input_samples="$($BCFTOOLS_BIN query -l "$input_vcf" | wc -l)"
    if [[ "$input_samples" -ne "$EXPECTED_SAMPLE_COUNT" ]]; then
      status="failed_sample_count"
    else
      printf '[%(%F %T)T ] Stage 1/4: normalize variants, retain only biallelic SNPs, remove old phase fields, and perform genotype/site quality control.\n' -1 >> "$chromosome_log"

      if "$BCFTOOLS_BIN" norm --threads "$THREADS" -f "$REFERENCE_FASTA" -m -any -Ou "$input_vcf" 2>> "$chromosome_log" \
        | "$BCFTOOLS_BIN" view --threads "$THREADS" -v snps -m2 -M2 -Ou 2>> "$chromosome_log" \
        | "$BCFTOOLS_BIN" annotate --threads "$THREADS" -x FORMAT/PGT,FORMAT/PID,FORMAT/PS -Ou 2>> "$chromosome_log" \
        | "$BCFTOOLS_BIN" +setGT -Ou -- -t a -n u 2>> "$chromosome_log" \
        | "$BCFTOOLS_BIN" +setGT -Ou -- -t q -n . -i "FMT/DP<${MIN_DP} || FMT/GQ<${MIN_GQ}" 2>> "$chromosome_log" \
        | "$BCFTOOLS_BIN" +fill-tags -Ou -- -t AC,AN,AF,F_MISSING 2>> "$chromosome_log" \
        | "$BCFTOOLS_BIN" filter --threads "$THREADS" -s HardFilter -e "$SITE_HARD_FILTER" -Oz -o "$temporary_hard_filtered" 2>> "$chromosome_log"; then

        printf '[%(%F %T)T ] stage 2/4: create hard-filtered VCF index.\n' -1 >> "$chromosome_log"
        if "$BCFTOOLS_BIN" index --threads "$THREADS" -t "$temporary_hard_filtered" >> "$chromosome_log" 2>&1; then
          printf '[%(%F %T)T ] Stage 3/4: retain SNPs with PASS status, call rate>=0.90, and remaining polymorphism in dog11.\n' -1 >> "$chromosome_log"
          if "$BCFTOOLS_BIN" view --threads "$THREADS" -f PASS \
              -i "INFO/F_MISSING<=${MAX_MISSING_RATE} && INFO/AC>0 && INFO/AC<INFO/AN" \
              -Oz -o "$temporary_phase_input" "$temporary_hard_filtered" >> "$chromosome_log" 2>&1 \
            && "$BCFTOOLS_BIN" index --threads "$THREADS" -t "$temporary_phase_input" >> "$chromosome_log" 2>&1; then

            printf '[%(%F %T)T ] stage 4/4: check final sample count, record count and unphased GT contract.\n' -1 >> "$chromosome_log"
            output_samples="$($BCFTOOLS_BIN query -l "$temporary_phase_input" | wc -l)"
            output_records="$($BCFTOOLS_BIN view -H "$temporary_phase_input" | wc -l)"
            phased_gt_count="$($BCFTOOLS_BIN query -f '[%GT\n]' "$temporary_phase_input" | awk 'index($0,"|") {n++} END {print n+0}')"

            if [[ "$output_samples" -ne "$EXPECTED_SAMPLE_COUNT" ]]; then
              status="failed_output_sample_count"
            elif [[ "$output_records" -le 0 ]]; then
              status="failed_empty_output"
            elif [[ "$phased_gt_count" -ne 0 ]]; then
              status="failed_phased_gt_remains"
            elif "$BCFTOOLS_BIN" stats "$temporary_phase_input" > "$temporary_stats" 2>> "$chromosome_log"; then
              mv "$temporary_hard_filtered" "$hard_filtered_vcf"
              mv "${temporary_hard_filtered}.tbi" "$hard_filtered_tbi"
              mv "$temporary_phase_input" "$phase_input_vcf"
              mv "${temporary_phase_input}.tbi" "$phase_input_tbi"
              mv "$temporary_stats" "$stats_file"
              status="run_success"
            else
              status="failed_stats_exit"
            fi
          else
            status="failed_phase_input_filter_exit"
          fi
        else
          status="failed_hard_filtered_index_exit"
        fi
      else
        status="failed_normalize_filter_exit"
      fi
    fi
  fi

  {
    printf '[%(%F %T)T ] chromosome status: %s; input_samples=%s; output_samples=%s; output_records=%s\n' \
      -1 "$status" "$input_samples" "$output_samples" "$output_records"
    printf '[%(%F %T)T ] single-chromosome task completed / done\n' -1
  } >> "$chromosome_log"

  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$chromosome" "$input_vcf" "$hard_filtered_vcf" "$phase_input_vcf" \
    "$output_samples" "$output_records" "$status" > "$chromosome_state"
}

if [[ $# -ne 0 ]]; then
  printf 'ERROR: this script does not accept arguments.\n' >&2
  exit 1
fi

mkdir -p "$OUT_DIR" "$HARD_FILTERED_ROOT" "$PHASE_INPUT_ROOT" "$STATS_ROOT" "$LOG_ROOT" "$TMP_ROOT" "$STATE_DIR"
: > "$LOG_FILE"

exec 9> "$LOCK_FILE"
if ! /usr/bin/flock -n 9; then
  fail "existing another items 1.5 task holds run lock: ${LOCK_FILE}"
fi

log_stage "start prepare dog11 all chromosome phase input."
log_stage "script path: ${SCRIPT_PATH}"
log_stage "input joint genotyping success table: ${SUCCESS_TSV}"
log_stage "output directory: ${OUT_DIR}"
log_stage "parallel chromosome count=${PARALLEL_CHROMOSOMES}; per-task bcftools threads=${THREADS}"
log_stage "genotype threshold: DP>=${MIN_DP}, GQ>=${MIN_GQ}; site call rate>=0.90"
log_stage "site hard filter expression: ${SITE_HARD_FILTER}"

require_file "$SUCCESS_TSV" "Step 4 joint genotyping success table"
require_file "$REFERENCE_FASTA" "CanFam3.1 FASTA"
require_file "$REFERENCE_FAI" "CanFam3.1 FASTA index"
[[ -x "$BCFTOOLS_BIN" ]] || fail "bcftools does not exist or is not executable: ${BCFTOOLS_BIN}"

expected_success_header=$'chromosome\tjoint_vcf\tstatus'
read -r success_header < "$SUCCESS_TSV"
[[ "$success_header" == "$expected_success_header" ]] || fail "Step 4 success table has an unexpected header."

rm -f "${STATE_DIR}"/*.status.tsv

declare -A INPUT_VCF_BY_CHROM=()
INPUT_COUNT=0
while IFS=$'\t' read -r chromosome joint_vcf status; do
  [[ "$status" == "run_success" || "$status" == "already_done" ]] || \
    fail "Step 4 success table contains a non-success status: ${chromosome}/${status}"
  [[ "$chromosome" =~ ^chr([1-9]|[12][0-9]|3[0-8])$ ]] || fail "chromosome name does not match chr1-chr38: ${chromosome}"
  [[ -z "${INPUT_VCF_BY_CHROM[$chromosome]+x}" ]] || fail "success table contains duplicate chromosome: ${chromosome}"
  require_file "$joint_vcf" "${chromosome} joint VCF"
  require_file "${joint_vcf}.tbi" "${chromosome} joint VCF index"
  INPUT_VCF_BY_CHROM["$chromosome"]="$joint_vcf"
  INPUT_COUNT=$((INPUT_COUNT + 1))
done < <(tail -n +2 "$SUCCESS_TSV")

[[ "$INPUT_COUNT" -eq "$EXPECTED_CHROMOSOME_COUNT" ]] || \
  fail "joint VCF count is not expected of ${EXPECTED_CHROMOSOME_COUNT}: actual=${INPUT_COUNT}"

mapfile -t CHROMOSOMES < <(awk -F '\t' '$1 ~ /^chr([1-9]|[12][0-9]|3[0-8])$/ {print $1}' "$REFERENCE_FAI")
[[ "${#CHROMOSOMES[@]}" -eq "$EXPECTED_CHROMOSOME_COUNT" ]] || \
  fail "reference FAI in chr1-chr38 unexpected count: actual=${#CHROMOSOMES[@]}"

SUBMITTED_COUNT=0
for chromosome in "${CHROMOSOMES[@]}"; do
  [[ -n "${INPUT_VCF_BY_CHROM[$chromosome]+x}" ]] || fail "success table missing chromosome: ${chromosome}"
  wait_for_slot
  SUBMITTED_COUNT=$((SUBMITTED_COUNT + 1))
  log_stage "submit chromosome ${SUBMITTED_COUNT}/${EXPECTED_CHROMOSOME_COUNT}: ${chromosome}"
  run_one_chromosome "$chromosome" "${INPUT_VCF_BY_CHROM[$chromosome]}" &
done

log_stage "all chromosome jobs submitted, waiting for background tasks to finish."
wait
log_stage "all chromosome tasks have returned, starting the summary status."

printf 'chromosome\tinput_joint_vcf\thard_filtered_vcf\tphase_input_vcf\toutput_samples\toutput_records\tstatus\n' > "$STATUS_TSV"
printf 'chromosome\tphase_input_vcf\toutput_records\tstatus\n' > "$SUCCESS_CHROM_TSV"
printf 'chromosome\tphase_input_vcf\tstatus\n' > "$FAILED_CHROM_TSV"

ALREADY_DONE_COUNT=0
RUN_SUCCESS_COUNT=0
FAILED_COUNT=0
STATUS_ROW_COUNT=0
TOTAL_OUTPUT_RECORDS=0

for chromosome in "${CHROMOSOMES[@]}"; do
  chromosome_state="${STATE_DIR}/${chromosome}.status.tsv"
  if [[ ! -s "$chromosome_state" ]]; then
    FAILED_COUNT=$((FAILED_COUNT + 1))
    printf '%s\t%s\t%s\n' "$chromosome" "${PHASE_INPUT_ROOT}/${chromosome}.dog11.biallelic_snps.pass.callrate90.vcf.gz" "failed_missing_state" >> "$FAILED_CHROM_TSV"
    continue
  fi

  IFS=$'\t' read -r state_chromosome input_vcf hard_filtered_vcf phase_input_vcf output_samples output_records status < "$chromosome_state"
  STATUS_ROW_COUNT=$((STATUS_ROW_COUNT + 1))
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$state_chromosome" "$input_vcf" "$hard_filtered_vcf" "$phase_input_vcf" "$output_samples" "$output_records" "$status" >> "$STATUS_TSV"

  case "$status" in
    already_done)
      ALREADY_DONE_COUNT=$((ALREADY_DONE_COUNT + 1))
      TOTAL_OUTPUT_RECORDS=$((TOTAL_OUTPUT_RECORDS + output_records))
      printf '%s\t%s\t%s\t%s\n' "$state_chromosome" "$phase_input_vcf" "$output_records" "$status" >> "$SUCCESS_CHROM_TSV"
      ;;
    run_success)
      RUN_SUCCESS_COUNT=$((RUN_SUCCESS_COUNT + 1))
      TOTAL_OUTPUT_RECORDS=$((TOTAL_OUTPUT_RECORDS + output_records))
      printf '%s\t%s\t%s\t%s\n' "$state_chromosome" "$phase_input_vcf" "$output_records" "$status" >> "$SUCCESS_CHROM_TSV"
      ;;
    *)
      FAILED_COUNT=$((FAILED_COUNT + 1))
      printf '%s\t%s\t%s\n' "$state_chromosome" "$phase_input_vcf" "$status" >> "$FAILED_CHROM_TSV"
      ;;
  esac
done

{
  printf 'key\tvalue\n'
  printf 'step\t%s\n' "1.5_prepare_dog11_phase_input"
  printf 'script_path\t%s\n' "$SCRIPT_PATH"
  printf 'input_success_tsv\t%s\n' "$SUCCESS_TSV"
  printf 'reference_fasta\t%s\n' "$REFERENCE_FASTA"
  printf 'output_dir\t%s\n' "$OUT_DIR"
  printf 'hard_filtered_root\t%s\n' "$HARD_FILTERED_ROOT"
  printf 'phase_input_root\t%s\n' "$PHASE_INPUT_ROOT"
  printf 'chromosome_count\t%s\n' "$EXPECTED_CHROMOSOME_COUNT"
  printf 'sample_count\t%s\n' "$EXPECTED_SAMPLE_COUNT"
  printf 'status_row_count\t%s\n' "$STATUS_ROW_COUNT"
  printf 'already_done_count\t%s\n' "$ALREADY_DONE_COUNT"
  printf 'run_success_count\t%s\n' "$RUN_SUCCESS_COUNT"
  printf 'failed_count\t%s\n' "$FAILED_COUNT"
  printf 'total_output_records\t%s\n' "$TOTAL_OUTPUT_RECORDS"
  printf 'parallel_chromosomes\t%s\n' "$PARALLEL_CHROMOSOMES"
  printf 'threads_per_chromosome\t%s\n' "$THREADS"
  printf 'min_dp\t%s\n' "$MIN_DP"
  printf 'min_gq\t%s\n' "$MIN_GQ"
  printf 'max_missing_rate\t%s\n' "$MAX_MISSING_RATE"
  printf 'site_hard_filter\t%s\n' "$SITE_HARD_FILTER"
  printf 'bcftools_bin\t%s\n' "$BCFTOOLS_BIN"
  printf 'lock_file\t%s\n' "$LOCK_FILE"
} > "$RUN_SUMMARY_TSV"

log_stage "write the status table: ${STATUS_TSV}"
log_stage "write success table: ${SUCCESS_CHROM_TSV}"
log_stage "write the failure table: ${FAILED_CHROM_TSV}"
log_stage "total chromosomes=${EXPECTED_CHROMOSOME_COUNT}; status_rows=${STATUS_ROW_COUNT}; already_done=${ALREADY_DONE_COUNT}; run_success=${RUN_SUCCESS_COUNT}; failed=${FAILED_COUNT}; final_records=${TOTAL_OUTPUT_RECORDS}"

if [[ "$FAILED_COUNT" -gt 0 || "$STATUS_ROW_COUNT" -ne "$EXPECTED_CHROMOSOME_COUNT" ]]; then
  log_stage "Some chromosome tasks failed or are missing; see ${FAILED_CHROM_TSV}"
  exit 1
fi

log_stage "task completed / done"
