#!/usr/bin/env bash


set -euo pipefail
shopt -s nullglob

SCRIPT_PATH="/path/to/repository/scripts/05_joint_genotyping_phasing/dog11/1.4_joint_genotype_dog11/1.4_joint_genotype_dog11.sh"
STEP3_OUT_DIR="/path/to/project/analysis_20260709/5-6dogs/1.ASE-processing/1.3_run_dog11_gvcf_calling"
SUCCESS_TSV="${STEP3_OUT_DIR}/dog11_gvcf_batch_success.tsv"

REFERENCE_FASTA="/path/to/project/legacy/reference/GCF_000002285.3_CanFam3.1_genomic_chr.fna"
REFERENCE_FAI="/path/to/project/legacy/reference/GCF_000002285.3_CanFam3.1_genomic_chr.fna.fai"
REFERENCE_DICT="/path/to/project/legacy/reference/GCF_000002285.3_CanFam3.1_genomic_chr.dict"

OUT_DIR="/path/to/project/analysis_20260709/5-6dogs/1.ASE-processing/1.4_joint_genotype_dog11"
WORKSPACE_ROOT="${OUT_DIR}/genomicsdb"
VCF_ROOT="${OUT_DIR}/joint_vcf"
LOG_ROOT="${OUT_DIR}/logs"
TMP_ROOT="${OUT_DIR}/tmp"
STATE_DIR="${OUT_DIR}/state"
SAMPLE_MAP="${OUT_DIR}/dog11_gvcf_sample_map.tsv"
STATUS_TSV="${OUT_DIR}/dog11_joint_genotype_status.tsv"
SUCCESS_CHROM_TSV="${OUT_DIR}/dog11_joint_genotype_success.tsv"
FAILED_CHROM_TSV="${OUT_DIR}/dog11_joint_genotype_failed.tsv"
RUN_SUMMARY_TSV="${OUT_DIR}/run_summary.tsv"
LOG_FILE="${OUT_DIR}/1.4_joint_genotype_dog11.log"
LOCK_FILE="${OUT_DIR}/1.4_joint_genotype_dog11.lock"

GATK_ENV_BIN="/path/to/software/miniconda3/envs/dogs-harmony-latest/bin"
GATK_BIN="${GATK_ENV_BIN}/gatk"
BCFTOOLS_BIN="${GATK_ENV_BIN}/bcftools"
JAVA_BIN="${GATK_ENV_BIN}/java"
export PATH="${GATK_ENV_BIN}:${PATH}"

EXPECTED_SAMPLE_COUNT=11
EXPECTED_CHROMOSOME_COUNT=38
PARALLEL_CHROMOSOMES=20
READER_THREADS=5
JAVA_MEM="24g"

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
    sleep 30
  done
}

run_one_chromosome() {
  local chromosome="$1"
  local chromosome_log="${LOG_ROOT}/${chromosome}.joint_genotype.log"
  local chromosome_tmp="${TMP_ROOT}/${chromosome}"
  local building_workspace="${chromosome_tmp}/${chromosome}.genomicsdb.building"
  local final_workspace="${WORKSPACE_ROOT}/${chromosome}.genomicsdb"
  local temporary_vcf="${chromosome_tmp}/${chromosome}.dog11.joint.vcf.gz"
  local temporary_tbi="${temporary_vcf}.tbi"
  local final_vcf="${VCF_ROOT}/${chromosome}.dog11.joint.vcf.gz"
  local final_tbi="${final_vcf}.tbi"
  local chromosome_state="${STATE_DIR}/${chromosome}.status.tsv"
  local status="pending"

  mkdir -p "$chromosome_tmp"
  : > "$chromosome_log"

  {
    printf '[%(%F %T)T ] start chromosome: %s\n' -1 "$chromosome"
    printf '[%(%F %T)T ] input sample map: %s\n' -1 "$SAMPLE_MAP"
    printf '[%(%F %T)T ] output joint VCF: %s\n' -1 "$final_vcf"
  } >> "$chromosome_log"

  if [[ -s "$final_vcf" && -s "$final_tbi" && -d "$final_workspace" ]]; then
    status="already_done"
  elif [[ -e "$final_vcf" || -e "$final_tbi" || -e "$final_workspace" ]]; then
    status="failed_incomplete_existing_output"
    printf '[%(%F %T)T ] ERROR: final output is incomplete, please check manually this chromosome.\n' -1 >> "$chromosome_log"
  elif [[ -e "$building_workspace" || -e "$temporary_vcf" || -e "$temporary_tbi" ]]; then
    status="failed_stale_temporary_output"
    printf '[%(%F %T)T ] ERROR: temporary results from an interrupted previous run were found, please check manually.\n' -1 >> "$chromosome_log"
  elif "$GATK_BIN" --java-options "-Xmx${JAVA_MEM} -Djava.io.tmpdir=${chromosome_tmp}" GenomicsDBImport \
      --genomicsdb-workspace-path "$building_workspace" \
      --sample-name-map "$SAMPLE_MAP" \
      --intervals "$chromosome" \
      --reader-threads "$READER_THREADS" \
      --tmp-dir "$chromosome_tmp" >> "$chromosome_log" 2>&1; then
    printf '[%(%F %T)T ] GenomicsDBImport completed, start GenotypeGVCFs.\n' -1 >> "$chromosome_log"

    if "$GATK_BIN" --java-options "-Xmx${JAVA_MEM} -Djava.io.tmpdir=${chromosome_tmp}" GenotypeGVCFs \
        -R "$REFERENCE_FASTA" \
        -V "gendb://${building_workspace}" \
        -L "$chromosome" \
        -O "$temporary_vcf" \
        --create-output-variant-index true >> "$chromosome_log" 2>&1; then
      if [[ -s "$temporary_vcf" && -s "$temporary_tbi" && -d "$building_workspace" ]]; then
        mv "$building_workspace" "$final_workspace"
        mv "$temporary_vcf" "$final_vcf"
        mv "$temporary_tbi" "$final_tbi"
        status="run_success"
      else
        status="failed_missing_output_after_gatk"
      fi
    else
      status="failed_genotypegvcfs_exit"
    fi
  else
    status="failed_genomicsdbimport_exit"
  fi

  {
    printf '[%(%F %T)T ] chromosome status: %s\n' -1 "$status"
    printf '[%(%F %T)T ] single-chromosome task completed / done\n' -1
  } >> "$chromosome_log"

  printf '%s\t%s\t%s\t%s\t%s\n' \
    "$chromosome" "$final_workspace" "$final_vcf" "$final_tbi" "$status" > "$chromosome_state"
}

if [[ $# -ne 0 ]]; then
  printf 'ERROR: this script does not accept arguments.\n' >&2
  exit 1
fi

mkdir -p "$OUT_DIR" "$WORKSPACE_ROOT" "$VCF_ROOT" "$LOG_ROOT" "$TMP_ROOT" "$STATE_DIR"
: > "$LOG_FILE"

exec 9> "$LOCK_FILE"
if ! /usr/bin/flock -n 9; then
  fail "existing another items 1.4 task holds run lock: ${LOCK_FILE}"
fi

log_stage "start dog11 joint genotyping."
log_stage "script path: ${SCRIPT_PATH}"
log_stage "input gVCF success table: ${SUCCESS_TSV}"
log_stage "reference genome: ${REFERENCE_FASTA}"
log_stage "output directory: ${OUT_DIR}"
log_stage "parallel chromosome count: ${PARALLEL_CHROMOSOMES}; GenomicsDB reader threads: ${READER_THREADS}; single task Java memory: ${JAVA_MEM}"

require_file "$SUCCESS_TSV" "Step 3 gVCF success table"
require_file "$REFERENCE_FASTA" "CanFam3.1 FASTA"
require_file "$REFERENCE_FAI" "CanFam3.1 FASTA index"
require_file "$REFERENCE_DICT" "CanFam3.1 sequence dictionary"
[[ -x "$GATK_BIN" ]] || fail "GATK does not exist or is not executable: ${GATK_BIN}"
[[ -x "$BCFTOOLS_BIN" ]] || fail "bcftools does not exist or is not executable: ${BCFTOOLS_BIN}"
[[ -x "$JAVA_BIN" ]] || fail "Java does not exist or is not executable: ${JAVA_BIN}"

expected_success_header=$'project\tsample_id\toutput_gvcf\tbatch_status'
read -r success_header < "$SUCCESS_TSV"
[[ "$success_header" == "$expected_success_header" ]] || fail "Step 3 success table has an unexpected header."

rm -f "${STATE_DIR}"/*.status.tsv
: > "$SAMPLE_MAP"

SAMPLE_COUNT=0
declare -A SEEN_SAMPLE=()
while IFS=$'\t' read -r project sample_id output_gvcf batch_status; do
  [[ "$batch_status" == "run_success" || "$batch_status" == "already_done" ]] || \
    fail "success table contains non-success status: ${sample_id}/${batch_status}"
  [[ -z "${SEEN_SAMPLE[$sample_id]+x}" ]] || fail "success table contains duplicate sample: ${sample_id}"
  require_file "$output_gvcf" "${sample_id} gVCF"
  require_file "${output_gvcf}.tbi" "${sample_id} gVCF index"

  actual_sample="$($BCFTOOLS_BIN query -l "$output_gvcf")"
  [[ "$actual_sample" == "$sample_id" ]] || \
    fail "gVCF sample names do not match: in the table=${sample_id}; VCF in=${actual_sample}; file=${output_gvcf}"

  printf '%s\t%s\n' "$sample_id" "$output_gvcf" >> "$SAMPLE_MAP"
  SEEN_SAMPLE["$sample_id"]=1
  SAMPLE_COUNT=$((SAMPLE_COUNT + 1))
done < <(tail -n +2 "$SUCCESS_TSV")

[[ "$SAMPLE_COUNT" -eq "$EXPECTED_SAMPLE_COUNT" ]] || \
  fail "gVCF sample count differs from the expected ${EXPECTED_SAMPLE_COUNT}: actual=${SAMPLE_COUNT}"
log_stage "input check passed, total ${SAMPLE_COUNT} per-sample files gVCF; sample map: ${SAMPLE_MAP}"

mapfile -t CHROMOSOMES < <(
  awk -F '\t' '$1 ~ /^chr([1-9]|[12][0-9]|3[0-8])$/ {print $1}' "$REFERENCE_FAI"
)
[[ "${#CHROMOSOMES[@]}" -eq "$EXPECTED_CHROMOSOME_COUNT" ]] || \
  fail "reference FAI in chr1-chr38 unexpected count: actual=${#CHROMOSOMES[@]}"

SUBMITTED_COUNT=0
for chromosome in "${CHROMOSOMES[@]}"; do
  wait_for_slot
  SUBMITTED_COUNT=$((SUBMITTED_COUNT + 1))
  log_stage "submit chromosome ${SUBMITTED_COUNT}/${EXPECTED_CHROMOSOME_COUNT}: ${chromosome}"
  run_one_chromosome "$chromosome" &
done

log_stage "all chromosome jobs submitted, waiting for background tasks to finish."
wait
log_stage "all chromosome tasks have returned, starting the summary status."

printf 'chromosome\tgenomicsdb_workspace\tjoint_vcf\tjoint_vcf_tbi\tstatus\n' > "$STATUS_TSV"
printf 'chromosome\tjoint_vcf\tstatus\n' > "$SUCCESS_CHROM_TSV"
printf 'chromosome\tjoint_vcf\tstatus\n' > "$FAILED_CHROM_TSV"

ALREADY_DONE_COUNT=0
RUN_SUCCESS_COUNT=0
FAILED_COUNT=0
STATUS_ROW_COUNT=0

for chromosome in "${CHROMOSOMES[@]}"; do
  chromosome_state="${STATE_DIR}/${chromosome}.status.tsv"
  if [[ ! -s "$chromosome_state" ]]; then
    FAILED_COUNT=$((FAILED_COUNT + 1))
    printf '%s\t%s\t%s\n' "$chromosome" "${VCF_ROOT}/${chromosome}.dog11.joint.vcf.gz" "failed_missing_state" >> "$FAILED_CHROM_TSV"
    continue
  fi

  IFS=$'\t' read -r state_chromosome final_workspace final_vcf final_tbi status < "$chromosome_state"
  STATUS_ROW_COUNT=$((STATUS_ROW_COUNT + 1))
  printf '%s\t%s\t%s\t%s\t%s\n' "$state_chromosome" "$final_workspace" "$final_vcf" "$final_tbi" "$status" >> "$STATUS_TSV"

  case "$status" in
    already_done)
      ALREADY_DONE_COUNT=$((ALREADY_DONE_COUNT + 1))
      printf '%s\t%s\t%s\n' "$state_chromosome" "$final_vcf" "$status" >> "$SUCCESS_CHROM_TSV"
      ;;
    run_success)
      RUN_SUCCESS_COUNT=$((RUN_SUCCESS_COUNT + 1))
      printf '%s\t%s\t%s\n' "$state_chromosome" "$final_vcf" "$status" >> "$SUCCESS_CHROM_TSV"
      ;;
    *)
      FAILED_COUNT=$((FAILED_COUNT + 1))
      printf '%s\t%s\t%s\n' "$state_chromosome" "$final_vcf" "$status" >> "$FAILED_CHROM_TSV"
      ;;
  esac
done

{
  printf 'key\tvalue\n'
  printf 'step\t%s\n' "1.4_joint_genotype_dog11"
  printf 'script_path\t%s\n' "$SCRIPT_PATH"
  printf 'input_success_tsv\t%s\n' "$SUCCESS_TSV"
  printf 'reference_fasta\t%s\n' "$REFERENCE_FASTA"
  printf 'sample_map\t%s\n' "$SAMPLE_MAP"
  printf 'output_dir\t%s\n' "$OUT_DIR"
  printf 'genomicsdb_root\t%s\n' "$WORKSPACE_ROOT"
  printf 'joint_vcf_root\t%s\n' "$VCF_ROOT"
  printf 'sample_count\t%s\n' "$SAMPLE_COUNT"
  printf 'chromosome_count\t%s\n' "$EXPECTED_CHROMOSOME_COUNT"
  printf 'status_row_count\t%s\n' "$STATUS_ROW_COUNT"
  printf 'already_done_count\t%s\n' "$ALREADY_DONE_COUNT"
  printf 'run_success_count\t%s\n' "$RUN_SUCCESS_COUNT"
  printf 'failed_count\t%s\n' "$FAILED_COUNT"
  printf 'parallel_chromosomes\t%s\n' "$PARALLEL_CHROMOSOMES"
  printf 'reader_threads_per_chromosome\t%s\n' "$READER_THREADS"
  printf 'java_mem_per_chromosome\t%s\n' "$JAVA_MEM"
  printf 'gatk_bin\t%s\n' "$GATK_BIN"
  printf 'lock_file\t%s\n' "$LOCK_FILE"
} > "$RUN_SUMMARY_TSV"

log_stage "write the status table: ${STATUS_TSV}"
log_stage "write success table: ${SUCCESS_CHROM_TSV}"
log_stage "write the failure table: ${FAILED_CHROM_TSV}"
log_stage "total chromosomes=${EXPECTED_CHROMOSOME_COUNT}; status_rows=${STATUS_ROW_COUNT}; already_done=${ALREADY_DONE_COUNT}; run_success=${RUN_SUCCESS_COUNT}; failed=${FAILED_COUNT}"

if [[ "$FAILED_COUNT" -gt 0 || "$STATUS_ROW_COUNT" -ne "$EXPECTED_CHROMOSOME_COUNT" ]]; then
  log_stage "Some chromosome tasks failed or are missing; see ${FAILED_CHROM_TSV}"
  exit 1
fi

log_stage "task completed / done"
