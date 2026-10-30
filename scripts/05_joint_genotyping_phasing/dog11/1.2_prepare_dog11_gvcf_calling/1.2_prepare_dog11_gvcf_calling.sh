#!/usr/bin/env bash


set -euo pipefail

SCRIPT_PATH="/path/to/repository/scripts/05_joint_genotyping_phasing/dog11/1.2_prepare_dog11_gvcf_calling/1.2_prepare_dog11_gvcf_calling.sh"
STEP1_OUT_DIR="/path/to/project/analysis_20260709/5-6dogs/1.ASE-processing/1.1_build_dog11_wgs_manifest"
MANIFEST_TSV="${STEP1_OUT_DIR}/dog11_wgs_manifest.tsv"
SUMMARY_TSV="${STEP1_OUT_DIR}/dog11_wgs_manifest_summary.tsv"
PREFLIGHT_TSV="${STEP1_OUT_DIR}/dog11_wgs_manifest_preflight.tsv"

OUT_DIR="/path/to/project/analysis_20260709/5-6dogs/1.ASE-processing/1.2_prepare_dog11_gvcf_calling"
INPUT_TARGETS_TSV="${OUT_DIR}/dog11_gvcf_input_targets.tsv"
RUN_SUMMARY_TSV="${OUT_DIR}/run_summary.tsv"
LOG_FILE="${OUT_DIR}/1.2_prepare_dog11_gvcf_calling.log"

REFERENCE_DICT="/path/to/project/legacy/reference/GCF_000002285.3_CanFam3.1_genomic_chr.dict"

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

mkdir -p "$OUT_DIR"
: > "$LOG_FILE"

log_stage "start prepare dog11 GVCF calling input table."
log_stage "script path: ${SCRIPT_PATH}"
log_stage "input manifest: ${MANIFEST_TSV}"
log_stage "input summary: ${SUMMARY_TSV}"
log_stage "input preflight: ${PREFLIGHT_TSV}"
log_stage "output directory: ${OUT_DIR}"

require_file "$MANIFEST_TSV" "Step 1 manifest"
require_file "$SUMMARY_TSV" "Step 1 summary"
require_file "$PREFLIGHT_TSV" "Step 1 preflight"
require_file "$REFERENCE_DICT" "reference dict"

expected_manifest_header=$'project\tsample_id\treference_fasta\tmarkdup_bam\tmarkdup_bai\tmarkdup_metrics\tsingle_sample_vcf\tsingle_sample_vcf_tbi\twhathap_vcf\twhathap_vcf_tbi\twhathap_stats\twhathap_log\trecommended_role'
expected_summary_header=$'project\tsample_id\treference_fasta_exists\treference_fai_exists\tmarkdup_bam_exists\tmarkdup_bai_exists\tmarkdup_metrics_exists\tsingle_sample_vcf_exists\tsingle_sample_vcf_tbi_exists\twhathap_vcf_exists\twhathap_vcf_tbi_exists\twhathap_stats_exists\twhathap_log_exists\tready_for_gvcf\trecommended_next_step'
expected_preflight_header=$'project\tsample_id\treference_fasta\treference_fai\treference_dict\tmarkdup_bam\tmarkdup_bai\treference_fasta_exists\treference_fai_exists\treference_dict_exists\tmarkdup_bam_exists\tmarkdup_bai_exists\tsamtools_available\tbam_quickcheck_pass\tbam_index_readable\tbam_rg_sm\tbam_rg_sm_matches_sample\tbam_sq_md5\tbam_sq_matches_reference\treference_sq_md5\tpreflight_status\tready_for_gvcf'

read -r manifest_header < "$MANIFEST_TSV"
read -r summary_header < "$SUMMARY_TSV"
read -r preflight_header < "$PREFLIGHT_TSV"

if [[ "$manifest_header" != "$expected_manifest_header" ]]; then
  fail "manifest unexpected header."
fi

if [[ "$summary_header" != "$expected_summary_header" ]]; then
  fail "summary unexpected header."
fi

if [[ "$preflight_header" != "$expected_preflight_header" ]]; then
  fail "preflight unexpected header."
fi

declare -A READY_BY_SAMPLE
declare -A SUMMARY_PROJECT_BY_SAMPLE
declare -A PREFLIGHT_PROJECT_BY_SAMPLE
declare -A PREFLIGHT_STATUS_BY_SAMPLE
declare -A PREFLIGHT_READY_BY_SAMPLE

while IFS=$'\t' read -r project sample_id reference_fasta_exists reference_fai_exists markdup_bam_exists markdup_bai_exists markdup_metrics_exists single_sample_vcf_exists single_sample_vcf_tbi_exists whatshap_vcf_exists whatshap_vcf_tbi_exists whatshap_stats_exists whatshap_log_exists ready_for_gvcf recommended_next_step; do
  READY_BY_SAMPLE["$sample_id"]="$ready_for_gvcf"
  SUMMARY_PROJECT_BY_SAMPLE["$sample_id"]="$project"
done < <(tail -n +2 "$SUMMARY_TSV")

while IFS=$'\t' read -r preflight_project preflight_sample_id preflight_reference_fasta preflight_reference_fai preflight_reference_dict preflight_markdup_bam preflight_markdup_bai preflight_reference_fasta_exists preflight_reference_fai_exists preflight_reference_dict_exists preflight_markdup_bam_exists preflight_markdup_bai_exists preflight_samtools_available preflight_bam_quickcheck_pass preflight_bam_index_readable preflight_bam_rg_sm preflight_bam_rg_sm_matches_sample preflight_bam_sq_md5 preflight_bam_sq_matches_reference preflight_reference_sq_md5 preflight_status preflight_ready_for_gvcf; do
  PREFLIGHT_PROJECT_BY_SAMPLE["$preflight_sample_id"]="$preflight_project"
  PREFLIGHT_STATUS_BY_SAMPLE["$preflight_sample_id"]="$preflight_status"
  PREFLIGHT_READY_BY_SAMPLE["$preflight_sample_id"]="$preflight_ready_for_gvcf"
done < <(tail -n +2 "$PREFLIGHT_TSV")

printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
  "project" "sample_id" "reference_fasta" "reference_fai" "reference_dict" \
  "input_bam" "input_bai" "ready_for_gvcf" "preflight_status" "status" > "$INPUT_TARGETS_TSV"

TOTAL_SAMPLE_COUNT=0
READY_SAMPLE_COUNT=0
SKIPPED_SAMPLE_COUNT=0

while IFS=$'\t' read -r project sample_id reference_fasta markdup_bam markdup_bai markdup_metrics single_sample_vcf single_sample_vcf_tbi whatshap_vcf whatshap_vcf_tbi whatshap_stats whatshap_log recommended_role; do
  TOTAL_SAMPLE_COUNT=$((TOTAL_SAMPLE_COUNT + 1))

  reference_fai="${reference_fasta}.fai"
  ready_for_gvcf="${READY_BY_SAMPLE[$sample_id]:-no}"
  preflight_status="${PREFLIGHT_STATUS_BY_SAMPLE[$sample_id]:-missing}"
  status="ready"

  if [[ "${SUMMARY_PROJECT_BY_SAMPLE[$sample_id]:-}" != "$project" ]]; then
    status="summary_project_mismatch"
  elif [[ "${PREFLIGHT_PROJECT_BY_SAMPLE[$sample_id]:-}" != "$project" ]]; then
    status="preflight_project_mismatch"
  elif [[ "$ready_for_gvcf" != "yes" ]]; then
    status="not_ready_in_step1_summary"
  elif [[ "$preflight_status" != "pass" || "${PREFLIGHT_READY_BY_SAMPLE[$sample_id]:-}" != "yes" ]]; then
    status="not_passed_step1_preflight"
  elif [[ ! -s "$reference_fasta" || ! -s "$reference_fai" || ! -s "$REFERENCE_DICT" || ! -s "$markdup_bam" || ! -s "$markdup_bai" ]]; then
    status="missing_required_input"
  fi

  if [[ "$status" == "ready" ]]; then
    READY_SAMPLE_COUNT=$((READY_SAMPLE_COUNT + 1))
  else
    SKIPPED_SAMPLE_COUNT=$((SKIPPED_SAMPLE_COUNT + 1))
  fi

  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$project" "$sample_id" "$reference_fasta" "$reference_fai" "$REFERENCE_DICT" \
    "$markdup_bam" "$markdup_bai" "$ready_for_gvcf" "$preflight_status" "$status" >> "$INPUT_TARGETS_TSV"
done < <(tail -n +2 "$MANIFEST_TSV")

{
  printf 'key\tvalue\n'
  printf 'step\t%s\n' "1.2_prepare_dog11_gvcf_calling"
  printf 'script_path\t%s\n' "$SCRIPT_PATH"
  printf 'input_manifest_tsv\t%s\n' "$MANIFEST_TSV"
  printf 'input_summary_tsv\t%s\n' "$SUMMARY_TSV"
  printf 'input_preflight_tsv\t%s\n' "$PREFLIGHT_TSV"
  printf 'output_dir\t%s\n' "$OUT_DIR"
  printf 'input_targets_tsv\t%s\n' "$INPUT_TARGETS_TSV"
  printf 'total_sample_count\t%s\n' "$TOTAL_SAMPLE_COUNT"
  printf 'ready_sample_count\t%s\n' "$READY_SAMPLE_COUNT"
  printf 'skipped_sample_count\t%s\n' "$SKIPPED_SAMPLE_COUNT"
  printf 'does_not_run_gatk\t%s\n' "yes"
} > "$RUN_SUMMARY_TSV"

log_stage "write GVCF calling input table: ${INPUT_TARGETS_TSV}"
log_stage "total samples=${TOTAL_SAMPLE_COUNT}; ready=${READY_SAMPLE_COUNT}; skipped=${SKIPPED_SAMPLE_COUNT}"
log_stage "this step only prepare input table, does not define downstream execution logic, and does not run GATK."
log_stage "task completed / done"
