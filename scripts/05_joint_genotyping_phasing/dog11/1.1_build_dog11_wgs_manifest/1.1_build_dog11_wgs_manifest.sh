#!/usr/bin/env bash


set -euo pipefail

SCRIPT_PATH="/path/to/repository/scripts/05_joint_genotyping_phasing/dog11/1.1_build_dog11_wgs_manifest/1.1_build_dog11_wgs_manifest.sh"
OUT_DIR="/path/to/project/analysis_20260709/5-6dogs/1.ASE-processing/1.1_build_dog11_wgs_manifest"
MANIFEST_TSV="${OUT_DIR}/dog11_wgs_manifest.tsv"
SUMMARY_TSV="${OUT_DIR}/dog11_wgs_manifest_summary.tsv"
PREFLIGHT_TSV="${OUT_DIR}/dog11_wgs_manifest_preflight.tsv"
RUN_SUMMARY_TSV="${OUT_DIR}/run_summary.tsv"
LOG_FILE="${OUT_DIR}/1.1_build_dog11_wgs_manifest.log"

REFERENCE_FASTA="/path/to/project/legacy/reference/GCF_000002285.3_CanFam3.1_genomic_chr.fna"
REFERENCE_FAI="${REFERENCE_FASTA}.fai"
REFERENCE_DICT="/path/to/project/legacy/reference/GCF_000002285.3_CanFam3.1_genomic_chr.dict"

FIVE_WGS_GATK_ROOT="/path/to/project/legacy/five-dogs/data/WGS/WGS_gatk"
FIVE_SINGLE_VCF_ROOT="/path/to/project/legacy/five-dogs/data/WGS/step5-WGS_gatk_norm-filter"
FIVE_WHATHAP_ROOT="/path/to/project/legacy/five-dogs/4.ASE-process/batch-v7/4.7haplotype-phasing"

SIX_WGS_GATK_ROOT="/path/to/project/doga_six_dogs/2.wgs-fastq/2.3WGS_gatk"
SIX_SINGLE_VCF_ROOT="/path/to/project/doga_six_dogs/2.wgs-fastq/2.5WGS_gatk_norm_filter"
SIX_WHATHAP_ROOT="/path/to/project/doga_six_dogs/3.ASE/8.haplotype-phasing"

FIVE_SAMPLES=(
  "Crak_358336"
  "CodRou_346627"
  "SanMcg_349891"
  "Dozer_386519"
  "Ace_396226"
)

SIX_SAMPLES=(
  "PI153"
  "BC1642"
  "RW589"
  "SLK1048"
  "RW356"
  "JPK344"
)

log_stage() {
  printf '[%(%F %T)T ] %s\n' -1 "$*" | tee -a "$LOG_FILE"
}

exists_flag() {
  local path="$1"

  if [[ -s "$path" ]]; then
    printf 'yes'
  else
    printf 'no'
  fi
}

command_exists_flag() {
  local command_name="$1"

  if command -v "$command_name" >/dev/null 2>&1; then
    printf 'yes'
  else
    printf 'no'
  fi
}

reference_sq_md5() {
  awk -F '\t' '$1 == "@SQ" { print $2 "\t" $3 }' "$REFERENCE_DICT" | md5sum | awk '{print $1}'
}

compute_bam_sq_md5() {
  local bam_path="$1"

  samtools view -H "$bam_path" 2>/dev/null | \
    awk -F '\t' '$1 == "@SQ" { print $2 "\t" $3 }' | \
    md5sum | awk '{print $1}'
}

extract_bam_rg_sm() {
  local bam_path="$1"

  samtools view -H "$bam_path" 2>/dev/null | \
    awk -F '\t' '$1 == "@RG" { for (i = 2; i <= NF; i++) if ($i ~ /^SM:/) print substr($i, 4) }' | \
    sort -u | paste -sd, -
}

samtools_quickcheck_flag() {
  local bam_path="$1"

  if samtools quickcheck -v "$bam_path" >/dev/null 2>&1; then
    printf 'yes'
  else
    printf 'no'
  fi
}

samtools_index_flag() {
  local bam_path="$1"

  if samtools idxstats "$bam_path" >/dev/null 2>&1; then
    printf 'yes'
  else
    printf 'no'
  fi
}

write_sample() {
  local project="$1"
  local sample_id="$2"
  local wgs_gatk_root="$3"
  local single_vcf_root="$4"
  local whatshap_root="$5"
  local recommended_role="gvcf_from_markdup_bam;old_vcf_and_whathap_for_qc_anchor_only"

  local markdup_bam="${wgs_gatk_root}/${sample_id}/${sample_id}.sorted.markdup.bam"
  local markdup_bai="${wgs_gatk_root}/${sample_id}/${sample_id}.sorted.markdup.bai"
  local markdup_metrics="${wgs_gatk_root}/${sample_id}/${sample_id}.markdup.metrics.txt"
  local single_sample_vcf="${single_vcf_root}/${sample_id}/${sample_id}.filtered.PASS.snp.norm.uniquePos.vcf.gz"
  local single_sample_vcf_tbi="${single_sample_vcf}.tbi"
  local whatshap_vcf="${whatshap_root}/${sample_id}/${sample_id}.phased.vcf.gz"
  local whatshap_vcf_tbi="${whatshap_vcf}.tbi"
  local whatshap_stats="${whatshap_root}/${sample_id}/${sample_id}.phasing.stats.txt"
  local whatshap_log="${whatshap_root}/${sample_id}/${sample_id}.phasing.log"

  local ref_exists
  local ref_fai_exists
  local ref_dict_exists
  local bam_exists
  local bai_exists
  local metrics_exists
  local single_vcf_exists
  local single_tbi_exists
  local whatshap_vcf_exists
  local whatshap_tbi_exists
  local whatshap_stats_exists
  local whatshap_log_exists
  local bam_quickcheck_pass="not_checked"
  local bam_index_readable="not_checked"
  local bam_rg_sm=""
  local bam_rg_sm_matches_sample="no"
  local bam_sq_md5=""
  local bam_sq_matches_reference="no"
  local preflight_status="failed_preflight"
  local ready_for_gvcf="no"
  local next_step="check_missing_markdup_or_reference_before_gvcf"

  ref_exists="$(exists_flag "$REFERENCE_FASTA")"
  ref_fai_exists="$(exists_flag "$REFERENCE_FAI")"
  ref_dict_exists="$(exists_flag "$REFERENCE_DICT")"
  bam_exists="$(exists_flag "$markdup_bam")"
  bai_exists="$(exists_flag "$markdup_bai")"
  metrics_exists="$(exists_flag "$markdup_metrics")"
  single_vcf_exists="$(exists_flag "$single_sample_vcf")"
  single_tbi_exists="$(exists_flag "$single_sample_vcf_tbi")"
  whatshap_vcf_exists="$(exists_flag "$whatshap_vcf")"
  whatshap_tbi_exists="$(exists_flag "$whatshap_vcf_tbi")"
  whatshap_stats_exists="$(exists_flag "$whatshap_stats")"
  whatshap_log_exists="$(exists_flag "$whatshap_log")"

  if [[ "$SAMTOOLS_AVAILABLE" == "yes" && "$bam_exists" == "yes" ]]; then
    bam_quickcheck_pass="$(samtools_quickcheck_flag "$markdup_bam")"
    bam_index_readable="$(samtools_index_flag "$markdup_bam")"
    bam_rg_sm="$(extract_bam_rg_sm "$markdup_bam" || true)"
    bam_sq_md5="$(compute_bam_sq_md5 "$markdup_bam" || true)"

    if [[ "$bam_rg_sm" == "$sample_id" ]]; then
      bam_rg_sm_matches_sample="yes"
    fi

    if [[ -n "$REFERENCE_SQ_MD5" && "$bam_sq_md5" == "$REFERENCE_SQ_MD5" ]]; then
      bam_sq_matches_reference="yes"
    fi
  fi

  if [[ "$ref_exists" == "yes" && "$ref_fai_exists" == "yes" && "$ref_dict_exists" == "yes" && \
    "$bam_exists" == "yes" && "$bai_exists" == "yes" && "$SAMTOOLS_AVAILABLE" == "yes" && \
    "$bam_quickcheck_pass" == "yes" && "$bam_index_readable" == "yes" && \
    "$bam_rg_sm_matches_sample" == "yes" && "$bam_sq_matches_reference" == "yes" ]]; then
    preflight_status="pass"
    ready_for_gvcf="yes"
    next_step="run_per_sample_HaplotypeCaller_ERC_GVCF_from_markdup_bam"
    READY_FOR_GVCF_COUNT=$((READY_FOR_GVCF_COUNT + 1))
  fi

  TOTAL_SAMPLE_COUNT=$((TOTAL_SAMPLE_COUNT + 1))

  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$project" "$sample_id" "$REFERENCE_FASTA" "$markdup_bam" "$markdup_bai" \
    "$markdup_metrics" "$single_sample_vcf" "$single_sample_vcf_tbi" \
    "$whatshap_vcf" "$whatshap_vcf_tbi" "$whatshap_stats" "$whatshap_log" \
    "$recommended_role" >> "$MANIFEST_TSV"

  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$project" "$sample_id" "$REFERENCE_FASTA" "$REFERENCE_FAI" "$REFERENCE_DICT" \
    "$markdup_bam" "$markdup_bai" "$ref_exists" "$ref_fai_exists" "$ref_dict_exists" \
    "$bam_exists" "$bai_exists" "$SAMTOOLS_AVAILABLE" "$bam_quickcheck_pass" \
    "$bam_index_readable" "$bam_rg_sm" "$bam_rg_sm_matches_sample" "$bam_sq_md5" \
    "$bam_sq_matches_reference" "$REFERENCE_SQ_MD5" "$preflight_status" "$ready_for_gvcf" >> "$PREFLIGHT_TSV"

  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$project" "$sample_id" "$ref_exists" "$ref_fai_exists" "$bam_exists" \
    "$bai_exists" "$metrics_exists" "$single_vcf_exists" "$single_tbi_exists" \
    "$whatshap_vcf_exists" "$whatshap_tbi_exists" "$whatshap_stats_exists" \
    "$whatshap_log_exists" "$ready_for_gvcf" "$next_step" >> "$SUMMARY_TSV"
}

mkdir -p "$OUT_DIR"
: > "$LOG_FILE"

log_stage "start dog11 WGS/WhatsHap manifest step one."
log_stage "script path: ${SCRIPT_PATH}"
log_stage "output directory: ${OUT_DIR}"

printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
  "project" "sample_id" "reference_fasta" "markdup_bam" "markdup_bai" \
  "markdup_metrics" "single_sample_vcf" "single_sample_vcf_tbi" \
  "whathap_vcf" "whathap_vcf_tbi" "whathap_stats" "whathap_log" \
  "recommended_role" > "$MANIFEST_TSV"

printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
  "project" "sample_id" "reference_fasta_exists" "reference_fai_exists" \
  "markdup_bam_exists" "markdup_bai_exists" "markdup_metrics_exists" \
  "single_sample_vcf_exists" "single_sample_vcf_tbi_exists" \
  "whathap_vcf_exists" "whathap_vcf_tbi_exists" "whathap_stats_exists" \
  "whathap_log_exists" "ready_for_gvcf" "recommended_next_step" > "$SUMMARY_TSV"

printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
  "project" "sample_id" "reference_fasta" "reference_fai" "reference_dict" \
  "markdup_bam" "markdup_bai" "reference_fasta_exists" "reference_fai_exists" "reference_dict_exists" \
  "markdup_bam_exists" "markdup_bai_exists" "samtools_available" "bam_quickcheck_pass" \
  "bam_index_readable" "bam_rg_sm" "bam_rg_sm_matches_sample" "bam_sq_md5" \
  "bam_sq_matches_reference" "reference_sq_md5" "preflight_status" "ready_for_gvcf" > "$PREFLIGHT_TSV"

TOTAL_SAMPLE_COUNT=0
READY_FOR_GVCF_COUNT=0
SAMTOOLS_AVAILABLE="$(command_exists_flag samtools)"
REFERENCE_SQ_MD5=""

if [[ -s "$REFERENCE_DICT" ]]; then
  REFERENCE_SQ_MD5="$(reference_sq_md5 || true)"
fi

log_stage "samtools available: ${SAMTOOLS_AVAILABLE}; reference sequence dictionary MD5: ${REFERENCE_SQ_MD5:-not_available}"

log_stage "write BarkBase 5 dogs path."
for sample_id in "${FIVE_SAMPLES[@]}"; do
  write_sample "BarkBase5" "$sample_id" "$FIVE_WGS_GATK_ROOT" "$FIVE_SINGLE_VCF_ROOT" "$FIVE_WHATHAP_ROOT"
done

log_stage "write DoGA 6 dogs path."
for sample_id in "${SIX_SAMPLES[@]}"; do
  write_sample "DoGA6" "$sample_id" "$SIX_WGS_GATK_ROOT" "$SIX_SINGLE_VCF_ROOT" "$SIX_WHATHAP_ROOT"
done

MISSING_FOR_GVCF_COUNT=$((TOTAL_SAMPLE_COUNT - READY_FOR_GVCF_COUNT))

{
  printf 'key\tvalue\n'
  printf 'step\t%s\n' "1.1_build_dog11_wgs_manifest"
  printf 'script_path\t%s\n' "$SCRIPT_PATH"
  printf 'output_dir\t%s\n' "$OUT_DIR"
  printf 'manifest_tsv\t%s\n' "$MANIFEST_TSV"
  printf 'summary_tsv\t%s\n' "$SUMMARY_TSV"
  printf 'preflight_tsv\t%s\n' "$PREFLIGHT_TSV"
  printf 'reference_fasta\t%s\n' "$REFERENCE_FASTA"
  printf 'reference_fai\t%s\n' "$REFERENCE_FAI"
  printf 'reference_dict\t%s\n' "$REFERENCE_DICT"
  printf 'reference_sq_md5\t%s\n' "$REFERENCE_SQ_MD5"
  printf 'samtools_available\t%s\n' "$SAMTOOLS_AVAILABLE"
  printf 'barkbase_sample_count\t%s\n' "${#FIVE_SAMPLES[@]}"
  printf 'doga_sample_count\t%s\n' "${#SIX_SAMPLES[@]}"
  printf 'total_sample_count\t%s\n' "$TOTAL_SAMPLE_COUNT"
  printf 'ready_for_gvcf_count\t%s\n' "$READY_FOR_GVCF_COUNT"
  printf 'missing_for_gvcf_count\t%s\n' "$MISSING_FOR_GVCF_COUNT"
} > "$RUN_SUMMARY_TSV"

log_stage "write manifest: ${MANIFEST_TSV}"
log_stage "write summary: ${SUMMARY_TSV}"
log_stage "write preflight: ${PREFLIGHT_TSV}"
log_stage "write run summary: ${RUN_SUMMARY_TSV}"
log_stage "total samples=${TOTAL_SAMPLE_COUNT}; ready_for_gvcf=${READY_FOR_GVCF_COUNT}; missing_for_gvcf=${MISSING_FOR_GVCF_COUNT}"
log_stage "task completed / done"
