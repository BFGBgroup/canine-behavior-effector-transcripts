#!/usr/bin/env bash


set -euo pipefail
shopt -s nullglob

if (( $# != 0 )); then
  printf 'ERROR: this script does not accept arguments, please direct run.\n' >&2
  exit 2
fi

SCRIPT_PATH="/path/to/repository/scripts/05_joint_genotyping_phasing/dog11/1.10_reconcile_dog11_whatsHap_shapeit5/1.10_reconcile_dog11_whatsHap_shapeit5.sh"
PYTHON_SCRIPT="/path/to/repository/scripts/05_joint_genotyping_phasing/dog11/1.10_reconcile_dog11_whatsHap_shapeit5/1.10_reconcile_dog11_whatsHap_shapeit5.py"

WHATSHAP_SUCCESS="/path/to/project/analysis_20260709/5-6dogs/1.ASE-processing/1.6_rerun_dog11_whatsHap/dog11_whatsHap_success.tsv"
WHATSHAP_ROOT="/path/to/project/analysis_20260709/5-6dogs/1.ASE-processing/1.6_rerun_dog11_whatsHap/phased_vcf"
SHAPEIT_SUCCESS="/path/to/project/analysis_20260709/5-6dogs/1.ASE-processing/1.9_phase_rare_dog11/dog11_phase_rare_success.tsv"
SHAPEIT_ROOT="/path/to/project/analysis_20260709/5-6dogs/1.ASE-processing/1.9_phase_rare_dog11/phased_common_rare_bcf"

OUT_DIR="/path/to/project/analysis_20260709/5-6dogs/1.ASE-processing/1.10_reconcile_dog11_whatsHap_shapeit5"
BCF_ROOT="${OUT_DIR}/reconciled_bcf"
BLOCK_QC_ROOT="${OUT_DIR}/block_qc"
CORRECTION_ROOT="${OUT_DIR}/corrections"
SAMPLE_SUMMARY_ROOT="${OUT_DIR}/sample_summaries"
METRICS_ROOT="${OUT_DIR}/metrics"
LOG_ROOT="${OUT_DIR}/logs"
TMP_ROOT="${OUT_DIR}/tmp"
STATE_DIR="${OUT_DIR}/state"
STATUS_TSV="${OUT_DIR}/dog11_whatsHap_shapeit5_reconcile_status.tsv"
SUCCESS_TSV="${OUT_DIR}/dog11_whatsHap_shapeit5_reconcile_success.tsv"
FAILED_TSV="${OUT_DIR}/dog11_whatsHap_shapeit5_reconcile_failed.tsv"
SAMPLE_SUMMARY_TSV="${OUT_DIR}/dog11_sample_chromosome_reconcile_summary.tsv"
RUN_SUMMARY_TSV="${OUT_DIR}/run_summary.tsv"
LOG_FILE="${OUT_DIR}/1.10_reconcile_dog11_whatsHap_shapeit5.log"
LOCK_FILE="${OUT_DIR}/1.10_reconcile_dog11_whatsHap_shapeit5.lock"

PYTHON_BIN="/path/to/software/miniconda3/envs/dogs-harmony-latest/bin/python3.12"
BCFTOOLS_BIN="/path/to/software/miniconda3/envs/dogs-harmony-latest/bin/bcftools"

EXPECTED_CHROMOSOMES=38
EXPECTED_SAMPLES=11
EXPECTED_WHATSHAP_TASKS=418
EXPECTED_TOTAL_RECORDS=9145272
PARALLEL_CHROMOSOMES=12

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

count_indexed_records() {
  local bcf="$1"
  "$BCFTOOLS_BIN" index -s "$bcf" \
    | awk -F '\t' '
        $3 ~ /^[0-9]+$/ { total += $3; contigs++ }
        END {
          if (contigs == 0) exit 1
          print total + 0
        }'
}

count_unphased_heterozygotes() {
  local bcf="$1"
  "$BCFTOOLS_BIN" query -f '[%GT\t]\n' "$bcf" \
    | awk '
        {
          for (i = 1; i <= NF; i++) {
            if ($i ~ /^[0-9]+\/[0-9]+$/) {
              split($i, allele, "/")
              if (allele[1] != allele[2]) count++
            }
          }
        }
        END { print count + 0 }'
}

write_state() {
  local state_file="$1"
  local chromosome="$2"
  local shapeit_bcf="$3"
  local output_bcf="$4"
  local expected_records="$5"
  local output_records="$6"
  local output_samples="$7"
  local unphased_hets="$8"
  local status="$9"
  local metrics="${10}"
  local block_qc="${11}"
  local corrections="${12}"

  {
    printf 'chromosome\tshapeit_bcf\toutput_bcf\texpected_records\toutput_records\toutput_samples\tunphased_heterozygotes\tmetrics_tsv\tblock_qc_tsv_gz\tcorrections_tsv_gz\tstatus\n'
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      "$chromosome" "$shapeit_bcf" "$output_bcf" "$expected_records" "$output_records" \
      "$output_samples" "$unphased_hets" "$metrics" "$block_qc" "$corrections" "$status"
  } > "${state_file}.tmp"
  mv "${state_file}.tmp" "$state_file"
}

run_chromosome() {
  local chromosome="$1"
  local chromosome_number="${chromosome#chr}"
  local shapeit_bcf="${SHAPEIT_ROOT}/${chromosome}.dog11.ref1519.common_rare.phased.bcf"
  local output_bcf="${BCF_ROOT}/${chromosome}.dog11.ref1519.whatsHap_reconciled.phased.bcf"
  local block_qc="${BLOCK_QC_ROOT}/${chromosome}.whatsHap_shapeit5.block_qc.tsv.gz"
  local corrections="${CORRECTION_ROOT}/${chromosome}.whatsHap_shapeit5.corrections.tsv.gz"
  local sample_summary="${SAMPLE_SUMMARY_ROOT}/${chromosome}.sample_reconcile_summary.tsv"
  local metrics="${METRICS_ROOT}/${chromosome}.reconcile_metrics.tsv"
  local runner_log="${LOG_ROOT}/${chromosome}.runner.log"
  local state_file="${STATE_DIR}/${chromosome}.status.tsv"
  local tmp_dir="${TMP_ROOT}/${chromosome}"
  local tmp_bcf="${tmp_dir}/${chromosome}.reconciled.tmp.bcf"
  local tmp_block_qc="${tmp_dir}/${chromosome}.block_qc.tmp.tsv.gz"
  local tmp_corrections="${tmp_dir}/${chromosome}.corrections.tmp.tsv.gz"
  local tmp_sample_summary="${tmp_dir}/${chromosome}.sample_summary.tmp.tsv"
  local tmp_metrics="${tmp_dir}/${chromosome}.metrics.tmp.tsv"
  local expected_records output_records=0 output_samples=0 unphased_hets=-1 status="pending"

  mkdir -p "$tmp_dir"
  : > "$runner_log"
  expected_records="$(awk -F '\t' -v chromosome="$chromosome" '
    NR > 1 && $1 == chromosome && $8 == "run_success" { print $6 }
  ' "$SHAPEIT_SUCCESS")"
  if [[ ! "$expected_records" =~ ^[1-9][0-9]*$ ]]; then
    status="failed_expected_records"
  elif [[ -s "$output_bcf" && -s "${output_bcf}.csi" && -s "$block_qc" && -s "$corrections" && -s "$sample_summary" && -s "$metrics" ]]; then
    output_records="$(count_indexed_records "$output_bcf" || printf '0')"
    output_samples="$("$BCFTOOLS_BIN" query -l "$output_bcf" | wc -l)"
    unphased_hets="$(count_unphased_heterozygotes "$output_bcf")"
    if [[ "$output_records" == "$expected_records" && "$output_samples" == "$EXPECTED_SAMPLES" && "$unphased_hets" == "0" ]]; then
      status="already_done"
    else
      status="failed_existing_output_contract"
    fi
  elif [[ -e "$output_bcf" || -e "${output_bcf}.csi" || -e "$block_qc" || -e "$corrections" || -e "$sample_summary" || -e "$metrics" ]]; then
    status="failed_incomplete_existing_output"
  else
    rm -f "$tmp_bcf" "${tmp_bcf}.csi" "$tmp_block_qc" "$tmp_corrections" "$tmp_sample_summary" "$tmp_metrics"
    printf '[%(%F %T)T ] stage 1/4: read%s of 11 items WhatsHap VCF and orient PS block.\n' -1 "$chromosome" >> "$runner_log"
    if ! "$PYTHON_BIN" "$PYTHON_SCRIPT" \
      --chromosome "$chromosome" \
      --shapeit-bcf "$shapeit_bcf" \
      --whatshap-root "$WHATSHAP_ROOT" \
      --output-bcf "$tmp_bcf" \
      --block-qc "$tmp_block_qc" \
      --corrections "$tmp_corrections" \
      --sample-summary "$tmp_sample_summary" \
      --metrics "$tmp_metrics" >> "$runner_log" 2>&1; then
      status="failed_reconcile"
    fi

    if [[ "$status" == "pending" ]]; then
      printf '[%(%F %T)T ] stage 2/4: create temporary BCF index.\n' -1 >> "$runner_log"
      if ! "$BCFTOOLS_BIN" index --threads 1 -f "$tmp_bcf" >> "$runner_log" 2>&1; then
        status="failed_index"
      fi
    fi

    if [[ "$status" == "pending" ]]; then
      printf '[%(%F %T)T ] stage 3/4: check record count, sample count and unphased heterozygous GT.\n' -1 >> "$runner_log"
      output_records="$(count_indexed_records "$tmp_bcf" || printf '0')"
      output_samples="$("$BCFTOOLS_BIN" query -l "$tmp_bcf" | wc -l)"
      unphased_hets="$(count_unphased_heterozygotes "$tmp_bcf")"
      if [[ "$output_records" != "$expected_records" ]]; then
        status="failed_output_records"
      elif [[ "$output_samples" != "$EXPECTED_SAMPLES" ]]; then
        status="failed_output_samples"
      elif [[ "$unphased_hets" != "0" ]]; then
        status="failed_unphased_heterozygotes"
      fi
    fi

    if [[ "$status" == "pending" ]]; then
      printf '[%(%F %T)T ] stage 4/4: release BCF and QC result.\n' -1 >> "$runner_log"
      mv "$tmp_block_qc" "$block_qc"
      mv "$tmp_corrections" "$corrections"
      mv "$tmp_sample_summary" "$sample_summary"
      mv "$tmp_metrics" "$metrics"
      mv "$tmp_bcf" "$output_bcf"
      mv "${tmp_bcf}.csi" "${output_bcf}.csi"
      status="run_success"
      printf '[%(%F %T)T ] %s harmonize completed / done.\n' -1 "$chromosome" >> "$runner_log"
    fi
  fi

  write_state "$state_file" "$chromosome" "$shapeit_bcf" "$output_bcf" "$expected_records" \
    "$output_records" "$output_samples" "$unphased_hets" "$status" "$metrics" "$block_qc" "$corrections"
  [[ "$status" == "run_success" || "$status" == "already_done" ]]
}

mkdir -p "$OUT_DIR" "$BCF_ROOT" "$BLOCK_QC_ROOT" "$CORRECTION_ROOT" \
  "$SAMPLE_SUMMARY_ROOT" "$METRICS_ROOT" "$LOG_ROOT" "$TMP_ROOT" "$STATE_DIR"

exec 9> "$LOCK_FILE"
flock -n 9 || fail "detected another items No. 1.10 step task positive in run: ${LOCK_FILE}"
: > "$LOG_FILE"

log_stage "Starting correction of SHAPEIT5 autosomal phasing using local WhatsHap PS blocks."
log_stage "script path: ${SCRIPT_PATH}"
log_stage "Python tool: ${PYTHON_SCRIPT}"
log_stage "WhatsHap input directory: ${WHATSHAP_ROOT}"
log_stage "SHAPEIT5 input directory: ${SHAPEIT_ROOT}"
log_stage "output directory: ${OUT_DIR}"
log_stage "parallel chromosome upper limit=${PARALLEL_CHROMOSOMES}"

require_file "$PYTHON_SCRIPT" "Python tool"
require_file "$WHATSHAP_SUCCESS" "No. 1.6 step success table"
require_file "$SHAPEIT_SUCCESS" "No. 1.9 step success table"
[[ -x "$PYTHON_BIN" ]] || fail "Python not executable: ${PYTHON_BIN}"
[[ -x "$BCFTOOLS_BIN" ]] || fail "bcftools not executable: ${BCFTOOLS_BIN}"

whatshap_success_count="$(awk -F '\t' 'NR > 1 && $6 == "run_success" { count++ } END { print count + 0 }' "$WHATSHAP_SUCCESS")"
shapeit_success_count="$(awk -F '\t' 'NR > 1 && $8 == "run_success" { count++ } END { print count + 0 }' "$SHAPEIT_SUCCESS")"
[[ "$whatshap_success_count" == "$EXPECTED_WHATSHAP_TASKS" ]] \
  || fail "Step 1.6 successful-task count is not ${EXPECTED_WHATSHAP_TASKS}: ${whatshap_success_count}"
[[ "$shapeit_success_count" == "$EXPECTED_CHROMOSOMES" ]] \
  || fail "Step 1.9 successful-chromosome count is not ${EXPECTED_CHROMOSOMES}: ${shapeit_success_count}"

for chromosome_number in $(seq 1 "$EXPECTED_CHROMOSOMES"); do
  chromosome="chr${chromosome_number}"
  shapeit_bcf="${SHAPEIT_ROOT}/${chromosome}.dog11.ref1519.common_rare.phased.bcf"
  require_file "$shapeit_bcf" "${chromosome} SHAPEIT5 BCF"
  require_file "${shapeit_bcf}.csi" "${chromosome} SHAPEIT5 BCF index"
  sample_rows="$(awk -F '\t' -v chromosome="$chromosome" '
    NR > 1 && $2 == chromosome && $6 == "run_success" { count++ }
    END { print count + 0 }
  ' "$WHATSHAP_SUCCESS")"
  [[ "$sample_rows" == "$EXPECTED_SAMPLES" ]] \
    || fail "${chromosome} of WhatsHap success sample count is not ${EXPECTED_SAMPLES}: ${sample_rows}"
done
log_stage "input preflight check completed: 38 records autosome, each 11 items WhatsHap sample all has success record."

active_jobs=0
for chromosome_number in $(seq 1 "$EXPECTED_CHROMOSOMES"); do
  chromosome="chr${chromosome_number}"
  log_stage "submit chromosome ${chromosome_number}/${EXPECTED_CHROMOSOMES}: ${chromosome}"
  run_chromosome "$chromosome" &
  active_jobs=$((active_jobs + 1))
  if (( active_jobs >= PARALLEL_CHROMOSOMES )); then
    wait -n || true
    active_jobs=$((active_jobs - 1))
  fi
done
while (( active_jobs > 0 )); do
  wait -n || true
  active_jobs=$((active_jobs - 1))
done
log_stage "all chromosome tasks have returned, starting the summary."

{
  printf 'chromosome\tshapeit_bcf\toutput_bcf\texpected_records\toutput_records\toutput_samples\tunphased_heterozygotes\tmetrics_tsv\tblock_qc_tsv_gz\tcorrections_tsv_gz\tstatus\n'
  for chromosome_number in $(seq 1 "$EXPECTED_CHROMOSOMES"); do
    state_file="${STATE_DIR}/chr${chromosome_number}.status.tsv"
    if [[ -s "$state_file" ]]; then
      tail -n 1 "$state_file"
    else
      printf 'chr%s\tNA\tNA\t0\t0\t0\t-1\tNA\tNA\tNA\tfailed_missing_state\n' "$chromosome_number"
    fi
  done
} > "${STATUS_TSV}.tmp"
mv "${STATUS_TSV}.tmp" "$STATUS_TSV"

awk -F '\t' 'BEGIN { OFS = FS } NR == 1 || $11 == "run_success" || $11 == "already_done"' \
  "$STATUS_TSV" > "${SUCCESS_TSV}.tmp"
mv "${SUCCESS_TSV}.tmp" "$SUCCESS_TSV"
awk -F '\t' 'BEGIN { OFS = FS } NR == 1 || ($11 != "run_success" && $11 != "already_done")' \
  "$STATUS_TSV" > "${FAILED_TSV}.tmp"
mv "${FAILED_TSV}.tmp" "$FAILED_TSV"

failed_count="$(awk -F '\t' 'NR > 1 && $11 != "run_success" && $11 != "already_done" { count++ } END { print count + 0 }' "$STATUS_TSV")"
run_success_count="$(awk -F '\t' 'NR > 1 && $11 == "run_success" { count++ } END { print count + 0 }' "$STATUS_TSV")"
already_done_count="$(awk -F '\t' 'NR > 1 && $11 == "already_done" { count++ } END { print count + 0 }' "$STATUS_TSV")"
if [[ "$failed_count" != "0" ]]; then
  fail "has ${failed_count} chromosomes harmonize failed, see: ${FAILED_TSV}"
fi

{
  head -n 1 "${SAMPLE_SUMMARY_ROOT}/chr1.sample_reconcile_summary.tsv"
  for chromosome_number in $(seq 1 "$EXPECTED_CHROMOSOMES"); do
    tail -n +2 "${SAMPLE_SUMMARY_ROOT}/chr${chromosome_number}.sample_reconcile_summary.tsv"
  done
} > "${SAMPLE_SUMMARY_TSV}.tmp"
mv "${SAMPLE_SUMMARY_TSV}.tmp" "$SAMPLE_SUMMARY_TSV"

total_output_records="$(awk -F '\t' 'NR > 1 { total += $3 } END { print total + 0 }' "${METRICS_ROOT}"/*.reconcile_metrics.tsv)"
total_whatshap_hets="$(awk -F '\t' 'NR > 1 { total += $5 } END { print total + 0 }' "${METRICS_ROOT}"/*.reconcile_metrics.tsv)"
total_ps_blocks="$(awk -F '\t' 'NR > 1 { total += $6 } END { print total + 0 }' "${METRICS_ROOT}"/*.reconcile_metrics.tsv)"
total_matched_hets="$(awk -F '\t' 'NR > 1 { total += $7 } END { print total + 0 }' "${METRICS_ROOT}"/*.reconcile_metrics.tsv)"
total_genotype_mismatches="$(awk -F '\t' 'NR > 1 { total += $8 } END { print total + 0 }' "${METRICS_ROOT}"/*.reconcile_metrics.tsv)"
total_oriented_blocks="$(awk -F '\t' 'NR > 1 { total += $9 } END { print total + 0 }' "${METRICS_ROOT}"/*.reconcile_metrics.tsv)"
total_ambiguous_blocks="$(awk -F '\t' 'NR > 1 { total += $12 } END { print total + 0 }' "${METRICS_ROOT}"/*.reconcile_metrics.tsv)"
total_no_anchor_blocks="$(awk -F '\t' 'NR > 1 { total += $13 } END { print total + 0 }' "${METRICS_ROOT}"/*.reconcile_metrics.tsv)"
total_discordant_sites="$(awk -F '\t' 'NR > 1 { total += $14 } END { print total + 0 }' "${METRICS_ROOT}"/*.reconcile_metrics.tsv)"
total_ps_annotated_sites="$(awk -F '\t' 'NR > 1 { total += $15 } END { print total + 0 }' "${METRICS_ROOT}"/*.reconcile_metrics.tsv)"
total_corrected_sites="$(awk -F '\t' 'NR > 1 { total += $16 } END { print total + 0 }' "${METRICS_ROOT}"/*.reconcile_metrics.tsv)"
[[ "$total_output_records" == "$EXPECTED_TOTAL_RECORDS" ]] \
  || fail "final record total count is not ${EXPECTED_TOTAL_RECORDS}: ${total_output_records}"

{
  printf 'metric\tvalue\n'
  printf 'expected_chromosomes\t%s\n' "$EXPECTED_CHROMOSOMES"
  printf 'run_success_chromosomes\t%s\n' "$run_success_count"
  printf 'already_done_chromosomes\t%s\n' "$already_done_count"
  printf 'failed_chromosomes\t%s\n' "$failed_count"
  printf 'output_records\t%s\n' "$total_output_records"
  printf 'whatshap_phased_heterozygotes\t%s\n' "$total_whatshap_hets"
  printf 'whatshap_ps_blocks\t%s\n' "$total_ps_blocks"
  printf 'matched_heterozygotes\t%s\n' "$total_matched_hets"
  printf 'genotype_mismatches\t%s\n' "$total_genotype_mismatches"
  printf 'oriented_blocks\t%s\n' "$total_oriented_blocks"
  printf 'ambiguous_tie_blocks\t%s\n' "$total_ambiguous_blocks"
  printf 'no_anchor_blocks\t%s\n' "$total_no_anchor_blocks"
  printf 'discordant_sites_after_block_orientation\t%s\n' "$total_discordant_sites"
  printf 'ps_annotated_sites\t%s\n' "$total_ps_annotated_sites"
  printf 'corrected_sites\t%s\n' "$total_corrected_sites"
  printf 'parallel_chromosomes\t%s\n' "$PARALLEL_CHROMOSOMES"
} > "${RUN_SUMMARY_TSV}.tmp"
mv "${RUN_SUMMARY_TSV}.tmp" "$RUN_SUMMARY_TSV"

log_stage "all 38 records autosome harmonize completed."
log_stage "total output sites=${total_output_records}; WhatsHap PS block=${total_ps_blocks}; orientable blocks=${total_oriented_blocks}"
log_stage "conflicting sites corrected within blocks=${total_corrected_sites}; tied blocks=${total_ambiguous_blocks}; blocks without anchors=${total_no_anchor_blocks}"
log_stage "1.10_reconcile_dog11_whatsHap_shapeit5 done"
