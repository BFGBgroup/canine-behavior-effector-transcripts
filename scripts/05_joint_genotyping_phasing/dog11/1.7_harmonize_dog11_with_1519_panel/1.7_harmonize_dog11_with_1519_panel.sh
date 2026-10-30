#!/usr/bin/env bash


set -euo pipefail
shopt -s nullglob

SCRIPT_PATH="/path/to/repository/scripts/05_joint_genotyping_phasing/dog11/1.7_harmonize_dog11_with_1519_panel/1.7_harmonize_dog11_with_1519_panel.sh"
STEP5_SUCCESS="/path/to/project/analysis_20260709/5-6dogs/1.ASE-processing/1.5_prepare_dog11_phase_input/dog11_phase_input_success.tsv"
PANEL_DIR="/path/to/project/analysis_20260709/5-6dogs/1.ASE-processing/reference/1,519-CanFam3.1-WGS-26.7.14"
REFERENCE_FAI="/path/to/project/legacy/reference/GCF_000002285.3_CanFam3.1_genomic_chr.fna.fai"

OUT_DIR="/path/to/project/analysis_20260709/5-6dogs/1.ASE-processing/1.7_harmonize_dog11_with_1519_panel"
TARGET_ALL_ROOT="${OUT_DIR}/target_all_unique_bcf"
TARGET_SHARED_ROOT="${OUT_DIR}/target_shared_exact_bcf"
REFERENCE_SHARED_ROOT="${OUT_DIR}/reference_shared_exact_bcf"
DUPLICATE_ROOT="${OUT_DIR}/duplicate_positions"
METRICS_ROOT="${OUT_DIR}/metrics"
LOG_ROOT="${OUT_DIR}/logs"
TMP_ROOT="${OUT_DIR}/tmp"
STATE_DIR="${OUT_DIR}/state"
STATUS_TSV="${OUT_DIR}/dog11_1519_harmonization_status.tsv"
SUCCESS_TSV="${OUT_DIR}/dog11_1519_harmonization_success.tsv"
FAILED_TSV="${OUT_DIR}/dog11_1519_harmonization_failed.tsv"
RUN_SUMMARY_TSV="${OUT_DIR}/run_summary.tsv"
LOG_FILE="${OUT_DIR}/1.7_harmonize_dog11_with_1519_panel.log"
LOCK_FILE="${OUT_DIR}/1.7_harmonize_dog11_with_1519_panel.lock"

BCFTOOLS_BIN="/path/to/software/miniconda3/envs/dogs-harmony-latest/bin/bcftools"

EXPECTED_CHROMOSOME_COUNT=38
EXPECTED_TARGET_SAMPLES=11
EXPECTED_PANEL_SAMPLES=1519
PARALLEL_CHROMOSOMES=8
THREADS=2

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


write_duplicate_positions() {
  local input_vcf="$1"
  local output_tsv="$2"
  "$BCFTOOLS_BIN" query -f '%CHROM\t%POS\n' "$input_vcf" \
    | awk 'BEGIN {FS=OFS="\t"}
      NR==1 {previous_chrom=$1; previous_pos=$2; count=1; next}
      $1==previous_chrom && $2==previous_pos {count++; next}
      {if (count>1) print previous_chrom, previous_pos; previous_chrom=$1; previous_pos=$2; count=1}
      END {if (NR>0 && count>1) print previous_chrom, previous_pos}' > "$output_tsv"
}

write_unique_position_bcf() {
  local input_vcf="$1"
  local duplicate_tsv="$2"
  local output_bcf="$3"
  if [[ -s "$duplicate_tsv" ]]; then
    "$BCFTOOLS_BIN" view --threads "$THREADS" -T "^${duplicate_tsv}" -Ob -o "$output_bcf" "$input_vcf"
  else
    "$BCFTOOLS_BIN" view --threads "$THREADS" -Ob -o "$output_bcf" "$input_vcf"
  fi
  "$BCFTOOLS_BIN" index --threads "$THREADS" -f "$output_bcf"
}

run_one_chromosome() {
  local chromosome="$1"
  local target_vcf="$2"
  local panel_vcf="${PANEL_DIR}/ref-panel_${chromosome}_sample-snp_filltags_filter.phased.vcf.gz"
  local target_all_bcf="${TARGET_ALL_ROOT}/${chromosome}.dog11.unique_pos.bcf"
  local target_shared_bcf="${TARGET_SHARED_ROOT}/${chromosome}.dog11.shared_exact.bcf"
  local reference_shared_bcf="${REFERENCE_SHARED_ROOT}/${chromosome}.ref1519.shared_exact.bcf"
  local target_duplicate_tsv="${DUPLICATE_ROOT}/${chromosome}.dog11.duplicate_positions.tsv"
  local panel_duplicate_tsv="${DUPLICATE_ROOT}/${chromosome}.ref1519.duplicate_positions.tsv"
  local metrics_tsv="${METRICS_ROOT}/${chromosome}.harmonization_metrics.tsv"
  local chromosome_log="${LOG_ROOT}/${chromosome}.harmonization.log"
  local chromosome_tmp="${TMP_ROOT}/${chromosome}"
  local state_tsv="${STATE_DIR}/${chromosome}.status.tsv"
  local tmp_target_dup="${chromosome_tmp}/target.duplicate_positions.tsv"
  local tmp_panel_dup="${chromosome_tmp}/panel.duplicate_positions.tsv"
  local tmp_target_unique="${chromosome_tmp}/target.unique_pos.bcf"
  local tmp_panel_unique="${chromosome_tmp}/panel.unique_pos.bcf"
  local tmp_target_shared="${chromosome_tmp}/target.shared_exact.bcf"
  local tmp_panel_shared="${chromosome_tmp}/panel.shared_exact.bcf"
  local tmp_target_samepos="${chromosome_tmp}/target.shared_position.bcf"
  local tmp_metrics="${chromosome_tmp}/harmonization_metrics.tsv"
  local panel_match_input=""
  local target_original=0 target_duplicate=0 target_unique=0
  local panel_original=0 panel_duplicate=0 panel_unique=0
  local same_position=0 exact_shared=0 allele_conflict=0 target_only=0 panel_only=0
  local reference_shared_count=0 output_target_samples=0 output_panel_samples=0
  local target_samples=0 panel_samples=0 status="pending"

  mkdir -p "$chromosome_tmp"
  : > "$chromosome_log"
  {
    printf '[%(%F %T)T ] start chromosome: %s\n' -1 "$chromosome"
    printf '[%(%F %T)T ] dog11 input: %s\n' -1 "$target_vcf"
    printf '[%(%F %T)T ] 1,519 panel input: %s\n' -1 "$panel_vcf"
  } >> "$chromosome_log"

  if [[ -s "$target_all_bcf" && -s "${target_all_bcf}.csi" \
        && -s "$target_shared_bcf" && -s "${target_shared_bcf}.csi" \
        && -s "$reference_shared_bcf" && -s "${reference_shared_bcf}.csi" \
        && -e "$target_duplicate_tsv" && -e "$panel_duplicate_tsv" && -s "$metrics_tsv" ]]; then
    IFS=$'\t' read -r _ target_original target_duplicate target_unique panel_original panel_duplicate panel_unique \
      same_position exact_shared allele_conflict target_only panel_only < <(tail -n 1 "$metrics_tsv")
    status="already_done"
  elif [[ -e "$target_all_bcf" || -e "${target_all_bcf}.csi" || -e "$target_shared_bcf" \
        || -e "${target_shared_bcf}.csi" || -e "$reference_shared_bcf" || -e "${reference_shared_bcf}.csi" \
        || -e "$target_duplicate_tsv" || -e "$panel_duplicate_tsv" || -e "$metrics_tsv" ]]; then
    status="failed_incomplete_existing_output"
  elif [[ ! -s "$target_vcf" || ! -s "${target_vcf}.tbi" || ! -s "$panel_vcf" || ! -s "${panel_vcf}.tbi" ]]; then
    status="failed_missing_input_or_index"
  else
    target_samples="$($BCFTOOLS_BIN query -l "$target_vcf" | wc -l)"
    panel_samples="$($BCFTOOLS_BIN query -l "$panel_vcf" | wc -l)"
    if [[ "$target_samples" -ne "$EXPECTED_TARGET_SAMPLES" ]]; then
      status="failed_target_sample_count"
    elif [[ "$panel_samples" -ne "$EXPECTED_PANEL_SAMPLES" ]]; then
      status="failed_panel_sample_count"
    else
      printf '[%(%F %T)T ] stage 1/5: scan and record target with panel in duplicate coordinates.\n' -1 >> "$chromosome_log"
      if write_duplicate_positions "$target_vcf" "$tmp_target_dup" 2>> "$chromosome_log" \
        && write_duplicate_positions "$panel_vcf" "$tmp_panel_dup" 2>> "$chromosome_log"; then
        target_duplicate="$(wc -l < "$tmp_target_dup")"
        panel_duplicate="$(wc -l < "$tmp_panel_dup")"
        target_original="$($BCFTOOLS_BIN index -n "$target_vcf")"
        panel_original="$($BCFTOOLS_BIN index -n "$panel_vcf")"

        printf '[%(%F %T)T ] stage 2/5: entire site exclude duplicate coordinates and create unique position BCF.\n' -1 >> "$chromosome_log"
        if write_unique_position_bcf "$target_vcf" "$tmp_target_dup" "$tmp_target_unique" >> "$chromosome_log" 2>&1; then
          target_unique="$($BCFTOOLS_BIN index -n "$tmp_target_unique")"


          if [[ -s "$tmp_panel_dup" ]]; then
            if ! write_unique_position_bcf "$panel_vcf" "$tmp_panel_dup" "$tmp_panel_unique" >> "$chromosome_log" 2>&1; then
              status="failed_panel_unique_position_build_exit"
            else
              panel_match_input="$tmp_panel_unique"
              panel_unique="$($BCFTOOLS_BIN index -n "$tmp_panel_unique")"
            fi
          else
            panel_match_input="$panel_vcf"
            panel_unique="$panel_original"
            printf '[%(%F %T)T ] the reference panel has no duplicate coordinates, directly use the original panel VCF.\n' -1 >> "$chromosome_log"
          fi

          printf '[%(%F %T)T ] stage 3/5: by CHROM/POS/REF/ALT exact extract total has record.\n' -1 >> "$chromosome_log"
          if [[ -z "$panel_match_input" ]]; then
            :
          elif "$BCFTOOLS_BIN" isec --threads "$THREADS" -n=2 -w1 -Ob -o "$tmp_target_shared" \
              "$tmp_target_unique" "$panel_match_input" >> "$chromosome_log" 2>&1 \
            && "$BCFTOOLS_BIN" isec --threads "$THREADS" -n=2 -w1 -Ob -o "$tmp_panel_shared" \
              "$panel_match_input" "$tmp_target_unique" >> "$chromosome_log" 2>&1 \
            && "$BCFTOOLS_BIN" index --threads "$THREADS" -f "$tmp_target_shared" >> "$chromosome_log" 2>&1 \
            && "$BCFTOOLS_BIN" index --threads "$THREADS" -f "$tmp_panel_shared" >> "$chromosome_log" 2>&1; then
            exact_shared="$($BCFTOOLS_BIN index -n "$tmp_target_shared")"
            reference_shared_count="$($BCFTOOLS_BIN index -n "$tmp_panel_shared")"

            printf '[%(%F %T)T ] stage 4/5: by coordinates statistics total has, conflict and respective specific site.\n' -1 >> "$chromosome_log"
            if "$BCFTOOLS_BIN" isec --threads "$THREADS" -c all -n=2 -w1 -Ob -o "$tmp_target_samepos" \
                "$tmp_target_unique" "$panel_match_input" >> "$chromosome_log" 2>&1; then
              same_position="$($BCFTOOLS_BIN view -H "$tmp_target_samepos" | wc -l)"
              allele_conflict=$((same_position - exact_shared))
              target_only=$((target_unique - same_position))
              panel_only=$((panel_unique - same_position))

              printf '[%(%F %T)T ] stage 5/5: check output contract and release result.\n' -1 >> "$chromosome_log"
              output_target_samples="$($BCFTOOLS_BIN query -l "$tmp_target_shared" | wc -l)"
              output_panel_samples="$($BCFTOOLS_BIN query -l "$tmp_panel_shared" | wc -l)"
              if [[ "$exact_shared" -le 0 ]]; then
                status="failed_no_exact_shared_site"
              elif [[ "$exact_shared" -ne "$reference_shared_count" ]]; then
                status="failed_shared_record_count_mismatch"
              elif [[ "$output_target_samples" -ne "$EXPECTED_TARGET_SAMPLES" || "$output_panel_samples" -ne "$EXPECTED_PANEL_SAMPLES" ]]; then
                status="failed_output_sample_count"
              elif (( allele_conflict < 0 || target_only < 0 || panel_only < 0 )); then
                status="failed_negative_metric"
              else
                {
                  printf 'chromosome\ttarget_original_records\ttarget_duplicate_positions\ttarget_unique_records\tpanel_original_records\tpanel_duplicate_positions\tpanel_unique_records\tsame_position_records\texact_shared_records\tallele_conflict_positions\ttarget_only_positions\tpanel_only_positions\n'
                  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
                    "$chromosome" "$target_original" "$target_duplicate" "$target_unique" \
                    "$panel_original" "$panel_duplicate" "$panel_unique" "$same_position" \
                    "$exact_shared" "$allele_conflict" "$target_only" "$panel_only"
                } > "$tmp_metrics"
                mv "$tmp_target_unique" "$target_all_bcf"
                mv "${tmp_target_unique}.csi" "${target_all_bcf}.csi"
                mv "$tmp_target_shared" "$target_shared_bcf"
                mv "${tmp_target_shared}.csi" "${target_shared_bcf}.csi"
                mv "$tmp_panel_shared" "$reference_shared_bcf"
                mv "${tmp_panel_shared}.csi" "${reference_shared_bcf}.csi"
                mv "$tmp_target_dup" "$target_duplicate_tsv"
                mv "$tmp_panel_dup" "$panel_duplicate_tsv"
                mv "$tmp_metrics" "$metrics_tsv"
                rm -f "$tmp_panel_unique" "${tmp_panel_unique}.csi" "$tmp_target_samepos"
                status="run_success"
              fi
            else
              status="failed_position_intersection_exit"
            fi
          else
            status="failed_exact_intersection_exit"
          fi
        else
          status="failed_unique_position_build_exit"
        fi
      else
        status="failed_duplicate_position_scan_exit"
      fi
    fi
  fi

  {
    printf '[%(%F %T)T ] status=%s; target_unique=%s; panel_unique=%s; exact_shared=%s; allele_conflict=%s\n' \
      -1 "$status" "$target_unique" "$panel_unique" "$exact_shared" "$allele_conflict"
    printf '[%(%F %T)T ] single-chromosome task completed / done\n' -1
  } >> "$chromosome_log"
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$chromosome" "$target_original" "$target_duplicate" "$target_unique" "$panel_original" "$panel_duplicate" \
    "$panel_unique" "$same_position" "$exact_shared" "$allele_conflict" "$target_only" "$panel_only" \
    "$target_all_bcf" "$target_shared_bcf" "$reference_shared_bcf" "$metrics_tsv" "$status" > "$state_tsv"
}

if [[ $# -ne 0 ]]; then
  printf 'ERROR: this script does not accept arguments.\n' >&2
  exit 1
fi

mkdir -p "$OUT_DIR" "$TARGET_ALL_ROOT" "$TARGET_SHARED_ROOT" "$REFERENCE_SHARED_ROOT" \
  "$DUPLICATE_ROOT" "$METRICS_ROOT" "$LOG_ROOT" "$TMP_ROOT" "$STATE_DIR"
: > "$LOG_FILE"

exec 9> "$LOCK_FILE"
if ! /usr/bin/flock -n 9; then
  fail "existing another items 1.7 task holds run lock: ${LOCK_FILE}"
fi

log_stage "start harmonize dog11 with 1,519 dogs CanFam3.1 reference panel."
log_stage "script path: ${SCRIPT_PATH}"
log_stage "dog11 input success table: ${STEP5_SUCCESS}"
log_stage "reference panel directory: ${PANEL_DIR}"
log_stage "output directory: ${OUT_DIR}"
log_stage "parallel chromosome count=${PARALLEL_CHROMOSOMES}; per-task bcftools threads=${THREADS}"
log_stage "this step only prepare SHAPEIT5 input, does not run phasing, and does not read step 1.6 step WhatsHap result."

require_file "$STEP5_SUCCESS" "No. 1.5 step phase-input success table"
require_file "$REFERENCE_FAI" "CanFam3.1 FASTA index"
[[ -x "$BCFTOOLS_BIN" ]] || fail "bcftools does not exist or is not executable: ${BCFTOOLS_BIN}"

expected_header=$'chromosome\tphase_input_vcf\toutput_records\tstatus'
read -r observed_header < "$STEP5_SUCCESS"
[[ "$observed_header" == "$expected_header" ]] || fail "No. 1.5 step success table has an unexpected header."

rm -f "${STATE_DIR}"/*.status.tsv

declare -A TARGET_VCF_BY_CHROM=()
input_count=0
while IFS=$'\t' read -r chromosome target_vcf output_records input_status; do
  [[ "$input_status" == "run_success" || "$input_status" == "already_done" ]] || \
    fail "No. 1.5 step success table contains a non-success status: ${chromosome}/${input_status}"
  [[ "$chromosome" =~ ^chr([1-9]|[12][0-9]|3[0-8])$ ]] || fail "chromosome name does not match chr1-chr38: ${chromosome}"
  [[ -z "${TARGET_VCF_BY_CHROM[$chromosome]+x}" ]] || fail "No. 1.5 step success table contains duplicate chromosomes: ${chromosome}"
  require_file "$target_vcf" "${chromosome} dog11 phase-input VCF"
  require_file "${target_vcf}.tbi" "${chromosome} dog11 phase-input VCF index"
  TARGET_VCF_BY_CHROM["$chromosome"]="$target_vcf"
  input_count=$((input_count + 1))
done < <(tail -n +2 "$STEP5_SUCCESS")
[[ "$input_count" -eq "$EXPECTED_CHROMOSOME_COUNT" ]] || fail "dog11 input unexpected chromosome count: ${input_count}"

mapfile -t CHROMOSOMES < <(awk -F '\t' '$1 ~ /^chr([1-9]|[12][0-9]|3[0-8])$/ {print $1}' "$REFERENCE_FAI")
[[ "${#CHROMOSOMES[@]}" -eq "$EXPECTED_CHROMOSOME_COUNT" ]] || fail "reference FAI in chr1-chr38 unexpected count."

submitted=0
for chromosome in "${CHROMOSOMES[@]}"; do
  [[ -n "${TARGET_VCF_BY_CHROM[$chromosome]+x}" ]] || fail "dog11 success table missing chromosome: ${chromosome}"
  panel_vcf="${PANEL_DIR}/ref-panel_${chromosome}_sample-snp_filltags_filter.phased.vcf.gz"
  require_file "$panel_vcf" "${chromosome} 1,519 reference panel VCF"
  require_file "${panel_vcf}.tbi" "${chromosome} 1,519 reference panel VCF index"
  wait_for_slot
  submitted=$((submitted + 1))
  log_stage "submit chromosome ${submitted}/${EXPECTED_CHROMOSOME_COUNT}: ${chromosome}"
  run_one_chromosome "$chromosome" "${TARGET_VCF_BY_CHROM[$chromosome]}" &
done

log_stage "all chromosome jobs submitted, waiting for background tasks to finish."
wait
log_stage "all chromosome tasks have returned, starting the summary."

header=$'chromosome\ttarget_original_records\ttarget_duplicate_positions\ttarget_unique_records\tpanel_original_records\tpanel_duplicate_positions\tpanel_unique_records\tsame_position_records\texact_shared_records\tallele_conflict_positions\ttarget_only_positions\tpanel_only_positions\ttarget_all_unique_bcf\ttarget_shared_exact_bcf\treference_shared_exact_bcf\tmetrics_tsv\tstatus'
printf '%s\n' "$header" > "$STATUS_TSV"
printf 'chromosome\ttarget_all_unique_bcf\ttarget_shared_exact_bcf\treference_shared_exact_bcf\texact_shared_records\tstatus\n' > "$SUCCESS_TSV"
printf 'chromosome\tstatus\n' > "$FAILED_TSV"

run_success=0 already_done=0 failed=0 total_exact_shared=0 total_conflict=0
for chromosome in "${CHROMOSOMES[@]}"; do
  state_tsv="${STATE_DIR}/${chromosome}.status.tsv"
  if [[ ! -s "$state_tsv" ]]; then
    failed=$((failed + 1))
    printf '%s\tfailed_missing_state\n' "$chromosome" >> "$FAILED_TSV"
    continue
  fi
  state_line="$(<"$state_tsv")"
  printf '%s\n' "$state_line" >> "$STATUS_TSV"
  IFS=$'\t' read -r state_chr _ _ _ _ _ _ _ exact_shared allele_conflict _ _ target_all target_shared reference_shared _ status <<< "$state_line"
  case "$status" in
    run_success)
      run_success=$((run_success + 1))
      total_exact_shared=$((total_exact_shared + exact_shared))
      total_conflict=$((total_conflict + allele_conflict))
      printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$state_chr" "$target_all" "$target_shared" "$reference_shared" "$exact_shared" "$status" >> "$SUCCESS_TSV"
      ;;
    already_done)
      already_done=$((already_done + 1))
      total_exact_shared=$((total_exact_shared + exact_shared))
      total_conflict=$((total_conflict + allele_conflict))
      printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$state_chr" "$target_all" "$target_shared" "$reference_shared" "$exact_shared" "$status" >> "$SUCCESS_TSV"
      ;;
    *)
      failed=$((failed + 1))
      printf '%s\t%s\n' "$state_chr" "$status" >> "$FAILED_TSV"
      ;;
  esac
done

{
  printf 'metric\tvalue\n'
  printf 'expected_chromosomes\t%s\n' "$EXPECTED_CHROMOSOME_COUNT"
  printf 'run_success_chromosomes\t%s\n' "$run_success"
  printf 'already_done_chromosomes\t%s\n' "$already_done"
  printf 'failed_chromosomes\t%s\n' "$failed"
  printf 'total_exact_shared_records\t%s\n' "$total_exact_shared"
  printf 'total_allele_conflict_positions\t%s\n' "$total_conflict"
} > "$RUN_SUMMARY_TSV"

if [[ "$failed" -ne 0 || $((run_success + already_done)) -ne "$EXPECTED_CHROMOSOME_COUNT" ]]; then
  fail "harmonization did not succeed for all inputs: run_success=${run_success}; already_done=${already_done}; failed=${failed}"
fi

log_stage "all ${EXPECTED_CHROMOSOME_COUNT} records autosome harmonize completed."
log_stage "Total exactly shared records=${total_exact_shared}; Total positions with allele conflicts=${total_conflict}"
log_stage "1.7_harmonize_dog11_with_1519_panel done"
