#!/usr/bin/env bash


set -euo pipefail
shopt -s nullglob

SCRIPT_PATH="/path/to/repository/scripts/05_joint_genotyping_phasing/dog11/1.8_phase_common_dog11_with_1519_panel/1.8_phase_common_dog11_with_1519_panel.sh"
HARMONIZATION_SUCCESS="/path/to/project/analysis_20260709/5-6dogs/1.ASE-processing/1.7_harmonize_dog11_with_1519_panel/dog11_1519_harmonization_success.tsv"
MAP_DIR="/path/to/project/analysis_20260709/5-6dogs/1.ASE-processing/reference/Campbell2016_CanFam3.1_genetic_map/1.convert_Campbell2016_to_SHAPEIT5/shapeit5_maps"

OUT_DIR="/path/to/project/analysis_20260709/5-6dogs/1.ASE-processing/1.8_phase_common_dog11_with_1519_panel"
COMMON_TARGET_ROOT="${OUT_DIR}/common_target_bcf"
COMMON_REFERENCE_ROOT="${OUT_DIR}/common_reference_bcf"
PHASED_ROOT="${OUT_DIR}/phased_common_bcf"
METRICS_ROOT="${OUT_DIR}/metrics"
LOG_ROOT="${OUT_DIR}/logs"
TMP_ROOT="${OUT_DIR}/tmp"
STATE_DIR="${OUT_DIR}/state"
STATUS_TSV="${OUT_DIR}/dog11_ref1519_phase_common_status.tsv"
SUCCESS_TSV="${OUT_DIR}/dog11_ref1519_phase_common_success.tsv"
FAILED_TSV="${OUT_DIR}/dog11_ref1519_phase_common_failed.tsv"
RUN_SUMMARY_TSV="${OUT_DIR}/run_summary.tsv"
LOG_FILE="${OUT_DIR}/1.8_phase_common_dog11_with_1519_panel.log"
LOCK_FILE="${OUT_DIR}/1.8_phase_common_dog11_with_1519_panel.lock"

SHAPEIT5_BIN="/path/to/software/environments/shapeit5/bin/SHAPEIT5_phase_common"
BCFTOOLS_BIN="/path/to/software/miniconda3/envs/dogs-harmony-latest/bin/bcftools"

EXPECTED_CHROMOSOMES=38
EXPECTED_TARGET_SAMPLES=11
EXPECTED_REFERENCE_SAMPLES=1519
COMMON_MAF="0.001"
COMMON_MAF_LABEL="${COMMON_MAF//./p}"
PARALLEL_CHROMOSOMES=40
THREADS=1
BCFTOOLS_THREADS=0
SEED=15052011

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

wait_for_slot() {
  while [[ "$(jobs -pr | wc -l)" -ge "$PARALLEL_CHROMOSOMES" ]]; do
    sleep 20
  done
}

validate_map() {
  local map_file="$1"
  awk '
    BEGIN { FS = "\t" }
    NR == 1 {
      if (NF != 3 || $1 != "position" || $2 != "COMBINED_rate(cM/Mb)" || $3 != "Genetic_Map(cM)") exit 10
      next
    }
    {
      if (NF != 3 || $1 !~ /^[0-9]+$/) exit 11
      if ($2 !~ /^[0-9]+([.][0-9]+)?([eE][-+]?[0-9]+)?$/ || $3 !~ /^[0-9]+([.][0-9]+)?([eE][-+]?[0-9]+)?$/) exit 12
      if ($2 < 0 || $3 < 0) exit 13
      if (seen && ($1 <= previous_position || $3 < previous_cm)) exit 14
      previous_position = $1
      previous_cm = $3
      seen = 1
      rows++
    }
    END { if (rows < 2) exit 15 }
  ' "$map_file"
}

has_info_tag() {
  local bcf="$1"
  local tag="$2"
  "$BCFTOOLS_BIN" view -h "$bcf" | grep -F "##INFO=<ID=${tag}," > /dev/null
}

count_unphased_heterozygotes() {
  local bcf="$1"
  "$BCFTOOLS_BIN" query -f '[%GT\t]\n' "$bcf" \
    | awk '{
        for (i = 1; i <= NF; i++) {
          if ($i ~ /^[0-9]+\/[0-9]+$/) {
            split($i, allele, "/")
            if (allele[1] != allele[2]) count++
          }
        }
      }
      END { print count + 0 }'
}

run_one_chromosome() {
  local chromosome="$1"
  local target_bcf="$2"
  local reference_bcf="$3"
  local harmonized_shared_records="$4"
  local chromosome_number="${chromosome#chr}"
  local map_file="${MAP_DIR}/chr${chromosome_number}.Campbell2016_CanFam3.1.average.SHAPEIT5.map"
  local common_target_bcf="${COMMON_TARGET_ROOT}/${chromosome}.dog11.ref1519.maf${COMMON_MAF_LABEL}.common.bcf"
  local common_reference_bcf="${COMMON_REFERENCE_ROOT}/${chromosome}.ref1519.maf${COMMON_MAF_LABEL}.common.bcf"
  local output_bcf="${PHASED_ROOT}/${chromosome}.dog11.ref1519.phase_common.bcf"
  local metrics_tsv="${METRICS_ROOT}/${chromosome}.phase_common_metrics.tsv"
  local runner_log="${LOG_ROOT}/${chromosome}.runner.log"
  local shapeit_log="${LOG_ROOT}/${chromosome}.shapeit5.log"
  local state_tsv="${STATE_DIR}/${chromosome}.status.tsv"
  local chromosome_tmp="${TMP_ROOT}/${chromosome}"
  local tmp_common_target="${chromosome_tmp}/${chromosome}.dog11.common.tmp.bcf"
  local tmp_common_reference="${chromosome_tmp}/${chromosome}.ref1519.common.tmp.bcf"
  local tmp_output="${chromosome_tmp}/${chromosome}.phase_common.tmp.bcf"
  local tmp_shapeit_log="${chromosome_tmp}/${chromosome}.shapeit5.tmp.log"
  local tmp_metrics="${chromosome_tmp}/${chromosome}.metrics.tmp.tsv"
  local common_target_records=0 common_reference_records=0
  local output_records=0 output_samples=0 unphased_hets=-1
  local observed_unphased_hets=-1
  local status="pending"
  local -a final_artifacts=(
    "$common_target_bcf" "${common_target_bcf}.csi"
    "$common_reference_bcf" "${common_reference_bcf}.csi"
    "$output_bcf" "${output_bcf}.csi" "$metrics_tsv" "$shapeit_log"
  )

  mkdir -p "$chromosome_tmp"
  : > "$runner_log"
  {
    printf '[%(%F %T)T ] start chromosome: %s\n' -1 "$chromosome"
    printf '[%(%F %T)T ] dog11 harmonize input: %s\n' -1 "$target_bcf"
    printf '[%(%F %T)T ] 1,519 reference harmonize input: %s\n' -1 "$reference_bcf"
    printf '[%(%F %T)T ] reference panel common SNP threshold: MAF>=%s\n' -1 "$COMMON_MAF"
    printf '[%(%F %T)T ] genetic map: %s\n' -1 "$map_file"
  } >> "$runner_log"

  local complete_artifacts=0 existing_artifacts=0 artifact
  for artifact in "${final_artifacts[@]}"; do
    [[ -e "$artifact" ]] && existing_artifacts=$((existing_artifacts + 1))
    [[ -s "$artifact" ]] && complete_artifacts=$((complete_artifacts + 1))
  done

  if [[ "$complete_artifacts" -eq "${#final_artifacts[@]}" ]]; then
    IFS=$'\t' read -r _ metrics_shared metrics_maf common_target_records common_reference_records \
      output_records output_samples unphased_hets _ _ _ _ < <(tail -n 1 "$metrics_tsv")
    observed_unphased_hets="$(count_unphased_heterozygotes "$output_bcf")"
    if [[ "$metrics_shared" -eq "$harmonized_shared_records" && "$metrics_maf" == "$COMMON_MAF" \
          && "$common_target_records" -gt 0 && "$common_target_records" -eq "$common_reference_records" \
          && "$common_target_records" -eq "$(count_indexed_records "$common_target_bcf")" \
          && "$common_reference_records" -eq "$(count_indexed_records "$common_reference_bcf")" \
          && "$output_records" -eq "$common_target_records" \
          && "$output_records" -eq "$(count_indexed_records "$output_bcf")" \
          && "$output_samples" -eq "$EXPECTED_TARGET_SAMPLES" \
          && "$output_samples" -eq "$($BCFTOOLS_BIN query -l "$output_bcf" | wc -l)" \
          && "$observed_unphased_hets" -eq 0 ]] \
      && cmp -s <("$BCFTOOLS_BIN" query -l "$target_bcf") <("$BCFTOOLS_BIN" query -l "$common_target_bcf") \
      && cmp -s <("$BCFTOOLS_BIN" query -l "$target_bcf") <("$BCFTOOLS_BIN" query -l "$output_bcf"); then
      unphased_hets="$observed_unphased_hets"
      status="already_done"
    else
      status="failed_existing_output_contract"
    fi
  elif [[ "$existing_artifacts" -gt 0 ]]; then
    status="failed_incomplete_existing_output"
  else
    rm -f "$tmp_common_target" "${tmp_common_target}.csi" \
      "$tmp_common_reference" "${tmp_common_reference}.csi" \
      "$tmp_output" "${tmp_output}.csi" "$tmp_shapeit_log" "$tmp_metrics"

    printf '[%(%F %T)T ] Stage 1/5: select common SNPs with MAF>=%s using AC/AN from the 1,519-dog reference panel.\n' -1 "$COMMON_MAF" >> "$runner_log"
    if "$BCFTOOLS_BIN" view --threads "$BCFTOOLS_THREADS" \
        -i "INFO/AN>0 && INFO/AC/INFO/AN>=${COMMON_MAF} && (INFO/AN-INFO/AC)/INFO/AN>=${COMMON_MAF}" \
        -Ob -o "$tmp_common_reference" "$reference_bcf" >> "$runner_log" 2>&1 \
      && "$BCFTOOLS_BIN" index --threads "$BCFTOOLS_THREADS" -f "$tmp_common_reference" >> "$runner_log" 2>&1; then
      common_reference_records="$(count_indexed_records "$tmp_common_reference")"
      if [[ "$common_reference_records" -le 0 || "$common_reference_records" -gt "$harmonized_shared_records" ]]; then
        status="failed_reference_common_record_count"
      else
        printf '[%(%F %T)T ] stage 2/5: from dog11 harmonize BCF extract the identical common SNP.\n' -1 >> "$runner_log"
        if "$BCFTOOLS_BIN" isec --threads "$BCFTOOLS_THREADS" -n=2 -w1 -Ob \
            -o "$tmp_common_target" "$target_bcf" "$tmp_common_reference" >> "$runner_log" 2>&1 \
          && "$BCFTOOLS_BIN" index --threads "$BCFTOOLS_THREADS" -f "$tmp_common_target" >> "$runner_log" 2>&1; then
          common_target_records="$(count_indexed_records "$tmp_common_target")"
          if [[ "$common_target_records" -ne "$common_reference_records" ]]; then
            status="failed_common_target_reference_mismatch"
          elif ! cmp -s <("$BCFTOOLS_BIN" query -l "$target_bcf") <("$BCFTOOLS_BIN" query -l "$tmp_common_target"); then
            status="failed_common_target_sample_order"
          else
            printf '[%(%F %T)T ] stage 3/5: run SHAPEIT5 phase_common create common SNP phase skeleton.\n' -1 >> "$runner_log"
            if "$SHAPEIT5_BIN" \
                --input "$tmp_common_target" \
                --reference "$tmp_common_reference" \
                --map "$map_file" \
                --region "$chromosome" \
                --filter-snp \
                --seed "$SEED" \
                --thread "$THREADS" \
                --output "$tmp_output" \
                --output-format bcf \
                --log "$tmp_shapeit_log" >> "$runner_log" 2>&1; then
              printf '[%(%F %T)T ] stage 4/5: create CSI and check sample, record count and sample order.\n' -1 >> "$runner_log"
              if "$BCFTOOLS_BIN" index --threads "$BCFTOOLS_THREADS" -f "$tmp_output" >> "$runner_log" 2>&1; then
                output_records="$(count_indexed_records "$tmp_output")"
                output_samples="$($BCFTOOLS_BIN query -l "$tmp_output" | wc -l)"
                if [[ "$output_samples" -ne "$EXPECTED_TARGET_SAMPLES" ]]; then
                  status="failed_output_sample_count"
                elif [[ "$output_records" -ne "$common_target_records" ]]; then
                  status="failed_output_record_count"
                elif ! cmp -s <("$BCFTOOLS_BIN" query -l "$target_bcf") <("$BCFTOOLS_BIN" query -l "$tmp_output"); then
                  status="failed_output_sample_order"
                else
                  printf '[%(%F %T)T ] stage 5/5: scan the output for heterozygous genotypes that remain unphased.\n' -1 >> "$runner_log"
                  unphased_hets="$(count_unphased_heterozygotes "$tmp_output")"
                  if [[ "$unphased_hets" -ne 0 ]]; then
                    status="failed_unphased_heterozygotes"
                  else
                    {
                      printf 'chromosome\tharmonized_shared_records\tcommon_maf\tcommon_target_records\tcommon_reference_records\toutput_records\toutput_samples\tunphased_heterozygotes\tcommon_target_bcf\tcommon_reference_bcf\toutput_bcf\tshapeit_log\n'
                      printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
                        "$chromosome" "$harmonized_shared_records" "$COMMON_MAF" \
                        "$common_target_records" "$common_reference_records" "$output_records" \
                        "$output_samples" "$unphased_hets" "$common_target_bcf" \
                        "$common_reference_bcf" "$output_bcf" "$shapeit_log"
                    } > "$tmp_metrics"
                    mv "$tmp_common_target" "$common_target_bcf"
                    mv "${tmp_common_target}.csi" "${common_target_bcf}.csi"
                    mv "$tmp_common_reference" "$common_reference_bcf"
                    mv "${tmp_common_reference}.csi" "${common_reference_bcf}.csi"
                    mv "$tmp_output" "$output_bcf"
                    mv "${tmp_output}.csi" "${output_bcf}.csi"
                    mv "$tmp_shapeit_log" "$shapeit_log"
                    mv "$tmp_metrics" "$metrics_tsv"
                    status="run_success"
                  fi
                fi
              else
                status="failed_output_index"
              fi
            else
              status="failed_shapeit5_exit"
            fi
          fi
        else
          status="failed_target_common_extraction"
        fi
      fi
    else
      status="failed_reference_common_extraction"
    fi
  fi

  {
    printf '[%(%F %T)T ] status=%s; harmonized_shared=%s; common_target=%s; common_reference=%s; output=%s; unphased_hets=%s\n' \
      -1 "$status" "$harmonized_shared_records" "$common_target_records" \
      "$common_reference_records" "$output_records" "$unphased_hets"
    printf '[%(%F %T)T ] single-chromosome task completed / done\n' -1
  } >> "$runner_log"
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$chromosome" "$target_bcf" "$reference_bcf" "$map_file" "$harmonized_shared_records" \
    "$COMMON_MAF" "$common_target_bcf" "$common_reference_bcf" "$common_target_records" \
    "$common_reference_records" "$output_bcf" "$output_records" "$unphased_hets" "$shapeit_log" "$status" > "$state_tsv"
}

if [[ $# -ne 0 ]]; then
  printf 'ERROR: this script does not accept arguments.\n' >&2
  exit 1
fi

mkdir -p "$OUT_DIR" "$COMMON_TARGET_ROOT" "$COMMON_REFERENCE_ROOT" "$PHASED_ROOT" \
  "$METRICS_ROOT" "$LOG_ROOT" "$TMP_ROOT" "$STATE_DIR"
: > "$LOG_FILE"

exec 9> "$LOCK_FILE"
if ! /usr/bin/flock -n 9; then
  fail "existing another items 1.8 task holds run lock: ${LOCK_FILE}"
fi

log_stage "Starting SHAPEIT5 phase_common for common dog11 SNPs using the 1,519-dog reference panel."
log_stage "script path: ${SCRIPT_PATH}"
log_stage "No. 1.7 step success table: ${HARMONIZATION_SUCCESS}"
log_stage "genetic map directory: ${MAP_DIR}"
log_stage "output directory: ${OUT_DIR}"
log_stage "common SNP definition=1,519 reference panel MAF>=${COMMON_MAF}"
log_stage "maximum parallel tasks=${PARALLEL_CHROMOSOMES}; per-task SHAPEIT5 threads=${THREADS}; bcftools additional threads=${BCFTOOLS_THREADS}; seed=${SEED}"
log_stage "this step does not read WhatsHap result; phase_common only create common SNP skeleton, rare variant leave for downstream phase_rare."

require_file "$HARMONIZATION_SUCCESS" "No. 1.7 step harmonize success table"
[[ -x "$SHAPEIT5_BIN" ]] || fail "SHAPEIT5_phase_common does not exist or is not executable: ${SHAPEIT5_BIN}"
[[ -x "$BCFTOOLS_BIN" ]] || fail "bcftools does not exist or is not executable: ${BCFTOOLS_BIN}"

expected_header=$'chromosome\ttarget_all_unique_bcf\ttarget_shared_exact_bcf\treference_shared_exact_bcf\texact_shared_records\tstatus'
read -r observed_header < "$HARMONIZATION_SUCCESS"
[[ "$observed_header" == "$expected_header" ]] || fail "No. 1.7 step success table has an unexpected header."

declare -A TARGET_BY_CHROM=()
declare -A REFERENCE_BY_CHROM=()
declare -A RECORDS_BY_CHROM=()
input_count=0
while IFS=$'\t' read -r chromosome _ target_bcf reference_bcf exact_shared input_status; do
  [[ "$input_status" == "run_success" || "$input_status" == "already_done" ]] || \
    fail "No. 1.7 step success table contains a non-success status: ${chromosome}/${input_status}"
  [[ "$chromosome" =~ ^chr([1-9]|[12][0-9]|3[0-8])$ ]] || fail "chromosome name does not match chr1-chr38: ${chromosome}"
  [[ -z "${TARGET_BY_CHROM[$chromosome]+x}" ]] || fail "No. 1.7 step success table contains duplicate chromosomes: ${chromosome}"
  [[ "$exact_shared" =~ ^[1-9][0-9]*$ ]] || fail "${chromosome} exact total has record count invalid: ${exact_shared}"
  require_file "$target_bcf" "${chromosome} dog11 total has SNP BCF"
  require_file "${target_bcf}.csi" "${chromosome} dog11 total has SNP CSI"
  require_file "$reference_bcf" "${chromosome} 1,519 reference total has SNP BCF"
  require_file "${reference_bcf}.csi" "${chromosome} 1,519 reference total has SNP CSI"
  TARGET_BY_CHROM["$chromosome"]="$target_bcf"
  REFERENCE_BY_CHROM["$chromosome"]="$reference_bcf"
  RECORDS_BY_CHROM["$chromosome"]="$exact_shared"
  input_count=$((input_count + 1))
done < <(tail -n +2 "$HARMONIZATION_SUCCESS")
[[ "$input_count" -eq "$EXPECTED_CHROMOSOMES" ]] || fail "Unexpected successful-chromosome count in Step 1.7: ${input_count}"

mapfile -t CHROMOSOMES < <(printf 'chr%s\n' {1..38})
for chromosome in "${CHROMOSOMES[@]}"; do
  [[ -n "${TARGET_BY_CHROM[$chromosome]+x}" ]] || fail "No. 1.7 step success table missing chromosome: ${chromosome}"
  target_bcf="${TARGET_BY_CHROM[$chromosome]}"
  reference_bcf="${REFERENCE_BY_CHROM[$chromosome]}"
  expected_records="${RECORDS_BY_CHROM[$chromosome]}"
  chromosome_number="${chromosome#chr}"
  map_file="${MAP_DIR}/chr${chromosome_number}.Campbell2016_CanFam3.1.average.SHAPEIT5.map"
  require_file "$map_file" "${chromosome} SHAPEIT5 genetic map"
  validate_map "$map_file" || fail "${chromosome} genetic map format, numeric value or order does not satisfy SHAPEIT5 contract: ${map_file}"

  target_records="$(count_indexed_records "$target_bcf")"
  reference_records="$(count_indexed_records "$reference_bcf")"
  target_samples="$($BCFTOOLS_BIN query -l "$target_bcf" | wc -l)"
  reference_samples="$($BCFTOOLS_BIN query -l "$reference_bcf" | wc -l)"
  sample_overlap="$(comm -12 <("$BCFTOOLS_BIN" query -l "$target_bcf" | sort) <("$BCFTOOLS_BIN" query -l "$reference_bcf" | sort) | wc -l)"
  [[ "$target_records" -eq "$expected_records" ]] || fail "${chromosome} dog11 record count differs from the 1.7 success table."
  [[ "$reference_records" -eq "$expected_records" ]] || fail "${chromosome} reference panel record count differs from dog11."
  [[ "$target_samples" -eq "$EXPECTED_TARGET_SAMPLES" ]] || fail "${chromosome} dog11 unexpected sample count: ${target_samples}"
  [[ "$reference_samples" -eq "$EXPECTED_REFERENCE_SAMPLES" ]] || fail "${chromosome} unexpected reference-panel sample count: ${reference_samples}"
  [[ "$sample_overlap" -eq 0 ]] || fail "${chromosome} dog11 and the reference panel share sample names: ${sample_overlap}"
  has_info_tag "$reference_bcf" "AC" || fail "${chromosome} reference panel BCF missing INFO/AC, cannot filter by reference-panel MAF."
  has_info_tag "$reference_bcf" "AN" || fail "${chromosome} reference panel BCF missing INFO/AN, cannot filter by reference-panel MAF."
done

rm -f "${STATE_DIR}"/*.status.tsv
submitted=0
for chromosome in "${CHROMOSOMES[@]}"; do
  wait_for_slot
  submitted=$((submitted + 1))
  log_stage "submit chromosome ${submitted}/${EXPECTED_CHROMOSOMES}: ${chromosome}"
  run_one_chromosome "$chromosome" "${TARGET_BY_CHROM[$chromosome]}" \
    "${REFERENCE_BY_CHROM[$chromosome]}" "${RECORDS_BY_CHROM[$chromosome]}" &
done

log_stage "all chromosome jobs submitted, waiting for background tasks to finish."
wait || true
log_stage "all chromosome tasks have returned, starting the summary."

header=$'chromosome\ttarget_shared_exact_bcf\treference_shared_exact_bcf\tgenetic_map\tharmonized_shared_records\tcommon_maf\tcommon_target_bcf\tcommon_reference_bcf\tcommon_target_records\tcommon_reference_records\toutput_bcf\toutput_records\tunphased_heterozygotes\tshapeit_log\tstatus'
printf '%s\n' "$header" > "$STATUS_TSV"
printf 'chromosome\tcommon_target_bcf\tcommon_reference_bcf\tcommon_records\toutput_bcf\toutput_records\tunphased_heterozygotes\tstatus\n' > "$SUCCESS_TSV"
printf 'chromosome\tstatus\n' > "$FAILED_TSV"

run_success=0
already_done=0
failed=0
total_harmonized_shared_records=0
total_common_records=0
total_output_records=0
for chromosome in "${CHROMOSOMES[@]}"; do
  state_tsv="${STATE_DIR}/${chromosome}.status.tsv"
  if [[ ! -s "$state_tsv" ]]; then
    failed=$((failed + 1))
    printf '%s\tfailed_missing_state\n' "$chromosome" >> "$FAILED_TSV"
    continue
  fi
  state_line="$(<"$state_tsv")"
  printf '%s\n' "$state_line" >> "$STATUS_TSV"
  IFS=$'\t' read -r state_chr _ _ _ harmonized_shared _ common_target_bcf common_reference_bcf \
    common_target_records common_reference_records output_bcf output_records unphased_hets _ status <<< "$state_line"
  case "$status" in
    run_success|already_done)
      if [[ "$status" == "run_success" ]]; then
        run_success=$((run_success + 1))
      else
        already_done=$((already_done + 1))
      fi
      total_harmonized_shared_records=$((total_harmonized_shared_records + harmonized_shared))
      total_common_records=$((total_common_records + common_target_records))
      total_output_records=$((total_output_records + output_records))
      printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$state_chr" "$common_target_bcf" "$common_reference_bcf" "$common_target_records" \
        "$output_bcf" "$output_records" "$unphased_hets" "$status" >> "$SUCCESS_TSV"
      ;;
    *)
      failed=$((failed + 1))
      printf '%s\t%s\n' "$state_chr" "$status" >> "$FAILED_TSV"
      ;;
  esac
done

{
  printf 'metric\tvalue\n'
  printf 'expected_chromosomes\t%s\n' "$EXPECTED_CHROMOSOMES"
  printf 'run_success_chromosomes\t%s\n' "$run_success"
  printf 'already_done_chromosomes\t%s\n' "$already_done"
  printf 'failed_chromosomes\t%s\n' "$failed"
  printf 'reference_panel_common_maf\t%s\n' "$COMMON_MAF"
  printf 'total_harmonized_shared_records\t%s\n' "$total_harmonized_shared_records"
  printf 'total_common_records\t%s\n' "$total_common_records"
  printf 'total_output_records\t%s\n' "$total_output_records"
  printf 'parallel_chromosomes\t%s\n' "$PARALLEL_CHROMOSOMES"
  printf 'threads_per_chromosome\t%s\n' "$THREADS"
  printf 'bcftools_extra_threads\t%s\n' "$BCFTOOLS_THREADS"
  printf 'seed\t%s\n' "$SEED"
} > "$RUN_SUMMARY_TSV"

if [[ "$failed" -ne 0 || $((run_success + already_done)) -ne "$EXPECTED_CHROMOSOMES" ]]; then
  fail "phase_common did not succeed for all tasks: run_success=${run_success}; already_done=${already_done}; failed=${failed}"
fi

log_stage "all ${EXPECTED_CHROMOSOMES} records autosome phase_common common SNP skeleton completed."
log_stage "No. 1.7 step harmonize total has SNP total count=${total_harmonized_shared_records}"
log_stage "reference panel MAF>=${COMMON_MAF} common SNP and phasing output total count=${total_output_records}"
log_stage "1.8_phase_common_dog11_with_1519_panel done"
