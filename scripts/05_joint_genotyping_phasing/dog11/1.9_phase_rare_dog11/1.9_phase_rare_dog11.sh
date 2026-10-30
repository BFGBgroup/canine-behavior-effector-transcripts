#!/usr/bin/env bash


set -euo pipefail
shopt -s nullglob

SCRIPT_PATH="/path/to/repository/scripts/05_joint_genotyping_phasing/dog11/1.9_phase_rare_dog11/1.9_phase_rare_dog11.sh"
HARMONIZATION_SUCCESS="/path/to/project/analysis_20260709/5-6dogs/1.ASE-processing/1.7_harmonize_dog11_with_1519_panel/dog11_1519_harmonization_success.tsv"
PHASE_COMMON_SUCCESS="/path/to/project/analysis_20260709/5-6dogs/1.ASE-processing/1.8_phase_common_dog11_with_1519_panel/dog11_ref1519_phase_common_success.tsv"
MAP_DIR="/path/to/project/analysis_20260709/5-6dogs/1.ASE-processing/reference/Campbell2016_CanFam3.1_genetic_map/1.convert_Campbell2016_to_SHAPEIT5/shapeit5_maps"

OUT_DIR="/path/to/project/analysis_20260709/5-6dogs/1.ASE-processing/1.9_phase_rare_dog11"
PHASED_ROOT="${OUT_DIR}/phased_common_rare_bcf"
CHUNK_ROOT="${OUT_DIR}/chunks"
CHUNK_LOG_ROOT="${OUT_DIR}/chunk_logs"
CHUNK_MANIFEST_ROOT="${OUT_DIR}/chunk_manifests"
METRICS_ROOT="${OUT_DIR}/metrics"
LOG_ROOT="${OUT_DIR}/logs"
TMP_ROOT="${OUT_DIR}/tmp"
STATE_DIR="${OUT_DIR}/state"
STATUS_TSV="${OUT_DIR}/dog11_phase_rare_status.tsv"
SUCCESS_TSV="${OUT_DIR}/dog11_phase_rare_success.tsv"
FAILED_TSV="${OUT_DIR}/dog11_phase_rare_failed.tsv"
RUN_SUMMARY_TSV="${OUT_DIR}/run_summary.tsv"
LOG_FILE="${OUT_DIR}/1.9_phase_rare_dog11.log"
LOCK_FILE="${OUT_DIR}/1.9_phase_rare_dog11.lock"

SHAPEIT5_RARE_BIN="/path/to/software/environments/shapeit5/bin/SHAPEIT5_phase_rare"
BCFTOOLS_BIN="/path/to/software/miniconda3/envs/dogs-harmony-latest/bin/bcftools"

EXPECTED_CHROMOSOMES=38
EXPECTED_SAMPLES=11
EXPECTED_TOTAL_INPUT_RECORDS=9145272
EXPECTED_TOTAL_SCAFFOLD_RECORDS=7690454
CHUNK_SIZE_BP=5000000
SCAFFOLD_BUFFER_BP=500000
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

count_region_records() {
  local bcf="$1"
  local region="$2"
  "$BCFTOOLS_BIN" view -r "$region" -H "$bcf" | wc -l
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

count_non_pass_records() {
  local bcf="$1"
  "$BCFTOOLS_BIN" view -e 'FILTER="PASS"' -H "$bcf" | wc -l
}

count_non_biallelic_snp_records() {
  local bcf="$1"
  "$BCFTOOLS_BIN" view -e 'TYPE="snp" && N_ALT=1' -H "$bcf" | wc -l
}

count_scaffold_not_in_input() {
  local scaffold_bcf="$1"
  local input_bcf="$2"
  "$BCFTOOLS_BIN" isec -C -w1 -Ou "$scaffold_bcf" "$input_bcf" \
    | "$BCFTOOLS_BIN" view -H | wc -l
}

has_info_tag() {
  local bcf="$1"
  local tag="$2"
  "$BCFTOOLS_BIN" view -h "$bcf" | grep -F "##INFO=<ID=${tag}," > /dev/null
}

get_contig_length() {
  local bcf="$1"
  local chromosome="$2"
  "$BCFTOOLS_BIN" view -h "$bcf" \
    | awk -v chromosome="$chromosome" '
        $0 ~ "^##contig=<ID=" chromosome "," {
          if (match($0, /length=[0-9]+/)) {
            print substr($0, RSTART + 7, RLENGTH - 7)
          }
        }'
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

common_genotypes_match() {
  local scaffold_bcf="$1"
  local output_bcf="$2"
  cmp -s \
    <("$BCFTOOLS_BIN" isec -n=2 -w1 -Ou "$scaffold_bcf" "$output_bcf" \
      | "$BCFTOOLS_BIN" query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n') \
    <("$BCFTOOLS_BIN" isec -n=2 -w2 -Ou "$scaffold_bcf" "$output_bcf" \
      | "$BCFTOOLS_BIN" query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n')
}

wait_for_slot() {
  while [[ "$(jobs -pr | wc -l)" -ge "$PARALLEL_CHROMOSOMES" ]]; do
    sleep 20
  done
}

run_one_chromosome() {
  local chromosome="$1"
  local input_bcf="$2"
  local scaffold_bcf="$3"
  local expected_input_records="$4"
  local expected_scaffold_records="$5"
  local chromosome_number="${chromosome#chr}"
  local map_file="${MAP_DIR}/chr${chromosome_number}.Campbell2016_CanFam3.1.average.SHAPEIT5.map"
  local output_bcf="${PHASED_ROOT}/${chromosome}.dog11.ref1519.common_rare.phased.bcf"
  local output_csi="${output_bcf}.csi"
  local metrics_tsv="${METRICS_ROOT}/${chromosome}.phase_rare_metrics.tsv"
  local runner_log="${LOG_ROOT}/${chromosome}.runner.log"
  local state_tsv="${STATE_DIR}/${chromosome}.status.tsv"
  local chunk_dir="${CHUNK_ROOT}/${chromosome}"
  local chunk_log_dir="${CHUNK_LOG_ROOT}/${chromosome}"
  local coordinates_tsv="${CHUNK_MANIFEST_ROOT}/${chromosome}.phase_rare_chunks.tsv"
  local chunk_list="${CHUNK_MANIFEST_ROOT}/${chromosome}.phase_rare_chunk_files.txt"
  local chromosome_tmp="${TMP_ROOT}/${chromosome}"
  local tmp_coordinates="${chromosome_tmp}/${chromosome}.coordinates.tmp.tsv"
  local tmp_chunk_list="${chromosome_tmp}/${chromosome}.chunk_files.tmp.txt"
  local tmp_output="${chromosome_tmp}/${chromosome}.common_rare.tmp.bcf"
  local tmp_metrics="${chromosome_tmp}/${chromosome}.metrics.tmp.tsv"
  local contig_length input_records scaffold_records output_records=0 output_samples=0
  local unphased_hets=-1 chunk_total=0 chunk_phased=0 chunk_skipped=0
  local expected_added_records=$((expected_input_records - expected_scaffold_records))
  local status="pending"

  mkdir -p "$chunk_dir" "$chunk_log_dir" "$chromosome_tmp"
  : > "$runner_log"
  {
    printf '[%(%F %T)T ] start chromosome: %s\n' -1 "$chromosome"
    printf '[%(%F %T)T ] complete dog11 input: %s\n' -1 "$input_bcf"
    printf '[%(%F %T)T ] common scaffold: %s\n' -1 "$scaffold_bcf"
    printf '[%(%F %T)T ] genetic map: %s\n' -1 "$map_file"
  } >> "$runner_log"

  if [[ -s "$output_bcf" && -s "$output_csi" && -s "$metrics_tsv" ]]; then
    IFS=$'\t' read -r _ input_records scaffold_records _ output_records output_samples \
      unphased_hets chunk_total chunk_phased chunk_skipped _ < <(tail -n 1 "$metrics_tsv")
    if [[ "$input_records" -eq "$expected_input_records" \
          && "$scaffold_records" -eq "$expected_scaffold_records" \
          && "$output_records" -eq "$expected_input_records" \
          && "$output_records" -eq "$(count_indexed_records "$output_bcf")" \
          && "$output_samples" -eq "$EXPECTED_SAMPLES" \
          && "$output_samples" -eq "$($BCFTOOLS_BIN query -l "$output_bcf" | wc -l)" \
          && "$(count_unphased_heterozygotes "$output_bcf")" -eq 0 ]] \
      && cmp -s <("$BCFTOOLS_BIN" query -l "$input_bcf") <("$BCFTOOLS_BIN" query -l "$output_bcf") \
      && common_genotypes_match "$scaffold_bcf" "$output_bcf"; then
      unphased_hets=0
      status="already_done"
    else
      status="failed_existing_output_contract"
    fi
  elif [[ -e "$output_bcf" || -e "$output_csi" || -e "$metrics_tsv" ]]; then
    status="failed_incomplete_existing_output"
  else
    rm -f "$tmp_coordinates" "$tmp_chunk_list" "$tmp_output" "${tmp_output}.csi" "$tmp_metrics"
    contig_length="$(get_contig_length "$input_bcf" "$chromosome")"
    if [[ ! "$contig_length" =~ ^[1-9][0-9]*$ ]]; then
      status="failed_contig_length"
    else
      printf '[%(%F %T)T ] stage 1/6: generate 5 Mb input and with 0.5 Mb buffered scaffold block.\n' -1 >> "$runner_log"
      awk -v chromosome="$chromosome" -v contig_length="$contig_length" \
        -v size="$CHUNK_SIZE_BP" -v buffer="$SCAFFOLD_BUFFER_BP" '
        BEGIN {
          OFS = "\t"
          print "chunk_id", "input_region", "scaffold_region"
          chunk = 1
          for (start = 1; start <= contig_length; start += size) {
            end = start + size - 1
            if (end > contig_length) end = contig_length
            scaffold_start = start - buffer
            if (scaffold_start < 1) scaffold_start = 1
            scaffold_end = end + buffer
            if (scaffold_end > contig_length) scaffold_end = contig_length
            printf "%04d\t%s:%d-%d\t%s:%d-%d\n", chunk, chromosome, start, end, chromosome, scaffold_start, scaffold_end
            chunk++
          }
        }' > "$tmp_coordinates"

      if [[ -s "$coordinates_tsv" ]]; then
        if cmp -s "$tmp_coordinates" "$coordinates_tsv"; then
          rm -f "$tmp_coordinates"
        else
          status="failed_existing_chunk_coordinates"
        fi
      else
        mv "$tmp_coordinates" "$coordinates_tsv"
      fi

      if [[ "$status" == "pending" ]]; then
        chunk_total=$(( $(wc -l < "$coordinates_tsv") - 1 ))
        if [[ "$chunk_total" -le 0 ]]; then
          status="failed_empty_chunk_plan"
        else
          printf '[%(%F %T)T ] Stage 2/6: run or reuse SHAPEIT5 phase_rare for each chunk; total chunks: %s.\n' -1 "$chunk_total" >> "$runner_log"
          : > "$tmp_chunk_list"
          while IFS=$'\t' read -r chunk_id input_region scaffold_region; do
            local chunk_bcf="${chunk_dir}/${chromosome}.phase_rare.chunk${chunk_id}.bcf"
            local chunk_csi="${chunk_bcf}.csi"
            local chunk_log="${chunk_log_dir}/${chromosome}.phase_rare.chunk${chunk_id}.log"
            local tmp_chunk_bcf="${chromosome_tmp}/${chromosome}.chunk${chunk_id}.tmp.bcf"
            local tmp_chunk_log="${chromosome_tmp}/${chromosome}.chunk${chunk_id}.tmp.log"
            local chunk_input_records scaffold_region_records chunk_output_records chunk_output_samples

            chunk_input_records="$(count_region_records "$input_bcf" "$input_region")"
            if [[ "$chunk_input_records" -eq 0 ]]; then
              chunk_skipped=$((chunk_skipped + 1))
              printf '[%(%F %T)T ] skip empty block %s/%s: %s\n' -1 "$chunk_id" "$chunk_total" "$input_region" >> "$runner_log"
              continue
            fi
            scaffold_region_records="$(count_region_records "$scaffold_bcf" "$scaffold_region")"
            if [[ "$scaffold_region_records" -eq 0 ]]; then
              status="failed_empty_scaffold_region"
              break
            fi

            if [[ -s "$chunk_bcf" && -s "$chunk_csi" && -s "$chunk_log" ]]; then
              chunk_output_records="$(count_indexed_records "$chunk_bcf")"
              chunk_output_samples="$($BCFTOOLS_BIN query -l "$chunk_bcf" | wc -l)"
              if [[ "$chunk_output_records" -ne "$chunk_input_records" \
                    || "$chunk_output_samples" -ne "$EXPECTED_SAMPLES" ]] \
                || ! cmp -s <("$BCFTOOLS_BIN" query -l "$input_bcf") <("$BCFTOOLS_BIN" query -l "$chunk_bcf"); then
                status="failed_existing_chunk_contract"
                break
              fi
              printf '[%(%F %T)T ] reuse chunk %s/%s: %s\n' -1 "$chunk_id" "$chunk_total" "$chunk_bcf" >> "$runner_log"
            elif [[ -e "$chunk_bcf" || -e "$chunk_csi" || -e "$chunk_log" ]]; then
              status="failed_incomplete_existing_chunk"
              break
            else
              rm -f "$tmp_chunk_bcf" "${tmp_chunk_bcf}.csi" "$tmp_chunk_log"
              printf '[%(%F %T)T ] run chunk %s/%s: input=%s; scaffold=%s\n' \
                -1 "$chunk_id" "$chunk_total" "$input_region" "$scaffold_region" >> "$runner_log"
              if "$SHAPEIT5_RARE_BIN" \
                  --input "$input_bcf" \
                  --input-region "$input_region" \
                  --scaffold "$scaffold_bcf" \
                  --scaffold-region "$scaffold_region" \
                  --map "$map_file" \
                  --seed "$SEED" \
                  --thread "$THREADS" \
                  --output "$tmp_chunk_bcf" \
                  --log "$tmp_chunk_log" >> "$runner_log" 2>&1 \
                && "$BCFTOOLS_BIN" index --threads "$BCFTOOLS_THREADS" -f "$tmp_chunk_bcf" >> "$runner_log" 2>&1; then
                chunk_output_records="$(count_indexed_records "$tmp_chunk_bcf")"
                chunk_output_samples="$($BCFTOOLS_BIN query -l "$tmp_chunk_bcf" | wc -l)"
                if [[ "$chunk_output_records" -ne "$chunk_input_records" ]]; then
                  status="failed_chunk_record_count"
                  break
                elif [[ "$chunk_output_samples" -ne "$EXPECTED_SAMPLES" ]]; then
                  status="failed_chunk_sample_count"
                  break
                elif ! cmp -s <("$BCFTOOLS_BIN" query -l "$input_bcf") <("$BCFTOOLS_BIN" query -l "$tmp_chunk_bcf"); then
                  status="failed_chunk_sample_order"
                  break
                else
                  mv "$tmp_chunk_bcf" "$chunk_bcf"
                  mv "${tmp_chunk_bcf}.csi" "$chunk_csi"
                  mv "$tmp_chunk_log" "$chunk_log"
                fi
              else
                status="failed_phase_rare_chunk"
                break
              fi
            fi
            printf '%s\n' "$chunk_bcf" >> "$tmp_chunk_list"
            chunk_phased=$((chunk_phased + 1))
          done < <(tail -n +2 "$coordinates_tsv")
        fi
      fi

      if [[ "$status" == "pending" ]]; then
        if [[ "$chunk_phased" -le 0 || $((chunk_phased + chunk_skipped)) -ne "$chunk_total" ]]; then
          status="failed_chunk_completion_count"
        else
          printf '[%(%F %T)T ] Stage 3/6: concatenate %s nonempty common+rare chunks without loss.\n' -1 "$chunk_phased" >> "$runner_log"
          if "$BCFTOOLS_BIN" concat --naive -Ob -o "$tmp_output" -f "$tmp_chunk_list" >> "$runner_log" 2>&1; then
            printf '[%(%F %T)T ] stage 4/6: build whole-chromosome CSI and check record count and sample order.\n' -1 >> "$runner_log"
            if "$BCFTOOLS_BIN" index --threads "$BCFTOOLS_THREADS" -f "$tmp_output" >> "$runner_log" 2>&1; then
              input_records="$(count_indexed_records "$input_bcf")"
              scaffold_records="$(count_indexed_records "$scaffold_bcf")"
              output_records="$(count_indexed_records "$tmp_output")"
              output_samples="$($BCFTOOLS_BIN query -l "$tmp_output" | wc -l)"
              if [[ "$output_records" -ne "$expected_input_records" ]]; then
                status="failed_output_record_count"
              elif [[ "$output_samples" -ne "$EXPECTED_SAMPLES" ]]; then
                status="failed_output_sample_count"
              elif ! cmp -s <("$BCFTOOLS_BIN" query -l "$input_bcf") <("$BCFTOOLS_BIN" query -l "$tmp_output"); then
                status="failed_output_sample_order"
              else
                printf '[%(%F %T)T ] stage 5/6: check common scaffold orientation, and scan unphased heterozygous GT.\n' -1 >> "$runner_log"
                if ! common_genotypes_match "$scaffold_bcf" "$tmp_output"; then
                  status="failed_common_scaffold_genotype_mismatch"
                else
                  unphased_hets="$(count_unphased_heterozygotes "$tmp_output")"
                  if [[ "$unphased_hets" -ne 0 ]]; then
                    status="failed_unphased_heterozygotes"
                  else
                    printf '[%(%F %T)T ] stage 6/6: write QC and atomic release final result.\n' -1 >> "$runner_log"
                    {
                      printf 'chromosome\tinput_records\tscaffold_records\tadded_noncommon_records\toutput_records\toutput_samples\tunphased_heterozygotes\tchunk_total\tchunk_phased\tchunk_skipped\toutput_bcf\n'
                      printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
                        "$chromosome" "$input_records" "$scaffold_records" "$expected_added_records" \
                        "$output_records" "$output_samples" "$unphased_hets" "$chunk_total" \
                        "$chunk_phased" "$chunk_skipped" "$output_bcf"
                    } > "$tmp_metrics"
                    mv "$tmp_output" "$output_bcf"
                    mv "${tmp_output}.csi" "$output_csi"
                    mv "$tmp_chunk_list" "$chunk_list"
                    mv "$tmp_metrics" "$metrics_tsv"
                    status="run_success"
                  fi
                fi
              fi
            else
              status="failed_output_index"
            fi
          else
            status="failed_chunk_concat"
          fi
        fi
      fi
    fi
  fi

  {
    printf '[%(%F %T)T ] status=%s; input=%s; scaffold=%s; added=%s; output=%s; unphased_hets=%s; chunks=%s/%s\n' \
      -1 "$status" "$expected_input_records" "$expected_scaffold_records" "$expected_added_records" \
      "$output_records" "$unphased_hets" "$chunk_phased" "$chunk_total"
    printf '[%(%F %T)T ] single-chromosome task completed / done\n' -1
  } >> "$runner_log"
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$chromosome" "$input_bcf" "$scaffold_bcf" "$map_file" "$expected_input_records" \
    "$expected_scaffold_records" "$expected_added_records" "$output_bcf" "$output_records" \
    "$output_samples" "$unphased_hets" "$chunk_total" "$chunk_phased" "$chunk_skipped" \
    "$metrics_tsv" "$status" > "$state_tsv"
}

if [[ $# -ne 0 ]]; then
  printf 'ERROR: this script does not accept arguments.\n' >&2
  exit 1
fi

mkdir -p "$OUT_DIR" "$PHASED_ROOT" "$CHUNK_ROOT" "$CHUNK_LOG_ROOT" \
  "$CHUNK_MANIFEST_ROOT" "$METRICS_ROOT" "$LOG_ROOT" "$TMP_ROOT" "$STATE_DIR"
: > "$LOG_FILE"

exec 9> "$LOCK_FILE"
if ! /usr/bin/flock -n 9; then
  fail "existing another items 1.9 task holds run lock: ${LOCK_FILE}"
fi

log_stage "start set dog11 complete SNP phasing to No. 1.8 step common scaffold."
log_stage "script path: ${SCRIPT_PATH}"
log_stage "No. 1.7 step success table: ${HARMONIZATION_SUCCESS}"
log_stage "No. 1.8 step success table: ${PHASE_COMMON_SUCCESS}"
log_stage "genetic map directory: ${MAP_DIR}"
log_stage "output directory: ${OUT_DIR}"
log_stage "block: input=${CHUNK_SIZE_BP} bp; scaffold both sides buffer=${SCAFFOLD_BUFFER_BP} bp"
log_stage "maximum parallel tasks=${PARALLEL_CHROMOSOMES}; per-task phase_rare threads=${THREADS}; bcftools additional threads=${BCFTOOLS_THREADS}; seed=${SEED}"

require_file "$HARMONIZATION_SUCCESS" "No. 1.7 step harmonize success table"
require_file "$PHASE_COMMON_SUCCESS" "No. 1.8 step phase_common success table"
[[ -x "$SHAPEIT5_RARE_BIN" ]] || fail "SHAPEIT5_phase_rare does not exist or is not executable: ${SHAPEIT5_RARE_BIN}"
[[ -x "$BCFTOOLS_BIN" ]] || fail "bcftools does not exist or is not executable: ${BCFTOOLS_BIN}"

harmonization_header=$'chromosome\ttarget_all_unique_bcf\ttarget_shared_exact_bcf\treference_shared_exact_bcf\texact_shared_records\tstatus'
read -r observed_header < "$HARMONIZATION_SUCCESS"
[[ "$observed_header" == "$harmonization_header" ]] || fail "No. 1.7 step success table has an unexpected header."

phase_common_header=$'chromosome\tcommon_target_bcf\tcommon_reference_bcf\tcommon_records\toutput_bcf\toutput_records\tunphased_heterozygotes\tstatus'
read -r observed_header < "$PHASE_COMMON_SUCCESS"
[[ "$observed_header" == "$phase_common_header" ]] || fail "No. 1.8 step success table has an unexpected header."

declare -A INPUT_BY_CHROM=()
declare -A INPUT_RECORDS_BY_CHROM=()
declare -A SCAFFOLD_BY_CHROM=()
declare -A SCAFFOLD_RECORDS_BY_CHROM=()

while IFS=$'\t' read -r chromosome input_bcf _ _ _ input_status; do
  [[ "$input_status" == "run_success" || "$input_status" == "already_done" ]] || \
    fail "No. 1.7 step success table contains a non-success status: ${chromosome}/${input_status}"
  [[ "$chromosome" =~ ^chr([1-9]|[12][0-9]|3[0-8])$ ]] || fail "chromosome name does not match chr1-chr38: ${chromosome}"
  [[ -z "${INPUT_BY_CHROM[$chromosome]+x}" ]] || fail "No. 1.7 step success table contains duplicate chromosomes: ${chromosome}"
  require_file "$input_bcf" "${chromosome} dog11 complete BCF"
  require_file "${input_bcf}.csi" "${chromosome} dog11 complete CSI"
  INPUT_BY_CHROM["$chromosome"]="$input_bcf"
  INPUT_RECORDS_BY_CHROM["$chromosome"]="$(count_indexed_records "$input_bcf")"
done < <(tail -n +2 "$HARMONIZATION_SUCCESS")

while IFS=$'\t' read -r chromosome _ _ common_records scaffold_bcf output_records unphased_hets phase_status; do
  [[ "$phase_status" == "run_success" || "$phase_status" == "already_done" ]] || \
    fail "No. 1.8 step success table contains a non-success status: ${chromosome}/${phase_status}"
  [[ "$common_records" =~ ^[1-9][0-9]*$ && "$output_records" -eq "$common_records" ]] || \
    fail "${chromosome} No. 1.8 step record count invalid or mismatch."
  [[ "$unphased_hets" -eq 0 ]] || fail "${chromosome} No. 1.8 step still contains unphased heterozygous GT: ${unphased_hets}"
  [[ -z "${SCAFFOLD_BY_CHROM[$chromosome]+x}" ]] || fail "No. 1.8 step success table contains duplicate chromosomes: ${chromosome}"
  require_file "$scaffold_bcf" "${chromosome} common scaffold BCF"
  require_file "${scaffold_bcf}.csi" "${chromosome} common scaffold CSI"
  SCAFFOLD_BY_CHROM["$chromosome"]="$scaffold_bcf"
  SCAFFOLD_RECORDS_BY_CHROM["$chromosome"]="$common_records"
done < <(tail -n +2 "$PHASE_COMMON_SUCCESS")

mapfile -t CHROMOSOMES < <(printf 'chr%s\n' {1..38})
total_input_records=0
total_scaffold_records=0
for chromosome in "${CHROMOSOMES[@]}"; do
  [[ -n "${INPUT_BY_CHROM[$chromosome]+x}" ]] || fail "No. 1.7 step success table missing chromosome: ${chromosome}"
  [[ -n "${SCAFFOLD_BY_CHROM[$chromosome]+x}" ]] || fail "No. 1.8 step success table missing chromosome: ${chromosome}"
  input_bcf="${INPUT_BY_CHROM[$chromosome]}"
  scaffold_bcf="${SCAFFOLD_BY_CHROM[$chromosome]}"
  input_records="${INPUT_RECORDS_BY_CHROM[$chromosome]}"
  scaffold_records="${SCAFFOLD_RECORDS_BY_CHROM[$chromosome]}"
  chromosome_number="${chromosome#chr}"
  map_file="${MAP_DIR}/chr${chromosome_number}.Campbell2016_CanFam3.1.average.SHAPEIT5.map"
  require_file "$map_file" "${chromosome} SHAPEIT5 genetic map"
  validate_map "$map_file" || fail "${chromosome} genetic map format, numeric value or order does not satisfy SHAPEIT5 contract."
  [[ "$(count_indexed_records "$scaffold_bcf")" -eq "$scaffold_records" ]] || \
    fail "${chromosome} scaffold record count differs from the Step 1.8 success table."
  [[ "$input_records" -gt "$scaffold_records" ]] || fail "${chromosome} complete input record count is not greater than scaffold."
  [[ "$($BCFTOOLS_BIN query -l "$input_bcf" | wc -l)" -eq "$EXPECTED_SAMPLES" ]] || \
    fail "${chromosome} complete input unexpected sample count."
  [[ "$($BCFTOOLS_BIN query -l "$scaffold_bcf" | wc -l)" -eq "$EXPECTED_SAMPLES" ]] || \
    fail "${chromosome} scaffold unexpected sample count."
  cmp -s <("$BCFTOOLS_BIN" query -l "$input_bcf") <("$BCFTOOLS_BIN" query -l "$scaffold_bcf") || \
    fail "${chromosome} complete input and scaffold differ in sample identity or order."
  [[ "$(count_non_pass_records "$input_bcf")" -eq 0 ]] || fail "${chromosome} complete input contains non-PASS record."
  [[ "$(count_non_biallelic_snp_records "$input_bcf")" -eq 0 ]] || fail "${chromosome} The complete input contains a variant other than a biallelic SNP."
  [[ "$(count_scaffold_not_in_input "$scaffold_bcf" "$input_bcf")" -eq 0 ]] || \
    fail "${chromosome} scaffold contains entries without an exact match in the complete input REF/ALT."
  has_info_tag "$input_bcf" "AC" || fail "${chromosome} complete input missing INFO/AC."
  has_info_tag "$input_bcf" "AN" || fail "${chromosome} complete input missing INFO/AN."
  total_input_records=$((total_input_records + input_records))
  total_scaffold_records=$((total_scaffold_records + scaffold_records))
done

[[ "$total_input_records" -eq "$EXPECTED_TOTAL_INPUT_RECORDS" ]] || \
  fail "complete input unexpected total record count: expected=${EXPECTED_TOTAL_INPUT_RECORDS}; observed=${total_input_records}"
[[ "$total_scaffold_records" -eq "$EXPECTED_TOTAL_SCAFFOLD_RECORDS" ]] || \
  fail "scaffold unexpected total record count: expected=${EXPECTED_TOTAL_SCAFFOLD_RECORDS}; observed=${total_scaffold_records}"

rm -f "${STATE_DIR}"/*.status.tsv
submitted=0
for chromosome in "${CHROMOSOMES[@]}"; do
  wait_for_slot
  submitted=$((submitted + 1))
  log_stage "submit chromosome ${submitted}/${EXPECTED_CHROMOSOMES}: ${chromosome}"
  run_one_chromosome "$chromosome" "${INPUT_BY_CHROM[$chromosome]}" \
    "${SCAFFOLD_BY_CHROM[$chromosome]}" "${INPUT_RECORDS_BY_CHROM[$chromosome]}" \
    "${SCAFFOLD_RECORDS_BY_CHROM[$chromosome]}" &
done

log_stage "all chromosome jobs submitted, waiting for background tasks to finish."
wait || true
log_stage "all chromosome tasks have returned, starting the summary."

header=$'chromosome\tinput_bcf\tscaffold_bcf\tgenetic_map\tinput_records\tscaffold_records\tadded_noncommon_records\toutput_bcf\toutput_records\toutput_samples\tunphased_heterozygotes\tchunk_total\tchunk_phased\tchunk_skipped\tmetrics_tsv\tstatus'
printf '%s\n' "$header" > "$STATUS_TSV"
printf 'chromosome\toutput_bcf\tinput_records\tscaffold_records\tadded_noncommon_records\toutput_records\tunphased_heterozygotes\tstatus\n' > "$SUCCESS_TSV"
printf 'chromosome\tstatus\n' > "$FAILED_TSV"

run_success=0
already_done=0
failed=0
total_output_records=0
total_unphased_hets=0
total_chunks=0
for chromosome in "${CHROMOSOMES[@]}"; do
  state_tsv="${STATE_DIR}/${chromosome}.status.tsv"
  if [[ ! -s "$state_tsv" ]]; then
    failed=$((failed + 1))
    printf '%s\tfailed_missing_state\n' "$chromosome" >> "$FAILED_TSV"
    continue
  fi
  state_line="$(<"$state_tsv")"
  printf '%s\n' "$state_line" >> "$STATUS_TSV"
  IFS=$'\t' read -r state_chr _ _ _ input_records scaffold_records added_records output_bcf \
    output_records _ unphased_hets chunk_total _ _ _ status <<< "$state_line"
  case "$status" in
    run_success|already_done)
      if [[ "$status" == "run_success" ]]; then
        run_success=$((run_success + 1))
      else
        already_done=$((already_done + 1))
      fi
      total_output_records=$((total_output_records + output_records))
      total_unphased_hets=$((total_unphased_hets + unphased_hets))
      total_chunks=$((total_chunks + chunk_total))
      printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$state_chr" "$output_bcf" "$input_records" "$scaffold_records" "$added_records" \
        "$output_records" "$unphased_hets" "$status" >> "$SUCCESS_TSV"
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
  printf 'total_input_records\t%s\n' "$total_input_records"
  printf 'total_scaffold_records\t%s\n' "$total_scaffold_records"
  printf 'total_added_noncommon_records\t%s\n' "$((total_input_records - total_scaffold_records))"
  printf 'total_output_records\t%s\n' "$total_output_records"
  printf 'total_unphased_heterozygotes\t%s\n' "$total_unphased_hets"
  printf 'total_chunks\t%s\n' "$total_chunks"
  printf 'chunk_size_bp\t%s\n' "$CHUNK_SIZE_BP"
  printf 'scaffold_buffer_bp\t%s\n' "$SCAFFOLD_BUFFER_BP"
  printf 'parallel_chromosomes\t%s\n' "$PARALLEL_CHROMOSOMES"
  printf 'threads_per_chromosome\t%s\n' "$THREADS"
  printf 'bcftools_extra_threads\t%s\n' "$BCFTOOLS_THREADS"
  printf 'seed\t%s\n' "$SEED"
} > "$RUN_SUMMARY_TSV"

if [[ "$failed" -ne 0 || $((run_success + already_done)) -ne "$EXPECTED_CHROMOSOMES" \
      || "$total_output_records" -ne "$total_input_records" || "$total_unphased_hets" -ne 0 ]]; then
  fail "phase_rare did not succeed for all tasks: run_success=${run_success}; already_done=${already_done}; failed=${failed}"
fi

log_stage "all ${EXPECTED_CHROMOSOMES} records autosome phase_rare completed."
log_stage "common scaffold site total count=${total_scaffold_records}"
log_stage "supplement non-common site total count=$((total_input_records - total_scaffold_records))"
log_stage "final common+rare phase site total count=${total_output_records}; unphased heterozygous GT=${total_unphased_hets}"
log_stage "1.9_phase_rare_dog11 done"
