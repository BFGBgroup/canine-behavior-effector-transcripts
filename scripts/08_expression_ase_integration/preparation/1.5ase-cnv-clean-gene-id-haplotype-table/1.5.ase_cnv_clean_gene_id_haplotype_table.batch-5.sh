#!/usr/bin/bash
set -euo pipefail


PARALLEL=19
THREADS_PER_TASK=1
EXPECTED_TASKS=19
EXPECTED_CLEAN_MODE="known_cnv_removed"
EXPECTED_CLEAN_POLICY="exclude_only_single_state_consensus_cnv_at_7_of_9_scales"
EXPECTED_CONSENSUS_FRACTION="0.75"
EXPECTED_PHASE_SET_SCOPE="chromosome"
EXPECTED_NOT_EVALUABLE_POLICY="excluded_in_strict;retained_with_flag_in_known_cnv_removed"
EXPECTED_UNANNOTATED_CONTRACT="GeneID/GeneID_base/GeneSymbol/AllGeneIDs=(NA)<Variant> when GeneHitN=0"

TASK_TSV="/path/to/project/legacy/five-dogs/data/statistics/five_dogs_fastq_pair_stats.tsv"
INPUT_14_RUN_ID="20260821_154914"
INPUT_ROOT="/path/to/project/analysis_20260817/5-dogs/10.finemaping/preparation/1.4ase-cnv-qc-filter/${INPUT_14_RUN_ID}"
INPUT_COMPLETION_MARKER="${INPUT_ROOT}/execution_completion.marker"
SCRIPT_DIR="/path/to/repository/scripts/08_expression_ase_integration/preparation/1.5ase-cnv-clean-gene-id-haplotype-table"
SINGLE_SCRIPT="${SCRIPT_DIR}/1.5.ase_cnv_clean_gene_id_haplotype_table.single-5.py"
OUTPUT_ROOT="/path/to/project/analysis_20260817/5-dogs/10.finemaping/preparation/1.5ase-cnv-clean-gene-id-haplotype-table"
LOG_DIR="${SCRIPT_DIR}/log"

BASH_BIN="/usr/bin/bash"
PYTHON_BIN="/path/to/software/miniconda3/envs/dogs-harmony-latest/bin/python3"
AWK_BIN="/usr/bin/awk"
CAT_BIN="/usr/bin/cat"
DATE_BIN="/usr/bin/date"
FLOCK_BIN="/usr/bin/flock"
GREP_BIN="/usr/bin/grep"
MKDIR_BIN="/usr/bin/mkdir"
MKTEMP_BIN="/usr/bin/mktemp"
RM_BIN="/usr/bin/rm"
XARGS_BIN="/usr/bin/xargs"

RUN_ID="${RUN_ID:-$($DATE_BIN '+%Y%m%d_%H%M%S')}"
if [[ ! "$RUN_ID" =~ ^[0-9]{8}_[0-9]{6}$ ]]; then
  printf '[ERROR] RUN_ID must be meets YYYYMMDD_HHMMSS: %s\n' "$RUN_ID" >&2
  exit 1
fi
RUN_DIR="${OUTPUT_ROOT}/${RUN_ID}"
MASTER_LOG="${LOG_DIR}/1.5.ase_cnv_clean_gene_id_haplotype_table.${RUN_ID}.log"
SUCCESS_LIST="${LOG_DIR}/1.5.ase_cnv_clean_gene_id_haplotype_table.${RUN_ID}.success.tsv"
FAILED_LIST="${LOG_DIR}/1.5.ase_cnv_clean_gene_id_haplotype_table.${RUN_ID}.failed.tsv"
COMPLETION_MARKER="${RUN_DIR}/execution_completion.marker"

for tool in "$BASH_BIN" "$PYTHON_BIN" "$AWK_BIN" "$CAT_BIN" "$DATE_BIN" \
  "$FLOCK_BIN" "$GREP_BIN" "$MKDIR_BIN" "$MKTEMP_BIN" "$RM_BIN" "$XARGS_BIN"; do
  if [[ ! -x "$tool" ]]; then
    printf '[ERROR] required tool not executable: %s\n' "$tool" >&2
    exit 1
  fi
done
for required_file in "$TASK_TSV" "$SINGLE_SCRIPT" "$INPUT_COMPLETION_MARKER"; do
  if [[ ! -s "$required_file" ]]; then
    printf '[ERROR] required file missing or is empty: %s\n' "$required_file" >&2
    exit 1
  fi
done
if [[ ! -d "$INPUT_ROOT" ]]; then
  printf '[ERROR] 1.4 input directory does not exist: %s\n' "$INPUT_ROOT" >&2
  exit 1
fi
if ! "$GREP_BIN" -F -x -q 'execution_status=success' "$INPUT_COMPLETION_MARKER"; then
  printf '[ERROR] 1.4 the batch is not marked as fully successful: %s\n' "$INPUT_COMPLETION_MARKER" >&2
  exit 1
fi
if [[ -e "$RUN_DIR" ]]; then
  printf '[ERROR] timestamp result directory already exists, refuse to overwrite: %s\n' "$RUN_DIR" >&2
  exit 1
fi

mapfile -t TASKS < <(
  "$AWK_BIN" -F'\t' '
    NR == 1 {
      if ($1 != "animal_id" || $2 != "term") exit 10
      next
    }
    $1 != "" && $2 != "" {
      gsub(/\r/, "", $1); gsub(/\r/, "", $2)
      key = $1 "\t" $2
      if (!seen[key]++) print key
    }
  ' "$TASK_TSV"
)
TOTAL="${#TASKS[@]}"
if [[ "$TOTAL" -ne "$EXPECTED_TASKS" ]]; then
  printf '[ERROR] the task table should contain%s unique task, observed%s: %s\n' \
    "$EXPECTED_TASKS" "$TOTAL" "$TASK_TSV" >&2
  exit 1
fi

qc_metric_equals() {
  local qc_file="$1"
  local metric="$2"
  local expected="$3"
  "$AWK_BIN" -F'\t' -v metric="$metric" -v expected="$expected" '
    $1 == metric && $2 == expected {found=1}
    END {exit !found}
  ' "$qc_file"
}

PREFLIGHT_ERROR=0
for task in "${TASKS[@]}"; do
  IFS=$'\t' read -r animal tissue <<< "$task"
  sample_id="${animal}_${tissue}"
  clean_tsv="${INPUT_ROOT}/${animal}/${tissue}/${sample_id}.ase_haplotype.cnv_clean.tsv"
  qc_tsv="${INPUT_ROOT}/${animal}/${tissue}/${sample_id}.ase_haplotype.cnv_qc_summary.tsv"
  if [[ ! -s "$clean_tsv" ]]; then
    printf '[ERROR] 1.4 CNV-clean input missing or is empty: %s\n' "$clean_tsv" >&2
    PREFLIGHT_ERROR=1
  fi
  if [[ ! -s "$qc_tsv" ]]; then
    printf '[ERROR] 1.4 QC input missing or is empty: %s\n' "$qc_tsv" >&2
    PREFLIGHT_ERROR=1
    continue
  fi
  if ! "$AWK_BIN" -F'\t' '
      $1 == "cnv_result_status" &&
      ($2 == "OK" || $2 == "NO_STEP7_SEGMENT_RESULT" ||
       $2 == "NO_STEP7_SAMPLE_RESULT" || $2 == "INCOMPLETE_STEP7_SCALE_SET") {found=1}
      END {exit !found}
    ' "$qc_tsv"; then
    printf '[ERROR] 1.4 cnv_result_status does not belong to allowed status: %s\n' "$qc_tsv" >&2
    PREFLIGHT_ERROR=1
  fi
  for contract in \
    "clean_mode=${EXPECTED_CLEAN_MODE}" \
    "clean_policy=${EXPECTED_CLEAN_POLICY}" \
    "consensus_fraction=${EXPECTED_CONSENSUS_FRACTION}" \
    "phase_set_scope=${EXPECTED_PHASE_SET_SCOPE}" \
    "flag_not_evaluable_policy=${EXPECTED_NOT_EVALUABLE_POLICY}" \
    "unannotated_id_contract=${EXPECTED_UNANNOTATED_CONTRACT}"; do
    metric="${contract%%=*}"
    expected="${contract#*=}"
    if ! qc_metric_equals "$qc_tsv" "$metric" "$expected"; then
      printf '[ERROR] 1.4 QC contract does not match: %s expected%s: %s\n' "$metric" "$expected" "$qc_tsv" >&2
      PREFLIGHT_ERROR=1
    fi
  done
done
if [[ "$PREFLIGHT_ERROR" -ne 0 ]]; then
  printf '[ERROR] input preflight check failed, result directory was not created, no task was started\n' >&2
  exit 1
fi

"$MKDIR_BIN" -p "$LOG_DIR" "$OUTPUT_ROOT"
if [[ -e "$MASTER_LOG" || -e "$SUCCESS_LIST" || -e "$FAILED_LIST" ]]; then
  printf '[ERROR] this run RUN_ID of log or status file already exists, refuse to overwrite: %s\n' "$RUN_ID" >&2
  exit 1
fi
"$MKDIR_BIN" "$RUN_DIR"
: > "$MASTER_LOG"
: > "$SUCCESS_LIST"
: > "$FAILED_LIST"
exec >> "$MASTER_LOG" 2>&1

printf 'run_id\tanimal_id\ttissue\tstatus\texit_code\tstarted_at\tfinished_at\toutput_tsv\tqc_tsv\ttask_log\treason\n' > "$SUCCESS_LIST"
printf 'run_id\tanimal_id\ttissue\tstatus\texit_code\tstarted_at\tfinished_at\toutput_tsv\tqc_tsv\ttask_log\treason\n' > "$FAILED_LIST"

STATE_DIR="$($MKTEMP_BIN -d -t ase_gene_id_haplotype_batch_XXXXXX)"
COMPLETED_FILE="${STATE_DIR}/completed.count"
SUCCESS_COUNT_FILE="${STATE_DIR}/success.count"
FAIL_COUNT_FILE="${STATE_DIR}/fail.count"
LOCK_FILE="${STATE_DIR}/lock"
printf '0\n' > "$COMPLETED_FILE"
printf '0\n' > "$SUCCESS_COUNT_FILE"
printf '0\n' > "$FAIL_COUNT_FILE"
cleanup() { "$RM_BIN" -rf -- "$STATE_DIR"; }
trap cleanup EXIT

export OMP_NUM_THREADS="$THREADS_PER_TASK"
export OPENBLAS_NUM_THREADS="$THREADS_PER_TASK"
export MKL_NUM_THREADS="$THREADS_PER_TASK"
export NUMEXPR_NUM_THREADS="$THREADS_PER_TASK"
export VECLIB_MAXIMUM_THREADS="$THREADS_PER_TASK"

printf '[%s ] start 1.5 GeneID haplotype-count aggregation at level; run_id=%s\n' \
  "$($DATE_BIN '+%F %T')" "$RUN_ID"
printf '[batch] TASK_TSV=%s\n' "$TASK_TSV"
printf '[batch] INPUT_ROOT=%s\n' "$INPUT_ROOT"
printf '[batch] SINGLE_SCRIPT=%s\n' "$SINGLE_SCRIPT"
printf '[batch] RUN_DIR=%s\n' "$RUN_DIR"
printf '[batch] PARALLEL=%s THREADS_PER_TASK=%s TOTAL=%s\n' "$PARALLEL" "$THREADS_PER_TASK" "$TOTAL"
printf '[batch] AGGREGATION_KEY=GeneID+CHROM+SAMPLE+TISSUE\n'

run_one() {
  local animal="$1"
  local tissue="$2"
  local sample_id="${animal}_${tissue}"
  local input_tsv="${INPUT_ROOT}/${animal}/${tissue}/${sample_id}.ase_haplotype.cnv_clean.tsv"
  local out_dir="${RUN_DIR}/${animal}/${tissue}"
  local out_tsv="${out_dir}/${sample_id}.cnv_clean.gene_id_haplotype.tsv"
  local qc_tsv="${out_dir}/${sample_id}.cnv_clean.gene_id_haplotype.qc.tsv"
  local task_log_dir="${out_dir}/log"
  local task_log="${task_log_dir}/${sample_id}.ase_cnv_clean_gene_id_haplotype_table.log"
  local started_at finished_at rc

  "$MKDIR_BIN" -p "$task_log_dir"
  started_at="$($DATE_BIN '+%F %T')"
  {
    printf '[task] %s START animal_id=%s tissue=%s\n' "$started_at" "$animal" "$tissue"
    printf '[task] input_tsv=%s\n' "$input_tsv"
    printf '[task] output_tsv=%s\n' "$out_tsv"
    printf '[task] qc_tsv=%s\n' "$qc_tsv"
  } > "$task_log"

  set +e
  "$PYTHON_BIN" "$SINGLE_SCRIPT" "$animal" "$tissue" \
    --input-tsv "$input_tsv" \
    --out-dir "$out_dir" >> "$task_log" 2>&1
  rc=$?
  set -e
  finished_at="$($DATE_BIN '+%F %T')"

  (
    "$FLOCK_BIN" -x 200
    printf '%s\n' "$(( $("$CAT_BIN" "$COMPLETED_FILE") + 1 ))" > "$COMPLETED_FILE"
    if [[ "$rc" -eq 0 ]]; then
      printf '%s\n' "$(( $("$CAT_BIN" "$SUCCESS_COUNT_FILE") + 1 ))" > "$SUCCESS_COUNT_FILE"
      printf '%s\t%s\t%s\tsuccess\t0\t%s\t%s\t%s\t%s\t%s\tNA\n' \
        "$RUN_ID" "$animal" "$tissue" "$started_at" "$finished_at" \
        "$out_tsv" "$qc_tsv" "$task_log" >> "$SUCCESS_LIST"
    else
      printf '%s\n' "$(( $("$CAT_BIN" "$FAIL_COUNT_FILE") + 1 ))" > "$FAIL_COUNT_FILE"
      printf '%s\t%s\t%s\tfailed\t%s\t%s\t%s\t%s\t%s\t%s\tworker_exit_nonzero\n' \
        "$RUN_ID" "$animal" "$tissue" "$rc" "$started_at" "$finished_at" \
        "$out_tsv" "$qc_tsv" "$task_log" >> "$FAILED_LIST"
    fi
    printf '[%s ] progress completed=%s/%s success=%s failed=%s last=%s_%s rc=%s\n' \
      "$($DATE_BIN '+%F %T')" "$("$CAT_BIN" "$COMPLETED_FILE")" "$TOTAL" \
      "$("$CAT_BIN" "$SUCCESS_COUNT_FILE")" "$("$CAT_BIN" "$FAIL_COUNT_FILE")" \
      "$animal" "$tissue" "$rc"
  ) 200>"$LOCK_FILE"
  printf '[task] %s END animal_id=%s tissue=%s exit_code=%s\n' \
    "$finished_at" "$animal" "$tissue" "$rc" >> "$task_log"
  return "$rc"
}

export -f run_one
export RUN_ID INPUT_ROOT RUN_DIR PYTHON_BIN SINGLE_SCRIPT MKDIR_BIN DATE_BIN FLOCK_BIN CAT_BIN
export COMPLETED_FILE SUCCESS_COUNT_FILE FAIL_COUNT_FILE LOCK_FILE SUCCESS_LIST FAILED_LIST TOTAL

set +e
printf '%s\n' "${TASKS[@]}" | "$XARGS_BIN" -I '{}' -P "$PARALLEL" \
  "$BASH_BIN" --noprofile --norc -c 'IFS=$'"'"'\t'"'"' read -r animal tissue <<< "$1"; run_one "$animal" "$tissue"' _ '{}'
XARGS_RC=$?
set -e

COMPLETED="$("$CAT_BIN" "$COMPLETED_FILE")"
SUCCESS="$("$CAT_BIN" "$SUCCESS_COUNT_FILE")"
FAILED="$("$CAT_BIN" "$FAIL_COUNT_FILE")"
printf '[%s ] 1.5 batch processing finished; completed=%s/%s success=%s/%s failed=%s/%s xargs_rc=%s\n' \
  "$($DATE_BIN '+%F %T')" "$COMPLETED" "$TOTAL" "$SUCCESS" "$TOTAL" "$FAILED" "$TOTAL" "$XARGS_RC"
printf '[batch] MASTER_LOG=%s\n[batch] SUCCESS_LIST=%s\n[batch] FAILED_LIST=%s\n' \
  "$MASTER_LOG" "$SUCCESS_LIST" "$FAILED_LIST"

if [[ "$COMPLETED" -eq "$TOTAL" && "$SUCCESS" -eq "$TOTAL" && "$FAILED" -eq 0 && "$XARGS_RC" -eq 0 ]]; then
  {
    printf 'run_id=%s\n' "$RUN_ID"
    printf 'execution_status=success\n'
    printf 'completed_tasks=%s\n' "$SUCCESS"
    printf 'completed_at=%s\n' "$($DATE_BIN '+%F %T')"
  } > "$COMPLETION_MARKER"
  printf '[batch] completion_marker=%s\n' "$COMPLETION_MARKER"
  printf '[%s ] 1.5 GeneID haplotype-count aggregation task completed at level / done\n' "$($DATE_BIN '+%F %T')"
  printf '[%s ] reminder: execution success does not imply result completeness, format, counts and biological quality passes requirements, still requires separately check.\n' \
    "$($DATE_BIN '+%F %T')"
  printf 'execution_status=success\n'
  exit 0
fi

printf '[%s ] 1.5 GeneID haplotype-count aggregation failed at level, completion marker was not written.\n' "$($DATE_BIN '+%F %T')"
printf 'execution_status=failed\n'
exit 1
