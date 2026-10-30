#!/usr/bin/bash
set -euo pipefail


PARALLEL=19
THREADS_PER_TASK=1
CONSENSUS_FRACTION=0.75
CLEAN_MODE="known_cnv_removed"
PROGRESS_EVERY=50000

TASK_TSV="/path/to/project/legacy/five-dogs/data/statistics/five_dogs_fastq_pair_stats.tsv"
INPUT_13_RUN_ID="20260821_151317"
ASE_ROOT="/path/to/project/analysis_20260817/5-dogs/10.finemaping/preparation/1.3ase-haplotype-integration/${INPUT_13_RUN_ID}"
INPUT_13_COMPLETION_MARKER="${ASE_ROOT}/execution_completion.marker"
STEP7_ROOT="/path/to/project/analysis_20260609/5-dogs/4.cnv-rna/8.casper_dog_keep_locus_genes-v2"

PYTHON_BIN="/path/to/software/miniconda3/envs/dogs-harmony-latest/bin/python3"
BASH_BIN="/usr/bin/bash"
XARGS_BIN="/usr/bin/xargs"
FLOCK_BIN="/usr/bin/flock"
AWK_BIN="/usr/bin/awk"
GREP_BIN="/usr/bin/grep"
DATE_BIN="/usr/bin/date"
MKDIR_BIN="/usr/bin/mkdir"
MKTEMP_BIN="/usr/bin/mktemp"
RM_BIN="/usr/bin/rm"
CAT_BIN="/usr/bin/cat"

SCRIPT_DIR="/path/to/repository/scripts/08_expression_ase_integration/preparation/1.4ase-cnv-qc-filter"
SINGLE_SCRIPT="${SCRIPT_DIR}/1.4.variant_cnv_qc_filter.single.py"
OUTPUT_ROOT="/path/to/project/analysis_20260817/5-dogs/10.finemaping/preparation/1.4ase-cnv-qc-filter"
RUN_ID="${RUN_ID:-$(${DATE_BIN} '+%Y%m%d_%H%M%S')}"
if [[ ! "$RUN_ID" =~ ^[0-9]{8}_[0-9]{6}$ ]]; then
  echo "[ERROR] RUN_ID must use YYYYMMDD_HHMMSS format: $RUN_ID" >&2
  exit 1
fi
RUN_DIR="${OUTPUT_ROOT}/${RUN_ID}"

LOG_DIR="${SCRIPT_DIR}/log"
MASTER_LOG="${LOG_DIR}/1.4.ase_cnv_qc_filter.${RUN_ID}.log"
SUCCESS_LIST="${LOG_DIR}/1.4.ase_cnv_qc_filter.${RUN_ID}.success.tsv"
FAILED_LIST="${LOG_DIR}/1.4.ase_cnv_qc_filter.${RUN_ID}.failed.tsv"
COMPLETION_MARKER="${RUN_DIR}/execution_completion.marker"

UNANNOTATED_ID_CONTRACT=$'unannotated_id_contract\tGeneID/GeneID_base/GeneSymbol=(NA)<Variant> when GeneHitN=0'
UNANNOTATED_METADATA_CONTRACT=$'unannotated_metadata_contract\tNCBI_GeneID/SourceGTFGeneIDs=NA; GeneFeatureType/GeneBiotype=unannotated; SourceGTFGeneIDN=0'
GENE_ANNOTATION_CONTRACT=$'gene_annotation_contract\t3.7 known genes retain NCBI metadata; manual genes use GeneFeatureType/GeneBiotype=manual_new_gene and NCBI_GeneID=NA'

timestamp() {
  "${DATE_BIN}" '+[%F %T ]'
}

"${MKDIR_BIN}" -p "$LOG_DIR"
if [[ -e "$MASTER_LOG" || -e "$SUCCESS_LIST" || -e "$FAILED_LIST" ]]; then
  echo "$(timestamp) [ERROR] current RUN_ID of log or status file already exists: $RUN_ID" >&2
  exit 1
fi
: > "$MASTER_LOG"
exec >> "$MASTER_LOG" 2>&1
printf 'run_id\tanimal_id\tterm\tstatus\texit_code\tstarted_at\tfinished_at\tfull_tsv\tclean_tsv\tsummary_tsv\ttask_log\treason\n' > "$SUCCESS_LIST"
printf 'run_id\tanimal_id\tterm\tstatus\texit_code\tstarted_at\tfinished_at\tfull_tsv\tclean_tsv\tsummary_tsv\ttask_log\treason\n' > "$FAILED_LIST"

echo "$(timestamp) [batch] START"
echo "$(timestamp) [batch] run_id=$RUN_ID"
echo "$(timestamp) [batch] task_tsv=$TASK_TSV"
echo "$(timestamp) [batch] ase_root=$ASE_ROOT"
echo "$(timestamp) [batch] step7_root=$STEP7_ROOT"
echo "$(timestamp) [batch] single_script=$SINGLE_SCRIPT"
echo "$(timestamp) [batch] run_dir=$RUN_DIR"
echo "$(timestamp) [batch] consensus_fraction=$CONSENSUS_FRACTION"
echo "$(timestamp) [batch] clean_mode=$CLEAN_MODE"
echo "$(timestamp) [batch] concurrent_tasks=$PARALLEL"
echo "$(timestamp) [batch] threads_per_task=$THREADS_PER_TASK"

for executable in \
  "$PYTHON_BIN" "$BASH_BIN" "$XARGS_BIN" "$FLOCK_BIN" "$AWK_BIN" \
  "$GREP_BIN" "$DATE_BIN" "$MKDIR_BIN" "$MKTEMP_BIN" "$RM_BIN" "$CAT_BIN"; do
  if [[ ! -x "$executable" ]]; then
    echo "$(timestamp) [ERROR] required program not executable: $executable"
    exit 1
  fi
done

for required_file in "$TASK_TSV" "$SINGLE_SCRIPT" "$INPUT_13_COMPLETION_MARKER"; do
  if [[ ! -s "$required_file" ]]; then
    echo "$(timestamp) [ERROR] input file missing or is empty: $required_file"
    exit 1
  fi
done

if ! "$AWK_BIN" 'BEGIN{ok=0} $0=="execution_status=success"{ok=1} END{exit(ok?0:1)}' \
  "$INPUT_13_COMPLETION_MARKER"; then
  echo "$(timestamp) [ERROR] 1.3 completion marker does not contain execution_status=success: $INPUT_13_COMPLETION_MARKER"
  exit 1
fi

if [[ ! -d "$ASE_ROOT" || ! -d "$STEP7_ROOT" ]]; then
  echo "$(timestamp) [ERROR] 1.3 ASE or Step7 input directory does not exist"
  exit 1
fi
if [[ -e "$RUN_DIR" ]]; then
  echo "$(timestamp) [ERROR] refuse to overwrite existing run directory: $RUN_DIR"
  exit 1
fi

mapfile -t TASKS < <(
  "$AWK_BIN" -F'\t' '
    NR == 1 {
      if ($1 != "animal_id" || $2 != "term" || $3 != "raw_data_name") exit 41
      next
    }
    $1 != "" && $2 != "" {
      gsub(/\r/, "", $1)
      gsub(/\r/, "", $2)
      key = $1 FS $2
      if (!seen[key]++) print $1 FS $2
    }
  ' "$TASK_TSV"
)

TOTAL="${#TASKS[@]}"
if [[ "$TOTAL" -ne 19 ]]; then
  echo "$(timestamp) [ERROR] the task table should contain 19 unique animal_id+term, observed ${TOTAL}"
  exit 1
fi

PREFLIGHT_ERROR=0
for term in cerebellum cortex_frontal cortex_occipital pituitary_gland; do
  segment_tsv="${STEP7_ROOT}/${term}/segments.tsv"
  if [[ ! -s "$segment_tsv" ]]; then
    echo "$(timestamp) [WARN] Step7 segments missing/is empty, corresponding sample will retain and flag flag_not_evaluable: $segment_tsv"
  fi
done

for task in "${TASKS[@]}"; do
  animal="${task%%$'\t'*}"
  term="${task#*$'\t'}"
  sample_id="${animal}_${term}"
  ase_tsv="${ASE_ROOT}/${animal}/${term}/${sample_id}.ase_haplotype.tsv"
  ase_qc="${ASE_ROOT}/${animal}/${term}/${sample_id}.ase_haplotype.qc.tsv"

  if [[ ! -s "$ase_tsv" ]]; then
    echo "$(timestamp) [ERROR] 1.3 ASE input missing/is empty: $ase_tsv"
    PREFLIGHT_ERROR=1
  fi
  if [[ ! -s "$ase_qc" ]]; then
    echo "$(timestamp) [ERROR] 1.3 QC input missing/is empty: $ase_qc"
    PREFLIGHT_ERROR=1
  elif ! "$GREP_BIN" -Fxq "$UNANNOTATED_ID_CONTRACT" "$ase_qc"; then
    echo "$(timestamp) [ERROR] 1.3 QC missing unannotated synthetic ID contract: $ase_qc"
    PREFLIGHT_ERROR=1
  elif ! "$GREP_BIN" -Fxq "$UNANNOTATED_METADATA_CONTRACT" "$ase_qc"; then
    echo "$(timestamp) [ERROR] 1.3 QC missing the unannotated metadata contract: $ase_qc"
    PREFLIGHT_ERROR=1
  elif ! "$GREP_BIN" -Fxq "$GENE_ANNOTATION_CONTRACT" "$ase_qc"; then
    echo "$(timestamp) [ERROR] 1.3 QC missing known/manual gene annotation contract: $ase_qc"
    PREFLIGHT_ERROR=1
  fi
done

if [[ "$PREFLIGHT_ERROR" -ne 0 ]]; then
  echo "$(timestamp) [ERROR] run before input verify failed, no task was started"
  exit 1
fi

"$MKDIR_BIN" -p "$OUTPUT_ROOT"
if ! "$MKDIR_BIN" "$RUN_DIR"; then
  echo "$(timestamp) [ERROR] cannot create unique run directory: $RUN_DIR"
  exit 1
fi

STATE_DIR="$("$MKTEMP_BIN" -d -t ase_cnv_qc_batch_XXXXXX)"
trap '"$RM_BIN" -rf "$STATE_DIR"' EXIT
COMPLETED_FILE="${STATE_DIR}/completed.count"
SUCCESS_COUNT_FILE="${STATE_DIR}/success.count"
FAIL_COUNT_FILE="${STATE_DIR}/fail.count"
LOCK_FILE="${STATE_DIR}/lock"

printf '0\n' > "$COMPLETED_FILE"
printf '0\n' > "$SUCCESS_COUNT_FILE"
printf '0\n' > "$FAIL_COUNT_FILE"


export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

run_one() {
  local animal="$1"
  local term="$2"
  local sample_id="${animal}_${term}"
  local outdir="${RUN_DIR}/${animal}/${term}"
  local full_tsv="${outdir}/${sample_id}.ase_haplotype.cnv_annotated.tsv"
  local clean_tsv="${outdir}/${sample_id}.ase_haplotype.cnv_clean.tsv"
  local summary_tsv="${outdir}/${sample_id}.ase_haplotype.cnv_qc_summary.tsv"
  local task_log="${outdir}/log/${sample_id}.variant_cnv_qc_filter.log"
  local started_at finished_at rc status reason

  "$MKDIR_BIN" -p "${outdir}/log"
  started_at="$("$DATE_BIN" '+%F %T')"

  (
    "$FLOCK_BIN" -x 200
    echo "$(timestamp) [start] animal_id=$animal term=$term task_log=$task_log"
  ) 200>"$LOCK_FILE"

  set +e
  "$PYTHON_BIN" "$SINGLE_SCRIPT" "$animal" "$term" \
    --ase-tsv "${ASE_ROOT}/${animal}/${term}/${sample_id}.ase_haplotype.tsv" \
    --out-dir "$outdir" \
    --step7-root "$STEP7_ROOT" \
    --consensus-fraction "$CONSENSUS_FRACTION" \
    --clean-mode "$CLEAN_MODE" \
    --progress-every "$PROGRESS_EVERY" > "$task_log" 2>&1
  rc=$?
  set -e

  finished_at="$("$DATE_BIN" '+%F %T')"
  if [[ "$rc" -eq 0 ]]; then
    status="success"
    reason="worker_exit_0"
  else
    status="failed"
    reason="worker_nonzero_exit"
  fi

  (
    "$FLOCK_BIN" -x 200
    printf '%s\n' "$(( $("$CAT_BIN" "$COMPLETED_FILE") + 1 ))" > "$COMPLETED_FILE"
    if [[ "$rc" -eq 0 ]]; then
      printf '%s\n' "$(( $("$CAT_BIN" "$SUCCESS_COUNT_FILE") + 1 ))" > "$SUCCESS_COUNT_FILE"
      printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$RUN_ID" "$animal" "$term" "$status" "$rc" "$started_at" "$finished_at" \
        "$full_tsv" "$clean_tsv" "$summary_tsv" "$task_log" "$reason" >> "$SUCCESS_LIST"
    else
      printf '%s\n' "$(( $("$CAT_BIN" "$FAIL_COUNT_FILE") + 1 ))" > "$FAIL_COUNT_FILE"
      printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$RUN_ID" "$animal" "$term" "$status" "$rc" "$started_at" "$finished_at" \
        "$full_tsv" "$clean_tsv" "$summary_tsv" "$task_log" "$reason" >> "$FAILED_LIST"
    fi
    echo "$(timestamp) [progress] completed $("$CAT_BIN" "$COMPLETED_FILE")/${TOTAL} | success $("$CAT_BIN" "$SUCCESS_COUNT_FILE")/${TOTAL} | failed $("$CAT_BIN" "$FAIL_COUNT_FILE")/${TOTAL} | last=${animal} ${term} rc=${rc}"
  ) 200>"$LOCK_FILE"

  return "$rc"
}

export -f timestamp run_one
export PYTHON_BIN SINGLE_SCRIPT ASE_ROOT STEP7_ROOT RUN_DIR RUN_ID TOTAL
export DATE_BIN MKDIR_BIN FLOCK_BIN CAT_BIN
export COMPLETED_FILE SUCCESS_COUNT_FILE FAIL_COUNT_FILE LOCK_FILE SUCCESS_LIST FAILED_LIST
export CONSENSUS_FRACTION CLEAN_MODE PROGRESS_EVERY

echo "$(timestamp) [stage 1/2] preflight=PASS; dispatching ${TOTAL} tasks"
set +e
printf '%s\n' "${TASKS[@]}" \
  | "$XARGS_BIN" -I {} -P "$PARALLEL" "$BASH_BIN" --noprofile --norc -c '
      line="$1"
      animal="${line%%$'\''\t'\''*}"
      term="${line#*$'\''\t'\''}"
      run_one "$animal" "$term"
    ' _ {}
XARGS_RC=$?
set -e

COMPLETED="$("$CAT_BIN" "$COMPLETED_FILE")"
SUCCESS="$("$CAT_BIN" "$SUCCESS_COUNT_FILE")"
FAIL="$("$CAT_BIN" "$FAIL_COUNT_FILE")"

echo "$(timestamp) [stage 2/2] dispatch finished; completed=${COMPLETED}/${TOTAL}; success=${SUCCESS}/${TOTAL}; failed=${FAIL}/${TOTAL}; xargs_rc=${XARGS_RC}"
echo "$(timestamp) [batch] success_list=$SUCCESS_LIST"
echo "$(timestamp) [batch] failed_list=$FAILED_LIST"

if [[ "$COMPLETED" -ne "$TOTAL" || $((SUCCESS + FAIL)) -ne "$TOTAL" || "$FAIL" -ne 0 || "$XARGS_RC" -ne 0 ]]; then
  echo "$(timestamp) [batch] ERROR: Execution states are incomplete or a worker failed; no completion marker will be written"
  exit 1
fi

printf 'run_id=%s\nexecution_status=success\ncompleted_tasks=%s\ncompleted_at=%s\n' \
  "$RUN_ID" "$COMPLETED" "$("$DATE_BIN" '+%F %T')" > "$COMPLETION_MARKER"
echo "$(timestamp) [batch] task completed / done; execution_status=success"
echo "$(timestamp) [batch] reminder: execution completion does not establish result completeness, format, count and biological quality passes requirements, still requires separately check."
exit 0
