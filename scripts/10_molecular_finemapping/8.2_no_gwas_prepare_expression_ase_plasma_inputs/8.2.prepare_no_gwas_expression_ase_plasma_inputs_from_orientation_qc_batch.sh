#!/usr/bin/bash
set -euo pipefail


SCRIPT_DIR="/path/to/repository/scripts/10_molecular_finemapping/8.2_no_gwas_prepare_expression_ase_plasma_inputs"
PY_SCRIPT="${SCRIPT_DIR}/8.2.prepare_no_gwas_expression_ase_plasma_inputs_from_orientation_qc.py"
PYTHON_BIN="/path/to/software/miniconda3/envs/dogs-harmony-latest/bin/python3"
DATE_BIN="/usr/bin/date"
MKDIR_BIN="/usr/bin/mkdir"

OUTPUT_ROOT="/path/to/project/analysis_20260817/5-dogs/10.finemaping/formal/8.2_no_gwas_prepare_expression_ase_plasma_inputs"
RUN_ID="${RUN_ID:-$(${DATE_BIN} '+%Y%m%d_%H%M%S')}"
if [[ ! "$RUN_ID" =~ ^[0-9]{8}_[0-9]{6}$ ]]; then
  printf '[ERROR] RUN_ID must be meets YYYYMMDD_HHMMSS format: %s\n' "$RUN_ID" >&2
  exit 1
fi
RUN_DIR="${OUTPUT_ROOT}/${RUN_ID}"

LOG_DIR="${SCRIPT_DIR}/log"
MASTER_LOG="${LOG_DIR}/8.2.prepare_no_gwas_expression_ase_plasma_inputs.${RUN_ID}.log"
SUCCESS_LIST="${LOG_DIR}/8.2.prepare_no_gwas_expression_ase_plasma_inputs.${RUN_ID}.success.tsv"
FAILED_LIST="${LOG_DIR}/8.2.prepare_no_gwas_expression_ase_plasma_inputs.${RUN_ID}.failed.tsv"
COMPLETION_MARKER="${RUN_DIR}/execution_completion.marker"

timestamp() {
  "${DATE_BIN}" '+[%F %T ]'
}

for executable in "$PYTHON_BIN" "$DATE_BIN" "$MKDIR_BIN"; do
  if [[ ! -x "$executable" ]]; then
    printf '[ERROR] required program not executable: %s\n' "$executable" >&2
    exit 1
  fi
done
if [[ ! -s "$PY_SCRIPT" ]]; then
  printf '[ERROR] Python main program does not exist or is empty: %s\n' "$PY_SCRIPT" >&2
  exit 1
fi

"${MKDIR_BIN}" -p "$LOG_DIR" "$OUTPUT_ROOT"
for path in "$RUN_DIR" "$MASTER_LOG" "$SUCCESS_LIST" "$FAILED_LIST"; do
  if [[ -e "$path" ]]; then
    printf '[ERROR] current RUN_ID corresponding path already exists, refuse to overwrite: %s\n' "$path" >&2
    exit 1
  fi
done
"${MKDIR_BIN}" "$RUN_DIR"
: > "$MASTER_LOG"
printf 'run_id\ttask_id\tstatus\texit_code\tstarted_at\tfinished_at\toutput_dir\ttask_log\treason\n' > "$SUCCESS_LIST"
printf 'run_id\ttask_id\tstatus\texit_code\tstarted_at\tfinished_at\toutput_dir\ttask_log\treason\n' > "$FAILED_LIST"
exec >> "$MASTER_LOG" 2>&1

started_at="$(${DATE_BIN} '+%F %T')"
echo "$(timestamp) [start] five dogs no-GWAS PLASMA expression+ASE joint input prepare"
echo "$(timestamp) [config] run_id=$RUN_ID"
echo "$(timestamp) [config] python=$PYTHON_BIN"
echo "$(timestamp) [config] worker=$PY_SCRIPT"
echo "$(timestamp) [config] run_dir=$RUN_DIR"

set +e
"${PYTHON_BIN}" "$PY_SCRIPT" --out-dir "$RUN_DIR"
exit_code=$?
set -e
finished_at="$(${DATE_BIN} '+%F %T')"

if [[ "$exit_code" -ne 0 ]]; then
  printf '%s\t%s\tfailed\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$RUN_ID" "prepare_joint_inputs" "$exit_code" "$started_at" "$finished_at" \
    "$RUN_DIR" "$MASTER_LOG" "Python main program returned a nonzero exit code" >> "$FAILED_LIST"
  echo "$(timestamp) [failed] Python main program exit code=$exit_code"
  echo "$(timestamp) [summary] dispatched=1 success=0 failed=1"
  echo "$(timestamp) execution_status=failed"
  exit "$exit_code"
fi

printf '%s\t%s\tsuccess\t0\t%s\t%s\t%s\t%s\t%s\n' \
  "$RUN_ID" "prepare_joint_inputs" "$started_at" "$finished_at" \
  "$RUN_DIR" "$MASTER_LOG" "" >> "$SUCCESS_LIST"
{
  printf 'run_id=%s\n' "$RUN_ID"
  printf 'execution_status=success\n'
  printf 'finished_at=%s\n' "$finished_at"
} > "$COMPLETION_MARKER"

echo "$(timestamp) task completed / done"
echo "$(timestamp) [summary] dispatched=1 success=1 failed=0"
echo "$(timestamp) execution_status=success"
echo "$(timestamp) [reminder] execution completion does not establish output completeness, format, count or biological quality passes requirements, please check separately result."
