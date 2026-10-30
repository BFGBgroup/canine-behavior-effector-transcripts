#!/bin/bash
set -euo pipefail


PARALLEL=1

TSV="/path/to/project/legacy/five-dogs/data/statistics/five_dogs_fastq_pair_stats.tsv"
SINGLE_SCRIPT="/path/to/repository/scripts/04_wgs_preprocessing/barkbase5/1.fastp/step1.fastp.single-wgs.sh"

BATCH_OUT_DIR="/path/to/project/legacy/five-dogs/scripts/2.wgs-fastq/1.fastp"
MASTER_LOG="${BATCH_OUT_DIR}/step1.fastp.batch-wgs.master.log"
SUCCESS_LIST="${BATCH_OUT_DIR}/step1.fastp.success.samples-wgs.tsv"
FAILED_LIST="${BATCH_OUT_DIR}/step1.fastp.failed.samples-wgs.tsv"

mkdir -p "$BATCH_OUT_DIR"

if [[ ! -f "$TSV" ]]; then
  echo "[ERROR] TSV not found: $TSV" >&2
  exit 1
fi

if [[ ! -f "$SINGLE_SCRIPT" ]]; then
  echo "[ERROR] single script not found: $SINGLE_SCRIPT" >&2
  exit 1
fi

mapfile -t TASKS < <(
  awk -F'\t' '
    NR>1 && $1!="" {
      gsub(/\r/, "", $1)
      print $1
    }
  ' "$TSV" | awk '!seen[$0]++'
)

TOTAL="${#TASKS[@]}"

if [[ "$TOTAL" -eq 0 ]]; then
  echo "[ERROR] No tasks found in TSV (after header)." >&2
  exit 1
fi

STATE_DIR="$(mktemp -d -t step1_fastp_batch_wgs_XXXXXX)"
DONE_FILE="${STATE_DIR}/done.count"
FAIL_FILE="${STATE_DIR}/fail.count"
LOCK_FILE="${STATE_DIR}/lock"

echo 0 > "$DONE_FILE"
echo 0 > "$FAIL_FILE"

echo -e "animal_id" > "$SUCCESS_LIST"
echo -e "animal_id" > "$FAILED_LIST"

{
  echo "[batch] $(date '+%F %T') START"
  echo "[batch] TSV=$TSV"
  echo "[batch] SINGLE_SCRIPT=$SINGLE_SCRIPT"
  echo "[batch] PARALLEL=$PARALLEL"
  echo "[batch] TOTAL=$TOTAL"
  echo "[batch] SUCCESS_LIST=$SUCCESS_LIST"
  echo "[batch] FAILED_LIST=$FAILED_LIST"
  echo "[batch] STATE_DIR=$STATE_DIR"
} | tee -a "$MASTER_LOG"

progress_print() {
  local animal="$1"
  local rc="$2"

  local done fail
  done="$(cat "$DONE_FILE")"
  fail="$(cat "$FAIL_FILE")"

  echo "[progress] done ${done}/${TOTAL} | failed ${fail}/${TOTAL} | last=${animal} rc=${rc}"
  echo "[progress] $(date '+%F %T') done ${done}/${TOTAL} | failed ${fail}/${TOTAL} | last=${animal} rc=${rc}" >> "$MASTER_LOG"
}

run_one() {
  local animal="$1"

  local logdir="${BATCH_OUT_DIR}/logs/${animal}"
  mkdir -p "$logdir"
  local logfile="${logdir}/${animal}.fastp.log"

  echo "[start] $(date '+%F %T') animal_id=$animal log=$logfile" >> "$MASTER_LOG"

  set +e
  bash "$SINGLE_SCRIPT" "$animal" >"$logfile" 2>&1
  local rc=$?
  set -e

  (
    flock -x 200
    if [[ "$rc" -eq 0 ]]; then
      echo $(( $(cat "$DONE_FILE") + 1 )) > "$DONE_FILE"
      printf "%s\n" "$animal" >> "$SUCCESS_LIST"
    else
      echo $(( $(cat "$FAIL_FILE") + 1 )) > "$FAIL_FILE"
      printf "%s\n" "$animal" >> "$FAILED_LIST"
    fi
  ) 200>"$LOCK_FILE"

  progress_print "$animal" "$rc"
  return "$((rc != 0))"
}

export -f run_one progress_print
export SINGLE_SCRIPT BATCH_OUT_DIR MASTER_LOG DONE_FILE FAIL_FILE LOCK_FILE TOTAL SUCCESS_LIST FAILED_LIST

set +e
printf "%s\n" "${TASKS[@]}" \
| xargs -n 1 -P "$PARALLEL" bash -lc 'run_one "$1"' _
XARGS_RC=$?
set -e

DONE="$(cat "$DONE_FILE")"
FAIL="$(cat "$FAIL_FILE")"

echo "[batch] $(date '+%F %T') FINISH done ${DONE}/${TOTAL} | failed ${FAIL}/${TOTAL} | xargs_rc=${XARGS_RC}" | tee -a "$MASTER_LOG"
echo "[batch] SUCCESS_LIST=$SUCCESS_LIST" | tee -a "$MASTER_LOG"
echo "[batch] FAILED_LIST=$FAILED_LIST"   | tee -a "$MASTER_LOG"

echo "==================== SUCCESS SAMPLES (done ${DONE}/${TOTAL}) ===================="
if [[ -s "$SUCCESS_LIST" ]]; then
  cat "$SUCCESS_LIST"
else
  echo "(none)"
fi

echo "==================== FAILED SAMPLES (failed ${FAIL}/${TOTAL}) ==================="
if [[ -s "$FAILED_LIST" ]]; then
  cat "$FAILED_LIST"
else
  echo "(none)"
fi

echo "[batch] STATE_DIR=$STATE_DIR (kept for debugging)" | tee -a "$MASTER_LOG"

if (( XARGS_RC != 0 || FAIL != 0 || DONE + FAIL != TOTAL )); then
  exit 1
fi
exit 0
