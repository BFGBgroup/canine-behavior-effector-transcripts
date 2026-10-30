#!/bin/bash
set -euo pipefail


SCRIPT="/path/to/repository/scripts/03_candidate_catalog/1.6.5_feelnc_transcript_annotation/1.6.5b_prepare_query_transcripts/1.6.5b_prepare_query_transcripts.py"
PYTHON_BIN="/usr/bin/python3"

ANNOTATED_GTF="/path/to/project/legacy/five-dogs/9.new-gene/1.6_gffcompare_stringtie3_merged_transcriptome-v6/gffcmp.stringtie3.annotated.gtf"
CLASSIFICATION_TSV="/path/to/project/legacy/five-dogs/9.new-gene/1.6_gffcompare_stringtie3_merged_transcriptome-v6/1.6_gffcompare_transcript_classification.tsv"
OUTPUT_ROOT="/path/to/project/legacy/five-dogs/9.new-gene/1.6.5_feelnc_transcript_annotation/1.6.5b_prepare_query_transcripts"
QUERY_DIR="${OUTPUT_ROOT}/query_transcripts"
LOG_DIR="${OUTPUT_ROOT}/logs-v1"

OUT_QUERY_GTF="${QUERY_DIR}/stringtie3_gffcompare_all_query_transcripts.gtf"
OUT_BASIC_INFO="${QUERY_DIR}/stringtie3_gffcompare_all_query_transcript_basic_info.tsv"
OUT_CLASS_SUMMARY="${QUERY_DIR}/stringtie3_gffcompare_all_query_class_code_summary.tsv"
OUT_ISSUE_EXAMPLES="${QUERY_DIR}/stringtie3_gffcompare_all_query_issue_examples.tsv"
QC_SUMMARY="${QUERY_DIR}/1.6.5b_prepare_query_transcripts.qc_summary.tsv"
MANIFEST="${QUERY_DIR}/1.6.5b_prepare_query_transcripts.manifest.tsv"

OVERWRITE=0
MAX_ISSUE_EXAMPLES=200

if [[ "$#" -ne 0 ]]; then
  echo "[ERROR] This batch does not accept positional arguments." >&2
  echo "Usage: bash $0" >&2
  exit 2
fi

if [[ "$OVERWRITE" != "0" && "$OVERWRITE" != "1" ]]; then
  echo "[ERROR] OVERWRITE must be 0 or 1, got: $OVERWRITE" >&2
  exit 2
fi

if [[ ! -f "$SCRIPT" ]]; then
  echo "[ERROR] required Python script not found: $SCRIPT" >&2
  exit 1
fi

if [[ ! -x "$PYTHON_BIN" ]]; then
  echo "[ERROR] python not executable: $PYTHON_BIN" >&2
  exit 1
fi

if [[ ! -s "$ANNOTATED_GTF" ]]; then
  echo "[ERROR] annotated GTF missing or empty: $ANNOTATED_GTF" >&2
  exit 1
fi

if [[ ! -s "$CLASSIFICATION_TSV" ]]; then
  echo "[ERROR] classification TSV missing or empty: $CLASSIFICATION_TSV" >&2
  exit 1
fi

mkdir -p "$QUERY_DIR" "$LOG_DIR"

MASTER_LOG="${LOG_DIR}/1.6.5b_prepare_query_transcripts.batch.master.log"
RUN_LOG="${LOG_DIR}/1.6.5b_prepare_query_transcripts.run.log"
SUCCESS_LIST="${LOG_DIR}/1.6.5b_prepare_query_transcripts.success.tsv"
FAILED_LIST="${LOG_DIR}/1.6.5b_prepare_query_transcripts.failed.tsv"

printf "task\toutdir\trc\trun_log\n" > "$SUCCESS_LIST"
printf "task\toutdir\trc\trun_log\n" > "$FAILED_LIST"

CMD=(
  "$PYTHON_BIN"
  "$SCRIPT"
  --annotated-gtf "$ANNOTATED_GTF"
  --classification-tsv "$CLASSIFICATION_TSV"
  --out-query-gtf "$OUT_QUERY_GTF"
  --out-basic-info "$OUT_BASIC_INFO"
  --out-class-summary "$OUT_CLASS_SUMMARY"
  --out-issue-examples "$OUT_ISSUE_EXAMPLES"
  --qc-summary "$QC_SUMMARY"
  --manifest "$MANIFEST"
  --max-issue-examples "$MAX_ISSUE_EXAMPLES"
)

if [[ "$OVERWRITE" == "1" ]]; then
  CMD+=(--overwrite)
fi

{
  echo "[batch] $(date '+%F %T') START"
  echo "[batch] SCRIPT=$SCRIPT"
  echo "[batch] PYTHON_BIN=$PYTHON_BIN"
  echo "[batch] ANNOTATED_GTF=$ANNOTATED_GTF"
  echo "[batch] CLASSIFICATION_TSV=$CLASSIFICATION_TSV"
  echo "[batch] OUTPUT_ROOT=$OUTPUT_ROOT"
  echo "[batch] QUERY_DIR=$QUERY_DIR"
  echo "[batch] OUT_QUERY_GTF=$OUT_QUERY_GTF"
  echo "[batch] OUT_BASIC_INFO=$OUT_BASIC_INFO"
  echo "[batch] OUT_CLASS_SUMMARY=$OUT_CLASS_SUMMARY"
  echo "[batch] OUT_ISSUE_EXAMPLES=$OUT_ISSUE_EXAMPLES"
  echo "[batch] QC_SUMMARY=$QC_SUMMARY"
  echo "[batch] MANIFEST=$MANIFEST"
  echo "[batch] LOG_DIR=$LOG_DIR"
  echo "[batch] OVERWRITE=$OVERWRITE"
  echo "[batch] MAX_ISSUE_EXAMPLES=$MAX_ISSUE_EXAMPLES"
  printf "[batch] CMD="
  printf " %q" "${CMD[@]}"
  printf "\n"
  echo "[batch] RUN_LOG=$RUN_LOG"
} > "$MASTER_LOG"

set +e
"${CMD[@]}" > "$RUN_LOG" 2>&1
RC=$?
set -e

if [[ "$RC" -eq 0 ]]; then
  printf "prepare_query_transcripts_1.6.5b\t%s\t%s\t%s\n" "$QUERY_DIR" "$RC" "$RUN_LOG" >> "$SUCCESS_LIST"
  echo "[batch] $(date '+%F %T') FINISH success rc=$RC" | tee -a "$MASTER_LOG"
  echo "[batch] OUT_QUERY_GTF=$OUT_QUERY_GTF" | tee -a "$MASTER_LOG"
  echo "[batch] OUT_BASIC_INFO=$OUT_BASIC_INFO" | tee -a "$MASTER_LOG"
  echo "[batch] OUT_ISSUE_EXAMPLES=$OUT_ISSUE_EXAMPLES" | tee -a "$MASTER_LOG"
  echo "[batch] QC_SUMMARY=$QC_SUMMARY" | tee -a "$MASTER_LOG"
else
  printf "prepare_query_transcripts_1.6.5b\t%s\t%s\t%s\n" "$QUERY_DIR" "$RC" "$RUN_LOG" >> "$FAILED_LIST"
  echo "[batch] $(date '+%F %T') FINISH failed rc=$RC" | tee -a "$MASTER_LOG"
  echo "[batch] RUN_LOG=$RUN_LOG" | tee -a "$MASTER_LOG"
fi

exit "$RC"
