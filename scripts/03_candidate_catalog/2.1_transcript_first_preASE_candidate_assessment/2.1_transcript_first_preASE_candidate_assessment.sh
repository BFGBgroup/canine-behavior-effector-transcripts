#!/bin/bash
set -euo pipefail


SCRIPT="/path/to/repository/scripts/03_candidate_catalog/2.1_transcript_first_preASE_candidate_assessment/2.1_transcript_first_preASE_candidate_assessment.py"
PYTHON_BIN="/usr/bin/python3"

CLASSIFICATION_TSV="/path/to/project/legacy/five-dogs/9.new-gene/1.6_gffcompare_stringtie3_merged_transcriptome-v6/1.6_gffcompare_transcript_classification.tsv"
TRANSCRIPT_EXON_BED="/path/to/project/legacy/five-dogs/9.new-gene/1.6_gffcompare_stringtie3_merged_transcriptome-v6/transcript_exon.bed"
FEELNC_STATUS_TSV="/path/to/project/legacy/five-dogs/9.new-gene/1.6.5_feelnc_transcript_annotation/1.6.5c_run_feelnc_on_all_query-v4/qc/1.6.5c_all_query_feelnc_module_status.tsv"
REFERENCE_GTF="/path/to/project/legacy/reference/GCF_000002285.3_CanFam3.1_genomic_chr.gtf"

OUTPUT_ROOT="/path/to/project/analysis_20260709/5-dogs/9.new-gene/2.1_transcript_first_preASE_candidate_assessment"
LOG_DIR="/path/to/repository/scripts/03_candidate_catalog/2.1_transcript_first_preASE_candidate_assessment/logs"
CONTAINMENT_OVERLAP_FRACTION="0.95"
OVERWRITE=0

log_msg() {
  printf '[%s ] %s\n' "$(date '+%F %T')" "$*"
}

if [[ "$#" -ne 0 ]]; then
  echo "[ERROR] this script does not accept positional arguments." >&2
  echo "Usage: bash $0" >&2
  exit 2
fi

if [[ ! -f "$SCRIPT" ]]; then
  echo "[ERROR] Python script does not exist: $SCRIPT" >&2
  exit 1
fi

if [[ ! -x "$PYTHON_BIN" ]]; then
  echo "[ERROR] Python not executable: $PYTHON_BIN" >&2
  exit 1
fi

for required_file in \
  "$CLASSIFICATION_TSV" \
  "$TRANSCRIPT_EXON_BED" \
  "$FEELNC_STATUS_TSV" \
  "$REFERENCE_GTF"; do
  if [[ ! -s "$required_file" ]]; then
    echo "[ERROR] required input does not exist or is empty: $required_file" >&2
    exit 1
  fi
done

if [[ "$OVERWRITE" != "0" && "$OVERWRITE" != "1" ]]; then
  echo "[ERROR] OVERWRITE must be 0 or 1, current: $OVERWRITE" >&2
  exit 2
fi

EXPECTED_TRANSCRIPT_COUNT="$(( $(wc -l < "$CLASSIFICATION_TSV") - 1 ))"
if [[ "$EXPECTED_TRANSCRIPT_COUNT" -le 0 ]]; then
  echo "[ERROR] classification no valid data rows: $CLASSIFICATION_TSV" >&2
  exit 1
fi

if [[ -d "$OUTPUT_ROOT" && "$OVERWRITE" == "0" ]]; then
  shopt -s nullglob
  existing_outputs=("${OUTPUT_ROOT}/2.1_"*.tsv)
  shopt -u nullglob
  if [[ "${#existing_outputs[@]}" -gt 0 ]]; then
    echo "[ERROR] output already exists, do not overwrite by default: $OUTPUT_ROOT" >&2
    printf '[ERROR] existing: %s\n' "${existing_outputs[@]}" >&2
    exit 1
  fi
fi

mkdir -p "$LOG_DIR"
MASTER_LOG="${LOG_DIR}/2.1_transcript_first_preASE_candidate_assessment.master.log"
RUN_LOG="${LOG_DIR}/2.1_transcript_first_preASE_candidate_assessment.run.log"
SUCCESS_LIST="${LOG_DIR}/2.1_transcript_first_preASE_candidate_assessment.success.tsv"
FAILED_LIST="${LOG_DIR}/2.1_transcript_first_preASE_candidate_assessment.failed.tsv"

printf 'task\toutdir\trc\trun_log\n' > "$SUCCESS_LIST"
printf 'task\toutdir\trc\trun_log\n' > "$FAILED_LIST"

CMD=(
  "$PYTHON_BIN"
  "$SCRIPT"
  --classification-tsv "$CLASSIFICATION_TSV"
  --transcript-exon-bed "$TRANSCRIPT_EXON_BED"
  --feelnc-status-tsv "$FEELNC_STATUS_TSV"
  --reference-gtf "$REFERENCE_GTF"
  --outdir "$OUTPUT_ROOT"
  --containment-overlap-fraction "$CONTAINMENT_OVERLAP_FRACTION"
  --expected-transcript-count "$EXPECTED_TRANSCRIPT_COUNT"
)

if [[ "$OVERWRITE" == "1" ]]; then
  CMD+=(--overwrite)
fi

{
  log_msg "START"
  log_msg "SCRIPT=$SCRIPT"
  log_msg "CLASSIFICATION_TSV=$CLASSIFICATION_TSV"
  log_msg "TRANSCRIPT_EXON_BED=$TRANSCRIPT_EXON_BED"
  log_msg "FEELNC_STATUS_TSV=$FEELNC_STATUS_TSV"
  log_msg "REFERENCE_GTF=$REFERENCE_GTF"
  log_msg "OUTPUT_ROOT=$OUTPUT_ROOT"
  log_msg "CONTAINMENT_OVERLAP_FRACTION=$CONTAINMENT_OVERLAP_FRACTION"
  log_msg "EXPECTED_TRANSCRIPT_COUNT=$EXPECTED_TRANSCRIPT_COUNT"
  log_msg "OVERWRITE=$OVERWRITE"
  printf '[%s ] CMD=' "$(date '+%F %T')"
  printf ' %q' "${CMD[@]}"
  printf '\n'
} > "$MASTER_LOG"

set +e
"${CMD[@]}" > "$RUN_LOG" 2>&1
RC=$?
set -e

if [[ "$RC" -eq 0 ]]; then
  printf 'transcript_first_preASE_candidate_assessment\t%s\t%s\t%s\n' "$OUTPUT_ROOT" "$RC" "$RUN_LOG" >> "$SUCCESS_LIST"
  log_msg "FINISH success rc=$RC, task completed / done" | tee -a "$MASTER_LOG"
else
  printf 'transcript_first_preASE_candidate_assessment\t%s\t%s\t%s\n' "$OUTPUT_ROOT" "$RC" "$RUN_LOG" >> "$FAILED_LIST"
  log_msg "FINISH failed rc=$RC" | tee -a "$MASTER_LOG"
fi

log_msg "RUN_LOG=$RUN_LOG" | tee -a "$MASTER_LOG"
log_msg "SUCCESS_LIST=$SUCCESS_LIST" | tee -a "$MASTER_LOG"
log_msg "FAILED_LIST=$FAILED_LIST" | tee -a "$MASTER_LOG"
exit "$RC"
