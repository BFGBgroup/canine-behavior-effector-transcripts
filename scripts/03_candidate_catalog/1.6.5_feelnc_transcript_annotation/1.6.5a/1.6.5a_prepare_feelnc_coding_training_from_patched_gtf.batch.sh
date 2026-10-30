#!/bin/bash
set -euo pipefail


SCRIPT="/path/to/repository/scripts/03_candidate_catalog/1.6.5_feelnc_transcript_annotation/1.6.5a/1.6.5a_prepare_feelnc_coding_training_from_patched_gtf.py"
PYTHON_BIN="/usr/bin/python3"

PATCHED_GTF="/path/to/project/legacy/five-dogs/5.cnv/0.gtf/GCF_000002285.3_CanFam3.1_genomic_chr.gene_id_patched.gtf"
OUTPUT_ROOT="/path/to/project/legacy/five-dogs/9.new-gene/1.6.5_feelnc_transcript_annotation/1.6.5a"
TRAINING_DIR="${OUTPUT_ROOT}/training_sets"
LOG_DIR="${OUTPUT_ROOT}/logs-v1"

OUT_TRANSCRIPT_IDS="${TRAINING_DIR}/CanFam3.1_refseq_protein_coding_transcript_ids.from_patched_CDS.txt"
OUT_TRAINING_GTF="${TRAINING_DIR}/CanFam3.1_refseq_protein_coding_transcripts.from_patched_CDS.gtf"
QC_SUMMARY="${TRAINING_DIR}/CanFam3.1_refseq_protein_coding_training_from_patched_CDS.qc_summary.tsv"
MANIFEST="${TRAINING_DIR}/CanFam3.1_refseq_protein_coding_training_from_patched_CDS.manifest.tsv"

OVERWRITE=0
MIN_CODING_TRANSCRIPTS=1000

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

if [[ ! -s "$PATCHED_GTF" ]]; then
  echo "[ERROR] patched GTF missing or empty: $PATCHED_GTF" >&2
  exit 1
fi

mkdir -p "$TRAINING_DIR" "$LOG_DIR"

MASTER_LOG="${LOG_DIR}/1.6.5a_prepare_feelnc_coding_training.batch.master.log"
RUN_LOG="${LOG_DIR}/1.6.5a_prepare_feelnc_coding_training.run.log"
SUCCESS_LIST="${LOG_DIR}/1.6.5a_prepare_feelnc_coding_training.success.tsv"
FAILED_LIST="${LOG_DIR}/1.6.5a_prepare_feelnc_coding_training.failed.tsv"

printf "task\toutdir\trc\trun_log\n" > "$SUCCESS_LIST"
printf "task\toutdir\trc\trun_log\n" > "$FAILED_LIST"

CMD=(
  "$PYTHON_BIN"
  "$SCRIPT"
  --patched-gtf "$PATCHED_GTF"
  --out-training-gtf "$OUT_TRAINING_GTF"
  --out-transcript-ids "$OUT_TRANSCRIPT_IDS"
  --qc-summary "$QC_SUMMARY"
  --manifest "$MANIFEST"
  --min-coding-transcripts "$MIN_CODING_TRANSCRIPTS"
)

if [[ "$OVERWRITE" == "1" ]]; then
  CMD+=(--overwrite)
fi

{
  echo "[batch] $(date '+%F %T') START"
  echo "[batch] SCRIPT=$SCRIPT"
  echo "[batch] PYTHON_BIN=$PYTHON_BIN"
  echo "[batch] PATCHED_GTF=$PATCHED_GTF"
  echo "[batch] OUTPUT_ROOT=$OUTPUT_ROOT"
  echo "[batch] TRAINING_DIR=$TRAINING_DIR"
  echo "[batch] OUT_TRANSCRIPT_IDS=$OUT_TRANSCRIPT_IDS"
  echo "[batch] OUT_TRAINING_GTF=$OUT_TRAINING_GTF"
  echo "[batch] QC_SUMMARY=$QC_SUMMARY"
  echo "[batch] MANIFEST=$MANIFEST"
  echo "[batch] LOG_DIR=$LOG_DIR"
  echo "[batch] OVERWRITE=$OVERWRITE"
  echo "[batch] MIN_CODING_TRANSCRIPTS=$MIN_CODING_TRANSCRIPTS"
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
  printf "prepare_feelnc_coding_training_1.6.5a\t%s\t%s\t%s\n" "$TRAINING_DIR" "$RC" "$RUN_LOG" >> "$SUCCESS_LIST"
  echo "[batch] $(date '+%F %T') FINISH success rc=$RC" | tee -a "$MASTER_LOG"
  echo "[batch] OUT_TRAINING_GTF=$OUT_TRAINING_GTF" | tee -a "$MASTER_LOG"
  echo "[batch] OUT_TRANSCRIPT_IDS=$OUT_TRANSCRIPT_IDS" | tee -a "$MASTER_LOG"
  echo "[batch] QC_SUMMARY=$QC_SUMMARY" | tee -a "$MASTER_LOG"
else
  printf "prepare_feelnc_coding_training_1.6.5a\t%s\t%s\t%s\n" "$TRAINING_DIR" "$RC" "$RUN_LOG" >> "$FAILED_LIST"
  echo "[batch] $(date '+%F %T') FINISH failed rc=$RC" | tee -a "$MASTER_LOG"
  echo "[batch] RUN_LOG=$RUN_LOG" | tee -a "$MASTER_LOG"
fi

exit "$RC"
